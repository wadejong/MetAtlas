#
# Base code for calculating most likely protonation sites
#
# Using OpenBabel as the toolkit to generate protonated and deprotonated structures
#

import openbabel, uuid, csv, subprocess, pymongo, bson
  
obConversion = openbabel.OBConversion()
obConversion.SetInAndOutFormats('inchi','xyz')
obBuilder = openbabel.OBBuilder()
obMolecule = openbabel.OBMol()
obElementTable = openbabel.OBElementTable()
obForceField = openbabel.OBForceField.FindForceField('UFF')

#
# Goal is to store the data in database, here trying MongoDb
#

mongoClient = pymongo.MongoClient()
mongoDB = mongoClient.test    
mongoDocument = {}


def mySplit(myString, numberOfVariables):
  return (myString.split() + [None] * numberOfVariables)[:numberOfVariables]


def numberOfElectrons(myMolecule):
  electronCount = 0
  for myAtom in openbabel.OBMolAtomIter(myMolecule):
    electronCount += myAtom.GetAtomicNum()
  return(electronCount)

#
# From structure and charge, generate an Orca input
#

def createInput(myMolecule):
  myCharge = myMolecule.GetTotalCharge()
  myMult = (1 if (numberOfElectrons(myMolecule)+myCharge) % 2 == 0 else 2)
  inputFileName = str(uuid.uuid4())+'.inp'
  inputFile = open(inputFileName,'w')
  inputFile.write('!PM3 Opt \n%coords \n  CTyp xyz\n')
  inputFile.write(' Charge '+str(myCharge)+'\n Mult '+str(myMult)+' \n coords\n')
  for myAtom in openbabel.OBMolAtomIter(myMolecule):
    inputFile.write('  '+obElementTable.GetSymbol(myAtom.GetAtomicNum())+' '+str(myAtom.GetX())+' '+str(myAtom.GetY())+' '+str(myAtom.GetZ())+' \n')
  inputFile.write(' end\nend\n')
  inputFile.write('%geom\n MaxIter 200\n end\n')
  inputFile.close()
  mongoDocument['calculationSetup'] = { "molecularSpinMultiplicity" : myMult, "charge" : myCharge, "numberOfElectrons" : numberOfElectrons(myMolecule), "waveFunctionTheory" : "PM3" } 
  return(inputFileName)

#
# Run Orca with generated input
#

def calculateEnergy(myInputFileName):
  myOutputFileName = myInputFileName.strip('inp')+'out'
  myOutput = open(myOutputFileName,'wb')
  subprocess.call(['../../orca_3_0_3/orca',myInputFileName],stdout=myOutput)
  return(myOutputFileName)

#
# Get energy for given structure and charge
#

def getEnergy(myMolecule):
  mongoDocument = {}
  myInputFileName = createInput(myMolecule)
  myOutputFileName = calculateEnergy(myInputFileName)
  myOutput = open(myOutputFileName,'r')
  outputLine = myOutput.readline()
  moleculeObj = {}
  atomsObj = []
  myEnergy = 0.0
  while outputLine:
    if outputLine.find('THE OPTIMIZATION HAS CONVERGED')>=0:
      while outputLine: 
        if outputLine.find('CARTESIAN COORDINATES')>=0:
          outputLine = myOutput.readline()
          outputLine = myOutput.readline()
          while outputLine:
            atomObj = {}
            atomName,xCoord,yCoord,zCoord = mySplit(outputLine,4)
            if atomName == None: 
              break
            else:
              atomObj['elementSymbol'] = atomName
              atomObj['cartesianCoordinates'] = { 'value' : [xCoord, yCoord, zCoord], 'units' : 'Angstrom' }
              atomsObj.append(atomObj)
            outputLine = myOutput.readline()
          while outputLine:
            if outputLine.find('Total Energy')>=0:
               myEnergy = outputLine.split()[3]
               break
            outputLine = myOutput.readline()
          break
        outputLine = myOutput.readline()
      break
    outputLine = myOutput.readline()
#
#  Store data in MongoDb
#
  mongoDocument['_id'] = myInputFileName.strip('.inp')
  mongoDocument['parentInchi'] = obMolecule.GetTitle(inchiString)
  mongoDocument['childInchi'] = obConversion.WriteString(myMolecule)
  obConversion.SetOutFormat('smi')
  mongoDocument['childSmiles'] = obConversion.WriteString(myMolecule)
  obConversion.SetOutFormat('inchi')
  mongoDocument['molecularFormula'] = myMolecule.GetFormula()
  moleculeObj['atoms'] = atomsObj
  mongoDocument['molecule'] = moleculeObj
  mongoDocument['totalEnergy'] = {"value":myEnergy, "units":"Hartree"}
  insertResult = mongoDB.metatlas.insert_one(mongoDocument)
  return(myEnergy)

#
# Main routine reads metatlas inchi and converts to xyz
# Optimize structure
# Find deprotonation sites
# Find protonation sites
#

with open('metatlas_inchi_inchikey.csv',newline='') as csvFile:
 csvReader = csv.reader(csvFile)
 for csvRow in csvReader:
   speciesNumber,inchiString,inchiKey = csvRow[0],csvRow[1],csvRow[2]
   print(inchiString)
   obConversion.SetInAndOutFormats('smi','xyz')
   obConversion.ReadString(obMolecule,inchiString)
# I need this line to have the correct number of hydrogens
   moleculeFormula = obMolecule.GetFormula()
   print(obMolecule.GetFormula())
   print(obMolecule.NumAtoms())
   print(obMolecule.NumConformers())
   print(obConversion.WriteString(obMolecule))
   obMolecule.AddHydrogens()
   obBuilder.Build(obMolecule)
   print(obMolecule.GetFormula())
   print(obMolecule.NumAtoms())
   print(obMolecule.NumConformers())
   print(obMolecule.Has2D())
   print(obMolecule.Has3D())
   obConversion.WriteFile(obMolecule,'test.xyz')
#   obForceField.Setup(obMolecule)
#   obForceField.ConjugateGradients(1000,1.0e-6)
#   obForceField.UpdateCoordinates(obMolecule)
   obMolecule.SetTitle(inchiString)
   neutralEnergy = getEnergy(obMolecule)
   protonatedEnergy, deprotonatedEnergy = 0.0, 0.0
   protonatedAtom, deprotonatedAtom = 0, 0
#
# Deprotonate
#
   for obAtom in openbabel.OBMolAtomIter(obMolecule):
     if obAtom.IsHydrogen:
       modifiedMolecule = openbabel.OBMol(obMolecule)
       modifiedMolecule.DeleteAtom(obAtom)
       modifiedMolecule.SetTotalCharge(obMolecule.GetTotalCharge()-1)
       calculatedEnergy = getEnergy(modifiedMolecule)
       if calculatedEnergy < deprotonationEnergy:
         deprotonationEnergy = calculatedEnergy
         for connectedAtom in openbabel.OBAtomAtomIter(obAtom):
           deprotonatedAtom = connectedAtom.GetIdx()
       print('deprotonation',calculatedEnergy)
     else:
#
# Protonate, still needs work to make this go
#
       obAtom.IncrementImplicitValence()
       obMolecule.AddHydrogens(obAtom)
       totalCharge = obMolecule.GetTotalCharge()
       obMolecule.SetTotalCharge(totalCharge+1)
       calculatedEnergy = getEnergy(obMolecule)
       if calculatedEnergy < protonationEnergy:
         protonationEnergy = calculatedEnergy
         protonatedAtom = obAtom.GetIdx()
       atomToDelete = obMolecule.getAtom(obMolecule.NumAtoms())
       obMolecule.DeleteAtom(atomToDelete)
       obMolecule.SetTotalCharge(totalCharge)
       print('protonation',calculatedEnergy)
   print('lowest (de)protonation',deprotonationEnergy,protonationEnergy)
