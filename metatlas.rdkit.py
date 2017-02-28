#
# Base code for calculating most likely protonation sites
#
# Using RDKit as the toolkit to generate protonated and deprotonated structures
#

import csv, subprocess, pymongo, bson, uuid
from rdkit import Chem
from rdkit.Chem import AllChem, rdqueries, Descriptors
  
isHydrogen = rdqueries.AtomNumEqualsQueryAtom(1)
isUnderSaturated = rdqueries.IsUnsaturatedQueryAtom()

#
# Goal is to store the data in database, here trying MongoDb
#

mongoClient = pymongo.MongoClient()
mongoDB = mongoClient.test    
mongoDocument = {}


def mySplit(myString, numberOfVariables):
  return (myString.split() + [None] * numberOfVariables)[:numberOfVariables]


def numberOfElectrons(myMolecule):
  tbl = Chem.GetPeriodicTable()
  return(sum(atom.GetAtomicNum() for atom in myMolecule.GetAtoms()))

#
# From structure and charge, generate an Orca input
#

def createInput(myMolecule,myCharge):
  myMult = (1 if (numberOfElectrons(myMolecule)+myCharge) % 2 == 0 else 2)
  inputFileName = str(uuid.uuid4())+'.inp'
  inputFile = open(inputFileName,'w')
  inputFile.write('!PM3 Opt \n%coords \n  CTyp xyz\n')
  inputFile.write(' Charge '+str(myCharge)+'\n Mult '+str(myMult)+' \n coords\n')
  for myAtom in myMolecule.GetAtoms():
   index = myAtom.GetIdx()
   inputFile.write('  '+myAtom.GetSymbol()+' '+str(myMolecule.GetConformer().GetAtomPosition(index).x)+' '+str(myMolecule.GetConformer().GetAtomPosition(index).y)+' '+str(myMolecule.GetConformer().GetAtomPosition(index).z)+' \n')
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

def getEnergy(myMolecule,myCharge):
  mongoDocument = {}
  myInputFileName = createInput(myMolecule,myCharge)
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
               myEnergy = float(outputLine.split()[3])
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
  mongoDocument['parentInchi'] = myMolecule.GetProp('inchiString')
  mongoDocument['childInchi'] = AllChem.MolToInchi(myMolecule)
  mongoDocument['childSmiles'] = AllChem.MolToSmiles(myMolecule)
  mongoDocument['molecularFormula'] = AllChem.CalcMolFormula(myMolecule)
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
 csvRow = next(csvReader)
 for csvRow in csvReader:
   speciesNumber,inchiString,inchiKey = csvRow[0],csvRow[1],csvRow[2]
   myMolecule = AllChem.AddHs(AllChem.MolFromInchi(inchiString))
   myMolecule.SetProp('inchiString',inchiString)
   AllChem.EmbedMolecule(myMolecule)
   AllChem.UFFOptimizeMolecule(myMolecule)
   neutralEnergy = getEnergy(myMolecule,AllChem.GetFormalCharge(myMolecule))
   protonatedEnergy, deprotonatedEnergy = 0.0, 0.0
   protonatedAtom, deprotonatedAtom = 0, 0
#
# Deprotonate
#
   for myAtom in myMolecule.GetAtomsMatchingQuery(isHydrogen):
     tempMol = AllChem.EditableMol(myMolecule) 
     tempMol.RemoveAtom(myAtom.GetIdx())
     modifiedMolecule = tempMol.GetMol()
     calculatedEnergy = getEnergy(modifiedMolecule,AllChem.GetFormalCharge(myMolecule)-1)
     if calculatedEnergy < deprotonatedEnergy:
       deprotonatedEnergy = calculatedEnergy
       deprotonatedAtom = myAtom.GetBonds()[0]
     print('deprotonation',calculatedEnergy)
#
# Protonate, still needs work to make this go
#
   for myAtom is myMolecule.GetAtomsMatchingQuery(isUnderSaturated):
     tempMol = Chem.EditableMol(myMolecule)
     tempMol.GetAtomWithIdx(myAtom.GetIdx()).SetFormalCharge(1)
     AllChem.AddHs(tempMol,addCoords=True)
     modifiedMolecule = tempMol.GetMol()
     calculatedEnergy = getEnergy(modifiedMolecule,AllChem.GetFormalCharge(modifiedMolecule))
     if calculatedEnergy < protonationEnergy:
       protonationEnergy = calculatedEnergy
       protonatedAtom = myAtom.GetIdx()
      print('protonation',calculatedEnergy)
   print('lowest (de)protonation',deprotonatedEnergy,protonatedEnergy)
