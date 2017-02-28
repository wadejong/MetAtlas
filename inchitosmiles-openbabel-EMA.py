#
# Using Openbabel to convert inchi strings to smiles used by Harvard group
# Reading Harvard data set
#

import csv, openbabel
  
csvfileOut = open('SmilesData', 'w')
csvFile = open('EMA_0p1.txt','r') 
obMolecule = openbabel.OBMol()
obConversion = openbabel.OBConversion()
obConversion.SetInAndOutFormats('inchi','smi')
csvReader = csv.reader(csvFile)
speciesNumber = 0
for csvRow in csvReader:
  inchiString = csvRow[0]
  speciesNumber += 1
  print('Entry',speciesNumber)
  obConversion.ReadString(obMolecule,inchiString)
  mySmiles = obConversion.WriteString(obMolecule).rstrip()
  csvfileOut.write(mySmiles+'  '+str(speciesNumber)+'\n')
