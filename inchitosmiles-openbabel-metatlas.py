#
# Using Openbabel to convert inchi strings to smiles used by Harvard group
# Reading MetAtlas data set
#

import csv, openbabel
  
csvfileOut = open('SmilesData', 'w')
csvFile = open('metatlas_inchi_inchikey.csv.org','r') 
obMolecule = openbabel.OBMol()
obConversion = openbabel.OBConversion()
obConversion.SetInAndOutFormats('inchi','smi')
csvReader = csv.reader(csvFile)
speciesNumber = 0
for csvRow in csvReader:
  speciesNumber,inchiString,inchiKey=csvRow[0],csvRow[1],csvRow[2]
  print('Entry',speciesNumber)
  obConversion.ReadString(obMolecule,inchiString)
  mySmiles = obConversion.WriteString(obMolecule).rstrip()
  csvfileOut.write(mySmiles+'  '+str(speciesNumber)+'\n')
