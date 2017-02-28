#
# Using RDKit to convert inchi strings to smiles used by Harvard group
#

import csv
from rdkit.Chem import AllChem
  
csvfileOut = open('SmilesData', 'w')
speciesNumber = 0

#
# Reading from MetAtlas file
#

with open('metatlas_inchi_inchikey.csv.org',newline='') as csvFile:
 csvReader = csv.reader(csvFile)
 for csvRow in csvReader:
   speciesNumber,inchiString,inchiKey = csvRow[0],csvRow[1],csvRow[2]
   print('Entry',speciesNumber,inchiString)
   myMolecule = AllChem.MolFromInchi(inchiString)
#  mySmiles  = AllChem.MolToSmiles(myMolecule, isomericSmiles=True, canonical=False )
   mySmiles  = AllChem.MolToSmiles(myMolecule)
   csvfileOut.write(mySmiles+'  '+str(speciesNumber)+'\n')
