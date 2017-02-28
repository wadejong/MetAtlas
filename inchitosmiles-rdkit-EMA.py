#
# Using RDKit to convert inchi strings to smiles used by Harvard group
#

import csv
from rdkit.Chem import AllChem
  
csvfileOut = open('SmilesData', 'w')
speciesNumber = 0

#
# Reading from Harvard EMA set 
#

with open('EMA_0p1.txt',newline='') as csvFile:
 csvReader = csv.reader(csvFile)
 for csvRow in csvReader:
   inchiString = csvRow[0]
   speciesNumber += 1
   print('Entry',speciesNumber,inchiString)
   myMolecule = AllChem.MolFromInchi(inchiString)
#  mySmiles  = AllChem.MolToSmiles(myMolecule, isomericSmiles=True, canonical=False )
   mySmiles  = AllChem.MolToSmiles(myMolecule)
   csvfileOut.write(mySmiles+'  '+str(speciesNumber)+'\n')
