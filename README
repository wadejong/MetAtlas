Python code and data to start generating protonation and deprotonation energies.
--------------------------------------------------------------------------------
1. Generate data for data sets.

2. Once we have that data, we can correlate this through ML to structural features.
   Want to use the Harvard autoencoder. It is based on sentence recognition, but each work is 
   representing one letter in smiles string. I want to look at defining chemical groups as words.

3. Try and predict protonation/deprotonation sites and energies of molecules outside the data set.


Second part of project is to predict fragmentation spectra. 
This is a step up from protonation/deprotonation.
-----------------------------------------------------------
1. Generate energy costs for fragmenting bonds in molecules. We'll need to do a subset as this is combinatorial.

2. Correlate autoencoded structures (with and without protonation information) to experimentally available mass spectra.

3. From a structure outside the testset, predict a mass spectrum. i

4. If we can do this, we can start to enumerate predicted mass spectra and we can train a neural net to assign
   experimental spectra.


Files for first part:
---------------------
metatlas_inchi_inchikey.csv - Inchi data from LBNL's Ben Bowen
metatlas.openbabel.py - Main code using OpenBabel
metatlas.rdkit.py - Main code using RDKit
pubchempy - Notes on investigating using PubChem to generate 3d structure (can only do 2d through REST API)

Tools to be installed: OpenBabel, RDKit (needs Conda or MiniConda), Orca (open-source download)


Files for second part:
----------------------
EMA_0p1.txt - Inchi data set from Harvard group used for autoencoder
SmilesData* - Smiles generated from EMA with OpenBabel or RDKit (latter also generated canonical smiles)
inchitosmiles* - Python for converting Inchi to smiles for both data sets and both toolkits
