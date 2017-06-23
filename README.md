# Goals
1. Generate MetAtlas dataset containing minimized structures and energy for each
   molecule in neutral, protonated, and deprotonated states.  
2. Correlate dataset to structural features through machine learning.
    - Want to use the Harvard autoencoder. It is based on sentence recognition,
      but each work is representing one letter in smiles string. I want to look
      at defining chemical groups as words.  
3. Try and predict protonation/deprotonation sites and energies of molecules outside the data set.  
4. Predict fragmentation spectra
    1. Generate energy costs for fragmenting bonds in molecules. We'll need to
       do a subset as this is combinatorial.  
    2. Correlate autoencoded structures (with and without protonation
       information) to experimentally available mass spectra.  
    3. From a structure outside the testset, predict a mass spectrum.  
    4. If we can do this, we can start to enumerate predicted mass spectra and
       we can train a neural net to assign
       experimental spectra.

# Getting started
It's a tad non-trivial to get this whole thing setup
TODO: setup.py/requirements.txt 

## Requirements
1. [https://orcaforum.cec.mpg.de/](Orca) v3
2. [https://github.com/materialsproject/fireworks](Fireworks)
3. Python modules
    a. mendeleev
    b. pybel
4. Database config file

# Generating data
A pre-populated CSV file containing the metatlas database is included, with 
molecules stored in the 'inchi' format (metatlas_inchi_inchikey.csv). 

The [https://github.com/materialsproject/fireworks](Fireworks) framework
(developed by Jain, et. al, here at LBL), is used as an interface between the
slurm queueing system on NERSC and the the mongodb database that is used to
catalog dataset information. 

## main.py
This is the main script to read in the CSV file, process each molecule, 
and add it to the mongodb database as a Fireworks task to be run on (likely)
edison later. 

## metatlas.py
This is the application that defines special Fireworks tasks for the
optimization, collection, and addition to the database. This is required
for Fireworks to know how to submit to queue, what jobs to run, etc. 

# Database
The dataset is stored in a monogodb database created/maintained/hosted by NERSC.
I (brandon) currently am the only user with the credentials to access the
databse. Send me an email if you'd like to access.
