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
A conda environment can be made using the provided env.yaml file. It has 
a lot of dependencies due to bloat-y modules such as openbabel and psi4. The 
primary mode of package installation is conda, with pip being the fallback 
in cases where no conda package exists.

## Requirements
1. [Orca](https://orcaforum.cec.mpg.de/) v3
2. [Fireworks](https://github.com/materialsproject/fireworks)
3. Python modules
    a. mendeleev
    b. pybel
4. Database config file

# Generating data
A pre-populated CSV file containing the metatlas database is included, with 
molecules stored in the 'inchi' format (metatlas_inchi_inchikey.csv). 

The [Fireworks](https://github.com/materialsproject/fireworks)framework
(developed by Jain, et. al, here at LBL), is used as an interface between the
slurm queueing system on NERSC and the the mongodb database that is used to
catalog dataset information.

There is currently a mongodb instance managed and running on NERSC; I have the
credentials and can share them with whomever requires them.

## Defining Fireworks Tasks
In order to create a firework task that runs a simulation of some kind, one must
define a new _class_ that derives from FiretaskBase and has a method _run_task_
associated with it. Example tasks can be seen in `metatlas.py`. `metatlas.py**
contains all of the firetasks that have been written thus far.

## Using the newly-defined Fireworks Tasks
**In order for the conda-installed fireworks package to use your newly defined
Firetasks** that are defined in metatlas.py, you must make a soft link in the
package's user_objects directory:

```
cd /path/to/anaconda/lib/python2.7/site-packages/fireworks/user_objects
ln -s /path/to/MetAtlas/metatlas.py 
```

Once a new Task has been defined, the next thing is to get them into the
Fireworks mongodb via the queueing system. This requires a simple script that
reads in a new set of molecules (likely in SMILES format), converting them to an
input that can be fed into a Task, and then launched into the queue. `main.py`
has a simple example of doing this.
# Running all the calculations on Edison

Turns out this is kinda non-trivial and I've clearly forgotten how to do it. So
let me refresh my memory and keep track of what's going on so this doesn't
happen again.

## Actually submitting jobs to the queue
From the MOM node of a compute facility with a queuing system, I basically just
left a running 'qlaunch' instance that continually keeps the queue with m number
of jobs.

```qlaunch rapidfire -m 50 --nlaunches infinite```

I'm not sure this is sufficient because this doesn't  necessarily make sure the
calculations are running in the correct place.

Currently the configuration does look like it works. I'm able to successfully
submit a job to the SLURM queue, it creates a directory, waits, gets resources,
and attempts to run but doesn't seem to actually do the calculation
