# Nanoparticle Shape Lab

## Install instructions on BigRed2
* First, git clone the project:
```git clone https://github.com/softmaterialslab/np-shape-lab.git```
* Then, load the required modules using following command:
```module swap PrgEnv-cray PrgEnv-gnu && module load boost/1.65.0 && module load gsl```
* Next, go to the root directory:
 ```cd np-shape-lab```
* Then, install the project:
```make cluster-install```
* Next, submit a test job:
```make cluster-test-submit```
* Then, clean the datafiles from the test job:
```make dataclean```
* Fianlly, submit the job:
```makecluster-submit```
* All outputs from the simulation will be stored in the bin folder when the simulation is completed.
* Check and compare files (ex: energy_nanomembrane.dat) inside the ```bin/outfiles``` directory.

## Install instructions on Local computer
* Load the necessary modules:
```module load gsl && module load openmpi/3.0.1 && module load boost/1_67_0```
* Next, go to the root directory:
 ```cd np-shape-lab```
* Then, install the project:
```make all```
