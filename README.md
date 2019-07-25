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
* Finally, submit the job:
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
* This should create np_shape_lab executable in the bin directory

## Testing
### Homogeneously-charged Disc Formation:
* A reference set of parameters for testing homogeneously charged disc formation are provided below in the complete executable command:
* time ./np_shape_lab -R 10 -q 600 -c 0.005 -t 1 -b 40 -s 40 -S 25000 -D 4
* Respectively, these are the (radius (in nm), net charge (in e), salt concentration (in Molar), surface tension (in dyn/cm), bending rigidity (in kBT), stretching rigidity (in kBT), net number of steps, and discretization parameter).
* After about 10 mintues, this should produce a disc of final reduced area (A = 15.69), local potential (U = 2625.28 kB T), and conserved total energy (E = 2689.97 kB T).
* Note that minor changes on order of a percent are expected due to shuffling of the initial charge distributions dependent on different machines' random seed.

### Inhomogeneously-charged Hemisphere Formation:
* The same parameters may be used to test hemisphere formation with two added parameters (N) and (p):
* time ./np_shape_lab -R 10 -q 600 -N 2 -p 0.5 -c 0.005 -t 1 -b 40 -s 40 -S 25000 -D 4
* Respectively, the new parameters specify the collective number of stripes (N) and the fractional area of the charged patch (p), such that (p = 0.5) is a standard Janus particle.  The charge (q) now specifies the charge were it homogeneously charged, effectively specifying a charge density in the charged region.
* After about 10 minutes, this should produce a hemisphere of unchanged final reduced area (A = 12.51), local potential (U = 1224.04 kB T), and conserved total energy (E = 1457.35).
