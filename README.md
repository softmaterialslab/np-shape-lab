# Nanoparticle Shape Lab

## Install instructions on BigRed3
* First, git clone the project:
```git clone https://github.com/softmaterialslab/np-shape-lab.git```
* Then, load the required modules using following command:
```module swap PrgEnv-cray PrgEnv-gnu && module load boost/1.65.0 && module load gsl```
* Next, go to the root directory:
 ```cd np-shape-lab```
* Then, install the project:
```make cluster-install```
* Finally, submit the job:
```make cluster-submit```
* All outputs from the simulation will be stored in the bin folder when the simulation is completed.
* Check and compare files (ex: energy_nanomembrane.dat) inside the ```bin/outfiles``` directory.
* Clean the datafiles if desired:
```make dataclean```

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
```time ./np_shape_lab -R 10 -q 600 -c 0.005 -t 1 -v 1 -b 40 -s 40 -S 25000 -D 4 -F n```
* Respectively, these are the (radius (in nm), net charge (in e), salt concentration (in Molar), surface tension (in dyn/cm), volume tension (in atm/nm^3), bending rigidity (in kBT), stretching rigidity (in kBT), net number of steps, and discretization parameter).
* After a few (3 - 10) minutes, this should produce a disc of final reduced area (A = 15.675), local potential (U = 2620.98 kB T), and conserved total energy (E = 2699.97 kB T).
* Note that minor changes on order of a percent are expected due to shuffling of the initial charge distributions dependent on different machines' random seed.
* For testing a higher resolution grid (D = 8), use:
```time ./np_shape_lab -R 10 -q 600 -c 0.005 -t 1 -v 1 -b 40 -s 40 -S 25000 -D 8 -F n```

### Homogeneously-charged Rod Formation:
* Simply increasing the salt concentration (c) allows for testing homogeneously charged rod formation using the command below:
```time ./np_shape_lab -R 10 -q 600 -c 0.01 -t 1 -v 1 -b 40 -s 40 -S 25000 -D 4 -F n```
* After a few minutes, this should produce a rod of final reduced area (A = 13.259), local potential (U = 1919.14 kB T), and conserved total energy (E = 1939.96 kB T).

### Inhomogeneously-charged Hemisphere Formation:
* The same parameters as in the disc example above may be used to test hemisphere formation, with two added parameters (N) and (p):
```time ./np_shape_lab -R 10 -q 600 -N 2 -p 0.5 -c 0.005 -t 1 -v 1 -b 40 -s 40 -S 25000 -D 4 -F n```
* Respectively, the new parameters specify the collective number of stripes (N) and the fractional area of the charged patch (p), such that (p = 0.5) is a standard Janus particle.  The charge (q) now specifies the charge were it homogeneously charged, effectively specifying a charge density in the charged region.
* After a few minutes, this should produce a hemisphere of unchanged final reduced area (A = 12.53), local potential (U = 1224 kB T), and conserved total energy (E = 1437).

###  Uncharged, Icosahedrally Buckled Control:
* The above three examples use the default setting that does not induce spontaneous elastic buckling.:
```time ./np_shape_lab -R 10 -q 0 -c 0.005 -t 0 -v 0 -b 1 -s 1000 -S 25000 -D 4 -F n -B y```
* The added parameter is a buckling flag ("-B").
* After a few (3 - 10) minutes, this should produce an icosahedron of final reduced area (A = 12.3795), local potential (U = 26.86 kB T), and conserved total energy (E = 127.69 kB T).

###  Yin-yang Patterns:
* The same parameters as in the Inhomogeneously-charged Hemisphere example, with one additional parameter (H):
```time ./np_shape_lab -R 10 -q 600 -N 2 -p 0.5 -c 0.005 -t 1 -v 1 -b 40 -s 40 -S 25000 -D 4 -F n -H y```
* The added parameter is a yinyang function flag ("-H").
* After a few minutes, this should produce an yinyang-sphere of unchanged final reduced area (A = 12.5108), local potential (U = 1218.2 kB T), and conserved total energy (E = 1433.21 kB T).