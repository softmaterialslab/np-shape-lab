# Nanoparticle Shape Lab

* load the necessary modules; module load gsl

* make the code by ```make```; then ```make clean```

* run the simulation for the following set of parameters using ```echo 900 0.01 3 30 30 | ./oncluster_simulate_membrane```. These are charge, salt_conc, surface tension, bending, stretching.

* quick check the results by comparing the following in outfiles/ folder
        * energy_nanomembrane.dat with new energy file. col 2 is kinetic, col 3 is potential, and col 5 is total. Total is conserved outside of annealing zones
        * compare area.dat (col 2) and volume.dat (col 2) with new files.
        * check simulation.out for the nanoparticle system configuration
