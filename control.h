// This is a header file to set control

#ifndef _CONTROL_H
#define _CONTROL_H

class CONTROL
{
  public:
    double timestep;		// timestep used in molecular dynamics
    int steps;			    // number of steps in molecular dynamics
    char anneal;		    // if you want to anneal, applies only to fmd
    int annealfreq;		  // frequency of annealing
    int annealDuration; // NB added to anneal gradually over this many steps.  Value ascribed in main.
    double TAnnealFac;  // NB added.  Fractional decrement to temperature per annealing procedure.
    double QAnnealFac;  // NB added.  Same for the thermostat fake mass.\
    double stepStart;   // NB added.  Dictates the step at which the MD should start.  Useful if loading off-file.

    double fakemass;		// fictitious mass used in molecular dynamics
    int hiteqm; 		    // wait till this step
    int freq; 			    // frequency of sampling
    char on;			      // y if fmd needs to turn on
    int extra_compute;	// energy computed after these many steps
    int verify;			    // verification with correct answer for ind. density done after these many steps
    int writeverify;		// write the verification files
    int writedensity; 	// write the density files
    int writedata; 		  // write the data files
};

#endif

