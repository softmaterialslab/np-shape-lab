// This is a header file to set control

#ifndef _CONTROL_H
#define _CONTROL_H

class CONTROL {
private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive &ar, const unsigned int version) {
        ar & timestep;
        ar & steps;
        ar & anneal;
        ar & annealfreq;
        ar & annealDuration;
        ar & TAnnealFac;
        ar & QAnnealFac;
        ar & moviestart;
        ar & offfreq;
        ar & moviefreq;
        ar & povfreq;
        ar & writedata;
    }

public:
    double timestep;        // timestep used in molecular dynamics
    int steps;                // number of steps in molecular dynamics
    char anneal;            // if you want to anneal, applies only to fmd
    int annealfreq;          // frequency of annealing
    int annealDuration; // NB added to anneal gradually over this many steps.  Value ascribed in main.
    double TAnnealFac;  // NB added.  Fractional decrement to temperature per annealing procedure.
    double QAnnealFac;  // NB added.  Same for the thermostat fake mass.
    
    //double stepStart;   // NB added.  Dictates the step at which the MD should start.  Useful if loading off-file.

    int moviestart;         //movie strat point
    int offfreq;            // wait till this step
    int moviefreq;                // frequency of sampling
    int povfreq;    // energy computed after these many steps
    int writedata;          // write the data files
};

#endif

