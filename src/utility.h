// This is header file for the UTILITY class.
// This file includes standard library files and gsl functions that are utilized in the code
// This file also has useful constant parameters for the problem at hand

#ifndef _UTILITY_H
#define _UTILITY_H

#include<iostream>
#include<iomanip>
#include<fstream>
#include<cmath>
#include<vector>
#include<cassert>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<gsl/gsl_sf_legendre.h>
#include<gsl/gsl_poly.h>

//OPENMP
#include <omp.h>
//BOOST MPI
#include <boost/mpi.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/base_object.hpp>

namespace mpi = boost::mpi;

using namespace std;



// constants
const double pi = 3.1415926535897932384626433832795028841971693993751;		// Pi, NB added to significant digits
const double dcut = 1.12246204830937298142;					// Cutoff distance in Lennard-Jones = 2 ^ (1/6)
const double dcut2 = 1.25992104989487316476;
const double lB_water = 0.714;							// Bjerrum length in water, in nanometers
const double epsilon_water = 80;						// Dielectric constant of water
const double room_temperature = 298;						// Room temperature in Kelvin

//const double unit_radius_sphere = 10;						// NB moved to specification in main.cpp and forwarded as arguments
const double unitenergy = 1.3807 * pow(10.0,-16) * room_temperature;		// Unit of energy (thermal energy at room temperature in CGS)
const double unitmass = 23 * 1.67 * pow(10.0, -24);				// Unit of mass (mass of sodium ion in CGS)
const double kB = 1;								// Boltzmann constant in reduced units

//const double scalefactor = epsilon_water * lB_water / unit_radius_sphere;	//NOTE scaling of electrostatic interaction found erroneous (below) before

// const double scalefactor = lB_water / unit_radius_sphere;			// Reduced units lead to this scale factor for Coloumb interaction

class UTILITY
{
public:
    long seed;
    const gsl_rng_type * T;
    gsl_rng * r;

    UTILITY()
    {
        gsl_rng_env_setup();
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);

//       srand((time(0)));                	// srand & time are built-in
//       unsigned long int s = random();  	// gsl_rng_uniform will eventually
//       gsl_rng_set(r,s); 		// seed the random number generator;

    }

    ~UTILITY()
    {
        gsl_rng_free(r);
    }
};

#endif