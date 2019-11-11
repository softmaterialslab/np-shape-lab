//
// Created by nbrunk on 10/16/2019.
//

// This is particle class

#ifndef _PARTICLE_H
#define _PARTICLE_H

#include "utility.h"
#include "vector3d.h"
#include "thermostat.h"

class PARTICLE
{
private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        ar & id;
        ar & diameter;
        ar & valency;
        ar & q;
        ar & m;
        ar & epsilon;
        ar & posvec;
        ar & velvec;
        ar & forvec;
        ar & pe;
        ar & electrostaticPE;
        ar & ke;
        ar & energy;
    }

public:

    // members
    int id;		// id of the particle
    double diameter;	// diameter of the particle
    int valency;		// valency of the ion
    double q;		// charge of the particle
    double m; 		// mass of the particle
    double epsilon;	// dielectric constant of the medium
    VECTOR3D posvec;	// position vector of the particle
    VECTOR3D velvec;	// velocity vector of the particle
    VECTOR3D forvec;	// force vector on the particle
    double pe;		// potential energy
    double electrostaticPE;  // the ES component of the PE
    long double ke;	// kinetic energy
    double energy;	// energy

    // member functions

    // make a particle
    PARTICLE(int get_id = 0, double get_diameter = 0, int get_valency = 0, double get_charge = 0, double get_mass = 0, VECTOR3D get_position = VECTOR3D(0,0,0))
    {
        id = get_id;
        diameter = get_diameter;
        valency = get_valency;
        q = get_charge;
        m = get_mass;
        posvec = get_position;
    }

    // update position of the particle
    void update_real_position(double dt) {
        posvec = (posvec + (velvec ^ dt));
        return;
    }

    // update velocity of the particle
    void update_velocity(double dt)
    {
        velvec = ( velvec + ( forvec ^ ( 0.5 * dt ) ) );
        return;
    }

    void update_real_velocity(double dt, THERMOSTAT ion_bath, long double expfac)
    {
        velvec = ((velvec ^ (expfac)) + (forvec ^ (0.5 * dt * sqrt(expfac))));
        return;
    }

    // calculate kinetic energy of a particle
    void real_kinetic_energy()
    {
        ke = 0.5 * m * (velvec * velvec);
        return;
    }
};

#endif