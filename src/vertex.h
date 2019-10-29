// This is vertex class

#ifndef _VERTEX_H
#define _VERTEX_H

#include "utility.h"
#include "vector3d.h"
#include "thermostat.h"
#include "edge.h"
#include "face.h"
#include "interface.h"

class EDGE;

class FACE;

class INTERFACE;

class VERTEX {

private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive &ar, const unsigned int version) {
        ar & index;
        ar & m;
        ar & q;
        ar & posvec;
        ar & velvec;
        ar & itsE;
        ar & itsF;
        ar & itsnormal;
        ar & itsarea;
        ar & itsrefarea;
        ar & aforce;
        ar & vforce;
        ar & lforce;
        ar & forvec;
        ar & bforce;
        ar & sforce;
        ar & TForce;
        ar & VolTForce; // NB added
        ar & neg_GradVi;
        ar & neg_GradAi;
        ar & neg_GradConi;
        ar & ke;
        ar & r;
        ar & theta;
        ar & phi;
    }

public:

    // members
    unsigned int index;            // index
    double m;                // mass
    double q;                                // charge
    VECTOR3D posvec;            // position vector
    VECTOR3D velvec;            // velocity vector
    vector<EDGE *> itsE;            // its edges
    vector<FACE *> itsF;            // its faces

    VECTOR3D itsnormal;
    double itsarea;
    double itsrefarea;

    VECTOR3D aforce;
    VECTOR3D vforce;
    VECTOR3D lforce;
    VECTOR3D forvec;
    VECTOR3D bforce;
    VECTOR3D sforce;
    VECTOR3D TForce;
    VECTOR3D VolTForce;

    VECTOR3D neg_GradVi;      // NB modified.  Precise name for what was previously "dummy_vforce".
    VECTOR3D neg_GradAi;      // NB added.
    VECTOR3D neg_GradConi;    // NB added.  Neg constraint gradient on vertex, differs from above if quadratic.

    double ke;
    double r, theta, phi;            //  Returns the polar coordinates of the vertex.

    void compute_area_normal();

    void computearea();

    void
    compute_neg_VGrads(INTERFACE */*, char*/); // 2017.09.16 NB added:  Compute given vertex net area & volume gradient.
    void stretching_forces(INTERFACE *, char);

    void tension_forces(INTERFACE *);           // 2017.09.16 NB added:  Compute the tension forces on a given vertex.
    void volume_tension_forces(INTERFACE *);    // 2019.08.16 NB added:  Comptue volume tension forces on a given vertex.

    //void constraint_forces(INTERFACE*);
    //void volume_constraint_gradient(INTERFACE*);

    // member functions

    // make a vertex
    VERTEX(VECTOR3D position = VECTOR3D(0, 0, 0)) : posvec(position) {
    }

    // update velocity of the vertex
    void update_real_velocity(double dt, THERMOSTAT main_bath, long double expfac) {
        velvec = ((velvec ^ (expfac)) + (forvec ^ (0.5 * dt * sqrt(expfac))));
        return;
    }

    // compute kinetic energy
    void real_kinetic_energy() {
        ke = 0.5 * m * (velvec * velvec);
        return;
    }

    // update position of the vertex
    void update_real_position(double dt) {
        posvec = (posvec + (velvec ^ dt));
        return;
    }

    /*Unused Functions:
    // get the polar coordinates for the vertex
    void get_polar() {
        r = posvec.GetMagnitude();
        theta = acos(posvec.z / r);
        phi = posvec.y > 0 ? acos(posvec.x / sqrt(r * r - posvec.z * posvec.z)) : 2 * pi - acos(posvec.x / sqrt(r * r -
                                                                                                                posvec.z *
                                                                                                                posvec.z));
        return;
    }*/
};

#endif
