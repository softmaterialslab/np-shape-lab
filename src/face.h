// This is vertex class

#ifndef _FACE_H
#define _FACE_H

#include "utility.h"
#include "vector3d.h"
#include "thermostat.h"
#include "vertex.h"
#include "edge.h"

class VERTEX;

class EDGE;

class FACE {

private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive &ar, const unsigned int version) {
        ar & index;
        ar & itsV;
        ar & itsE;
        ar & itsnormal;
        ar & itsdirection;
        ar & itscenter;
        ar & itsarea;
        ar & subtends_volume;
    }

public:

    // members
    unsigned int index;
    vector<VERTEX *> itsV;            // three vertices define a face
    vector<EDGE *> itsE;            // three edges define a face

    VECTOR3D itsnormal;            // normal to the face
    VECTOR3D itsdirection;            // normal to the face
    VECTOR3D itscenter;            // center of the face
    double itsarea;
    double subtends_volume;

    // member functions
    EDGE *across(VERTEX *);

    VERTEX *across(EDGE *);

    void compute_area_normal();

    void computecenter();

    void computearea();

    void computevolume();

    void computegradients(VERTEX *, VECTOR3D *, VECTOR3D *);
};

#endif
