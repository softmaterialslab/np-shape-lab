// This is vertex class

#ifndef _EDGE_H
#define _EDGE_H

#include "utility.h"
#include "vector3d.h"
#include "thermostat.h"
#include "interface.h"
#include <map>
#include "vertex.h"
#include "face.h"

class VERTEX;
class FACE;
class INTERFACE;

class EDGE 
{
  public:

  // members
  unsigned int index;
  vector<VERTEX*> itsV;			// two vertices define an edge
  vector<FACE*> itsF;
  double itsS;
  int isOnBoundary;
  
  // member function
  VERTEX* opposite(VERTEX*);
  FACE*   opposite(FACE*);
  double lengthSquared();
  double length();
  double l0;

  double compute_S();
  void     bending_forces(INTERFACE* boundary);
  VECTOR3D compute_gradS(VERTEX* wrt);
  double crossingLengthSquared();
  void    flipIfFavorable();
  void   flip();
};


#endif

