// Face cpp file

#include "face.h"
#include "functions.h"

VERTEX* FACE::across(EDGE* curE)
{
  VERTEX* result = NULL;
  if (curE == itsE[0] || curE == itsE[1] || curE == itsE[2])
    for (unsigned int i=0; i < itsV.size(); i++)
      if (itsV[i] != curE->itsV[0] && itsV[i] != curE->itsV[1])
	result = itsV[i];
  return result;
}

void FACE::compute_area_normal()
{
  assert(itsV.size()==3);	// make sure that there are three vertices to the face
  
  //VECTOR3D r10 = itsV[1]->posvec - itsV[0]->posvec;   // NB renamed for consistency (were just a & b).
  //VECTOR3D r21 = itsV[2]->posvec - itsV[1]->posvec;
  itsdirection = (itsV[1]->posvec - itsV[0]->posvec)%(itsV[2]->posvec - itsV[1]->posvec);
  itsarea = 0.5*(itsdirection.GetMagnitude());
  itsnormal = itsdirection^(0.5/itsarea);		// normalize
  return;
}

void FACE::computearea()
{
  assert(itsV.size()==3);	// make sure that there are three vertices to the face

  //VECTOR3D r10 = itsV[1]->posvec - itsV[0]->posvec; // NB renamed for consistency (were just a & b).
  //VECTOR3D r21 = itsV[2]->posvec - itsV[1]->posvec;
  itsarea = 0.5*((itsV[1]->posvec - itsV[0]->posvec)%(itsV[2]->posvec - itsV[1]->posvec)).GetMagnitude();
  return;
}

void FACE::computevolume()
{
  assert(itsV.size()==3);	// make sure that there are three vertices to the face
  
  subtends_volume = (1./6.)*((itsV[0]->posvec)*((itsV[1]->posvec)%(itsV[2]->posvec)));
  
  return;
}

void FACE::computegradients(VERTEX* wrtv, VECTOR3D* gradarea, VECTOR3D* gradvolume)
{
  assert(itsV.size()==3);	// make sure that there are three vertices to the face
  
  VERTEX* nV;
  VERTEX* nnV;
  if (itsV[0] == wrtv)
  {
    nV = (itsV[1]);
    nnV = (itsV[2]);
  }
  else if (itsV[1] == wrtv)
  {
    nV = itsV[2];
    nnV = itsV[0];
  }
  else if (itsV[2] == wrtv)
  {
    nV = itsV[0];
    nnV = itsV[1];
  }
  else
  {
     (*gradarea) = VECTOR3D(0,0,0);
     (*gradvolume) = VECTOR3D(0,0,0);
     return;
  }
  
  (*gradvolume) = VECTOR3D((1./6.)*(nV->posvec.y*nnV->posvec.z - nV->posvec.z*nnV->posvec.y),
			  (1./6.)*(nV->posvec.z*nnV->posvec.x - nV->posvec.x*nnV->posvec.z),
			  (1./6.)*(nV->posvec.x*nnV->posvec.y - nV->posvec.y*nnV->posvec.x)); // NB verified.
  
  VECTOR3D r10 = nV->posvec - wrtv->posvec; //NB changed from v10 name to r10, consistency.
  VECTOR3D r20 = nnV->posvec - wrtv->posvec;
  VECTOR3D cp = r10%r20;
  
  double xcomp = (nV->posvec.y - nnV->posvec.y)*(cp.z) - (nV->posvec.z - nnV->posvec.z)*(cp.y);
  double ycomp = (nV->posvec.z - nnV->posvec.z)*(cp.x) - (nV->posvec.x - nnV->posvec.x)*(cp.z);
  double zcomp = (nV->posvec.x - nnV->posvec.x)*(cp.y) - (nV->posvec.y - nnV->posvec.y)*(cp.x); //NB verified.
  
  double denominator = 2 * cp.GetMagnitude(); // NB verified: 2 * |cp| = 4 * A_k .
  
  (*gradarea) = VECTOR3D(xcomp/denominator, ycomp/denominator, zcomp/denominator);
  
}

//  %%% Unused Functons: %%%

void FACE::computecenter()
{
    assert(itsV.size() == 3);
    itscenter = (itsV[0]->posvec + itsV[1]->posvec + itsV[2]->posvec)^(1.0/3.0);
}

EDGE* FACE::across(VERTEX* curV)
{
    EDGE* result = NULL;
    if (curV == itsV[0] || curV == itsV[1] || curV == itsV[2])
        for (unsigned int i=0; i < itsE.size(); i++)
            if (itsE[i]->itsV[0] != curV && itsE[i]->itsV[1] != curV)
                result = itsE[i];
    return result;
}