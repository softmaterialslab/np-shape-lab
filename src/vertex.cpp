// Vertex cpp file

#include "vertex.h"
#include "functions.h"

void VERTEX::computenormal()
{
  itsnormal = VECTOR3D(0,0,0);
  for (unsigned int i = 0; i < itsF.size(); i++)
    itsnormal += ((itsF[i]->itsnormal)^(itsF[i]->itsarea));	// area weighted
  itsnormal = (itsnormal^(1/itsnormal.GetMagnitude()));
}

//  Compute the area associated with each vertex (1/3 the sum of its faces):
void VERTEX::computearea()
{
  itsarea = 0;
  for (unsigned int i = 0; i < itsF.size(); i++)
    itsarea = itsarea + (itsF[i]->itsarea);
  itsarea = itsarea * (1./3.);
}

// Compute the negative of the geometric gradients associated with a given vertex:
/*2017.09.17 NB modified:  Renamed function from "constraint_gradients(~)".  Functionally identical.
 * Just confusing to have the tension forces rely upon a call to "constraint_gradients(~)" and its return
 * defining the negative gradient as "dummy_vforce", when in reality it is the (-gradient*lambda/m) = force.
 * The SHAKE & RATTLE functions reflect this and also have a flag/conversion to quadratic constraint if desired. */
void VERTEX::compute_neg_VGrads(INTERFACE* boundary/*, char constraintFlag*/)
{
  VECTOR3D temp_GradAki;        // Volume gradient of F_k  w.r.t.  v_i .
  VECTOR3D temp_GradVki;        // Area gradient of F_k  w.r.t.  v_i .      NB added area portions.
  neg_GradVi = VECTOR3D(0,0,0); // Net volume gradient w.r.t. v_i .
  neg_GradAi = VECTOR3D(0,0,0);	// Net area gradient w.r.t. v_i .
  
  for (unsigned int i = 0; i < itsF.size(); i++)
  {
    itsF[i]->computegradients(this, &temp_GradAki, &temp_GradVki);
    (boundary->areaGradient[make_pair(this,itsF[i])]) = temp_GradAki;
    neg_GradVi -= (temp_GradVki); // NB:  Scaling by (1/3) not needed as one would expect from A_i = Sum_k[(1/3)*A_k].
    neg_GradAi -= (temp_GradAki); // As Grad(A_i) = {Sum_k[(1/3)*Grad(A_k_i)] / A_i} where A_i above has same prefactor.
  }                               // In English, each face contributes 1/3 the rate of change but also 1/3 the net area.
/*  if(constraintFlag == 'Q')
  {
    neg_GradVi *= (2)*(boundary->total_volume - boundary->ref_volume); // NB added:  for a quadratic constraint (V - V_0)^2 = 0.
    neg_GradAi *= (2)*(boundary->total_area - boundary->ref_area);  // NB added:  for a quadratic constraint (A - A_0)^2 = 0.
  }*/
}

// Compute the force on a given vertex due to the stretching energy penalty:
void VERTEX::stretching_forces(INTERFACE* boundary)
{
  VECTOR3D temp_sforce;       // The stretching force due to one edge on the vertex.
  sforce = VECTOR3D(0,0,0);   // The net stretching force due to all eges on the vertex (iterative sum of the above).
  
  for (unsigned int i = 0; i < itsE.size(); i++)
  {
    temp_sforce = ( (itsE[i]->length() - itsE[i]->l0) / itsE[i]->length() ) ^ (posvec - itsE[i]->opposite(this)->posvec);
    
    temp_sforce = temp_sforce ^ (1.0 / (itsE[i]->l0 * itsE[i]->l0) );
    
    sforce += temp_sforce;
  }
  
  sforce *= (-1) * boundary->sconstant * boundary->avg_edge_length * boundary->avg_edge_length;
}

// (2017.09.17 NB added.)  Compute the force on a given vertex due to the surface tension energy penalty:
void VERTEX::tension_forces(INTERFACE* boundary)
{
    //  The net surface tension force on a vertex is (sigma_a) * neg_GradAi:
    //TForce = ((2.0 * boundary->sigma_a * ((boundary->total_Area_Vertices / boundary->ref_Area_Vertices) - 1)) ^ neg_GradAi);  // Quadratic area difference from sphere.
    //TForce = ((2.0 * boundary->sigma_a * boundary->total_Area_Vertices / boundary->ref_Area_Vertices) ^ neg_GradAi);          // Quadratic in absolute area.
    TForce = (boundary->sigma_a ^ neg_GradAi);                                                                                  // Linear in absolute area.
}

// Compute the force on a given vertex due to the volume tension energy penalty:
void VERTEX::volume_tension_forces(INTERFACE* boundary)
{   //  The net volume tension force on a vertex:
    VolTForce = ((2.0 * boundary->sigma_v * (boundary->total_volume - boundary->ref_volume)) ^ neg_GradVi);  // Quadratic volume difference from sphere.
    //VolTForce = (2.0 * boundary->sigma_v * (boundary->total_volume - boundary->ref_volume) ^ neg_GradVi);          // Quadratic in absolute area.
    //VolTForce = (boundary->sigma_v ^ neg_GradVi); // Linear in absolute volume.
}

// ### Unused functions: ###
/*// gradient of the volume constraint V - V0 = 0 wrt the current vertex
void VERTEX::volume_constraint_gradient(INTERFACE* boundary)
{
  VECTOR3D temp_vforce;
  VECTOR3D temp_aforce; // redundant here
  neg_GradVi = VECTOR3D(0,0,0);
  for (unsigned int i = 0; i < itsF.size(); i++)
  {
//    itsF[i]->computegradients(this, &temp_aforce, &temp_vforce);
    neg_GradVi += temp_vforce;
  }
  neg_GradVi *= (-1);	// we need negative gradient in shake implementation
}
// total area is used in area constraint force
void VERTEX::constraint_forces(INTERFACE* boundary)
{
  VECTOR3D temp_aforce;
  VECTOR3D temp_vforce;
  aforce = VECTOR3D(0,0,0);
  vforce = VECTOR3D(0,0,0);
  lforce = VECTOR3D(0,0,0);
  
  for (unsigned int i = 0; i < itsF.size(); i++)
  {
    itsF[i]->computegradients(this, &temp_aforce, &temp_vforce);
    (boundary->areaGradient[make_pair(this,itsF[i])]) = temp_aforce;
    aforce += temp_aforce;
    vforce += temp_vforce;
  }
  // global volume constraint
//   vforce *= (-2) * boundary->lambda_v * (boundary->total_volume - boundary->ref_volume);
  // local volume constraint
  vforce *= (-1) * boundary->lambda_v;
  
  // dummy_vforce used in shake
//   dummy_vforce = (-1) ^ vforce;
  
  // global area constraint								// added NB
//  aforce *= (-2) * boundary->lambda_a * (boundary->total_area - boundary->ref_area);	// added NB
  // linear area constraint (global = local)
  aforce *= (-1) * boundary->lambda_a;
  
//  for (unsigned int i = 0; i < itsE.size(); i++)
//    if (itsE[i]->isOnBoundary)
//      lforce -= (posvec-itsE[i]->opposite(this)->posvec) ^ 1./itsE[i]->length();
//  lforce *= boundary->lambda_l;
  
  for (unsigned int i = 0; i < itsF.size(); i++)
  {
    itsF[i]->computecenter();
    //lforce -= itsF[i]->boundary_force();
    for (unsigned int j = 0; j < itsF[i]->itsE.size(); j++)
      if (itsF[i]->itsE[j]->isOnBoundary)
      {
	VECTOR3D tmp = (itsF[i]->itscenter
	                - itsF[i]->itsE[j]->opposite(itsF[i])->itscenter);
	tmp.normalize();
	lforce -= tmp;
      }
  }
  lforce *= (1.0/3.0)*boundary->lambda_l;

  forvec = aforce + vforce + lforce;
}*/
