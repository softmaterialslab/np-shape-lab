#include "newforces.h"
#include "functions.h"

VECTOR3D lj_vertex_vertex(VERTEX& v1, VERTEX& v2, double d, double elj)
{
  VECTOR3D r_vec = v1.posvec - v2.posvec;
  double r2 = r_vec.GetMagnitudeSquared();
  double d2 = d * d;
  
  if (r2 < dcut2 * d2)
  {
    double r6 = r2 * r2 * r2;
    double r12 = r6 * r6;
    double d6 = d2 * d2 * d2;
    double d12 = d6 * d6;
    return(r_vec ^ ( 48 * elj * ( (d12 / r12) - 0.5 * (d6 / r6) ) * (1 / r2) ));
  } 
  else
    return VECTOR3D(0,0,0);
}

VECTOR3D es_vertex_vertex(VERTEX& v1, VERTEX& v2, double em, double inv_kappa)
{
  // salt-free interaction
//   return((Grad(v1.posvec, v2.posvec)) ^ (-scalefactor * v1.q * v2.q / em));
  
  // presence of salt
  return( (-scalefactor * v1.q * v2.q / em) * (exp(-(1.0/inv_kappa) * (v1.posvec - v2.posvec).GetMagnitude())) ^ 
  (Grad(v1.posvec, v2.posvec) + ( ( (-1.0/inv_kappa) * (1 / (v1.posvec - v2.posvec).GetMagnitudeSquared()) ) ^ (v1.posvec - v2.posvec) ) ) );
}
