#include "newenergies.h"

double energy_lj_vertex_vertex(VERTEX& v1, VERTEX& v2, double d, double elj)
{
  VECTOR3D r_vec = v1.posvec - v2.posvec;
  double r2 = r_vec.GetMagnitudeSquared();
  double d2 = d * d;
  
  if (r2 < dcut2 * d2)
  {
    double r6 = r2 * r2 * r2;
    double d6 = d2 * d2 * d2;
    return(4 * elj * (d6 / r6) * ( ( d6 / r6 ) - 1 ) + elj);
  } 
  else
    return 0;
}

double energy_es_vertex_vertex(VERTEX& v1, VERTEX& v2, double em, double inv_kappa, const double scalefactor)
{
  double r = (v1.posvec - v2.posvec).GetMagnitude();
  return(scalefactor * v1.q * v2.q * ( exp(-(1.0/inv_kappa) * r) ) / (em * r));
}
