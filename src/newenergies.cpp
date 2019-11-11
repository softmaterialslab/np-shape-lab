#include "newenergies.h"

// Common LJ function for all pairs:
double energy_lj_pair(VECTOR3D &pos1, VECTOR3D &pos2, double d, double elj)
{
  VECTOR3D r_vec = pos1 - pos2;
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

// ### ES Pair-specific functions (overloaded for mesh and particle classes, but entities always referred to as 'vN'):
    // For vertex-vertex pairs:
double energy_es_pair(VERTEX &v1, VERTEX &v2, double em, double inv_kappa, const double scalefactor)
{
  double r = (v1.posvec - v2.posvec).GetMagnitude();
  return(scalefactor * v1.q * v2.q * ( exp(-(1.0/inv_kappa) * r) ) / (em * r));
}
    //  For vertex-particle pairs (both permutations, for safety):
double energy_es_pair(VERTEX &v1, PARTICLE &v2, double em, double inv_kappa, const double scalefactor)
{
    double r = (v1.posvec - v2.posvec).GetMagnitude();
    return(scalefactor * v1.q * v2.q * ( exp(-(1.0/inv_kappa) * r) ) / (em * r));
}
double energy_es_pair(PARTICLE &v1, VERTEX &v2, double em, double inv_kappa, const double scalefactor)
{
    double r = (v1.posvec - v2.posvec).GetMagnitude();
    return(scalefactor * v1.q * v2.q * ( exp(-(1.0/inv_kappa) * r) ) / (em * r));
}
    // For particle-particle pairs:
double energy_es_pair(PARTICLE &v1, PARTICLE &v2, double em, double inv_kappa, const double scalefactor)
{
    double r = (v1.posvec - v2.posvec).GetMagnitude();
    return(scalefactor * v1.q * v2.q * ( exp(-(1.0/inv_kappa) * r) ) / (em * r));
}