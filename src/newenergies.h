#ifndef _NEWENERGIES_H
#define _NEWENERGIES_H

#include "vertex.h"
#include "thermostat.h"

// Computes the kinetic energy of the vertices & returns the net kinetic energy for assignment, if desired:
inline long double compute_kinetic_energy(vector<VERTEX>& V)
{
  for (unsigned int i = 0; i < V.size(); i++)
    V[i].real_kinetic_energy();
  long double kinetic_energy = 0.0;
  for (unsigned int i = 0; i < V.size(); i++)
    kinetic_energy += V[i].ke;
  return kinetic_energy;
}
// Overloads to also compute the KE of the ions & return their net kinetic energy for assignment, if desired:
inline long double compute_kinetic_energy(vector<PARTICLE>& counterions)
{
    for (unsigned int i = 0; i < counterions.size(); i++)
        counterions[i].real_kinetic_energy();
    long double kinetic_energy = 0.0;
    for (unsigned int i = 0; i < counterions.size(); i++)
        kinetic_energy += counterions[i].ke;
    return kinetic_energy;
}

inline double bath_kinetic_energy(vector<THERMOSTAT>& bath)
{
  for (unsigned int j = 0; j < bath.size(); j++)
    bath[j].kinetic_energy();
  double kinetic_energy = 0.0;
  for (unsigned int j = 0; j < bath.size(); j++)
    kinetic_energy += bath[j].ke;
  return kinetic_energy;
}

inline double bath_potential_energy(vector<THERMOSTAT>& bath)
{
  for (unsigned int j = 0; j < bath.size(); j++)
    bath[j].potential_energy();
  double potential_energy = 0.0;
  for (unsigned int j = 0; j < bath.size(); j++)
    potential_energy += bath[j].pe;
  return potential_energy;
}

double energy_lj_pair(VECTOR3D &pos1, VECTOR3D &pos2, double d, double elj);
double energy_es_pair(VERTEX &v1, VERTEX &v2, double em, double inv_kappa, const double scalefactor);
double energy_es_pair(VERTEX &v1, PARTICLE &v2, double em, double inv_kappa, const double scalefactor);
double energy_es_pair(PARTICLE &v1, VERTEX &v2, double em, double inv_kappa, const double scalefactor);
double energy_es_pair(PARTICLE &v1, PARTICLE &v2, double em, double inv_kappa, const double scalefactor);

#endif
