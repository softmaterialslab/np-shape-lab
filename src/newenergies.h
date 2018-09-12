#ifndef _NEWENERGIES_H
#define _NEWENERGIES_H

#include "vertex.h"
#include "thermostat.h"

// Computes the kinetic energy of the vertices && returns the net kinetic energy for assignment, if desired:
inline long double vertex_kinetic_energy(vector<VERTEX>& V)
{
  for (unsigned int i = 0; i < V.size(); i++)
    V[i].real_kinetic_energy();
  long double kinetic_energy = 0.0;
  for (unsigned int i = 0; i < V.size(); i++)
    kinetic_energy += V[i].realke;
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

double energy_lj_vertex_vertex(VERTEX&, VERTEX&, double, double);
double energy_es_vertex_vertex(VERTEX&, VERTEX&, double, double);

#endif
