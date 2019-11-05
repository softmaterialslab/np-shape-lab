#ifndef _NEWFORCES_H
#define _NEWFORCES_H

#include "vertex.h"

void force_calculation_init(INTERFACE &, vector<PARTICLE> &, const double, char bucklingFlag);	// calculate the forces
void force_calculation(INTERFACE &, vector<PARTICLE> &, const double, char bucklingFlag);	// calculate the forces

#endif 
