// This is header file for the UTILITY class.
// This file includes standard library files and gsl functions that are utilized in the code
// This file also has useful constant parameters for the problem at hand

#ifndef _MPI_UTILITY_H
#define _MPI_UTILITY_H

extern mpi::environment env;
extern mpi::communicator world;

//MPI boundary parameters
extern unsigned int lowerBoundMesh;
extern unsigned int upperBoundMesh;
extern int lowerBoundIons;
extern int upperBoundIons;
extern unsigned int sizFVecMesh;
extern int sizFVecIons;
//extern unsigned int extraElementsMesh;
//extern unsigned int extraElementsIons;

#endif
