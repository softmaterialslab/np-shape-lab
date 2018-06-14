#include "vector3d.h"

VECTOR3D operator^(const VECTOR3D& vec, long double scalar)	// product of a vector and a scalar
{
  return VECTOR3D(vec.x*scalar, vec.y*scalar, vec.z*scalar);					
}

VECTOR3D operator^(long double scalar, const VECTOR3D& vec)	// product of a vector and a scalar
{
      return VECTOR3D(vec.x*scalar, vec.y*scalar, vec.z*scalar);					
//   return vec^scalar;					
}
