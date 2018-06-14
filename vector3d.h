// This is 3d vector class
// VECTOR3D is a vector in cartesian coordinate system

#ifndef _VECTOR3D_H
#define _VECTOR3D_H

#include "utility.h"

class VECTOR3D 
{
  public:
    long double x, y, z;									// component along each axis (cartesian)
    VECTOR3D(long double xx = 0.0, long double yy = 0.0, long double zz = 0.0) : x(xx), y(yy), z(zz) 	// make a 3d vector
    {
    } 
    long double GetMagnitude()			// magnitude of the vector
    {
      return sqrt(x*x + y*y +z*z);							
    }

    void normalize()
    {
      long double mag = GetMagnitude();
      x /= mag;
      y /= mag;
      z /= mag;
    }

    long double GetMagnitudeSquared()			// magnitude of the vector
    {
      return (x*x + y*y + z*z);							
    }

    VECTOR3D operator-(const VECTOR3D& vec)	// subtract two vectors
    {
      return VECTOR3D(x - vec.x, y - vec.y, z - vec.z);					
    }

    long double operator*(const VECTOR3D& vec)	// dot product of two vectors
    {
      return x*vec.x + y*vec.y + z*vec.z;						
    }

//     VECTOR3D operator^(long double scalar)	// product of a vector and a scalar
//     {
//       return VECTOR3D(x*scalar, y*scalar, z*scalar);					
//     }
    VECTOR3D operator+(const VECTOR3D& vec)	// add two vectors
    {
      return VECTOR3D(x + vec.x, y + vec.y, z + vec.z);					
    }

    bool operator==(const VECTOR3D& vec)		// compare two vectors
    {
      if (x==vec.x && y==vec.y && z==vec.z) 
	return true;
      return false;
    }

    VECTOR3D operator%(const VECTOR3D& vec)
    {
      return(VECTOR3D(y*vec.z - z*vec.y, z*vec.x - x*vec.z, x*vec.y - y*vec.x));
    }

    // overloading *=
    void operator*=(const long double factor)	// const used
    {
      x *= factor;
      y *= factor;
      z *= factor;
    }

    // reluctantly overloading +=
    void operator+=(const VECTOR3D& vec)
    {
      x += vec.x;
      y += vec.y;
      z += vec.z;
    }

    // reluctantly overloading -=
    void operator-=(const VECTOR3D& vec)
    {
      x -= vec.x;
      y -= vec.y;
      z -= vec.z;
    }
};

VECTOR3D operator^(const VECTOR3D& vec, long double scalar);	// product of a vector and a scalar
VECTOR3D operator^(long double scalar, const VECTOR3D& vec);	// product of a vector and a scalar

// VECTOR3D operator^(const VECTOR3D& vec, long double scalar)	// product of a vector and a scalar
// {
//   return VECTOR3D(vec.x*scalar, vec.y*scalar, vec.z*scalar);					
// }
// 
// VECTOR3D operator^(long double scalar, const VECTOR3D& vec)	// product of a vector and a scalar
// {
//       return VECTOR3D(vec.x*scalar, vec.y*scalar, vec.z*scalar);					
// //   return vec^scalar;					
// }
#endif

