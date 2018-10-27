//****************************************************
//* quaternion.h                                     *
//*                                                  *
//* Implementaion for a generalized quaternion class *   
//*                                                  *
//* Written 1.25.00 by Angela Bennett                *
//****************************************************


#ifndef _QUATERNION_H_
#define _QUATERNION_H_

#include <iostream>
#include <math.h>

#ifdef SHOEMAKE
#include "EulerAngles.h"
#endif

template<class DataType> 
class Quaternion
{

 public:
  
  //Quaternion
  // -default constructor
  // -creates a new quaternion with all parts equal to zero
  Quaternion(void);
  
  //Quaternion
  // -constructor
  // -parametes : w, x, y, z elements of the quaternion
  // -creates a new quaternion based on the elements passed in
  Quaternion(DataType wi, DataType xi, DataType yi, DataType zi);
  
  //Quaternion
  // -constructor
  // -parameters : 4D vector
  // -creates a new quaternion based on the elements passed in
  Quaternion(DataType v[4]);

  std::array<DataType,4> get() const { 
       std::array<DataType,4> ref;
       ref[0] = w;
       ref[1] = x;
       ref[2] = y;
       ref[3] = z;
       return ref;
  }

  DataType get(int i) const { 
       DataType val = 0.0;
       switch(i)
       {
         case 0:
              val = w;
              break;
         case 1:
              val = x;
              break;
         case 2:
              val = y;
              break;
         case 3:
              val = z;
              break;
       }
       return val;
   }

  //Quaternion
  // -copy constructor
  // -parameters : const quaternion q
  // -creates a new quaternion based on the quaternion passed in
  Quaternion(const Quaternion<DataType>& q); 

#ifdef SHOEMAKE
  //Quaternion
  // -constructor
  // -parameters : yaw, pitch, and roll of an Euler angle
  // -creates a new quaternion based on the Euler elements passed in
  // -used with Shoemakes code
  Quaternion(DataType e[3], int order);
#endif  

  //~Quaternion
  // -default destructor
  ~Quaternion();
  
  //operator=
  // -parameters : q1- Quaternion object
  // -return values : Quaternion
  // -when called on quaternion q2 sets q2 to be an object of  q3 
  Quaternion<DataType> operator = (const Quaternion<DataType>& q);
 
  //operator+
  // -parameters : q1 - Quaternion object
  // -return value : Quaternion 
  // -when called on quaternion q2 adds q1 + q2 and returns the sum in a new quaternion
  Quaternion<DataType> operator + (const Quaternion<DataType>& q);
  
  //operator-
  // -parameters : q1- Quaternion object
  // -return values : Quaternion 
  // -when called on q1 subtracts q1 - q2 and returns the difference as a new quaternion
  Quaternion<DataType> operator - (const Quaternion<DataType>& q);

  //operator*
  // -parameters : q1 - Quaternion object
  // -return values : Quaternion 
  // -when called on a quaternion q2, multiplies q2 *q1  and returns the product in a new quaternion 
  Quaternion<DataType> operator * (const Quaternion<DataType>& q);
  
  //operator/
  // -parameters : q1 and q2- Quaternion objects
  // -return values : Quaternion 
  // -divide q1 by q2 and returns the quotient as q1 
  Quaternion<DataType> operator / (Quaternion<DataType>& q);
  
  //operator+=
  // -parameters : q1- Quaternion object
  // -return values : Quaternion 
  // -when called on quaternion q3 adds q1 and q3 and returns the sum as q3 
  Quaternion<DataType>& operator += (const Quaternion<DataType>& q);
  
  //operator-=
  // -parameters : q1- Quaternion object
  // -return values : Quaternion 
  // -when called on quaternion q3, subtracts q1 from q3 and returns the difference as q3
  Quaternion<DataType>& operator -= (const Quaternion<DataType>& q);
 
  //operator*=
  // -parameters : q1- Quaternion object
  // -return values : Quaternion 
  // -when called on quaternion q3, multiplies q3 by q1 and returns the product as q3
  Quaternion<DataType>& operator *= (const Quaternion<DataType>& q);
 
  //operator/=
  // -parameters : q1- Quaternion object
  // -return values : quaternion
  // -when called on quaternion q3, divides q3 by q1 and returns the quotient as q3
  Quaternion<DataType>& operator /= (Quaternion<DataType>& q);
  
  //operator<<
  // -parameters : ostream o, quaternion q
  // -return values :
  // -prints out a quaternion by it's components
  friend inline std::ostream& operator << (std::ostream& output, const Quaternion<DataType>& q)
    {
      output << "[" << q.w << ", " << "(" << q.x << ", " << q.y << ", " << q.z << ")]";
      return output; 
    }
  
  //operator!=
  // -parameters : q1 and q2- Quaternion objects
  // -return value : bool
  // -determines if q1 and q2 and equal
  bool operator != (const Quaternion<DataType>& q);
  
  //operator==
  // -parameters : q1 and q2- Quaternion objects
  // -return value : bool
  // -determines if q1 and q2 and equal
  bool operator == (const Quaternion<DataType>& q);  
    
  //other methods: norm, inverse, conjugate, toEuler
  
  //norm
  // -parameters : none
  // -return value : DataType
  // -when called on a quaternion object q, returns the norm of q
  DataType norm();
  
  //magnitude
  // -parameters : none
  // -return value : DataType
  // -when called on a quaternion object q, returns the magnitude q
  DataType magnitude();
  
  //scale
  // -parameters :  s- a value to scale q1 by
  // -return value: quaternion
  // -returns the original quaternion with each part, w,x,y,z, multiplied by some scalar s
  Quaternion<DataType> scale(DataType s);
  
  //inverse
  // -parameters : none
  // -return value : quaternion
  // -when called on a quaternion object q, returns the inverse of q
  Quaternion<DataType> inverse();
  
  //conjugate
  // -parameters : none
  // -return value : quaternion
  // -when called on a quaternion object q, returns the conjugate of q
  Quaternion<DataType> conjugate();
  
  //UnitQuaternion
  // -parameters : none
  // -return value : quaternion
  // -when called on quaterion q, takes q and returns the unit quaternion of q
  Quaternion<DataType> UnitQuaternion();
  
  // -parameters : 3D vector of type DataType
  // -return value : void
  // -when given a  3D vector, v, rotates v by the quaternion
  void QuatRotation(DataType v[3]);
  
#ifdef SHOEMAKE
  // -parameters : empty 3D vector, rotation order
  // -return : void
  // - converts this quaternion into Euler angles
  void toEuler(DataType e[3], int order);
#endif

 private:
  // [w, (x, y, z)]
  DataType w, x, y, z;
  
};

#include "Quaternion.cpp"

#endif

