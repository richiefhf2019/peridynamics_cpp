#ifndef VEC_H
#define VEC_H

#include "realtypes.h"

namespace periDynamics {
  
class Vec{
 public:
  Vec():x(0),y(0),z(0) {}
  Vec(REAL d):x(d),y(d),z(d) {}
  Vec(REAL _x, REAL _y, REAL _z):x(_x),y(_y),z(_z) {}
  
  REAL getx() const {return x;}
  REAL gety() const {return y;}
  REAL getz() const {return z;}
  void setx(REAL _x) {x = _x;}
  void sety(REAL _y) {y = _y;}
  void setz(REAL _z) {z = _z;}
  void set(REAL _x, REAL _y, REAL _z) {x = _x; y = _y; z = _z;}
  void set(Vec v) {x = v.getx(); y = v.gety(); z = v.getz();}
  
  bool operator==(const Vec v);
  bool operator==(const REAL d);   
  bool operator!=(const Vec v); 
  void operator+=(const Vec v);
  void operator-=(const Vec v);
  void operator*=(REAL d);
  void operator/=(REAL d);
  Vec  operator+(Vec v) const;
  Vec  operator-(Vec v) const;
  Vec  operator*(Vec p) const;   // cross product of this vector and p
  Vec  operator*(REAL d) const;
  REAL operator%(Vec p) const;   // dot product of this and p
  void print() const;
  
 private:
  REAL x;
  REAL y;
  REAL z;
};
 
// Non-member functions
Vec operator*(REAL d, Vec v);
Vec operator/(Vec v, REAL d);
Vec operator-(Vec v);
REAL vfabs(Vec v);
Vec vcos(Vec v);
Vec vacos(Vec v);
Vec rotateVec(Vec v, Vec alf);    // find the exact vector after v is rotated alf in space
Vec normalize(Vec v);
/*calculate the angle between v1 and v2 if rotating v1 in the plane
  composed of v1 and v2 from itself to v2, the angle could be 0<alf<360
  norm specify that the rotation must be around norm according to right hand rule,
  even if 180<alf<360
*/
REAL angle(Vec v1, Vec v2, Vec norm); 
 
} // namespace periDynamics

#endif
