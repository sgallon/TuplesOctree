#include "point.hpp"

namespace octree{
Point::Point():
  X(0.),Y(0.),Z(0.)
{}

Point::Point(double X_, double Y_, double Z_):
  X(X_),Y(Y_),Z(Z_)
{}

Point Point::operator+(const Point& p2) const {
  return Point(this-> X + p2.X, this->Y + p2.Y , this->Z + p2.Z  );
}
Point& Point::operator+=(const Point& p2){
  this-> X += p2.X;
  this-> Y += p2.Y;
  this-> Z += p2.Z;
  return *this;
}


Point Point::operator-(const Point& p2) const {
  return Point(this-> X - p2.X, this->Y - p2.Y , this->Z - p2.Z  );
}
Point& Point::operator-=(const Point& p2){
  this-> X -= p2.X;
  this-> Y -= p2.Y;
  this-> Z -= p2.Z;
  return *this;
}


Point Point::operator*( double alpha) const {
  return Point(this-> X * alpha, this->Y * alpha , this->Z * alpha );
}
Point& Point::operator*=(double alpha) {
  this->X *= alpha;
  this->Y *= alpha;
  this->Z *= alpha;
  return *this;
}


Point Point::operator/( const double alpha) const {
  return Point(this-> X / alpha, this->Y / alpha, this->Z/alpha);
}

bool Point::operator==(const Point& p2) const {
  return (this->X == p2.X && this->Y == p2.Y && this->Z == p2.Z);
}

bool Point::operator!=(const Point& p2) const {
  return !(*this == p2);
}

double& Point::operator[](const int i) {
  if(i==0){
    return this->X;
  }
  else if (i==1){
    return this->Y;
  }
  else if (i==2){
    return this->Z;
  }
  else{
    std::cerr << " Invalid argument: " <<  i << "for octree::Point[]" << std::endl;
    throw std::runtime_error("InvalidArgument");
  }
}

double Point::operator[](const int i) const{
  if(i==0){
    return this->X;
  }
  else if (i==1){
    return this->Y;
  }
  else if (i==2){
    return this->Z;
  }
  else{
    std::cerr << " Invalid argument: " <<  i << "for octree::Point[]" << std::endl;
    return 0.;
  }
}



std::ostream& operator<<(std::ostream& os,const Point& p1){
  return os << "{ "<< p1.X << " " << p1.Y << " " << p1.Z << " }" ;
}


Point elementwise_absolute(const Point& p1){
  // Elementwise abs
  return Point(std::abs(p1.X),std::abs(p1.Y),std::abs(p1.Z));
}
double norm_difference(const Point& p1 , const Point& p2){
  double dx = (p2.X - p1.X);
  double dy = (p2.Y - p1.Y);
  double dz = (p2.Z - p1.Z);
  return sqrt(dx*dx + dy*dy + dz*dz);
}
double norm_difference(const Point* p1 , const Point* p2){
  double dx = (p2->X - p1->X);
  double dy = (p2->Y - p1->Y);
  double dz = (p2->Z - p1->Z);
  return sqrt(dx*dx + dy*dy + dz*dz);
}
double norm(const Point& p1){
  return sqrt(p1.X*p1.X + p1.Y *p1.Y + p1.Z *p1.Z);
}
double norm(Point* p1){
  return sqrt(p1->X*p1->X + p1->Y *p1->Y + p1->Z *p1->Z);
}
double norm(const double x, const double y, const double z){
  return sqrt(x*x + y*y + z*z);
}

double squared_norm(const Point& p1){
  return p1.X * p1.X + p1.Y * p1.Y + p1.Z * p1.Z;
}
double squared_norm(const Point* p1){
  return p1->X * p1->X + p1->Y * p1->Y + p1->Z * p1->Z;
}

double squared_norm_difference(const Point& p1 , const Point& p2){
  double dx = (p2.X - p1.X);
  double dy = (p2.Y - p1.Y);
  double dz = (p2.Z - p1.Z);
  return (dx*dx + dy*dy + dz*dz);
}

double squared_norm_difference(const Point* p1 , const Point* p2){
  double dx = (p2->X - p1->X);
  double dy = (p2->Y - p1->Y);
  double dz = (p2->Z - p1->Z);
  return (dx*dx + dy*dy + dz*dz);
}

double cosAngleToZAxis(const Point& p1, const Point& p2){
  double dx = p2.X - p1.X;
  double dy = p2.Y - p1.Y;
  double dz = std::abs(p2.Z - p1.Z);
  double norm = sqrt(dx*dx + dy*dy + dz*dz);
  return dz / norm;
}

double cosAngleToZAxis(const Point* p1, const Point* p2){
  double dx = p2->X - p1->X;
  double dy = p2->Y - p1->Y;
  double dz = std::abs(p2->Z - p1->Z);
  double norm = sqrt(dx*dx + dy*dy + dz*dz);
  return dz / norm;
}



double dotproduct(Point p1, Point p2){
  return p1.X * p2.X + p1.Y *p2.Y + p1.Z * p2.Z;
}
Point crossproduct(Point p1, Point p2){
  return Point(p1.Y * p2.Z - p1.Z * p2.Y,
               p1.Z * p2.X - p1.X * p2.Z,
               p1.X * p2.Y - p1.Y * p2.X);
}
} // namespace
