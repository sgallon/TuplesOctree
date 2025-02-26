#ifndef POINT_HPP
#define POINT_HPP

#include <iostream>

#include <math.h>
#include <vector>


namespace octree{
class Point
{
  public:
    Point();
    Point(double X_, double Y_, double Z_);
    double X;
    double Y;
    double Z;


    Point operator+( const Point& p2) const;
    Point& operator+=(const Point& p2);
    Point operator-( const Point& p2) const;
    Point& operator-=(const Point& p2);
    Point operator*( double alpha) const;
    Point& operator*=(double alpha);
    Point operator/( double alpha) const;

    bool operator==(const Point& p2) const;
    bool operator!=(const Point& p2) const;
    
    double operator[](const int i) const;
    double& operator[](const int i);
    friend std::ostream& operator<<(std::ostream& os,const Point& p1 );


};
using PointVec = std::vector< Point > ;

Point elementwise_absolute(const Point& p1);

double norm(const Point& p1);
double norm(const Point* p1);
double norm_difference(const Point& p1 ,const Point& p2);
double norm_difference(const Point* p1 ,const Point* p2);
double squared_norm(const Point& p1);
double squared_norm(const Point* p1);
double squared_norm_difference(const Point& p1 , const Point& p2);
double squared_norm_difference(const Point& p1 , const Point& p2);

double cosAngleToZAxis(const Point& p1, const Point& p2);
double cosAngleToZAxis(const Point* p1, const Point* p2);

double dotproduct(Point p1, Point p2);
Point crossproduct(Point p1, Point p2);


} // namespace
#endif //POINT_HPP
