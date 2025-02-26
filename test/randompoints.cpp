#ifndef RANDOMPOINTS
#define RANDOMPOINTS
#include <random>
#include <point.hpp>

std::random_device rd;
// std::mt19937 gen( rd() );
std::mt19937 gen(42);
octree::PointVec generateRandomPointArr( size_t N)
{
  std::uniform_real_distribution<> dis(0., 2.);

  octree::PointVec temp(N);
  for( size_t idx = 0 ;  idx<N ; idx++)
  {
    temp[idx] = { dis(gen), dis(gen),dis(gen)};
  }
  return temp;
}

#endif