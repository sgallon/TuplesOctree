#include <iostream>
#include "randompoints.cpp"
#include "testtuples.cpp"

#include "octree.hpp"
#include "tuple.hpp"

#include <algorithm>
#include <chrono>




int triangle_test_random( size_t N,size_t maxTriangle , int depth, double distance, double tolerance, bool doDirectSearch = true)
{
  auto randomPoints = generateRandomPointArr(N);

  octree::IndexTripletVec tuples_direct;
  if(doDirectSearch){
    tuples_direct = direct_search_triangles(randomPoints,distance,tolerance);
    if(tuples_direct.size()==0){
      std::cout << "No triangles found in direct search. Update parameter" << std::endl;
      return 1;
    }
  }
  auto tree = octree::TOctree(randomPoints,{1.,1.,1.},2.,depth);
  auto tuples = tree.search_triangles(distance,tolerance,maxTriangle);
  
  if((doDirectSearch) && (tuples.size() < tuples_direct.size())){
    std::cout << "Number of triangles in direct search: " << tuples_direct.size() << " vs " << tuples.size() << std::endl;
  }
  if (tuples.size() == 0) {
    std::cout << "No triangles found" << std::endl;
    return 1;
  }
  std::cout << "Searched for " << maxTriangle << " triangles, found " << tuples.size() << std::endl;
  auto mismatch = test_triangles(randomPoints,tuples,distance,tolerance);
  if(mismatch != 0){
    std::cout << "number of mismatches: " << mismatch << std::endl;
    return mismatch;
  }else{
    std::cout << "no mismatches found" << std::endl;
  }
  auto nDuplicates = count_duplicates(tuples);
  if(nDuplicates != 0){
    std::cout << "number of duplicates: " << nDuplicates << std::endl;
    return nDuplicates;
  }else{
    std::cout << "no duplicates found" << std::endl;
  }
  if((doDirectSearch) && (tuples_direct.size() > maxTriangle) && (tuples.size() < 0.9 *maxTriangle)){
    std::cout << "Tree search found less than 90% of triangles than maxTriangle, but direct search found more than maxTriangle" << std::endl;
    return 1;
  }
  auto maxTuplePerP1 = count_max_tuplePerParticle1(tuples);
  if((doDirectSearch) && 
     (tuples_direct.size() > 5 * maxTriangle) &&
     (maxTuplePerP1 > std::max(1.,(double) tuples_direct.size()/N)) ){
    std::cout << "Unsufficent randomization of tuples" << std::endl;
    return 1;
  }
  return 0;
}


int triangle_test( size_t N, int depth, double distance, double tolerance)
{
  auto randomPoints = generateRandomPointArr(N);
  auto tree = octree::TOctree(randomPoints,{1.,1.,1.},2.,depth);

  auto tuples = tree.search_triangles(distance,tolerance);
  std::cout << "Number of triangles: " << tuples.size() << std::endl;

  // test for wrong pairs
  auto mismatch = test_triangles(randomPoints,tuples,distance,tolerance);
  if(mismatch != 0){
    std::cout << "number of mismatches: " << mismatch << std::endl;
    return mismatch;
  }else{
    std::cout << "no mismatches found" << std::endl;
  }
  // test for duplicates
  auto nDuplicates = count_duplicates(tuples);
  if(nDuplicates != 0){
    std::cout << "number of duplicates: " << nDuplicates << std::endl;
    return nDuplicates;
  }else{
    std::cout << "no duplicates found" << std::endl;
  }
  auto triangles_direct = direct_search_triangles(randomPoints,distance,tolerance);
  if(tuples.size() != triangles_direct.size()){
    std::cout << "Number of triangles does not match: " << tuples.size() << " vs " << triangles_direct.size() << std::endl;
    sort_triplets(triangles_direct);
    size_t iMax = std::min(tuples.size(),triangles_direct.size());
    for(size_t i = 0; i<iMax; ++i){
      if(tuples[i] != triangles_direct[i]){
        std::cout << "Triangle " << i << " does not match: " << tuples[i][0] << " " << tuples[i][1] << " " << tuples[i][2] << " vs " << triangles_direct[i][0] << " " << triangles_direct[i][1] << " " << triangles_direct[i][2] << std::endl;
      }
    }
    return 1;
  }
  return 0;
}


double test_count( size_t N, int depth, double distance, double tolerance)
{
  auto randomPoints = generateRandomPointArr(N);
  auto tree = octree::TOctree(randomPoints,{1.,1.,1.},2.,depth);

  auto tuples = tree.search_triangles(distance,tolerance);
  auto nTuples = tree.count_triangles(distance,tolerance);
  if(nTuples != tuples.size()){
    std::cout << "Number of triangles does not match: " << nTuples << " vs " << tuples.size() << std::endl;
    return 1;
  }
  return 0;
}


int main( int argc, char *argv[] )
{
  int returnVal = 0;
  returnVal += triangle_test(50000,3,0.2,1e-2);
  if(returnVal > 0){
    std::cout << "test 1 failed" << std::endl;
    return returnVal ;
  }
  returnVal += triangle_test(50000,5,0.2,1e-2);
  if(returnVal > 0){
    std::cout << "test 2 failed" << std::endl;
    return returnVal ;
  }

  
  // test maxTriangle << nPairs
  returnVal += triangle_test_random(50000,1000,3,0.2,1e-2,false);
  if(returnVal > 0){
    std::cout << "test 3 failed" << std::endl;
    return returnVal ;
  }
  returnVal += triangle_test_random(50000,1000,5,0.2,1e-2,false);
  if(returnVal > 0){
    std::cout << "test 4 failed" << std::endl;
    return returnVal ;
  }

  returnVal += triangle_test_random(50000,8000,3,0.2,1e-2,false);
  if(returnVal > 0){
    std::cout << "test 5 failed" << std::endl;
    return returnVal ;
  }
  returnVal += triangle_test_random(50000,8000,5,0.2,1e-2,false);
  if(returnVal > 0){
    std::cout << "test 6 failed" << std::endl;
    return returnVal ;
  }

  returnVal += triangle_test_random(50000,10000,3,0.2,1e-2,true);
  if(returnVal > 0){
    std::cout << "test 7 failed" << std::endl;
    return returnVal ;
  }
  returnVal += triangle_test_random(50000,10000,6,0.2,1e-2,true);
  if(returnVal > 0){
    std::cout << "test 8 failed" << std::endl;
    return returnVal ;
  }
  
  return returnVal;
}