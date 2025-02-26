#include <iostream>
#include "randompoints.cpp"

#include <octree.hpp>
#include <tuple.hpp>

#include <algorithm>
#include <chrono>
#include "testtuples.cpp"

int tetrahedron_test_random(size_t N, size_t maxTetrahedron, int depth, double distance, double tolerance, bool doDirectSearch = true) {
  auto randomPoints = generateRandomPointArr(N);

  octree::IndexQuadrupletVec tetrahedra_direct;
  if(doDirectSearch){
      tetrahedra_direct = direct_search_tetrahedra(randomPoints, distance, tolerance);
      if(tetrahedra_direct.size() == 0){
          std::cout << "No tetrahedra found in direct search. Update parameter" << std::endl;
          return 1;
      }
  }

  auto tree = octree::TOctree(randomPoints, {1., 1., 1.}, 2., depth);
  auto tetrahedra = tree.search_tetrahedra(distance, tolerance, maxTetrahedron);

  if((doDirectSearch) && (tetrahedra.size() < tetrahedra_direct.size())){
      std::cout << "Number of tetrahedra in direct search: " << tetrahedra_direct.size() << " vs " << tetrahedra.size() << std::endl;
  }
  if (tetrahedra.size() == 0) {
      std::cout << "No tetrahedra found" << std::endl;
      return 1;
  }
  std::cout << "Searched for " << maxTetrahedron << " tetrahedra, found " << tetrahedra.size() << std::endl;
  auto mismatch = test_tetrahedra(randomPoints, tetrahedra, distance, tolerance);
  if(mismatch != 0){
      std::cout << "number of mismatches: " << mismatch << std::endl;
      return mismatch;
  }else{
      std::cout << "no mismatches found" << std::endl;
  }

  auto nDuplicates = count_duplicates(tetrahedra);
  if(nDuplicates != 0){
      std::cout << "number of duplicates: " << nDuplicates << std::endl;
      return nDuplicates;
  }else{
      std::cout << "no duplicates found" << std::endl;
  }
  if((doDirectSearch) && (tetrahedra_direct.size() > maxTetrahedron) && (tetrahedra.size() < 0.8 *maxTetrahedron)){
      std::cout << "Tree search found less than 80% of tetrahedra than maxTetrahedron, but direct search found more than maxTetrahedron" << std::endl;
      return 1;
  }
  auto maxTuplePerP1 = count_max_tuplePerParticle1(tetrahedra);
  if((doDirectSearch) &&
      (tetrahedra_direct.size() > 5 * maxTetrahedron) &&
      (maxTuplePerP1 > std::max(1.,(double) 3 * tetrahedra_direct.size()/N)) ){
      std::cout << "Unsufficent randomization of tuples" << std::endl;
      std::cout << "maxTuplePerP1: " << maxTuplePerP1 << std::endl;
      return 1;
  }
  return 0;
}

int tetrahedra_test( size_t N, int depth, double distance, double tolerance)
{
  auto randomPoints = generateRandomPointArr(N);
  auto tree = octree::TOctree(randomPoints,{1.,1.,1.},2.,depth);

  auto tuples = tree.search_tetrahedra(distance,tolerance);
  std::cout << "Number of tetrahedra: " << tuples.size() << std::endl;
  auto mismatch = test_tetrahedra(randomPoints,tuples,distance,tolerance);
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
  auto tetrahedra_direct = direct_search_tetrahedra(randomPoints,distance,tolerance);
  if (tetrahedra_direct.size() == 0) {
    std::cout << "No tetrahedra found in direct search. Update parameter" << std::endl;
    return 1;
  }
  if(tuples.size() != tetrahedra_direct.size()){
    std::cout << "Number of tetrahedra does not match: " << tuples.size() << " vs " << tetrahedra_direct.size() << std::endl;
    sort_quadruplets(tuples);
    size_t iMax = std::min(tuples.size(),tetrahedra_direct.size());
    for(size_t i = 0; i<iMax; ++i){
      if(tuples[i] != tetrahedra_direct[i]){
        std::cout << "Tetrahedron " << i << " does not match: " << tuples[i][0] << " " << tuples[i][1] << " " << tuples[i][2] << " " << tuples[i][3] << " vs " << tetrahedra_direct[i][0] << " " << tetrahedra_direct[i][1] << " " << tetrahedra_direct[i][2] << " " << tetrahedra_direct[i][3] << std::endl;
      }
    }
    return 1;
  } 
  return 0;
}

double test_count(size_t N, int depth, double distance, double tolerance){
  // generate random points
  auto randomPoints = generateRandomPointArr(N);
  // create octree
  auto tree = octree::TOctree(randomPoints,{1.,1.,1.},2.,depth);
  // count tetrahedra
  auto tuples = tree.search_tetrahedra(distance,tolerance);
  auto nTuples = tree.count_tetrahedra(distance,tolerance);
  if(nTuples != tuples.size()){
    std::cout << "Number of tetrahedra does not match: " << nTuples << " vs " << tuples.size() << std::endl;
    return 1;
  }
  return 0;
}

int main( int argc, char *argv[] )
{
  
  int returnVal = 0;
  returnVal += tetrahedra_test(10000,3,0.2,0.1);
  if(returnVal > 0){
    std::cout << "test 0 failed" << std::endl;
    return returnVal ;
  }
  returnVal += tetrahedra_test(10000,5,0.2,0.1);
  if(returnVal > 0){
    std::cout << "test 1 failed" << std::endl;
    return returnVal ;
  }
  
  
  returnVal += tetrahedron_test_random(20000, 1000, 3, 0.2, 0.1,false);
  if(returnVal > 0){
    std::cout << "test 2 failed" << std::endl;
    return returnVal ;
  }
  returnVal += tetrahedron_test_random(20000, 1000, 5, 0.2, 0.1,false);
  if(returnVal > 0){
    std::cout << "test 3 failed" << std::endl;
    return returnVal ;
  }

  returnVal += tetrahedron_test_random(20000, 50000, 3, 0.2, 0.1,false);
  if(returnVal > 0){
    std::cout << "test 4 failed" << std::endl;
    return returnVal ;
  }
  returnVal += tetrahedron_test_random(20000, 50000, 5, 0.2, 0.1,false);
  if(returnVal > 0){
    std::cout << "test 5 failed" << std::endl;
    return returnVal ;
  }

  returnVal += tetrahedron_test_random(20000, 70000, 3, 0.2, 0.1,true);
  if(returnVal > 0){
    std::cout << "test 6 failed" << std::endl;
    return returnVal ;
  }
  returnVal += tetrahedron_test_random(20000, 70000, 5, 0.2, 0.1,true);
  if(returnVal > 0){
    std::cout << "test 7 failed" << std::endl;
    return returnVal ;
  }


  return returnVal;
}