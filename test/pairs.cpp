#include <iostream>
#include "randompoints.cpp"

#include <octree.hpp>
#include <tuple.hpp>

#include <algorithm>
#include <chrono>
#include "testtuples.cpp"


int pair_test_random( size_t N,size_t maxPairs , int depth, double distance, double tolerance,  bool doDirectSearch = true)
{
  auto randomPoints = generateRandomPointArr(N);
  
  octree::IndexPairVec pairs_direct;
  if(doDirectSearch){
    pairs_direct = direct_search_pairs(randomPoints,distance,tolerance);
    if(pairs_direct.size()==0){
      std::cout << "No pairs found in direct search. Update parameter" << std::endl;
      return 1;
    }
  }

  auto tree = octree::TOctree(randomPoints,{1.,1.,1.},2.,depth);
  auto pairs = tree.search_pairs(distance,tolerance,maxPairs);

  if((doDirectSearch) && (pairs.size() < pairs_direct.size())){
    std::cout << "Number of pairs in direct search: " << pairs_direct.size() << " vs " << pairs.size() << std::endl;
  }
  if (pairs.size() == 0) {
    std::cout << "No triangles found" << std::endl;
    return 1;
  }  
  std::cout << "Searched for " << maxPairs << " pairs, found " << pairs.size() << std::endl;

  auto mismatch = test_pairs(randomPoints,pairs,distance,tolerance);
  if(mismatch != 0){
    std::cout << "number of mismatches: " << mismatch << std::endl;
    return mismatch;
  }else{
    std::cout << "no mismatches found" << std::endl;
  }
  auto nDuplicates = count_duplicates(pairs);
  if(nDuplicates != 0){
    std::cout << "number of duplicates: " << nDuplicates << std::endl;
    return nDuplicates;
  }else{
    std::cout << "no duplicates found" << std::endl;
  }
  if((doDirectSearch) && (pairs_direct.size() > maxPairs) && (pairs.size() < 0.9 *maxPairs)){
    std::cout << "Tree search found less than 90% of pairs than maxPairs, but direct search found more than maxPairs" << std::endl;
    return 1;
  }
  auto maxTuplePerP1 = count_max_tuplePerParticle1(pairs);
  if((doDirectSearch) && 
     (pairs_direct.size() > 5 * maxPairs) &&
     (maxTuplePerP1 > std::max(1.,(double) pairs_direct.size()/N)) ){
    std::cout << "Unsufficent randomization of tuples" << std::endl;
    std::cout << "maxTuplePerP1: " << maxTuplePerP1 << std::endl;
    return 1;
  }
  return 0;
}

int pair_test( size_t N, int depth, double distance, double tolerance)
{
  auto randomPoints = generateRandomPointArr(N);
  auto tree = octree::TOctree(randomPoints,{1.,1.,1.},2.,depth);

  auto pairs = tree.search_pairs(distance,tolerance);
  std::cout << "Number of pairs: " << pairs.size() << std::endl;

  // test for wrong pairs
  auto mismatch = test_pairs(randomPoints,pairs,distance,tolerance);
  if(mismatch != 0){
    std::cout << "number of mismatches: " << mismatch << std::endl;
    return mismatch;
  }else{
    std::cout << "no mismatches found" << std::endl;
  }
  // test for duplicates
  auto nDuplicates = count_duplicates(pairs);
  if(nDuplicates != 0){
    std::cout << "number of duplicates: " << nDuplicates << std::endl;
    return nDuplicates;
  }else{
    std::cout << "no duplicates found" << std::endl;
  }
  // test for missing pairs
  auto pairs_direct = direct_search_pairs(randomPoints,distance,tolerance);
  if(pairs.size() != pairs_direct.size()){
    std::cout << "Number of pairs does not match: " << pairs.size() << " vs " << pairs_direct.size() << std::endl;
    sort_pairs(pairs);
    size_t iMax = std::min(pairs.size(),pairs_direct.size());
    for(size_t i = 0; i<iMax ; ++i){
      if(pairs[i] != pairs_direct[i]){
        std::cout << "Pair " << i << " does not match: " << pairs[i][0] << " " << pairs[i][1] << " vs " << pairs_direct[i][0] << " " << pairs_direct[i][1] << std::endl;
      }
    }
    return 1;
  }
  return 0;
}

double test_count(size_t N, int depth, double distance, double tolerance){
  // generate random points
  auto randomPoints = generateRandomPointArr(N);
  // generate octree
  auto tree = octree::TOctree(randomPoints,{1.,1.,1.},2.,depth);
  // count pairs
  auto pairs = tree.search_pairs(distance,tolerance);
  // count pairs
  auto nPairs = tree.count_pairs(distance,tolerance);
  // compare
  if(nPairs != pairs.size()){
    std::cout << "Number of pairs does not match: " << nPairs << " vs " << pairs.size() << std::endl;
    return 1;
  }
  return 0;
}


int main( int argc, char *argv[] )
{
  int returnVal = 0;
  returnVal += pair_test(10000, 3, 0.2, 1e-3);
  if(returnVal > 0){
    std::cout << "test 0 failed" << std::endl;
    return returnVal ;
  }
  returnVal += pair_test(10000, 5, 0.2, 1e-3);
  if(returnVal > 0){
    std::cout << "test 1 failed" << std::endl;
    return returnVal ;
  }

  returnVal += pair_test_random(100000, 1000, 3, 0.2, 1e-3,false);
  if(returnVal > 0){
    std::cout << "test 2 failed" << std::endl;
    return returnVal ;
  }
  returnVal += pair_test_random(100000, 1000, 5, 0.2, 1e-3,false);
  if(returnVal > 0){
    std::cout << "test 3 failed" << std::endl;
    return returnVal ;
  }

  returnVal += pair_test_random(100000, 70000, 3, 0.2, 1e-3,false);
  if(returnVal > 0){
    std::cout << "test 4 failed" << std::endl;
    return returnVal ;
  }
  returnVal += pair_test_random(100000, 70000, 5, 0.2, 1e-3,false);
  if(returnVal > 0){
    std::cout << "test 5 failed" << std::endl;
    return returnVal ;
  }

  returnVal += pair_test_random(100000, 108000, 3, 0.2, 1e-3,true);
  if(returnVal > 0){
    std::cout << "test 6 failed" << std::endl;
    return returnVal ;
  }
  returnVal += pair_test_random(100000, 108000, 5, 0.2, 1e-3,true);
  if(returnVal > 0){
    std::cout << "test 7 failed" << std::endl;
    return returnVal ;
  }
  

  
  return returnVal;
}