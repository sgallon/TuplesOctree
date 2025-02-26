#ifndef TESTUPLE
#define TESTUPLE
#include <algorithm>

#include <tuple.hpp>
#include <point.hpp>


bool particles_at_searched_distance(const octree::PointVec& data,const size_t particle1, const size_t particle2,const double dist_min, const double dist_max){
  auto distanceP1P2 = octree::norm_difference(data[particle1], data[particle2] );
  return ( distanceP1P2 > dist_min && distanceP1P2 < dist_max);
}

bool particles_in_z_cone(const octree::PointVec& data,const size_t particle1, const size_t particle2,const double cosAngle_min, const double cosAngle_max){
  auto cosAngle = octree::cosAngleToZAxis(data[particle1], data[particle2] );
  return ( cosAngle > cosAngle_min && cosAngle < cosAngle_max);
}

octree::IndexPairVec direct_search_pairs(const octree::PointVec& data,double distance, double tolerance){
  auto nPoints = data.size();
  auto dist_min = distance * (1-tolerance);
  auto dist_max = distance * (1+tolerance);
  octree::IndexPairVec pairs;
  for(size_t i = 0; i<nPoints; ++i){
    for(size_t j = i+1; j<nPoints; ++j){
      if (! particles_at_searched_distance(data,i,j,dist_min,dist_max)) continue;
      pairs.push_back({i,j});
    }
  }
  return pairs;
}

octree::IndexPairVec direct_search_pairs(const octree::PointVec& data, const double distance, const double tolerance, const double angleToZ, const double deltaAngle){
  auto nPoints = data.size();
  auto dist_min = distance * (1-tolerance);
  auto dist_max = distance * (1+tolerance);
  auto cosAngle_max = cos(std::max(0., angleToZ - deltaAngle));
  auto cosAngle_min = cos(std::min(M_PI*0.5,angleToZ + deltaAngle));
  octree::IndexPairVec pairs;
  for(size_t i = 0; i<nPoints; ++i){
    for(size_t j = i+1; j<nPoints; ++j){
      if (! particles_at_searched_distance(data,i,j,dist_min,dist_max)) continue;
      if (! particles_in_z_cone(data,i,j,cosAngle_min,cosAngle_max)) continue;
      pairs.push_back({i,j});
    }
  }
  return pairs;
}

octree::IndexTripletVec direct_search_triangles(const octree::PointVec& data,double distance, double tolerance){
  auto nPoints = data.size();
  auto dist_min = distance * (1-tolerance);
  auto dist_max = distance * (1+tolerance);
  octree::IndexTripletVec tuple;

  for(size_t ip1 = 0; ip1<nPoints; ++ip1){
    for(size_t ip2 = ip1+1; ip2<nPoints; ++ip2){
      if (! particles_at_searched_distance(data,ip1,ip2,dist_min,dist_max)) continue;
      // loop over ip3 
      for(size_t ip3 = ip2+1; ip3<nPoints; ++ip3){
        if (! particles_at_searched_distance(data,ip1,ip3,dist_min,dist_max)) continue;
        if (! particles_at_searched_distance(data,ip2,ip3,dist_min,dist_max)) continue;
        tuple.push_back({ip1,ip2,ip3});
      }
    }
  }
  return tuple;
}
octree::IndexQuadrupletVec direct_search_tetrahedra(const octree::PointVec& data,double distance, double tolerance){
  auto nPoints = data.size();
  auto dist_min = distance * (1-tolerance);
  auto dist_max = distance * (1+tolerance);
  octree::IndexQuadrupletVec tuple;

  for(size_t ip1 = 0; ip1<nPoints; ++ip1){
    for(size_t ip2 = ip1+1; ip2<nPoints; ++ip2){
      if (! particles_at_searched_distance(data,ip1,ip2,dist_min,dist_max)) continue;
      // loop over ip3 
      for(size_t ip3 = ip2+1; ip3<nPoints; ++ip3){
        if (! particles_at_searched_distance(data,ip1,ip3,dist_min,dist_max)) continue;
        if (! particles_at_searched_distance(data,ip2,ip3,dist_min,dist_max)) continue;
        // loop over ip4
        for(size_t ip4 = ip3+1; ip4<nPoints; ++ip4){
          if (! particles_at_searched_distance(data,ip1,ip4,dist_min,dist_max)) continue;
          if (! particles_at_searched_distance(data,ip2,ip4,dist_min,dist_max)) continue;
          if (! particles_at_searched_distance(data,ip3,ip4,dist_min,dist_max)) continue;
          tuple.push_back({ip1,ip2,ip3,ip4});
        }
      }
    }
  }
  return tuple;
}

int test_pairs(const octree::PointVec& data,const octree::IndexPairVec& pairs,double distance, double tolerance){
  int mismatches = 0;
  auto dist_min = distance * (1-tolerance);
  auto dist_max = distance * (1+tolerance);
  for(size_t i = 0; i<pairs.size(); ++i){
    if(particles_at_searched_distance(data,pairs[i][0],pairs[i][1],dist_min,dist_max)) continue;
    auto dist = octree::norm_difference(data[pairs[i][0]],data[pairs[i][1]]);
    std::cout << "Found a mismatch at pair " << i << " : " << dist/distance << std::endl;
    mismatches++;
  }
  return mismatches;
}


int test_pairs(const octree::PointVec& data,const octree::IndexPairVec& pairs,double distance, double tolerance, double angleToZ, double deltaAngle){
  int mismatches = 0;
  auto dist_min = distance * (1-tolerance);
  auto dist_max = distance * (1+tolerance);
  auto cosAngle_min = std::cos(std::min(0.5*M_PI, angleToZ + deltaAngle));
  auto cosAngle_max = std::cos(std::max(0.,angleToZ - deltaAngle));
  for(size_t i = 0; i<pairs.size(); ++i){
    if(particles_at_searched_distance(data,pairs[i][0],pairs[i][1],dist_min,dist_max) &&
       particles_in_z_cone(data,pairs[i][0],pairs[i][1],cosAngle_min,cosAngle_max)){
        continue;
    }
    std::cout << "Found a mismatch at pair " << i << " : ";
    if(! particles_at_searched_distance(data,pairs[i][0],pairs[i][1],dist_min,dist_max)){
      auto dist = octree::norm_difference(data[pairs[i][0]],data[pairs[i][1]]);
      std::cout << "distance " << dist/distance  << std::endl; 
    }
    if(! particles_in_z_cone(data,pairs[i][0],pairs[i][1],cosAngle_min,cosAngle_max)){
      auto cosAngle = octree::cosAngleToZAxis(data[pairs[i][0]],data[pairs[i][1]]);
      std::cout << "angle " << std::acos(cosAngle) * 180/M_PI << std::endl; 
    }
    mismatches++;
  }
  return mismatches;
}

int test_triangles(const octree::PointVec& data,const octree::IndexTripletVec& triplets,double distance, double tolerance){
  int mismatches = 0;
  auto dist_min = distance * (1-tolerance);
  auto dist_max = distance * (1+tolerance);
  for(size_t i = 0; i<triplets.size(); ++i){
    if(particles_at_searched_distance(data,triplets[i][0],triplets[i][1],dist_min,dist_max)) continue;
    if(particles_at_searched_distance(data,triplets[i][0],triplets[i][2],dist_min,dist_max)) continue;
    if(particles_at_searched_distance(data,triplets[i][1],triplets[i][2],dist_min,dist_max)) continue;
    auto dist12 = octree::norm_difference(data[triplets[i][0]],data[triplets[i][1]]);
    auto dist13 = octree::norm_difference(data[triplets[i][0]],data[triplets[i][2]]);
    auto dist23 = octree::norm_difference(data[triplets[i][1]],data[triplets[i][2]]);
    std::cout << "Found a mismatch at triplet " << i << " : " << dist12/distance << " " << dist13/distance << " " << dist23/distance << std::endl;
    mismatches++;
  }
  return mismatches;
}

int test_tetrahedra(const octree::PointVec& data,const octree::IndexQuadrupletVec& quadruplets,double distance, double tolerance){
  int mismatches = 0;
  auto dist_min = distance * (1-tolerance);
  auto dist_max = distance * (1+tolerance);
  for(size_t i = 0; i<quadruplets.size(); ++i){
    if(particles_at_searched_distance(data,quadruplets[i][0],quadruplets[i][1],dist_min,dist_max)) continue;
    if(particles_at_searched_distance(data,quadruplets[i][0],quadruplets[i][2],dist_min,dist_max)) continue;
    if(particles_at_searched_distance(data,quadruplets[i][0],quadruplets[i][3],dist_min,dist_max)) continue;
    if(particles_at_searched_distance(data,quadruplets[i][1],quadruplets[i][2],dist_min,dist_max)) continue;
    if(particles_at_searched_distance(data,quadruplets[i][1],quadruplets[i][3],dist_min,dist_max)) continue;
    if(particles_at_searched_distance(data,quadruplets[i][2],quadruplets[i][3],dist_min,dist_max)) continue;
    // calculate all pairwise distances
    auto dist12 = octree::norm_difference(data[quadruplets[i][0]],data[quadruplets[i][1]]);
    auto dist13 = octree::norm_difference(data[quadruplets[i][0]],data[quadruplets[i][2]]);
    auto dist14 = octree::norm_difference(data[quadruplets[i][0]],data[quadruplets[i][3]]);
    auto dist23 = octree::norm_difference(data[quadruplets[i][1]],data[quadruplets[i][2]]);
    auto dist24 = octree::norm_difference(data[quadruplets[i][1]],data[quadruplets[i][3]]);
    auto dist34 = octree::norm_difference(data[quadruplets[i][2]],data[quadruplets[i][3]]);
    std::cout << "Found a mismatch at quadruplet " << i << " : " << dist12/distance << " " << dist13/distance << " " << dist14/distance << " " << dist23/distance << " " << dist24/distance << " " << dist34/distance << std::endl;
    mismatches++;
  }
  return mismatches;
}

void sort_pair(octree::IndexPair& pair){
  if(pair[0] > pair[1]){
    auto temp = pair[0];
    pair[0] = pair[1];
    pair[1] = temp;
  }
}

void sort_pairs(octree::IndexPairVec& pairs){
  auto nPairs = pairs.size();
  // first sort each pair individually
  for(size_t i = 0; i<nPairs; ++i){
    sort_pair(pairs[i]);
  }
  // then sort the whole vector
  std::sort(pairs.begin(),pairs.end());
}

void sort_triplet(octree::IndexTriplet& triplet){
  if(triplet[0] > triplet[1]){
    auto temp = triplet[0];
    triplet[0] = triplet[1];
    triplet[1] = temp;
  }
  if(triplet[0] > triplet[2]){
    auto temp = triplet[0];
    triplet[0] = triplet[2];
    triplet[2] = temp;
  }
  if(triplet[1] > triplet[2]){
    auto temp = triplet[1];
    triplet[1] = triplet[2];
    triplet[2] = temp;
  }
}

void sort_triplets(octree::IndexTripletVec& triplets){
  auto nTriplets = triplets.size();
  // first sort each triplet individually
  for(size_t i = 0; i<nTriplets; ++i){
    sort_triplet(triplets[i]);
  }
  // then sort the whole vector
  std::sort(triplets.begin(),triplets.end());
}

void sort_quadruplet(octree::IndexQuadruplet& quadruplet){
  if(quadruplet[0] > quadruplet[1]){
    auto temp = quadruplet[0];
    quadruplet[0] = quadruplet[1];
    quadruplet[1] = temp;
  }
  if(quadruplet[0] > quadruplet[2]){
    auto temp = quadruplet[0];
    quadruplet[0] = quadruplet[2];
    quadruplet[2] = temp;
  }
  if(quadruplet[0] > quadruplet[3]){
    auto temp = quadruplet[0];
    quadruplet[0] = quadruplet[3];
    quadruplet[3] = temp;
  }
  if(quadruplet[1] > quadruplet[2]){
    auto temp = quadruplet[1];
    quadruplet[1] = quadruplet[2];
    quadruplet[2] = temp;
  }
  if(quadruplet[1] > quadruplet[3]){
    auto temp = quadruplet[1];
    quadruplet[1] = quadruplet[3];
    quadruplet[3] = temp;
  }
  if(quadruplet[2] > quadruplet[3]){
    auto temp = quadruplet[2];
    quadruplet[2] = quadruplet[3];
    quadruplet[3] = temp;
  }
}

void sort_quadruplets(octree::IndexQuadrupletVec& quadruplets){
  auto nQuadrupletts = quadruplets.size();
  // first sort each quadruplet individually
  for(size_t i = 0; i<nQuadrupletts; ++i){
    sort_quadruplet(quadruplets[i]);
  }
  // then sort the whole vector
  std::sort(quadruplets.begin(),quadruplets.end());
}

template <size_t tupleType>
void sort_tuple(octree::IndexNTuple<tupleType>& tuple){
  if(tupleType == 2){
    sort_pair(tuple);
  }else if(tupleType == 3){
    sort_triplet(tuple);
  }else if(tupleType == 4){
    sort_quadruplet(tuple);
  }
}

template <size_t tupleType>
int count_duplicates (octree::IndexNTupleVec<tupleType>& tuples){
  auto nTuples = tuples.size();
  if (nTuples == 0) return 0;
  int nDuplicates = 0;
  for(size_t i = 0; i<nTuples-1; ++i){
    for(size_t j = i+1; j<nTuples; ++j){
      if(tuples[i] == tuples[j]){
        std::cout << "Found duplicate at " << i << " and " << j << std::endl;
        nDuplicates++;
      }
    }
  }
  return nDuplicates;
}

template <size_t tupleType>
int count_max_tuplePerParticle1(octree::IndexNTupleVec<tupleType>& tuples){
  auto nTuples = tuples.size();
  size_t maxTuplesPerParticle1 = 0;
  for(size_t i = 0; i<nTuples-1; ++i){
    size_t nTuplesPerParticle1 = 0;
    for(size_t j = i; j<nTuples; ++j){
      if(tuples[i][0] == tuples[j][0]){
        nTuplesPerParticle1++;
      }
    }
    if(nTuplesPerParticle1 > maxTuplesPerParticle1){
      maxTuplesPerParticle1 = nTuplesPerParticle1;
    }
  }
  return maxTuplesPerParticle1;
}

#endif