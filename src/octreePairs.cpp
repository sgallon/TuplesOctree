#include "octree.hpp"

namespace octree
{
  /* ##############################################################################
                                searching
                            public functions
  ############################################################################## */
IndexPairVec TOctree::search_pairs(const  double distance, const double tolerance){
  double nPairsEstimation = pairs_estimation(distance, tolerance);
  auto pairs = IndexPairVec();
  pairs.reserve( (size_t)(1.5 * nPairsEstimation) ); //reserve a bit more space, to avoid reallocation

  const double dist_min = distance * (1. - tolerance);
  const double dist_max = distance * (1. + tolerance);

  size_t level = std::min(depth,(int)trunc(log(boxSize/dist_max)/log(2)) ); // level of the tree at which we start searching
  if ( dist_max < boxSize * pow(2,-depth) ){ // for small distances, we search only in the leaf nodes
    for(auto& nodePtr : nodesAtLevel[level]){
      direct_search_pairs(nodePtr,dist_min,dist_max, pairs);
    }
    return pairs;
  }
  for(auto& nodePtr : nodesAtLevel[level]){
    tree_search_pairs(nodePtr,dist_min,dist_max, pairs);
  }
  return pairs;
}

IndexPairVec TOctree::search_pairs(const  double distance, const double tolerance,const size_t maxPairs){
// returns at most maxPairs pairs
// searching the whole tree
// To avoid multiple counting, we only return tetrahedra with the ordering of indices:
// particle1 < particle2

double nPairsEstimation = pairs_estimation(distance, tolerance);
if (10 * nPairsEstimation < maxPairs ){
  // If we expect less than 10% of the maxPairs, we can search all pairs
  return search_pairs(distance,tolerance);
}

auto dist_min = distance * (1. - tolerance);
auto dist_max = distance * (1. + tolerance);

auto pairs = IndexPairVec();
pairs.reserve( (size_t)(std::min(nPairsEstimation, (double) maxPairs)) );

IndexVec particle1Vec(nParticles);
std::iota(particle1Vec.begin(), particle1Vec.end(), 0);


size_t pairsPerParticle = 0;
if (  0.75 * (size_t)(nPairsEstimation) < maxPairs) {
    // We are expecting roughly the same amount of Pairs as desired
    // Hence, Pairs are still "rare"
    // Therefore, we need to search in a semi-random way: returning all Pairs given a particle1
    // stopping only if we have found enough Pairs or if we have searched all particles
  pairsPerParticle =  maxPairs ;
}
else if( 0.5 > maxPairs/(double) nParticles){
  // We are expecting more Pairs than desired
  // but we are still expecting less than 50% of the maxPairs per particle
  pairsPerParticle = std::max((size_t) 1, (size_t) (1.5 * round( (double) maxPairs / (double) nParticles ) )) ;
}
else{
  // We are expecting more Pairs than desired
  // and we are expecting more than 50% of the maxPairs per particle
  pairsPerParticle = std::max((size_t) 1 ,(size_t) round( 1.5 * (double) nPairsEstimation / (double) nParticles ) ) ;
  // return slightly more PairsPerParticle than expected, to avoid finding to few, as for high indices, there are fewer Pairs with strictly increasing indices
}
if ( distance < boxSize * pow(2,-depth) ) {
  for(auto& particle1 : particle1Vec){
    if(pairs.size() == maxPairs) break;
    direct_search_pairs_randomized(particle1 , dist_min, dist_max, pairsPerParticle, maxPairs, pairs);
  }
  return pairs;
}
for(auto& particle1 : particle1Vec){
  if(pairs.size() == maxPairs) break;
  tree_search_pairs_randomized(particle1,dist_min, dist_max, pairsPerParticle,maxPairs, pairs); 
}
return pairs;
}

size_t TOctree::count_pairs(const double distance, const double tolerance){
  const double dist_min = distance * (1. - tolerance);
  const double dist_max = distance * (1. + tolerance);
  double nPairs=0;
  size_t level = std::min(depth,(int)trunc(log(boxSize/dist_max)/log(2)) ); // level of the tree at which we start searching
  if ( dist_max < boxSize * pow(2,-depth) ){ // for small distances, we search only in the leaf nodes
    for(auto& nodePtr : nodesAtLevel[level]){
      nPairs += count_pairs_direct(nodePtr, dist_min, dist_max);
    }
    return nPairs;
  }
  for(auto& nodePtr : nodesAtLevel[level]){
    nPairs += count_pairs_tree(nodePtr,dist_min,dist_max);
  }
  return nPairs;
}


/* ##############################################################################
                              searching
                          private functions
############################################################################## */


void TOctree::direct_search_pairs(Node* nodePtr, const double dist_min,  const double dist_max, IndexPairVec& pairs){
  // search all pairs with particle1 in nodePtr
  // and particle2 in nodePtr or any of its neighbors
  for(auto& particle1 : nodePtr->particles){
    for (auto& neighbor : nodePtr->neighbors){
      if(! node_is_in_searchspace(neighbor,particle1, dist_min, dist_max)) continue;
      for(size_t iP2=neighbor->upper_bound(particle1) ; iP2 < neighbor->nParticles; ++iP2 ){
        auto particle2 = *(neighbor -> particlesPtr + iP2) ;
        if( ! particles_at_searched_distance(particle1,particle2, dist_min, dist_max) ) continue;
        pairs.push_back({particle1,particle2});
      }
    }
  }
  return;
}

size_t TOctree::count_pairs_direct(Node* nodePtr, const double dist_min,  const double dist_max){
  size_t nPairsNode = 0;
  for(auto& particle1 : nodePtr->particles){
    for (auto& neighbor : nodePtr->neighbors){
      if(! node_is_in_searchspace(neighbor,particle1, dist_min, dist_max)) continue;
      for(size_t iP2=neighbor->upper_bound(particle1) ; iP2 < neighbor->nParticles; ++iP2 ){
        auto particle2 = *(neighbor -> particlesPtr + iP2) ;
        if( ! particles_at_searched_distance(particle1,particle2, dist_min, dist_max) ) continue;
        nPairsNode++;
      }
    }
  }
  return nPairsNode;

}


void TOctree::tree_search_pairs(Node* nodePtr, const double distance_min, const double distance_max, IndexPairVec& pairs){

  
  for(auto& particle1 : nodePtr->particles){
    // for (auto& neighbor : neighbors){
    for (auto& neighbor : nodePtr->neighbors){
      auto currentNodeP2 = neighbor;
      do{
        if(! node_is_in_searchspace(currentNodeP2,particle1, distance_min, distance_max)){
          currentNodeP2 = currentNodeP2 -> nextNode(nodePtr->level);
          continue;
        }
        if (! currentNodeP2->isLeaf){
          currentNodeP2 = currentNodeP2->children[0].get();
          continue;
        }
        // currentNode is in searchspace and leaf node -> check particles
        for(size_t iP2=currentNodeP2->upper_bound(particle1) ; iP2 < currentNodeP2->nParticles; ++iP2 ){
          auto particle2 = *(currentNodeP2 -> particlesPtr + iP2) ;
          if( ! particles_at_searched_distance(particle1,particle2, distance_min, distance_max) ) continue;
          pairs.push_back({particle1,particle2});
        }
        currentNodeP2 = currentNodeP2 -> nextNode(nodePtr->level);
      }while (currentNodeP2->level > nodePtr->level);
    }
  }
  return ;
}

size_t TOctree::count_pairs_tree(Node* nodePtr, const double distance_min, const double distance_max){
  size_t nPairs = 0;
  for(auto& particle1 : nodePtr->particles){
    // for (auto& neighbor : neighbors){
    for (auto& neighbor : nodePtr->neighbors){
      auto currentNodeP2 = neighbor;
      do{
        if(! node_is_in_searchspace(currentNodeP2,particle1, distance_min, distance_max)){
          currentNodeP2 = currentNodeP2 -> nextNode(nodePtr->level);
          continue;
        }
        if (! currentNodeP2->isLeaf){
          currentNodeP2 = currentNodeP2->children[0].get();
          continue;
        }
        // currentNode is in searchspace and leaf node -> check particles
        for(size_t iP2=currentNodeP2->upper_bound(particle1) ; iP2 < currentNodeP2->nParticles; ++iP2 ){
          auto particle2 = *(currentNodeP2 -> particlesPtr + iP2) ;
          if( ! particles_at_searched_distance(particle1,particle2, distance_min, distance_max) ) continue;
          nPairs++;
        }
        currentNodeP2 = currentNodeP2 -> nextNode(nodePtr->level);
      }while (currentNodeP2->level > nodePtr->level);
    }
  }
  return nPairs ;
}


void TOctree::direct_search_pairs_randomized(const size_t particle1, const double dist_min, const double dist_max, const size_t maxReturn, const size_t maxPairs, IndexPairVec& pairs){
  size_t nPairsFound = 0;

  auto nodeP1 = get_starting_node(particle1,dist_max);

  auto searchSpace = get_searchspace_large_scale(particle1,nodeP1 -> neighbors,dist_min,dist_max);

  auto searchSpaceParticles = get_all_particles(searchSpace);
  std::shuffle(searchSpaceParticles.begin(), searchSpaceParticles.end(), randomNumberGenerator);

  for(auto& particle2 : searchSpaceParticles){
    if (particle2 <= particle1) continue;
    if( ! particles_at_searched_distance(particle1,particle2, dist_min, dist_max) ) continue;
    pairs.push_back({particle1,particle2});
    if (++nPairsFound == maxReturn) break;
    if (pairs.size() == maxPairs) break;
  }
  return ;
}

void TOctree::tree_search_pairs_randomized(const size_t particle1,  const double dist_min, const double dist_max, const size_t maxReturn, const size_t maxPairs, IndexPairVec& pairs){
  size_t nPairsFound = 0;

  auto nodeP1 = get_starting_node(particle1,dist_max);



  auto searchSpaceL = get_searchspace_large_scale(particle1,nodeP1 ->neighbors,dist_min,dist_max);
  auto searchSpaceS = get_searchspace_pair(particle1,searchSpaceL,dist_min,dist_max) ;

  auto searchSpaceParticles = get_all_particles(searchSpaceS);
  std::shuffle(searchSpaceParticles.begin(),searchSpaceParticles.end(),randomNumberGenerator);


  for(auto& particle2 : searchSpaceParticles){
    if (particle2 <= particle1) continue;
    if (! particles_at_searched_distance(particle1,particle2, dist_min,dist_max) ) continue;
    pairs.push_back({particle1,particle2});
    if (++nPairsFound == maxReturn) break;
    if (pairs.size() == maxPairs) break;
  }
  return ;
}
}
