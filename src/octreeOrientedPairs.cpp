#include "octree.hpp"

namespace octree
{
  /* ##############################################################################
                                searching
                            public functions
  ############################################################################## */
IndexPairVec TOctree::search_oriented_pairs(const  double distance, const double tolerance, const double angleToZ, const double deltaAngle){
  double nPairsEstimation = oriented_pairs_estimation(distance, tolerance,angleToZ,deltaAngle);
  auto pairs = IndexPairVec();
  pairs.reserve(nPairsEstimation);

  const double dist_min = distance * (1. - tolerance);
  const double dist_max = distance * (1. + tolerance);
  const double cosAngle_min = cos(std::min(M_PI*0.5,angleToZ + deltaAngle));
  const double cosAngle_max = cos(std::max(0., angleToZ - deltaAngle));

  size_t level = std::min(depth,(int)trunc(log(boxSize/dist_max)/log(2)) );
  if ( dist_max < boxSize * pow(2,-depth) ){
    for(auto& node: nodesAtLevel[level]){
      direct_search_oriented_pairs(node, dist_min, dist_max, cosAngle_min, cosAngle_max, pairs); 
    }
    return pairs;
  }
  for( auto& node : nodesAtLevel[level]){
    tree_search_oriented_pairs(node,dist_min, dist_max, cosAngle_min, cosAngle_max, pairs);
  }
  return pairs;
}



IndexPairVec TOctree::search_oriented_pairs(const double distance, const double tolerance, const double angleToZ, const double deltaAngle, const size_t maxPairs){
  // returns at most maxPairs pairs
  // the pairs hereby are not isotropic, as the vector between p1 and p2 has an angle
  // angleToZ +/-deltaAngle to the z axis

  // searching the whole tree
  // To avoid multiple counting, we only return tetrahedra with the ordering of indices:
  // particle1 < particle2
  double nPairsEstimation = oriented_pairs_estimation(distance, tolerance,angleToZ,deltaAngle);
  auto pairs = IndexPairVec();

  if (10 * nPairsEstimation < maxPairs ){
  // It is not likely, that we will find more pairs than maxPairs, we thus can search systematically and not randomized at all
    pairs = search_oriented_pairs(distance, tolerance, angleToZ, deltaAngle);
    return pairs;
  }
  const double dist_min = distance * (1. - tolerance);
  const double dist_max = distance * (1. + tolerance);
  const double cosAngle_min = cos(std::min(M_PI*0.5,angleToZ + deltaAngle));
  const double cosAngle_max = cos(std::max(0., angleToZ - deltaAngle));

  pairs.reserve( std::min(maxPairs, (size_t)(nPairsEstimation)) );

  IndexVec particle1Vec(nParticles); //vector of particle1
  std::iota(particle1Vec.begin(),particle1Vec.end(),0); //fill up particle1Vec with ascending numbers
  std::shuffle(particle1Vec.begin(),particle1Vec.end(),randomNumberGenerator);

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

  size_t nPairsFound = 0;
  size_t level = std::min(depth,(int)trunc(log(boxSize/dist_max)/log(2)) );
  for(size_t particle1Idx = 0 ; particle1Idx < nParticles; ++particle1Idx ){
    if(nPairsFound == maxPairs) break;
    tree_search_oriented_pairs_randomized(particle1Vec[particle1Idx],dist_min, dist_max, cosAngle_min, cosAngle_max, level, pairsPerParticle, maxPairs, pairs, nPairsFound);
  }
  return pairs;
}

size_t TOctree::count_oriented_pairs(const double distance, const double tolerance, const double angleToZ, const double deltaAngle){
  size_t nPairs = 0;

  const double dist_min = distance * (1. - tolerance);
  const double dist_max = distance * (1. + tolerance);
  const double cosAngle_min = cos(std::min(M_PI*0.5,angleToZ + deltaAngle));
  const double cosAngle_max = cos(std::max(0., angleToZ - deltaAngle));

  size_t level = std::min(depth,(int)trunc(log(boxSize/dist_max)/log(2)) );
  if ( dist_max < boxSize * pow(2,-depth) ){
    for(auto& node: nodesAtLevel[level]){
    // for(size_t iNode=0; iNode < nNodes; ++iNode){
      nPairs += count_oriented_pairs_direct(node, dist_min, dist_max, cosAngle_min, cosAngle_max); 
    }
    return nPairs;
  }
  for( auto& node : nodesAtLevel[level]){
  // for(size_t iNode=0; iNode < nNodes; ++iNode){
    nPairs += count_oriented_pairs_tree(node,dist_min, dist_max, cosAngle_min, cosAngle_max);
  }
  return nPairs;
}
/* ##############################################################################
                              searching
                          private functions
############################################################################## */


void TOctree::direct_search_oriented_pairs(Node* nodePtr, const double dist_min, const double dist_max,  const double cosAngle_min, const double cosAngle_max, IndexPairVec& pairs){
  // return all pairs with all p1 in Node
  // searches particle wise all relevant Nodes


  // # calculate minimal and maximal distance


  for(auto& particle1 : nodePtr -> get_particles() ){
    for(auto& neighbor: nodePtr -> neighbors ){
      if(! node_is_in_searchspace(neighbor,particle1,dist_min,dist_max)) continue;
      if(! node_is_in_cone(neighbor, particle1, cosAngle_min,cosAngle_max)) continue;
      for(size_t iP2=neighbor->upper_bound(particle1) ; iP2 < neighbor->nParticles; ++iP2 ){
        auto particle2 = *(neighbor -> particlesPtr + iP2) ;
        if( ! particles_at_searched_distance(particle1,particle2, dist_min, dist_max) ) continue;
        if( ! particles_in_z_cone(particle1,particle2, cosAngle_min, cosAngle_max)) continue;
        pairs.push_back({particle1,particle2});
      }  
    }
  }
  return;
}

size_t TOctree::count_oriented_pairs_direct(Node* nodePtr, const double dist_min, const double dist_max,  const double cosAngle_min, const double cosAngle_max){
  size_t nPairs = 0;
  for(auto& particle1 : nodePtr -> get_particles() ){
    for(auto& neighbor: nodePtr -> neighbors ){
      if(! node_is_in_searchspace(neighbor,particle1,dist_min,dist_max)) continue;
      if(! node_is_in_cone(neighbor, particle1, cosAngle_min,cosAngle_max)) continue;
      for(size_t iP2=neighbor->upper_bound(particle1) ; iP2 < neighbor->nParticles; ++iP2 ){
        auto particle2 = *(neighbor -> particlesPtr + iP2) ;
        if( ! particles_at_searched_distance(particle1,particle2, dist_min, dist_max) ) continue;
        if( ! particles_in_z_cone(particle1,particle2, cosAngle_min, cosAngle_max)) continue;
        nPairs++;
      }  
    }
  }
  return nPairs;
}
void TOctree::tree_search_oriented_pairs(Node* nodePtr, const double dist_min, const double dist_max, const double cosAngle_min, const double cosAngle_max, IndexPairVec& pairs){
  // systemic search for all oriented pairs with the first particle in given node.

  auto neighbors =  nodePtr -> neighbors;
   for(auto& particle1 : nodePtr -> get_particles() ){
    for(auto& neighbor: neighbors ){
      auto currentNodeP2 = neighbor;
      do{
        if (! node_is_in_searchspace(currentNodeP2,particle1, dist_min, dist_max) ||
            ! node_is_in_cone(currentNodeP2,particle1,cosAngle_min,cosAngle_max)){
          currentNodeP2 = currentNodeP2 -> nextNode(nodePtr->level);
          continue;
        }
        if (! currentNodeP2->isLeaf){
          currentNodeP2 = currentNodeP2->children[0].get();
          continue;
        }
        // currentNode is in searchspace and LeafNode -> check particles
        for(size_t iP2=currentNodeP2->upper_bound(particle1) ; iP2 < currentNodeP2->nParticles; ++iP2 ){
          auto particle2 = *(currentNodeP2 -> particlesPtr + iP2) ;
          if( ! particles_at_searched_distance(particle1,particle2, dist_min, dist_max) ) continue;
          if( ! particles_in_z_cone(particle1, particle2, cosAngle_min, cosAngle_max)) continue;
          pairs.push_back({particle1,particle2});
        }
        currentNodeP2 = currentNodeP2 -> nextNode(nodePtr->level);
      }while (currentNodeP2->level > nodePtr->level);
    }
  }
  return ;
}

size_t TOctree::count_oriented_pairs_tree(Node* nodePtr, const double dist_min, const double dist_max,  const double cosAngle_min, const double cosAngle_max){
  size_t nPairs = 0;
  auto neighbors =  nodePtr -> neighbors;
   for(auto& particle1 : nodePtr -> get_particles() ){
    for(auto& neighbor: neighbors ){
      auto currentNodeP2 = neighbor;
      do{
        if (! node_is_in_searchspace(currentNodeP2,particle1, dist_min, dist_max) ||
            ! node_is_in_cone(currentNodeP2,particle1,cosAngle_min,cosAngle_max)){
          currentNodeP2 = currentNodeP2 -> nextNode(nodePtr->level);
          continue;
        }
        if (! currentNodeP2->isLeaf){
          currentNodeP2 = currentNodeP2->children[0].get();
          continue;
        }
        // currentNode is in searchspace and LeafNode -> check particles
        for(size_t iP2=currentNodeP2->upper_bound(particle1) ; iP2 < currentNodeP2->nParticles; ++iP2 ){
          auto particle2 = *(currentNodeP2 -> particlesPtr + iP2) ;
          if( ! particles_at_searched_distance(particle1,particle2, dist_min, dist_max) ) continue;
          if( ! particles_in_z_cone(particle1, particle2, cosAngle_min, cosAngle_max)) continue;
          nPairs++;
        }
        currentNodeP2 = currentNodeP2 -> nextNode(nodePtr->level);
      }while (currentNodeP2->level > nodePtr->level);
    }
  }
  return nPairs;
}


void TOctree::tree_search_oriented_pairs_randomized(const size_t particle1, const double dist_min, const double dist_max, const double cosAngle_min, const double cosAngle_max, const size_t level, const size_t maxReturn, const size_t maxPairs, IndexPairVec& pairs, size_t& pairsFound){

  size_t p1PairsFound=0;
  auto startNode = get_search_start_node(particle1, level);
  auto neighbors =  startNode -> neighbors;
  auto searchSpaceL = get_searchspace_large_scale(particle1, neighbors, dist_min,dist_max);
  auto searchSpaceS = get_searchspace_oriented_pair(particle1, searchSpaceL, dist_min,dist_max, cosAngle_min, cosAngle_max) ;

  auto searchSpaceParticles = get_all_particles(searchSpaceS);
  std::shuffle(searchSpaceParticles.begin(), searchSpaceParticles.end(), randomNumberGenerator);
  for(auto& particle2 : searchSpaceParticles){
    if (particle2 <= particle1) continue;
    if (!particles_at_searched_distance(particle1 ,particle2, dist_min, dist_max)) continue;
    if (!particles_in_z_cone(particle1 ,particle2 , cosAngle_min, cosAngle_max)) continue;
    pairs.push_back({particle1,particle2});
    if (++pairsFound   == maxPairs ) return;
    if (++p1PairsFound == maxReturn) return; 
  }
  return ;
}
}
