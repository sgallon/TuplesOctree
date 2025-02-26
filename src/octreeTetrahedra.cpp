#include "octree.hpp"

namespace octree
{
  /* ##############################################################################
                                searching
                            public functions
  ############################################################################## */
IndexQuadrupletVec TOctree::search_tetrahedra(const double distance,const double tolerance){
  // returns all tetrahedra
  // searching the whole tree
  // To avoid multiple counting, we only return tetrahedra with the ordering of indices:
  // particle1 < particle2 < particle3 < particle4
  double nTetrahedraEstimation = tetrahedra_estimation(distance, tolerance);
  auto tetrahedra = IndexQuadrupletVec();
  tetrahedra.reserve( (size_t)(nTetrahedraEstimation) );
  const double dist_min = distance * (1. - tolerance);
  const double dist_max = distance * (1. + tolerance);
  size_t level = std::min(depth,(int)trunc(log(boxSize/dist_max)/log(2)) );
  if ( dist_max < boxSize * pow(2,-depth) ){
    // for small distances, we search only in the leaf nodes
    for(auto& node : nodesAtLevel[level]){
      direct_search_tetrahedra(node, dist_min, dist_max, tetrahedra);
    }
    return tetrahedra;
  }
  for(auto& node : nodesAtLevel[level]){
    tree_search_tetrahedra(node,dist_min,dist_max, tetrahedra);
  }  
  return tetrahedra;
}

  IndexQuadrupletVec TOctree::search_tetrahedra(const double distance,const double tolerance, const size_t maxTetrahedra){
  // returns at most maxTetrahedra tetrahedra
  // searching the whole tree
  // To avoid multiple counting, we only return tetrahedra with the ordering of indices:
  // particle1 < particle2 < particle3 < particle4
  double nTetrahedraEstimation = tetrahedra_estimation(distance, tolerance);
  auto tetrahedra = IndexQuadrupletVec();


  if (10 * nTetrahedraEstimation < maxTetrahedra ){
    return search_tetrahedra(distance,tolerance);
  }
  // it is not-unlikely that we find more tetrahedra than desired
  // the algorithm needs to be randomized to obtain a stochastically representative subset
  // instead of systematically searching the Nodes by their position (smaller computational costs, as eg. neighborhood query needs to be called once per Node and not once per particle), we search systematically by varying particle1
  tetrahedra.reserve( std::min(maxTetrahedra, (size_t)(nTetrahedraEstimation)) );

  auto distMin = distance * (1. - tolerance);
  auto distMax = distance * (1. + tolerance);

  IndexVec particle1Vec(nParticles); //vector of particle1
  std::iota(particle1Vec.begin(),particle1Vec.end(),0); //fill up particle1Vec with ascending numbers
  std::shuffle(particle1Vec.begin(),particle1Vec.end(),randomNumberGenerator);

  // This is not necessary if we can assume that particles are sufficiently mixed and that the position of a particle is statistically independent of its index

  // variables for progress bar

  size_t tetrahedraPerParticle = 0;
  if ( 0.75* nTetrahedraEstimation < maxTetrahedra) {
    tetrahedraPerParticle =  maxTetrahedra ;
  }
  else if( 0.5  < maxTetrahedra/(double) nParticles){
    tetrahedraPerParticle = std::max( 1,(int) (3* round((double) maxTetrahedra/(double) nParticles)) );
  }
  else{
    tetrahedraPerParticle = std::max(1, (int)(2 * (double) nTetrahedraEstimation / (double) nParticles ) ) ;
  }
  if ( distance < boxSize * pow(2,-depth) ) {
    for(auto& particle1 : particle1Vec){
      if(tetrahedra.size() == maxTetrahedra) break;
      direct_search_tetrahedra_randomized(particle1,distMin,distMax,tetrahedraPerParticle, maxTetrahedra, tetrahedra);
    }
    return tetrahedra;
  }
  for(auto& particle1 : particle1Vec){
    if(tetrahedra.size() == maxTetrahedra) break;
    tree_search_tetrahedra_randomized(particle1,distMin,distMax,tetrahedraPerParticle, maxTetrahedra, tetrahedra);
  }
  return tetrahedra;
}

size_t TOctree::count_tetrahedra( const double distance, const double tolerance){
  size_t nTetrahedra = 0;
  const double dist_min = distance * (1. - tolerance);
  const double dist_max = distance * (1. + tolerance);
  size_t level = std::min(depth,(int)trunc(log(boxSize/dist_max)/log(2)) );
  if ( dist_max < boxSize * pow(2,-depth) ){
    // for small distances, we search only in the leaf nodes
    for(auto& node : nodesAtLevel[level]){
      count_tetrahedra_direct(node, dist_min, dist_max);
    }
    return nTetrahedra;
  }
  for(auto& node : nodesAtLevel[level]){
    count_tetrahedra_tree(node,dist_min,dist_max);
  }  
  return nTetrahedra;
}
/* ##############################################################################
                              searching
                          private functions
############################################################################## */
void TOctree::direct_search_tetrahedra(Node* nodePtr, const double distMin, const double distMax, IndexQuadrupletVec& tetrahedra){
  // return all tetrahedra with all p1 in Node
  // searches particle wise all relevant Nodes
  for(auto& particle1 : nodePtr -> get_particles()){
    for (auto& neighbor : nodePtr -> neighbors){
      if(! node_is_in_searchspace(neighbor, particle1, distMin, distMax) )  continue;
      auto iP2start = neighbor->upper_bound(particle1);
      for(size_t iP2=iP2start ; iP2 < neighbor->nParticles; ++iP2 ){
        auto particle2 = *(neighbor -> particlesPtr + iP2) ;
        if( ! particles_at_searched_distance(particle1,particle2, distMin, distMax) ) continue;
        // particle 1, particle2 is a pair
        for(auto& neighbor3 : nodePtr -> neighbors){
          if(! node_is_in_searchspace(neighbor3, particle2, distMin, distMax) )  continue;
          if(! node_is_in_searchspace(neighbor3, particle1, distMin, distMax) )  continue;
          auto iP3start = neighbor3->upper_bound(particle2);
          for( auto iP3 = iP3start ; iP3 < neighbor3 -> nParticles; ++iP3){
            auto particle3 = *(neighbor3 -> particlesPtr + iP3) ;
            if( ! particles_at_searched_distance(particle2,particle3, distMin, distMax) ) continue;
            if( ! particles_at_searched_distance(particle1,particle3, distMin, distMax) ) continue;
            // particle1, particle2, particle3 is a regular triangle
            for(auto& neighbor4 : nodePtr -> neighbors){
              if(! node_is_in_searchspace(neighbor4, particle3, distMin, distMax) )  continue;
              if(! node_is_in_searchspace(neighbor4, particle2, distMin, distMax) )  continue;
              if(! node_is_in_searchspace(neighbor4, particle1, distMin, distMax) )  continue;
              auto iP4start = neighbor4->upper_bound(particle3);
              for( auto iP4 = iP4start ; iP4 < neighbor4 -> nParticles; ++iP4){
                auto particle4 = *(neighbor4 -> particlesPtr + iP4) ;
                if( ! particles_at_searched_distance(particle3,particle4, distMin, distMax) ) continue;
                if( ! particles_at_searched_distance(particle2,particle4, distMin, distMax) ) continue;
                if( ! particles_at_searched_distance(particle1,particle4, distMin, distMax) ) continue;
                // particle1, particle2, particle3, particle4 is a tetrahedron
                tetrahedra.push_back({particle1,particle2,particle3,particle4});
              }
            }
          }
        }
      }
    }
  }
  return ;
}

size_t TOctree::count_tetrahedra_direct(Node* nodePtr, const double distMin, const double distMax){
  // return all tetrahedra with all p1 in Node
  // searches particle wise all relevant Nodes
  size_t nTetrahedra= 0;
  for(auto& particle1 : nodePtr -> get_particles()){
    for (auto& neighbor : nodePtr -> neighbors){
      if(! node_is_in_searchspace(neighbor, particle1, distMin, distMax) )  continue;
      auto iP2start = neighbor->upper_bound(particle1);
      for(size_t iP2=iP2start ; iP2 < neighbor->nParticles; ++iP2 ){
        auto particle2 = *(neighbor -> particlesPtr + iP2) ;
        if( ! particles_at_searched_distance(particle1,particle2, distMin, distMax) ) continue;
        // particle 1, particle2 is a pair
        for(auto& neighbor3 : nodePtr -> neighbors){
          if(! node_is_in_searchspace(neighbor3, particle2, distMin, distMax) )  continue;
          if(! node_is_in_searchspace(neighbor3, particle1, distMin, distMax) )  continue;
          auto iP3start = neighbor3->upper_bound(particle2);
          for( auto iP3 = iP3start ; iP3 < neighbor3 -> nParticles; ++iP3){
            auto particle3 = *(neighbor3 -> particlesPtr + iP3) ;
            if( ! particles_at_searched_distance(particle2,particle3, distMin, distMax) ) continue;
            if( ! particles_at_searched_distance(particle1,particle3, distMin, distMax) ) continue;
            // particle1, particle2, particle3 is a regular triangle
            for(auto& neighbor4 : nodePtr -> neighbors){
              if(! node_is_in_searchspace(neighbor4, particle3, distMin, distMax) )  continue;
              if(! node_is_in_searchspace(neighbor4, particle2, distMin, distMax) )  continue;
              if(! node_is_in_searchspace(neighbor4, particle1, distMin, distMax) )  continue;
              auto iP4start = neighbor4->upper_bound(particle3);
              for( auto iP4 = iP4start ; iP4 < neighbor4 -> nParticles; ++iP4){
                auto particle4 = *(neighbor4 -> particlesPtr + iP4) ;
                if( ! particles_at_searched_distance(particle3,particle4, distMin, distMax) ) continue;
                if( ! particles_at_searched_distance(particle2,particle4, distMin, distMax) ) continue;
                if( ! particles_at_searched_distance(particle1,particle4, distMin, distMax) ) continue;
                // particle1, particle2, particle3, particle4 is a tetrahedron
                nTetrahedra++;
              }
            }
          }
        }
      }
    }
  }
  return nTetrahedra ;
}

void TOctree::tree_search_tetrahedra(Node* nodePtr, const double distMin, const double distMax, IndexQuadrupletVec& tetrahedra){
  // return all tetrahedra with all p1 in Node
  // searches particle wise all relevant Nodes
  for(auto& particle1 : nodePtr->particles){
    for (auto& neighbor : nodePtr->neighbors){
      auto currentNodeP2 = neighbor;
      do{
        if(! node_is_in_searchspace(currentNodeP2,particle1, distMin, distMax)){
          currentNodeP2 = currentNodeP2 -> nextNode(nodePtr->level);
          continue;
        }
        if (! currentNodeP2->isLeaf){
          currentNodeP2 = currentNodeP2->children[0].get();
          continue;
        }
        auto iP2start = currentNodeP2->upper_bound(particle1);
        for(size_t iP2=iP2start ; iP2 < currentNodeP2->nParticles; ++iP2 ){
          auto particle2 = *(currentNodeP2 -> particlesPtr + iP2) ;
          if( ! particles_at_searched_distance(particle1,particle2, distMin, distMax) ) continue;
          // particle1, particle2 is a pair‚
          for(auto& neighbor3 : nodePtr->neighbors){
            auto currentNodeP3 = neighbor3;
            do{
              if ((! node_is_in_searchspace(currentNodeP3,particle2, distMin, distMax)) ||
                  (! node_is_in_searchspace(currentNodeP3,particle1, distMin, distMax)))
              {
                currentNodeP3 = currentNodeP3 -> nextNode(nodePtr->level);
                continue;
              }
              if (! currentNodeP3->isLeaf){
                currentNodeP3 = currentNodeP3->children[0].get();
                continue;
              }
              // currentNode3 is leaf and in searchspace of particle1 and particle2
              // search all particles in currentNode
              auto iP3start = currentNodeP3->upper_bound(particle2);
              for( auto iP3 = iP3start ; iP3 < currentNodeP3 -> nParticles; ++iP3){
                auto particle3 = *(currentNodeP3 -> particlesPtr + iP3) ;
                if( ! particles_at_searched_distance(particle2,particle3, distMin, distMax) ) continue;
                if( ! particles_at_searched_distance(particle1,particle3, distMin, distMax) ) continue;
                // particle1, particle2, particle3 is a regular triangle
                for(auto& neighbor4 : nodePtr->neighbors){
                  auto currentNodeP4 = neighbor4;
                  do{
                    if ((! node_is_in_searchspace(currentNodeP4,particle3, distMin, distMax)) ||
                        (! node_is_in_searchspace(currentNodeP4,particle2, distMin, distMax)) ||
                        (! node_is_in_searchspace(currentNodeP4,particle1, distMin, distMax)))
                    {
                      currentNodeP4 = currentNodeP4 -> nextNode(nodePtr->level);
                      continue;
                    }
                    if (! currentNodeP4->isLeaf){
                      currentNodeP4 = currentNodeP4->children[0].get();
                      continue;
                    }
                    // currentNode4 is leaf and in searchspace of particle1, particle2 and particle3
                    // search all particles in currentNode
                    auto iP4start = currentNodeP4->upper_bound(particle3);
                    for( auto iP4 = iP4start ; iP4 < currentNodeP4 -> nParticles; ++iP4){
                      auto particle4 = *(currentNodeP4 -> particlesPtr + iP4) ;
                      if( ! particles_at_searched_distance(particle3,particle4, distMin, distMax) ) continue;
                      if( ! particles_at_searched_distance(particle2,particle4, distMin, distMax) ) continue;
                      if( ! particles_at_searched_distance(particle1,particle4, distMin, distMax) ) continue;
                      // particle1, particle2, particle3, particle4 is a tetrahedron
                      tetrahedra.push_back({particle1,particle2,particle3,particle4});
                    }
                    currentNodeP4 = currentNodeP4 -> nextNode(nodePtr->level);
                  }while (currentNodeP4->level > nodePtr->level);
                }
              }
              currentNodeP3 = currentNodeP3 -> nextNode(nodePtr->level);
            }while (currentNodeP3->level > nodePtr->level);
          }
        }
        currentNodeP2 = currentNodeP2 -> nextNode(nodePtr->level);
      }while(currentNodeP2->level > nodePtr->level);
    }
  }
  return;
}



size_t TOctree::count_tetrahedra_tree(Node* nodePtr, const double distMin, const double distMax){
  // return all tetrahedra with all p1 in Node
  // searches particle wise all relevant Nodes
  size_t nTetrahedra = 0;
  for(auto& particle1 : nodePtr->particles){
    for (auto& neighbor : nodePtr->neighbors){
      auto currentNodeP2 = neighbor;
      do{
        if(! node_is_in_searchspace(currentNodeP2,particle1, distMin, distMax)){
          currentNodeP2 = currentNodeP2 -> nextNode(nodePtr->level);
          continue;
        }
        if (! currentNodeP2->isLeaf){
          currentNodeP2 = currentNodeP2->children[0].get();
          continue;
        }
        auto iP2start = currentNodeP2->upper_bound(particle1);
        for(size_t iP2=iP2start ; iP2 < currentNodeP2->nParticles; ++iP2 ){
          auto particle2 = *(currentNodeP2 -> particlesPtr + iP2) ;
          if( ! particles_at_searched_distance(particle1,particle2, distMin, distMax) ) continue;
          // particle1, particle2 is a pair‚
          for(auto& neighbor3 : nodePtr->neighbors){
            auto currentNodeP3 = neighbor3;
            do{
              if ((! node_is_in_searchspace(currentNodeP3,particle2, distMin, distMax)) ||
                  (! node_is_in_searchspace(currentNodeP3,particle1, distMin, distMax)))
              {
                currentNodeP3 = currentNodeP3 -> nextNode(nodePtr->level);
                continue;
              }
              if (! currentNodeP3->isLeaf){
                currentNodeP3 = currentNodeP3->children[0].get();
                continue;
              }
              // currentNode3 is leaf and in searchspace of particle1 and particle2
              // search all particles in currentNode
              auto iP3start = currentNodeP3->upper_bound(particle2);
              for( auto iP3 = iP3start ; iP3 < currentNodeP3 -> nParticles; ++iP3){
                auto particle3 = *(currentNodeP3 -> particlesPtr + iP3) ;
                if( ! particles_at_searched_distance(particle2,particle3, distMin, distMax) ) continue;
                if( ! particles_at_searched_distance(particle1,particle3, distMin, distMax) ) continue;
                // particle1, particle2, particle3 is a regular triangle
                for(auto& neighbor4 : nodePtr->neighbors){
                  auto currentNodeP4 = neighbor4;
                  do{
                    if ((! node_is_in_searchspace(currentNodeP4,particle3, distMin, distMax)) ||
                        (! node_is_in_searchspace(currentNodeP4,particle2, distMin, distMax)) ||
                        (! node_is_in_searchspace(currentNodeP4,particle1, distMin, distMax)))
                    {
                      currentNodeP4 = currentNodeP4 -> nextNode(nodePtr->level);
                      continue;
                    }
                    if (! currentNodeP4->isLeaf){
                      currentNodeP4 = currentNodeP4->children[0].get();
                      continue;
                    }
                    // currentNode4 is leaf and in searchspace of particle1, particle2 and particle3
                    // search all particles in currentNode
                    auto iP4start = currentNodeP4->upper_bound(particle3);
                    for( auto iP4 = iP4start ; iP4 < currentNodeP4 -> nParticles; ++iP4){
                      auto particle4 = *(currentNodeP4 -> particlesPtr + iP4) ;
                      if( ! particles_at_searched_distance(particle3,particle4, distMin, distMax) ) continue;
                      if( ! particles_at_searched_distance(particle2,particle4, distMin, distMax) ) continue;
                      if( ! particles_at_searched_distance(particle1,particle4, distMin, distMax) ) continue;
                      // particle1, particle2, particle3, particle4 is a tetrahedron
                      nTetrahedra++;
                    }
                    currentNodeP4 = currentNodeP4 -> nextNode(nodePtr->level);
                  }while (currentNodeP4->level > nodePtr->level);
                }
              }
              currentNodeP3 = currentNodeP3 -> nextNode(nodePtr->level);
            }while (currentNodeP3->level > nodePtr->level);
          }
        }
        currentNodeP2 = currentNodeP2 -> nextNode(nodePtr->level);
      }while(currentNodeP2->level > nodePtr->level);
    }
  }
  return nTetrahedra;
}

void TOctree::direct_search_tetrahedra_randomized(const size_t particle1, const double distMin, const double distMax, const size_t maxReturn, const size_t maxTetrahedra, IndexQuadrupletVec& tetrahedra){
  
  size_t tetrahedraFound = 0;
  auto node = get_starting_node(particle1,distMax);
  auto searchSpace = get_searchspace_large_scale(particle1,node -> neighbors, distMin, distMax);
  auto searchSpaceParticles = get_all_particles(searchSpace);

  std::shuffle(searchSpaceParticles.begin(), searchSpaceParticles.end(), randomNumberGenerator);
  size_t maxP2Tetrahedra = std::max(1, (int)(3 * (double) maxReturn / (double) searchSpaceParticles.size() ) ) ;
  for(auto& particle2: searchSpaceParticles){
    if (particle2 <= particle1 ) continue;
    if (! particles_at_searched_distance(particle1,particle2, distMin,distMax) ) continue;
    // found particle2 candidate
    size_t p2Tetrahedra = 0;
    auto searchSpace3 = get_searchspace_triplet(particle1,particle2,searchSpace,distMin,distMax);
    auto searchSpaceParticles3 = get_all_particles(searchSpace3);
    std::shuffle(searchSpaceParticles3.begin(), searchSpaceParticles3.end(), randomNumberGenerator);
    for(auto& particle3: searchSpaceParticles3){
      if (particle3 <= particle2 ) continue;
      if (! particles_at_searched_distance(particle2,particle3, distMin,distMax) ) continue;
      if (! particles_at_searched_distance(particle1,particle3, distMin,distMax) ) continue;
      // found particle3 candidate
      auto searchSpace4 = get_searchspace_quadruplet(particle1,particle2,particle3,searchSpace,distMin,distMax);
      auto searchSpaceParticles4 = get_all_particles(searchSpace4);
      std::shuffle(searchSpaceParticles4.begin(), searchSpaceParticles4.end(), randomNumberGenerator);
      for(auto& particle4: searchSpaceParticles4){
        if (particle4 <= particle3 ) continue;
        if (! particles_at_searched_distance(particle3,particle4, distMin,distMax) ) continue;
        if (! particles_at_searched_distance(particle2,particle4, distMin,distMax) ) continue;
        if (! particles_at_searched_distance(particle1,particle4, distMin,distMax) ) continue;
        // found particle4 candidate
        if(tetrahedra.size() == maxTetrahedra) return;
        tetrahedra.push_back({particle1,particle2,particle3,particle4});
        tetrahedraFound +=1 ;
        p2Tetrahedra +=1 ;
        if (tetrahedraFound >= maxReturn) return;
        if (p2Tetrahedra >= maxP2Tetrahedra) break;
      }
      if (p2Tetrahedra >= maxP2Tetrahedra) break;
    } 
  }
  return;
}

void TOctree::tree_search_tetrahedra_randomized(const size_t particle1, const double distMin, const double distMax, const size_t maxReturn, const size_t maxTetrahedra, IndexQuadrupletVec& tetrahedra){

  size_t tetrahedraFound = 0;
  auto node = get_starting_node(particle1,distMax);
  auto searchSpaceL = get_searchspace_large_scale(particle1,node->neighbors, distMin, distMax);
  auto searchSpace = get_searchspace_pair(particle1,searchSpaceL,distMin,distMax);
  auto searchSpaceParticles = get_all_particles(searchSpace);

  std::shuffle(searchSpaceParticles.begin(), searchSpaceParticles.end(), randomNumberGenerator);
  size_t maxP2Tetrahedra = std::max(1, (int)(3 * (double) maxReturn / (double) searchSpaceParticles.size() ) ) ;
  for(auto& particle2: searchSpaceParticles){
    if (particle2 <= particle1 ) continue;
    if (! particles_at_searched_distance(particle1,particle2, distMin,distMax) ) continue;
    // found particle2 candidate
    size_t p2Tetrahedra = 0;
    auto searchSpace3 = get_searchspace_triplet(particle1,particle2,searchSpaceL,distMin,distMax);
    auto searchSpaceParticles3 = get_all_particles(searchSpace3);
    std::shuffle(searchSpaceParticles3.begin(), searchSpaceParticles3.end(), randomNumberGenerator);
    for(auto& particle3: searchSpaceParticles3){
      if (particle3 <= particle2 ) continue;
      if (! particles_at_searched_distance(particle2,particle3, distMin,distMax) ) continue;
      if (! particles_at_searched_distance(particle1,particle3, distMin,distMax) ) continue;
      // found particle3 candidate
      auto searchSpace4 = get_searchspace_quadruplet(particle1,particle2,particle3,searchSpaceL,distMin,distMax);
      auto searchSpaceParticles4 = get_all_particles(searchSpace4);
      std::shuffle(searchSpaceParticles4.begin(), searchSpaceParticles4.end(), randomNumberGenerator);
      for(auto& particle4: searchSpaceParticles4){
        if (particle4 <= particle3 ) continue;
        if (! particles_at_searched_distance(particle3,particle4, distMin,distMax) ) continue;
        if (! particles_at_searched_distance(particle2,particle4, distMin,distMax) ) continue;
        if (! particles_at_searched_distance(particle1,particle4, distMin,distMax) ) continue;
        // found particle4 candidate
        if(tetrahedra.size() == maxTetrahedra) return;
        tetrahedra.push_back({particle1,particle2,particle3,particle4});
        tetrahedraFound +=1 ;
        p2Tetrahedra +=1 ;
        if (tetrahedraFound >= maxReturn) return;
        if (p2Tetrahedra >= maxP2Tetrahedra) break;
      }
      if (p2Tetrahedra >= maxP2Tetrahedra) break;
    } 
  }
  return;
}

}
