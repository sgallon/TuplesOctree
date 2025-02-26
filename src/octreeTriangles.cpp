#include "octree.hpp"
namespace octree
{
  /* ##############################################################################
                                searching
                            public functions
  ############################################################################## */
IndexTripletVec TOctree::search_triangles(const double distance,const double tolerance){
  // returns all triangles
  // searching the whole tree
  // To avoid multiple counting, we only return triangles with the ordering of indices:
  // particle1 < particle2 < particle3 < particle4
  auto triangles = IndexTripletVec();
  const double dist_min = distance * (1. - tolerance);
  const double dist_max = distance * (1. + tolerance);
  size_t level = std::min(depth,(int)trunc(log(boxSize/dist_max)/log(2)) );
  if ( dist_max < boxSize * pow(2,-depth) ){
    // for small distances, we search only in the leaf nodes
    for(auto& node : nodesAtLevel[level]){
      direct_search_triangles(node, dist_min, dist_max, triangles);
    }
    return triangles;
  }
  for(auto& node : nodesAtLevel[level]){
    tree_search_triangles(node,dist_min,dist_max, triangles);
  }  
  return triangles;
}

IndexTripletVec TOctree::search_triangles(const double distance,const double tolerance, const size_t maxTriangles){
  // returns at most maxTriangles triangles
  // searching the whole tree
  // To avoid multiple counting, we only return triangles with the ordering of indices:
  // particle1 < particle2 < particle3
  double nTrianglesEstimation = triangles_estimation(distance, tolerance);

  if (5 * nTrianglesEstimation < maxTriangles ){
    // It is not likely, that we will find more triangles than maxTriangles, we thus can search systematically and not randomized at all
    return search_triangles(distance,tolerance);
  }
  // it is not-unlikely that we find more triangles than desired
  // the algorithm needs to be randomized to obtain a stochastically representative subset
  // instead of systematically searching the Nodes by their position (smaller computational costs, as eg. neighborhood query needs to be called once per Node and not once per particle), we search systematically by varying particle1

  IndexTripletVec triangles;
  triangles.reserve( std::min(maxTriangles, (size_t)(nTrianglesEstimation)) );

  double dist_min = distance * (1. - tolerance);
  double dist_max = distance * (1. + tolerance);

  IndexVec particle1Vec(nParticles); //vector of particle1
  std::iota(particle1Vec.begin(),particle1Vec.end(),0); //fill up particle1Vec with ascending numbers
  std::shuffle(particle1Vec.begin(),particle1Vec.end(),randomNumberGenerator);

  // This is not necessary if we can assume that particles are sufficiently mixed and that the position of a particle is statistically independent of its index

  size_t trianglesPerParticle = 0;
  if (  0.75* nTrianglesEstimation < maxTriangles) {
    // We are expecting few triangles, hence triangles are "rare"
    // Therefore, we need to search in a semi-random way: returning all triangles given a particle1
    // stopping only if we have found enough triangles or if we have searched all particles
    trianglesPerParticle =  maxTriangles ;
  } 
  else if( 0.5 > maxTriangles/(double) nParticles){
    trianglesPerParticle = std::max((size_t) 1, (size_t) (1.5 * round( (double) maxTriangles / (double) nParticles ) )) ;
  } else{
    // triangles are not "rare"
    // Therefore, we need to search in a more randomized way
    // return only a few triangles per particle1, in a randomized way
    trianglesPerParticle = std::max(1 , (int)(2* (double) nTrianglesEstimation / (double) nParticles ) ) ;
    // return slightly more trianglesPerParticle than expected, to avoid finding to few, as for high indices, there are fewer triangles with strictly increasing indices
    // if maxTriangles < nParticles, return 1
  }
  if ( distance < boxSize * pow(2,-depth) ) {
    for(auto& particle1 : particle1Vec){
      if(triangles.size() == maxTriangles) break;
      direct_search_triangles_randomized(particle1,dist_min,dist_max,trianglesPerParticle, maxTriangles, triangles);
    }
    return triangles;
  }
  for(auto& particle1 : particle1Vec){
    if(triangles.size() == maxTriangles) break;
    tree_search_triangles_randomized(particle1,dist_min,dist_max,trianglesPerParticle,maxTriangles, triangles);
  }

  return triangles;
  }

size_t TOctree::count_triangles(const double distance,const double tolerance){
  size_t nTriangles = 0;
  const double dist_min = distance * (1. - tolerance);
  const double dist_max = distance * (1. + tolerance);
  size_t level = std::min(depth,(int)trunc(log(boxSize/dist_max)/log(2)) );
  if ( dist_max < boxSize * pow(2,-depth) ){
    // for small distances, we search only in the leaf nodes
    for(auto& node : nodesAtLevel[level]){
      nTriangles += count_triangles_direct(node, dist_min, dist_max);
    }
    return nTriangles;
  }
  for(auto& node : nodesAtLevel[level]){
    nTriangles += count_triangles_tree(node,dist_min,dist_max);
  }  
  return nTriangles;
}
/* ##############################################################################
                              searching
                          private functions
############################################################################## */
void TOctree::direct_search_triangles(Node* nodePtr, const double distMin, const double distMax, IndexTripletVec& triangles){
  // return all triangles with all p1 in Node
  // searches particle wise all relevant Nodes
  for(auto& particle1 : nodePtr->particles){
    for (auto& neighbor : nodePtr -> neighbors){
      if(! node_is_in_searchspace(neighbor, particle1, distMin, distMax) )  continue;
      auto iP2start = neighbor->upper_bound(particle1);
      for(size_t iP2=iP2start ; iP2 < neighbor->nParticles; ++iP2 ){
        auto particle2 = *(neighbor -> particlesPtr + iP2) ;
        if( ! particles_at_searched_distance(particle1,particle2, distMin, distMax) ) continue;
        // particle1, particle2 is a pair‚
        // search for particle3
        for(auto& neighbor3 : nodePtr -> neighbors){
          if(! node_is_in_searchspace(neighbor3, particle2, distMin, distMax) )  continue;
          if(! node_is_in_searchspace(neighbor3, particle1, distMin, distMax) )  continue;
          auto iP3start = neighbor3->upper_bound(particle2);
          for( auto iP3 = iP3start ; iP3 < neighbor3 -> nParticles; ++iP3){
            auto particle3 = *(neighbor3 -> particlesPtr + iP3) ;
            if( ! particles_at_searched_distance(particle2,particle3, distMin, distMax) ) continue;
            if( ! particles_at_searched_distance(particle1,particle3, distMin, distMax) ) continue;
            // particle1, particle2, particle3 is a regular triangle
            triangles.push_back({particle1,particle2,particle3});
          }
        }
      }
    }
  }
  return ;
}

size_t TOctree::count_triangles_direct(Node* nodePtr, const double distMin, const double distMax){
  // return all triangles with all p1 in Node
  // searches particle wise all relevant Nodes
  size_t nTriangles = 0;
  for(auto& particle1 : nodePtr->particles){
    for (auto& neighbor : nodePtr -> neighbors){
      if(! node_is_in_searchspace(neighbor, particle1, distMin, distMax) )  continue;
      auto iP2start = neighbor->upper_bound(particle1);
      for(size_t iP2=iP2start ; iP2 < neighbor->nParticles; ++iP2 ){
        auto particle2 = *(neighbor -> particlesPtr + iP2) ;
        if( ! particles_at_searched_distance(particle1,particle2, distMin, distMax) ) continue;
        // particle1, particle2 is a pair‚
        // search for particle3
        for(auto& neighbor3 : nodePtr -> neighbors){
          if(! node_is_in_searchspace(neighbor3, particle2, distMin, distMax) )  continue;
          if(! node_is_in_searchspace(neighbor3, particle1, distMin, distMax) )  continue;
          auto iP3start = neighbor3->upper_bound(particle2);
          for( auto iP3 = iP3start ; iP3 < neighbor3 -> nParticles; ++iP3){
            auto particle3 = *(neighbor3 -> particlesPtr + iP3) ;
            if( ! particles_at_searched_distance(particle2,particle3, distMin, distMax) ) continue;
            if( ! particles_at_searched_distance(particle1,particle3, distMin, distMax) ) continue;
            // particle1, particle2, particle3 is a regular triangle
            nTriangles++;
          }
        }
      }
    }
  }
  return nTriangles;
}

void TOctree::tree_search_triangles(Node* nodePtr, const double distMin, const double distMax, IndexTripletVec& triangles){
  // return all triangles with all p1 in Node
  // searches particle wise all relevant Nodes
  for(auto& particle1 : nodePtr->particles){
    for (auto& neighbor : nodePtr -> neighbors){
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
        // currentNode is leaf and in searchspace
        // search all particles in currentNode
        auto iP2start = currentNodeP2->upper_bound(particle1);
        for(size_t iP2=iP2start ; iP2 < currentNodeP2->nParticles; ++iP2 ){
          auto particle2 = *(currentNodeP2 -> particlesPtr + iP2) ;
          if( ! particles_at_searched_distance(particle1,particle2, distMin, distMax) ) continue;
          // particle1, particle2 is a pair‚
          // search for particle3
          for(auto& neighbor3 : nodePtr -> neighbors){
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
                triangles.push_back({particle1,particle2,particle3});
              }
              currentNodeP3 = currentNodeP3 -> nextNode(nodePtr->level);
            }while (currentNodeP3->level > nodePtr->level);
          }
        }
        currentNodeP2 = currentNodeP2 -> nextNode(nodePtr->level);
      }while(currentNodeP2->level > nodePtr->level);
    }
  }
  return ;
}


size_t TOctree::count_triangles_tree(Node* nodePtr, const double distMin, const double distMax){
  size_t nTriangles = 0;
  // return all triangles with all p1 in Node
  // searches particle wise all relevant Nodes
  for(auto& particle1 : nodePtr->particles){
    for (auto& neighbor : nodePtr -> neighbors){
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
        // currentNode is leaf and in searchspace
        // search all particles in currentNode
        auto iP2start = currentNodeP2->upper_bound(particle1);
        for(size_t iP2=iP2start ; iP2 < currentNodeP2->nParticles; ++iP2 ){
          auto particle2 = *(currentNodeP2 -> particlesPtr + iP2) ;
          if( ! particles_at_searched_distance(particle1,particle2, distMin, distMax) ) continue;
          // particle1, particle2 is a pair‚
          // search for particle3
          for(auto& neighbor3 : nodePtr -> neighbors){
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
                nTriangles++; 
              }
              currentNodeP3 = currentNodeP3 -> nextNode(nodePtr->level);
            }while (currentNodeP3->level > nodePtr->level);
          }
        }
        currentNodeP2 = currentNodeP2 -> nextNode(nodePtr->level);
      }while(currentNodeP2->level > nodePtr->level);
    }
  }
  return nTriangles;
}

void TOctree::direct_search_triangles_randomized(const size_t particle1, const double distMin, const double distMax, const size_t maxReturn, const size_t maxTriangles, IndexTripletVec& triangles){
  size_t trianglesFound = 0;
  auto node = get_starting_node(particle1,distMax);
  
  auto searchSpace = get_searchspace_large_scale(particle1,node -> neighbors, distMin, distMax);
  auto searchSpaceParticles = get_all_particles(searchSpace);

  std::shuffle(searchSpaceParticles.begin(), searchSpaceParticles.end(), randomNumberGenerator);
  size_t maxP2Triangles = std::max(1, (int)(2 * (double) maxReturn / (double) searchSpaceParticles.size() ) ) ;


  for(auto& particle2: searchSpaceParticles){
    if (particle2 <= particle1 ) continue;
    if (! particles_at_searched_distance(particle1,particle2, distMin,distMax) ) continue;
    // found particle2 candidate
    size_t p2Triangles = 0;
    auto searchSpace3 = get_searchspace_triplet(particle1,particle2,searchSpace,distMin,distMax);
    auto searchSpaceParticles3 = get_all_particles(searchSpace3);
    std::shuffle(searchSpaceParticles3.begin(), searchSpaceParticles3.end(), randomNumberGenerator);
    
    for(auto& particle3: searchSpaceParticles3){
      if (particle3 <= particle2 ) continue;
      if (! particles_at_searched_distance(particle2,particle3, distMin,distMax) ) continue;
      if (! particles_at_searched_distance(particle1,particle3, distMin,distMax) ) continue;
      // found particle3 candidate
      triangles.push_back({particle1,particle2,particle3});
      trianglesFound +=1 ;
      p2Triangles += 1;
      if ((trianglesFound == maxReturn) || (triangles.size() == maxTriangles )){
        return;
      }
      if (p2Triangles == maxP2Triangles){
        break;
      }
    } 
  }

  
  return ;
}

void TOctree::tree_search_triangles_randomized(const size_t particle1, const double distMin, const double distMax, const size_t maxReturn, const size_t maxTriangles, IndexTripletVec& triangles){
  size_t trianglesFound = 0;
  auto node = get_starting_node(particle1,distMax);
  auto searchSpaceL = get_searchspace_large_scale(particle1,node -> neighbors, distMin, distMax);

  auto searchSpaceParticle2 = get_searchspace_pair(particle1, searchSpaceL, distMin,distMax) ;
  auto searchSpaceParticles2 = get_all_particles(searchSpaceParticle2);
  std::shuffle(searchSpaceParticles2.begin(),searchSpaceParticles2.end(),randomNumberGenerator);
  size_t maxP2Triangles = std::max(1, (int)(2 * (double) maxReturn / (double) searchSpaceParticles2.size() ) ) ;

  for(auto& particle2: searchSpaceParticles2){
    if (particle2 <= particle1 ) continue;
    if (! particles_at_searched_distance(particle1,particle2, distMin,distMax) ) continue;
    auto searchSpaceParticle3 = get_searchspace_triplet(particle1,particle2,searchSpaceL,distMin,distMax);
    auto searchSpaceParticles3 = get_all_particles(searchSpaceParticle3);
    size_t p2Triangles = 0;
    std::shuffle(searchSpaceParticles3.begin(),searchSpaceParticles3.end(),randomNumberGenerator);
    for(auto& particle3: searchSpaceParticles3){
      if (particle3 <= particle2 ) continue;
      if (! particles_at_searched_distance(particle2,particle3, distMin,distMax) ) continue;
      if (! particles_at_searched_distance(particle1,particle3, distMin,distMax) ) continue;
      triangles.push_back({particle1,particle2,particle3});
      trianglesFound +=1 ;
      p2Triangles += 1;
      if ((trianglesFound == maxReturn) || (triangles.size() == maxTriangles )){
        return;
      }
      if (p2Triangles == maxP2Triangles){
        break;
      }
    }
  }
return;
}
}
