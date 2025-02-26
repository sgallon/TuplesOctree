#include "octree.hpp"

namespace octree
{
  // ##############################################################################
  // Constructor & building of tree
  // ##############################################################################
TOctree::TOctree(const PointVec dataIn,const Point center_,const double boxSize_,const int depth_)
  : data(dataIn),
    nParticles(dataIn.size()),
    center(center_),
    boxSize(boxSize_),
    depth(depth_)
{
  std::random_device rd;
  randomNumberGenerator  = std::mt19937(rd());

  IndexVec particles(nParticles);
  std::iota(particles.begin(),particles.end(),0);
    
  dataPtr = data.data();
  if (depth == 0){
    //Tree of depth 0, is just the domain. All search algorithms should reduce to direct search
    rootNode = std::make_unique<Node>( center,0, boxSize_,true,0, nullptr,particles);
    return;
  }
  rootNode  = std::make_unique< Node >(center, 0 ,boxSize_,false, 0, nullptr, particles);
  nodesAtLevel.resize(depth+1);
  for(int level=0 ; level < depth+1; ++level){
    nodesAtLevel[level].reserve(pow(8,level));
  }
  nodesAtLevel[0].push_back(rootNode.get());
  build_tree(rootNode.get(), particles); // build the tree recursively, also fills nodesAtLevel

  // store the neighbors of each node

  // loop over all levels
  for (int level = 0; level < depth+1; ++level){
    // loop over all nodes at this level
    for (size_t iNode = 0; iNode < nodesAtLevel[level].size(); ++iNode){
      nodesAtLevel[level][iNode]->store_neighbors();
    }
  }
}

TOctree::TOctree():
  data(PointVec()),
  dataPtr(nullptr),
  nParticles(0),
  center(Point()),
  boxSize(0.),
  depth(0.)
  {};

void TOctree::build_tree(Node* nodePtr,const IndexVec& inputParticles){
  if (!(nodePtr -> level < depth)){
    std::cerr << "Error: node level is larger than tree depth" << std::endl;
    throw std::runtime_error("TreeError");
    return;
  }
  
  std::array<IndexVec,8> children_particles;
  for (size_t i = 0; i < 8 ; i++){
  // reserve rough estimate of particle array size
    children_particles[i].reserve(inputParticles.size()/8);
  }
  // loop over particles to put them into boxes
  for (size_t i = 0; i < inputParticles.size() ; i++){
    children_particles[get_sub_octant(nodePtr -> center , inputParticles[i] )].push_back(inputParticles[i]);
  }

  create_children_of_node(nodePtr,(nodePtr -> level == depth - 1), children_particles);
  if (nodePtr -> level == depth - 1) return;
  for( size_t i = 0; i < 8; ++i ){
    build_tree(nodePtr -> children[i].get(), children_particles[i]);
  }
  return;
}

void TOctree::create_children_of_node(Node* nodePtr, const bool areLeaf,  std::array<IndexVec,8>& inputParticles){
  // xz plane low y        // xz plane high y
  // (front)               // (back)
  // 0  1                  // 4  5
  // 2  3                  // 6  7
  for( size_t i = 0; i < 8; ++i ){
    nodePtr -> children[i] =std::make_unique< Node >( get_child_center(nodePtr , i),nodePtr -> level + 1, nodePtr -> sideLength*0.5,areLeaf, i
    , nodePtr, inputParticles[i]);
  }
  for( size_t i = 0; i < 8; ++i ){
    nodePtr -> children[i] -> nextSibling = nodePtr -> children[(i+1)%8].get();
    nodesAtLevel[nodePtr -> level + 1].push_back(nodePtr -> children[i].get());
  }

  return;
}

 /* ##############################################################################
                               Estimation of tuples
 ############################################################################## */
double TOctree::tetrahedra_estimation(const  double distance, const double tolerance){
  return 128./9. *  M_PI * M_PI *  pow((double)nParticles,4)/pow(boxSize*boxSize*boxSize,3)
    * (3. * tolerance + pow(tolerance,3)) * pow(distance,9) * pow(tolerance,5);
}
double TOctree::pairs_estimation(const  double distance, const double tolerance){
  return 4./3. * M_PI * pow((double)nParticles,2)/(boxSize*boxSize*boxSize) *  pow(distance,3) * (3. * tolerance + pow(tolerance,3));
}

double TOctree::oriented_pairs_estimation(const  double distance, const double tolerance, const double angleToZ, const double deltaAngle ){
  auto angle_max = angleToZ + deltaAngle;
  if (angle_max > M_PI*0.5) {
    angle_max = M_PI*0.5;
  }
  auto angle_min = angleToZ - deltaAngle;
  if (angle_min < 0.) {
    angle_min = 0.;
  }
  return pairs_estimation(distance,tolerance) * (cos(angle_min)-cos(angle_max));
}

double TOctree::triangles_estimation(const  double distance, const double tolerance){
  return 32./9. * M_PI * M_PI * pow((double)nParticles,3)/pow(boxSize*boxSize*boxSize,2) * (3. * tolerance + pow(tolerance,3)) * pow(distance,6) * pow(tolerance,2);
}
/* ##############################################################################
                              searching
                          private functions
############################################################################## */

//Search space functions
bool TOctree::node_is_in_searchspace(Node* nodePtr,const size_t particle, const double dist_min,const double dist_max){
  if (nodePtr -> nParticles ==0) return false;
  if (nodePtr -> particleIdxMax < particle) return false;
  if ( min_dist_to_node(nodePtr,particle) > dist_max ){
    return false;
  }if (max_dist_to_node(nodePtr,particle) <  dist_min){
    return false;
  }
  return true;
}

Node* TOctree::get_starting_node(const size_t particle1, const double dist_max){
  const size_t searchLevel = std::min(depth,(int)trunc(log(boxSize/dist_max)/log(2)) );
  Node* currentNode = rootNode.get();

  for(size_t level = 0; level < searchLevel; level++){
    auto octant = get_sub_octant(currentNode,particle1);
    currentNode = currentNode->children[octant].get();
  }
  return currentNode;
}
bool TOctree::node_is_in_cone(Node* nodePtr, const size_t particle, const double cosAngle_min, const double cosAngle_max ){
  auto pointPtr = dataPtr + particle;
  // calculate the minimal and maximal cos that define the cone
  
  // define variables for the cosine of the minimal angle and the maximal angle between the particle and the node
  double cosMaxAngle = 0.; 
  double cosMinAngle = 0.;

  double rminH2 = 0.;
  double rmaxH2 = 0.;
  double rminV2 = 0.;
  double rmaxV2 = 0.;
  double temp = 0;

  // calculate the maximum angle to z between particle and node
  // calculate squared minimal horizontal distance from particle to node
  temp = std::abs(nodePtr->center.X - pointPtr->X) - nodePtr->halfSideLength;
  rminH2 += (temp > 0 ? temp*temp : 0);
  temp = std::abs(nodePtr->center.Y - pointPtr->Y) - nodePtr->halfSideLength;
  rminH2  += (temp > 0 ? temp*temp : 0);  
  if (rminH2 ==0.){
    cosMinAngle = 1.;
  }else{
    // calculate squared maximal vertical distance  from particle to node
    temp = nodePtr->center.Z + (nodePtr->center.Z - pointPtr->Z < 0 ? -1 : 1 ) * nodePtr->halfSideLength - pointPtr->Z;    
    rmaxV2 = temp*temp;
    cosMinAngle = sqrt(rmaxV2)/sqrt(rmaxV2+rminH2);
  }
  if (cosMinAngle < cosAngle_min){
    // minimal angle is larger than the maximal angle of the cone
    return false;
  }
  // calculate the minimum angle to z between particle and node
  // calculate squared minimal vertical distance from particle to node
  temp = std::abs(nodePtr->center.Z - pointPtr->Z) - nodePtr->halfSideLength;
  rminV2 += (temp > 0 ? temp*temp : 0);
  if (rminV2 ==0.){
    // vertical distance is smaller than halfSideLength
    // maximal angle is 90°
    cosMaxAngle = 0.;
  }else{
    // calculate squared maximal horizontal distance from particle to node
    temp = nodePtr->center.X + (nodePtr->center.X - pointPtr->X < 0 ? -1 : 1 ) * nodePtr->halfSideLength - pointPtr->X;    
    rmaxH2 = temp*temp;
    temp = nodePtr->center.Y + (nodePtr->center.Y - pointPtr->Y < 0 ? -1 : 1 ) * nodePtr->halfSideLength - pointPtr->Y;
    rmaxH2 += temp*temp;

    cosMaxAngle = sqrt(rminV2)/sqrt(rminV2+rmaxH2);
  }
  if (cosMaxAngle > cosAngle_max) {
    // maximal angle is smaller than the minimal angle of the cone
    return false;
  }
  return true;
}

NodePtrVec TOctree::get_searchspace_large_scale(const size_t particle1, const NodePtrVec& neighbors, const double dist_min, const double dist_max){
  auto searchSpaceL = NodePtrVec() ;
  searchSpaceL.reserve(27);
  for(auto& neighbor: neighbors){
    if (node_is_in_searchspace(neighbor, particle1, dist_min, dist_max) ){
      searchSpaceL.push_back(neighbor);
    }
  }
  return searchSpaceL;
}

NodePtrVec TOctree::get_searchspace_pair(const size_t particle1,
  const NodePtrVec& searchSpaceLargeScale, const double dist_min,
  const double dist_max){
  auto searchSpaceS = NodePtrVec() ;
  
  for(auto& candidate : searchSpaceLargeScale){
    auto currentNodeP2 = candidate;
    auto level = currentNodeP2->level;
    do{
      if(! node_is_in_searchspace(currentNodeP2,particle1, dist_min, dist_max)){
        currentNodeP2 = currentNodeP2->nextNode(level);
        continue;
      }
      if (! currentNodeP2->isLeaf){
        currentNodeP2 = currentNodeP2 -> children[0].get();
        continue;
      }
      searchSpaceS.push_back(currentNodeP2);
      currentNodeP2 = currentNodeP2->nextNode(level);
    } while (currentNodeP2 -> level > level);
  }
  return searchSpaceS;
}

NodePtrVec TOctree::get_searchspace_oriented_pair(const size_t particle1,
  const NodePtrVec& searchSpaceLargeScale, const double dist_min,
  const double dist_max, const double cosAngle_min, const double cosAngle_max){


  auto searchSpaceS = NodePtrVec() ;
  for(auto& candidate : searchSpaceLargeScale){
    auto currentNodeP2 = candidate;
    auto level = currentNodeP2->level;
    do{
      if (
        ! node_is_in_searchspace(currentNodeP2, particle1, dist_min, dist_max) ||
        ! node_is_in_cone(currentNodeP2,particle1, cosAngle_min, cosAngle_max ) 
      ){  
        currentNodeP2 = currentNodeP2->nextNode(level);
        continue;
      }
      if (! currentNodeP2->isLeaf){
        currentNodeP2 = currentNodeP2 -> children[0].get();
        continue;
      }
      searchSpaceS.push_back(currentNodeP2);
      currentNodeP2 = currentNodeP2->nextNode(level);
    } while (currentNodeP2 -> level > level);
  }
  return searchSpaceS;
}

NodePtrVec TOctree::get_searchspace_triplet(const size_t particle1,
  const size_t particle2 ,  const NodePtrVec& searchSpaceLargeScale,
  const double dist_min, const double dist_max){
/*
Return search space of a third particle, given two particles
*/
  auto searchSpaceS = NodePtrVec() ;
  for(auto& candidate : searchSpaceLargeScale){
    auto currentNode = candidate;
    auto level = currentNode->level;
    do{
      if(! node_is_in_searchspace(currentNode,particle2, dist_min, dist_max) ||
         ! node_is_in_searchspace(currentNode,particle1, dist_min, dist_max))
      {
        currentNode = currentNode->nextNode(level);
        continue;
      }
      if (! currentNode->isLeaf){
        currentNode = currentNode -> children[0].get();
        continue;
      }
      searchSpaceS.push_back(currentNode);
      currentNode = currentNode->nextNode(level);
    } while (currentNode -> level > level);
  }
  return searchSpaceS;
}

NodePtrVec TOctree::get_searchspace_quadruplet(const size_t particle1, const size_t particle2 , const size_t particle3,
  const NodePtrVec& searchSpaceLargeScale,  const double dist_min, const double dist_max){
  auto searchSpaceS = NodePtrVec() ;
  for(auto& candidate : searchSpaceLargeScale){
    auto currentNode = candidate;
    auto level = currentNode->level;
    do{
      if (! node_is_in_searchspace(currentNode,particle3, dist_min, dist_max) ||
          ! node_is_in_searchspace(currentNode,particle2, dist_min, dist_max) ||
          ! node_is_in_searchspace(currentNode,particle1, dist_min, dist_max) ){
        currentNode = currentNode->nextNode(level);
        continue;
      }
      if (! currentNode->isLeaf){
        currentNode = currentNode -> children[0].get();
        continue;
      }
      searchSpaceS.push_back(currentNode);
      currentNode = currentNode->nextNode(level);
    } while (currentNode -> level > level);
  }
  return searchSpaceS;
}

// ##############################################################################
// Small helper functions
// ##############################################################################

Node* TOctree::get_search_start_node(const size_t particle1, const size_t level){
  Node* currentNodeP1 = rootNode.get();
  size_t currentLevel=0;
  while (currentLevel < level)
  {
    currentNodeP1 = currentNodeP1->children[get_sub_octant(currentNodeP1,particle1)].get();
    currentLevel++;
  }
  
  return  currentNodeP1;
}

size_t TOctree::get_sub_octant(Node* nodePtr, const size_t particle){
  auto pointPtr = dataPtr + particle;
  return ( pointPtr->X < nodePtr->center.X ? 0 : 1 ) + (pointPtr->Y < nodePtr->center.Y ? 0 : 4 ) + (pointPtr->Z < nodePtr->center.Z ? 2 : 0 );
}
inline size_t TOctree::get_sub_octant(const Point& nodeCenter, const Point querypoint){
  return (querypoint.X < nodeCenter.X ? 0 : 1 ) + (querypoint.Y < nodeCenter.Y ? 0 : 4 ) + (querypoint.Z < nodeCenter.Z ? 2 : 0 );
}
inline size_t TOctree::get_sub_octant(const Point& nodeCenter, const size_t particle){
  auto pointPtr = dataPtr + particle;
  return ( pointPtr->X < nodeCenter.X ? 0 : 1 ) + (pointPtr->Y < nodeCenter.Y ? 0 : 4 ) + (pointPtr->Z < nodeCenter.Z ? 2 : 0 );
}

inline Point TOctree::get_child_center(Node* nodePtr, const size_t octant){
  return Point (
    nodePtr->center.X +  ( (octant % 2) == 0 ? -1. : +1. )  * 0.25 * nodePtr -> sideLength,
    nodePtr->center.Y +  ( (octant / 4) == 0 ? -1. : +1. )  * 0.25 * nodePtr -> sideLength,
    nodePtr->center.Z +  ( (octant / 2)%2==0 ? +1. : -1. )  * 0.25 * nodePtr -> sideLength
  );
}
double TOctree::min_dist_to_node(Node* nodePtr, const size_t particle ){
  double dist2 = 0.;
  auto pointPtr = dataPtr+particle;
  double halfSideLength = nodePtr -> sideLength*0.5;
  auto d = std::abs(nodePtr->center.X - pointPtr->X) - halfSideLength;
  dist2 += (d > 0 ? d*d : 0);
  d = std::abs(nodePtr->center.Y - pointPtr->Y) - halfSideLength;
  dist2 += (d > 0 ? d*d : 0);
  d = std::abs(nodePtr->center.Z - pointPtr->Z )- halfSideLength;
  dist2 += (d > 0 ? d*d : 0);
  return sqrt(dist2);
}

double TOctree::max_dist_to_node(Node* nodePtr, const size_t particle){
  auto pointPtr = dataPtr + particle;
  auto halfSideLength = nodePtr -> sideLength*0.5;
  double dist2 = 0.;
  double temp = 0;
  temp = nodePtr->center.X + (nodePtr->center.X - pointPtr->X < 0 ? -1 : 1 ) * halfSideLength - pointPtr->X;
  dist2 += temp*temp;
  temp = nodePtr->center.Y + (nodePtr->center.Y - pointPtr->Y < 0 ? -1 : 1 ) * halfSideLength - pointPtr->Y;
  dist2 += temp*temp;
  temp = nodePtr->center.Z + (nodePtr->center.Z - pointPtr->Z < 0 ? -1 : 1 ) * halfSideLength - pointPtr->Z;
  dist2 += temp*temp;
  return sqrt(dist2);
}

double TOctree::min_angle_z_axis(Node* nodePtr, const size_t particle ){
  auto pointPtr = dataPtr+particle;
  auto halfSideLength = nodePtr -> sideLength * 0.5;
  if ( (std::abs(nodePtr->center.X - pointPtr->X) < halfSideLength) &&
       (std::abs(nodePtr->center.Y - pointPtr->Y) < halfSideLength) )
  {
    return 0;
  }

  auto d = std::abs(nodePtr->center.X - pointPtr->X) - halfSideLength;
  double distPerpendicular = (d > 0 ? d*d : 0.);
  d = std::abs(nodePtr->center.Y - pointPtr->Y) - halfSideLength;
  distPerpendicular += (d > 0 ? d*d : 0.);
  distPerpendicular = sqrt(distPerpendicular);
  auto distParallel = std::abs(center.Z - pointPtr->Z) + halfSideLength;
  return atan2(distPerpendicular, distParallel);
}
inline double TOctree::get_distance(const Point& point1, const Point& point2){
      return norm( (point2 - point1) ) ;
}
inline double TOctree::get_distance(const size_t particle1, const size_t particle2){
  return norm_difference((dataPtr + particle1) , (dataPtr + particle2) );
}

bool TOctree::particles_at_searched_distance(const size_t particle1, const size_t particle2,
  const double distance_min, const double distance_max){
    auto distanceP1P2 = get_distance(particle1, particle2);
    if (distanceP1P2 > distance_max) return false;
    if (distanceP1P2 < distance_min) return false;
    return true;
}
bool TOctree::particles_in_z_cone(const size_t particle1, const size_t particle2, const double cosAngle_min, const double cosAngle_max){
  auto cosAngleToZ = cosAngleToZAxis(dataPtr + particle1, dataPtr + particle2);
  if (cosAngleToZ < cosAngle_min) return false;
  if (cosAngleToZ > cosAngle_max) return false;
  return true;
}
} // namespace
