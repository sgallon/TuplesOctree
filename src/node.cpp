#include "node.hpp"

namespace octree{

Node::Node(const Point& center_,const int level_ , const double sideLength_, bool isLeaf_,
    const size_t octant_,  Node* parent_, const IndexVec particles_)
: center(center_),
  level(level_),
  sideLength(sideLength_),
  halfSideLength(0.5 * sideLength_),
  isLeaf(isLeaf_),
  octant(octant_),
  parent(parent_),
  particles(particles_),
  nParticles(particles.size()),
  particleIdxMin((particles.size() == 0 ? 0 :*min_element(particles.begin(),particles.end() ) )),
  particleIdxMax((particles.size() == 0 ? 0 :*max_element(particles.begin(),particles.end() ) )),
  particlesPtr(particles.data())
{

}

void Node::store_neighbors(){
  neighbors.reserve(27);
  for(size_t i=0; i < 27; i++){
    // loop over all directions and store neighbors
    auto direction = Point( -1 + int(i < 9 ? 0 : i/9) ,-1 + fmod(i/3,3),-1 + fmod(i,3));
    auto temp = get_neighbor_node(this , direction);
    if (temp){
      neighbors.push_back(temp);
      neighborArr[i] = temp;
    }
    else{
      neighborArr[i] = nullptr;
    }
  }
}


Node* Node::neighbor(const Point& direction){
   size_t directionIdx = (direction.X < 0 ? 0 : 9 * (1 + direction.X)) + (direction.Y < 0 ? 0 : 3 * (1 + direction.Y)) + (direction.Z < 0 ? 0 : 1 * (1 + direction.Z));
  return neighborArr[directionIdx];
}

Node* get_neighbor_node(Node* nodePtr, const Point& direction){
  /*
  The neighbor finding algorithm is based on the algorithm described in :

  Kunio Aizawa, Koyo Motomura, Shintaro Kimura, Ryosuke Kadowaki and Jia Fan
  "Constant time neighbor finding in quadtrees: An experimental result,"
  2008 3rd International Symposium on Communications, Control and Signal Processing,
  2008, pp. 505-510, doi: 10.1109/ISCCSP.2008.4537278.

  It was implemented here as a recursive algorithm.
  It returns the same-level nodes, if they exist given any direction D
  fulfilling D.X,D.Y and D.Z in {-1,0,1}.


  In each step, the direction are treated separately and iteratively

  recall the labeling of nodes:
  // xz plane low y        // xz plane high y
  // (front)               // (back)
  // 0  1                  // 4  5
  // 2  3                  // 6  7

  Example:
  Find a neighbor in direction (1,1,0) for quadrant 2

   1. find the quadrant number of the neighbor in x-direction
        x-neighbor of 2 is 3
   2. Is the neighbor in same parent node?
        neighbor of 2 in pos. x direction is in same parent node as 2
   3. find the quadrant number of the neighbor in y-direction
        y-neighbor of 3 is 7
   4. Is the neighbor in same parent node?
        neighbor of 3 in pos. y direction is in same  parent node as 3
   5. find the quadrant number of the neighbor in z-direction
        We are not searching in z direction

  return sibling with octant number 6

  Example:
  Find a neighbor in direction (-1,1,-1) for quadrant 4

  1. find the quadrant number of the neighbor in x-direction
        x-neighbor of 4 is 5
  2. Is the neighbor in same parent node?
        neighbor of 4 in neg. x direction is not in same Node as 5
        remember for later: we need to search in neg. x direction of parent
  3. find the quadrant number of the neighbor in y-direction
        y-neighbor of 5 is 1
  4. Is the neighbor in same parent node?
        neighbor of 5 in pos. y direction is not in same Node:
        remember for later: we need to search in pos. y direction of parent
  5. find the quadrant number of the neighbor in z-direction
        z-neighbor of 1 is 3
  4. Is the neighbor in same parent node?
        neighbor of 1 in negative z direction is in same parent node as 1

  Find Neighbor of parent-node in direction (-1,1,0)
  return its child with octant number 3

  */

  if (nodePtr->parent){
    // if parent node exists, nodePtr is not the rootNode

    if (direction.X == 0 && direction.Y == 0 && direction.Z == 0){
      // neighbor in direction 0 is node itself
      return nodePtr;
    }
    auto nextDirection = Point(0.,0.,0.);
    auto nextOctant = nodePtr->octant;
    if (std::abs(direction.X)== 1.){
      switch(nextOctant){
        case 0 : case 2: case 4: case 6: {
          nextDirection.X = (direction.X < 0 ? -1. : 0 );
          break;
        }
        case 1 : case 3: case 5: case 7: {
          nextDirection.X = (direction.X > 0 ? +1. : 0 );
          break;
        }
      }
      // nextOctant_old 0 1 2 3 4 5 6 7
      // nextOctant_new 1 0 3 2 5 4 7 6
      nextOctant = std::fmod(8 + nextOctant + 1 - 2 * std::fmod(nextOctant, 2)  ,8);
    }
    if (std::abs(direction.Y)== 1.){
      switch(nextOctant){
        case 0 : case 1: case 2: case 3: {
          nextDirection.Y = (direction.Y < 0 ? -1. : 0 );
          break;
        }
        case 4 : case 5: case 6: case 7: {
          nextDirection.Y = (direction.Y > 0 ? +1. : 0 );
          break;
        }
      }
      // nextOctant_old 0 1 2 3 4 5 6 7
      // nextOctant_new 4 5 6 7 0 1 2 3
      nextOctant = std::fmod(8 + 4 + nextOctant ,8);
    }
    if (std::abs(direction.Z)== 1.){
      switch(nextOctant){
        case 2 : case 3: case 6: case 7: {
          nextDirection.Z = (direction.Z < 0 ? -1. : 0 );
          break;
        }
        case 0 : case 1: case 4: case 5: {
          nextDirection.Z = (direction.Z > 0 ? +1. : 0 );
          break;
        }
      }
      // nextOctant_old 0 1 2 3 4 5 6 7
      // nextOctant_new 2 3 0 1 6 7 4 5
      nextOctant = std::fmod(8 + nextOctant   + 2 - 4 * std::fmod((nextOctant / 2), 2), 8);
    }
    if (nextDirection.X == 0 && nextDirection.Y == 0 && nextDirection.Z == 0){
      // nextDirection solely 0 <=> searched neighbor is sibling
      return nodePtr ->parent -> children[nextOctant].get();
    }
    else{
      // searching parent neighbor in direction nextDirection
      // if parent neighbor exists:
      // return its child with octant number nextOctant
      auto temp = nodePtr->parent->neighbor(nextDirection);
      return (temp ? temp -> children[nextOctant].get() : nullptr);
    }
  }
  else{
    // nodePtr is root, neighbors therefore do not exist
    return nullptr;
  }
}

IndexVec Node::get_particles(){
  return particles;
}

bool Node::pos_is_in_node(Point pos){
  auto temp = elementwise_absolute(pos - center) / (sideLength/2);
  return ( (temp.X<= 1) && (temp.Y <= 1) && (temp.Z<= 1) );
}

bool Node::pos_is_in_node(const Point* pos){
  return ( (std::abs(pos->X - center.X)<= 0.5*sideLength) &&
           (std::abs(pos->Y - center.Y)<= 0.5*sideLength) &&
           (std::abs(pos->Z - center.Z)<= 0.5*sideLength) ) ;
}

std::ostream& operator<<(std::ostream& os,const Node* nodePtr ){
  return os << "[Node with center " << nodePtr -> center << " and sideLength " << nodePtr -> sideLength << "]" ;
}

std::ostream& operator<<(std::ostream& os,const Node& node ){
  return os << "[Node with center " << node.center << " and sideLength " << node.sideLength << "]";
}

size_t Node::upper_bound(const size_t idx){
  if(idx < particleIdxMin) return 0;
  size_t count = nParticles;
  size_t step = 0;
  size_t it = 0;
  size_t first = 0;
  
  while (count > 0)
  {
    it = first;
    step = count * 0.5;
    it += step;
    if (!( idx < *(particlesPtr + it) )){
      first = ++it;
      count -= step + 1;
      continue;
    }
    count = step;
  }
  return first; 
}
Node* Node::nextNode(int minLevel){
  if (level == minLevel) return this;
  if(octant < 7) return nextSibling;
  return parent->nextNode(minLevel);
}


IndexVec get_all_particles(NodePtrVec nodes){
  size_t nParticlesInNodes = 0 ;
  for(auto& node: nodes){
    nParticlesInNodes += node -> nParticles;
  }

  IndexVec allParticles;
  allParticles.reserve(nParticlesInNodes);
  for(auto& node: nodes){
    allParticles.insert(allParticles.end(), node->particles.begin(), node->particles.end() );
  }
  return allParticles;
}


} // namespace
