#ifndef NODE_HPP
#define NODE_HPP

#include "point.hpp"

#include <array>
#include <iostream>
#include <math.h>
#include <memory>
#include <vector>
#include <algorithm> //max_element

namespace octree
{
using IndexVec = std::vector< size_t >;

class Node
{
public:
  // Constructor
  Node(const Point& center_,const int level_ , const double sideLength_, bool isLeaf_,
      const size_t octant_,  Node* parent_, const IndexVec particles_);

  IndexVec get_particles();
  void store_neighbors();

  std::vector< Node* > neighbors;
  Node* neighbor(const Point& direction);
  Node* nextNode(int minLevel);
  const Point center;
  const int level;
  const double sideLength;
  const double halfSideLength;
  const bool isLeaf;
  const size_t octant;
  Node* parent;
  Node* nextSibling;
  const IndexVec particles;
  const size_t nParticles;
  const size_t particleIdxMin;
  const size_t particleIdxMax;
  const size_t* particlesPtr;

  std::array< std::unique_ptr<Node> ,8> children;

  bool pos_is_in_node(Point pos);
  bool pos_is_in_node(const Point* pos);

  size_t upper_bound(const size_t idx);

  friend std::ostream& operator<<(std::ostream& os,const Node* nodePtr );
  friend std::ostream& operator<<(std::ostream& os,const Node& node );

private:
  // std::vector< Node* > neighbors;
  std::array< Node* ,27> neighborArr; // only used when a neighbor in a certain direction is requested
};
Node* get_neighbor_node(Node* nodePtr, const Point& direction);
using NodePtrVec = std::vector<Node*>;
IndexVec get_all_particles(NodePtrVec Nodes);

} // namespace


#endif // NODE_HPP
