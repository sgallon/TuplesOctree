#ifndef OCTREE_HPP
#define OCTREE_HPP

#include "point.hpp"
#include "node.hpp"
#include "tuple.hpp"

#include <array>
#include <algorithm>
#include <iostream>
#include <math.h>
#include <memory>
#include <numeric>
#include <random>
#include <vector>
#include <algorithm> //max_element


namespace octree
{
class TOctree {
  private:
    std::unique_ptr<Node> rootNode;
    std::mt19937 randomNumberGenerator;

    void build_tree(Node* nodePtr, const IndexVec& inputParticles);
    void create_children_of_node(Node* nodePtr,const bool areLeaf,std::array< IndexVec,8>& inputParticles);

    void store_neighbors(Node* nodePtr);
    
    // helper functions
    Node * get_search_start_node(const size_t point1, const size_t level);

    inline double get_distance(const Point& point1, const Point& point2);
    inline double get_distance(const size_t particle1, const size_t particle2);

    bool particles_at_searched_distance(const size_t particle1, const size_t particle2, const double dist_min, const double dist_max);
    bool particles_in_z_cone(const size_t particle1, const size_t particle2, const double cosAngle_min, const double cosAngle_max);
    
    size_t get_sub_octant(Node* nodePtr, const size_t particle);
    inline size_t get_sub_octant(const Point& center, const Point querypoint);
    inline size_t get_sub_octant(const Point& center, const size_t particle);
    inline Point get_child_center(Node* nodePtr,const size_t octant);



    double min_dist_to_node(Node* nodePtr, const size_t particle );
    double max_dist_to_node(Node* nodePtr, const size_t particle );
    double min_angle_z_axis(Node* nodePtr, const size_t particle );
    //Searching

    Node* get_starting_node(const size_t particle1, const double dist_max);
    bool node_is_in_searchspace(Node* nodePtr,const size_t particle,const double dist_min,const double dist_max);
    bool node_is_in_cone(Node* nodePtr, const size_t particle, const double cosAngle_min, const double cosAngle_max  );

    NodePtrVec get_searchspace_large_scale(const size_t particle1, const NodePtrVec& neighbors, const double dist_min, const double dist_max);

    NodePtrVec get_searchspace_pair(const size_t particle1, const NodePtrVec& searchSpaceLargeScale,
      const double dist_min, const double dist_max);

    NodePtrVec get_searchspace_oriented_pair(const size_t particle1, const NodePtrVec& searchSpaceLargeScale, const double dist_min,const double dist_max, const double cosAngle_min, const double cosAngle_max);

    NodePtrVec get_searchspace_triplet(const size_t particle1,const size_t particle2 ,
      const NodePtrVec& searchSpaceLargeScale, const double dist_min, const double dist_max);

    NodePtrVec get_searchspace_quadruplet(const size_t particle1, const size_t particle2 , const size_t particle3 , const NodePtrVec& searchSpaceLargeScale, const double dist_min, const double dist_max);

    //implemented in src/octreePairs.cpp
    void direct_search_pairs(Node* nodePtr, const double dist_min,  const double dist_max, IndexPairVec& pairs);
    void direct_search_pairs_randomized(const size_t particle1, const double dist_min, const double dist_max, const size_t maxReturn, const size_t maxPairs, IndexPairVec& pairs);

    void tree_search_pairs(Node* nodePtr, const double distance_min, const double distance_max, IndexPairVec& pairs);
    void tree_search_pairs_randomized(const size_t particle1, const double dist_min, const double dist_max, const size_t maxReturn, const size_t maxPairs , IndexPairVec& pairs);

    size_t count_pairs_direct(Node* nodePtr, const double dist_min,  const double dist_max);
    size_t count_pairs_tree(Node* nodePtr, const double distance_min, const double distance_max);

    //implemented in src/octreeAnisotropicPairs.cpp
    void direct_search_oriented_pairs(Node* nodePtr, const double dist_min, const double dist_max,  const double cosAngle_min, const double cosAngle_max, IndexPairVec& pairs);        
    void tree_search_oriented_pairs(Node* nodePtr, const double dist_min, const double dist_max, const double cosAngle_min, const double cosAngle_max, IndexPairVec& pairs);
    void tree_search_oriented_pairs_randomized(const size_t particle1, const double dist_min, const double dist_max, const double cosAngle_min, const double cosAngle_max, const size_t level, const size_t maxReturn, const size_t maxPairs, IndexPairVec& pairs, size_t& pairsFound);  

    size_t count_oriented_pairs_direct(Node* nodePtr, const double dist_min, const double dist_max,  const double cosAngle_min, const double cosAngle_max);
    size_t count_oriented_pairs_tree(Node* nodePtr, const double dist_min, const double dist_max,  const double cosAngle_min, const double cosAngle_max);
    // implemented in src/octreeTriangles.cpp:
    void direct_search_triangles(Node* nodePtr, const double distMin, const double distMax, IndexTripletVec& triangles);
    void tree_search_triangles(Node* nodePtr, const double distMin, const double distMax, IndexTripletVec& triangles);

    void direct_search_triangles_randomized(const size_t particle1, const double distMin, const double distMax, const size_t maxReturn,const size_t maxTriangles, IndexTripletVec& triangles);
    void tree_search_triangles_randomized(const size_t particle1, const double distMin, const double distMax, const size_t maxReturn, const size_t maxTriangles, IndexTripletVec& triangles);

    size_t count_triangles_direct(Node* nodePtr, const double dist_min,  const double dist_max);
    size_t count_triangles_tree(Node* nodePtr, const double distance_min, const double distance_max);


    //implemented in src/octreeTetrahedra.cpp
    void direct_search_tetrahedra(Node* nodePtr, const double distMin, const double distMax, IndexQuadrupletVec& tetrahedra);
    void tree_search_tetrahedra(Node* nodePtr, const double distMin, const double distMax, IndexQuadrupletVec& tetrahedra);

    void direct_search_tetrahedra_randomized(const size_t particle1, const double distMin, const double distMax, const size_t maxReturn, const size_t maxTetrahedra, IndexQuadrupletVec& tetrahedra);
    void tree_search_tetrahedra_randomized(const size_t particle1, const double distMin, const double distMax, const size_t maxReturn, const size_t maxTetrahedra, IndexQuadrupletVec& tetrahedra);

    size_t count_tetrahedra_direct(Node* nodePtr, const double dist_min,  const double dist_max);
    size_t count_tetrahedra_tree(Node* nodePtr, const double distance_min, const double distance_max);

    const PointVec data;
    const Point * dataPtr;
    std::vector< NodePtrVec > nodesAtLevel;

  public:
    //Constructor
    TOctree();
    TOctree(const PointVec data_,const  Point center_,const double box_size_,const int depth);
    const size_t nParticles;
    const Point center;
    const double boxSize;
    const int depth;

    double tetrahedra_estimation(const  double distance, const double tolerance);
    double pairs_estimation(const  double distance, const double tolerance);
    double oriented_pairs_estimation(const  double distance, const double tolerance, const double angleToZ, const double deltaAngle );
    double triangles_estimation(const  double distance, const double tolerance);
    
    //implemented in src/octreePairs.cpp:
    IndexPairVec search_pairs(const  double distance, const double tolerance);
    IndexPairVec search_pairs(const  double distance, const double tolerance,const size_t maxPairs);
    
    //implemented in src/octreeAnisotropicPairs.cpp:
    IndexPairVec search_oriented_pairs(const  double distance, const double tolerance, const double angleToZ, const double deltaAngle);
    IndexPairVec search_oriented_pairs(const  double distance, const double tolerance, const double angleToZ, const double deltaAngle, const size_t maxPairs);

    // implemented in src/octreeTriangles.cpp:
    IndexTripletVec search_triangles(const double distance, const double tolerance);
    IndexTripletVec search_triangles(const double distance, const double tolerance, const size_t maxTriangles);
    //implemented in src/octreeTetrahedra.cpp:
    IndexQuadrupletVec search_tetrahedra(const double distance,const double tolerance);
    IndexQuadrupletVec search_tetrahedra(const double distance,const double tolerance,const size_t maxTetrahedra);

    size_t count_pairs(const double distance, const double tolerance);
    size_t count_oriented_pairs(const double distance, const double tolerance, const double angleToZ, const double deltaAngle);
    size_t count_triangles(const double distance, const double tolerance);
    size_t count_tetrahedra(const double distance, const double tolerance);

};
} // namespace
#endif //OCTOTREE_HPP
