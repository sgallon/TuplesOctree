#ifndef TUPLE_HPP
#define TUPLE_HPP

#include "point.hpp"

#include <array>
#include <vector>
#include <string>
#include <fstream>



namespace octree{
template<size_t N> using IndexNTuple = std::array<size_t,N>;
template<size_t N> using PointNTuple = std::array<Point,N>;
template<size_t N> using IndexNTupleVec = std::vector<IndexNTuple<N> >;
template<size_t N> using PointNTupleVec = std::vector<PointNTuple<N> >;

using IndexPair  = IndexNTuple<2>;
using PointPair  = PointNTuple<2>;
using IndexPairVec   = std::vector<IndexPair>;
using PointPairVec   = std::vector<PointPair>;

using IndexTriplet = IndexNTuple<3>;
using PointTriplet = PointNTuple<3>;
using IndexTripletVec = std::vector<IndexTriplet>;
using PointTripletVec = std::vector<PointTriplet>;

using IndexQuadruplet  = IndexNTuple<4>;
using PointQuadruplet  = PointNTuple<4>;
using IndexQuadrupletVec = std::vector<IndexQuadruplet>;
using PointQuadrupletVec = std::vector<PointQuadruplet>;

template <size_t N> bool writeTuples(std::string filename, std::vector<IndexNTuple<N>> tuples){
auto size = tuples.size();
std::fstream file;
file.open (filename, std::fstream::binary | std::fstream::out);
if (!file){
  std::cerr << "could not open file \""<< filename <<"\" to write" << std::endl;
  throw std::runtime_error("FileWriterError");
}
else{
  for (size_t i = 0; i < size; ++i){
    for(size_t j = 0; j < N ; ++j){
      size_t x = tuples[i][j];
      file.write(reinterpret_cast<char*>(&x),sizeof( size_t ));

    }
  }
}
return true;
};
template <size_t N> std::vector<IndexNTuple<N>> readTuples(std::string filename){
  std::vector<IndexNTuple<N>> tuples;
  std::fstream file(filename,std::fstream::binary | std::fstream::in);
  if (!file){
    std::cerr << "could not read file \""<< filename <<"\"" << std::endl;
    throw std::runtime_error("FileReaderError");
  }
  else{
    while(! file.eof() ){
      IndexNTuple<N> tuple ;
      for(size_t j = 0; j < N ; ++j){
        size_t idx ;
        file.read(reinterpret_cast<char*>(&idx),sizeof(size_t));
        if(file.eof()){
          return tuples;
        }
        else{
          tuple[j] = idx;
        }
      }
      tuples.push_back(tuple);
    }
    return tuples;
  }
};
}
#endif //TUPLE_HPP
