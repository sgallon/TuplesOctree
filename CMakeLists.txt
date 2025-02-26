cmake_minimum_required(VERSION 3.18)

project(TuplesTOctree 
        VERSION 1.0
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)



include_directories(include)

add_library(
    ${PROJECT_NAME} STATIC
    src/node.cpp
    src/point.cpp
    src/octree.cpp
    src/octreeOrientedPairs.cpp
    src/octreePairs.cpp
    src/octreeTriangles.cpp
    src/octreeTetrahedra.cpp
)
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)

target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -Wpedantic)

set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR})
install(TARGETS ${PROJECT_NAME} ARCHIVE DESTINATION lib)

enable_testing()
add_subdirectory(test)