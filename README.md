# TuplesOctree
TuplesOctree is a static C++ library for the efficient search of regular tuples (ensembles of points that are pairwise equidistant) using an octree data structure.
Implemented are algorithms for the search of
1. Pairs
2. Pairs with an angle to the z axis
3. Triangles
4. Tetrahedra 
## Literature
For more information see:

Sebastian Gallon. "Multiparticle Dispersion and Irreversibility in Geophysical Turbulent Flows" PhD thesis. Ecole normale supérieure de Lyon, 2024. 
https://theses.hal.science/tel-04963641

## Getting started
### Build TuplesOctree
Clone the repository into a local folder of your choice, here we use /exemplary/path/to/TuplesOctree
```console
git clone https://github.com/sgallon/TuplesOctree.git /exemplary/path/to/TuplesOctree
```
run the following commands from the terminal
```console
cd /exemplary/path/to/TuplesOctree
cmake -B build/
cd build
make 
make install
```
In the last step, the library is installed to the directory:
```console
/exemplary/path/to/TuplesOctree/lib
```
A different install location can be set in CMakeLists.txt
### Include TuplesOctree in a CMake Project
Lets say we want to use the installed library in a CMakeProject called MyCoolProject that builds an executable called 
```console
fancyProgram.out
```
Simply add the following lines to the '''CMakeLists.txt''' of MyCoolProject:
```
target_link_directories(fancyProgram.out PRIVATE  /path/to/TuplesOctree/lib)
target_include_directories(fancyProgram.out PRIVATE /path/to/TuplesOctree/include)
target_link_libraries(fancyProgram.out libTuplesTOctree.a)
```
If MyCoolProject already uses some other libraries, say libExample1.a and libExample2.a, the last line might look like this:
```
target_link_libraries(fancyProgram.out libExample1.a libTuplesTOctree.a libExample2.a)
```
## Examples
The test executables can be seen as examples to use the code.
## Acknowledgements
I am grateful for the many discussions on this code I have had with my PhD advisor Alain Pumir.