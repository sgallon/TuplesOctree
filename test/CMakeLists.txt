set(test_executables
    testpairs.out
    testorientedpairs.out
    testtriangles.out
    testtetrahedra.out
    )


set(testpairs.out_SOURCES 
    pairs.cpp 
    )

set(testorientedpairs.out_SOURCES
    orientedpairs.cpp
    )

set(testtriangles.out_SOURCES
    triangles.cpp
    )

set(testtetrahedra.out_SOURCES
    tetrahedra.cpp
    )


foreach(test_executable ${test_executables})
    add_executable(${test_executable} ${${test_executable}_SOURCES})
    target_link_libraries(${test_executable} ${PROJECT_NAME} )
    add_test(NAME ${test_executable} COMMAND ${test_executable} )   
endforeach()