cmake_minimum_required(VERSION 3.14)

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall -Wno-sign-compare -O3 -g3")

SET(TEST_CLASSES
    testArc/src/testArc.cc
    testGeodesic/src/testGeodesic.cc
    testLocus/src/testLocus.cc
    testLLPoint/src/testLLPoint.cc 
    testShape/src/testShape.cc
    testVector/src/testVector.cc
    )

SET(TEST_SRC
        ${TEST_SRC_DIR}/testGeolib.cc
        ${TEST_SRC_DIR}/testUtil.cc
        ${TEST_CLASSES}
        )

add_executable(testGeolib ${TEST_SRC})
target_include_directories(testGeolib PUBLIC ${geolib_INCLUDE_DIRS} ${geolib_test_INCLUDE_DIRS})
target_link_libraries(testGeolib geolib)
target_compile_definitions(testGeolib 
        PRIVATE 
        FILEROOT=".")
set_target_properties(testGeolib 
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin
        )
