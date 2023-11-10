# Use this to build an examples binary
cmake_minimum_required(VERSION 3.14)
project(GEOLIB_EXAMPLE VERSION 0.0.1)

set(CMAKE_CXX_STANDARD 20)

include(${PROJECT_SOURCE_DIR}/../../.cmake/CPM.cmake)

CPMAddPackage(
        NAME geolib_library
        GITHUB_REPOSITORY ddrake-mitre/geodetic_library
        VERSION 3.2.7-SNAPSHOT
        GIT_TAG 2a8ee52
)

# Define source files for this executable
set(EXAMPLES_MAIN_SRC
        ${PROJECT_SOURCE_DIR}/main.cpp)

add_executable(examples ${EXAMPLES_MAIN_SRC})
target_include_directories(examples PUBLIC ${geolib_library_SOURCE_DIR}/include)
target_link_libraries(examples geolib)
set_target_properties(examples PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin
        )
