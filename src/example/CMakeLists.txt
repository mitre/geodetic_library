# Use this to build an examples library of source
cmake_minimum_required(VERSION 3.14)

set(CMAKE_CXX_STANDARD 11)

include(${PROJECT_SOURCE_DIR}/.cmake/CPM.cmake)

CPMAddPackage(
        NAME geolib_library
        GIT_REPOSITORY https://github.com/mitre/geodetic_library.git
        VERSION 3.2.7-SNAPSHOT
        GIT_TAG main
)

# Define source files for this executable
set(EXAMPLES_MAIN_SRC
        ${EXAMPLES_SRC_DIR}/main.cpp)

add_executable(examples ${EXAMPLES_MAIN_SRC})
target_link_libraries(examples geolib_library)
set_target_properties(examples PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin
        )
