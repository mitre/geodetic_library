cmake_minimum_required(VERSION 3.14)
project(GEOLIB VERSION 3.2.10.0) # the nano value is a boolean. 1 == SNAPSHOT, 0 == release

set (CMAKE_CXX_STANDARD 11)

string(TOLOWER "${PROJECT_NAME}" PROJECT_NAME_LOWER)

option(BUILD_TESTING "Build the testGeolib target" ON)
set(SOLVER_TYPE "VINCENTY" CACHE STRING "Solver type desired for the build: VINCENTY or KARNEY")

# Include a suite of helper functions for output and building the aaesim package
include(.cmake/functions.cmake)

# Append -SNAPSHOT to the version name if this is a pre-release version
if (${PROJECT_VERSION_TWEAK})
    set(geolib_VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}-SNAPSHOT)
else()
    set(geolib_VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH})
endif()
message(STATUS "Generating build for ${PROJECT_NAME} version ${geolib_VERSION}")

# Add stuff to turn on code coverage only if a Debug build
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    include(${PROJECT_SOURCE_DIR}/.cmake/CodeCoverage.cmake)
    APPEND_COVERAGE_COMPILER_FLAGS()
endif()

# Gather git info if git exists
find_package(Git)
if(GIT_FOUND)
    message("Git found: ${GIT_EXECUTABLE}")
    include (.cmake/GetGitRevisionDescription.cmake)
    get_git_head_revision(GIT_BRANCH GIT_HASH)
    git_local_changes(GIT_LOCAL_CHANGES)
endif()

# Write the build_info
set(USERN $ENV{USER})
set(HOSTNAME $ENV{HOSTNAME})
string(TIMESTAMP BUILDDATETIME %Y/%m/%d-%H:%M:%S) # example: 20160101.010101
configure_file(${PROJECT_SOURCE_DIR}/.cmake/build_info.h.template ${PROJECT_SOURCE_DIR}/include/geolib/geolib_build_info.h)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-strict-aliasing")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wparentheses -Wno-sign-compare -O0 -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall -Wparentheses -Wno-sign-compare -O3 -g3")

# Create some common paths
set (GEOLIB_SRC_DIR
        ${PROJECT_SOURCE_DIR}/src/c)
set (geolib_INCLUDE_DIRS
        ${PROJECT_SOURCE_DIR}/include/geolib)

add_subdirectory(${GEOLIB_SRC_DIR})
if (BUILD_TESTING)
    set (geolib_test_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/include/test)
    set (TEST_SRC_DIR ${PROJECT_SOURCE_DIR}/test)
    add_subdirectory(${TEST_SRC_DIR})
else()
    message(STATUS "Skipping tests for ${PROJECT_NAME}")
endif()

# Install and Export for downstream users of the library
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME_LOWER}-config-version.cmake"
    COMPATIBILITY AnyNewerVersion
)
configure_package_config_file(
    ".cmake/${PROJECT_NAME_LOWER}-config.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME_LOWER}-config.cmake"
    INSTALL_DESTINATION lib/cmake/${PROJECT_NAME_LOWER}
)
install(
    TARGETS geolib
    EXPORT ${PROJECT_NAME_LOWER}-targets
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    PUBLIC_HEADER DESTINATION include
)
install(
    EXPORT ${PROJECT_NAME_LOWER}-targets
    FILE ${PROJECT_NAME_LOWER}-targets.cmake
    NAMESPACE mitre::
    DESTINATION lib/cmake/${PROJECT_NAME_LOWER}
)
