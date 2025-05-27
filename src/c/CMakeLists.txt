cmake_minimum_required(VERSION 3.14)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")

# -----------
# Select one solver algorithm that will be used by the code.
# For AAESim, we will use Vincenty.
#
# NOTE: Vincenty is considered less accurate but faster.
# 		Karney is considered more accurate but slower.
if(NOT SOLVER_TYPE MATCHES "^(VINCENTY|KARNEY)$")
    message(FATAL_ERROR "Invalid SOLVER_TYPE: ${SOLVER_TYPE}. Must be either VINCENTY or KARNEY.")
endif()
message(STATUS "Using solver type: ${SOLVER_TYPE}")
add_compile_definitions(${SOLVER_TYPE})
# -----------

set(GEOLIB_SOURCE_FILES
    # NOTE: the file extension of source is important to CMAKE. It makes decisions
    #       based on them, not well documented. Be careful about making changes to the file extensions.
    Display.cc
    ErrorCodes.cc
    Geodesic.cc
    LLPoint.cc
    Util.cc
    VincentyAlgorithms.cc

# Uncomment these as needed to expand library functionality.
    Arc.cxx              # COMPILES. UNTESTED
#	Boundary.cxx         HAS COMPILE ERRORS
#	ComplexBoundary.cxx  HAS COMPILE ERRORS
#	Configuration.cxx    COMPILES. UNTESTED
    Constructions.cxx    # COMPILES. UNTESTED
#	insideOutside.cxx    COMPILES. UNTESTED
    Intersections.cxx    # COMPILES. UNTESTED
#	KarneyAlgorithms.cxx COMPILES. UNTESTED
    Locus.cxx            # COMPILES. UNTESTED
#	NearBoundary.cxx     COMPILES. UNTESTED
#	NGSVincentyAlgorithms.cxx    COMPILES. UNTESTED
    Projections.cxx      # COMPILES. UNTESTED
    Sphere.cxx           # COMPILES. UNTESTED
    Spiral.cxx           # COMPILES. UNTESTED
    Vector.cxx           # COMPILES. UNTESTED
    )

set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

add_library(geolib STATIC ${GEOLIB_SOURCE_FILES})
add_library(mitre::geolib ALIAS geolib)

target_include_directories(geolib PUBLIC 
    $<BUILD_INTERFACE:${geolib_INCLUDE_DIRS}>
    $<INSTALL_INTERFACE:include>
)

target_compile_definitions(geolib PUBLIC -D_LINUX_ -Dunix)
