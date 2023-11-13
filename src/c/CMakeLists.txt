cmake_minimum_required(VERSION 3.14)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")

# -----------
# Select one solver algorithm that will be used by the code.
# For AAESim, we will use Vincenty.
#
# NOTE: Vincenty is considered less accurate but faster.
# 		Karney is considered more accurate but slower.
add_compile_definitions(VINCENTY)
#add_compile_definitions(KARNEY)
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

target_include_directories(geolib PUBLIC ${geolib_INCLUDE_DIRS})

target_compile_definitions(geolib PUBLIC -D_LINUX_ -Dunix)
