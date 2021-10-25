/*
 * Copyright 2007-2011 The MITRE Corporation.  All Rights reserved.
 *
 * This is the copyright work of The MITRE Corporation and was produced
 * for the U.S. Government under Contract Number DTFAWA-10-C-00080
 * and is subject to Federal Aviation Administration Acquisition
 * Management System Clause 3.5-13, Rights in Data-General, Alt. III
 * and Alt. IV (Oct. 1996).  No other use than that granted to the U.S.
 * Government, or to those acting on behalf of the U.S. Government, under
 * that Clause is authorized without the express written permission of
 * The MITRE Corporation.  
 * For further information, please contact The MITRE Corporation,
 * Contracts Office, 7515 Colshire Drive, McLean, VA 22102 (703) 983-6000.
 */


/* A NOTE ABOUT PARAMETERS AND UNITS */
/* Unless otherwise noted below, all distances are in nautical miles and
 * all angles, latitudes, longitudes, courses, and azimuths are in radians.
 * Standard conversions are:
 *   nautical miles to feet: multiply by 1852/0.3048 or use FEET_PER_NMI macro defined
 *   radians to degrees: multiply by 180.0/pi.
 */


#include "testUtil.h"

#ifndef EPS
#define EPS 1.0e-20
#endif

#ifndef TOL
#define TOL 1.37e-9
#endif

/* Print test cases that passed to the console */
#ifndef PRINT_PASSED_CASES
#define PRINT_PASSED_CASES 0  //0 = Do not print to console, 1 = Print to Console
#endif

/*Solves issue of strcasecmp giving errors in VS*/
#ifdef _MSC_VER
#define strcasecmp    _stricmp
#define strncasecmp   _strnicmp
#else
#include "strings.h"
#endif

namespace geolib_idealab {
/*
 * geolibTestData function declarations
 */

double randLat();

double randLon();

double randAzimuth();

double randDist();

double randDouble(double minDouble, double maxDouble);

double randSignedDist();

double randSlope();

/*
 * testUtilityFunctions declarations
 */

TestSet testFindSetMaxAndMin_Set1();

TestSuite testFindSetMaxAndMin_AllSets();

/*
 * testShape function declarations
 */

TestSet testAddPointToLLPointSet_Set1();

/*
 * testLLPoint function declarations
 */

TestSet testPtsAreSame_Set1();

TestSuite testPtsAreSame_AllSets();

TestSet testSphereInverse_Set1();

TestSuite testSphereInverse_AllSets();

TestSet testsphereInvDist_Set1();

TestSuite testsphereInvDist_AllSets();

TestSet testsphereInvCrs_Set1();

TestSuite testsphereInvCrs_AllSets();

TestSet testInverse_Set1();

TestSuite testInverse_AllSets();

TestSet testInvCrs_Set1();

TestSuite testInvCrs_AllSets();

TestSet testInvDist_Set1();

TestSuite testInvDist_AllSets();

TestSet testDirect_Set1();

TestSuite testDirect_AllSets();

TestSet testDirectLat_Set1();

TestSuite testDirectLat_AllSets();

TestSet testDirectLon_Set1();

TestSuite testDirectLon_AllSets();

TestSet testDirectInverseConsistency();

TestSet testDirectInverseMathematicaData();

/*
 * testArc function declarations
 */

TestSet testInitArcIntx_Set1();

TestSuite testInitArcIntx_AllSets();

TestSet testPtIsOnArc_Set1();

TestSuite testPtIsOnArc_AllSets();

TestSet testPtIsInsideArc_Set1();

TestSuite testPtIsInsideArc_AllSets();

TestSet testArcIntx_Set1();

TestSet testArcIntx_Set2();

TestSet testArcIntx_Set3();

TestSuite testArcIntx_AllSets();

TestSet testArcLength_Set1();

TestSuite testArcLength_AllSets();

TestSet testArcFromLength_Set1();

TestSuite testArcFromLength_AllSets();

TestSet testGetArcExtent_Set1();

TestSuite testGetArcExtent_AllSets();

TestSet testPtsOnArcOnTanThruPt_Set1();

TestSuite testPtsOnArcOnTanThruPt_AllSets();

TestSet testRefSoftGetArcExtent_Set1();

TestSuite testRefSoftGetArcExtent_AllSets();

TestSet testRefSoftPtIsInsideArc_Set1();

TestSuite testRefSoftPtIsInsideArc_AllSets();

TestSet testRefSoftPtIsOnArc_Set1();

TestSuite testRefSoftPtIsOnArc_AllSets();

TestSet testGeoTanToTwoCircles_Set1();

TestSuite testGeoTanToTwoCircles_AllSets();

TestSet testComputeSubtendedAngle_Set1();

TestSuite testComputeSubtendedAngle_AllSets();

TestSet testGeoTanToArcAtAngleToGeo_Set1();

TestSet testGeoTanToArcAtAngleToGeo_Set2();

TestSuite testGeoTanToArcAtAngleToGeo_AllSets();

TestSet testArcTanToArcAndGeo_Set1();

TestSuite testArcTanToArcAndGeo_AllSets();

TestSet testArcFromStartAndEnd_Set1();

TestSuite testArcFromStartAndEnd_AllSets();

TestSet testArcTanToCrs_Set1();

TestSuite testArcTanToCrs_AllSets();

TestSet testArcEndFromStartAndRadius_Set1();

TestSuite testArcEndFromStartAndRadius_AllSets();

TestSet testArcEndFromStartAndCenter_Set1();

TestSuite testArcEndFromStartAndCenter_AllSets();

/*
 * testGeodesic function declarations
 */

TestSet testMinSubtendedAngle_Set1();

TestSuite testMinSubtendedAngle_AllSets();

TestSet testCrsIntx_Set1();

TestSuite testCrsIntx_AllSets();

TestSet testGeoIntx_Set1();

TestSuite testGeoIntx_AllSets();

TestSet testProjectToGeo_Set1();

TestSuite testProjectToGeo_AllSets();

TestSet testGeoArcIntx_Set1();

TestSet testGeoArcIntx_Set2();

TestSet testGeoArcIntx_Set3();

TestSuite testGeoArcIntx_AllSets();

TestSet testArcTanToTwoGeos_Set1();

TestSet testArcTanToTwoGeos_Set2();

TestSuite testArcTanToTwoGeos_AllSets();

TestSet testPtIsOnGeo_Set1();

TestSuite testPtIsOnGeo_AllSets();

TestSet testGeoCrs_Set1();

TestSuite testGeoCrs_AllSets();

TestSet testProjectArcTanPtsToGeo_Set1();

TestSuite testProjectArcTanPtsToGeo_AllSets();

TestSet testLocusGeoIntx_Set1();

TestSuite testLocusGeoIntx_AllSets();

TestSet testProjectToGeoAtAngle_Set1();

TestSuite testProjectToGeoAtAngle_AllSets();

TestSet testPtIsOnCrs_Set1();

TestSuite testPtIsOnCrs_AllSets();

/*
 * testLocus function declarations
 */

TestSet testCreateLocus_Set1();

TestSuite testCreateLocus_AllSets();

TestSet testDistToLocusFromGeoDist_Set1();

TestSuite testDistToLocusFromGeoDist_AllSets();

TestSet testDistToLocusFromGeoPt_Set1();

TestSuite testDistToLocusFromGeoPt_AllSets();

TestSet testPtOnLocusFromGeoPt_Set1();

TestSuite testPtOnLocusFromGeoPt_AllSets();

TestSet testPtIsOnLocus_Set1();

TestSuite testPtIsOnLocus_AllSets();

TestSet testLocusArcIntx_Set1();

TestSet testLocusArcIntx_Set2();

TestSuite testLocusArcIntx_AllSets();

TestSet testLocusIntx_Set1();

TestSet testLocusIntx_Set2();

TestSuite testLocusIntx_AllSets();

TestSet testArcTanToTwoLoci_Set1();

TestSuite testArcTanToTwoLoci_AllSets();

TestSet testLocusCrsAtPt_Set1();

TestSuite testLocusCrsAtPt_AllSets();

TestSet testProjectToLocus_Set1();

TestSuite testProjectToLocus_AllSets();

TestSet testLociCoincide_Set1();

TestSuite testLociCoincide_AllSets();


/*
 * testBoundary function declarations
 */

TestSet testBndryCircleIntx_Set1();

TestSet testBndryCircleIntx_Set2();

TestSet testBndryCircleIntx_Set3();

TestSuite testBndryCircleIntx_AllSets();

TestSet testBndryCircleIntxExists_Set1();

TestSet testBndryCircleIntxExists_Set2();

TestSuite testBndryCircleIntxExists_AllSets();

TestSet testBndryCircleIntxExists_Set1();

TestSuite testBndryIntxExists_AllSets();

TestSet testInsertionSort_Set1();

TestSuite testInsertionSort_AllSets();

TestSet testArcsCoincide_Set1();

TestSuite testArcsCoincide_AllSets();

TestSet testPtIsInsideBndry_Set1();

TestSet testPtIsInsideBndry_Set2();

TestSet testPtIsInsideBndry_Set3();

TestSet testPtIsInsideBndry_Set4();

TestSuite testPtIsInsideBndry_AllSets();

TestSet testPtsAreInsideBndry_Set1();

TestSet testPtsAreInsideBndry_Set2();

TestSet testPtsAreInsideBndry_Set3();

TestSet testPtsAreInsideBndry_Set4();

TestSuite testPtsAreInsideBndry_AllSets();

TestSet testBndryGeoIntx_Set1();

TestSet testBndryGeoIntx_Set2();

TestSet testBndryGeoIntx_Set3();

TestSet testBndryGeoIntx_Set4();

TestSuite testBndryGeoIntx_AllSets();

TestSet testBndryLocusIntx_Set1();

TestSet testBndryLocusIntx_Set2();

TestSuite testBndryLocusIntx_AllSets();

TestSet testBndryArcIntx_Set2();

TestSet testBndryArcIntx_Set3();

TestSet testBndryArcIntx_Set4();

TestSuite testBndryArcIntx_AllSets();

TestSet testOrderBndry_Set1();

TestSuite testOrderBndry_AllSets();

TestSet testSeparateBndry_Set1();

TestSet testSeparateBndry_Set2();

TestSuite testSeparateBndry_AllSets();

ErrorSet constructLocusArcBoundary(LLPoint geoStart, double geocrs, double geolen, double locusDist, Boundary *b);

ErrorSet constructGeoBoundary(LLPoint center, double radius, double crs1, double crs2, double crs3, double crs4, Boundary *b);

ErrorSet construct2Arc2LocusBoundary(LLPoint center, double innerRadius, double outerRadius, double startSubAngle, double endSubAngle, double length, Boundary *b);

ErrorSet createLocusSpiralBndry(LLPoint sp, double lineAz, double len, double startRad, double endRad, Boundary *b); 

/*
 * testComplexBoundary function declarations
 */

TestSet testNewComplexBoundary_Set1();

TestSuite testNewComplexBoundary_AllSets();

TestSet testlongitudinallyPartitionBoundary_Set1();

TestSet testComplexBoundaryCircleIntersectionExists_Set1();

TestSet testBoundaryCircleListIntersections_Set1();

TestSuite testSpiralMidChord_AllSets();

TestSet testSpiralMidChord_Set1();

TestSuite testSpiralGeoIntx_AllSets();

TestSet testSpiralGeoIntx_Set1();

TestSuite testSpiralLocusIntx_AllSets();

TestSet testSpiralLocusIntx_Set1();

TestSuite testPtsOnSpiralOnTanThruPt_AllSets();

TestSet testPtsOnSpiralOnTanThruPt_Set1();

TestSuite testGeoTanToSpiralAtAngleToGeo_AllSets();

TestSet testGeoTanToSpiralAtAngleToGeo_Set1();

TestSuite testProjectToSpiral_AllSets();

TestSet testProjectToSpiral_Set1();

TestSuite testGeoTanToTwoSpirals_AllSets();

TestSet testGeoTanToTwoSpirals_Set1();

void testSpiralArcIntx_Set1();

void testSpiralIntx_Set1();

void testBug33911();

} // namespace
