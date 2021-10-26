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

#include <string>
#include "Util.h"
#include "Geolib.h"

namespace geolib_idealab {

#ifndef TESTTOL
#define TESTTOL 5.0e-7
#endif

#ifndef FILEROOT
#define FILEROOT "."
#endif

#ifndef TESTTYPE
typedef enum {
	SET = 0,
	SUITE = 1
}
TestType;
#define TESTTYPE
#endif

#ifndef TESTSET
typedef struct
{
	std::string name;
    int pass;
    int fail;
    int errors;
    int unverified;
    int setupFailures;
    int testCases;
    int testingError;
}
TestSet;
#define TESTSET
#endif

#ifndef TESTELEMENT
typedef struct
{
	void* ptr;
	TestType type;
}
TestElement;
#define TESTELEMENT
#endif

#ifndef TESTSUITE
typedef struct
{
	std::string name;
	int length;
	TestElement* elements;
    int pass;
    int fail;
    int errors;
    int unverified;
    int setupFailures;
    int testCases;
	int testingSuiteError;
}
TestSuite;
#define TESTSUITE
#endif

TestSet newTestSet(std::string name);

TestSuite newTestSuite(std::string name);

void displayTestSet(TestSet set);

void displayTestSuite(TestSuite suite);

void outputTestSuiteMetrics(TestSuite suite, char* outputFile);

void growSuite(TestSuite* suite);

void addTestSet(TestSet set, TestSuite* suite);

void addTestSuite(TestSuite childSuite, TestSuite* parentSuite);

int checkForNullAndDuplicatePointers(void* firstP, ...);

int isEqualLLPoint(LLPoint pt1, LLPoint pt2);

ShapeType convertStringToShapeType(char* str, char* nullDelim, ShapeType dflt);

LineType convertStringToLineType(char* str, char* nullDelim, LineType dflt);

int quadrant(LLPoint pt, LLPoint ref);

//TODO This function is a temporary duplicate of the displayPt function.
//Once the refactoring of libWGS84.c is complete and displayPt is public
//all references to printLLPoint should be replaced with displayPt.
void printLLPoint(std::string ptName, LLPoint pt);

void newflush();

void newpause();

//TODO This function is a temporary duplicate of the displayLocus function.
//Once the refactoring of libWGS84.c is complete and displayLocus is public
//all references to printLocus should be replaced with displayLocus.
void printLocus(char* locusName, Locus locus);

double crsZeroToTwoPI(double crs);

//TODO refactor and redesign this method so that EPS and TOL are passed in arguments
//TODO refactor to handle error codes properly; use masks and this should probably return the error code(s) generated
void createLocusThroughPoint(Locus *locusPtr, LLPoint pt, double crsLoc,
		double slopeLocAngle, double d1, double d2, double d3);

int isEqualLLPointWithinBox(LLPoint pt1, LLPoint pt2);

int isEqualLatitude(double lat1, double lat2);

int isEqualLongitude(double lon1, double lon2);

void display(LLPoint p);

void displayOld(LLPoint p);

void displayOld2(LLPoint p);

void displayOld3(LLPoint p);

void parseInputFindArcLength (char* str, LLPoint* center, double* radius, double* startAz,
                double* endAz, int* dir);

void parseInputDiscretizedArcLength(char* str, LLPoint* center, double* radius, double* startAz,
                double* endAz, int* dir);

void parseInputInterceptAtAngle(char* str, LLPoint* pt1, double* crs, LLPoint* pt2, double* angle);

double parseDMSFormat(const char* text);

void displayBoundaryArc(Arc a, int displayRadians);

void displayBoundaryLocus(Locus l, int displayRadians);

void displayBoundaryGeodesic(Geodesic g, int displayRadians);

void displayBoundaryLLPoint(LLPoint p, int displayRadians);

int compareArcs(Arc arc1, Arc* arc2, double tol);

ErrorSet createGeoTestBndry(LLPoint center, double radius, Boundary* geoTestBndry);

ErrorSet createLocusTestBndry(LLPoint center, double radius, Boundary* locusTestBndry);

ErrorSet createArcTestBndry(LLPoint center, double radius, Boundary* arcTestBndry);

ErrorSet createSpiralTestBndry(LLPoint center, double radius1, double radius2, Boundary* spiralTestBndry);

} // namespace
