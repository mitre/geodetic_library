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

/** @file Util.h */

#pragma once

#include "ErrorCodes.h"
#include "Shape.h"
#include "Constants.h"

#ifndef FILEROOT
#define FILEROOT "."
#endif

/* Checks for use of MS compiler. Converts UNIX specific string commands to MS specific
 * string commands. */
#ifdef _MSC_VER
#define strcasecmp    _stricmp
#define strncasecmp   _strnicmp
#else
#include "strings.h"
#endif

#ifndef DATA_FORMATS
/* These format strings are used to read data from files that were written using the
 * createPtString, createGeoString, createArcString, or createLocusString functions
 * in DATA_OUT mode.  The writing formats employed in these functions are essentially the same, with
 * additional precision modifiers. */
#define LLPOINT_DATA_FORMAT "%lf%c%lf%c"
#define GEODESIC_DATA_FORMAT LLPOINT_DATA_FORMAT LLPOINT_DATA_FORMAT "%d"
#define ARC_DATA_FORMAT LLPOINT_DATA_FORMAT LLPOINT_DATA_FORMAT LLPOINT_DATA_FORMAT "%d"
#define LOCUS_DATA_FORMAT LLPOINT_DATA_FORMAT LLPOINT_DATA_FORMAT "%lf%lf%d"
#define DATA_FORMATS
#endif

#ifndef MAX_DATA_LINE_LENGTH
#define MAX_DATA_LINE_LENGTH 1024
#endif

#ifndef OUTPUTMODE
typedef enum {
    /* Values given explicitly for clarity */
    SYSTEM_OUT = 0,
    MATLAB_OUT = 1,
    JAVA_OUT = 2,
    DATA_OUT = 3
} OutputMode;
#define OUTPUTMODE
#endif

namespace geolib_idealab {
/** The signum function
 * @param x Number to determine sign of (double)
 * @return Returns the sign of the input number
 * @retval 1 Positive number
 * @retval -1 Negative number
 */

int sgn(double x);

/** Calculates the reciprocal of a number
 * @param x Number to find the reciprocal of (double)
 * @return Returns the reciprocal of the input number
 * @retval 0 <= return <= 2*PI
 */
double reciprocal(double x);

/** Calculate the angle subtended by an arc, taking into account its orientation
 * @param startCrs Initial azimuth (double)
 * @param endCrs Final azimuth (double)
 * @param orient Arc orientation (ArcDirection)
 * @return Returns the angle subtended between the two courses in the direction of the arc
 * @retval 0 <= return <= \f$2\pi\f$
 */
double computeSubtendedAngle(double startCrs, double endCrs, ArcDirection orient);

/** Calculate the angle subtended between two azimuths
* @param crs1 First azimuth in radians (double)
* @param crs2 Second azimuth in radians (double)
* @param angle Pointer to double that gets updated with value
*               of angle between courses (double*)
* @return Returns Error code that indicates success or cause of failure and updates
* 				pointer with calculated angle.
* @retval SUCCESS is always returned (no errors currently possible)
*/
ErrorSet minSubtendedAngle(double crs1, double crs2, double *angle);

/**
 * Compute radius of curvature in the prime vertical at latitude and longitude of input location.
 * @param p LLPoint location
 * @returns Curvature in prime vertical (\f$N\f$)
 */
double findN(LLPoint p);

/**
 * Compute meriodional radius of curvature at latitude and longitude of input location
 * @param p LLPoint location
 * @returns Meriodional curvature (\f$M\f$)
 */
double findM(LLPoint p);

/**
 * Based on linear approximation to error function evaluated at two locations \f$y(x_{n-1})\f$, \f$y(x_n)\f$, extrapolate to find root of 
 * approximation \f$x_{n+1}\f$ where \f$y(x_{n+1})=0\f$.
 * @param x  two-element double array Two element array containing two \f$x\f$ values, \f$x_{n-1}\f$, \f$x_n\f$(double*)
 * @param y  two-element array array containing value of error function \f$y\f$ at corresponding \f$x\f$ values
 * @param err  pointer to ErrorSet (double*)
 * @returns new estimate of error function root
 */
double findRootSecantMethod(double *x, double *y, ErrorSet *err);

/** Converts geocentric latitude to geodetic latitude
 * @param lat Geocentric latitude (double)
 * @return Returns the converted geodetic latitude
 */
double geodeticLat(double lat);

/** Convert geodetic latitude to geocentric latitude
 * @param lat Geodetic latitude (double)
 * @return Returns the converted geocentric latitude
 */
double geocentricLat(double lat);

/** Converts a LLPoint from geodetic to geocentric coordinates
 * @param pt Point in geodetic coordinates (LLPoint)
 * @return Returns the converted LLPoint
 * @retval LLPoint in geocentric coordinates
 */
LLPoint geodeticToGeocentric(LLPoint pt);

/** Converts an LLPoint from geocentric to geodetic coordinates
 * @param pt  (LLPoint)
 * @return Returns the converted LLPoint
 * @retval LLPoint in geocentric coordinates
 */
LLPoint geocentricToGeodetic(LLPoint pt);

/** Converts degrees minutes seconds to radians
 * @param d Degrees (int)
 * @param m Minutes (int)
 * @param s Seconds (double)
 * @return Returns the radian value of the input degree/minute/second
 */
double dms2rad(int d, int m, double s);

/** Converts degrees minutes seconds to radians.  Takes a hemisphere input to correct
 * the radian value accordingly.
 * @param d Degrees (double)
 * @param m Minutes (double)
 * @param s Seconds (double)
 * @param hemi Hemisphere.  Possible values are 'N', 'E', 'W', 'S' (char*)
 * @return Returns the radian value of the input degree/minute/second
 */
double dms2radHemi(double d, double m, double s, char *hemi);

/** Converts a radian value to degrees/minutes/seconds
 * @param crs Course in radians (double)
 * @param *d Pointer to an int to be updated with degrees (int)
 * @param *m Pointer to an int to be updated with minutes (int)
 * @param *s Pointer to an int to be updated with seconds (double)
 * @return Updates the three input pointers with calculated values
 * @retval Nothing
 */
void rad2dms(double crs, int *d, int *m, double *s);

/** Converts the input course to degrees/minutes/seconds
 * @param crs  (double)
 * @param *d  (int)
 * @param *m  (int)
 * @param *s  (double)
 * @return Updates the three input pointers with calculated values
 * @retval Nothing
 */
void crsintrad2dms(double crs, int *d, int *m, double *s);

/** More accurate version of the signum function
 * @param x  (double)
 * @return Returns the sign of the input number.  Handles zero case correctly.
 */
double crsintsgn(double x);

/** Determines the lower of the two input values
 * @param value1 (double)
 * @param value2 (double)
 * @return Returns the lower of the two input values
 */
double minimum(double value1, double value2);

/** Determines the greater of the two input values
 * @param value1 (double)
 * @param value2 (double)
 * @return Returns the greater of the two input values
 */
double maximum(double value1, double value2);

/** Steps through an array of values and finds the minimum and maximum values
 * @param values[] Array of values to check (double)
 * @param size Length of input array (int)
 * @param max Pointer to double to be updated with maximum value from list (double*)
 * @param min Pointer to double to be updated with minimum value from list (double*)
 * @return Returns Error code that indicates success or cause of failure and updates
 * 				pointers with minimum and maximum values in the list
 * @retval SUCCESS is always returned (no errors currently possible)
 */
ErrorSet findSetMaxAndMin(double values[], int size, double *max, double *min);

/** Calculates the modulus of the input course
 * @param crs (double)
 * @return Returns the modulus of the input course
 */
double modcrs(double crs);

/** Map azimuth or angle values from [0,2*pi] to [-pi,pi]
 * @param x (double)
 * @return Returns calculated latitude
 */
double modlat(double x);

/** Map azimuth or angle values from [0,2*pi] to [-pi,pi]
 * @param x (double)
 * @return Returns calculated value
 */
double modlon(double x);

/** Return positive remainder of x / y
 * @param x (double)
 * @param y (double)
 * @return Returns calculated value
 */
double modpos(double x, double y);

/* Functions in Vector.c */

/** Converts an LLPoint object to a vector in the ECEF coordinate system
 * @param geo  (LLPoint)
 * @return Returns a vector representing the LL position of the LLPoint
 */
Vector geodeticToECEF(LLPoint geo);

/** Calculates the cross product of two vectors
 * @param v1 (Vector)
 * @param v2 (Vector)
 * @return Returns the cross product of two input vectors
 * @retval Cross product
 */
Vector cross(Vector v1, Vector v2);

/** Calculates the dot product of two vectors
 * @param v1 (Vector)
 * @param v2 (Vector)
 * @return Returns the dot product of two input vectors
 * @retval Dot product
 */
double dot(Vector v1, Vector v2);

/** Calculates the length of a vector
 * @param v (Vector)
 * @return Returns the calculated length of the vector
 */
double norm(Vector v);

/** Maps a point onto the spherical earth model
 * @param p (LLPoint)
 * @return Returns a vector representing the point on the spherical earth model
 */
Vector mapToUnitSphere(LLPoint p);

/** Normalizes a vector such that it's length = 1
 * @param v (Vector*)
 * @return Sets the vector's properties such the direction is unchanged and
 * 				length = 1
 * @retval Nothing
 */
void normalize(Vector *v);

/** Multiplies a vector by a scalar
 * @param v (Vector*)
 * @param scal (Vector*)
 * @return Sets the vector's properties according to recalculated length
 * @retval Nothing
 */
void scalarMultiply(Vector *v, double scal);

/** Adds two vectors
 * @param v1  (Vector)
 * @param v2  (Vector)
 * @return returns the sum of the two input vectors
 */
Vector vectorAdd(Vector v1, Vector v2);

/** Subtracts two vectors
 * @param v1 (Vector)
 * @param v2 (Vector)
 * @return Returns the difference of the two input vectors
 */
Vector vectorSubtract(Vector v1, Vector v2);

/** Maps a vector onto the spherical earth model
 * @param v  (Vector)
 * @return Returns an LLPoint representing the vector on the spherical earth model
 */
LLPoint mapVectorToSphere(Vector v);

/* Functions in Display.c */

/** Creates a string describing the properties of an LLPoint in the desired format
 * @param p Point to be described (LLPoint)
 * @param pointName Name of point (optional) (char*)
 * @param mode Output mode.  Options are SYSTEM_OUT and MATLAB_OUT (OutputMode)
 * @param displayRadians Flag to determine output in degrees (0) or radians (1) (int)
 * @return Returns a string containing properties of the input LLPoint
 */
char *createPtString(LLPoint p, char *pointName, OutputMode mode, int displayRadians);

/** Creates a string describing the properties of a Geodesic in the desired format
 * @param g Geodesic to be described (Geodesic)
 * @param geoName Name of geodesic (optional) (char*)
 * @param mode Output mode.  Options are SYSTEM_OUT and MATLAB_OUT (OutputMode)
 * @param displayRadians Flag to determine output in degrees (0) or radians (1) (int)
 * @return Returns a string containing properties of the input Geodesic
 */
char *createGeoString(Geodesic g, char *geoName, OutputMode mode, int displayRadians);

/** Creates a string describing the properties of a Locus in the desired format
 * @param l Locus to be described (Locus)
 * @param locusName Name of locus (optional) (char*)
 * @param mode Output mode.  Options are SYSTEM_OUT and MATLAB_OUT (OutputMode)
 * @param displayRadians Flag to determine output in degrees (0) or radians (1) (int)
 * @return Returns a string containing properties of the input Locus
 */
char *createLocusString(Locus l, char *locusName, OutputMode mode, int displayRadians);

/** Creates a string describing the properties of an Arc in the desired format
 * @param a Arc to be described (Arc)
 * @param arcName Name of the arc (optional) (char*)
 * @param mode Output mode.  Options are SYSTEM_OUT and MATLAB_OUT (OutputMode)
 * @param displayRadians Flag to determine output in degrees (0) or radians (1) (int)
 * @return Returns a string containing properties of the input Arc
 */
char *createArcString(Arc a, char *arcName, OutputMode mode, int displayRadians);

/** Creates a string describing the properties of a Spiral in the desired format
 * @param a Spiral to be described (Spiral)
 * @param spiralName Name of the spiral (optional) (char*)
 * @param mode Output mode.  Options are SYSTEM_OUT and MATLAB_OUT (OutputMode)
 * @param displayRadians Flag to determine output in degrees (0) or radians (1) (int)
 * @return Returns a string containing properties of the input Spiral
 */
char *createSpiralString(Spiral s, char *spiralName, OutputMode mode, int displayRadians);

/** Creates a string describing the properties of a Boundary in the desired format
 * @param b Boundary to be described (Boundary)
 * @param boundaryName Name of the boundary (optional) (char*)
 * @param mode Output mode.  Options are SYSTEM_OUT and MATLAB_OUT (OutputMode)
 * @param displayRadians Flag to determine output in degrees (0) or radians (1) (int)
 * @return Returns a string containing properties of the input Boundary
 */
char *createBndryString(Boundary b, char *boundaryName, OutputMode mode, int displayRadians);

/** Creates a string describing the properties of a Complex Boundary in the desired format
 * @param c Complex Boundary to be described (ComplexBoundary)
 * @param ComplexBoundaryName Name of the complex boundary (optional) (char*)
 * @param mode Output mode.  Options are SYSTEM_OUT and MATLAB_OUT (OutputMode)
 * @param displayRadians Flag to determine output in degrees (0) or radians (1) (int)
 * @return Returns a string containing properties of the input Complex Boundary
 */
char *createComplexBndryString(ComplexBoundary c, char *ComplexBoundaryName, OutputMode mode, int displayRadians);

/** Prints the properties of an LLPoint to the console
 * @param p LLPoint to be described (LLPoint)
 * @param pointName Name of the LLPoint (optional) (char*)
 * @param displayRadians Flag to determine output in degrees (0) or radians (1) (int)
 * @return Prints the properties of the LLPoint to the system console
 * @retval Nothing
 */
void displayPt(LLPoint p, char *pointName, int displayRadians);

/** Prints the properties of a Geodesic to the console
 * @param g Geodesic to be described (Geodesic)
 * @param geoName Name of the Geodesic (optional) (char*)
 * @param displayRadians Flag to determine output in degrees (0) or radians (1) (int)
 * @return Prints the properties of the Geodesic to the system console
 * @retval Nothing
 */
void displayGeo(Geodesic g, char *geoName, int displayRadians);

/** Prints the properties of a Locus to the console
 * @param l Locus to be described (Locus)
 * @param locusName Name of the locus (optional) (char*)
 * @param displayRadians Flag to determine output in degrees (0) or radians (1) (int)
 * @return Prints the properties of the Locus to the system console
 * @retval Nothing
 */
void displayLocus(Locus l, char *locusName, int displayRadians);

/** Prints the properties of an Arc to the console
 * @param a Arc to be described (Arc)
 * @param arcName Name of the arc (optional) (char*)
 * @param displayRadians displayRadians Flag to determine output in degrees (0) or radians (1) (int)
 * @return Prints the properties of the Arc to the system console
 * @retval Nothing
 */
void displayArc(Arc a, char *arcName, int displayRadians);

/** Prints the properties of a Spiral to the console
 * @param a Spiral to be described (Spiral)
 * @param spiralName Name of the spiral (optional) (char*)
 * @param displayRadians displayRadians Flag to determine output in degrees (0) or radians (1) (int)
 * @return Prints the properties of the Spiral to the system console
 * @retval Nothing
 */
void displaySpiral(Spiral s, char *spiralName, int displayRadians);

/** Prints the properties of a Boundary to the console
 * @param b Boundary to be described (Boundary)
 * @param boundaryName Name of the boundary (char*)
 * @param displayRadians displayRadians Flag to determine output in degrees (0) or radians (1) (int)
 * @return Prints the properties of the Boundary to the system console
 * @retval Nothing
 */
void displayBndry(Boundary b, char *boundaryName, int displayRadians);

/** Prints the properties of a Complex Boundary to the console
 * @param c Complex Boundary to be described (ComplexBoundary)
 * @param complexBoundaryName Name of the complex boundary (optional) (char*)
 * @param displayRadians displayRadians Flag to determine output in degrees (0) or radians (1) (int)
 * @return Prints the properties of the Complex Boundary to the system console
 * @retval Nothing
 */
void displayComplexBndry(ComplexBoundary c, char *complexBoundaryName, int displayRadians);

/** Prints the properties of an LLPoint to the console in a Matlab friendly format
 * @param p LLPoint to be described (LLPoint)
 * @param pointName Name of the LLPoint (optional) (char*)
 * @param displayRadians displayRadians Flag to determine output in degrees (0) or radians (1) (int)
 * @return Prints the properties of the LLPoint to the system console in the proper format
 * @retval Nothing
 */
void displayMatlabPt(LLPoint p, char *pointName, int displayRadians);

/** Prints the properties of a Geodesic to the console in a Matlab friendly format
 * @param g Geodesic to be described (Geodesic)
 * @param geoName Name of the Geodesic (optional) (char*)
 * @param displayRadians displayRadians Flag to determine output in degrees (0) or radians (1) (int)
 * @return Prints the properties of the Geodesic to the system console in the proper format
 * @retval Nothing
 */
void displayMatlabGeo(Geodesic g, char *geoName, int displayRadians);

/** Prints the properties of a Locus to the console in a Matlab friendly format
 * @param l Locus to be described (Locus)
 * @param locusName Name of the Locus (char*)
 * @param displayRadians displayRadians Flag to determine output in degrees (0) or radians (1) (int)
 * @return Prints the properties of the Locus to the system console in the proper format
 * @retval Nothing
 */
void displayMatlabLocus(Locus l, char *locusName, int displayRadians);

/** Prints the properties of an Arc to the console in a Matlab friendly format
 * @param a Arc to be described (Arc)
 * @param arcName Name of the Arc (optional) (char*)
 * @param displayRadians displayRadians Flag to determine output in degrees (0) or radians (1) (int)
 * @return Prints the properties of the Arc to the system console in the proper format
 * @retval Nothing
 */
void displayMatlabArc(Arc a, char *arcName, int displayRadians);

/** Prints the properties of a Spiral to the console in a Matlab friendly format
 * @param a Spiral to be described (Spiral)
 * @param spiralName Name of the Spiral (optional) (char*)
 * @param displayRadians displayRadians Flag to determine output in degrees (0) or radians (1) (int)
 * @return Prints the properties of the Spiral to the system console in the proper format
 * @retval Nothing
 */
void displayMatlabSpiral(Spiral s, char *spiralName, int displayRadians);

/** Prints the properties of a Boundary to the console in a Matlab friendly format
 * @param b Boundary to be described (Boundary)
 * @param boundaryName Name of the Boundary (optional) (char*)
 * @param displayRadians displayRadians Flag to determine output in degrees (0) or radians (1) (int)
 * @return Prints the properties of the Boundary to the system console in the proper format
 * @retval Nothing
 */
void displayMatlabBndry(Boundary b, char *boundaryName, int displayRadians);

/** Prints the properties of a Complex Boundary to the console in a Matlab friendly format
 * @param c Complex Boundary to be described (ComplexBoundary)
 * @param complexBoundaryName Name of the Complex Boundary (optional) (char*)
 * @param displayRadians displayRadians Flag to determine output in degrees (0) or radians (1) (int)
 * @return Prints the properties of the Complex Boundary to the system console in the proper format
 * @retval Nothing
 */
void displayMatlabComplexBndry(ComplexBoundary c, char *complexBoundaryName, int displayRadians);

/** Prints the properties of an LLPoint to the console in a Matlab friendly format
 * @param p LLPoint to be described (LLPoint)
 * @param pointName Name of the LLPoint (optional) (char*)
 * @param displayRadians displayRadians Flag to determine output in degrees (0) or radians (1) (int)
 * @return Prints the properties of the LLPoint to the system console in the proper format
 * @retval Nothing
 */
void displayDataPt(LLPoint p, char *pointName, int displayRadians);

/** Prints the properties of a Geodesic to the console in a Data friendly format
 * @param g Geodesic to be described (Geodesic)
 * @param geoName Name of the Geodesic (optional) (char*)
 * @param displayRadians displayRadians Flag to determine output in degrees (0) or radians (1) (int)
 * @return Prints the properties of the Geodesic to the system console in the proper format
 * @retval Nothing
 */
void displayDataGeo(Geodesic g, char *geoName, int displayRadians);

/** Prints the properties of a Locus to the console in a Data friendly format
 * @param l Locus to be described (Locus)
 * @param locusName Name of the Locus (char*)
 * @param displayRadians displayRadians Flag to determine output in degrees (0) or radians (1) (int)
 * @return Prints the properties of the Locus to the system console in the proper format
 * @retval Nothing
 */
void displayDataLocus(Locus l, char *locusName, int displayRadians);

/** Prints the properties of an Arc to the console in a Data friendly format
 * @param a Arc to be described (Arc)
 * @param arcName Name of the Arc (optional) (char*)
 * @param displayRadians displayRadians Flag to determine output in degrees (0) or radians (1) (int)
 * @return Prints the properties of the Arc to the system console in the proper format
 * @retval Nothing
 */
void displayDataArc(Arc a, char *arcName, int displayRadians);

/** Prints the properties of a Spiral to the console in a Data friendly format
 * @param a Spiral to be described (Spiral)
 * @param spiralName Name of the Spiral (optional) (char*)
 * @param displayRadians displayRadians Flag to determine output in degrees (0) or radians (1) (int)
 * @return Prints the properties of the Spiral to the system console in the proper format
 * @retval Nothing
 */
void displayDataSpiral(Spiral s, char *spiralName, int displayRadians);

/** Prints the properties of a Boundary to the console in a Data friendly format
 * @param b Boundary to be described (Boundary)
 * @param boundaryName Name of the Boundary (optional) (char*)
 * @param displayRadians displayRadians Flag to determine output in degrees (0) or radians (1) (int)
 * @return Prints the properties of the Boundary to the system console in the proper format
 * @retval Nothing
 */
void displayDataBndry(Boundary b, char *boundaryName, int displayRadians);

/** Prints the properties of a Complex Boundary to the console in a Data friendly format
 * @param c Complex Boundary to be described (ComplexBoundary)
 * @param complexBoundaryName Name of the Complex Boundary (optional) (char*)
 * @param displayRadians displayRadians Flag to determine output in degrees (0) or radians (1) (int)
 * @return Prints the properties of the Complex Boundary to the system console in the proper format
 * @retval Nothing
 */
void displayDataComplexBndry(ComplexBoundary c, char *complexBoundaryName, int displayRadians);

/** Print out a list of constants that were set at compile time
 * @param fileName  (const char*)
 * @return Prints the constant list into the input file
 * @retval Nothing
 */
void reportConstants(const char *fileName);

void aboutGeolib(char *result, const unsigned int size);

} // end namespace