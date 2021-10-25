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

/*! \file Geolib.h
 *  \brief This file defines the function signatures used by geolib.
 *  \author Michael Mills, Richard Snow, Stuart Bowman, Juan Amezcua, John Landrigan
 */

/*
 * A NOTE ABOUT PARAMETERS AND UNITS
 * Unless otherwise noted below, all distances are in nautical miles and
 * all angles, latitudes, longitudes, courses, and azimuths are in radians.
 * Standard conversions are:
 *   nautical miles to feet: multiply by 1852/0.3048 or use FEET_PER_NMI macro defined below
 *   radians to degrees: multiply by 180.0/pi, or use DEG_PER_RAD macro defined below.
 *
 * 
 */

#pragma once

#include "Util.h"

#ifndef LIBWGS84_H_
#define LIBWGS84_H_
#endif /*LIBWGS84_H_*/

#ifdef WIN32
#define strdup _strdup
#endif

namespace geolib_idealab {
/* @struct Configuration Geolib.h "include/Geolib.h"
* @brief Defines numerical values used by most algorithms.
*
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @param semiMajorAxis The ellipsoid's equatorial radius (in NM)
* @param flattening The ratio of the difference between the equatorial and polar radii to
* 		the equatorial radius (dimensionless)
* @param lastErr Long parameter that stores the most recent error value
*/
#ifndef CONFIGURATION_STRUCT
typedef struct {
  double tol;
  double eps; // Leave this out?
  double semiMajorAxis;
  double flattening;
  double sphereRadius;
  double internalZero;
  double smallDistThreshold;
  int iterationCount;
  ErrorSet lastErr;
} Configuration;
#define CONFIGURATION_STRUCT
#endif

ErrorSet ptOnLocus2FromGeoDist(Locus loc, double secondDist, double geoDist, LLPoint *loc2Pt, double tol, double eps);
ErrorSet locus2ndOrderMeridianIntx(Locus2ndOrder locus2ndOrder,
                                   double longitude,
                                   LLPoint *intx,
                                   double tol,
                                   double eps);
ErrorSet areTerrainPtsNearPt(LLPoint pt,
                             double nearDist,
                             double *latList,
                             double *lonList,
                             int *near,
                             int numberOfPts,
                             double tol,
                             double eps);
ErrorSet areTerrainPtsNearGeoAbeam(Geodesic geo,
                                   double nearDist,
                                   double *latList,
                                   double *lonList,
                                   int *near,
                                   int numberOfPts,
                                   double tol,
                                   double eps);
ErrorSet areTerrainPtsNearLocusAbeam(Locus locus,
                                     double nearDist,
                                     double *latList,
                                     double *lonList,
                                     int *near,
                                     int numberOfPts,
                                     double tol,
                                     double eps);
ErrorSet areTerrainPtsNearArcAbeam(Arc arc,
                                   double nearDist,
                                   double *latList,
                                   double *lonList,
                                   int *near,
                                   int numberOfPts,
                                   double tol,
                                   double eps);
ErrorSet areTerrainPtsNearSpiralAbeam(Spiral spiral,
                                      double nearDist,
                                      double *latList,
                                      double *lonList,
                                      int *near,
                                      int numberOfPts,
                                      double tol,
                                      double eps);
ErrorSet areTerrainPtsNearBoundary(Boundary b,
                                   double nearDist,
                                   double *latList,
                                   double *lonList,
                                   int *near,
                                   int *inside,
                                   int numberOfPts,
                                   double tol,
                                   double eps);
ErrorSet distBetweenGeos(Geodesic line1, Geodesic line2, double *dist, double tol, double eps);
ErrorSet distBetweenArcGeo(Arc arc, Geodesic line, double *dist, double tol, double eps);
int azIsInArcExtent(Arc arc, double az);

/* API function prototypes */


/** Check if given testPt lies on geodesic defined by startPt and crs12.  It is
* considered to be on geodesic if distance from testPt to geodesic is less than tol.
* @param startPt Start point of geodesic (LLPoint)
* @param crs12 Azimuth of geodesic at startPt (double)
* @param testPt Point to test (LLPoint)
* @param crsTest1 (optional) Pointer to azimuth from startPt to testPt (double*)
* @param dist1Test (optional) Pointer to distance from startPt to testPt (double*)
* @param err Pointer to error set for returning error codes (long*)
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return Returns an int indicating whether the input point is on the input course.
* @retval 1 is returned if testPt is within tol of course
* @retval 0 is returned if testPt is a distance greater than tol from the course
*
*/
int ptIsOnCrs(LLPoint startPt, double crs12, LLPoint testPt,
              double *crsTest1, double *dist1Test, ErrorSet *err, double tol,
              double eps);

/** Check if two geodetic positions (latitude/longitude pairs, represented by LLPoints)
* are the same to with tol (nmi) distance
* @param p1 First point to compare (LLPoint)
* @param p2 Second point to compare (LLPoint)
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @return Returns an int indicating whether the two input points are within tolerance of each other.
* @retval 1 if positions are separated by no more than tol
* @retval 0 if positions are separated by at least tol
*
*/
int ptsAreSame(LLPoint p1, LLPoint p2, double tol);

/** Determines if two points are antipodal
*
* @param p1 First input point (LLPoint)
* @param p2 Second input point (LLPoint)
* @return Returns an ErrorSet containing information on whether points are anti-podal or not.
* @retval SUCCESS Indicates successful execution.
* @retval ANTIPODAL_POINTS_ERR Indicates that two LLPoints are on opposite sides of the earth ellipsoid.
 */
ErrorSet ptsAreAntipodal(LLPoint p1, LLPoint p2);

/** Calculate a geodetic position (lat/lon in radians) given starting position,
* course and distance
* @param origin Starting position lat/lon in radians (LLPoint)
* @param course Azimuth of geodesic at origin in radians (0=north, pi/2 = east, etc.) (double)
* @param distance Distance to desired point (in nmi) (double)
* @param dest Pointer to LLPoint that will be updated with lat/lon of destination (LLPoint*)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return Returns error code that indicates success or cause of failure.  Updates given memory address with
* 				calculated lat/long values.
* NOTE: If starting latitude is at either pole (where azimuth is not uniquely
*       defined), then the point returned will have longitude equal to the
*       input course. Latitude will be determined from input distance.
* @retval SUCCESS Indicates successful execution.
 */
ErrorSet direct(LLPoint origin, double course, double distance, LLPoint *dest,
                double eps);

/** Calculate just the destination latitude given starting position, course, distance.
*
* @param origin Starting position lat/lon in radians (LLPoint)
* @param course Azimuth of geodesic at origin in radians (0=north, pi/2 = east, etc.) (double)
* @param distance Distance to desired point (in nmi) (double)
* @param lat Pointer to a double that will be updated with latitude of destination (double*)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return Returns error code that indicates success or cause of failure; updates given memory
* address with calculated latitude.
* @retval SUCCESS Indicates successful execution.
 */
ErrorSet directLat(LLPoint origin, double course, double distance, double *lat,
                   double eps);

/** Calculate only the destination longitude given starting position, course,
* distance.
* @param origin Starting position lat/lon in radians (LLPoint)
* @param course Azimuth of geodesic at origin in radians (0=north, pi/2 = east, etc.) (double)
* @param distance Distance to desired point (in nmi) (double)
* @param lon Pointer to a double that will be updated with longitude of destination (double*)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return Returns error code that indicates success or cause of failure; updates given memory
* address with calculated longitude.
* @retval SUCCESS Indicates successful execution.
 */
ErrorSet directLon(LLPoint origin, double course, double distance, double *lon,
                   double eps);

/** Calculate distance and courses between two geodetic positions
* @param origin Starting position lat/lon in radians (LLPoint)
* @param dest Ending position lat/lon in radians (LLPoint)
* @param crs Pointer to double that will be updated with
*              azimuth at origin in radians (optional) (double*)
* @param bcrs Pointer to double that will be updated with
*               reciprocal azimuth at destination in radians (optional) (double*)
* @param dist Pointer to double that will be updated with
*               distance between points in nmi. (optional) (double*)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return Returns error code that indicates success or cause of failure; updates given memory
* addresses with calculated azimuths and distance.
* @retval SUCCESS Indicates successful execution.
* @retval UNEXPECTED_ERR Indicates that an unknown error has occurred.
 */
ErrorSet inverse(LLPoint origin, LLPoint dest, double *crs, double *bcrs,
                 double *dist, double eps);

/** Calculate courses (in radians) between two geodetic positions
* @param origin Starting position lat/lon in radians (LLPoint)
* @param dest Ending position lat/lon in radians (LLPoint)
* @param fcrs Pointer to double that will be updated with
*              azimuth at origin in radians (optional) (double*)
* @param bcrs Pointer to double that will be updated with
*               reciprocal azimuth at destination in radians (optional) (double*)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return Returns error code that indicates success or cause of failure; updates given memory
* addresses with calculated azimuths
* @retval SUCCESS Indicates successful execution.
 */
ErrorSet invCrs(LLPoint origin, LLPoint dest, double *fcrs, double *bcrs,
                double eps);

/** Calculate distance (in nmi) between two geodetic positions
*
* @param origin Starting position lat/lon in radians (LLPoint)
* @param dest Ending position lat/lon in radians (LLPoint)
* @param dist Pointer to double that will be updated with
*               distance between points in nmi. (optional) (double*)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return Returns error code that indicates success or cause of failure; updates given memory
* addresses with calculated distance.
* @retval SUCCESS Indicates successful execution.
 */
ErrorSet invDist(LLPoint origin, LLPoint dest, double *dist, double eps);

/** Calculate intersection of two geodesics.  Each geodesic is defined by its starting point
* and starting azimuth, and is assumed to be unbounded on the ellipsoid.
* @param pt1 Starting position of first geodesic (LLPoint)
* @param crs13 Starting azimuth of first geodesic, in radians (double)
* @param crs31 Pointer to reciprocal azimuth of 2nd geodesic at
*                 intersection point (result, in radians) (double*)
* @param dist13 (optional) Pointer to distance from pt1 to intersection point
*                  (result, in radians) (double*)
* @param pt2 Starting position of second geodesic (LLPoint)
* @param crs23 Starting azimuth of second geodesic (double)
* @param crs32 (optional) Pointer to reciprocal azimuth of 2nd geodesic at
*                 intersection point (result, in radians) (double*)
* @param dist23 (optional) Pointer to distance from pt2 to intersection point
*                 (result, in radians) (double*)
* @param intx Pointer to LLPoint struct that will be updated with calculated
*               position. (LLPoint*)
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return Returns error code that indicates success or cause of failure; updates given memory
* addresses with calculated values.
* @retval SUCCESS Indicates successful execution.
*/
ErrorSet crsIntx(LLPoint pt1, double crs13, double *crs31,
                 double *dist13, LLPoint pt2, double crs23,
                 double *crs32, double *dist23, LLPoint *intx,
                 double tol, double eps);

/** Calculate intersection of two geodesics.  Each geodesic is defined by starting and ending
* positions.  The extent of each geodesic is specified by its LineType parameter.  LineType == 0 implies the
* geodesic exists only between its start and end points; LineType == 1 implies the geodesic begins at its start point
* but extends beyond its endpoint; LineType == 2 implies the geodesic extends beyond both points.
* @param start1 Starting position of first geodesci (LLPoint)
* @param end1 Ending position of first geodesic (LLPoint)
* @param lineType1 Value from LineType enum that specifies extent of the first geodesic (LineType)
* @param crs31 Pointer to reciprocal azimuth of 2nd geodesic at
*                 intersection point (result, in radians) (double*)
* @param dist13 Pointer to distance from pt1 to intersection point
*                  (result, in radians) (optional) (double*)
* @param start2 Starting position of second geodesic (LLPoint)
* @param end2 Ending position of second geodesic (LLPoint)
* @param lineType2 Value from LineType enum that specifies extent of the second geodesic (LineType)
* @param crs32 Pointer to reciprocal azimuth of 2nd geodesic at
*                 intersection point (result, in radians) (optional) (double*)
* @param dist23 Pointer to distance from pt2 to intersection point
*                 (result, in radians) (optional) (double*)
* @param intx  Pointer to LLPoint struct that will be updated with calculated
*               position. (LLPoint*)
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return Returns error code that indicates success or cause of failure; updates given memory
* addresses with calculated values.
* @retval SUCCESS Indicates successful execution.
 */
ErrorSet geoIntx(LLPoint start1, LLPoint end1, LineType lineType1,
                 double *crs31, double *dist13, LLPoint start2,
                 LLPoint end2, LineType lineType2, double *crs32,
                 double *dist23, LLPoint *intx, double tol, double eps);

/** Calculate the point on geodesic nearest a given point.  In other words, project
* the given point to the geodesic.
* @param pt1 Starting point of geodesic (LLPoint)
* @param crs12 azimuth of geodesic at pt1 in radians (double)
* @param pt3 point to be projected to geodesic (LLPoint)
* @param pt2 pointer to LLPoint that will be updated with coordinates of
*              projected point. (LLPoint*)
* @param crsFromPoint pointer to azimuth of geodesic from pt3 to projected
*                        point, in radians (result) (double*)
* @param distFromPoint pointer to distance from pt3 to geodesic at projected
*                         point, in nmi (result) (double*)
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return Returns error code that indicates success or cause of failure; updates given memory
* addresses with calculated values.
* @retval SUCCESS Indicates successful execution.
* @retval ITERATION_MAX_REACHED_ERR Indicates that a method has looped more than the allowed iteration count.
* @retval ERROR_MAX_REACHED_ERR Indicates that an error value has grown larger than allowed.
 */
ErrorSet projectToGeo(LLPoint pt1, double crs12, LLPoint pt3, LLPoint *pt2,
                      double *crsFromPoint, double *distFromPoint,
                      double tol, double eps);

/** Given point pt3, the azimuth of the geodesic through pt1 not containing that point, and the desired
* intersection angle intAngle, find the point pt2 on geodesic such that angle between
* given geodesic and geodesic containing pt2 and pt3 is equal to intAngle.
* @param pt1 Starting point of geodesic (LLPoint)
* @param crs12 Azimuth of geodesic at pt1 in radians (double)
* @param pt3 Point to be projected to geodesic (LLPoint)
* @param intAngle Desired intersection angle in radians (!= 0, sign not important) (double)
* @param pt2 Pointer to LLPoint that will be updated with coordinates of
*				projected point. (LLPoint*)
* @param crsFromPoint Pointer to azimuth of geodesic from pt3 to projected
*              point, in radians (result) (double*)
* @param distFromPoint Pointer to distance from pt3 to geodesic at projected
*              point, in nmi (result) (double*)
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return Returns error code that indicates success or cause of failure; updates given memory
* addresses with calculated values.
* @retval SUCCESS Indicates successful execution.
* @retval INVALID_CRS_ERR Indicates that an azimuth value was out of range (e.g., desired course from pole other than \f$\pm\pi\f$).
* @retval ITERATION_MAX_REACHED_ERR Indicates that a method has looped more than the allowed iteration count.
* @retval ERROR_MAX_REACHED_ERR Indicates that an error value has grown larger than allowed.
 */
ErrorSet projectToGeoAtAngle(LLPoint pt1, double crs12, LLPoint pt3,
                             double intAngle, LLPoint *pt2, double *crsFromPoint,
                             double *distFromPoint, double tol, double eps);

/* Calculate points of intersection of two arcs.  Arcs are treated as full
* circles. Arc bounds must be applied by parent function, when applicable.
* @param center1 Center of first arc (LLPoint)
* @param r1 Radius of first arc in nmi. (double)
* @param center2 Center of second arc (LLPoint)
* @param r2 Radius of second arc in nmi. (double)
* @param intx Two-element array of LLPoint objects that will be updated with
*               intersections' coordinates. (LLPointPair)
* @param n Pointer to number of intersections found (result. 0, 1, or 2) (double*)
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return Returns error code that indicates success or cause of failure; updates given memory
* addresses with calculated values.
* @retval SUCCESS Indicates successful execution.
* @retval CONCENTRIC_CIRCLE_ERR Indicates that two arcs or circles either do not intersect or are identical. Used as a status code.
* @retval NO_INTERSECTION_ERR Status code indicates that no intersection was found in the case that no intersection point gets returned.
* @retval RADIUS_OUT_OF_RANGE_ERR Indicates that given radius does not meet algorithm requirement.
* @retval SEC_NOT_CONVERGED_ERR Indicates that the secant method failed to converge.
* @retval ERROR_MAX_REACHED_ERR Indicates that an error value has grown larger than allowed.
* @retval UNEXPECTED_ERR Indicates that an unknown error has occurred.
 */
ErrorSet arcIntx(LLPoint center1, double r1, LLPoint center2, double r2,
                 LLPointPair intx, int *n, double tol, double eps);

/* Calculate points of intersection of arc and geodesic.  The arc is treated as a full
* circle and geodesic is treated as unbounded.  Bounds of either object must
* be applied by the parent function, when applicable.
* @param pt1 Starting point of geodesic (LLPoint)
* @param crs1 Azimuth of geodesic at pt1, in radians (double)
* @param center Center of arc (LLPoint)
* @param radius Radius of arc, in nmi (double)
* @param intx Two-element array of LLPoint objects that will be updated with
*               intersections' coordinates. (LLPointPair*)
* @param n Pointer to number of intersections found (result. 0, 1, or 2) (double*)
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return Returns error code that indicates success or cause of failure; updates given memory
* addresses with calculated values.
* @retval SUCCESS Indicates successful execution.
* @retval CONCENTRIC_CIRCLE_ERR Indicates that two arcs or circles either do not intersect or are identical. Used as a status code.
* @retval TOL_TOO_SMALL_ERR Indicates that the requested tolerance cannot be met due to large requested Vincenty algorithm precision.
* @retval ITERATION_MAX_REACHED_ERR Indicates that a method has looped more than the allowed iteration count.
* @retval ERROR_MAX_REACHED_ERR Indicates that an error value has grown larger than allowed.
 */
ErrorSet geoArcIntx(LLPoint pt1, double crs1, LLPoint center,
                    double radius, LLPointPair intx, int *n, double tol,
                    double eps);

/** Calculates the arc of given radius that is tangent to two given geodesics.  The arc is defined by its center point,
* start point, and end point. The star point will be located on the first input geodesic and the end point will be
* located on the second input geodesic.
*
* @param pt1 Position of a point on 1st geodesic (LLPoint)
* @param crs1 Forward azimuth of 1st geodesic at pt1, in radians. (double)
* @param pt3 Position of a point on second geodesic (LLPoint)
* @param crs3 Forward azimuth of 2nd geodesic at pt3, in radians (double)
* @param radius Radius of desired arc, in nmi (double)
* @param centerPoint Pointer to an LLPoint struct that will contain arc's center (LLPoint*)
*                      coordinates
* @param startPoint  Pointer to LLPoint that will contain arc's start
*                     point coordinates (the point where the arc is tangent to the 1st geodesic) (LLPoint*)
* @param endPoint  Pointer to LLPoint that will contain arc's end point
*                   coordinates (the point where the arc is tangent to the 2nd geodesic) (LLPoint*)
* @param dir Pointer to ArcDirection enum that indicates the direction of computed arc (ArcDirection*)
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
*
* @return Returns error code that indicates success or cause of failure; updates given memory
* addresses with calculated values.  If no arc is found, then the latitudes and longitudes of
* the three arc points will be set equal to 0.0.
*
* @remarks This algorithm will fail to meet the required tolerance if the arc tangent
* points are sufficiently close to the middle of the two geodesics (i.e., the points where
* the geodesics reach their maximum separation, where they are nearly "parallel").  If this case
* is encountered, the RADIUS_OUT_OF_RANGE_ERR will be returned because the arc radius will be greater than
* 98% of half the maximum geodesic separation.
*
* @retval SUCCESS Indicates successful execution.
* @retval NO_TANGENT_ARC_ERR Status code indicates that no tangent arc could be found.
* @retval SEC_NOT_CONVERGED_ERR Indicates that the secant method failed to converge.
* @retval ITERATION_MAX_REACHED_ERR Indicates that a method has looped more than the allowed iteration count.
* @retval ERROR_MAX_REACHED_ERR Indicates that an error value has grown larger than allowed.
 */
ErrorSet arcTanToTwoGeos(LLPoint pt1, double crs1, LLPoint pt3,
                         double crs3, double radius,
                         LLPoint *centerPoint, LLPoint *startPoint,
                         LLPoint *endPoint, ArcDirection *dir,
                         double tol, double eps);

/** Constructs an arc beginning at given start point and ending at given end point. The given start course
* defines the geodesic that is tangent to the desired arc at its start point.  Fixing this direction
* forces the arc center point to lie on the same side of the tangent geodesic as the given end point and
* determines the direction of the arc.
*
* This function solves for the arc radius and center point.  It uses a spherical triangle solution to determine an approximate
* initial radius and places a test center point using that radius.  The distance from the test center point
* to the given end point is computed and compared to the the radius.  The radius is adjusted until the centerpoint
* ends up equidistant from the start and end point to within the tol distance.
*
* @image html ArcFromStartAndEnd_illustration.png
*
* @param arcStart Position of the arc start point. (LLPoint)
* @param startCrs Course or azimuth of the arc at the start point.  This is also the azimuth
* of the tangent geodesic at arcStart. (double)
* @param arcEnd Position of the arc end point. (LLPoint)
* @param newArc Pointer to Arc struct that will be updated with the computed arc (Arc*)
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
*
* @return Returns error code that indicates success or cause of failure; updates given memory address
* with the appropriate Arc struct.  If no arc is found, then start/end points of the Arc struct will
* be set to non-physical locations.
*
* @retval SUCCESS Indicates successful execution.
* @retval ITERATION_MAX_REACHED_ERR Indicates that a method has looped more than the allowed iteration count.
* @retval ERROR_MAX_REACHED_ERR Indicates that an error value has grown larger than allowed.
 */
ErrorSet arcFromStartAndEnd(LLPoint arcStart, double startCrs, LLPoint arcEnd, Arc *newArc,
                            double tol, double eps);

/** Constructs an arc beginning at arcStart with tangent course at arcStart equal to arcStartCrs that
* is tangent at its end point to the geodesic with course outCrs at outPoint.  The end point of the arc
* will not necessarily coincide with outPoint.  This algorithm solves for the arc radius, arc end point,
* and arc direction.  Two valid input geometries and their solutions are shown below:
*
* @image html ArcRadiusTangentToOutboundCourse_illustration.png
*
* @param arcStart Position of arc start point. (LLPoint)
* @param arcStartCrs Course or azimuth of the arc at the start point.  This is also the
* azimuth of the tangent geodesic at arcStart. (double)
* @param outPoint Point on tangent outbound geodesic. (LLPoint)
* @param outCrs Course or azimuth of tangent outbound geodesic at outPoint. (double)
* @param newArc Pointer to Arc struct that will be updated with computed arc parameters. (Arc*)
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
*
* @return Returns error code that indicates success or cause of failure; updates given memory address
* with the appropriate Arc struct.  If no arc is found, then radius of the Arc struct will
* be set to zero.
*
* @retval SUCCESS Indicates successful execution.
* @retval RADIUS_OUT_OF_RANGE_ERR Indicates that given radius does not meet algorithm requirement.
* @retval ITERATION_MAX_REACHED_ERR Indicates that a method has looped more than the allowed iteration count.
* @retval ERROR_MAX_REACHED_ERR Indicates that an error value has grown larger than allowed.
 */
ErrorSet arcTanToCrs(LLPoint arcStart, double arcStartCrs, LLPoint outPoint,
                     double outCrs, Arc *newArc, double tol, double eps);

/** Given arc and nearby point, calculate the points on arc where geodesics
* through given point are tangent to arc.  Arc is treated as complete
* circle.  Arc bounds must be applied by parent function, when applicable.
* @param point Given point (LLPoint)
* @param center Center point of arc (LLPoint)
* @param radius Radius of arc in nmi (double)
* @param tanPt Two-element array of LLPoint objects that will be updated with
*                tangent points' coordinates. (LLPointPair)
* @param n Pointer to number of tangent points found (result. 0, 1, or 2) (double*)
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return Returns error code that indicates success or cause of failure; updates the given memory addresses with
* the computed tangent points and number found. If no tangent points are found, then count will be set to zero.
* @retval SUCCESS is returned if no errors are encountered
* @retval NO_TANGENT_ARC_ERR is returned if outPoint is inside the arc
* @retval ITERATION_MAX_REACHED_ERR is returned if convergence is not achieved
* @retval ERROR_MAX_REACHED_ERR is returned if the error measure of an intermediate solution exceeds the limit
* @retval SUCCESS Indicates successful execution.
* @retval NO_TANGENT_ARC_ERR Status code indicates that no tangent arc could be found.
* @retval RADIUS_OUT_OF_RANGE_ERR Indicates that given radius does not meet algorithm requirement.
* @retval ITERATION_MAX_REACHED_ERR Indicates that a method has looped more than the allowed iteration count.
* @retval ERROR_MAX_REACHED_ERR Indicates that an error value has grown larger than allowed.
 */
ErrorSet ptsOnArcOnTanThruPt(LLPoint point, LLPoint center, double radius,
                             LLPointPair tanPt, int *n, double tol, double eps);

/** Given arc and nearby geodesic, find pairs of points on arc and geodesic
* which define two geodesics that are perpendicular to input geodesic and tangent to
* arc.  Arc and geodesic are treated as unbounded.  Bounds must be applied
* by parent function, when applicable.
* @param lineStart Point on geodesic (LLPoint)
* @param crs Azimuth of geodesic at lineStart in radians (double)
* @param center Center point of arc (LLPoint)
* @param radius Radius of arc in nmi (double)
* @param linePts Two element array containing computed points on
*                     geodesic (LLPointPair)
* @param tanPts Two element array containing computed points on
*                    arc (LLPointPair)
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return Returns error code that indicates success or cause of failure; updates the given memory addresses with
* the computed tangent points and number found. If no tangent points are found, then count will be set to zero.
* @retval SUCCESS is returned if no errors are encountered
* @retval INVALID_SHAPE_ERROR is returned if radius is less than the tolerance parameter
* @retval NO_PROJECTED_POINT_ERROR is returned if the arc center cannot be projected onto the input geodesic
* @retval ITERATION_MAX_REACHED_ERR is returned if convergence is not achieved
* @retval ERROR_MAX_REACHED_ERR is returned if the error measure of an intermediate solution exceeds the limit
* @retval SUCCESS Indicates successful execution.
* @retval NO_PROJECTED_POINT_ERR Indicates that an intermediate calculation point based on projecting to a geodesic could not be found
* @retval TOL_TOO_SMALL_ERR Indicates that the requested tolerance cannot be met due to large requested Vincenty algorithm precision.
* @retval INVALID_SHAPE_ERR Indicates that an invaled shape type has been passed.
* @retval ITERATION_MAX_REACHED_ERR Indicates that a method has looped more than the allowed iteration count.
* @retval ERROR_MAX_REACHED_ERR Indicates that an error value has grown larger than allowed.
 */
ErrorSet projectArcTanPtsToGeo(LLPoint lineStart, double crs, LLPoint center,
                               double radius, LLPointPair linePts,
                               LLPointPair tanPts, double tol, double eps);

/** Compute subtended angle of arc. By convention, the arc extent has the same sign as the arc
* orientation. Specifically, clockwise arcs have arc extent < 0.
* @param center Center point of the input arc (LLPoint)
* @param radius Radius of the input arc in nmi (double)
* @param startCrs Azimuth from arc center to arc start point, in radians (double)
* @param endCrs Azimuth from arc center to arc end point, in radians (double)
* @param orientation Direction of computed arc, +1 for counter-clockwise,
*                      -1 for clockwise (ArcDirection)
* @param arcExtent Value of the arc extent in radians (double*)
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param vincentyEps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return Returns error code that indicates success or cause of failure; updates the given memory address with
* the computed arc extent
*
* NOTE: This function treats startCrs and endCrs as if they were of infinite
*       precision.  If you want 2*PI to be returned whenever the start/end
*       points of the arcs are closer than your distance tolerance, then you
*       have to check that before calling this function.
* @retval SUCCESS Indicates successful execution.
* @retval UNEXPECTED_ERR Indicates that an unknown error has occurred.
 */
ErrorSet getArcExtent(LLPoint center, double radius, double startCrs, double endCrs,
                      ArcDirection orientation, double *arcExtent, double tol, double vincentyEps);

/* Test whether a point is on given geodesic.
*
* @param startPt Start point of geodesic (LLPoint)
* @param endPt End point of geodesic (LLPoint)
* @param testPt Point that is to be tested in relation to geodesic (LLPoint)
* @param length Code that indicates bounds of geodesic
*                 (0 = geodesic exists only between start/end,
*                  1 = geodesic extends beyong end
*                  2 = geodesic is not bounded by start/end points) (LineType)
* @param err Pointer to error code.  Will be set to non-zero value on
*              error. (ErrorSet*)
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return Integer indicating true (1) or false (0).
* @retval 1 if point is on geodesic
* @retval 0 if point is not on geodesic
* @retval SUCCESS Indicates successful execution.
 */
int ptIsOnGeo(LLPoint startPt, LLPoint endPt, LLPoint testPt,
              LineType length, ErrorSet *err, double tol, double eps);

/** Test whether point is on given bounded arc.
* @param center Center of arc (LLPoint)
* @param radius Radius of arc in nmi (double)
* @param startCrs Azimuth from center to arc start point in radians (double)
* @param endCrs Azimuth from center to arc end point in radians (double)
* @param orientation Direction of computed arc, +1 for counter-clockwise,
*                      -1 for clockwise (ArcDirection)
* @param testPt Point that is to be tested in relation to arc bounds (LLPoint)
* @param err Pointer to error code.  Will be set to non-zero value on
*              error. (ErrorSet*)
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return Integer indicating true (1) or false (0).
* @retval 1 if point is on arc
* @retval 0 if point is not on arc
* @retval SUCCESS Indicates successful execution.
 */
int ptIsOnArc(LLPoint center, double radius, double startCrs,
              double endCrs, ArcDirection orientation, LLPoint testPt,
              ErrorSet *err, double tol, double eps);

/** Test whether point is inside area bounded by given arc.  For this to return true,
* the test point must be no farther than radius+tol from the arc center point and must
* lie within the angle subtended by the arc.
*
* @param center Center of arc (LLPoint)
* @param radius Radius of arc in nmi (double)
* @param startCrs Azimuth from center to arc start point in radians (double)
* @param endCrs Azimuth from center to arc end point in radians (double)
* @param orientation Direction of computed arc, +1 for counter-clockwise,
*                      -1 for clockwise (ArcDirection)
* @param testPt Point that is to be tested in relation to arc bounds (LLPoint)
* @param err Pointer to error code.  Will be set to non-zero value on
*              error. (ErrorSet*)
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return Integer indicating true (1) or false (0).
* @retval 1 if point is within arc bounds
* @retval 0 if point is not within arc bounds
* @retval SUCCESS Indicates successful execution.
 */
int ptIsInsideArc(LLPoint center, double radius, double startCrs,
                  double endCrs, ArcDirection orientation, LLPoint testPt,
                  ErrorSet *err, double tol, double eps);

/* Compute length of fixed radius arc.  This function discretizes the
* ellipsoidal arc and approximates it by many small spherical arcs (pieces of small circles on the sphere).
* The number of sub-arcs in increased until the path length converges.
* @param center Center point of arc (LLPoint)
* @param radius Radius of arc in nmi (double)
* @param startCrs Azimuth from center to arc start point in radians (double)
* @param endCrs Azimuth from center to arc end point in radians (double)
* @param orient Direction of computed arc, +1 for counter-clockwise,
*                 -1 for clockwise (int)
* @param steps Pointer to number of discretized arcs (int*)
* @param err Pointer to error code.  Will be set to non-zero value on
*              error.*        tol = accuracy tolerance in nmi for intersection calculation
*              (max distance from found intersection to true intersection) (ErrorSet*)
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return Path length of arc in NM.
* @retval SUCCESS Indicates successful execution.
* @retval ITERATION_MAX_REACHED_ERR Indicates that a method has looped more than the allowed iteration count.
* @retval ERROR_MAX_REACHED_ERR Indicates that an error value has grown larger than allowed.
 */
double arcLength(LLPoint center, double radius,
                 double startCrs, double endCrs, ArcDirection orient,
                 int *steps, ErrorSet *err, double tol, double eps);

/** Quickly calculate approximation to arc length.  This algorithm is
* non-iterative, so its accuracy degrades as the arc radius increases.
* @param center Center point of the Arc (LLPoint)
* @param radius Radius of the arc in nmi (double)
* @param startCrs Azimuth from center to arc start point in radians (double)
* @param endCrs Azimuth from center to arc end point in radians (double)
* @param orient Direction of computed arc (int)
* @param meanR Best fit ROC for spherical approximation of ellipsoid (double*)
* @return Returns the approximate arc length.
*/
double approxArcLength(LLPoint center, double radius, double startCrs,
                       double endCrs, ArcDirection orient, double *meanR);

/* Find a point on an arc that is a given distance from another point on the
* arc.  If the given distance is positive, it is measured along the arc in the same
* direction as the arc orientation.  If the distance is negative, then the new point is
* computed by moving opposite the arc orientation.
* @param center Center point of arc (LLPoint)
* @param givenPoint Point from which distance is to be measured (LLPoint)
* @param dir Direction or orientation of arc (ArcDirection)
* @param arcDist Distance from givenPoint (in nmi).  Use negative value
*                  if you want to measure in direction opposite to arc's
*                  orientation. (double)
* @param endPoint Calculated point (LLPoint*)
* @param subtendedAngle Signed difference in azimuth (double*)
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return Returns error code that indicates success or cause of failure; updates the given memory address with
* 				the calculated values
* @retval SUCCESS Indicates successful execution.
* @retval RADIUS_OUT_OF_RANGE_ERR Indicates that given radius does not meet algorithm requirement.
* @retval ITERATION_MAX_REACHED_ERR Indicates that a method has looped more than the allowed iteration count.
* @retval ERROR_MAX_REACHED_ERR Indicates that an error value has grown larger than allowed.
 */
ErrorSet arcFromLength(LLPoint center, LLPoint givenPoint, ArcDirection dir, double arcDist,
                       LLPoint *endPoint, double *subtendedAngle, double tol, double eps);

/** Calculate distance from defining geodesic to locus at point on the geodesic
* that is given distance from geodesic's start point
* @param loc Locus structure (Locus)
* @param distance Distance in nmi from geodesic start point to point
*                   of interest (double)
* @return distance from geodesic to locus in NM. Distance > 0 if locus is to right of geodesic
*               distance < 0 if locus is to left of geodesic
***/
double distToLocusFromGeoDist(Locus loc, double distance);

/* Calculate distance from defining geodesic to locus at a given point on the geodesic. The
* boundedness of the locus specified by the loc.lineType attribute is ignored by this function so
* it effectively treats the input locus as unbounded.
* @param loc Locus structure (Locus)
* @param geopt Point on locus's defining geodesic (LLPoint)
* @param faz Pointer to forward azimuth of geodesic at geopt (double*)
* @param err Pointer to error code.  Set to non-zero value on error. (ErrorSet*)
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return Signed distance to locus; err is set to 0.
* @retval SUCCESS Indicates successful execution.
* @retval RADIUS_OUT_OF_RANGE_ERR Indicates that given radius does not meet algorithm requirement.
 */
double distToLocusFromGeoPt(Locus loc, LLPoint geopt, double *faz, ErrorSet *err,
                            double tol, double eps);

/** Determine point on locus abeam given point on defining geodesic. The
* boundedness of the locus specified by the loc.lineType attribute is ignored by this function so
* it effectively treats the input locus as unbounded.
* @param loc Locus structure (Locus)
* @param geoPt Point on locus's defining geodesic (LLPoint)
* @param ptOnLoc Pointer to LLPoint, updated with coordinates of point
*                  on locus abeam given point. (LLPoint*)
* @param perpCrs Pointer to double, updated with azimuth from point
*                  on geodesic to point on locus (double*)
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return an error code that indicates success or failure; on success, ptOnLoc is updated
* with the locus point's coordinates.
* @retval SUCCESS Indicates successful execution.
* @retval RADIUS_OUT_OF_RANGE_ERR Indicates that given radius does not meet algorithm requirement.
 */
ErrorSet ptOnLocusFromGeoPt(Locus loc, LLPoint geoPt, LLPoint *ptOnLoc,
                            double *perpCrs, double tol, double eps);

/** Determine whether a given point is on a locus.  If point is on Locus,
* then pointer to projected point on geodesic is returned; otherwise, NULL is
* returned.  Returning pointer to projected point allows use of that point in
* other calculations, possibly saving steps.
* @param loc Locus structure (Locus)
* @param testPt Point on defining geodesic (LLPoint)
* @param ptOnGeo Pointer to LLPoint, updated with point on defining geodesic
*              abeam the given point on the locus (LLPoint*)
* @param err Pointer to error code.  Set to non-zero value on error. (ErrorSet*)
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return Integer indicating true (1) or false (0).
* @retval 1 if point is on arc
* @retval 0 if point is not on arc
* @retval SUCCESS Indicates successful execution.
 */
int ptIsOnLocus(Locus loc, LLPoint testPt, LLPoint *ptOnGeo, ErrorSet *err,
                double tol, double eps);

/** Creates a locus from an input geodesic and the locus start and end points.
* @param loc Pointer to a locus structure to be filled by the function (Locus*)
* @param geoStart Start point of the input geodesic (LLPoint)
* @param geoStartAz Start azimuth of the input geodesic (double)
* @param locStart Start point of the output locus (LLPoint)
* @param locEnd End point of the output locus (LLPoint)
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return Returns error code that indicates success or cause of failure; updates the given memory address with
* 				the constructed locus
* @retval SUCCESS Indicates successful execution.
 */
ErrorSet locusFromGeoAndPts(Locus *loc, LLPoint geoStart, double geoStartAz,
                            LLPoint locStart, LLPoint locEnd, double tol, double eps);

/** Find points where arc intersects with locus.  Arc is treated as complete
* circle, however locus bounds are checked and applied.  Arc bounds must be
* applied by parent function, when applicable.
* @param loc Locus structure (Locus)
* @param center Center of arc (LLPoint)
* @param radius Radius of arc, in nmi (double)
* @param intx Two element array of LLPoints that will be updated with
*               intersection coordinates. (LLPointPair)
* @param n Pointer to number of intersections found (result. 0, 1, or 2) (int*)
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return Returns error code that indicates success or failure; on success, intx is updated
* with n intersection coordinate pairs.
* @retval SUCCESS Indicates successful execution.
* @retval SEC_NOT_CONVERGED_ERR Indicates that the secant method failed to converge.
* @retval ERROR_MAX_REACHED_ERR Indicates that an error value has grown larger than allowed.
 */
ErrorSet locusArcIntx(Locus loc, LLPoint center, double radius,
                      LLPointPair intx, int *n, double tol, double eps);

/** Find points where geodesic and locus intersect.  Bounds of locus are
* applied, bounds of geodesic are not.
* KNOWN LIMITATION:  Only one intersection will be returned.  It is possible for
* nearly-parallel locus and geodesic to intersect twice, but this algorithm will
* only return (at most) one intersection.
* @param geost Start point of geodesic (LLPoint)
* @param geoend End point of geodesic (LLPoint)
* @param loc Locus structure (Locus)
* @param intx Pointer to LLPoint that will be updated with intersection
*               coordinates. (LLPoint*)
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return Returns error code that indicates success or failure; on success, intx is updated
* if an intersection is found.
* @retval SUCCESS Indicates successful execution.
* @retval NO_INTERSECTION_ERR Status code indicates that no intersection was found in the case that no intersection point gets returned.
* @retval RADIUS_OUT_OF_RANGE_ERR Indicates that given radius does not meet algorithm requirement.
* @retval INVALID_SHAPE_ERR Indicates that an invaled shape type has been passed.
* @retval SEC_NOT_CONVERGED_ERR Indicates that the secant method failed to converge.
* @retval ERROR_MAX_REACHED_ERR Indicates that an error value has grown larger than allowed.
* @retval UNEXPECTED_ERR Indicates that an unknown error has occurred.
 */
ErrorSet locusGeoIntx(LLPoint geost, LLPoint geoend, Locus loc,
                      LLPoint *intx, double tol, double eps);

/* Find points where two loci intersect.  The bounds of both loci are applied.
* KNOWN LIMITATION:  Only one intersection will be returned.  It is possible for
* nearly-parallel loci to intersect twice, but this algorithm will only return
* (at most) one intersection.
* @param loc1 First input locus (Locus)
* @param loc2 Second input locus (Locus)
* @param intx Pointer to LLPoint that will be updated with intersection
*               coordinates. (LLPoint*)
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return an error code that indicates success or failure; on success, intx is updated
* if an intersection is found.
* @retval SUCCESS Indicates successful execution.
* @retval NO_INTERSECTION_ERR Status code indicates that no intersection was found in the case that no intersection point gets returned.
* @retval INVALID_TYPE_ERR Indicates that the incorrect variable type has been passed.
* @retval SHAPE_NOT_DEFINED_ERR Indicates that a shape object has not been defined.
* @retval SEC_NOT_CONVERGED_ERR Indicates that the secant method failed to converge.
* @retval ERROR_MAX_REACHED_ERR Indicates that an error value has grown larger than allowed.
 */
ErrorSet locusIntx(Locus loc1, Locus loc2, LLPoint *intx, double tol,
                   double eps);

/** Determine arc of given radius that is tangent to two given loci.
* @param loc1 First input locus (Locus)
* @param loc2 Second input locus (Locus)
* @param radius Radius of desired arc in nmi. (double)
* @param centerPoint Pointer to LLPoint that will contain arc's center
*              coordinates. (LLPoint*)
* @param startPoint Pointer to LLPoint that will contain arc's start
*              point coordinates. (LLPoint*)
* @param endPoint Pointer to LLPoint that will contain arc's end point
*              coordinates. (LLPoint*)
* @param dir Pointer to integer direction of computed arc. Enum type,
*              valid values are CLOCKWISE and COUNTERCLOCKWISE. (ArcDirection*)
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return Returns error code that indicates success or failure; on success, the arc center, start, and end
* 				points are updated with the appropriate coordinates.
* @retval SUCCESS Indicates successful execution.
* @retval INVALID_CRS_ERR Indicates that an azimuth value was out of range (e.g., desired course from pole other than \f$\pm\pi\f$).
* @retval ITERATION_MAX_REACHED_ERR Indicates that a method has looped more than the allowed iteration count.
* @retval ERROR_MAX_REACHED_ERR Indicates that an error value has grown larger than allowed.
 */
ErrorSet arcTanToTwoLoci(Locus loc1, Locus loc2, double radius,
                         LLPoint *centerPoint, LLPoint *startPoint,
                         LLPoint *endPoint, ArcDirection *dir,
                         double tol, double eps);

/** Given geodesic defined by startPt and endPt, along with testPt on geodesic,
* return course of geodesic at testPt.
* @param geo Geodesic struct containing the geodesic parameters (Geodesic)
* @param testPt Point at which azimuth of geo is to be computed (LLPoint)
* @param startCrs Azimuth of geodesic from start point to end point (double*)
* @param revCrs Azimuth of geodesic from end point to start point (double*)
* @param distToPt Distance from startPt to testPt in nmi (double*)
* @param err Pointer to error code. Set to non-zero value on error. (ErrorSet*)
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return Returns azimuth in radians of geodesic at test point
* @retval SUCCESS Indicates successful execution.
* @retval UNEXPECTED_ERR Indicates that an unknown error has occurred.
 */
double geoCrs(Geodesic geo, LLPoint testPt,
              double *startCrs, double *revCrs,
              double *distToPt, ErrorSet *err, double tol,
              double eps);

/** Given locus structure and testPt on locus, return azimuth of locus at
* testPt.  pt must lie on locus.
* @param loc Locus structure (Locus)
* @param testPt Point at which azimuth of geodesic is to be determined (LLPoint)
* @param geoPt Point on locus's defining geodesic corresponding to testPt (LLPoint)
* @param perpCrs Azimuth from testPt to geoPt, in radians (double*)
* @param err Pointer to integer error code.  Set to non-zero value on error. (ErrorSet*)
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return azimuth in radians of locus at test point.
* @retval SUCCESS Indicates successful execution.
* @retval POINT_NOT_ON_LINE_ERR Status code indicates that a point was not on the expected geodesic or locus
 */
double locusCrsAtPt(Locus loc, LLPoint testPt, LLPoint *geoPt,
                    double *perpCrs, ErrorSet *err, double tol, double eps);

/** Finds the course on the locus at a point abeam the geodesic at a distance from the geodesic start point
* @param loc Input locus on which to find the course (Locus)
* @param distFromGeoStart Distance from the geodesic start point at which to calculate the locus course (double)
* @param locPt Point on the locus at which to find the course (LLPoint*)
* @param perpCrs Course perpendicular to the geodesic passing through locPt (double*)
* @param err Errors returned during calculation (ErrorSet*)
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return Returns locus course abeam the geodesic at the input distance
* @retval SUCCESS Indicates successful execution.
 */
double locusCrsAtGeoDist(Locus loc, double distFromGeoStart, LLPoint *locPt,
                         double *perpCrs, ErrorSet *err, double tol, double eps);

/** Calculate point on locus nearest a given point.  In other words, project
* the given point to the locus
* @param loc Locus structure (Locus)
* @param pt3 Point to be projected to geodesic (LLPoint)
* @param perpPt Pointer to LLPoint that will be updated with coordinates of
*                 projected point. (LLPoint*)
* @param crsFromPoint Pointer to azimuth of geodesic from pt3 to projected
*                        point, in radians (double*)
* @param distFromPoint Pointer to distance from pt3 to locus at projected
*                         point, in nmi (double*)
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return Returns error code that indicates success or failure; on success, perpPt is updated
* with the appropriate coordinates.
* @retval SUCCESS Indicates successful execution.
* @retval NO_PROJECTED_POINT_ERR Indicates that an intermediate calculation point based on projecting to a geodesic could not be found
* @retval ITERATION_MAX_REACHED_ERR Indicates that a method has looped more than the allowed iteration count.
* @retval ERROR_MAX_REACHED_ERR Indicates that an error value has grown larger than allowed.
 */
ErrorSet projectToLocus(Locus loc, LLPoint pt3, LLPoint *perpPt,
                        double *crsFromPoint, double *distFromPoint,
                        double tol, double eps);

/* Given two circles, find geodesic that connects tangent points on both circles
* @param center1 Center point of the first circle (LLPoint)
* @param r1 Radius of the first circle (double)
* @param dir1 Turn direction of the first circle (ArcDirection)
* @param center2 Center point of the second circle (LLPoint)
* @param r2 Radius of the second circle (double)
* @param dir2 Turn direction of the second circle (ArcDirection)
* @param tanLines Two geodesics tangent to the input arcs (Geodesic)
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return Returns error code that indicates success or cause of failure; updates the given memory address with
* 				the calculated geodesics.
* @retval SUCCESS Indicates successful execution.
* @retval CIRCLE_INSIDE_CIRCLE_ERR Status code indicates that no arc-arc intersection was found because one arc lies entirely inside the other.
* @retval ITERATION_MAX_REACHED_ERR Indicates that a method has looped more than the allowed iteration count.
* @retval ERROR_MAX_REACHED_ERR Indicates that an error value has grown larger than allowed.
 */
ErrorSet geoTanToTwoCircles(LLPoint center1, double r1, ArcDirection dir1,
                            LLPoint center2, double r2, ArcDirection dir2,
                            Geodesic tanLines[2], double tol, double eps);

/* The following functions reside in Spiral.c.
*  */
/** Finds the radius of the input spiral at the input azimuth
** @param sp Input spiral (Spiral)
** @param az Azimuth at which to find the radius of sp (double)
** @param radius Spiral radius, updated and returned after computation (double)
** @return Returns error code that indicates success or cause of failure; updates given memory
** 				address with calculated radius.
* @retval SUCCESS Indicates successful execution.
* @retval NO_MEMORY_ALLOCATED_ERR Indicates that a NULL pointer was passed for a required reference
 */
ErrorSet spiralRadius(Spiral sp, double az, double *radius);

/** Finds the point on the input spiral at the input azimuth
* @param sp Input spiral (Spiral)
* @param az Azimuth at which to find the point on sp (double)
* @param pt Pointer to LLPoint structure updated by function (LLPoint)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return Returns error code that indicates success or cause of failure; updates given memory
* 				address with calculated point.
* @retval SUCCESS is returned if no error is encountered
* ERROR RESULT: Error code > 0 is returned at err address.  Return value is undefined.
* @retval SUCCESS Indicates successful execution.
 */
ErrorSet ptOnSpiral(Spiral sp, double az, LLPoint *pt, double eps);

/* Determines if the input point lies on the input spiral
* @param sp Input spiral (Spiral)
* @param point Input point (LLPoint)
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return returns 1 if the spiral contains the point, 0 if not */
int ptIsOnSpiral(Spiral sp, LLPoint point, double tol, double eps);

/* Truncates the input spiral's extent such that the new end point is the input point
* @param sp Input spiral
* @param pt New end point of the input spiral
* @param newSp Pointer to Spiral structure to be updated by function with new spiral
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return Returns error code that indicates success or cause of failure; updates given memory
* 				address with new spiral.
* @retval SUCCESS Indicates successful execution.
* @retval NO_MEMORY_ALLOCATED_ERR Indicates that a NULL pointer was passed for a required reference
* @retval POINT_NOT_ON_ARC_ERR Indicates that start/end points used to define an arc are not equidistant from arc center.
 */
ErrorSet moveSpiralEndToPt(Spiral sp, LLPoint pt, Spiral *newSp, double tol, double eps);

/* Truncates the input spiral's extent such that the new start point is the input point
* @param sp Input spiral (Spiral)
* @param pt New start point of the input spiral (LLPoint)
* @param newSp Pointer to Spiral structure to be updated by function with new spiral (Spiral*)
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return Returns error code that indicates success or cause of failure; updates given memory
* 				address with new spiral.
* @retval SUCCESS Indicates successful execution.
* @retval NO_MEMORY_ALLOCATED_ERR Indicates that a NULL pointer was passed for a required reference
* @retval POINT_NOT_ON_ARC_ERR Indicates that start/end points used to define an arc are not equidistant from arc center.
 */
ErrorSet moveSpiralStartToPt(Spiral sp, LLPoint pt, Spiral *newSp, double tol, double eps);

/** Finds the intersection points of the input spiral and input line
* @param sp Input spiral (Spiral)
* @param line Input geodesic (Geodesic)
* @param pts Array to hold intersections points (LLPointSet*)
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return Returns error code that indicates success or cause of failure; updates given memory
* 				address with intersection points.
* @retval SUCCESS Indicates successful execution.
 */
ErrorSet spiralGeoIntx(Spiral sp, Geodesic line, LLPointSet *pts, double tol, double eps);

/** Finds the intersection points of the input spiral and input locus
* @param sp Input spiral (Spiral)
* @param locus Input locus (Locus)
* @param pts Array to hold intersections points (LLPointSet*)
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return Returns error code that indicates success or cause of failure; updates given memory
* 				address with intersection points.
* @retval SUCCESS Indicates successful execution.
* @retval NO_INTERSECTION_ERR Status code indicates that no intersection was found in the case that no intersection point gets returned.
 */
ErrorSet spiralLocusIntx(Spiral sp, Locus locus, LLPointSet *pts, double tol, double eps);

/** Finds the intersection points of the input spiral and input arc
* @param sp Input spiral (Spiral)
* @param arc Input arc (Arc)
* @param pts Array to hold intersections points (LLPointSet*)
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return Returns error code that indicates success or cause of failure; updates given memory
* 				address with intersection points.
* ERROR RESULT: Error code > 0 is returned at err address.  Return value is undefined.
* @retval SUCCESS Indicates successful execution.
 */
ErrorSet spiralArcIntx(Spiral sp, Arc arc, LLPointSet *pts, double tol, double eps);

/** Finds the intersection points of the input spiral and a second input spiral
* @param sp1 First Input Spiral (Spiral)
* @param sp2 Second Input Spiral (Spiral)
* @param pts Array to hold intersections points (LLPointSet)
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return Returns error code that indicates success or cause of failure; updates given memory
* 				address with intersection points.
* @retval SUCCESS Indicates successful execution.
 */
ErrorSet spiralIntx(Spiral sp1, Spiral sp2, LLPointSet *pts, double tol, double eps);

/** Finds the points on the input spiral that create tangent lines when paired with input point
* @param sp Input spiral (Spiral)
* @param pt Input point (LLPoint)
* @param pts Array to hold tangent points (LLPointSet*)
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return Returns error code that indicates success or cause of failure; updates given memory
* 				address with tangent points.
* @retval SUCCESS Indicates successful execution.
* @retval NO_INTERSECTION_ERR Status code indicates that no intersection was found in the case that no intersection point gets returned.
 */
ErrorSet ptsOnSpiralOnTanThruPt(Spiral sp, LLPoint pt, LLPointSet *pts, double tol, double eps);

/** Finds the lines tangent to both input spirals
* @param sp1 First input spiral (Spiral)
* @param sp2 Second input spiral (Spiral)
* @param pts Array to hold tangent points (LLPointSet*)
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return Returns error code that indicates success or cause of failure; updates given memory
* 				address with tangent points.
* @retval SUCCESS Indicates successful execution.
 */
ErrorSet geoTanToTwoSpirals(Spiral sp1, Spiral sp2, LLPointSet *pts, double tol, double eps);

/** Finds the perpendicular projection of the input point on the input spiral
* @param sp Input spiral (Spiral)
* @param pt Input point (LLPoint)
* @param perpPt Perpendicular projection of point onto spiral (LLPoint*)
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return Returns error code that indicates success or cause of failure; updates given memory
* 				address with projected point.
* @retval SUCCESS Indicates successful execution.
* @retval NO_INTERSECTION_ERR Status code indicates that no intersection was found in the case that no intersection point gets returned.
 */
ErrorSet projectToSpiral(Spiral sp, LLPoint pt, LLPoint *perpPt, double tol, double eps);

/** Finds the locus that is an angle from the input geodesic that is tangential to the input spiral
* @param sp Input spiral (Spiral)
* @param geo Input geodesic (Geodesic)
* @param angle Course difference between the line and the locus (double)
* @param tanGeo Output geodesic (Geodesic*)
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return Returns error code that indicates success or cause of failure; updates given memory
* 				address with output geodesic.
* @retval SUCCESS Indicates successful execution.
 */
ErrorSet geoTanToSpiralAtAngleToGeo(Spiral sp, Geodesic geo, double angle, Geodesic *tanGeo, double tol, double eps);

/** Finds the points on the input spiral at which the tangent course at the point is equal to the course of the line at a point perpendicular to the point on the spiral
* @param sp Input spiral (Spiral)
* @param line Input line (Geodesic)
* @param pts Output points (LLPointSet*)
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return Returns error code that indicates success or cause of failure; updates given memory
* 				address with output locus.
* @retval SUCCESS Indicates successful execution.
 */
ErrorSet spiralTanPts(Spiral sp, Geodesic line, LLPointSet *pts, double tol, double eps);

/** Returns an integer that indicates whether a testPt is at a pole and which one.
* @param testPt Point in question (LLPoint)
* @param err Pointer used to update the error code value (ErrorSet*)
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return -1 if at South Pole, 0 if NOT at a pole, 1 it at North Pole
* @retval SUCCESS Indicates successful execution.
* @retval ANTIPODAL_POINTS_ERR Indicates that two LLPoints are on opposite sides of the earth ellipsoid.
 */
int ptIsAtPole(LLPoint testPt, ErrorSet *err, double tol, double eps);

/*
* Arc Functions
*/

/**
*
* Find the geodesic that is tangent to the given input geodesic and arc
* such that the tangent geodesic forms the given input angle with respect
* to the input geodesic.
*
* Interior angle is formed by the given geodesic and the found
* tangent line.
*
* NOTE - Up to four distinct solutions exists for any given input angle,
* 		  geodesic and arc.  It is possible for no solution to exist as well.
* 		  In addition to checking the error code after this function is called,
* 		  the tangentLineLocation argument should also be checked to see if a solution exists.
*
* 		  Furthermore, if a solution exists, only one of the distinct tangent line
* 		  solutions is returned.  If the returned tangent line is not of the desired
* 		  orientation, changing either the direction of the input geodesic and/or
* 		  the orientation of the input arc will result in a different tangent line
*        from the original solution set to be returned.  Therefore,
*        running through the four combinations of input geo directions and arc orientations
*        will generate the entire solution set of tangent lines.
*
*
* @param arc The input arc to which we want to find a tangential line to (Arc)
* @param geo The input geodesic to which we want to find a tangential line to (Geodesic)
* @param angle The inner angle, between 0 <= angle <= pi/2 radians, that should be formed by the
* 				  tangential line with respect to the input geodesic. (double)
* @param tangentLine pointer for the output tangential line, if it exists. (Geodesic*)
* @param tangentLineLocation A pointer for the location of the tangent line with respect to the input geo:
* 			                    -1 = tangent line strictly to the left of the geodesic, 0 = no solution exists,
*                              1 = tangent line strictly to the right the geodesic (int*)
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @returns err - Error code zero indicates the function successfully completed.
* 			      Otherwise an error occurred.
* @retval SUCCESS Indicates successful execution.
* @retval ITERATION_MAX_REACHED_ERR Indicates that a method has looped more than the allowed iteration count.
 */
ErrorSet geoTanToArcAtAngleToGeo(Arc arc,
                                 Geodesic geo,
                                 double angle,
                                 Geodesic *tangentLine,
                                 int *tangentLineLocation,
                                 double tol,
                                 double eps);
/* 
*
*
*                lineEndPoint
*                    x
*                    |\
*       input line   | \  Line tangent to arc
*               ---> |  \  <---
*                    |   \
*                    |    \
*                    |   ^ \
*                    |   |  \
*                    x angle \     *   *
*       lineStartPoint        \ *         *
*         tangent point to arc x           *
*                               \*       *
*                                \  *  *
*                                 \
*
*/

/**
*
* Find the arc that is tangent to the given input geodesic and arc
* such that the tangent arc is of the desired radius.
*
*
*
* NOTE - Up to four distinct solutions exists for any given input
* 		  geodesic and arc.  It is possible for no solution to exist as well.
* 		  In addition to checking the error code after this function is called,
* 		  the tangentArcLocation argument should also be checked to see if a solution exists.
*
* 		  Furthermore, if a solution exists, only one of the distinct tangent arc
* 		  solutions is returned.  If the returned tangent arc is not of the desired
* 		  orientation, changing either the direction of the input geodesic and/or
* 		  the orientation of the input arc will result in a different tangent arc
*        from the original solution set to be returned.  Therefore,
*        running through the four combinations of input geo directions and arc orientations
*        will generate the entire solution set of tangent arcs.
*
*
* @param arc The input arc to which we want to find a tangential arc to (Arc)
* @param geo The input geodesic to which we want to find a tangential arc to (Geodesic)
* @param radius The desired radius of the tangent arc in nautical miles (double)
* @param tangentArc A pointer for the output tangential arc, if it exists. (Arc*)
* @param tangentArcLocation A pointer for the location of the tangent arc with respect to the input geo:
* 			                    -1 = tangent arc strictly to the left of the geodesic, 0 = no solution exists,
*                              1 = tangent arc strictly to the right the geodesic (int*)
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return Returns error code that indicates success or cause of failure; updates given memory
* 				address with output arc and arc status.
* @retval SUCCESS Indicates successful execution.
* @retval ITERATION_MAX_REACHED_ERR Indicates that a method has looped more than the allowed iteration count.
* @retval ERROR_MAX_REACHED_ERR Indicates that an error value has grown larger than allowed.
* @retval UNEXPECTED_ERR Indicates that an unknown error has occurred.
 */
ErrorSet arcTanToArcAndGeo(Arc arc,
                           Geodesic geo,
                           double radius,
                           Arc *tangentArc,
                           int *tangentArcLocation,
                           double tol,
                           double eps);

/**
* Construct arc from start point, start course, and radius.  The arc's
* end point is computed such that the course from the end point to the
* given nextEndPoint is tangent to the arc at the end point.
* @param arcStart Start point of the arc (LLPoint)
* @param startCrs Course of the arc at the start point of the arc (double)
* @param radius Radius of the arc (double)
* @param arcDir Turn direction of the arc (ArcDirection)
* @param nextEndPoint Point through which the tangent line to the arc at
*     the arc's end point passes through (LLPoint)
* @param newArc Pointer to an arc structure to hold the output arc (Arc*)
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return Returns error code that indicates success or cause of failure; updates given memory
* 				address with output arc.
* @retval SUCCESS Indicates successful execution.
 */
ErrorSet arcEndFromStartAndRadius(LLPoint arcStart, double startCrs, double radius, ArcDirection arcDir,
                                  LLPoint nextEndPoint, Arc *newArc, double tol, double eps);

/*************************************************************************
* Construct arc from start point, start course, and center.  The arc's
* end point is computed such that the course from the end point to the
* given nextEndPoint is tangent to the arc at the end point.'
* @param arcStart Start point of the arc (LLPoint)
* @param startCrs Course of the arc at the start point of the arc (double)
* @param arcCenter Center point of the arc (LLPoint)
* @param arcDir Turn direction of the arc (ArcDirection)
* @param nextEndPoint Point through which the tangent line to the arc at
*     the arc's end point passes through (LLPoint)
* @param newArc Pointer to an arc structure to hold the output arc (Arc*)
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return Returns error code that indicates success or cause of failure; updates given memory
* 				address with output arc.
* @retval SUCCESS Indicates successful execution.
 */
ErrorSet arcEndFromStartAndCenter(LLPoint arcStart, double startCrs, LLPoint arcCenter,
                                  ArcDirection arcDir, LLPoint nextEndPoint, Arc *newArc,
                                  double tol, double eps);

/**
* Test whether a given locus lies on other given locus
*
* Note - It is assumed that the geodesics used to define the input loci are
* 		  well-defined, i.e. the geodesic start/end points are not identical
* 		  (with respect to tol) for a given locus
*
* @param locus1 First input locus (Locus)
* @param locus2 Second input locus (Locus)
* @param commonShape If locus1 and locus2 coincide then either the common locus segment or the common point is returned (Shape*)
* @param err The pointer to the error variable (ErrorSet*)
* @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
* @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
* @return Returns 0 if loci do not coincide or if error resulted.  Return 1 otherwise.
* @retval SUCCESS Indicates successful execution.
* @retval NO_MEMORY_ALLOCATED_ERR Indicates that a NULL pointer was passed for a required reference
* @retval INVALID_SHAPE_ERR Indicates that an invaled shape type has been passed.
* @retval UNEXPECTED_ERR Indicates that an unknown error has occurred.
 */
int lociCoincide(Locus locus1, Locus locus2, Shape *commonShape, ErrorSet *err, double tol, double eps);

/**
*  Find distance and azimuth of great circle between two points on spherical
* earth model
* @param org Start point of great circle (LLPoint)
* @param dest End point of great circle (LLPoint)
* @param crs Pointer to double that will be updated with the azimuth of great circle at org in radians (double*)
* @param dist Point to double that will be updated with the distance between org and dest in nmi (double*)
* @param userROC Pointer to double used to pass desired sphere radius
*                   (NULL is valid and revert to best fit or default) (double*)
* @param eps Convergence parameter for Vincenty forward/inverse algorithms (double)
* @return Updates given memory addresses with calculated course and distance.
* @retval Nothing
*/
void sphereInverse(LLPoint org, LLPoint dest, double *crs, double *dist,
                   double *userROC, double eps);

/**
* Find distance between points on great circle (using spherical earth model).
* @param org Start point of great circle (LLPoint)
* @param dest End point of great circle (LLPoint)
* @param userROC Pointer to double used to pass desired sphere radius
*                   (NULL is valid and revert to best fit or default) (double*)
* @return Returns the calculated spherical distance.
*/
double sphereInvDist(LLPoint org, LLPoint dest, double *userROC);

/* Find azimuth of great circle between two points on spherical
* earth model
* @param org Start point of great circle (LLPoint)
* @param dest End point of great circle (LLPoint)
* @param eps Convergence parameter for Vincenty forward/inverse algorithms (double)
* @return Returns the calculated course on the earth sphere.
*/
double sphereInvCrs(LLPoint org, LLPoint dest, double eps);

/* Calculate a spherical position (lat/lon in radians) given starting position,
* course and distance
* @param org Start point of great circle (LLPoint)
* @param course Course to destination point (double)
* @param dist Distance to destination point (double)
* @param dest Pointer to LLPoint struct to be updated by algorithm
* @return Returns error code that indicates success or cause of failure; updates given memory
* 				address with destination point.
* @retval SUCCESS Indicates successful execution.
 */
ErrorSet sphereDirect(LLPoint org, double course, double dist,
                      LLPoint *dest);

/* Find the intersection of two fixed-radius arcs on spherical earth model.
* Bounds of arcs are not applied.
* @param center1 Center of the first arc (LLPoint)
* @param r1 Radius of the first arc (double)
* @param center2 Center of the second arc (LLPoint)
* @param r2 Radius of the second arc (double)
* @param intx LLPointPair to hold intersection points (LLPointPair)
* @param n Number of intersection points found (int*)
* @param bestFitROC Pointer to double with value used as
* 			sphere radius within function (double*)
* @param eps convergence parameter for Vincenty forward/inverse algorithms (double)
* @return Returns error code that indicates success or cause of failure; updates given memory
* 				addresses with number of intersections and calculated points.
* @retval SUCCESS Indicates successful execution.
* @retval NO_INTERSECTION_ERR Status code indicates that no intersection was found in the case that no intersection point gets returned.
* @retval RADIUS_OUT_OF_RANGE_ERR Indicates that given radius does not meet algorithm requirement.
* @retval UNEXPECTED_ERR Indicates that an unknown error has occurred.
 */
ErrorSet initArcIntx(LLPoint center1, double r1, LLPoint center2, double r2,
                     LLPointPair intx, int *n, double *bestFitROC,
                     double eps);

#define BOUNDARY_ARRAY_INCREMENT 20

/** Initialize and return a Boundary object.
* @return Returns a Boundary with length = 0, elements pointer initialized and
* 				pointing to unitialized memory sufficiently large to hold BOUNDARY_ARRAY_INCREMENT Shape pointers. (Boundary)
*/
Boundary createBndry();

/** Deallocate memory reserved for Boundary struct. Memory pointed to by given pointer is freed.
* @param b Pointer to Boundary struct to be deallocated (Boundary*)
* @return Updates given memory address with cleared boundary struct.
* @retval Nothing */
void clearBndry(Boundary *b);

/** Add element to Boundary struct
* @param b Pointer to Boundary struct (Boundary*)
* @param element Pointer to element to be added (const void*)
* @param type of element being added (ShapeType)
* @return Returns error code that indicates success or failure (ErrorSet); updates given memory for Boundary with added element
* @retval SUCCESS Indicates successful execution.
* @retval INVALID_TYPE_ERR Indicates that the incorrect variable type has been passed.
 */
ErrorSet addElementToBndry(Boundary *b, const void *element, ShapeType type);

/** Add Geodesic to Boundary struct
* @param b Pointer to Boundary struct (Boundary*)
* @param g Pointer to Geodesic to be added (Geodesic*)
* @return Returns error code that indicates success or failure (ErrorSet); updates given memory for Boundary with added Geodesic
* @retval SUCCESS Indicates successful execution.
 */
ErrorSet addGeoToBndry(Boundary *b, const Geodesic *g);

/** Add Arc to Boundary struct
* @param b Pointer to Boundary struct (Boundary*)
* @param a Pointer to Arc to be added (Arc*)
* @return Returns error code that indicates success or failure (ErrorSet); updates given memory for Boundary with added Arc
* @retval SUCCESS Indicates successful execution.
 */
ErrorSet addArcToBndry(Boundary *b, const Arc *a);

/** Add Locus to Boundary struct
* @param b Pointer to Boundary struct (Boundary*)
* @param l Pointer to Locus to be added (Locus*)
* @return Returns error code that indicates success or failure (ErrorSet); updates given memory for Boundary with added Locus
* @retval SUCCESS Indicates successful execution.
 */
ErrorSet addLocusToBndry(Boundary *b, const Locus *l);

/** Add a Locus2ndOrder to Boundary struct
* @param b Pointer to Boundary struct (Boundary*)
* @param l Pointer to Locus2ndOrder to be added (Locus2ndOrder*)
* @return Returns error code that indicates success or failure (ErrorSet); updates given memory for Boundary with added Locus2ndOrder
* @retval SUCCESS Indicates successful execution.
 */
ErrorSet addLocus2ndOrderToBndry(Boundary *b, const Locus2ndOrder *l);

/** Add Spiral to Boundary struct
* @param b Pointer to Boundary struct (Boundary*)
* @param s Pointer to Spiral to be added (Spiral*)
* @return Returns error code that indicates success or failure (ErrorSet); updates given memory for Boundary with added Spiral
* @retval SUCCESS Indicates successful execution.
 */
ErrorSet addSpiralToBndry(Boundary *b, const Spiral *s);

/** Determines whether a given point is inside a Boundary
* @param b Boundary (Boundary)
* @param p point to be tested (LLPoint)
* @param err Error code that indicates success or cause of failure (ErrorSet*)
* @param tol Tolerance for function calculations (double)
* @param eps Convergence tolerance for Vincenty algorithm (double)
* @return Returns 1 if the point is in or on the Boundary and 0 if it is outside; updates given memory address with error code
* @retval SUCCESS Indicates successful execution.
* @retval SHAPE_NOT_DEFINED_ERR Indicates that a shape object has not been defined.
 */
int ptIsInsideBndry(Boundary b, LLPoint p, ErrorSet *err, double tol, double eps);

/** Determines whether a set of points are inside a Boundary
* @param b Boundary (Boundary)
* @param points set of LLPoints to be tested (LLPoint)
* @param inside array of 1's (inside) and 0's (outside) (int)
* @param numberOfPts the number of points to be tested (int)
* @param tol Tolerance for function calculations (double)
* @param eps Convergence tolerance for Vincenty algorithm (double)
* @return Returns Error code that indicates success or cause of failure (ErrorSet*)
 */
ErrorSet ptsAreInsideBndry(Boundary b,
                           double latList[],
                           double lonList[],
                           int inside[],
                           int numberOfPts,
                           double tol,
                           double eps);

void MergeSort(double *array, int left, int right);
void sortPoints(double latList[], double lonList[], int idx[], int len);
void sortArrayByIndex(int array[], int idx[], int len);

/** Order a simple connected boundary by the sequence in which the shapes connect
 * at their respective start/end points, starting with the shape at index zero
 * of the original boundary's shape list
 *
 * A simple boundary implies one that contains no holes in the boundary area.  This also
 * means that the boundary area is not partitioned (contains sub areas) in any way
 *
 * @param boundary Boundary to be ordered (Boundary)
 * @param orderedBoundary Pointer to ordered Boundary (Boundary*)
 * @param tol Tolerance for function calculations (double)
 * @param eps Convergence tolerance for Vincenty algorithm (double)
 * @return Returns error code that indicates success or failure; updates given memory address 
 * for orderedBoundary
 * @retval SUCCESS Indicates successful execution.
 * @retval NO_MEMORY_ALLOCATED_ERR Indicates that a NULL pointer was passed for a required reference
 * @retval UNEXPECTED_ERR Indicates that an unknown error has occurred.
 */
ErrorSet orderBndry(Boundary boundary, Boundary *orderedBoundary, double tol, double eps);
/*
 *                    not partitioned          partitioned
 *                       ---------              ---------
 *                       |       |              |   |   |
 *                       |       |              |   |   |
 *                       |       |              |   |   |
 *                       ---------              ---------
 *
 * Let unordered and ordered be Boundary structs and
 * let the numbers suffix indicate the array index at which the given boundary shape is located
 *
 *                       unordered               ordered
 *
 *                          S_0                    S_0
 *                       ---------              ---------
 *                       |       |              |       |
 *                   S_3 |       | S_2      S_3 |       | S_1
 *                       |       |              |       |
 *                       ---------              ---------
 *                          S_1                    S_2
 *
 * In this example, the sequence {S_0, S_1, S_2, S_3} represents the order in which the shapes are defined
 * in the shapes array of the unordered boundary, i.e. S_0 = unordered.elements[0], 
 * S_1 = unordered.elements[1], etc
 * In this sequence, there exists pairs of consecutive shapes in the array such that these consecutive shapes
 * are not always connected.  For example the shapes located at indexes 0 and 1 (S_0 and S_1) do not share a
 * common start/end point */


/** Separate a boundary into simple ordered boundaries
 *
 * A simple boundary implies one that contains no holes in the boundary area.  This also
 * means that the boundary area is not partitioned (contains sub areas) in any way
 *
 * @param boundary Boundary to be ordered (Boundary)
 * @param separateBoundaries Pointer to ordered Boundary (Boundary*)
 * @param tol Tolerance for function calculations (double)
 * @param eps Convergence tolerance for Vincenty algorithm (double)
 * @return Returns error code that indicates success or failure; updates given memory address
 * for separateBoundaries
 * @retval SUCCESS Indicates successful execution.
 * @retval NO_MEMORY_ALLOCATED_ERR Indicates that a NULL pointer was passed for a required reference
 * @retval UNEXPECTED_ERR Indicates that an unknown error has occurred.
 */
ErrorSet separateBndry(Boundary boundary,
                       Boundary separateBoundaries[],
                       int *numberOfBoundaries,
                       double tol,
                       double eps);

/** Finds non-tangent intersections of a Geodesic with a Boundary
* NOTE: Intersections are only found on forward direction of geodesic,
* starting at point p following azimuth az.  This is sufficient for the
* ptIsInside algorithm
* @param b Boundary being intersected (Boundary)
* @param p = Start point of intersecting geodesic (LLPoint)
* @param az = Start azimuth of Geodesic (double)
* @param intx = Array of pointers to intersection points (LLPoint**)
* @param intxCount Number of intersections found (int*)
* @param exitCode = Code used to return intersection status and alert
*     ptIsInside function if vertex is hit or if given point
*     lies on boundary
* @param tol Tolerance for function calculations (double)
* @param eps Convergence tolerance for Vincenty algorithm (double)
* @return Returns error code that indicates success or failure (ErrorSet); updates given memory addresses for intx, intxCount, and exitCode variables with the calculated values
* @retval SUCCESS Indicates successful execution.
* @retval CONCENTRIC_CIRCLE_ERR Indicates that two arcs or circles either do not intersect or are identical. Used as a status code.
* @retval NO_MEMORY_ALLOCATED_ERR Indicates that a NULL pointer was passed for a required reference
* @retval LINE_TOO_FAR_FROM_ARC_ERR Status code indicates that no geodesic/locus-arc intersection could be found because the geodesic/locus is too far from the arc.
* @retval SHAPE_NOT_DEFINED_ERR Indicates that a shape object has not been defined.
*/
ErrorSet bndryGeoIntx(Boundary b,
                      LLPoint p,
                      double az,
                      LLPoint **intx,
                      int *intxCount,
                      int *exitCode,
                      double tol,
                      double eps);

/** Find intersections of a Locus with a Boundary
* @param b Boundary being intersected (Boundary)
* @param loc = Intersecting Locus (Locus)
* @param intx = Array of pointers to intersection points (LLPoint**)
* @param intxCount Number of intersections found (int*)
* @param tol Tolerance for function calculations (double)
* @param eps Convergence tolerance for Vincenty algorithm (double)
* @return Returns error code that indicates success or failure (ErrorSet); updates given memory addresses for intx and intxCount variables with the calculated values
* @retval SUCCESS Indicates successful execution.
* @retval NO_INTERSECTION_ERR Status code indicates that no intersection was found in the case that no intersection point gets returned.
* @retval NO_MEMORY_ALLOCATED_ERR Indicates that a NULL pointer was passed for a required reference
* @retval SHAPE_NOT_DEFINED_ERR Indicates that a shape object has not been defined.
 */
ErrorSet bndryLocusIntx(Boundary b, Locus loc, LLPoint **intx, int *intxCount, double tol, double eps);

/** Find intersections of an arc (or a circle) with a Boundary
* @param b Boundary being intersected (Boundary)
* @param a = Intersecting Arc (Arc)
* @param intx = Array of pointers to intersection points (LLPoint**)
* @param intxCount Number of intersections found (int*)
* @param tol Tolerance for function calculations (double)
* @param eps Convergence tolerance for Vincenty algorithm (double)
* @return Returns error code that indicates success or failure (ErrorSet); updates given memory addresses for intx and intxCount variables with the calculated values
* @retval SUCCESS Indicates successful execution.
* @retval CONCENTRIC_CIRCLE_ERR Indicates that two arcs or circles either do not intersect or are identical. Used as a status code.
* @retval NO_MEMORY_ALLOCATED_ERR Indicates that a NULL pointer was passed for a required reference
* @retval SHAPE_NOT_DEFINED_ERR Indicates that a shape object has not been defined.
 */
ErrorSet bndryArcIntx(Boundary b, Arc a, LLPoint **intx,
                      int *intxCount, double tol, double eps);

/** Find the Boundary that is the intersection of a circle with a given Boundary
* @param boundary Boundary being intersected (Boundary)
* @param circle Intersecting circle (Arc)
* @param newBoundary Pointer to the intersection boundary (Boundary*)
* @param tol Tolerance for function calculations (double)
* @param eps Convergence tolerance for Vincenty algorithm (double)
* @return Returns error code that indicates success or failure (ErrorSet); updates given memory address for the newBoundary variable with the calculated value
* @retval SUCCESS Indicates successful execution.
* @retval CONCENTRIC_CIRCLE_ERR Indicates that two arcs or circles either do not intersect or are identical. Used as a status code.
* @retval NO_MEMORY_ALLOCATED_ERR Indicates that a NULL pointer was passed for a required reference
* @retval CIRCLE_INSIDE_CIRCLE_ERR Status code indicates that no arc-arc intersection was found because one arc lies entirely inside the other.
* @retval SUBTENDED_ANGLE_OUT_OF_RANGE_ERR Indicates that arc's subtended angle does not meet algorithm requirement.
* @retval SHAPE_NOT_DEFINED_ERR Indicates that a shape object has not been defined.
* @retval UNEXPECTED_ERR Indicates that an unknown error has occurred.
 */
ErrorSet bndryCircleIntx(Boundary boundary, Arc circle, Boundary *newBoundary, double tol, double eps);

/** Determine whether a circle intersects a Boundary
* @param boundary Boundary being intersected (Boundary)
* @param circle Intersecting circle (Arc)
* @param checkSurface Treat Boundary as a surface (int)
* Case 1: checkSurface is non-zero
* In this case, the methods treats the boundary as a surface
* and if the circle lies on the boundary or intersects the boundary then true is returned
* Otherwise false is returned
*
* Case 2: checkSurface is zero
* In this case, the method only evaluates against the shapes that comprise the boundary
* itself.  If any intersection occurs with any shape that comprises the boundary, then
* true is returned.  Otherwise false is returned
* @param intersectionExists Indicates whether an intersection exists (1 is yes and 0 is no).  (int*)
* @param tol Tolerance for function calculations (double)
* @param eps Convergence tolerance for Vincenty algorithm (double)
* @return Returns error code that indicates success or failure (ErrorSet); updates given memory address for intersectionExists
* @retval SUCCESS Indicates successful execution.
* @retval CONCENTRIC_CIRCLE_ERR Indicates that two arcs or circles either do not intersect or are identical. Used as a status code.
* @retval NO_MEMORY_ALLOCATED_ERR Indicates that a NULL pointer was passed for a required reference
* @retval CIRCLE_INSIDE_CIRCLE_ERR Status code indicates that no arc-arc intersection was found because one arc lies entirely inside the other.
* @retval SUBTENDED_ANGLE_OUT_OF_RANGE_ERR Indicates that arc's subtended angle does not meet algorithm requirement.
* @retval SHAPE_NOT_DEFINED_ERR Indicates that a shape object has not been defined.
 */
ErrorSet bndryCircleIntxExists(Boundary boundary,
                               Arc circle,
                               int checkSurface,
                               int *intersectionExists,
                               double tol,
                               double eps);

/** Determine whether a given arc lies on another given arc
* @param arc1 The first arc (Arc)
* @param arc2 The second arc (Arc)
* @param commonShapes Shapes that are common to both arcs (Shape)
* @param shapeCount Number of common shapes (int*)
* @param err Error code that indicates success or cause of failure (ErrorSet*)
* @param tol Tolerance for function calculations (double)
* @param eps Convergence tolerance for Vincenty algorithm (double)
* @return Returns 1 if the arcs coincide and 0 if they do not; updates given memory addresses for commonShapes, shapeCount and error code with the calculated values
* @retval SUCCESS Indicates successful execution.
* @retval NO_MEMORY_ALLOCATED_ERR Indicates that a NULL pointer was passed for a required reference
* @retval SHAPE_NOT_DEFINED_ERR Indicates that a shape object has not been defined.
 */
int arcsCoincide(Arc arc1, Arc arc2, Shape commonShapes[], int *shapeCount, ErrorSet *err,
                 double tol, double eps);

/** Determine if the given point is already in the pointset
 * @param pt Input point (LLPoint)
 * @param pts Input pointset (LLPointSet)
 * @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
 * @return Returns 1 if the input point is already in the point set, 0 if not
 */
int ptIsInSet(LLPoint pt, LLPointSet pts, double tol);

/** Determine the principal value (-PI < pangle <= PI) of an angle
 * @param angle Angle to be recalculated (double)
 * @return Returns the recalculated angle
 */
double principalAngle(double angle);

/** Calculates the angle between the input values in the correct orientation
 * @param startAngle Initial input angle (double)
 * @param endAngle Final input angle (double)
 * @param dir Direction to calculate subtended angle [1, -1] (int)
 * @return Returns the calculated subtended angle.  This value will never be negative
 */
double calculateSubtendedAngle(double startAngle, double endAngle, int dir);

/** Calculates the difference between two courses.  The return value will always be < PI and > -PI
 * @param az1 First input azimuth (double)
 * @param az2 Second input azimuth (double)
 * @return  Returns the difference in azimuths.  -PI < returned value < PI
 */
double azDifference(double az1, double az2);

/** Determines if an azimuth is within the extent of a spiral
 * @param sp Input spiral (Spiral)
 * @param az Azimuth to be checked (double)
 * @return Returns 1 if the azimuth is within the spiral extent, 0 if not
 */
int azIsInSpiralExtent(Spiral sp, double az);

/** Calculates the course of the spiral at a given azimuth from the centerpoint
 * @param sp Input spiral (Spiral)
 * @param az Azimuth at which to calculate the tangent course (double)
 * @param tanCrs Pointer to a double to be updated with calculated course (double*)
 * @return Updates input pointer with calculated course
 * @retval SUCCESS Indicates successful execution.
 */
ErrorSet spiralTanCrs(Spiral sp, double az, double *tanCrs);

/** Creates the spiral section with the same properties of the input spiral centered on the input azimuth and input radius
 * @param inputsp Input spiral to be recreated (Spiral)
 * @param az Azimuth to center the spiral section on (double)
 * @param rad Reference radius for constructing the spiral section (double)
 * @param sp Pointer to a spiral structure to be updated with new section (Spiral*)
 * @param eps Convergence tolerance for Vincenty algorithm (double)
 * @return Updates the input pointer with the constructed spiral section
 * @retval SUCCESS Indicates successful execution.
 */
ErrorSet createSpiralSection(Spiral inputsp, double az, double rad, Spiral *sp, double eps);

/** Determines a point on the input spiral:  If the geodesic doesn't intersect the spiral, the calculated point
 * is the closest point on the spiral to the geodesic.  If the geodesic intercepts the spiral, the calculated point
 * is the point on the spiral farthest from the geodesic and on the opposite side of the line as the spiral center point.
 *
 * @param sp Input spiral (Spiral)
 * @param line Input geodesic (Geodesic)
 * @param midChord Pointer to an LLPoint struct to be updated with calculated point (LLPoint*)
 * @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
 * @param eps Convergence tolerance for Vincenty algorithm (double)
 * @return Updates the input pointer with the calculated point
 * @retval SUCCESS Indicates successful execution.
 */
ErrorSet spiralMidChord(Spiral sp, Geodesic line, LLPoint *midChord, double tol, double eps);

/** Determines if the two input boundaries intersect.  Return value indicates the type of intersection:
 * intersectionExists = 0:  No intersection occurs between the boundaries
 * intersectionExists = 1:  The boundaries intersect
 * intersectionExists = 2:  Boundary 2 is fully inside Boundary 1
 * intersectionExists = 3:  Boundary 1 is fully inside Boundary 2
 *
 * @param B1 First input boundary (Boundary)
 * @param B2 Second input boundary (Boundary)
 * @param intersectionExists Return value indicating the type of intersection (int*)
 * @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
 * @param eps Convergence tolerance for Vincenty algorithm (double)
 * @return Updates the input pointer with the calculated point
 * @retval SUCCESS Indicates successful execution.
 */
ErrorSet bndryIntxExists(Boundary B1, Boundary B2, int *intersectionExists, double tol, double eps);

/** Projects a given point onto a boundary and returns the closest projection point if one exists.
 *
 * @param b Input boundary (Boundary)
 * @param pt Input point (LLPoint)
 * @param perpPt The projected point onto the boundary (LLPpoint*)
 * @param perpDist The distance between the input point and the projected point (double*)
 * @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
 * @param eps Convergence tolerance for Vincenty algorithm (double)
 * @return Updates the input pointer with the calculated point
 * @retval SUCCESS Indicates successful execution.
 */

ErrorSet projectToBndry(Boundary b, LLPoint pt, LLPoint *perpPt, double *perpDist, double tol, double eps);

} // namespace


