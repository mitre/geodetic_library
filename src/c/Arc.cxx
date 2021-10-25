/*
 * Copyright 2007-2021 The MITRE Corporation.  All Rights reserved.
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

/* \file Arc.c
 *  \brief This file defines the arc functions used by geolib.
 *  \author Michael Mills, Richard Snow, Stuart Bowman, Juan Amezcua, John Landrigan
 */

#include <sys/types.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#if REPLACE_WITH_AMDLIBM
#include "amdlibm.h"
#endif
#include "Geolib.h"

namespace geolib_idealab {

/*******************************************************************************
 * Return arc's subtended angle in radians
 *  The sign of the returned value is the same as the sign of the arc's
 * orientation.  In other words, if the orientation < 0 (clockwise),
 * then arc extent < 0.  This is counterintuitive.  If you want to calculate
 * the arc extent with the opposite sign convention, use computeSubtendedAngle.
 *
 * NOTE: This function treats startCrs and endCrs as if they were of infinite
 * precision.  If you want 2*PI to be returned in the start/end points of the
 * arcs are closer than your distance tolerance, then you have to check that
 * before calling this function.
 *  */

ErrorSet getArcExtent(LLPoint center, double radius, double startCrs, double endCrs,
        ArcDirection orientation, double *arcExtent, double tol, double vincentyEps)
{

	ErrorSet err = 0;
    double distToPoint, tempStartCrs, tempEndCrs;
    LLPoint startPoint, endPoint;
    int ptsAreClose = 0;

    //map courses to [0,2Pi)
    tempStartCrs = modpos(startCrs, M_2PI);
    tempEndCrs = modpos(endCrs, M_2PI);

    //Calculate the start/end points of the given arc
    err |= direct(center, tempStartCrs, radius, &startPoint, vincentyEps);

    err |= direct(center, tempEndCrs, radius, &endPoint, vincentyEps);

    //Check if the arc start and end points are within tolerance of each other
    err |= invDist(startPoint, endPoint, &distToPoint, vincentyEps);

    if (err) {
    	*arcExtent = 0;
    	return err;
    }
    else if (distToPoint <= tol) ptsAreClose = 1;

    //Compute the arc extent
    *arcExtent = computeSubtendedAngle(tempStartCrs, tempEndCrs, orientation);

    if ((tempStartCrs != tempEndCrs) && (ptsAreClose)){

    	/*
    	 * Since points are close, arc extent should either be really small or really close to 2Pi
    	 * No need to check to angular tolerance critera since the arc extent should be of several
    	 * orders of magnitude of difference with respect to Pi.
    	*/
    	if (fabs(*arcExtent) < M_PI) *arcExtent = 0.0;
    	else if (fabs(*arcExtent) > M_PI) *arcExtent = M_2PI * ((double)orientation);
    	else {
    		//if fabs(arc extent) ~= PI and points are close something is really wrong
    		err |= UNEXPECTED_ERR;
    		*arcExtent = 0;
    	}

    }

    return err;
}

int azIsInArcExtent(Arc arc, double az)

{
	double angle;
	int val = 0;

	angle = calculateSubtendedAngle(arc.startAz, az, arc.dir);
	if (fabs(arc.subtendedAngle) >= fabs(angle))
	{
		val = 1;
	}
	return val;
}

/*******************************************************************************
 * Test whether given point lies on arc between start and end courses.
 *
 * Input:
 *
 * Input:
 *   center      = The arc center point
 *   radius      = The arc radius in nautical miles
 *   startCrs    = The course from the arc center point to the arc start point, in radians
 *   endCrs      = The course from the arc center point to the arc end point, in radians
 *   orientation = +1 if arc is traversed counter-clockwise (right handed)
 *                 -1 if arc is traversed clockwise
 *   testPt      = The point to be tested
 *   err         = A pointer to store any errors
 *   tol         = The distance tolerance in nautical miles
 *   eps         = The convergence criteria for the Vincenty Inverse/Forward algorithms, in radians
 *
 * Output:
 * 		Return 1 if point is on the arc.  Return 0 otherwise.
 *
 */

int ptIsOnArc(LLPoint center, double radius, double startCrs, double endCrs,
        ArcDirection orientation, LLPoint testPt, ErrorSet* err, double tol,
        double vincentyEps)
{

	double distToPoint, crsToPoint, crsFromPoint;
	LLPoint startPoint, endPoint;
	double arcExtent, subExtent;
	ErrorSet newErr = 0;

	//Calculate the start/end points of the given arc
	newErr |= direct(center, startCrs, radius, &startPoint, vincentyEps);

	newErr |= direct(center, endCrs, radius, &endPoint, vincentyEps);

	//Check if the test point is within the neighborhood of the arc start point
	newErr |= invDist(startPoint, testPt, &distToPoint, vincentyEps);

	if (newErr) {
		*err |= newErr;
		return 0;
	}
	else if (distToPoint <= tol) return 1;

	//Check if the test point is within the neighborhood of the arc end point
	newErr |= invDist(endPoint, testPt, &distToPoint, vincentyEps);

	if (newErr) {
		*err |= newErr;
		return 0;
	}
	else if (distToPoint <= tol) return 1;

	//Get the forward/inverse course with respect to test point and arc center point
	//Get the distance in nautical miles from test point to arc center point
	newErr |= inverse(center,testPt,&crsToPoint,&crsFromPoint,&distToPoint,vincentyEps);

	//Return 0 and error value if Vincenty Inverse fails for any reason
	if (newErr) {
		*err |= newErr;
		return 0;
	}

	//Check if the test point is outside the neighborhood around the arc for a distance of tol
	//The arc extent may be assumed to be a full circle at this point
	if (fabs(radius-distToPoint) > tol)
	{
		//test point is outside neighborhood around arc and thus cannot be on the arc
		return 0;
	}

	//Get the actual arc extent of the given arc
	newErr |= getArcExtent(center, radius, startCrs, endCrs, orientation, &arcExtent, tol, vincentyEps);

	if (newErr){
		*err |= newErr;
		return 0;
	}

	//Check whether the arc is actually a full circle
	/*If so then the test point must be on the arc since we know it lies within its neighborhood
	* Note: The distance check of the start/end points occurs when getArcExtent is called.
	*/
	if (fabs(arcExtent) >= M_2PI) return 1;

	/* find extent of new arc with same startCrs but testPt as endCrs */
	*err |= getArcExtent(center, radius, startCrs, crsToPoint,orientation, &subExtent, tol, vincentyEps);

	if (newErr){
		*err |= newErr;
		return 0;
	}

	//Check whether the test point is in the arc extent of the given arc
	/*Note: If you have made it this far in the code then the test point is not within a distance of tol
	*      of the arc start/end point.  Hence the subextent must be greater than angleTol where the
	*      value of angleTol should be a function of the value of tol (due to the relationship d = r(theta)).
	*      Thus checking that the traversed arc extent difference is less than zero is sufficient as opposed
	*      to checking that it is less than angleTol
	*/
	if (fabs(subExtent) - fabs(arcExtent) < 0)
	{
	/* traversing arc from startPt, one would run into testPt before endPt */
		return 1;
	}

	return 0;

}

/*******************************************************************************
 * Test whether given point lies inside or on the boundary of the "pie slice"
 * created by the given arc start, end, and center points.
 *
 * Input:
 *   center      = The arc center point
 *   radius      = The arc radius in nautical miles
 *   startCrs    = The course from the arc center point to the arc start point, in radians
 *   endCrs      = The course from the arc center point to the arc end point, in radians
 *   orientation = +1 if arc is traversed counter-clockwise (right handed)
 *                 -1 if arc is traversed clockwise
 *   testPt      = The point to be tested
 *   err         = A pointer to store any errors
 *   tol         = The distance tolerance in nautical miles
 *   eps         = The convergence criteria for the Vincenty Inverse/Forward algorithms, in radians
 *
 * Output:
 * 		Return 1 if point is inside or on the boundary.  Return 0 otherwise.
 *
 */

int ptIsInsideArc(LLPoint center, double radius, double startCrs, double endCrs,
        ArcDirection orientation, LLPoint testPt, ErrorSet* err, double tol,
        double vincentyEps)
{

	double distToPoint, crsToPoint, crsFromPoint;
	double radiusMinusDistToPoint;
	LLPoint projectedTestPt;
	ErrorSet newErr = 0;
	int isOnArc;

	if(ptsAreSame(center,testPt,tol)){
		return 1;
	}

	//Find the start true course from the arc center point to the test point
	//Find the distance between the arc center point and the test point
	newErr |= inverse(center, testPt, &crsToPoint, &crsFromPoint, &distToPoint, vincentyEps);
	if (newErr) {
		*err |= newErr;
		return 0;
	}

	//Check if point is further from the center point than the radius
	if ( distToPoint > radius+tol)
	{
		return 0;
	}

	//Check if test point is within tol of the arc
	if (fabs(distToPoint - radius) <= tol) {

		//If so, no need to project the test point.
		projectedTestPt = testPt;

	} else {
		//Calculate the difference between the radius and the distance to the test point from the center
		radiusMinusDistToPoint = radius - distToPoint;

		//Project the test point a distance of radiusMinusDistToPoint along the crsToPoint
		//At this point, treat the arc as if it were a complete circle.  It may or may not be, but
		//assume it to be so for now
		//Hence the projected point must be on the same circle that the given arc is a part of
		newErr |= direct(testPt, modcrs(crsFromPoint + M_PI), radiusMinusDistToPoint, &projectedTestPt, vincentyEps);
		if (newErr) {
			*err |= newErr;
			return 0;
		}
	}

	//Check if the projected point is on the actual arc (no longer assuming it is a circle)
	//If the projected point is on the actual arc, then we must also have that the test point
	//is inside the "pie slice" defined by the given arc and it's corresponding center point
	//Note: this also handles the pie slice edge areas that result due to the bubble (neighborhood)
	//of distance tol around the end points of the arc (should it not be a complete circle).
	isOnArc = ptIsOnArc(center, radius, startCrs, endCrs, orientation, projectedTestPt, &newErr, tol, vincentyEps);

	if (newErr) {
		*err |= newErr;
		return 0;
	}

	return isOnArc;

}

/******************************************************************************
 * Discretizes the arc into small pie-slices.  Finds three points on
 * circ.  of each pie slice and calculates the radius of curvature and
 * central angle for that slice.  ROC*angle = length of slice.  Adds
 * up the lengths of the slices to get total arc length.
 */

double arcLength(LLPoint center, double radius,
                                 double startCrs, double endCrs, ArcDirection orient,
                                 int *n, ErrorSet* err, double tol, double eps)
{

    LLPoint p1, p2, p3;
    Vector v0, v1, v2, v3;
    Vector chord1, chord2;
    int i;
    int k = 0;

    double length, error;
    double meanR;
    double oldLength = approxArcLength(center, radius, startCrs, endCrs,
            orient, &meanR);
    double subtAngle;
    double theta, dtheta;

    double x1, x2, phi1, phi2, xi, psi, phibar;
    double arg, sigma, R, cosphi2;
    int npN = 5;

    ErrorSet newErr = 0;

    /* Assign local storage if optional pointers are not provided */
    if (n == NULL){
    	n = &npN;
    } else if (*n < 10)
        *n = 5;

    //Map arc center point to ECEF coordinate system
    v0 = geodeticToECEF(center);

    //Determine the subtended angle
    subtAngle = computeSubtendedAngle(startCrs, endCrs, orient);

    //Loop until the change in arc length between loops is within tolerance
    while ((k == 0) || ((fabs(error) > tol) && (k <= MAX_ITERATIONS) 
            ))
    {
    	//Increase the number of subsegments
        *n += 5;

        //Calculate the angle of the subarc
        dtheta = subtAngle / ((double) *n);

        length = 0.0;

        //Calculate the start point of the arc and map it to ECEF coordinate system
        newErr |= direct(center, startCrs, radius, &p1, eps);
        if (newErr)
        {
            *err |= newErr;
            return 0.0;
        }
        v1 = geodeticToECEF(p1);

        //Loop through n arcs and calculate the arc lengths
        for (i = 0; i < *n; i++)
        {

        	//Calculate points 2 and 3 and map them to ECEF coordinate system
            theta = startCrs + ((double) i) * dtheta;

            newErr = direct(center, theta + 0.5 * dtheta, radius, &p2, eps);
            if (newErr)
            {
                *err |= newErr;
                return 0.0;
            }
            v2 = geodeticToECEF(p2);
            newErr |= direct(center, theta + dtheta, radius, &p3, eps);
            if (newErr)
            {
                *err |= newErr;
                return 0.0;
            }
            v3 = geodeticToECEF(p3);

            //Calculate the length of the current sub-arc
            chord1 = vectorSubtract(v1, v2);
            chord2 = vectorSubtract(v3, v2);

            x1 = norm(chord1);
            x2 = norm(chord2);

            xi = dot(chord1, chord2) / x1 / x2; /* cos(phibar) */
            sigma = sqrt(1.0 - xi * xi); /* sin(phibar) */
            phibar = acos(xi);
            arg = x1 / x2 - xi;
            cosphi2 = sigma / sqrt(arg * arg + sigma * sigma);
            phi2 = acos(cosphi2);
            R = x2 / 2.0 / cosphi2;
            phi1 = phibar - phi2;

            psi = M_2PI - 2.0 * (phi1 + phi2);

            length += R * psi;

            p1 = p3;
            v1 = v3;

        }

        //Calculate difference between previous arc length and most recent arc length
        error = length - oldLength;

        oldLength = length;

        k = k + 1;

    }

    if (k >= MAX_ITERATIONS)
    {
        *err |= ITERATION_MAX_REACHED_ERR;

    }

    if (fabs(error) >= MAX_DISTANCE_ERROR)
    {
        *err |= ERROR_MAX_REACHED_ERR;
    }

    return length;

}

/******************************************************************************
 * This is a non-iterative approximation to arc length.  This function
 * calculates the best-fit radius of curvature for the ellipsoid in
 * the vicinity of the arc.  Then it treats the ellipsoidal arc as a
 * spherical arc with constant ROC.  Accuracy of this method degrades
 * as the size of the arc increases.
 */

double approxArcLength(LLPoint center, double radius, double startCrs,
                            double endCrs, ArcDirection orient, double* meanR)
{

    double f = FLATTENING;
    double ee = f * (2.0 - f); /* eccentricity squared */
    double a = SEMI_MAJOR_AXIS_NMI;
    double lat = center.latitude; /* geodetic latitude at center */
    double sinLat = sin(lat);
    double sinSq = sinLat * sinLat;
    double M = a * (1.0 - ee) / pow(1.0 - ee * sinSq, 1.5);
    double N = a / sqrt(1.0 - ee * sinSq);
    double subtAngle = computeSubtendedAngle(startCrs, endCrs, orient);
    double arcLength;
	double npMeanR;

    /* Assign local storage if optional pointers are not provided */
    if (NULL == meanR) meanR = &npMeanR;

    //Calculate the approximate arc length
    *meanR = sqrt(M * N);
    arcLength = fabs(subtAngle) * (*meanR) * sin(radius / (*meanR));

    return arcLength;

}

/* Find a point on an arc that is a given distance from another point on the
* arc.  If the given distance is positive, it is measured along the arc in the same
* direction as the arc orientation.  If the distance is negative, then the new point is
* computed by moving opposite the arc orientation. */
ErrorSet arcFromLength(LLPoint center, LLPoint givenPoint,
                             ArcDirection dir, double arcDist,
                             LLPoint* endPoint, double *subtendedAngle,
                             double tol, double eps)
{

    double radius;
    double actualDist;
    double startAz, tmpAz, endAz;
    double chordR;
    double error;

    ErrorSet err = 0;

    int k = 0;
    int n = 0;

	double npSubtendedAngle;
	LLPoint npEndPoint;

    /* Assign local storage if optional pointers are not provided */
    if (NULL == subtendedAngle) subtendedAngle = &npSubtendedAngle;
    if (NULL == endPoint) endPoint = &npEndPoint;

    //Zero distance given, return input point
    if (fabs(arcDist) < tol)
    {
        *subtendedAngle = 0;
        *endPoint = givenPoint;
        return err;
    }
    else if (arcDist < 0)
    //Input distance is negative, swap arc orientation
    {
        arcDist = -arcDist;
        if (dir == CLOCKWISE)
            dir = COUNTERCLOCKWISE;
        else
            dir = CLOCKWISE;
    }

    err |= inverse(center, givenPoint, &startAz, &tmpAz, &radius, eps);
    if (err)
        return err;

    /* Calculate radius of spherical small circle */
    if (radius > SPHERE_RADIUS)
    {
        err |= RADIUS_OUT_OF_RANGE_ERR;
        return err;
    }

    /* Spherical approximation is initial guess */
    chordR = SPHERE_RADIUS * sin(radius / SPHERE_RADIUS);
    *subtendedAngle = arcDist / chordR;

    //Iterate until arc length equals input distance
    while ((k == 0) || ((fabs(error) > tol) && (k < MAX_ITERATIONS) ))
    {
    	//Calculate the end azimuth, compute the arc length and check the result
        if (dir == CLOCKWISE)
            endAz = startAz + *subtendedAngle;
        else
            endAz = startAz - *subtendedAngle;

        n = 0;

        actualDist = arcLength(center, radius, startAz, endAz,
                dir, &n, &err, tol, eps);

        if (err)
            return err;

        error = actualDist - arcDist;

        //Calculate next guess
        *subtendedAngle = *subtendedAngle * arcDist / actualDist;

        k = k + 1;

    }

    if (k >= MAX_ITERATIONS)
    {
        err |= ITERATION_MAX_REACHED_ERR;
    }

    if (fabs(error) >= MAX_DISTANCE_ERROR)
    {
        err |= ERROR_MAX_REACHED_ERR;
    }

    err |= direct(center, endAz, radius, endPoint, eps);

    *subtendedAngle  = fmod(*subtendedAngle,M_2PI);

    return err;

}

/*******************************************************************************
 * Initialize a new Arc structure from given input parameters
 */
ErrorSet createArc(Arc* arc, LLPoint centerPoint, LLPoint startPoint,
                 LLPoint endPoint, ArcDirection dir, double tol, double eps)
{

    double radius = 0.0;
    double startRad = 0.0;
    double endRad = 0.0;
    double startAz = 0.0;
    double endAz = 0.0;
    double tmpAz = 0.0;

    ErrorSet err = 0;

    Arc npArc;
    /* Assign local storage if optional pointers are not provided */
    if (NULL == arc) arc = &npArc;


    err = inverse(centerPoint, startPoint, &startAz, &tmpAz, &startRad,
            eps);
    err |= inverse(centerPoint, endPoint, &endAz, &tmpAz, &endRad, eps);

    radius = 0.5 * (startRad + endRad);
    if (radius <= 0.0)
    {
        //Make sure radius is non-zero & positive
        err |= RADIUS_OUT_OF_RANGE_ERR;
        return err;
    }
    if ((fabs(startRad - radius) > tol) || (fabs(endRad - radius) > tol))
    {
        //Throw if start & end points are not on Arc
        err |= POINT_NOT_ON_ARC_ERR;
    }
    else
    {
        arc->centerPoint = centerPoint;
        arc->startPoint = startPoint;
        arc->endPoint = endPoint;
        arc->startAz = startAz;
        arc->endAz = endAz;
        arc->radius = radius;
        arc->dir = dir;
        arc->subtendedAngle = computeSubtendedAngle(startAz, endAz, dir);
    }

    return err;

}

/*******************************************************************************
 * Initialize a new Arc structure from given input parameters
 */
ErrorSet createArcFromCourses(Arc* arc, LLPoint centerPoint, double radius,
                 double startAz, double endAz, ArcDirection dir, double tol, double eps){

    ErrorSet err = 0;
    LLPoint startPoint;
    LLPoint endPoint;

    if (radius <= tol){
        //Make sure radius is non-zero & positive
        err |= RADIUS_OUT_OF_RANGE_ERR;
        return err;
    }

    err |= direct(centerPoint, startAz, radius, &startPoint, eps);
    err |= direct(centerPoint, endAz, radius, &endPoint, eps);

	arc->centerPoint = centerPoint;
	arc->startPoint = startPoint;
	arc->endPoint = endPoint;
	arc->startAz = startAz;
	arc->endAz = endAz;
	arc->radius = radius;
	arc->dir = dir;
	arc->subtendedAngle = computeSubtendedAngle(startAz, endAz, dir);

    return err;
}

/*******************************************************************************
 * Initialize a new Arc structure from given input parameters
 */
ErrorSet createArcCircle(Arc* arc, LLPoint centerPoint, double radius, ArcDirection dir, double tol, double eps){

    ErrorSet err = 0;
    LLPoint startPoint;
    LLPoint endPoint;

    if (radius <= tol){
        //Make sure radius is non-zero & positive
        err |= RADIUS_OUT_OF_RANGE_ERR;
        return err;
    }

    err |= direct(centerPoint, 0.0, radius, &startPoint, eps);
    err |= direct(centerPoint, 0.0, radius, &endPoint, eps);

	arc->centerPoint = centerPoint;
	arc->startPoint = startPoint;
	arc->endPoint = endPoint;
	arc->startAz = 0.0;
	arc->endAz = 0.0;
	arc->radius = radius;
	arc->dir = dir;
	arc->subtendedAngle = computeSubtendedAngle(0.0, 0.0, dir);

    return err;
}

/*************************************************************************
 * Construct arc from start point, start course, and radius.  The arc's
 * end point is computed such that the course from the end point to the
 * given nextEndPoint is tangent to the arc at the end point.
 */
ErrorSet arcEndFromStartAndRadius(LLPoint arcStart, double startCrs, double radius, ArcDirection arcDir,  LLPoint nextEndPoint, Arc* newArc, double tol, double eps)
{
  ErrorSet err = SUCCESS;
  LLPoint arcCenter, arcEnd, tempLLPoint;
  double perpCrs, dist, az1, az2, azTowardCenter;
  LLPointPair tanPt;
  int i, n;
  Arc npArc;

  /* Assign local storage if optional pointers are not provided */
  if (newArc == NULL) newArc = &npArc;

  if (arcDir == CLOCKWISE)
    perpCrs = modcrs(startCrs + M_PI_2);
  else
    perpCrs = modcrs(startCrs - M_PI_2);

  /* Find arc center and do error checking */
  err |= direct(arcStart, perpCrs, radius, &arcCenter, eps);
  err |= invDist(arcCenter, nextEndPoint, &dist, eps);

  if ((radius < tol) || (dist < radius + tol))
  {
    newArc->radius = 0.0;
	newArc->startPoint = arcStart;
	newArc->endPoint = arcStart;
	newArc->centerPoint = arcStart;
	newArc->dir = arcDir;
	newArc->startAz = 0.0;
	newArc->endAz = 0.0;
	return err;
  }

  /* Find tangent points */
  err |= ptsOnArcOnTanThruPt(nextEndPoint, arcCenter, radius, tanPt, &n, tol, eps);

  /* Select correct tangent point */
  for (i = 0; i < n; i++)
  {
    err |= invCrs(tanPt[i], nextEndPoint, &az1, &az2, eps);
	if (arcDir == CLOCKWISE)
	  azTowardCenter = az1 + M_PI_2;
    else
	  azTowardCenter = az1 - M_PI_2;

    err |= direct(tanPt[i], azTowardCenter, radius, &tempLLPoint, eps);
	err |= invDist(arcCenter, tempLLPoint, &dist, eps);
	if (dist < radius)
	{
	  arcEnd.latitude = tanPt[i].latitude;
	  arcEnd.longitude = tanPt[i].longitude;
	  break;
	}
  }

  /* Create new arc */
  err |= createArc(newArc, arcCenter, arcStart, arcEnd, arcDir, tol, eps);

  return err;
}

/*************************************************************************
 * Construct arc from start point, start course, and center.  The arc's
 * end point is computed such that the course from the end point to the
 * given nextEndPoint is tangent to the arc at the end point.
 */
ErrorSet arcEndFromStartAndCenter(LLPoint arcStart, double startCrs, LLPoint arcCenter, ArcDirection arcDir,  LLPoint nextEndPoint, Arc* newArc, double tol, double eps)
{
  ErrorSet err = SUCCESS;
  LLPoint arcEnd, tempLLPoint;
  double radius, dist, az1, az2, azTowardCenter;
  LLPointPair tanPt;
  int i, n;
  Arc npArc;

  /* Assign local storage if optional pointers are not provided */
  if (newArc == NULL) newArc = &npArc;

  /* Find arc radius and do error checking */
  err |= invDist(arcStart, arcCenter, &radius, eps);
  err |= invDist(arcCenter, nextEndPoint, &dist, eps);

  if ((radius < tol) || (dist < radius + tol))
  {
    newArc->radius = 0.0;
	newArc->startPoint = arcStart;
	newArc->endPoint = arcStart;
	newArc->centerPoint = arcStart;
	newArc->dir = arcDir;
	newArc->startAz = 0.0;
	newArc->endAz = 0.0;
	return err;
  }

  /* Find tangent points */
  err |= ptsOnArcOnTanThruPt(nextEndPoint, arcCenter, radius, tanPt, &n, tol, eps);

  /* Select correct tangent point */
  for (i = 0; i < n; i++)
  {
    err |= invCrs(tanPt[i], nextEndPoint, &az1, &az2, eps);
	if (arcDir == CLOCKWISE)
	  azTowardCenter = az1 + M_PI_2;
    else
	  azTowardCenter = az1 - M_PI_2;

    err |= direct(tanPt[i], azTowardCenter, radius, &tempLLPoint, eps);
	err |= invDist(arcCenter, tempLLPoint, &dist, eps);
	if (dist < radius)
	{
	  arcEnd.latitude = tanPt[i].latitude;
	  arcEnd.longitude = tanPt[i].longitude;
	  break;
	}
  }

  /* Create new arc */
  err |= createArc(newArc, arcCenter, arcStart, arcEnd, arcDir, tol, eps);

  return err;
}
} //namespace
