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
 * Arc Tangent to two Lines
 *
 */
ErrorSet arcTanToTwoGeos(LLPoint pt1, double crs12, LLPoint pt3,
                                double crs3, double radius,
                                LLPoint* centerPoint, LLPoint* startPoint,
                                LLPoint* endPoint, ArcDirection* dir,
                                double tol, double eps)
{
    /*
     * INPUT:
     *   pt1: start point of first line
     *   crs12: azimuth of line from pt1 to pt2 (measured at pt1)
     *   pt3: end point of second line
     *   crs3: azimuth of line from pt2 to pt3 (measured at pt3)
     *   radius: radius of desired tangent arc
     *   tol: maximum error in calculated point position
     *   eps: round-off error (also accuracy of forward and inverse calcs
     * OUTPUT:
     *   arcPts: array of three points: arc center, tangent point 1, and
     *           tangent point 2
     *   dir: If zero then no arc was found.  NEEDS TO BE CHECKED BY CALLING
     *        FUNCTION along with error code.
     *
     */

    ErrorSet err = 0;

    LLPoint pt2; /* vertex point */
    LLPoint iterationStart;

    double distArray[2] = { 0.0, 0.0 }, errArray[2] = { 0.0, 0.0 };
    double cca; /* course change angle, must be in range [-pi,pi] */
    double vertexAngle; /* interior angle at vertex pt2.  In range [0,pi) */
    double error = 0.0;
    ErrorSet tempErr = 0;
    double crs21, crs23, dist12, dist23;
    double crs1;
    double tmpCrs, perpCrs, perpDist, distToStart;
    double adjust = 0.0;
    double stepSize = 1.0e99;
    double B; /* value for spherical approx */
    //double  crsTest1,distTest1;
    int rht = 0; /* true for right hand turn */
    int k = 0;
    int revDirFlag = 0;
    double c12, c21, c32, c23;

#ifdef USE_BEST_FIT_ROC

    double sphereRad =
    lookUpROC(geocentricLat(0.5*fabs(pt1.latitude+pt3.latitude)));
#else

    double sphereRad = SPHERE_RADIUS;
#endif

    LLPoint npCenterPoint;
    LLPoint npStartPoint;
    LLPoint npEndPoint;
    ArcDirection npDir;

    /* Assign local storage if optional pointers are not provided */
    if (NULL == centerPoint) centerPoint = &npCenterPoint;
    if (NULL == startPoint) startPoint = &npStartPoint;
    if (NULL == endPoint) endPoint = &npEndPoint;
    if (NULL == dir) dir = &npDir;

    //*dir = 0; /* Initial direction parameter, indicates no arc yet found */

    /* There are NINE cases, depending on geometry */
    /* Case 1: pt1 and pt3 are the same */
    if (ptsAreSame(pt1, pt3, tol))
    {
        pt2 = pt1;
        crs1 = reciprocal(crs12);
        /* crs3 does not change */
    }
    /* Case 2: pt3 is on forward geodesic through pt1
     * NOTE: This should be the most common situation for aviation */
    else if (ptIsOnCrs(pt1, crs12, pt3, &crs21, &dist12, &err, tol, eps))
    {
        if (ptIsOnCrs(pt3, crs3, pt1, NULL, NULL, &tempErr, tol, eps)
                || ptIsOnCrs(pt3, crs3 + M_PI, pt1, NULL, NULL, &tempErr,
                        tol, eps))
        {
            /* Geodesics are collinear */
            err |= NO_TANGENT_ARC_ERR;
            return err;
        }
        pt2 = pt3;
        crs1 = crs21;
        /* crs3 does not change */
    }
    /* Case 3: pt3 is on reciprocal geodesic through pt1 */
    else if (ptIsOnCrs(pt1, crs12 + M_PI, pt3, &crs21, &dist12, &err, tol,
            eps))
    {
        if (ptIsOnCrs(pt3, crs3, pt1, NULL, NULL, &tempErr, tol, eps)
                || ptIsOnCrs(pt3, crs3 + M_PI, pt1, NULL, NULL, &tempErr,
                        tol, eps))
        {
            /* Geodesics are collinear */
            err |= NO_TANGENT_ARC_ERR;
            return err;
        }
        pt2 = pt3;
        crs1 = reciprocal(crs21);
        /* crs3 does not change */
    }
    /* Case 4: pt1 is on forward geodesic through pt3 */
    else if (ptIsOnCrs(pt3, crs3, pt1, &crs23, &dist23, &err, tol, eps))
    {
        pt2 = pt1;
        crs1 = reciprocal(crs12);
        crs3 = reciprocal(crs23);
    }
    /* Case 5: pt1 is on reciprocal geodesic through pt3 */
    else if (ptIsOnCrs(pt3, crs3 + M_PI, pt1, &crs23, &dist23, &err, tol,
            eps))
    {
        pt2 = pt1;
        crs1 = reciprocal(crs12);
        crs3 = crs23;
    }
    /* Remaining cases require finding the intersection of the two geodesics */
    else
    {
        /* Find nearest intersection of geodesics */
        err |= crsIntx(pt1, crs12, &crs21, &dist12, pt3, crs3,
                &crs23, &dist23, &pt2, tol, eps);

        /* Find the courses used to determine location of p2 w/ respect to p1 and p3 */
        err |= invCrs(pt1, pt2, &c12, &c21, eps);
        err |= invCrs(pt3, pt2, &c32, &c23, eps);

        if (err)
        {
            /* CANNOT CONTINUE if no intersection is found or unable to determine courses to
             * setup geometry of the remaining four cases*/
            return err;
        }

        /* Roughly compare azimuths of geodesics at given points with those at intersection
         * to determine which side of intersection the given points are on.*/
        /* Combination of the following two options yield the remaining 4 cases */
        if (fabs(modlon(c12 - crs12)) < M_PI_2)
            /* pt2 is on foward geodesic through pt1
             * (i.e. pt2 lies in front of pt1)
             */
            crs1 = crs21;
        else
            crs1 = reciprocal(crs21);

        if (fabs(modlon(c32 - crs3)) < M_PI_2)
            /* pt2 is on forward geodesic through pt3
             * (i.e. pt2 lies in front of pt3)
             * */
        	crs3 = reciprocal(crs23);
        else
        	crs3 = crs23;

    }

    /* Geometry is now set up:
     *   * pt2 is at intersection of the two geodesics
     *   * crs1 leads away from pt2 along first given geodesic
     *   * crs3 leads away from pt2 along second given geodesic
     */

    /* Now calculate the course change angle: */
    cca = modlon(crs3 - reciprocal(crs1));
    /* Positive cca indicates a right ha	nd turn */

    if (cca > 0)
    {
        rht = 1;
    }

    /* The interior angle between the two geodesics at the vertex */
    vertexAngle = M_PI - fabs(cca);

    /* If vertexAngle is very small, then it is possible that the input radius is
     * too large to fit anywhere between the geodesics.  The conditional here protects
     * against this case.
     * This will also handle the case of collinear geodesics, when pt1 is on the geodesic
     * through pt3 and vice versa. */
    //TODO This check relies on a spherical approximation. Can evaluate true in edge cases.
    B = 0.5 * vertexAngle;
    if (radius > sphereRad * B)
    {
        /* returned dir value will be 0.  Must be checked in calling function */
        err |= NO_TANGENT_ARC_ERR;
        return err;
    }

    /* Dist to start is distance from vertex (pt2) to start point of arc, along crs1 */
    /* Approximate formula from Napier's Rule of Circular Parts (via Janie Henry) */
    distToStart = sphereRad * asin(sin(radius / sphereRad) / sin(B));

    /* For small course change angles, where distToStart is small, azimuth computations are
     * not very accurate when the forward algorithm is initiated from the intersection point.
     * To handle these cases, we'll set a temporary point some larger distance
     * (at least 1 NM) from the intersection point along the incoming azimuth and project from
     * there to find the arc's start and center points.
     * For cases where the distToStart is larger than 10 NM, we'll use the intersection point as our
     * base   */

    if (distToStart < 10.0)
    {
        /* Move to new geodesic start point */
        err |= direct(pt2,crs1,distToStart+1.0,&iterationStart,eps);
        distToStart = 1.0;
        /* recompute azimuth from new start point to intersection point */
        err |= inverse(iterationStart,pt2,&crs1,NULL,NULL,eps);
        /* We'll be approaching arc start point from other direction, so have to switch
         * azimuth to center point */
        revDirFlag = 1;
    }
    else
    {
        iterationStart = pt2;
    }


    distArray[1] = distToStart;
    while ((k == 0) || (err == 0 && ((fabs(error) > tol) || (stepSize > tol))
            && (k < MAX_ITERATIONS) ))
    {

        /* Move to startPoint from pt2 */
        err |= direct(iterationStart, crs1, distToStart, startPoint, eps);
        err |= invCrs(*startPoint, iterationStart, &perpCrs, &tmpCrs, eps);

        /* perpCrs is crs from startPoint to vertex
         * Next we convert it to perpendicular crs by adding/subtracting pi/2
         * depending on turn direction and base point */
        if ((rht && !revDirFlag) || (!rht && revDirFlag))
        {
            perpCrs += M_PI_2;
        }
        else
        {
            perpCrs -= M_PI_2;
        }
        perpCrs = modcrs(perpCrs);

        /* place centerPoint by following perpendicular from startPoint */
        err |= direct(*startPoint, perpCrs, radius, centerPoint, eps);

        /* Project center point onto other geodesic */
        err |= projectToGeo(pt2, crs3, *centerPoint, endPoint, &perpCrs,
                &perpDist, tol, eps);

        error = radius - perpDist;
        if (k == 0)
        {
            if (fabs(error) <= tol)
            {
                /* For some geometries, spherical distToStart is good enough.
                 * In this case, we converge after one step.           */
                break;
            }

            errArray[1] = error;
            adjust = sin(vertexAngle);
            if (revDirFlag)
                adjust = -adjust;

            /* If vertexAngle << 1, then dividing by the sine causes problems */
            if (fabs(adjust) > error)
                distToStart = distToStart + error / adjust;
            else
                //TODO: Make this approximation smarter
                distToStart = distToStart + sgn(adjust)*error;
        }
        else
        {
            errArray[0] = errArray[1];
            distArray[0] = distArray[1];
            errArray[1] = error;
            distArray[1] = distToStart;
            distToStart = findRootSecantMethod(distArray, errArray, &err);
        }
        stepSize = fabs(distArray[1] - distArray[0]);

        k++;

    } //while

    /* Return direction of turn:
     * right hand turn: -1
     * left hand turn: +1       */
    if (rht)
    {
        *dir = CLOCKWISE;
    }
    else
    {
        *dir = COUNTERCLOCKWISE;
    }

    //Now check the iteration to see if it failed
    if (k >= MAX_ITERATIONS)
    {
        err |= ITERATION_MAX_REACHED_ERR;
//        printf("Error: ITERATION_MAX_REACHED in %s\n",__FUNCTION__);

    }
    if (fabs(error) > tol && fabs(error) < MAX_DISTANCE_ERROR)
    {
        err |= SEC_NOT_CONVERGED_ERR;
        return err;
    }
    if (fabs(error) >= MAX_DISTANCE_ERROR)
    {
        err |= ERROR_MAX_REACHED_ERR;
    }

    //    /* Check tangent points */
    //        printf("startPoint %s on crs12\n",ptIsOnCrs(pt2,crs1, *startPoint, NULL, NULL, &err,tol,eps)?"IS":"is NOT");
    //        printf("endPoint %s on crs32\n",ptIsOnCrs(pt2,crs3, *endPoint, NULL, NULL ,&err,tol,eps)?"IS":"is NOT");

    return err;

}

/********************************************************************************
 * Construct arc from start point, start course, and end point.  Solves for radius
 * and center point, returns complete arc struct.
 */
ErrorSet arcFromStartAndEnd(LLPoint arcStart, double startCrs, LLPoint arcEnd, Arc* newArc,
                                 double tol, double eps)
{

    LLPoint arcCenter;
    double radius = 0.0;  // Arc radius to solve for
    double centerEndDist = 0.0;  // Want this to converge to radius
    double dist12 = 0.0;  // Distance from arcStart to arcEnd
    double crs12  = 0.0;  // Course from arcStart to arcEnd
    double crs21  = 0.0;  // Course from arcEnd to arcStart
    double tmpCrs1;
    double angleAtP1 = 0.0; // Angle arcCenter->arcStart->arcEnd
    double turnAngle = 0.0; // Sign gives arc direction: turnAngle < 0 implies right hand turn and clockwise arc.
    double perpCrs   = 0.0; // Course from arcStart to arcCenter, pi/2 from startCrs
    double SPHERE_RAD = SPHERE_RADIUS_NMI;
    double distArray[2] = {0.0, 0.0};
    double errArray[2]  = {1.0e99, 1.0e99};
    double delta = 1.0e99;
    int k = 0;
    ArcDirection arcDir = CLOCKWISE;  // One value chosen for initialization; will not be kept
    ErrorSet err = SUCCESS;
    Arc npArc;

    /* Assign local storage if optional pointers are not provided */
    if (NULL == newArc) newArc = &npArc;

    /* Find course, distance from arcStart to arcEnd */
    err |= inverse(arcStart, arcEnd, &crs12, &crs21, &dist12, eps);

    /* Find magnitude and direction of turn, ignoring the arc and treating the chord
     * from arcStart to arcEnd as the next segment.
     */
    turnAngle = modlon(startCrs - crs12);

    /* Check for zero-radius arc:
     *  First conditional tests for arcEnd coincident with arcEnd
     *  Second conditional tests for extremely small turn or turn very close to 180 degrees */
    if ( (dist12 < tol) || (fabs(sin(turnAngle))*dist12 < tol) )
    {
        newArc->radius      = 0.0;
        newArc->startPoint  = arcStart;
        newArc->endPoint    = arcStart;
        newArc->centerPoint = arcStart;
        newArc->dir = arcDir;           // direction doesn't matter in this case
        newArc->startAz     = 0.0;
        newArc->endAz       = 0.0;
        return err;
    }

    if (turnAngle < 0.0)
    {
        /* Right hand turn */
        perpCrs = modcrs(startCrs + M_PI_2);
    }
    else
    {
        /* Left hand turn */
        arcDir = COUNTERCLOCKWISE;
        perpCrs = modcrs(startCrs - M_PI_2);
    }

    /* Find interior angle at arcStart of triangle formed by arcStart, arcCenter, arcEnd.
     * This angle will be in interval [0,pi/2]. */
    {
        double tmpDist = dist12/SPHERE_RAD;
        angleAtP1 = fabs(M_PI_2 - fabs(turnAngle));
        radius = SPHERE_RAD*atan( (1.0 - cos(tmpDist))/sin(tmpDist)/cos(angleAtP1) );
        /* Formula is bad for sin(tmpDist) =~ 0 or cos(angleAtP1) =~ 0.  These
         * conditions are check above so this formula is protected from them */
    }

    err |= direct(arcStart, perpCrs, radius, &arcCenter, eps);
    err |= inverse(arcCenter, arcEnd,&tmpCrs1,NULL,&centerEndDist,eps);

    distArray[1] = radius;
    errArray[1]  = radius - centerEndDist;  /* error is difference between actual distance and desired radius */

    /* Iterate to refine solution */
    while ( ( ((fabs(errArray[1]) > tol) || (delta > tol)) && (k < MAX_ITERATIONS)  ) || (k == 0) )
    {
        if (k == 0)
        {
            radius = 0.5*(radius+centerEndDist);
        }
        else
        {
            radius = findRootSecantMethod(distArray,errArray,&err);
        }
        err |= direct(arcStart,perpCrs,radius,&arcCenter,eps);
        err |= invDist(arcEnd,arcCenter,&centerEndDist,eps);
        errArray[0] = errArray[1];
        errArray[1] = radius - centerEndDist;
        distArray[0] = distArray[1];
        distArray[1] = radius;
        delta = fabs(distArray[1] - distArray[0]);
        k = k + 1;
    }

    if (fabs(errArray[1]) >= MAX_DISTANCE_ERROR)
    {
        err |= ERROR_MAX_REACHED_ERR;
        return err;
    }
    else if (k == MAX_ITERATIONS)
    {
        err |= ITERATION_MAX_REACHED_ERR;
        return err;
    }

    err |= createArc(newArc,arcCenter,arcStart,arcEnd,arcDir,tol,eps);

    return err;

}


/* Compute arc beginning at arcStart with tangent course at arcStart equal to arcStartCrs, that
 * is tangent at its end point to geodesic with course outCrs at outPoint.  The end point of the arc
 * will not necessarily coincide with outPoint.  This algorithm solves for the arc radius, arc end point,
 * and arc direction.
 */
ErrorSet arcTanToCrs(LLPoint arcStart, double arcStartCrs, LLPoint outPoint,
                                               double outCrs, Arc *newArc, double tol, double eps)
{

    ErrorSet err = SUCCESS;
    double SPHERE_RAD = SPHERE_RADIUS_NMI;
    int k = 0;
    double distArray[2] = {0.0, 0.0};
    double errArray[2]  = {1.0e99, 1.0e99};
    double revStartCrs = modcrs(arcStartCrs + M_PI);
    double revOutCrs   = modcrs(outCrs + M_PI);
    double distD = 0.0;
    double angleGamma = 0.0;
    double angleAlpha = 0.0;
    double startR = 0.0;
    double endR = 0.0;
    double crsToCenter = 0.0;
    double radDiff = 1.0e99;
    double stepSize = 9.0e99;
    ArcDirection turnDir = CLOCKWISE;
    int loopForm = 0;
    Arc npArc;

    /* Assign local storage if optional pointers are not provided */
    if (NULL == newArc) newArc = &npArc;

    /* Compute interior angles of spherical triangle that will be used for approx solution */
    err |= inverse(arcStart,outPoint,&angleGamma,&angleAlpha,&distD,eps);

    /* Check for zero-radius arc: No arc can be constructed if arc start point lies
     * on outbound geodesic. */
    //TODO: Check the following assumption: (MJM)
    /* NOTE: The case where arcStart == outPoint and startCrs == outCrs is also excluded
     * by this test even though two valid arcs could be constructed.  The problem is that
     * choosing which arc to construct requires that a direction be provided, but it is not.  It is
     * ASSUMED that this case is not needed for TERPS constructions. */
     if ( ptIsOnCrs(outPoint,outCrs,arcStart,NULL,NULL,&err,tol,eps) ||
          ptIsOnCrs(outPoint,modcrs(outCrs+M_PI),arcStart,NULL,NULL,&err,tol,eps) )
     {
         newArc->radius      = 0.0;
         newArc->startPoint  = arcStart;
         newArc->endPoint    = arcStart;
         newArc->centerPoint = arcStart;
         newArc->dir = turnDir;           // direction doesn't matter in this case
         newArc->startAz     = 0.0;
         newArc->endAz       = 0.0;
         return err;
     }

    angleGamma = modlon(angleGamma-revStartCrs);
    angleAlpha = modlon(angleAlpha-revOutCrs);

    /* Gamma > 0 => outLine.startPoint to right of this.startPoint
     * Gamma < 0 => outLine.startPoint to left of this.startPoint
     * sign(Alpha) = sign(Gamma) => RF turns directly toward outLine.startPoint (loopForm = 0)
     * sign(Alpha) = -sign(Gamma) => RF turns away from outLine.startPoint and  (loopForm = 1)
     *                              circles around                             */
    if (angleAlpha > 0) turnDir = COUNTERCLOCKWISE;
    if (sgn(angleGamma) == -sgn(angleAlpha)) loopForm = 1;

    distD = distD/SPHERE_RAD;
    /* Spherical Solution */
    {
        double cosDistD = cos(distD);
        double sinA = sin(fabs(angleAlpha));
        double cosA = cos(fabs(angleAlpha));
        double sinG = sin(fabs(angleGamma));
        double cosG = cos(fabs(angleGamma));
        double numer = 0.0;
        double denom = 0.0;
        double beta = 0.0, cosB = 0.0, sinL = 0.0;
        if (! loopForm)
        {
            numer = cosA*sinG-sinA*cosG*cosDistD;
            denom = 1.0 + cosA*cosG+sinA*sinG*cosDistD;
        }
        else
        {
            numer = cosA*sinG + sinA*cosG*cosDistD;
            denom = 1.0 + cosA*cosG - sinA*sinG*cosDistD;
        }
        beta = atan(numer/denom);
        cosB = cos(beta);
        sinL = sinA*sqrt(1.0-cosDistD*cosDistD)/cosB;
        startR = fabs(SPHERE_RAD*atan( (1.0 - sqrt(1.0-sinL*sinL))/sinL/cosB));
    }

    if (startR > MAX_ELLIPSOIDAL_ARC_RADIUS_NMI)
    {
        err |= RADIUS_OUT_OF_RANGE_ERR;
        newArc->radius = startR;
        return err;
    }

    crsToCenter = modcrs(arcStartCrs + (double) turnDir * M_PI_2);
    newArc->startPoint = arcStart;
    err |= direct(arcStart,crsToCenter,startR,&newArc->centerPoint,eps);
    err |= projectToGeo(outPoint,outCrs,newArc->centerPoint,&newArc->endPoint,NULL,
            &endR, tol, eps);
    radDiff = startR - endR;

    k = 0;
    distArray[1] = startR;
    errArray[1] =radDiff;
    while ( (((fabs(radDiff) > tol) || (stepSize > tol)) && (k < MAX_ITERATIONS) ) || (k < 2) ) /* at least 2 iterations */
    {
        if (k == 0)
            startR = 0.5*(startR+endR);
        else
            startR = findRootSecantMethod(distArray,errArray,&err);
        err |= direct(arcStart,crsToCenter,startR,&newArc->centerPoint,eps);
        err |= projectToGeo(outPoint,outCrs,newArc->centerPoint,&newArc->endPoint,NULL,
                &endR, tol, eps);
        radDiff = startR - endR;
        errArray[0] = errArray[1];
        errArray[1] = radDiff;
        distArray[0] = distArray[1];
        distArray[1] = startR;
        stepSize = fabs(distArray[1] - distArray[0]);

        k = k + 1;

    }

    if (radDiff >= MAX_DISTANCE_ERROR)
    {
        err |= ERROR_MAX_REACHED_ERR;
        return err;
    }

    if (k >= MAX_ITERATIONS)
    {
        err |= ITERATION_MAX_REACHED_ERR;
        /* Could still recover from this err (Arc might be valid) */
    }

    err |= createArc(newArc, newArc->centerPoint, newArc->startPoint, newArc->endPoint, turnDir, tol, eps);

    return err;

}

/*******************************************************************************
 * Find points on arc where lines from given point are tangent.
 *
 * INPUT:
 * point:
 *
 *
 *
 */

ErrorSet ptsOnArcOnTanThruPt(LLPoint point, LLPoint center, double radius,
                             LLPointPair tanPt, int* n, double tol, double eps)
{

    ErrorSet err = 0;

    double distToCenter;
    double crsToCenter, crsFromCenter;
    double radCrs, tanCrs;
    double a, b, sinb, C, temp;
    double distError, angularError, diff;

    int i, k;

#ifdef USE_BEST_FIT_ROC

    double sphereRad =
    lookUpROC(geocentricLat(0.5*fabs(center.latitude+point.latitude)));
#else

    double sphereRad = SPHERE_RADIUS;
#endif

	int npN;
	LLPointPair npTanPt;

    /* Assign local storage if optional pointers are not provided */
    if (NULL == n) n = &npN;
    if (NULL == tanPt) tanPt = npTanPt;

    *n = 0;

    //Make sure radius is non-zero & positive
    if (radius <= 0.0)
    {
        err |= RADIUS_OUT_OF_RANGE_ERR;
        return err;
    }

    if (err |= inverse(point, center, &crsToCenter, &crsFromCenter,
            &distToCenter, eps))
        return err;

    if (fabs(distToCenter - radius) < tol)
    {
        /* point is a tangent point -- return it */
        *n = 1;
        tanPt[0] = point;
        return err;
    }
    else if (distToCenter < radius)
    {
        /* point is inside arc--no tangents exist */
        /* returned value of n will be 0 */
        err |= NO_TANGENT_ARC_ERR;
        return err;
    }

    /* Find spherical approximation */
    /* Consider right spherical triangle formed by one point, center, and one tangent
     * point.  Then solve for angle between radius and
     *                            point
     *  C = 90 degrees             |\
     *  a = radius of arc          |B\
     *  c = distToCenter           |  \
     *                            c|   \a
     *                             |    \
     *                             |     \
     *                             |90___C\
     *                          tanPt  b  center
     */

    a = distToCenter / sphereRad; /* convert distance to angle */
    b = radius / sphereRad; /* convert distance to angle */
    sinb = sphereRad * sin(b);
    /* Formula from Napier's Rules of Circular Parts (via Janie Henry) */
    C = acos(tan(b) / tan(a));

    *n = 2;
    //  tanPt[0] = direct(center,modcrs(crsFromCenter+C),radius,eps);
    // Do following after finding tanPt[0]
    // tanPt[1] = direct(center,modcrs(crsFromCenter-C),radius,eps);

    /* Now we refine the approximations */

    for (i = 0; i < *n; i++)
    {

        k = 0;

        while ((k == 0) || ((fabs(distError) > tol) && (k < MAX_ITERATIONS) ))
        {

            //      printf("Finding Tan Points: i = %d, k = %d\n",i,k);

            if (i == 0)
            {
                if (err |= direct(center, modcrs(crsFromCenter + C), radius,
                        &tanPt[i], eps))
                    return err;
            }
            else
            {
                if (err |= direct(center, modcrs(crsFromCenter - C), radius,
                        &tanPt[i], eps))
                    return err;
            }

            if (err |= invCrs(tanPt[i], center, &radCrs, &temp, eps))
                return err;
            if (err |= invCrs(tanPt[i], point, &tanCrs, &temp, eps))
                return err;

            diff = modlon(radCrs - tanCrs); /* want this to equal 90 degrees */

            angularError = fabs(diff) - M_PI_2;

            C = C + angularError;

            distError = sinb * angularError;

            k++;

        }

        if (k >= MAX_ITERATIONS)
        {
            err |= ITERATION_MAX_REACHED_ERR;
//            printf("Error: ITERATION_MAX_REACHED in %s\n",__FUNCTION__);

        }

        if (fabs(distError) >= MAX_DISTANCE_ERROR)
        {
            err |= ERROR_MAX_REACHED_ERR;
        }

    }

    return err;

}

/*******************************************************************************
 * Find arc of given radius tangential to two loci
 * dir = -1 for counter clockwise arc
 * dir = +1 for clockwise arc
 */

ErrorSet arcTanToTwoLoci(Locus loc1, Locus loc2, double radius,
                                 LLPoint* centerPoint, LLPoint* startPoint,
                                 LLPoint* endPoint, ArcDirection* dir,
                                 double tol, double eps)
{
    ErrorSet err = 0;
    //int tmpTest;
    LLPoint intx[3];
    LLPoint geoPt1;
    LLPoint tempPt;
    LLPoint sphIntx0, sphIntx1, sphIntx2;

    double distbase = 0.0;
    double crs12, crs21, crs23, crs32, lcrs1;
    double crsFromPoint, distFromPoint;
    double gcrs1, geoLen1;
    double rcrs1, rcrs2, temp, vertexAngle;
    double locAngle;
    double error;
    double azDiff;
    double r2;
    double adjust = 0.0;
    double distarray[2] = { 0, 0 };
    double errarray[2] = { 0, 0 };
    double stepSize = 9.0e99;

    int k = 0;
#ifdef USE_BEST_FIT_ROC

    double sphereRad =
    lookUpROC(
            geocentricLat(
                    0.25*fabs(loc1.locusStart.latitude +
                            loc1.locusEnd.latitude +
                            loc2.locusStart.latitude +
                            loc2.locusEnd.latitude)
            )
    );
#else

    double sphereRad = SPHERE_RADIUS;
#endif

    LLPoint npCenterPoint;
    LLPoint npStartPoint;
    LLPoint npEndPoint;
    ArcDirection npDir;

    /* Assign local storage if optional pointers are not provided */
    if (NULL == centerPoint) centerPoint = &npCenterPoint;
    if (NULL == startPoint) startPoint = &npStartPoint;
    if (NULL == endPoint) endPoint = &npEndPoint;
    if (NULL == dir) dir = &npDir;

    err |= invCrs(loc1.locusStart, loc1.locusEnd, &crs12, &crs21, eps);
    err |= inverse(loc1.geoStart, loc1.geoEnd, &gcrs1, &temp, &geoLen1,
            eps);

    err |= invCrs(loc2.locusStart, loc2.locusEnd, &crs23, &crs32, eps);

    crs32 = modcrs(crs32 + M_PI);

    /* Use geodesic approximation of loci to calculate approximate arc */
    err |= arcTanToTwoGeos(loc1.locusStart, crs12, loc2.locusEnd,
            crs32, radius, &intx[0], &intx[1], &intx[2], dir, tol, eps);
    /* intx[0] will be approx centerPt
     * intx[1] will be near loc1
     * intx[2] will be near loc2 */

    /* If no approx arc found, then CANNOT CONTINUE, return null */
    if (err)
        return err;

    /* Calculate the angle at the vertex where loc1 and loc2 intersect
     * This is a spherical calculation based on geodesic approximation of loci */
    sphIntx0 = geodeticToGeocentric(intx[0]);
    sphIntx1 = geodeticToGeocentric(intx[1]);
    sphIntx2 = geodeticToGeocentric(intx[2]);
    rcrs1 = sphereInvCrs(sphIntx0, sphIntx1, eps);
    rcrs2 = sphereInvCrs(sphIntx0, sphIntx2, eps);

    vertexAngle = fabs(modlon(rcrs1 - rcrs2));
    vertexAngle = 2.0 * acos(sin(vertexAngle / 2.0) * cos(radius / sphereRad));

    /* Calculate inclination angle of loc1 relative to its geodesic */
    locAngle = atan((loc1.endDist - loc1.startDist) / geoLen1);

    //Calculate distance from geoStart to projection of intx[1] on loc1's geodesic
    if (err |= projectToGeo(loc1.geoStart, gcrs1, intx[1], &geoPt1,
            &crsFromPoint, &distFromPoint, tol, eps))
        return err;

    if (err |= invDist(loc1.geoStart, geoPt1, &distbase, eps))
        return err;

    err |= invCrs(loc1.geoStart, geoPt1, &crs12, &crs21, eps);

    err |= minSubtendedAngle(crs12, gcrs1, &azDiff);

    if (azDiff > M_PI/2) {
    	distbase = distbase * -1;
    }

    /* NOTE: error not set, relying on k==0 requirement to get into first iteration */
    while ((k == 0) || ((k < MAX_ITERATIONS) && ((fabs(error) > tol) || (stepSize > tol)) ))
    {
        if (k > 0)
        {
            /* after first iteration, must project point on geodesic */
            if (err |= direct(loc1.geoStart, gcrs1, distbase, &geoPt1, eps))
                return err;
        }

        if (err |= ptOnLocusFromGeoPt(loc1, geoPt1, &intx[1], &temp, tol, eps))
            return err;

        /* tempPt and &temp are required by function, but we don't need the data they'll
         * be storing */

        lcrs1 = locusCrsAtPt(loc1, intx[1], &tempPt, &temp, &err, tol,
                eps);
        if (lcrs1 < 0.0)
        {
            //locusCrsAtPt() failed
            return err |= INVALID_CRS_ERR;
        }

        /* Convert locus course to course toward arc center */
        lcrs1 = lcrs1 + (*dir) * M_PI_2;

        /* Updated arc center position */
        if (err |= direct(intx[1], lcrs1, radius, &intx[0], eps))
            return err;

        /* Project center point to locus2 */
        if (err |= projectToLocus(loc2, intx[0], &tempPt, &temp, &r2,
                tol, eps))
            return err;

        intx[2] = tempPt;
        /* Compare actual radius to desired radius */
        error = r2 - radius;

        if (fabs(error) <= tol)
            break; //No need to continue. Get out of loop now.

        distarray[0] = distarray[1];
        distarray[1] = distbase;
        errarray[0] = errarray[1];
        errarray[1] = error;

        stepSize = fabs(distarray[1] - distarray[0]);

        if (k == 0)
        {
            /* On first iteration, use geometric construction for next approximation */
            /* This is necessary because distarray and errarray are not full yet */
            adjust = cos(locAngle) / sin(vertexAngle);
            if (adjust < error)
            {
                distbase = distbase + error;
            }
            else
            {
                distbase = distbase + error * adjust;
            }
        }
        else
        {
            /* distarray and errarray are full, so use linear interpolation
             * to improve approx */
            distbase = findRootSecantMethod(distarray, errarray, &err);
        }

        k++;

    }

    if (fabs(error) >= MAX_DISTANCE_ERROR)
    {
        err |= ERROR_MAX_REACHED_ERR;
    }

    if (k >= MAX_ITERATIONS)
    {
        err |= ITERATION_MAX_REACHED_ERR;
    }

//    printf("plot(%f, %f, 'Marker', 'o', 'Color', 'b')\n", loc1.locusStart.longitude * 180 / M_PI, loc1.locusStart.latitude * 180 / M_PI);
//    printf("plot(%f, %f, 'Marker', 'o', 'Color', 'g')\n", loc2.locusStart.longitude * 180 / M_PI, loc2.locusStart.latitude * 180 / M_PI);
//
//    printf("plot(%f, %f, 'Marker', 'o', 'Color', 'k')\n", intx[0].longitude * 180 / M_PI, intx[0].latitude * 180 / M_PI);
//    printf("plot(%f, %f, 'Marker', 'o', 'Color', 'r')\n", intx[1].longitude * 180 / M_PI, intx[1].latitude * 180 / M_PI);
//    printf("plot(%f, %f, 'Marker', 'o', 'Color', 'r')\n", intx[2].longitude * 180 / M_PI, intx[2].latitude * 180 / M_PI);

	/* Copy solution to output addresses */
	*centerPoint = intx[0];
	*startPoint = intx[1];
	*endPoint = intx[2];


    return err;

}

static double outsideTanPtsApprox(double dist, double r1, double r2)
{
    /* Required: r1 >= r2 */
    /* Convert to angular distances */
    double delta = dist / SPHERE_RADIUS;
    double rho1 = r1 / SPHERE_RADIUS;
    double rho2 = r2 / SPHERE_RADIUS;
    double xi = atan(sin(delta) * sin(rho2) / (sin(rho1) - sin(rho2) * cos(
            delta)));
    /* return phi1: azimuth from center1 to approx. tan. point */
    return acos(tan(rho1) / tan(delta + xi));
}

static double crossingTanPtsApprox(double dist, double r1, double r2)
{
    /* Required: r1 >= r2 */
    /* Convert to angular distances */
    double delta = dist / SPHERE_RADIUS;
    double rho1 = r1 / SPHERE_RADIUS;
    double rho2 = r2 / SPHERE_RADIUS;
    double xi = atan(sin(delta) * sin(rho2) / (sin(rho1) + sin(rho2) * cos(
            delta)));
    /* return phi1: azimuth from center1 to approx. tan. point */
    return acos(tan(rho1) / tan(delta - xi));
}

/* MJM's implementation
 *
 * */
ErrorSet geoTanToTwoCircles(LLPoint center1, double r1, ArcDirection dir1,
                              LLPoint center2, double r2, ArcDirection dir2,
                              Geodesic tanLines[2], double tol, double eps)
{

    LLPoint tempCenter;
    LLPoint startPoints[2], endPoints[2];
    Geodesic tempGeo;

    double dist12, az12, az21; /* Distance, course from center1 to center2 */
    double tempR;
    ArcDirection tempDir;
    double error, errarray[2], azarray[2];
    double angle;
    double rightAngle = M_PI_2;
    double daz, newaz, az, baz;
    double azFromPoint, distFromPoint;
    double dist;

    ErrorSet err = 0;

    int swapped = 0;
    int i, k;

    Geodesic npTanLines[2];

    /* Assign local storage if optional pointers are not provided */
    if (NULL == tanLines) tanLines = npTanLines;

    /* Calculation assumes r1 >= r2.
     * If not true, swap circles before starting*/
    if (r1 < r2)
    {
        tempCenter = center1;
        tempR = r1;
        center1 = center2;
        r1 = r2;
        center2 = tempCenter;
        r2 = tempR;
        swapped = 1; /* Set flag so we remember swap has occurred */
        tempDir = dir1;
        dir1 = dir2;
        dir2 = tempDir;
    }

    /* Calculate distance, az from center1 to center2 */
    err |= inverse(center1, center2, &az12, &az21, &dist12, eps);
    if (err)
        return err;

    //    printf("az,dist = %f, %f\n",az12, dist12);

    /* Test for case when one circle is inside other: no tangent lines
     * NOTE: dist12 + r2 == r1 corresponds to tangent circles, we return
     * no solution in this case, since we're looking for tangent LINES
     * For opposite-orientation case, smaller circle must be completely
     * outside of larger circle */
    if (((dir1 == dir2) && (dist12 + r2 - r1 < tol)) || ((dir1 != dir2)
            && (dist12 - r2 - r1 < tol)))
    {
        err |= CIRCLE_INSIDE_CIRCLE_ERR;
        return err;
    }

    /* Solution depends on orientation of circles:
     * dir1 == dir2 => tangent lines do not cross center-center line
     * dir1 == -dir2 => tangent lines cross center-center line. */
    if (dir1 == dir2)
    {
        angle = outsideTanPtsApprox(dist12, r1, r2);
        //        printf("angle = %f\n",angle);
    }
    else
    {
        angle = crossingTanPtsApprox(dist12, r1, r2);
        //        printf("angle = %f\n",angle);
    }

    /* Set up angles so that first start point found is starting on circle 1
     * (based on swapped circles, if swapping has occurred) */
    if (dir1 == COUNTERCLOCKWISE)
    {
        //        printf("Negating angles\n");
        angle = -angle;
        rightAngle = -rightAngle;
    }

    /* Find starting point and then refine.  There are two possible
     * solutions, so we go through this loop twice. */

    azarray[1] = az12 - angle;
    for (i = 0; i < 2; i++)
    {

        k = 0;

        while ((k == 0) || ((fabs(error) > tol) && (k < MAX_ITERATIONS)
                ))
        {

            /* Plot estimate of startPoint on first circle */
            err |= direct(center1, azarray[1], r1, &startPoints[i], eps);
            if (err)
                return err;
            //            printf("startPoint[%d] = ",i);
            //            _display(startPoints[i]);


            /* Calculate azimuth from startPoint to center1 */
            err |= invCrs(startPoints[i], center1, &az, &baz, eps);
            if (err)
                return err;

            /* Calculate azimuth of tangent line at startPoint */
            az = az - rightAngle;

            /* Project center2 to tangent line.  This is approximately tangent point to second circle */
            err |= projectToGeo(startPoints[i], az, center2,
                    &endPoints[i], &azFromPoint, &distFromPoint, tol, eps);
            if (err)
                return err;

            //            printf("Approx end point = ");
            //            _display(endPoints[i]);

            /* Calculate error = distancd from tangent point to circle */
            error = distFromPoint - r2;

            //            printf("error = %e\n",error);

            errarray[1] = error;

            if (k == 0)
            {
                /* First improvement uses spherical triangle, since only one sample exists  */
                err |= invDist(endPoints[i], startPoints[i], &dist, eps);
                if (err)
                    return err;
                daz = asin(sin(error / SPHERE_RADIUS) / sin(dist
                        / SPHERE_RADIUS));
                newaz = azarray[1] - daz;
            }
            else
            {
                /* Once two approximations have been found, use secant method */
                newaz = findRootSecantMethod(azarray, errarray, &err);
            }

            errarray[0] = errarray[1];
            azarray[0] = azarray[1];

            azarray[1] = newaz;

            k++;

        }

        if (k >= MAX_ITERATIONS)
        {
            err |= ITERATION_MAX_REACHED_ERR;
//            printf("Error: ITERATION_MAX_REACHED in %s\n",__FUNCTION__);

        }

        if (fabs(error) >= MAX_DISTANCE_ERROR)
        {
            err |= ERROR_MAX_REACHED_ERR;
        }

        /* Create Geodesic for output */
        if (i == 0)
        {
            err |= createGeo(&tanLines[i], startPoints[i], endPoints[i],
                    SEGMENT, eps);
            if (err)
                return err;
        }
        else
        {
            err |= createGeo(&tanLines[i], endPoints[i], startPoints[i],
                    SEGMENT, eps);
            if (err)
                return err;
        }

        /* Make angle negative to find second solution */
        angle = az12 - azarray[1];
        azarray[1] = az12 + angle;

        rightAngle = -rightAngle;

    }

    if (swapped)
    {
        /* Swap lines so that first tan line returned always has start point on
         * first circle passed to this function */
        tempGeo = tanLines[0];
        tanLines[0] = tanLines[1];
        tanLines[1] = tempGeo;
    }

    return err;

}

static int computeSide(double crs1, double crs2){

	double angle;

	angle = computeSubtendedAngle(crs1, crs2, CLOCKWISE);

	if (angle < -M_PI){
		return 1; //Right
	} else if (angle < 0.0){
		return -1; //Left
	} else if (angle < M_PI){
		return 1; //Right
	} else {
		return -1;
	}

}

ErrorSet geoTanToArcAtAngleToGeo(Arc arc, Geodesic geo, double angle, Geodesic* tangentLine, int* tangentLineLocation, double tol, double eps){

	ErrorSet err = 0;

	LLPoint anglePt, tanPt;
	LLPoint centerProjToGeoPt;

	double crsFromGeoStartToCenter, crsFromGeoStartToTanPt, crsFromGeoStartToAnglePt;
	double crsFromCenterToGeo, crsFromCenterToTanPt;
	double crsFromTanPtToAnglePt, crsFromTanPtToCenter;
	double crsFromAnglePtToGeoStart, crsFromAnglePtToTanPt;
        double crsFromCenterToGeoEnd, crsFromGeoEndToCenter;

	double distFromCenterToTanPt, distFromProjCenterPt;

	double errorArray[2], distArray[2];
	double dist, centerToGeoDist;

	double minAcceptableAngle, maxAcceptableAngle;

	int LEFT = -1, RIGHT = 1;
	int centerSide, tanPtSide;
	int sgnErr;
	double delta = 9e99;

	int k = 0;

	LLPoint tempPt;
	Geodesic npTangentLine;
	int npTangentLineLocation;

	if(tangentLine == NULL) tangentLine = &npTangentLine;
	if(tangentLineLocation == NULL) tangentLineLocation = &npTangentLineLocation;


	/* Check bounds on angle such that solution exists */
	err |= projectToGeo(geo.startPoint, geo.startAz, arc.centerPoint, &centerProjToGeoPt, NULL, &distFromProjCenterPt, tol, eps);
	err |= inverse(arc.centerPoint, centerProjToGeoPt, &crsFromCenterToGeo, NULL, &centerToGeoDist, eps);

	err |= inverse(geo.startPoint, arc.centerPoint, &crsFromGeoStartToCenter, NULL, NULL, eps);
	centerSide = computeSide(geo.startAz, crsFromGeoStartToCenter);

	if( centerSide == LEFT ){
		minAcceptableAngle = (centerToGeoDist + arc.radius)/SPHERE_RADIUS_NMI;
		maxAcceptableAngle = M_PI - (fabs(centerToGeoDist - arc.radius))/SPHERE_RADIUS_NMI;
	} else {
		minAcceptableAngle = (fabs(centerToGeoDist - arc.radius))/SPHERE_RADIUS_NMI;
		maxAcceptableAngle = M_PI - (centerToGeoDist + arc.radius)/SPHERE_RADIUS_NMI;
	}
	if( !((angle > minAcceptableAngle) && (angle < maxAcceptableAngle)) ){
		/* No Solution Can Exist */
		*tangentLineLocation = 0;
		return err;
	}

        if (distFromProjCenterPt <= tol)
          centerSide = 0;

	//	Flat Earth Approximation
	//	err |= projectToGeo(geo.startPoint, geo.startAz, arc.centerPoint, &centerProjToGeoPt, &crsFromCenterToGeo, &d, tol, eps);
	//	distArray[0] = (arc.radius/cos(angle) + d)/tan(angle);

	// First Initial Guess

		// Find the geometric case
	if(arc.dir == CLOCKWISE && centerSide == LEFT){
		sgnErr = -1;
		crsFromCenterToTanPt = crsFromCenterToGeo + M_PI + angle;
		err |= direct(arc.centerPoint, crsFromCenterToTanPt, arc.radius, &tanPt, eps);
		err |= inverse(geo.startPoint, tanPt, &crsFromGeoStartToTanPt, NULL, NULL, eps);
		tanPtSide = computeSide(geo.startAz, crsFromGeoStartToTanPt);
		if(centerSide != tanPtSide){
			sgnErr *= -1;
			crsFromCenterToTanPt = crsFromCenterToGeo + M_PI - angle;
			err |= direct(arc.centerPoint, crsFromCenterToTanPt, arc.radius, &tanPt, eps);
		}
	} else if(arc.dir == CLOCKWISE && centerSide == RIGHT){
		sgnErr = 1;
		crsFromCenterToTanPt = crsFromCenterToGeo - angle;
		err |= direct(arc.centerPoint, crsFromCenterToTanPt, arc.radius, &tanPt, eps);
		err |= inverse(geo.startPoint, tanPt, &crsFromGeoStartToTanPt, NULL, NULL, eps);
		tanPtSide = computeSide(geo.startAz, crsFromGeoStartToTanPt);
		if(centerSide != tanPtSide){
			sgnErr *= -1;
			crsFromCenterToTanPt = crsFromCenterToGeo + angle;
			err |= direct(arc.centerPoint, crsFromCenterToTanPt, arc.radius, &tanPt, eps);
		}
	} else if(arc.dir == COUNTERCLOCKWISE && centerSide == LEFT){
		sgnErr = 1;
		crsFromCenterToTanPt = crsFromCenterToGeo + angle;
		err |= direct(arc.centerPoint, crsFromCenterToTanPt, arc.radius, &tanPt, eps);
		err |= inverse(geo.startPoint, tanPt, &crsFromGeoStartToTanPt, NULL, NULL, eps);
		tanPtSide = computeSide(geo.startAz, crsFromGeoStartToTanPt);
		if(centerSide != tanPtSide){
			sgnErr *= -1;
			crsFromCenterToTanPt = crsFromCenterToGeo - angle;
			err |= direct(arc.centerPoint, crsFromCenterToTanPt, arc.radius, &tanPt, eps);
		}
	} else if(arc.dir == COUNTERCLOCKWISE && centerSide == RIGHT){
		sgnErr = -1;
		crsFromCenterToTanPt = crsFromCenterToGeo + M_PI - angle;
		err |= direct(arc.centerPoint, crsFromCenterToTanPt, arc.radius, &tanPt, eps);
		err |= inverse(geo.startPoint, tanPt, &crsFromGeoStartToTanPt, NULL, NULL, eps);
		tanPtSide = computeSide(geo.startAz, crsFromGeoStartToTanPt);
		if(centerSide != tanPtSide){
			sgnErr *= -1;
			crsFromCenterToTanPt = crsFromCenterToGeo + M_PI + angle;
			err |= direct(arc.centerPoint, crsFromCenterToTanPt, arc.radius, &tanPt, eps);
		}
	}
        else if (centerSide == 0) {
          err |= invCrs(arc.centerPoint, geo.endPoint, &crsFromCenterToGeoEnd, &crsFromGeoEndToCenter, eps);
          if (arc.dir == CLOCKWISE) {
            sgnErr = -1;
            tanPtSide = LEFT;
            crsFromCenterToTanPt = modpos(crsFromCenterToGeoEnd - M_PI_2, M_2PI);
            crsFromCenterToTanPt = modpos(crsFromCenterToTanPt + angle, M_2PI);
          }
          else if (arc.dir == COUNTERCLOCKWISE) {
            sgnErr = 1;
            tanPtSide = RIGHT;
            crsFromCenterToTanPt = modpos(crsFromCenterToGeoEnd + M_PI_2, M_2PI);
            crsFromCenterToTanPt = modpos(crsFromCenterToTanPt - angle, M_2PI);
          }
          err |= direct(arc.centerPoint, crsFromCenterToTanPt, arc.radius, &tanPt, eps);
        }

	err |= inverse(arc.centerPoint, tanPt, NULL, &crsFromTanPtToCenter, NULL, eps);
	crsFromTanPtToAnglePt = crsFromTanPtToCenter - arc.dir*M_PI_2;

	err |= crsIntx(geo.startPoint, geo.startAz, NULL, NULL, tanPt, crsFromTanPtToAnglePt, NULL, NULL, &anglePt, tol, eps);

	err |= inverse(geo.startPoint, anglePt, &crsFromGeoStartToAnglePt, &crsFromAnglePtToGeoStart, &dist, eps);
	if(dist < 10.0){
		/* Move the geo.startPoint if it starts close to the anglePt */
		err |= direct(geo.startPoint, geo.startAz, 20.0, &tempPt, eps);
		err |= inverse(geo.startPoint, tempPt, NULL, &geo.startAz, NULL, eps);
		geo.startAz += M_PI;
		geo.startPoint = tempPt;
		err |= crsIntx(geo.startPoint, geo.startAz, NULL, NULL, tanPt, crsFromTanPtToAnglePt, NULL, NULL, &anglePt, tol, eps);
		err |= inverse(geo.startPoint, anglePt, &crsFromGeoStartToAnglePt, &crsFromAnglePtToGeoStart, &dist, eps);
	}

	// Give a sign to the distance in case the anglePt is behind geo.StartPt
	double angleDirTest;
	minSubtendedAngle(crsFromGeoStartToAnglePt,geo.startAz,&angleDirTest);
	if(angleDirTest > M_PI_2){
		distArray[0] = -1*dist;
	} else {
		distArray[0] = dist;
	}

	crsFromAnglePtToTanPt = crsFromAnglePtToGeoStart - tanPtSide*angle;
	double newCrsFromCenterToTanPt;
	err |= projectToGeo(anglePt, crsFromAnglePtToTanPt, arc.centerPoint, &tanPt, &newCrsFromCenterToTanPt, &distFromCenterToTanPt, tol, eps);
	minSubtendedAngle(newCrsFromCenterToTanPt,crsFromCenterToTanPt,&angleDirTest);
	if(angleDirTest > M_PI_2){
		errorArray[0] = -distFromCenterToTanPt - arc.radius;
	} else {
		errorArray[0] = distFromCenterToTanPt - arc.radius;
	}
	// Second Initial Guess
	distArray[1] = distArray[0] + sgnErr*errorArray[0]/sin(angle);

	err |= direct(geo.startPoint, geo.startAz, distArray[1], &anglePt, eps);
	err |= inverse(geo.startPoint, anglePt, NULL, &crsFromAnglePtToGeoStart, NULL, eps);
	crsFromAnglePtToTanPt = crsFromAnglePtToGeoStart - tanPtSide*angle;
	err |= projectToGeo(anglePt, crsFromAnglePtToTanPt, arc.centerPoint, &tanPt, &newCrsFromCenterToTanPt, &distFromCenterToTanPt, tol, eps);
	errorArray[1] = distFromCenterToTanPt - arc.radius;

	while((k == 0 || fabs(errorArray[1]) > tol  || delta > tol) && k < MAX_ITERATIONS){

		// Find new distance
		dist = findRootSecantMethod(distArray, errorArray, &err);
		distArray[0] = distArray[1];
		distArray[1] = dist;
		delta = fabs(distArray[1] - distArray[0]);

		// Find new error
		err |= direct(geo.startPoint, geo.startAz, distArray[1], &anglePt, eps);
		err |= inverse(geo.startPoint, anglePt, NULL, &crsFromAnglePtToGeoStart, NULL, eps);
		crsFromAnglePtToTanPt = crsFromAnglePtToGeoStart - tanPtSide*angle;
		err |= projectToGeo(anglePt, crsFromAnglePtToTanPt, arc.centerPoint, &tanPt, &newCrsFromCenterToTanPt, &distFromCenterToTanPt, tol, eps);

		errorArray[0] = errorArray[1];
		minSubtendedAngle(newCrsFromCenterToTanPt,crsFromCenterToTanPt,&angleDirTest);
		if(angleDirTest > M_PI_2){
			errorArray[1] = -distFromCenterToTanPt - arc.radius;
		} else {
			errorArray[1] = distFromCenterToTanPt - arc.radius;
		}
		k++;
	}

	createGeo(tangentLine,tanPt,anglePt,SEGMENT,eps);

	*tangentLineLocation = tanPtSide;

	if(k == MAX_ITERATIONS) err = ITERATION_MAX_REACHED_ERR;

	return err;
}

ErrorSet arcTanToArcAndGeo(Arc arc, Geodesic geo, double tangentArcRadius, Arc* tangentArc, int* tangentArcLocation, double tol, double eps)
{
	ErrorSet err = 0;

	LLPoint centerProjToGeo;
	int tanArcLoc = 0;//-1 = tangent arc strictly to the left of the geodesic, 0 = no solution exists, 1 = tangent arc strictly to the right the geodesic
    double distArray[2] = { 0.0, 0.0 };//This is the array for the x variable of the error function, where x represents the start true course from the circle center to the tangent arc center point.  Units are in radians
    double errorArray[2] = { 0.0, 0.0 };//This is the array that holds error function f(x) evaluated at a given x value.  Call the contents of this array y = f(x), where y represents the magnitude of the difference between the tangent arc tangentArcRadius and the distance from the the tangent arc center to the projection of the tangent arc center onto infinite geodesic.  Units are in nautical miles
    double error;//The magnitude of the difference between the tangent arc tangentArcRadius and the distance from the the tangent arc center to the projection of the tangent arc center onto infinite geodesic.  Units are in nautical miles
    double delta;
    double dist;
    double a, b, c, x; //Sides for spherical triangle
    double psi, dx;
    double R = SPHERE_RADIUS_NMI;
    double crsFromCtrProjToGeoToGeoStart, crsFromGeoStartToCtrProjToGeo;
    double distFromArcCtrToGeo;
    double angle;
    double crsFromGeoToArcCtr;
    int side;
    int LEFT = -1, RIGHT = 1, CENTER = 0;
    double crs, crsFromArcCtrToGeo;
    LLPoint tanArcProjOnGeo;

    LLPoint tanArcCenterPt;
    LLPoint tanArcStartPt, tanArcEndPt;
    ArcDirection tanArcDir;

    int k = 0;//iteration counter

    Arc npTangentArc;
    int npTangentArcLocation;

    /* Assign local storage if optional pointers are not provided */
    if (NULL == tangentArc) tangentArc = &npTangentArc;
    if (NULL == tangentArcLocation) tangentArcLocation = &npTangentArcLocation;

    err |= projectToGeo(geo.startPoint, geo.startAz, arc.centerPoint, &centerProjToGeo, &crsFromArcCtrToGeo, &distFromArcCtrToGeo, tol, eps);

    // Check for Solution Existence
    if((distFromArcCtrToGeo - (arc.radius + 2.0*tangentArcRadius)) > tol){
    	// No Solution Exists
    	*tangentArcLocation = 0;
    	return err;
    }

    // Move geo.startPoint to centerProjToGeo
    err |= invCrs(centerProjToGeo, geo.startPoint, &crsFromCtrProjToGeoToGeoStart, &crsFromGeoStartToCtrProjToGeo, eps);
    minSubtendedAngle(geo.startAz,crsFromGeoStartToCtrProjToGeo,&angle);
    if( angle < M_PI_2){
    	crs = modcrs(crsFromCtrProjToGeoToGeoStart + M_PI);
    } else if(angle > M_PI_2){
    	crs = crsFromCtrProjToGeoToGeoStart;
    }
    geo.startPoint = centerProjToGeo;
    geo.startAz = crs;

    // Determine which side the arc is on
    err |= invCrs(geo.startPoint, arc.centerPoint, &crsFromGeoToArcCtr, NULL, eps);
    if(distFromArcCtrToGeo < tol){
    	side = CENTER;
    } else {
    	side = computeSide(geo.startAz, crsFromGeoToArcCtr);
    }

    // Determine which side the tangent arc should go on
    if(arc.dir == CLOCKWISE && side == LEFT){
    	tanArcLoc = LEFT;
    }
    if(arc.dir == CLOCKWISE && side == CENTER){
    	tanArcLoc = LEFT;
    }
    if(arc.dir == CLOCKWISE && side == RIGHT && (distFromArcCtrToGeo - arc.radius < -tol)){
    	tanArcLoc = LEFT;
    }
    if(arc.dir == CLOCKWISE && side == RIGHT && (distFromArcCtrToGeo - arc.radius >= -tol)){
    	tanArcLoc = RIGHT;
    }
    if(arc.dir == COUNTERCLOCKWISE && side == RIGHT){
    	tanArcLoc = RIGHT;
    }
    if(arc.dir == COUNTERCLOCKWISE && side == CENTER){
    	tanArcLoc = RIGHT;
    }
    if(arc.dir == COUNTERCLOCKWISE && side == LEFT && (distFromArcCtrToGeo < arc.radius)){
    	tanArcLoc = RIGHT;
    }
    if(arc.dir == COUNTERCLOCKWISE && side == LEFT && (distFromArcCtrToGeo > arc.radius)){
    	tanArcLoc = LEFT;
    }

    // Determine the tangent arc direction. Should be consistent with geo direction.
    if(tanArcLoc == LEFT){
    	tanArcDir = COUNTERCLOCKWISE;
    } else {
    	tanArcDir = CLOCKWISE;
    }

	// Check for tangent Case
    if(fabs(distFromArcCtrToGeo - (arc.radius + 2.0*tangentArcRadius)) <= tol ){
    	err |= direct(arc.centerPoint, crsFromArcCtrToGeo, arc.radius, &tanArcStartPt, eps);
    	err |= direct(arc.centerPoint, crsFromArcCtrToGeo, (arc.radius + distFromArcCtrToGeo)/2.0, &tanArcCenterPt, eps);
    	tanArcEndPt = centerProjToGeo;
    	err |= createArc(tangentArc, tanArcCenterPt, tanArcStartPt, tanArcEndPt, tanArcDir, tol, eps);
    	*tangentArcLocation = tanArcLoc;
    	return err;
    }

    /* First Approximation */
    // Using the Spherical Law of Cosines
	a = distFromArcCtrToGeo;
	b = tangentArcRadius;
	if(tanArcLoc == side){
		a = M_PI_2 - a/R;
	} else {
		a = M_PI_2 + a/R;
	}
	b = M_PI_2 - b/R;
	c = arc.radius + tangentArcRadius;
	x = R*acos((cos(c/R) - cos(a)*cos(b))/(sin(a)*sin(b)));

	// Place tanArcCenterPt
    dist = x;
    err |= direct(geo.startPoint, geo.startAz, dist, &tanArcProjOnGeo, eps);
    err |= invCrs(geo.startPoint, tanArcProjOnGeo , NULL, &crs, eps);
    crs = fmod(crs + M_PI + tanArcLoc*M_PI_2, M_2PI);
    err |= direct(tanArcProjOnGeo, crs, tangentArcRadius, &tanArcCenterPt, eps);
    distArray[0] = dist;

    // Calculate Error
    err |= invDist(arc.centerPoint, tanArcCenterPt, &error, eps);
    error = error - c;
    errorArray[0] = error;

    /* Second Approximation */
    // Using the Spherical Law of Cosines
	psi = acos((cos(b) - cos(a)*cos(c/R))/(sin(a)*sin(c/R)));
	dx = R*atan(sin(psi)*tan(error/R));

	// Place tanArcCenterPt
    dist = dist + dx;
    err |= direct(geo.startPoint, geo.startAz, dist, &tanArcProjOnGeo, eps);
    err |= invCrs(geo.startPoint, tanArcProjOnGeo , NULL, &crs, eps);
    crs = fmod(crs + M_PI + tanArcLoc*M_PI_2, M_2PI);
    err |= direct(tanArcProjOnGeo, crs, tangentArcRadius, &tanArcCenterPt, eps);
    distArray[1] = dist;

    // Calculate Error
    err |= invDist(arc.centerPoint, tanArcCenterPt, &error, eps);
    error = error - c;
    errorArray[1] = error;

    // Calculate Step Size
    delta = distArray[1] - distArray[0];

    /* Iteration */
    while ( ((k == 0) || (fabs(error) > tol) || (fabs(delta) > tol))
    		&& (k < MAX_ITERATIONS)
           )
    {

        //solve for the root of the error function and use this as the new start true course from
        //the circle center to the tangent arc center point in the next iteration
        dist = findRootSecantMethod(distArray, errorArray, &err);

    	// Place tanArcCenterPt
        err |= direct(geo.startPoint, geo.startAz, dist, &tanArcProjOnGeo, eps);
        err |= invCrs(geo.startPoint, tanArcProjOnGeo , NULL, &crs, eps);
        crs = fmod(crs + M_PI + tanArcLoc*M_PI_2, M_2PI);
        err |= direct(tanArcProjOnGeo, crs, tangentArcRadius, &tanArcCenterPt, eps);

        // Calculate Error
        err |= invDist(arc.centerPoint, tanArcCenterPt, &error, eps);
        error = error - c;


        //update the error function arrays
        distArray[0] = distArray[1];
        errorArray[0] = errorArray[1];
        distArray[1] = dist;
        errorArray[1] = error;

        // Calculate Step Size
        delta = distArray[1] - distArray[0];

        k++;
	}

    if (k >= MAX_ITERATIONS)
    {
        err |= ITERATION_MAX_REACHED_ERR;
    }

    if (fabs(error) >= MAX_DISTANCE_ERROR)
    {
        err |= ERROR_MAX_REACHED_ERR;
    }

    // Create tangentArc
    err |= invCrs(tanArcCenterPt, arc.centerPoint, &crs, NULL, eps);
    err |= direct(tanArcCenterPt, crs, tangentArcRadius, &tanArcStartPt, eps);
    tanArcEndPt = tanArcProjOnGeo;
	err |= createArc(tangentArc, tanArcCenterPt, tanArcStartPt, tanArcEndPt, tanArcDir, tol, eps);

	*tangentArcLocation = tanArcLoc;

    return err;
}

//Finds the geodesic that is tangent to the spiral and intersects the input line at the input angle.
//This function is called by the higher level GeoAngleFromGeoTanToSpiral.The input spiral is created
//by GeoAngleFromGeoTanToSpiral such that no more than 1 point exists.
static ErrorSet initGeoTanToSpiralAtAngleToGeo(Spiral sp, Geodesic geo, double angle, Geodesic* tanGeo, double tol, double eps) {

double diff, spCrs, ptCrs, az12, az21, dist, crs31, crs32, dist31, dist32, perpCrs, perpDist, testAz, step, rad;
LLPoint perpPt, testPt, intPt;
Geodesic tempGeo;
ErrorSet err = 0;
int index = 0;

err  |= projectToGeo(geo.startPoint, geo.startAz, sp.centerPoint, &perpPt, &perpCrs, &perpDist, tol, eps);
ptCrs = geoCrs(geo, perpPt, &az12, &az21, &dist, &err, tol, eps);

testAz = fmod(perpCrs + sp.dir * M_PI / 2 - sp.dir * (M_PI / 2 + angle) + M_2PI, M_2PI);

step = 1;

while (1) {

	err |= ptOnSpiral(sp, testAz, &testPt, eps);
	err |= spiralTanCrs(sp, testAz, &spCrs);

	if (fabs(azDifference(spCrs, geo.startAz)) > M_PI / 2) {
		ptCrs = fmod(geo.startAz + M_PI, M_2PI);
	} else {
		ptCrs = geo.startAz;
	}

	err |= crsIntx(testPt, spCrs, &crs31, &dist31, geo.startPoint, ptCrs, &crs32, &dist32, &intPt, tol, eps);

	diff = azDifference(crs32, crs31);

	if ((fabs(azDifference(crs32, crs31) - angle) * dist31 < tol / 10) && (step < tol)) {
		err |= createGeo(&tempGeo, testPt, intPt, SEGMENT, eps);
		*tanGeo = tempGeo;
		break;
	}
	testAz = testAz - (diff - angle);
	err |= spiralRadius(sp, testAz, &rad);
	step = fabs(diff - angle) * rad;

	index++;

	if (index > 100) {
		err |= UNEXPECTED_ERR;
		break;
	}
}
return err;
}

//Finds the geodesic that is tangent to the spiral and intersects the input line at the input angle.
//This function calls initGeoTanToSpiralAtAngleToGeo to determine all candidate geodesics and then checks to determine
//which of them lie on the input spiral.
ErrorSet geoTanToSpiralAtAngleToGeo(Spiral sp, Geodesic geo, double angle, Geodesic* tanGeo, double tol, double eps) {

double perpDist, rad, dist, fCrs, bCrs, ptCrs, perpCrs, testAz;
int side;
ErrorSet err = 0;
LLPoint perpPt;
Spiral tempSp;
Geodesic testGeo;

err |= projectToGeo(geo.startPoint, geo.startAz, sp.centerPoint, &perpPt, &perpCrs, &perpDist, tol, eps);
ptCrs = geoCrs(geo, perpPt, &fCrs, &bCrs, &dist, &err, tol, eps);
err |= invCrs(perpPt, sp.centerPoint, &fCrs, &bCrs, eps);
side = 1;
if ((sp.dir * azDifference(fCrs, ptCrs)) < 0) {
	testAz = fmod(perpCrs - sp.dir * angle + 2*M_PI, M_2PI);
} else {
	testAz = fmod(perpCrs + M_PI + sp.dir * angle + 2*M_PI, M_2PI);
}

err |= spiralRadius(sp, testAz, &rad);

err |= createSpiralSection(sp, testAz, rad, &tempSp, eps);

err |= initGeoTanToSpiralAtAngleToGeo(sp, geo, angle, &testGeo, tol, eps);

if (ptIsOnSpiral(sp, testGeo.startPoint, tol, eps)) {
	*tanGeo = testGeo;
}
return err;
}

//Finds the points on the input spiral such that the tangent course (or reverse tangent course) at that point
//on the spiral passes through the input point.  This function is called by the higher level PtToSpiralTans.
//The input spiral is created by PtToSpiralTans such that no more than 2 points exist.
static ErrorSet initPtsOnSpiralOnTanThruPt(Spiral sp, LLPoint pt, LLPointPair tanPair, double tol, double eps) {

double ptAz, az21, ptDist, testAz, spCrs, rad, azDiff, ptCrs;
LLPoint spPt;
int sec, index = 0;
ErrorSet err = 0;

sec = 0;

err |= inverse(sp.centerPoint, pt, &ptAz, &az21, &ptDist, eps);
err |= spiralRadius(sp, ptAz, &rad);

if (rad < ptDist) {
	testAz = ptAz - sp.dir * acos(rad / ptDist);
} else {
	return err;
}

while (1) {

	err |= spiralTanCrs(sp, testAz, &spCrs);
	err |= ptOnSpiral(sp, testAz, &spPt, eps);

	err |= inverse(spPt, pt, &ptCrs, &az21, &ptDist, eps);

	if (sec) {
		azDiff = azDifference(ptCrs, fmod(spCrs + M_PI, 2*M_PI));
	} else {
		azDiff = azDifference(ptCrs, spCrs);
	}

	err |= spiralRadius(sp, testAz, &rad);

	if ((fabs(azDiff) * ptDist < tol) && (fabs(azDiff) * rad < tol)) {
		index = 0;
		if (sec) {
			tanPair[1] = spPt;
			return err;
		} else {
			err |= inverse(sp.centerPoint, pt, &ptAz, &az21, &ptDist, eps);
			err |= spiralRadius(sp, ptAz, &rad);
			testAz = ptAz + sp.dir * acos(rad / ptDist);
			tanPair[0] = spPt;
			sec = 1;
		}
	}

	if (index == 100) {
		err |= NO_INTERSECTION_ERR;
		return err;
	}
	index++;

	testAz = fmod(testAz - azDiff + 2*M_PI, 2*M_PI);
}

}

//Finds the points on the input spiral such that the tangent course (or reverse tangent course) at that point
//on the spiral passes through the input point.  This function calls findTanPts to determine
//all candidate points, and then checks to determine which of them lie on the input spiral.
ErrorSet ptsOnSpiralOnTanThruPt(Spiral sp, LLPoint pt, LLPointSet* pts, double tol, double eps) {
double initAz, az21, tempDist, rad;
Spiral tempSp;
LLPointPair temp1, temp2;
int i;
ErrorSet err = 0, temperr1 = 0, temperr2 = 0;

err |=  inverse(sp.centerPoint, pt, &initAz, &az21, &tempDist, eps);
err |= spiralRadius(sp, initAz, &rad);

err |= createSpiralSection(sp, initAz, rad, &tempSp, eps);

temperr1 |= initPtsOnSpiralOnTanThruPt(tempSp, pt, temp1, tol, eps);

if (rad + sp.growthRate * 2 * M_PI < tempDist) {
	err |= createSpiralSection(sp, initAz, rad + sp.growthRate * 2 * M_PI, &tempSp, eps);
} else {
	err |= createSpiralSection(sp, initAz, rad - sp.growthRate * 2 * M_PI, &tempSp, eps);
}

temperr2 |= initPtsOnSpiralOnTanThruPt(tempSp, pt, temp2, tol, eps);

if ((temperr1) && (temperr2)) {
	err |= NO_INTERSECTION_ERR;
}

for (i=0; i<=1; i++) {
	if (ptIsOnSpiral(sp, temp1[i], tol, eps)) {
		if (ptIsInSet(temp1[i], *pts, tol) == 0) {
			err |= addPtToPtSet(pts, &temp1[i]);
		}
	}
	if (ptIsOnSpiral(sp, temp2[i], tol, eps)) {
		if (ptIsInSet(temp2[i], *pts, tol) == 0) {
			err |= addPtToPtSet(pts, &temp2[i]);
		}
	}
}
return err;
}

//Finds the geodesic that is tangent to both input spirals.  The convention is that the line will have the same
//orientation as the first input spiral.  This function is called by the higher level GeoTanToTwoSpirals.
//The input spiral is created by GeoTanToTwoSpirals such that no more than 1 geodesic exists.
static ErrorSet initGeoTanToTwoSpirals(Spiral sp1, Spiral sp2, Geodesic* geo, double tol, double eps) {
LLPoint pt1, pt2, spPt1, spPt2, tempEnd, perpPt, geoStart, geoEnd;
double az12, az21, dist, testAz, spCrs, testCrs, azDiff;
Geodesic testGeo;
ErrorSet err = 0;

int index = 0;

err |= inverse(sp1.centerPoint, sp2.centerPoint, &az12, &az21,
	        	     &dist, eps);

testAz = fmod(az12 - sp1.dir * M_PI /2, 2*M_PI);
err |= direct(sp1.centerPoint, testAz, 5, &pt1, eps);

testAz = fmod(az21 + sp2.dir * M_PI /2, 2*M_PI);
err |= direct(sp2.centerPoint, testAz, 5, &pt2, eps);

err |= createGeo(&testGeo, pt1, pt2, INFINITE, eps);
err |= spiralMidChord(sp1, testGeo, &spPt1, tol, eps);
err |= invCrs(sp1.centerPoint, spPt1, &testAz, &az21, eps);

azDiff = 1;

if (err > 0) {
	err = 0;
}

while (1) {

	err |= ptOnSpiral(sp1, testAz, &spPt1, eps);
	err |= spiralTanCrs(sp1, testAz, &spCrs);

	err |= direct(spPt1, spCrs, 1, &tempEnd, eps);
	err |= createGeo(&testGeo, spPt1, tempEnd, INFINITE, eps);
	err |= spiralMidChord(sp2, testGeo, &spPt2, tol, eps);

	if (err == ITERATION_MAX_REACHED_ERR) {
		err = 0;
	}

	err |= inverse(spPt1, spPt2, &testCrs, &az21, &dist, eps);
	azDiff = azDifference(spCrs, testCrs);
	err |= projectToGeo(spPt1, spCrs, spPt2, &perpPt, &az12, &dist, tol, eps);

	if (dist < tol) {
		geoStart = spPt1;
		break;
	}

	if (index == 100) {
		err |= ITERATION_MAX_REACHED_ERR;
		return err;
	}
	testAz = testAz + azDiff;
	index++;
}

err |= createGeo(&testGeo, pt1, pt2, INFINITE, eps);
err |= spiralMidChord(sp2, testGeo, &spPt2, tol, eps);
err |= invCrs(sp2.centerPoint, spPt2, &testAz, &az21, eps);

index = 0;
azDiff = 1;

if (err > 0) {
	err = 0;
}

while (1) {

	err |= ptOnSpiral(sp2, testAz, &spPt2, eps);
	err |= spiralTanCrs(sp2, testAz, &spCrs);
	spCrs = fmod(spCrs + M_PI, 2*M_PI);

	err |= direct(spPt2, spCrs, 1, &tempEnd, eps);
	err |= createGeo(&testGeo, spPt2, tempEnd, INFINITE, eps);
	err |= spiralMidChord(sp1, testGeo, &spPt1, tol, eps);

	if (err == ITERATION_MAX_REACHED_ERR) {
		err = 0;
	}

	err |= inverse(spPt2, spPt1, &testCrs, &az21, &dist, eps);
	azDiff = azDifference(spCrs, testCrs);
	err |= projectToGeo(spPt2, spCrs, spPt1, &perpPt, &az12, &dist, tol, eps);

	if (dist < tol) {
		geoEnd = spPt2;
		break;
	}

	if (index == 100) {
		err |= ITERATION_MAX_REACHED_ERR;
		return err;
	}
	testAz = testAz + azDiff;
	index++;
}

err |= createGeo(geo, geoStart, geoEnd, SEGMENT, eps);

return err;
}

//Finds the geodesic that is tangent to both input spirals.  The convention is that the line will have the same
//orientation as the first input spiral.  This function calls initGeoTanToTwoSpirals to determine
//all candidate lines, and then checks to determine which of them lie on the both input spirals.
ErrorSet geoTanToTwoSpirals(Spiral sp1, Spiral sp2, LLPointSet* pts, double tol, double eps) {

Geodesic geo;
Spiral tempSp1, tempSp2;
double crs1, crs2, dist, initAz1, initAz2, rad1, rad2;

ErrorSet err = 0;

err |= inverse(sp1.centerPoint, sp2.centerPoint, &crs1, &crs2, &dist, eps);

if (sp1.dir == -1) {
		initAz1 = fmod(crs1 + M_PI / 2, M_2PI);
		initAz2 = fmod(crs2 - M_PI / 2, M_2PI);
} else {
	initAz1 = fmod(crs1 - M_PI / 2, M_2PI);
	initAz2 = fmod(crs2 + M_PI / 2, M_2PI);
}

err |= spiralRadius(sp1, initAz1, &rad1);
err |= createSpiralSection(sp1, initAz1, rad1, &tempSp1, eps);

err |= spiralRadius(sp2, initAz2, &rad2);
err |= createSpiralSection(sp2, initAz2, rad2, &tempSp2, eps);

err |= initGeoTanToTwoSpirals(tempSp1, tempSp2, &geo, tol, eps);


if (ptIsOnSpiral(sp1, geo.startPoint, tol, eps) && ptIsOnSpiral(sp2, geo.endPoint, tol, eps)) {
	err |= addPtToPtSet(pts, &geo.startPoint);
	err |= addPtToPtSet(pts, &geo.endPoint);
} else {
	err |= POINT_NOT_ON_ARC_ERR;
}

return err;

}
} //namespace
