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

/*  \file Locus.c
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
 * Allocate and initialize a new Locus structure from given input parameters
 */

ErrorSet createLocus(Locus* loc, LLPoint geoStart, LLPoint geoEnd,
                   double startDist, double endDist, LineType lineType,
                   double tol, double eps)
{

    ErrorSet err = 0;

    double crs12, crs21, length;
    double crsToLoc;

    LLPoint tempPt;

    /* Assign local storage if optional pointers are not provided */
    Locus npLoc;
    if (NULL == loc) loc = &npLoc;

    loc->geoStart = geoStart;
    loc->geoEnd = geoEnd;
    if (ptIsAtPole(geoStart, &err, tol, eps) != 0)
    {
        //geoStart is at a pole. Make sure geoStart has same longitude as geoEnd.
        loc->geoStart.longitude = geoEnd.longitude;
    }
    else if (ptIsAtPole(geoEnd, &err, tol, eps) != 0)
    {
        //geoEnd is at a pole. Make sure geoStart has same longitude as geoStart.
        loc->geoEnd.longitude = geoStart.longitude;
    }
    loc->startDist = startDist;
    loc->endDist = endDist;
    loc->lineType = lineType;

    err |= inverse(geoStart, geoEnd, &crs12, &crs21, &length, eps);
    if (err)
        return err;
    if (length < tol)
    {
        err |= INVALID_SHAPE_ERR;
        return err;
    }

    loc->geoLength = length;
    loc->geoAz = crs12;
    loc->geoRevAz = crs21;

    if (startDist > 0.0)
    {
        crsToLoc = crs12 + M_PI_2;
    }
    else
    {
        crsToLoc = crs12 - M_PI_2;
    }
    err |= direct(geoStart, crsToLoc, fabs(startDist), &tempPt, eps);
    if (err)
        return err;

    loc->locusStart = tempPt;

    if (endDist > 0.0)
    {
        crsToLoc = crs21 - M_PI_2;
    }
    else
    {
        crsToLoc = crs21 + M_PI_2;
    }
    err = direct(geoEnd, crsToLoc, fabs(endDist), &tempPt, eps);
    loc->locusEnd = tempPt;

    loc->slope = (endDist - startDist) / (loc->geoLength);

    return err;

}

static ErrorSet ptOnGeoAtDist(Geodesic geo, double distToStart, LLPoint* geoPt, double* geoAz, double tol, double eps)
 {
     ErrorSet err = SUCCESS;
     double npGeoAz;
     LLPoint npGeoPt;
     LLPoint newGeoStartPt = geo.startPoint;

     /* Assign local storage if optional pointers are not provided */
     if (NULL == geoPt) geoPt = &npGeoPt;
     if (NULL == geoAz) geoAz = &npGeoAz;

     /* Default values for case fabs(distToStart) <= tol */
     *geoAz = geo.startAz;
     *geoPt = geo.startPoint;

     if ( (distToStart > tol) || (distToStart < -tol) )
     {

         if ( (distToStart <= 10.0*SMALL_DIST_THRESHOLD) && (distToStart >= -10.0*SMALL_DIST_THRESHOLD) )
         {
             /* This block was added to fix Bug 24162 */
             /* Too close to geodesic start point for precise course, so back geoStart up 1000
              * meters and recalc course to start.  Also designed to improve precision for very short
              * geodesics */
             double fwdAccuracyDist = 1000.0/1852.0;  /* convert 1000 m to NM */
             double geoRevAz = modpos(geo.endAz + M_PI, M_2PI);

             /* Start at geodesic end point and project new start point */
             err |= direct(geo.endPoint, geoRevAz, geo.length + fwdAccuracyDist, &newGeoStartPt, eps);
             err |= inverse(newGeoStartPt, geo.endPoint, geoAz, NULL, NULL, eps);
             /* Adjust distToStart to account more movement of start point */
             distToStart = distToStart + fwdAccuracyDist;  /* Will be > 0 */
         }

         /* direct will handle distToStart<0 case */
         err |= direct(newGeoStartPt,*geoAz,distToStart,geoPt,eps);
         err |= inverse(*geoPt,newGeoStartPt,geoAz,NULL,NULL,eps);
         /* fcrs is now the azimuth from geoPt to newGeoStart. Two cases need to be handled:
          * If distToStart < 0, then fcrs is in same general direction as original geodesic's startCrs and fcrs+M_PI_2 will be to right of geodesic
          * If distToStart > 0, then fcrs is opposite from original geodesic's startCrs fcrs-M_PI_2 will be to right of geodesic
          */
         if (distToStart > 0)
         {
             *geoAz = modcrs(*geoAz + M_PI);
         }

     }
     return err;

 }

static ErrorSet ptAbeamGeoAtDist(Geodesic geo, double distToStart, double abmDist, LLPoint* abmPt,
                                              double* crsFromPoint, double tol, double eps)
 {

     ErrorSet err = SUCCESS;
     double npCrsFromPoint = 0.0;
     double geoAz;
     LLPoint geoPt, npAbmPt;

     /* Assign local storage if optional pointers are not provided */
     if (NULL == abmPt) abmPt = &npAbmPt;
     if (NULL == crsFromPoint) crsFromPoint = &npCrsFromPoint;

     err |= ptOnGeoAtDist(geo,distToStart,&geoPt,&geoAz,tol,eps);
     if (err) return err;

     if (abmDist > 0)
     {
         geoAz = modcrs(geoAz + M_PI_2);
     }
     else
     {
         geoAz = modcrs(geoAz - M_PI_2);
         abmDist = -abmDist;
     }

     if ( (abmDist <= tol) && (abmDist >= -tol) )
     {
         *abmPt = geoPt;
         *crsFromPoint = modcrs(geoAz + M_PI);
     }
     else if (abmDist <= 40.0*SMALL_DIST_THRESHOLD)
     {
         /* Point sought is very close to geodesic.  To improve precision, project
          * farther and then back to desired point.
          */
         LLPoint tmpPoint;
         double fwdAccuracyDist = 1000.0/1852.0; /* Convert 1000.0 m to NM */

         err |= direct(geoPt,geoAz,fwdAccuracyDist+abmDist,&tmpPoint,eps);
         err |= inverse(tmpPoint,geoPt,crsFromPoint,NULL,NULL,eps);
         err |= direct(tmpPoint,*crsFromPoint,fwdAccuracyDist,abmPt,eps);
         err |= inverse(*abmPt,tmpPoint,crsFromPoint,NULL,NULL,eps);
         *crsFromPoint = modcrs(*crsFromPoint+M_PI);
     }
     else
     {
         err |= direct(geoPt,geoAz,abmDist,abmPt,eps);
         err |= inverse(*abmPt,geoPt,crsFromPoint,NULL,NULL,eps);
     }

     return err;

 }

 /******************************************************************************
  * Find the course of locus at given point ON THE GEODESIC. Valid return values
  * are in range [0, 2*PI).
  *
  * @return Returns course (in radians) of locus at locus point abeam given
  * point on geodesic.
  *
  * @param loc Locus structure defining locus parameters
  * @param ptOnGeo LLPoint on locus's geodesic from which locus course will be computed
  * @param tol Required error tolerance (passed to subsequent function calls)
  * @param eps Fundamental forward/inverve computation accuracy
  *
  * @retval locPt LLPoint on locus at given ptOnGeo
  * @retval perpCrs Course from locPt to given ptOnGeo
  * @retval err Error code indicating exit status
  *
  */
double locusCrsAtGeoDist(Locus loc, double distFromGeoStart, LLPoint* locPt,
                             double* perpCrs, ErrorSet* err, double tol, double eps)
 {



     Geodesic geo;

     double locCrs = -1.0;
     double distToLoc = distToLocusFromGeoDist(loc,distFromGeoStart);
     LLPoint npLocPt;
     double npPerpCrs;

     /* Assign local storage if optional pointers are not provided */
     if (NULL == locPt) locPt = &npLocPt;
     if (NULL == perpCrs) perpCrs = &npPerpCrs;

 //    *err |= createGeo(&geo,loc.geoStart,loc.geoEnd,loc.lineType,eps);
     // Avoid recalculating geodesic parameters that are already stored in locus struct
     geo.startPoint = loc.geoStart;
     geo.endPoint = loc.geoEnd;
     geo.startAz = loc.geoAz;
     geo.length = loc.geoLength;
     geo.endAz = modpos(loc.geoRevAz+M_PI,M_2PI);
     geo.lineType = loc.lineType;

     *err |= ptAbeamGeoAtDist(geo,distFromGeoStart,distToLoc,locPt,&locCrs,tol,eps);
     /* locCrs is az from locus point to corresponding point on geodesic */
     if (distToLoc > 0)
     {
         /* Loc is to right of geodesic at this point, so rotate locCrs clockwise */
         locCrs = locCrs + M_PI_2;
     }
     else
     {
         locCrs = locCrs - M_PI_2;
     }

     locCrs = locCrs + atan(loc.slope);
     locCrs = modpos(locCrs, M_2PI);

     return locCrs;

 }

/*******************************************************************************
 * Find the distance to a locus given a distance along the geodesic
 * double distToLocusFromGeoDist (Locus loc, double distance, double eps);
 */

double distToLocusFromGeoDist(Locus loc, double distance)
{
    double distToLoc;

    if (loc.geoLength == 0.0)
    {
        printf("geodesic length is zero\n");
        distToLoc = 0.0;
    }
    /* Mike thinks we want to allow this case */
    //    else if (distance > geolen)
    //    {
    //        printf("distance is > length of geodesic\n");
    //        distToLoc = 0.0;
    //    }

    else
    {
        distToLoc = loc.startDist + distance * (loc.endDist - loc.startDist)
                / loc.geoLength;
    }
    return distToLoc;
}

/*******************************************************************************
 * Find the distance to a locus given a point on the geodesic double
 * distToLocusFromGeoPt (Locus loc, LLPoint geopt, double *faz, double tol, double
 * eps);
 * Pointer faz is used to return forward azimuth of geodesic at geopt.  This is
 * needed if geopt is not between geoStart and geoEnd.
 */
double distToLocusFromGeoPt(Locus loc, LLPoint geoPt, double* faz, ErrorSet* err,
                         double tol, double eps)
{
    LLPoint newPoint;
    double distToLoc;
    double crsFromPoint, distFromPoint, angle1;
    double crsFromStart, crsToStart, distToStart;
    int online;
    ErrorSet newErr = 0;
    double maxDist = SPHERE_RADIUS_NMI*M_PI/2.0;  /* Quarter circumference */
    double npFaz = 0.0;

    /* Assign local storage if optional pointers are not provided */
    if (NULL == faz) faz = &npFaz;

    online = ptIsOnGeo(loc.geoStart, loc.geoEnd, geoPt, loc.lineType,
            &newErr, tol, eps);

    if (newErr)
    {
        *err |= newErr;
        return 0.0;
    }

    if (!online)
    {
        newErr |= projectToGeo(loc.geoStart, loc.geoAz, geoPt, &newPoint,
                &crsFromPoint, &distFromPoint, tol, eps);
        if (newErr)
        {
            *err |= newErr;
            return 0.0;
        }
    }
    else
    {
        newPoint = geoPt;
    }

    /* Find azimuth, rev. azimuth, & distance from test point to start of geodesic */
    newErr |= inverse(newPoint, loc.geoStart, &crsToStart, &crsFromStart,
            &distToStart, eps);

    if (newErr)
    {
        *err |= newErr;
        return 0.0;
    }

    /* If geoPt is behind geoStart, then distance passed to WGSDistToLocusD
     * must be negative and direction from geoPt to ptonloc must be reversed */
    angle1 = fabs(modlon(loc.geoAz - crsFromStart));
    if (distToStart > tol)
    {
        if (angle1 > M_PI_2)
        { /* Case 1: geoPt lies behind geoStart */
            distToStart = -distToStart;
            *faz = crsToStart;
        }
        else
        {
            *faz = crsToStart + M_PI;
        }
        /* Find distance from geoPt to locus */
        distToLoc = distToLocusFromGeoDist(loc, distToStart);
    }
    else
    {
        distToLoc = loc.startDist;
        *faz = loc.geoAz;
    }

    return distToLoc;

}

/*******************************************************************************
 * For a given point on the geodesic of a locus, find the corresponding
 * point on the locus.
 */
ErrorSet ptOnLocusFromGeoPt(Locus loc, LLPoint geoPt, LLPoint* ptOnLoc,
                     double* perpCrs, double tol, double eps)
{

    double distp, crsToStart, crsFromStart, distToStart;
    double fcrs, angle1;
    ErrorSet err = 0;
    double maxDist = SPHERE_RADIUS_NMI*M_PI/2.0; /* Quarter circumference */
    double npPerpCrs;
    LLPoint npPtOnLoc;

    /* Assign local storage if optional pointers are not provided */
    if (NULL == perpCrs) perpCrs = &npPerpCrs;
    if (NULL == ptOnLoc) ptOnLoc = &npPtOnLoc;

    /* Find azimuth, rev. azimuth, & distance from test point to start of geodesic */
    err |= inverse(geoPt, loc.geoStart, &crsToStart, &crsFromStart,
            &distToStart, eps);

//    err |= minSubtendedAngle(crsToStart, loc.geoAz, &angle1);
//
//    if (angle1 < M_PI / 2) {
//    	err |= direct(loc.geoStart, loc.geoAz + M_PI, distToStart, &geoPt, eps);
//    } else {
//    	err |= direct(loc.geoStart, loc.geoAz, distToStart, &geoPt, eps);
//    }
//
//    err |= invCrs(geoPt, loc.geoStart, &crsToStart, &crsFromStart, eps);

    if (err) {
        return err;
    }

    /* If geoPt is behind geoStart, then distance passed to WGSDistToLocusD
     * must be negative and direction from geoPt to ptonloc must be reversed */
    if (distToStart <= tol)
    {
        /* if geoPt coincides with geoStart, then no calculation necessary */
        *ptOnLoc = loc.locusStart;
        if (loc.startDist >= 0.0)
            *perpCrs = modcrs(loc.geoAz + M_PI_2);
        else
            *perpCrs = modcrs(loc.geoAz - M_PI_2);
    }
    else
    {
        if (distToStart <= 100*SMALL_DIST_THRESHOLD)
        {
            /* This block was added to fix Bug 24162 */
            /* Too close to geodesic start point for precise course, so back geoStart up 100
             * meters and recalc course to start */
            LLPoint newGeoStart = loc.geoStart;
            err |= direct(loc.geoEnd, loc.geoRevAz, loc.geoLength + 100.0
                    / 1852.0, &newGeoStart, eps);
            err |= invCrs(geoPt, newGeoStart, &crsToStart, NULL,eps);
            /* Even if original point was behind geoStart, it could not have been farther
             * than SMALL_DIST_THRESHOLD, which is much less than the distance we moved geoStart.
             * Therefore, in any case, fcrs is now opposite crsToStart (fcrs is toward geoEnd).
             */
            fcrs = modcrs(crsToStart + M_PI);
            angle1 = fabs(modlon(loc.geoAz - crsFromStart));
            if (angle1 > M_PI_2)
            { /* Case 1: geoPt lies behind geoStart */
                distToStart = -distToStart;
            }

            if (err)
                return err; /* Can't continue if error has occurred here */
        }
        else
        {
            angle1 = fabs(modlon(loc.geoAz - crsFromStart));
            if (angle1 > M_PI_2)
            { /* Case 1: geoPt lies behind geoStart */
                distToStart = -distToStart;
                fcrs = crsToStart;
            }
            else
            {
                fcrs = modcrs(crsToStart + M_PI);
            }
        }

        /* Find distance from geoPt to locus */
        distp = distToLocusFromGeoDist(loc, distToStart);

        if (fabs(distp) <= tol)
        {
            /* locus start is same as geodesic start */
            *ptOnLoc = geoPt;
            return err;
        }

        if (distp > 0.0) // locus lies to right of geodesic
            *perpCrs = modcrs(fcrs + M_PI_2);
        else
            // locus lies to left of geodesic
            *perpCrs = modcrs(fcrs - M_PI_2);

        if (err |= direct(geoPt, *perpCrs, fabs(distp), ptOnLoc, eps))
            return err;
    }
    return err;

}

/*******************************************************************************
 * Determine whether a given point is on a locus.  If point is on Locus,
 * then pointer to projected point on geodesic is returned; otherwise, NULL is
 * returned.  Returning pointer to projected point allows use of that point in
 * other calculations, possibly saving steps.
 *
 *
 */
int ptIsOnLocus(Locus loc, LLPoint testPt, LLPoint* ptOnGeo, ErrorSet* err,
                     double tol, double eps)
{

    ErrorSet newErr = 0;

    double fcrs = loc.geoAz;
    double crsFromPoint, distFromPoint;
    double testCrs, testRevCrs, locDist;
    double maxDist = 0.0, minDist = 0.0, distError = 0.0;
    int online = 0;
    LLPoint npPtOnGeo;

    /* Assign local storage if optional pointers are not provided */
    if (ptOnGeo == NULL) ptOnGeo = &npPtOnGeo;

    if (ptsAreSame(testPt, loc.locusStart, tol))
    {
        *ptOnGeo = loc.geoStart;
        return 1;
    }
    else if (ptsAreSame(testPt, loc.locusEnd, tol))
    {
        *ptOnGeo = loc.geoEnd;
        return 1;
    }

    if (ptIsOnGeo(loc.geoStart, loc.geoEnd, testPt, loc.lineType, err,
            tol, eps) == 1)
    {
        //testPt is on the geodesic already
        //ptOnGeo = &testPt;
        *ptOnGeo = testPt;
        distFromPoint = 0.0;
        crsFromPoint = 0.0;
    }
    else
    {
        /* Project test point onto locus's geodesic.  This also returns the shortest
         * distance from the test point to the geodesic */

        newErr |= projectToGeo(loc.geoStart, fcrs, testPt, ptOnGeo,
                &crsFromPoint, &distFromPoint, tol, eps);

        if (newErr)
        {
            *err |= newErr;
            return 0;
        }

        /* Test whether projected point lies between end points that
         * define the geodesic. This will account for an infinite or
         * semi-infinite locus. */
        online = ptIsOnGeo(loc.geoStart, loc.geoEnd, *ptOnGeo,
                loc.lineType, &newErr, tol, eps);

        if ((!online) || (newErr))
        {
            /* If projected point is not on geodesic, then test point
             * cannot be on locus.*/
            *err |= newErr;
            return 0;
        }
    }

    /* Determine which side of locus the test point lies on */
    if (ptIsAtPole(loc.geoStart, err, tol, eps) == 0)
    {
        /* The point geoStart is NOT at a pole. Calculate course from geoStart to testPt */
        *err |= invCrs(loc.geoStart, testPt, &testCrs, &testRevCrs, eps);
    }
    else
    {
        /* The point geoStart is at a pole. We already have the course we need. */
        testCrs = -crsFromPoint;
    }
    testCrs = modlon(fcrs - testCrs);
    if (testCrs > 0)
    {
        /* Apply appropriate sign to distance */
        distFromPoint = -distFromPoint;
    }

    /* Calculate correct expected locus distance */
    locDist = distToLocusFromGeoPt(loc, *ptOnGeo, NULL, &newErr, tol, eps);
    if (newErr)
    {
        *err |= newErr;
        return 0;
    }

    /* ptOnGeo may be mislocated by up to tol.
     * This error in position translates to an error in locus distance, which may
     * be magnified if the locus slope is greater than 45 degrees */
    distError = fabs(tol * (loc.endDist - loc.startDist) / loc.geoLength);
    if (fabs(distError) < tol)
    {
        distError = tol;
    }
    minDist = locDist - distError;
    maxDist = locDist + distError;

    /* Return true if actual distance is within tol of expected distance */
    return ((distFromPoint > minDist) && (distFromPoint < maxDist));
    //    return (fabs(distFromPoint - locDist) <= tol);

}

/******************************************************************************
 * Find the course of locus at given point ON LOCUS. Valid return values are in
 * range [0, 2*PI).  Invalid return value is -1.0. If -1.0 is returned, then
 * err is also set to a value != 0.  The purpose of the -1.0 return value is to
 * help ensure that the caller of this function does not mistake the returned
 * course as being valid and correct.
 *
 * This function must project a point on the geodesic from the given locus point.
 * This requires a call to projectToGeo, which is subject to errors on the
 * order of tol.  These errors may propagate up through calling functions, leading
 * to compounding of errors.  Use locusCrsAtPt when possible.
 */

double locusCrsAtPt(Locus loc, LLPoint pt, LLPoint* geoPt,
                            double* perpCrs, ErrorSet* err, double tol, double eps)
{

    /* Calculated value perpCrs is course perpendicular to locus, toward geodesic. */
    Geodesic geo;
    double locCrs;
    //    double geoLen;
    double distToLoc;
    double slope;
    double temp;
    double faz;
    ErrorSet newErr = 0;
    int online = 0;
    LLPoint npGeoPt;
    double npPerpCrs;

    /* Assign local storage if optional pointers are not provided */
    if (NULL == geoPt) geoPt = &npGeoPt;
    if (NULL == perpCrs) perpCrs = &npPerpCrs;

    loc.lineType = INFINITE; //Bug Fix 15132 - SB - for this method, locus should be INFINITE always.
    online = ptIsOnLocus(loc, pt, geoPt, &newErr, tol, eps);
    if (!online)
    {
        newErr |= POINT_NOT_ON_LINE_ERR; //actually, not on the locus but we don't have POINT_NOT_ON_LOCUS_ERR and the enum is full
    }
    if (newErr)
    {
        *err |= newErr;
        return -1.0;
    }

    createGeo(&geo, loc.geoStart, loc.geoEnd, loc.lineType, eps);
    if (ptsAreSame(*geoPt, pt, tol))
    {
        /* pt and geoPt are the same (pt is on geodesic),
         * so can't find course from pt to geoPt
         * Calculate the course of the geodesic at geoPt and add
         * atan(slope) to it.
         */
        locCrs = geoCrs(geo, *geoPt, NULL, NULL, NULL,
        &newErr, tol, eps);
        if ((newErr != 0) || (locCrs == -1.0))
        {
            *err = newErr;
            return -1.0;
        }
        /* slope ratio */
        slope = (loc.endDist - loc.startDist) / geo.length;
        /* slope < 0 implies the locus is tending leftward */
        /* slope > 0 implies the locus is tending rightward */
        /* slope angle (in radians) */
        slope = atan(slope);
        //TODO: this choice for perpCrs always puts it to the right of the geodesic.  Is this correct?
        *perpCrs = modcrs(locCrs - M_PI_2);
        locCrs = locCrs + slope;
    }
    else
    {
        /* Calculate perpCrs: course from point on locus to corresponding point on geodesic */
        newErr |= invCrs(pt, *geoPt, perpCrs, &temp, eps);
        if (newErr)
        {
            *err |= newErr;
            return -1.0;
        }

        distToLoc = distToLocusFromGeoPt(loc, *geoPt, &faz, &newErr, tol, eps);
        if (newErr)
        {
            *err |= newErr;
            return -1.0;
        }

        //TODO:  Why not just use slope attribute of Locus structure? (MJM asks)
        /* slope ratio */
        slope = (loc.endDist - loc.startDist) / geo.length;
        /* slope < 0 implies the locus is tending leftward */
        /* slope > 0 implies the locus is tending rightward */

        /* slope angle (in radians) */
        slope = atan(slope);

        *perpCrs = *perpCrs + slope;

        /* course in direction of locus at given point */
        if (distToLoc < 0)
        {
            locCrs = *perpCrs - M_PI_2;
        }
        else
        {
            locCrs = *perpCrs + M_PI_2;
        }
    }

    locCrs = modpos(locCrs, M_2PI);

    return locCrs;

}



/* NOTE:  This function may return a very large slope in some cases */
ErrorSet locusFromGeoAndPts(Locus* loc, LLPoint geoStart, double geoStartAz,
                                            LLPoint locStart, LLPoint locEnd, double tol, double eps)
{

    ErrorSet err = SUCCESS;
    LLPoint startPtOnGeo, endPtOnGeo;
    double crsFromPt = 0.0, startDist = 0.0;
    double endDist = 0.0;
    double geoLength = 0.0;
    double geoRevAz = 0.0;
    double startAngle = 0.0, endAngle = 0.0;
    double startLateralAz = 0.0, endLateralAz = 0.0;
    Locus npLoc;

    /* Assign local storage if optional pointers are not provided */
    if (NULL == loc) loc = &npLoc;

    err |= projectToGeo(geoStart, geoStartAz, locStart,
            &startPtOnGeo,&crsFromPt,&startDist,tol,eps);
    err |= projectToGeo(geoStart, geoStartAz, locEnd,
            &endPtOnGeo,&crsFromPt,&endDist,tol,eps);

    err |= inverse(startPtOnGeo, endPtOnGeo, &geoStartAz, &geoRevAz, &geoLength,eps);
    err |= inverse(startPtOnGeo,locStart,&startLateralAz,NULL,&startDist,eps);
    err |= inverse(endPtOnGeo,locEnd,&endLateralAz,NULL,&endDist,eps);
    if (geoLength <= tol)
    {
        /* In this case, the locus is basically perpendicular to the
         * given geodesic and the slope is ~infinite and we can't really
         * discriminate between the start/end point (i.e., the geodesic is
         * so short we can't assign a direction to it).
         *
         * We'll assume the convention that if locus start/end points are on opposite sides of the geodesic, then
         * the startDist>0 and endDist<0.
         */
        if (fabs(modlon(startLateralAz-endLateralAz)) > M_PI_2) endDist = -endDist;
    } else  {

        /* Check which side of geo start and end points of locus lie.  Set distance
         * to negative values for those points to the left of geo */
        startAngle = modlon(startLateralAz-geoStartAz);
        if (startAngle < 0.0) startDist = -startDist;

        endAngle = modlon(geoRevAz-endLateralAz);
        if (endAngle < 0.0) endDist = -endDist;
    }

    if (geoLength < (endDist - startDist)*(1.0e-99))
        loc->slope = 1.0e99;
    else
        loc->slope = (endDist-startDist)/geoLength;

    loc->geoStart = startPtOnGeo;
    loc->geoEnd = endPtOnGeo;
    loc->geoLength = geoLength;
    loc->locusStart = locStart;
    loc->locusEnd = locEnd;
    loc->geoAz = geoStartAz;
    loc->geoRevAz = geoRevAz;
    loc->startDist = startDist;
    loc->endDist = endDist;
    loc->lineType = SEGMENT;

//    displayMatlabLocus(*loc,"tmpLoc",0);

    return err;
}

int lociCoincide(Locus locus1, Locus locus2, Shape* commonShape, ErrorSet* err, double tol, double eps)
{
	ErrorSet 			newErr = 0;
    LLPoint 		commonPt1;
    LLPoint			commonPt2;
    double			commonStartDist;
    double			commonEndDist;
    Locus 			tempLoc2;
    int 			loc1StartPtMatch = 0;
    int 			loc1StartEqualsLoc2Start = 0;
    int 			loc1StartEqualsLoc2End = 0;
    int 			loc1EndPtMatch = 0;
    int 			loc1EndEqualsLoc2Start = 0;
    int 			loc1EndEqualsLoc2End = 0;
    int 			loc1StartOnLoc2 = 0;
    int 			loc1EndOnLoc2 = 0;
    int 			loc2StartOnLoc1 = 0;
    int 			loc2EndOnLoc1 = 0;
    double 			distLoc1GeoStartToLoc2GeoStart = 0;
    double 			distIntX = 0;
    double                      angle, angle1, fcrs, bcrs;
    LLPoint 		intXGeoPt, intersection;
    int 			intXOnLoc1 = 0, intXOnLoc2 = 0;

    double dummy1;
    LLPoint dummy5;
    double MaxDist = SPHERE_RADIUS_NMI * M_PI_2; /* Quarter circumference */


	//null checks
    if (commonShape == NULL || err == NULL)
    {
    	if (err == NULL){
    		err = new ErrorSet();
    	}
    	*err |= NO_MEMORY_ALLOCATED_ERR;
    	return 0;
    }

    //check if loci are well defined
    if (	ptsAreSame(locus1.geoStart, locus1.geoEnd, tol) ||
    		ptsAreSame(locus2.geoStart, locus2.geoEnd, tol) ||
    		(locus1.startDist == 0 && locus1.endDist == 0) ||
    		(locus2.startDist == 0 && locus2.endDist == 0)
    	)
    {
    	*err |= INVALID_SHAPE_ERR;
    	return 0;
    }

    //check for collinear geodesi
    if ( !(ptIsOnGeo(locus1.geoStart, locus1.geoEnd, locus2.geoStart, INFINITE, &newErr, tol, eps) &&
    		ptIsOnGeo(locus1.geoStart, locus1.geoEnd, locus2.geoEnd, INFINITE, &newErr, tol, eps)) )
    {
    	*err |= newErr;
    	return 0;//loci cannot coincide if geodesi are not collinear
    }

    //orient locus 2 so that its start geo az is in the same direction as that of locus 1
        newErr |= minSubtendedAngle(locus2.geoAz, locus1.geoAz, &angle);
        if (angle < M_PI_4)
        {
          //This block was added to deal with loc1 near a Pole.
          newErr |= invCrs(locus1.geoStart, locus2.geoStart, &fcrs, &bcrs, eps);
          newErr |= minSubtendedAngle(locus1.geoAz, fcrs, &angle1);
          if (angle1 < M_PI_4)
            newErr |= minSubtendedAngle(locus2.geoAz, bcrs + M_PI, &angle);
        }
	if (angle < M_PI_4)// locus 1 and locus 2 point in same direction
	//if (locus2.geoAz - locus1.geoAz < M_PI_4)// locus 1 and locus 2 point in same direction
	{
		newErr |= createLocus(&tempLoc2, locus2.geoStart, locus2.geoEnd, locus2.startDist, locus2.endDist, locus2.lineType, tol, eps);
	}
	else if (angle > M_PI_4) // locus 1 and locus 2 point in opposite directions
	//else if (locus2.geoAz - locus1.geoAz > M_PI_4) // locus 1 and locus 2 point in opposite directions
	{
		newErr |= createLocus(&tempLoc2, locus2.geoEnd, locus2.geoStart, -locus2.endDist, -locus2.startDist, locus2.lineType, tol, eps);
	}
	else //contradiction here since the assumption is that geodesics are collinear - not good if this is happening because it means the geodesic are perpendicular
	{
		newErr |= UNEXPECTED_ERR;
		*err |= newErr;
		return 0;
	}


	//check end points
    loc1StartEqualsLoc2Start = ptsAreSame(locus1.locusStart, tempLoc2.locusStart, tol);
    loc1StartEqualsLoc2End = ptsAreSame(locus1.locusStart, tempLoc2.locusEnd, tol);
    loc1StartPtMatch = loc1StartEqualsLoc2Start || loc1StartEqualsLoc2End;

    loc1EndEqualsLoc2Start = ptsAreSame(locus1.locusEnd, tempLoc2.locusStart, tol);
    loc1EndEqualsLoc2End = ptsAreSame(locus1.locusEnd, tempLoc2.locusEnd, tol);
    loc1EndPtMatch = loc1EndEqualsLoc2Start || loc1EndEqualsLoc2End;

    loc1StartOnLoc2 = ptIsOnLocus(tempLoc2, locus1.locusStart, &dummy5, &newErr, tol, eps);
    loc1EndOnLoc2 = ptIsOnLocus(tempLoc2, locus1.locusEnd, &dummy5, &newErr, tol, eps);
    loc2StartOnLoc1 = ptIsOnLocus(locus1, tempLoc2.locusStart, &dummy5, &newErr, tol, eps);
    loc2EndOnLoc1 = ptIsOnLocus(locus1, tempLoc2.locusEnd, &dummy5, &newErr, tol, eps);


	/* if both loc1 start/end points match loc2 start/end
	 * points then loci must coincide, assuming both loci are well defined
	 * (i.e. neither has start/end points that are identical)
	 */
	if (loc1StartPtMatch && loc1EndPtMatch)
	{
		commonShape->type = LOCUS;
		commonShape->this_shape = malloc(sizeof(Locus));
		newErr |= createLocus((Locus*)commonShape->this_shape, locus1.geoStart, locus1.geoEnd, locus1.startDist, locus1.endDist, SEGMENT, tol, eps);
	}
	/* check the cases where only  one locus start or end point is common
	 * and the other end point lies on a locus
	 */
	else if ( loc1StartEqualsLoc2Start && ( loc1EndOnLoc2 || loc2EndOnLoc1 ) )
	{
		commonShape->type = LOCUS;
		commonShape->this_shape = malloc(sizeof(Locus));

		if (loc1EndOnLoc2)
		{
			newErr |= createPt(&commonPt2, locus1.geoEnd.latitude, locus1.geoEnd.longitude);
			commonEndDist = locus1.endDist;
		}
		else
		{
			newErr |= createPt(&commonPt2, locus2.geoEnd.latitude, locus2.geoEnd.longitude);
			commonEndDist = locus2.endDist;
		}

		newErr |= createLocus((Locus*)commonShape->this_shape, locus1.geoStart, commonPt2, locus1.startDist, commonEndDist, SEGMENT, tol, eps);
	}
	else if ( loc1EndEqualsLoc2End  && ( loc1StartOnLoc2 || loc2StartOnLoc1 ) )
	{
		commonShape->type = LOCUS;
		commonShape->this_shape = malloc(sizeof(Locus));

		if (loc1StartOnLoc2)
		{
			newErr |= createPt(&commonPt1, locus1.geoStart.latitude, locus1.geoStart.longitude);
			commonStartDist = locus1.startDist;
		}
		else
		{
			newErr |= createPt(&commonPt1, locus2.geoStart.latitude, locus2.geoStart.longitude);
			commonStartDist = locus2.startDist;
		}

		newErr |= createLocus((Locus*)commonShape->this_shape, commonPt1, locus1.geoEnd, commonStartDist, locus1.endDist, SEGMENT, tol, eps);
	}
	/* check the cases where no end points are common but some
	 * portion of loci coincide
	 */
	else if ( ( (loc1EndOnLoc2 && loc2StartOnLoc1) || (loc2EndOnLoc1 && loc1StartOnLoc2) ) && !loc1StartPtMatch && !loc1EndPtMatch )
	{
		commonShape->type = LOCUS;
		commonShape->this_shape = malloc(sizeof(Locus));

		if (loc1EndOnLoc2 && loc2StartOnLoc1)
		{
			newErr |= createPt(&commonPt1, locus2.geoStart.latitude, locus2.geoStart.longitude);
			commonStartDist = locus2.startDist;

			newErr |= createPt(&commonPt2, locus1.geoEnd.latitude, locus1.geoEnd.longitude);
			commonEndDist = locus1.endDist;
		}
		else
		{
			newErr |= createPt(&commonPt1, locus1.geoStart.latitude, locus1.geoStart.longitude);
			commonStartDist = locus1.startDist;

			newErr |= createPt(&commonPt2, locus2.geoEnd.latitude, locus2.geoEnd.longitude);
			commonEndDist = locus2.endDist;
		}

		newErr |= createLocus((Locus*)commonShape->this_shape, commonPt1, commonPt2, commonStartDist, commonEndDist, SEGMENT, tol, eps);
	}
	else if ( ( (loc2StartOnLoc1 && loc2EndOnLoc1) || (loc1StartOnLoc2 && loc1EndOnLoc2) ) && !loc1StartPtMatch && !loc1EndPtMatch )
	{
		commonShape->type = LOCUS;
		commonShape->this_shape = malloc(sizeof(Locus));

		if (loc2StartOnLoc1 && loc2EndOnLoc1)
		{
			newErr |= createPt(&commonPt1, locus2.geoStart.latitude, locus2.geoStart.longitude);
			commonStartDist = locus2.startDist;

			newErr |= createPt(&commonPt2, locus2.geoEnd.latitude, locus2.geoEnd.longitude);
			commonEndDist = locus2.endDist;
		}
		else
		{
			newErr |= createPt(&commonPt1, locus1.geoStart.latitude, locus1.geoStart.longitude);
			commonStartDist = locus1.startDist;

			newErr |= createPt(&commonPt2, locus1.geoEnd.latitude, locus1.geoEnd.longitude);
			commonEndDist = locus1.endDist;
		}

		newErr |= createLocus((Locus*)commonShape->this_shape, commonPt1, commonPt2, commonStartDist, commonEndDist, SEGMENT, tol, eps);
	}
	/* check the cases where one start/end point is common and loci do not coincide.
	 * then the common start/end point is the intersection
	 */
	else if ( loc1EndEqualsLoc2Start && !loc1StartOnLoc2 && !loc2EndOnLoc1 )
	{
		commonShape->type = LLPOINT;
		commonShape->this_shape = malloc(sizeof(LLPoint));
		newErr |= createPt((LLPoint*)commonShape->this_shape, locus1.locusEnd.latitude, locus1.locusEnd.longitude);
	}
	else if ( loc1StartEqualsLoc2End && !loc1EndOnLoc2 && !loc2StartOnLoc1 )
	{
		commonShape->type = LLPOINT;
		commonShape->this_shape = malloc(sizeof(LLPoint));
		newErr |= createPt((LLPoint*)commonShape->this_shape, locus1.locusStart.latitude, locus1.locusStart.longitude);
	}
	else if (loc1EndEqualsLoc2End && !loc1StartOnLoc2 && !loc2StartOnLoc1)
	{
		commonShape->type = LLPOINT;
		commonShape->this_shape = malloc(sizeof(LLPoint));
		newErr |= createPt((LLPoint*)commonShape->this_shape, locus1.locusEnd.latitude, locus1.locusEnd.longitude);
	}
	else if (loc1StartEqualsLoc2Start && !loc1EndOnLoc2 && !loc2EndOnLoc1)
	{
		commonShape->type = LLPOINT;
		commonShape->this_shape = malloc(sizeof(LLPoint));
		newErr |= createPt((LLPoint*)commonShape->this_shape, locus1.locusStart.latitude, locus1.locusStart.longitude);
	}
	/* check the cases where no end point is common, loci do not coincide, and
	 * one start/end point of only one loci lies on the other.
	 * then this start/end point is the intersection
	 */
	else if (loc1StartOnLoc2 && !loc1EndOnLoc2 && !loc2StartOnLoc1 && !loc2EndOnLoc1)
	{
		commonShape->type = LLPOINT;
		commonShape->this_shape = malloc(sizeof(LLPoint));
		newErr |= createPt((LLPoint*)commonShape->this_shape, locus1.locusStart.latitude, locus1.locusStart.longitude);
	}
	else if (!loc1StartOnLoc2 && loc1EndOnLoc2 && !loc2StartOnLoc1 && !loc2EndOnLoc1)
	{
		commonShape->type = LLPOINT;
		commonShape->this_shape = malloc(sizeof(LLPoint));
		newErr |= createPt((LLPoint*)commonShape->this_shape, locus1.locusEnd.latitude, locus1.locusEnd.longitude);
	}
	else if (!loc1StartOnLoc2 && !loc1EndOnLoc2 && loc2StartOnLoc1 && !loc2EndOnLoc1)
	{
		commonShape->type = LLPOINT;
		commonShape->this_shape = malloc(sizeof(LLPoint));
		newErr |= createPt((LLPoint*)commonShape->this_shape, tempLoc2.locusStart.latitude, tempLoc2.locusStart.longitude);
	}
	else if (!loc1StartOnLoc2 && !loc1EndOnLoc2 && !loc2StartOnLoc1 && loc2EndOnLoc1)
	{
		commonShape->type = LLPOINT;
		commonShape->this_shape = malloc(sizeof(LLPoint));
		newErr |= createPt((LLPoint*)commonShape->this_shape, tempLoc2.locusEnd.latitude, tempLoc2.locusEnd.longitude);
	}
	/*
	 * Finally, handle the cases where an intersection exists (such that the intersection
	 * point is not a start/end point of either locus) or no intersection exists.
	 * First find the expected intersection point.  If this intersection lies on both
	 * loci then we have an intersection.  Otherwise we have no intersection.
	 */
	else
	{
		/*
		 * The algorithm, in general, is to compute the standard line equation for locus 1,
		 * the standard line equation for locus 2, equate the two to find a solution (x).
		 * This will represent the distance from the locus 1 start point to the point on
		 * the geodesic that corresponds to the intersection point.  From there,
		 * the intersection point on the loci may be found.
		 * Treat the locus 1 geodesic start point as the origin of a Cartesian 2-D space.
		 */

		newErr |= invDist(locus1.geoStart, tempLoc2.geoStart, &distLoc1GeoStartToLoc2GeoStart, eps);

		/*
		 * At this point in the algorithm, cases where the slope difference is zero cannot
		 * have one intersection. Either no intersection exists or a common locus exists.
		 * IMPORTANT NOTE -
		 * This check is written in this way in order to handle the cases simultaneously where
		 * both slopes are zero and the cases where the slopes are the same in terms of magnitude
		 * and sign.  Do NOT use absolute values here.
		 */
		if (locus1.slope - tempLoc2.slope == 0)
			distIntX = 0;
		else
			distIntX = (tempLoc2.startDist - locus1.startDist - tempLoc2.slope * distLoc1GeoStartToLoc2GeoStart) / (locus1.slope - tempLoc2.slope);
                if (fabs(distIntX) > MaxDist)
		{
			//no intersection within a quarter circumference
			*err |= newErr;
			return 0;
		}

		newErr |= direct(locus1.geoStart, locus1.geoAz, distIntX, &intXGeoPt, eps);

		newErr |= ptOnLocusFromGeoPt(locus1, intXGeoPt, &intersection, &dummy1, tol, eps);

		intXOnLoc1 = ptIsOnLocus(locus1, intersection, &dummy5, &newErr, tol, eps);

		intXOnLoc2 = ptIsOnLocus(tempLoc2, intersection, &dummy5, &newErr, tol, eps);

		if (intXOnLoc1 && intXOnLoc2)
		{
			commonShape->type = LLPOINT;
			commonShape->this_shape = malloc(sizeof(LLPoint));
			newErr |= createPt((LLPoint*)commonShape->this_shape, intersection.latitude, intersection.longitude);
		}
		else
		{
			//no intersection found
			*err |= newErr;
			return 0;
		}
	}

	*err |= newErr;
   	return 1;

}
} //namespace
