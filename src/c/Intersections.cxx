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

/******************************************************************************
 * Once approximate line intersection is found by crsIntx or
 * geoIntx, then this function iterates to refine the solution.
 */
static ErrorSet iterateGeoIntx(LLPoint pt1, double crs13, double* crs31,
                                    double* dist13, LLPoint pt2, double crs23,
                                    double* crs32, double* dist23,
                                    LLPoint* intx, double* error, double tol,
                                    double eps)
{

    LLPoint newPt;
    LLPoint orgPt1;
    LLPoint orgPt2;

    ErrorSet err = 0;

    double dist12;
    double acrs13, acrs31; /* Course & dist from pt1 to approximate point newPt */
    double acrs23, acrs32; /* Course & dist from pt2 to approximate point newPt */
    //  double error;
    double minError = 1e99;
    double distarray[2] = { 0.0, 0.0 };
    double errarray[2] = { 0.0, 0.0 };
    double startPointShiftDist = 1.0; /* Dist to move start point if intersection is too close */

    int k = 0;

    int swapped = 0;
    int recalcAngles = 0;
    int start1Shifted = 0;
    int start2Shifted = 0;

    err |= invDist(pt1, *intx, dist13, eps);
    if (*dist13 <= 1.0 / 1852.0) /* Within 1 meter */
    {
        /* pt1 is near intersection */
        /* Move pt1 along back course of crs13 to make solution easier to find */
        /*    *intx = pt1; */
        err |= direct(pt1, modcrs(crs13 + M_PI), startPointShiftDist,
                &newPt, eps);

        err |= inverse(newPt, pt1, &acrs13, &acrs31, dist13, eps);

        orgPt1 = pt1;
        pt1 = newPt;
        crs13 = acrs13;
        start1Shifted = 1;
        recalcAngles = 1;
    }

    err |= invDist(pt2, *intx, dist23, eps);

    if (*dist23 <= 1.0 / 1852.0) /* Within 1 meter */
    {
        /* pt2 is near intersection */
        /*    *intx = pt2; */

        //        printf("Moving pt2 away from intersection \n");
        err |= direct(pt2, modcrs(crs23 + M_PI), startPointShiftDist,
                &newPt, eps);

        err |= inverse(newPt, pt2, &acrs23, &acrs32, dist23, eps);

        orgPt2 = pt2;
        pt2 = newPt;
        crs23 = acrs23;
        recalcAngles = 1;
        start2Shifted = 1;
    }

    if (fabs(*dist23) < fabs(*dist13))
    {
        /* swap points so that we track along the shorter course */
        newPt = pt1;
        pt1 = pt2;
        pt2 = newPt;
        acrs13 = crs13;
        crs13 = crs23;
        crs23 = acrs13;
        *dist13 = *dist23; /* don't need to keep adist23 */
        swapped = 1;
    }

    //Initialize the first distance & error guess
    distarray[0] = *dist13;
    err |= direct(pt1, crs13, distarray[0], intx, eps);
    err |= inverse(pt2, *intx, &acrs23, &acrs32, dist23, eps);
    errarray[0] = modlon(acrs23 - crs23);

    //Initialize the next distance & error guess
    distarray[1] = 1.01 * (*dist13);
    err |= direct(pt1, crs13, distarray[1], intx, eps);
    err |= inverse(pt2, *intx, &acrs23, &acrs32, dist23, eps);
    errarray[1] = modlon(acrs23 - crs23);
    if (err != SUCCESS)
        return err;

    while ((k == 0) || ((*error > tol) && (k < MAX_ITERATIONS)))
    {

        if ((fabs(distarray[0] - distarray[1]) < ONE_NANOMETER_NMI) || (fabs(
                errarray[0] - errarray[1]) < INTERNAL_ZERO))
        {
            /* If our linear approximation breaks down, exit loop */
            k = MAX_ITERATIONS; //TODO This will force an unhelpful error code to be raised. Need to change this to throw a useful error.
        }
        else
        {
            *dist13 = findRootSecantMethod(distarray, errarray, &err);

            /* Move newPt1 to line along course from pt1 */
            err |= direct(pt1, crs13, *dist13, intx, eps);
            if (err != SUCCESS)
                return err;

            /* calc distance/bearings from start of 2nd geodesic to approx intersection */
            err |= inverse(pt2, *intx, &acrs23, crs32, dist23, eps);
            if (err != SUCCESS)
                return err;

            /* place point on 2nd geodesic at same distance */
            err |= direct(pt2, crs23, *dist23, &newPt, eps);
            if (err != SUCCESS)
                return err;

            /* measure distance from approx point to similar point on 2nd geo,
             * this is the "error" in the current intersection estimate */
            err |= invDist(*intx, newPt, error, eps);
            if (err != SUCCESS)
                return err;

            //if (*error <= 1.0)
            //{
            //	//If intx and newPt are within 1 NM of each other,
            //	//start using a more accurate error prediction.
            //	*error = *error/fabs(sin(modlon(crs32-crs31)));  //developed by MJM
            //}

            if (fabs(*error) < minError)
                minError = *error;

            distarray[0] = distarray[1];
            errarray[0] = errarray[1];

            distarray[1] = *dist13;
            errarray[1] = modlon(acrs23 - crs23);
            //            printf("k = %d; dist = %20.15f; crserror = %2.15e; disterr = %2.15e\n",k,distarray[1],
            //                   errarray[1],*error);

            k++;

        }
    }

    if (k >= MAX_ITERATIONS)
    {
        //        printf("iterateGeoIntx did not converge to within %e tolerance. "
        //               "Minimum error was %e.\n",tol,*error);
        err |= SEC_NOT_CONVERGED_ERR;
    }

    if (fabs(*error) >= MAX_DISTANCE_ERROR)
    {
        err |= ERROR_MAX_REACHED_ERR;
    }

    // Force one more iteration. Testing has shown that the above
    // iteration sometimes fails to get within TOL of the actual
    // intersection point. But, testing has also shown that one
    // more iteration provides a highly accurate solution.
    if (err == SUCCESS)
    {
        //Force one more iteration
        *dist13 = findRootSecantMethod(distarray, errarray, &err);

        /* Move newPt1 to line along course from pt1 */
        err |= direct(pt1, crs13, *dist13, intx, eps);
        if (err != SUCCESS)
            return err;

        /* calc distance/bearings from start of 2nd geodesic to approx intersection */
        err |= inverse(pt2, *intx, &acrs23, crs32, dist23, eps);
        if (err != SUCCESS)
            return err;

        /* place point on 2nd geodesic at same distance */
        err |= direct(pt2, crs23, *dist23, &newPt, eps);
        if (err != SUCCESS)
            return err;

        /* measure distance from approx point to similar point on 2nd geo,
         * this is the "error" in the current intersection estimate */
        err |= invDist(*intx, newPt, error, eps);
        if (err != SUCCESS)
            return err;
    }
    //end sb chng

    err |= invCrs(pt1, *intx, &acrs13, crs31, eps);

    if (start1Shifted)
    {
        err |= invDist(orgPt1, *intx, dist13, eps);
    }
    if (start2Shifted)
    {
        err |= invDist(orgPt2, *intx, dist23, eps);
    }

    if (swapped)
    {
        /* unswap before reporting out */
        acrs31 = *crs31;
        *crs31 = *crs32;
        *crs32 = acrs31;
        newPt = pt1;
        pt1 = pt2;
        pt2 = newPt;
        //    acrs13 = crs13;
        //    crs13  = crs23;
        //    crs23  = acrs13;
        dist12 = *dist13;
        *dist13 = *dist23;
        *dist23 = dist12;
    }

    return err;

}

/**
 * recursiveCrsIntersect is a private method called only by itself (recursively)
 * and by crsIntx.
 *
 * See Geolib Bug 26599 for more details as to why
 */
static ErrorSet recursiveCrsIntx(LLPoint pt1, double crs13, double* crs31,
                       double* dist13, LLPoint pt2, double crs23,
                       double* crs32, double* dist23, int recursion, LLPoint* intx,
                       double tol, double eps)
{

    ErrorSet err = 0;
    int pt1IsOnGeo2 = 0, pt2IsOnGeo1 = 0;

    LLPoint newPt1, newPt2;
    //    LLPoint end1,end2;

    Vector pt1ECEF, intxECEF;
    LLPoint pt1TenNm, pt2TenNm; //sb
    double crs12, crs21, dist12;
    double angle1, angle2, angle3;
    double tmpCrs13 = crs13;
    double tmpCrs23 = crs23;
    double locCrs31 = 0.0, locDist13 = 0.0;
    double locCrs32 = 0.0, locDist23 = 0.0;
	double crs13Rotation = 0.0;
	double crs23Rotation = 0.0;
	double pushBackDist = 125.0;//nmi, IMPORTANT - Please see Geolib Bug 26599 before changing this value - jamezcua
    double error;
    double argument;
    double sinDist;
    int backedUp = 0;
	LLPoint perpPt2;
	double perpDist, crsPerpToPt2;
	double crsPt2ToPerp;
	double beta;
	double intDist;
	double crsPerpPt2ToPt1;

#ifdef USE_BEST_FIT_ROC

    double avgLat = 0.5*fabs(pt1.latitude + pt2.latitude);
    double sphereRad = lookUpROC(geocentricLat(avgLat));
#else

    double sphereRad = SPHERE_RADIUS;
#endif


    if (ptsAreSame(pt1, pt2, tol))
    {
        *intx = pt1;
    }
    else
    {
        //pt1IsOnGeo2 = ptIsOnCrs(pt2, crs23, pt1, &crs12, &dist12, &err,
        //        tol, eps);
        //pt2IsOnGeo1 = ptIsOnCrs(pt1, crs13, pt2, &crs21, &dist12, &err,
        //        tol, eps);
        //        err |= inverse(pt1,pt2,&crs12,&crs21,&dist12,eps);

        //sb -- use PtIsOnLine instead
        err |= direct(pt1, crs13, 10, &pt1TenNm, eps); //project out 10NM from pt1 to create a fake end point for this INFINITE line
        err |= direct(pt2, crs23, 10, &pt2TenNm, eps); //project out 10NM from pt2 to create a fake end point for this INFINITE line
        pt1IsOnGeo2 = ptIsOnGeo(pt2, pt2TenNm, pt1, INFINITE, &err, tol,
                eps);
        pt2IsOnGeo1 = ptIsOnGeo(pt1, pt1TenNm, pt2, INFINITE, &err, tol,
                eps);
        err |= inverse(pt1, pt2, &crs12, &crs21, &dist12, eps);
        if (err != SUCCESS)
            return err;
        //end sb chng

        /* Angles between pt1-pt2 geodesic and geodesics to be intersected */
        angle1 = modlon(crs13 - crs12);
        angle2 = modlon(crs21 - crs23);

//        /* The values of angle1 and angle2 can be used as discriminators for
//         * difficult geometries. In particular, geometries in which pt2 is nearly
//         * on geodesic 1 or in which pt1 is nearly on geodesic 2 need to be
//         * filtered out. This is numerically a hard geometry for the
//         * spherical intersection code to handle.
//         * Case 1
//         * if angle1 is very small, then the course along geodesic 1 is very close
//         * to the course from pt1 to pt2. Use pt2TenNm instead of pt2 to calculate
//         * the angles.
//         * Case 2
//         * if angle2 is very small, then the course along geodesic 2 is very close
//         * to the course from pt2 to pt1. Use pt1TenNm instead of pt1.
//         * Case 3
//         * if both angle1 AND angle 2 are both small, then use pt1TenNm and pt2TenNm.
//         */
//        if (fabs(angle1) < 1e-5 && fabs(angle2) > 1e-5)
//        {
//            // Case 1. Recalculate courses using pt2TenNm
//            err |= inverse(pt1, pt2TenNm, &crs12, &crs21, &dist12, eps);
//            angle1 = modlon(crs13 - crs12);
//            angle2 = modlon(crs21 - crs23);
//        }
//        else if (fabs(angle1) > 1e-5 && fabs(angle2) < 1e-5)
//        {
//            // Case 2. Recalculate courses using pt2TenNm
//            err |= inverse(pt1TenNm, pt2, &crs12, &crs21, &dist12, eps);
//            angle1 = modlon(crs13 - crs12);
//            angle2 = modlon(crs21 - crs23);
//        }
//        else if (fabs(angle1) < 1e-5 && fabs(angle2) < 1e-5)
//        {
//            // Case 3. Recalculate courses using pt1TenNm and pt2TenNm
//            err |= inverse(pt1TenNm, pt2TenNm, &crs12, &crs21, &dist12,
//                    eps);
//            angle1 = modlon(crs13 - crs12);
//            angle2 = modlon(crs21 - crs23);
//        }
//
//        sinDist = sin(dist12 / SPHERE_RADIUS);

        /* Place end points on lines */
        if (pt1IsOnGeo2 && pt2IsOnGeo1)
        {
            /* collinear courses */
            err |= COLLINEAR_COURSE_ERR;
            /* CANNOT CONTINUE */
            return err;
        }
        else if (pt1IsOnGeo2)
        {
            /* Intersection is pt1 */
            *intx = pt1;
            locDist23 = dist12;
            locCrs32 = crs12;
        }
        else if (pt2IsOnGeo1)
        {
            /* Intersection is pt2 */
            *intx = pt2;
            locCrs31 = crs21;
            locDist13 = dist12;
        }
        else
        {
        	//check for maximum allowed recursions
        	if(recursion < MAX_RECURSIONS){

				if (sin(angle1) * sin(angle2) < 0)
				{
					/* Given courses do not cross (they lie on opposite side of pt1-pt2 line)
					 * Use reciprocal courses so that nearest intersection is found */
					if (fabs(angle1) > fabs(angle2))
					{
						crs13 = modcrs(crs13 + M_PI);
						angle1 = modlon(crs13 - crs12);
					}
					else
					{
						crs23 = modcrs(crs23 + M_PI);
						angle2 = modlon(crs21 - crs23);
					}
				}
				/* Spherical solution for first approximation */
				angle1 = fabs(angle1);
				angle2 = fabs(angle2);
				angle3 = acos(-cos(angle1) * cos(angle2) + sin(angle1)
						* sin(angle2) * cos(dist12 / sphereRad));

				//Check for a divide by zero. If this occurs, one of the incoming
				//points is probably at a pole where crs is always PI.
				if (angle1 == 0.0 || angle3 == 0.0)
				{
					err |= UNEXPECTED_ERR;
					return err;
				}
				argument = (cos(angle2) + cos(angle1) * cos(angle3)) / (sin(angle1)
						* sin(angle3));

				if (fabs(argument) <= 0.9999)//IMPORTANT-Please do not change this value w/o consulting Mike Mills
				{
					locDist13 = sphereRad * acos(argument);
					err |= direct(pt1, crs13, locDist13, intx, eps);

					pt1ECEF = geodeticToECEF(pt1);
					intxECEF = geodeticToECEF(*intx);

					if (dot(pt1ECEF, intxECEF) < 0)
					{
						/* Nearest intersection is in other direction, so take
						 * supplement of courses and recalculate approximate
						 * intersection */

						crs13 = modcrs(crs13 + M_PI);
						crs23 = modcrs(crs23 + M_PI);
						angle1 = fabs(modlon(crs13 - crs12));
						angle2 = fabs(modlon(crs21 - crs23));

						angle3 = acos(-cos(angle1) * cos(angle2) + sin(angle1)
								* sin(angle2) * cos(dist12 / sphereRad));

						locDist13 = sphereRad * acos((cos(angle2) + cos(angle1)
								* cos(angle3)) / (sin(angle1) * sin(angle3)));

						err |= direct(pt1, crs13, locDist13, intx, eps);

					}

					/* Iterate to improve approximate solution */
					err |= iterateGeoIntx(pt1, crs13, &locCrs31,
							&locDist13, pt2, crs23, &locCrs32, &locDist23, intx,
							&error, tol, eps);

				}
				else if (fabs(argument) > 0.9999 && fabs(argument) < 1.0)//IMPORTANT-Please do not change these values w/o consulting Mike Mills or Juan Amezcua
				{
					/* Back up start points of both line segments to create more
					 * distance between starts and intersection point.  This appears
					 * to increase numerical stability. */
					//back up the start points in such a way that the distance between them grows
					//there are eight possible geometries, with respect to location and direction


					/* Check if the courses are relatively opposite from each other.
					 * Note, we don't need to handle the = (equals) case since that implies perfect
					 * collinearity, which should have been handled earlier.
					 */
					if ( fabs(crs13 - crs23) < M_PI_2)
					{
						//course are in relatively the same direction, invert one course
						crs23Rotation = M_PI;
					}


					/* Push back the start points by inverting the courses plus the addition
					 * of second rotation if necessary to force the push back courses to be in
					 * relatively opposite directions.
					 */
					err |= direct(pt1, crs13 + crs13Rotation + M_PI, pushBackDist, &newPt1, eps);
					err |= direct(pt2, crs23 + crs23Rotation + M_PI, pushBackDist, &newPt2, eps);
					err |= invCrs(newPt1, pt1, &tmpCrs13, &locCrs31, eps);
					err |= invCrs(newPt2, pt2, &tmpCrs23, &locCrs32, eps);
					/* !! Recursion !! */
					//increase the reursion count by 1
					err |= recursiveCrsIntx(newPt1, tmpCrs13, &locCrs31,
							&locDist13, newPt2, tmpCrs23, &locCrs32, &locDist23,
							recursion+1,intx, tol, eps);
					backedUp = 1;
				}
				else
				{
					/* Back up start points of both line segments to create more
					 * distance between starts and intersection point.  This appears
					 * to increase numerical stability. */

					//adjust the start points
					err |= direct(pt1, crs13 + M_PI, pushBackDist, &newPt1, eps);
					err |= direct(pt2, crs23 + M_PI, pushBackDist, &newPt2, eps);
					err |= invCrs(newPt1, pt1, &tmpCrs13, &locCrs31, eps);
					err |= invCrs(newPt2, pt2, &tmpCrs23, &locCrs32, eps);
					/* !! Recursion !! */
					//increase the recursion count by 1
					err |= recursiveCrsIntx(newPt1, tmpCrs13, &locCrs31,
							&locDist13, newPt2, tmpCrs23, &locCrs32, &locDist23,
							recursion+1,intx, tol, eps);
					backedUp = 1;
				}
        	}
        	else
        	{
        		//exceeded maximum number of recursions, use flat earth approximation

        		//using backed up points at this recursion, find perpendicular projection of one point onto the geodesic of the other
        		err |= projectToGeo(pt1, crs13, pt2, &perpPt2, &crsPerpToPt2, &perpDist, tol, eps);

        		//find course from pt2 to it's projection
        		err |= invCrs(pt2, perpPt2, &crsPt2ToPerp, NULL, eps);

				/* Check if the courses are relatively opposite from each other. */
        		crs23Rotation = 0.0;
				if ( fabs(crs13 - crs23) > M_PI_2)
				{
					//courses are in relatively the opposite directions,
					//invert one course to force them to be in the same relative direction
					crs23Rotation = M_PI;
				}

        		//find the magnitude of the smallest subtended angle between crs23 and crsPt2ToPerp,
				//making sure to use two courses that point in relatively the same direction
        		err |= minSubtendedAngle(crs23+crs23Rotation, crsPt2ToPerp, &beta);

        		//compute the flat earth intersection using trigonometry
        		intDist = perpDist / tan(M_PI_2 - beta);
         		err |= invCrs(perpPt2, pt1, &crsPerpToPt2, NULL, eps);
        		err |= direct(perpPt2, crsPerpPt2ToPt1, intDist, intx, eps);

        	}
        }
    }

    /* Copy local variables to external references, if available */
    /* If recursion was used to move input points, then need to calculate
     * correct return values for distance and azimuth. */
    if ((crs31 != NULL) || (dist13 != NULL))
    {
        if (backedUp)
        {
            err |= inverse(*intx, pt1, crs31, NULL, dist13, eps);
        }
        else
        {
            if (crs31)
                *crs31 = locCrs31;
            if (dist13)
                *dist13 = locDist13;
        }
    }

    if ((crs32 != NULL) || (dist23 != NULL))
    {
        if (backedUp)
        {
            err |= inverse(*intx, pt2, crs32, NULL, dist23, eps);
        }
        else
        {
            if (crs32)
                *crs32 = locCrs32;
            if (dist23)
                *dist23 = locDist23;
        }
    }
    return err;

}

/*****************************************************************************/
/*
 * Find intersection of two courses
 * Uses spherical solution as first approximation, then refines it until
 * courses from start points to intersection point match the courses
 * given.
 *
 * INPUT:
 *    pt1: start point of first geodesic
 *  crs13: azimuth from pt1 to intersection point
 *  crs31: reference to azimuth from intersection point to pt1
 * dist13: reference to distance from pt1 to intersection
 *    pt2: start point of second geodesic
 *  crs23: azimuth from pt2 to intersection point
 *  crs32: reference to azimuth from intersection to pt2
 * dist23: reference to distance from pt2 to intersection
 *    tol: required accuracy tolerance
 *    eps: forward/direct algorithm tolerance
 *         (soon to be deprecated and eventually removed)
 *
 * OUTPUT:
 *   intx: Pointer to intersection point
 *
 */

ErrorSet crsIntx(LLPoint pt1, double crs13, double* crs31,
                       double* dist13, LLPoint pt2, double crs23,
                       double* crs32, double* dist23, LLPoint* intx,
                       double tol, double eps)
{
	ErrorSet err = 0;
    double npCrs31, npDist13, npCrs32, npDist23;
    LLPoint npIntx;

    /* Assign local storage if optional pointers are not provided */
	if (NULL == crs31) crs31 = &npCrs31;
	if (NULL == dist13) dist13 = &npDist13;
	if (NULL == crs32) crs32 = &npCrs32;
	if (NULL == dist23) dist23 = &npDist23;
	if (NULL == intx) intx = &npIntx;

	err = recursiveCrsIntx(pt1, crs13, crs31, dist13, pt2, crs23, crs32, dist23, 0, intx, tol, eps);
	return err;
}


/*******************************************************************************/
/*
 * Find intersection of two geodesics, defined by start point and end point.
 * Calculates bearings and calls crsIntx to find radial intersection.
 * Tests whether intersection lies on given geodesics and between their end
 * points.
 *
 *  INPUT:
 * start1: start point of first geodesic
 *   end1: end point of first geodesic
 *  crs31: reference to azimuth from intersection point to pt1
 * dist13: reference to distance from pt1 to intersection
 * start2: start point of second geodesic
 *   end2: end point of second geodesic
 *  crs32: reference to azimuth from intersection to pt2
 * dist23: reference to distance from pt2 to intersection
 *   intx: Pointer to intersection point
 *    tol: required accuracy tolerance
 *    eps: forward/direct algorithm tolerance
 * OUTPUT:
 * 	  error code > 0 if no intersection is found
 */

ErrorSet geoIntx(LLPoint start1, LLPoint end1, LineType lineType1,
                        double* crs31, double* dist13, LLPoint start2,
                        LLPoint end2, LineType lineType2, double* crs32,
                        double* dist23, LLPoint* intx, double tol, double eps)
{

    ErrorSet err = 0;

    double crs13;
    double crs23;
    int endPointMatch = 0;

    double npCrs31, npDist13, npCrs32, npDist23;
    LLPoint npIntx;

    /* Assign local storage if optional pointers are not provided */
	if (NULL == crs31) crs31 = &npCrs31;
	if (NULL == dist13) dist13 = &npDist13;
	if (NULL == crs32) crs32 = &npCrs32;
	if (NULL == dist23) dist23 = &npDist23;
	if (NULL == intx) intx = &npIntx;

    /* Check whether start intersection occurs at start/end points */
    if (ptsAreSame(start1, start2, tol) || ptsAreSame(start1,
            end2, tol))
    {
        *intx = start1;
        endPointMatch = 1;
    }
    else if (ptsAreSame(end1, end2, tol) || ptsAreSame(end1,
            start2, tol))
    {
        *intx = end1;
        endPointMatch = 1;
    }

    if (endPointMatch)
    {
        err |= inverse(*intx, start1, crs31, &crs13, dist13, eps);
        err |= inverse(*intx, start2, crs32, &crs23, dist23, eps);
    }
    else
    {
        err |= inverse(start1, end1, &crs13, crs31, dist13, eps);
        err |= inverse(start2, end2, &crs23, crs32, dist23, eps);

        if (err |= crsIntx(start1, crs13, crs31, dist13,
                start2, crs23, crs32, dist23, intx, tol, eps))
            return err; /* Do not continue if no intersection found */
//        displayMatlabPt(*intx,"intx",0);

        /* If intx is not on both input lines in accordance with lineType,
         * then return error code.  Passed reference intx IS NOT set to NULL. */
        if (!(ptIsOnGeo(start1, end1, *intx, lineType1, &err, tol, eps)
                && ptIsOnGeo(start2, end2, *intx, lineType2, &err, tol,
                        eps)))
        {
            err |= NO_INTERSECTION_ERR;
            return err;
        }
    }

    return err;

}

/*******************************************************************************
 * Arc-Arc Intersection
 *
 * n = number of intersections (0, 1, 2).  If n = -1 then unable to determine intersection
 *
 */

ErrorSet arcIntx(LLPoint center1, double r1, LLPoint center2, double r2,
                       LLPointPair intx, int* n, double tol, double eps)
{

    ErrorSet err = 0;

    LLPoint pt, tempPt;

    double crs12, crs21, dist12, tempr;
    double crs1x, crsx1, crs2x, crsx2, dist1x, dist2x;
    double crsarray[2], errarray[2], longarray[2];
    int arc1IsAtPole = 0, arc2IsAtPole = 0;
    double error, y;
    double intersectionLatitude, newLongitude;
    Arc arcP, arcNP;

    int sn; /* Number of intersecting points found from sphere approx. */
    //    int narray; /* number of approximations stored in crsarray, errarray   */
    int k; /* Iteration number */
    int i;
    int c1EqualsC2 = -1;
	int c1OnArc2 = -1;
	int c2OnArc1 = -1;

	int npN;
	LLPointPair npIntx;

    /* Assign local storage if optional pointers are not provided */
    if (NULL == n) n = &npN;
    if (NULL == intx) intx = npIntx;


    /* initialize count of intersections */
    *n = -1;

    /* Arcs with radius greater than approximate half-circumference of the Earth
     * are not defined and can't be solved with this algorithm. */
    if ((r1 > MAX_ELLIPSOIDAL_ARC_RADIUS_NMI) ||
    		(r2 > MAX_ELLIPSOIDAL_ARC_RADIUS_NMI))
    {
        err |= RADIUS_OUT_OF_RANGE_ERR;
        return err;
    }
    /*The algorithm currently implemented cannot operate on
     *arcs which are centered at a pole. Dump if this occurs. */
    //else if (ptIsAtPole(center1,&err,tol,eps) != 0 ||
    //         ptIsAtPole(center2,&err,tol,eps) != 0)
    //{
    //    err |= UNEXPECTED_ERR;
    //    return err;
    //}

    //TODO update algorithm to not allow for an arc to be a singular point? - jamezcua
    //If allowed, this could lead to a logical inconsistency with respect to the geometry (including points) that is defined over a given mathematical space
    /*
     * Check for radius1 == 0, radius2 == 0, dist12 == 0
     */
    //handle arc(s) that are singular points
    if ((r1 == 0.0) && (r2 == 0.0))
    {
    	//Both arcs are actually singular points in space
    	//Use the arc center points as the representations of the arcs themselves

    	//check if the centers are identical
    	c1EqualsC2 = ptsAreSame(center1, center2, tol);

    	if (c1EqualsC2)
    	{
    		//centers are the same and both arcs are singular points
    		//therefore the intersection of the points are the points themselves
    		//pick center1 as the intersection point, with respect to the mathematical notion of neighborhoods
    		*n = 1;
    		err |= createPt(&intx[0], center1.latitude, center1.longitude);
    	}
    	else
    	{
    		//centers are the same and both arcs are singular points
    		//since the centers are distinct with respect to tol (and the mathematical notion of neighborhoods)
    		//no intersection may exist
    		*n = 0;
    	}

    	return err;
    }
    else if (r1 == 0.0)
    {
    	//Only arc1 is a singular point
    	//Use the center1 point as the representation of arc1

    	//check if center1 lies on arc2 (treating arc2 as a circle)
    	c1OnArc2 = ptIsOnArc(center2, r2, 0.0, 0.0, ArcDirection::COUNTERCLOCKWISE, center1, &err, tol, eps);

    	if (c1OnArc2)
    	{
    		//center1 lies on arc2 so center1 must be the intersection since center1 is also a singular point
    		*n = 1;
    		err |= createPt(&intx[0], center1.latitude, center1.longitude);
    	}
    	else
    	{
    		//center1 does not lie on arc2 and center1 is also a singular point
    		//no intersection may exist
    		*n = 0;
    	}

    	return err;
    }
    else if (r2 == 0.0)
    {
    	//Only arc2 is a singular point
    	//Use the center2 point as the representation of arc2

    	//check if center2 lies on arc1 (treating arc1 as a circle)
    	c2OnArc1 = ptIsOnArc(center1, r1, 0.0, 0.0, ArcDirection::COUNTERCLOCKWISE, center2, &err, tol, eps);

    	if (c2OnArc1)
    	{
    		//center2 lies on arc1 so center1 must be the intersection since center2 is also a singular point
    		*n = 1;
    		err |= createPt(&intx[0], center2.latitude, center2.longitude);
    	}
    	else
    	{
    		//center2 does not lie on arc1 and center2 is also a singular point
    		//no intersection may exist
    		*n = 0;
    	}

    	return err;
    }

    /*
     * determine the number intersections.
     * Check radii to figure out number of intersections
     */

     // Ensure that circle1 is not smaller than circle2
     // This step simplifies the check for non-intersecting cases later.
     if (r2 > r1)
     {
         tempr = r1;
         tempPt = center1;
         r1 = r2;
         center1 = center2;
         r2 = tempr;
         center2 = tempPt;
     }

     //find the distance between the arc centers
     err |= inverse(center1, center2, &crs12, &crs21, &dist12, eps);

     //at this point in the algorithm, it is guaranteed that r1 >= r2
    if (dist12 <= tol)
    {
    	//circles are concentric

    	if (fabs(r1 - r2) <= tol)
    	{
    		//circles are coincidental, infinite intersections
    		err |= CONCENTRIC_CIRCLE_ERR;
    		//allow n to = -1 indicating that no distinct intersections exist
    	}
    	else
    	{
    		//circles are concentric but not coincidental, no intersections
    		*n = 0;
    	}
    }
    /*
     * Please do not update this section of the code if it is unclear what is going on mathematically, see Mike or Juan for details
     */
    //Six "simple" cases - the use of tolerance complicates the problem
    if ( r1 + tol + r2 + tol < dist12)//see the figures/details in bug 18763 to understand this inequality - jamezcua
    {
        //circles are widely separated, no intersections
    	*n = 0;
    }
    else if (fabs((r1 - tol) - (r2 + tol)) > dist12)//see the figures/details in bug 18763 to understand this inequality - jamezcua
    {
    	//one circle lies within the other such that no intersection exists
    	*n = 0;
    }
    else if ( fabs(r1 + r2 - dist12) <= tol)//TODO need to verify the math on this inequality - jamezcua
    {
    	//The discs, represented by the circles, are disjoint other than at the point of tangency
    	*n = 1;
    }
    else if (fabs(fabs(r1 - r2) - dist12) <= tol)//TODO need to verify the math on this inequality - jamezcua
    {
    	//one circle lies within the other such that the inner circle is tangent to the outer circle
    	*n = 1;
    }
    else if (r1 + r2 > dist12)//TODO need to verify the math on this inequality - jamezcua
    {
    	//The discs, represented by the circles, overlap such that no center point overlaps either disc
    	*n = 2;
    }
    else if (fabs(r1 - r2) < dist12)//TODO need to verify the math on this inequality - jamezcua
    {
    	//The discs, represented by the circles, overlap such that at least one center point overlaps a disc
    	*n = 2;
    }




    //find the intersections if they exist
    if (*n == 1)
    {
    	err |= direct(center1, crs12, r1, &intx[0], eps);
    }
    else if (*n == 2)
    {

		// If both arcs are at a pole
		arc1IsAtPole = ptIsAtPole(center1, &err, tol, eps) != 0;
		arc2IsAtPole = ptIsAtPole(center2, &err, tol, eps) != 0;
		if (arc1IsAtPole && arc2IsAtPole)
		{
			//Umm...okay, what was the caller thinking?! Either no
			//intersections (if radii not equal) or infinite intersections.
			err |= UNEXPECTED_ERR;
			return err;
		}
		// If one arc is at a pole
		else if (arc1IsAtPole || arc2IsAtPole)
		{
			/* Algorithm Description:
			 We know the latitude at which the intersection must
			 occur (constant latitude of the arc at the pole). So,
			 we can just iterate along that constant latitude line
			 until we find a longitude that is on the non-pole arc.
			 Once one intersection is known, the other is simply on
			 the opposite side of the non-pole arc's center point.
			 */


			// Figure out which one is at the pole
			if (arc1IsAtPole)
			{
				arcP.centerPoint = center1;
				arcP.radius = r1;
				arcNP.centerPoint = center2;
				arcNP.radius = r2;
			}
			else
			{
				arcP.centerPoint = center2;
				arcP.radius = r2;
				arcNP.centerPoint = center1;
				arcNP.radius = r1;
			}

			// Calculate the intersection latitude from the pole arc
			err |= direct(arcP.centerPoint, 0, arcP.radius, &tempPt, eps);
			intersectionLatitude = tempPt.latitude;

			// Create a point on the latitude that we
			// know is correct for the intersection. Guess the longitude. This will help us
			// seed the iteration loop below.
			intx[0].latitude = intersectionLatitude; //initial guess for an intersection
			intx[0].longitude = arcNP.centerPoint.longitude;
			err |= invDist(arcNP.centerPoint, intx[0], &dist1x, eps); //distance btwn center and initial guess
			error = dist1x - arcNP.radius; //error of initial guess
			longarray[0] = intx[0].longitude; //store the longitude of initial guess
			errarray[0] = error; //store the error of initial guess

			// Create another point at a slightly different longitude. But,
			// the iteration below will fail if trying to find a solution across
			// the international date line (longitude of +/-pi). So, always
			// choose to search for the intersection point that is away from
			// that line.
			intx[0].latitude = intersectionLatitude; //next guess for intersection point
			intx[0].longitude = arcNP.centerPoint.longitude + (-sgn(
					arcNP.centerPoint.longitude)) * 0.1;
			err |= invDist(arcNP.centerPoint, intx[0], &dist1x, eps); //distance from center to intersection point
			error = dist1x - arcNP.radius; //error in this guess
			longarray[1] = intx[0].longitude; //store the longitude
			errarray[1] = error; //store the error

			// Now we're ready to iterate the longitude of our guess along a
			// line of constant latitude until the distance value reaches
			// the radius value (error = 0).
			k = 0;
			while ((k < MAX_ITERATIONS) && (fabs(error) > tol) && (fabs(error)
					< MAX_ELLIPSOIDAL_ARC_RADIUS_NMI))
			{

				k++;

				// Get a new guess for the intersection longitude
				newLongitude = findRootSecantMethod(longarray, errarray, &err);
				newLongitude = modlon(newLongitude); //force onto interval [-pi,pi]

				// Use newLongitude to create the new guess point
				intx[0].longitude = newLongitude;
				err |= invDist(arcNP.centerPoint, intx[0], &dist1x, eps); //distance from center to intersection point
				error = dist1x - arcNP.radius; //error in this guess

				// Now shift & store the values into the arrays
				longarray[0] = longarray[1];
				errarray[0] = errarray[1];
				longarray[1] = intx[0].longitude;
				errarray[1] = error;
			} //while

			// Now that we have one intersection point, we can calculate the other without difficulty
			// because it is opposite the known point.
			intx[1].latitude = intersectionLatitude;
			intx[1].longitude = modlon(arcNP.centerPoint.longitude
					+ (arcNP.centerPoint.longitude - intx[0].longitude));

		}
		// Neither arc is at a pole
		else
		{
			/* Spherical approx */
			err |= initArcIntx(geodeticToGeocentric(center1), r1,
					geodeticToGeocentric(center2), r2, intx, &sn, NULL, eps);
			/* Sphere arc intersect may return NO_INTERSECTION_ERR in cases where the
			 * circles are nearly tangent.  Non-intersecting cases are already handled by
			 * the midptErr check above, so we mask out the NO_INTERSECTION_ERR and proceed,
			 * pretending that the arcs are approximately tangent (number of approximate
			 * intersections will be 0).*/

			err = getMaskedError(err, NO_INTERSECTION_ERR);
			if (err)
			{
				return err;
			}

			if (sn < 2)
			{
				/* Spherical approximation doesn't exist, try more direct approximation *
				 * This case will occur when close to the tangent case.  Due to errors  *
				 * in the spherical approximation, the spherical solution may not exist *
				 * even though the ellipsoidal solution has at least one intersection.  */

				/* Flat-earth approximation of distance from midpoint to intersection */
				/* Distance from center1 to point along center1-center2 line abeam intersection */
				y = (r2 * r2 - r1 * r1 + dist12 * dist12) / (2.0 * dist12);

				/* angle between crs12 and crs to approximate intersection */
				crs1x = acos(y / r1);

				err |= direct(center1, crs12 + crs1x, r1, &intx[0], eps);
				err |= direct(center1, crs12 - crs1x, r1, &intx[1], eps);
				if (err)
				{
					return err;
				}
			}

			/* intx will point to array of two points now
			 * Refine position of each point until it lies on both circles.
			 */

			for (i = 0; i < 2; i++)
			{
				if (i == 0)
				{
					intx[i].latitude = geodeticLat(intx[i].latitude);
					pt = intx[i];
					/* Find course from center 2 to approx point */
					err |= inverse(center2, pt, &crs2x, &crsx2, &dist2x, eps);
				}
				else
				{
					/* use angle to first solution as starting point for
					 * second solution */
					crs2x = crs21 + modlon(crs21 - crs2x);
				}

				k = 0;

				/* Place point on circumference of circle 2 */

				err |= direct(center2, crs2x, r2, &pt, eps);
				/* Find distance from center 1 to pt */
				/* As pt approaches intersection, this value will approach r1 */
				err |= inverse(center1, pt, &crs1x, &crsx1, &dist1x, eps);

				error = r1 - dist1x;
				//        printf("initial error[%d] = %e => ",i,error);

				errarray[1] = error;
				crsarray[1] = crs1x;

				/* Iteratively improve solution */
				while ((k < MAX_ITERATIONS) && (fabs(error) > tol))
				{

					k++;

					/* Use calculated course to move point to circumference 1 */
					err |= direct(center1, crs1x, r1, &pt, eps);

					/* Find course from center 2 to new approx point */
					err |= inverse(center2, pt, &crs2x, &crsx2, &dist2x, eps);

					/* Move point along crs to circumference of circle 2 */
					err |= direct(center2, crs2x, r2, &pt, eps);

					/* Find distance from center 1 to pt */
					//      dist1x = invDist(center1,pt,eps);
					err |= inverse(center1, pt, &crs1x, &crsx1, &dist1x, eps);

					error = r1 - dist1x;

					crsarray[0] = crsarray[1];
					crsarray[1] = crs1x;
					errarray[0] = errarray[1];
					errarray[1] = error;
					/* The error is a function that depends on crs1x
					 * Use linear approx to this function to find new
					 * intersection estimate */
					crs1x = findRootSecantMethod(crsarray, errarray, &err);
				} /* end while */

				//        printf("k[%d] = %d, error[%d] = %e\n",i,k,i,error);
				intx[i] = pt;


			} /* end for */

		} /* if arc on pole */

		// Error check the iteration
		if (k >= MAX_ITERATIONS)
		{
			err |= SEC_NOT_CONVERGED_ERR;
		}

		if (fabs(error) >= MAX_DISTANCE_ERROR)
		{
			err |= ERROR_MAX_REACHED_ERR;
		}
    }

    return err;

}

/********************************************************************************
 * Line-Arc Intersections
 *
 */

ErrorSet geoArcIntx(LLPoint pt1, double crs1, LLPoint center,
                           double radius, LLPointPair intx, int* n, double tol,
                           double eps)
{

    ErrorSet err = 0;

    LLPoint perpPt;
    LLPoint newStart;
    int i, k, j = 0;

    double perpDist, perpCrs; // Distance & crs from center to perp. point
    double dist[2]; // Distance from perpendicular point to approx. point
    double crs[2], bcrs[2]; // Course and backcourse from perp. point to approx point
    double rcrs, brcrs; // Course and backcourse from intx to arc center point
    double distarray[2], errarray[2];
    double A, B, c, error, radDist;
    double tempCrs, newCrs; // throw-away course value
    double stepSize = 9.0e99;

#ifdef USE_BEST_FIT_ROC

    double sphereRad =
    lookUpROC(geocentricLat(0.5*fabs(pt1.latitude+center.latitude)));
#else

    double sphereRad = SPHERE_RADIUS;
#endif

	int npN;
	LLPointPair npIntx;

    /* Assign local storage if optional pointers are not provided */
    if (NULL == n) n = &npN;
    if (NULL == intx) intx = npIntx;

    *n = 0; /* initialize number of intersections found */

    /* eps (basic projection accuracy) must be less than
     * intersection accuracy tolerance                    */
    //TODO: With hard-coded eps, this check should not be necessary
    if (tol < eps)
    {
        err |= TOL_TOO_SMALL_ERR;
        return err;
    }

    /* Handle case if center and pt1 are essentially same point */
    /* Find dist from center to pt1 */
    err |= invDist(center, pt1, &radDist, eps);

    if (radDist < (tol - eps))
    {
        /* center and pt1 are indistinguishable */
        *n = 2;
        err |= direct(center, crs1, radius, &intx[0], eps);
        err |= direct(center, modcrs(crs1 + M_PI), radius, &intx[1], eps);
        return err;
    }
    else if (radDist < radius)
    {
        /* pt1 is inside circle, numerical accuracy may not be good */
        /* Move start of line outside circle for better accuracy */
        /* New start will be at least one mile from circle */
        err |= direct(pt1, crs1 + M_PI, 2.0 * radius + 10.0, &newStart, eps);

        err |= invCrs(newStart, pt1, &newCrs, &tempCrs, eps);
        pt1 = newStart;
        crs1 = newCrs;
    }

    /* project center of arc onto geodesic */
    if (err |= projectToGeo(pt1, crs1, center, &perpPt, &perpCrs,
            &perpDist, tol, eps))
        return err;
    /* calc distance & crs from pt1 to projected point */
    if (err |= invCrs(perpPt, pt1, &crs[0], &bcrs[0], eps))
        return err;
    crs[1] = modcrs(crs[0] + M_PI);

    if (fabs(perpDist - radius) < tol)
    {
        /* Line is tangent to circle */
        *n = 1;
        intx[0] = perpPt;
        return 0;
    }
    else if (perpDist > radius)
    {
        /* no intersections -- line too far from circle */
        *n = 0;
        return err;
    }
    else if (ptsAreSame(perpPt,center,tol))
    {
        /* Line goes through center of arc.  This is a special case where
         * we can project the intersection points directly. Find azimuth of geodesic
         * at arc center, then project out radius distance along this azimuth (in
         * both directions) to find two intersections with arc. */
        double azFromCenter, azToCenter;
        /* Find azimuth of geodesic at arc center */
        err |= inverse(center,pt1,&azFromCenter,&azToCenter,&radDist,eps);
        if (radDist < 100.0/1852.0)
        {
            /* line start is less than 100 meters from center point, so azimuth will not
             * be very precise.  Use precise given azimuth to move start of geodesic farther away
             * and then update pt1 to this location */
            if (fabs(modlon(crs1-azToCenter))<M_PI_2) {
                /* Arc center lies in crs1 direction from pt1, so move pt1 backward 1.0 NM */
                direct(pt1,crs1+M_PI,1.0,&pt1,eps);
            }
            else {
                /* Arc center lies "behind" pt1, so move pt1 forward 1 NM */
                direct(pt1,crs1,1.0,&pt1,eps);
            }
            /* Recompute azimuth */
            err |= inverse(center,pt1,&azFromCenter,&azToCenter,&radDist,eps);
        }

        *n = 2;
        err |= direct(center,azFromCenter,radius,&intx[0],eps);
        err |= direct(center,modcrs(azFromCenter+M_PI),radius,&intx[1],eps);

        return err;

    }

    if (cos(perpDist / sphereRad) > 0)
    {

        dist[0] = sphereRad * acos(cos(radius / sphereRad) / cos(perpDist
                / sphereRad));
    }
    else
    {
        /* Arc and geodesic describe the same great circle.
         * Intersection is entire arc. */
        //        err |= CONCENTRIC_CIRCLE_ERR;
        return err;
    }

    /* move first approximate point to line */
    if (err |= direct(perpPt, crs[0], dist[0], &intx[0], eps))
        return err;

    /* Iterate to improve approximations */

    for (i = 0; i < 2; i++)
    {

        if (i == 1)
        {
            // Use solution to first point to find approximation to second point
            dist[1] = dist[0];
            err |= direct(perpPt, crs[1], dist[1], &intx[1], eps);
        }

        k = 0;
        /* Calculate distance from center to approx. intersection point */
        if (err |= inverse(intx[i], center, &rcrs, &brcrs, &radDist, eps))
            return err;
        /* error in approximation is difference between approx. distance and
         * arc radius */
        error = radius - radDist;

        /* Preload array to enable linear extrapolation/interpolation of
         * solution */
        distarray[1] = dist[i];
        errarray[1] = error;

        /* calculate distance adjustment */
        while ((k == 0) || //force entry to ensure at least one iteration
                (((fabs(error) > tol) || (stepSize > tol)) &&
                        (k < MAX_ITERATIONS) ))
        {

            if (k == 0)
            {
                /* On first pass, improve approximation using triangle */
                if (err
                        |= invCrs(intx[i], perpPt, &bcrs[i], &tempCrs, eps))
                    return err;

                B = fabs(modlon(bcrs[i] - rcrs)); // into range [0,pi]

                /* C = 90 degrees             |\
                   a = error                  |A\
                   c = length adjustment      |  \
                                             b|   \c
                 |    \
                                              |     \
                                              |90___B\
                                                 a
                 */
                /* Formulae for spherical triangles */
                /* NOTE: this uses a great circle to approximate small portion of
                 * arc near the intersection (side b in diagram). */
                A = acos(sin(B) * cos(fabs(error) / sphereRad));
                if (fabs(sin(A)) < INTERNAL_ZERO)
                {
                    /* case where line is close to diameter */
                    c = error;
                }
                else if (fabs(A) < INTERNAL_ZERO)
                {
                    /* case where line is nearly tangent */
                    //TODO Is this a dead if condition?
                    c = error / cos(B);
                }
                else
                {
                    /* normal case */
                    c = sphereRad * asin(sin(error / sphereRad) / sin(A));
                }

                /* We move in different directions depending on which side of circle we are on */
                if (error > 0)
                {
                    /* if intx[i] is inside circle, move away from perpendicular point */
                    dist[i] += fabs(c);
                }
                else
                {
                    dist[i] -= fabs(c);
                }
            }
            else
            {
                /* Subsequent approximations use linear extrapolation/interpolation */
                /* Can't do this on first pass -- need two points for linear approx */
                dist[i] = findRootSecantMethod(distarray, errarray, &err);
                // DEBUG
                //                 printf("%20.15f, %20.15f, %20.15f, %20.15f, %20.15f \n", distarray[0],
                //                        errarray[0],distarray[1],errarray[1],dist[i]);
            }

            /* Place i-th approximate intersection point on geodesic */
            if (err |= direct(perpPt, crs[i], dist[i], &intx[i], eps))
                return err;
            /* Calculate distance from center of circle to approx. intersection */
            if (err |= inverse(intx[i], center, &rcrs, &brcrs, &radDist,
                    eps))
                return err;
            error = radius - radDist;

            /* This function uses the secant root finding method. If the secant method gets
             * stuck crossing back and forth over the root (see bug 33554) we switch to the
             * regula falsi method and insure that the two points used in the method have
             * opposite sign in the errarray. If that gets stuck we switch back to the
             * secant method.
             */
            if((sgn(errarray[0]) == sgn(errarray[1])) || (error == 0) || k < 5 || abs(j) == 2){
				distarray[0] = distarray[1];
				distarray[1] = dist[i];
				errarray[0] = errarray[1];
				errarray[1] = error;
				j = 0;
            } else if (sgn(error) == errarray[0]){
            	errarray[0] = error;
            	distarray[0] = dist[i];
            	j++;
            } else {
            	errarray[1] = error;
            	distarray[1] = dist[i];
            	j--;
            }

            stepSize = fabs(distarray[1] - distarray[0]);

            k++;

        } /* end while */

        if (k >= MAX_ITERATIONS)
        {
            err |= ITERATION_MAX_REACHED_ERR;
//            printf("Error: ITERATION_MAX_REACHED in %s\n",__FUNCTION__);

        }

        if (fabs(error) >= MAX_DISTANCE_ERROR)
        {
            err |= ERROR_MAX_REACHED_ERR;
        }

        *n = *n + 1;

    } /* end for(i) */

    return err;

}

/******************************************************************************
 * Find the intersections of a locus with an arc
 */
ErrorSet locusArcIntx(Locus orgloc, LLPoint center, double radius,
                            LLPointPair intx, int* n, double tol, double eps)
{
    ErrorSet err = 0;
    ErrorSet tmpErr;
    Locus loc = orgloc;
    LLPoint perpint;
    LLPoint locpt, pt1, locpt1, geoPt;
    double geodarray[2], errarray[2];
    double fcrs, bcrs, fcrs1, incr, crsFromPoint, distFromPoint, distbase;
    double iterationCrs, iterationInvCrs;
    double distcent, error, newdistbase, dist1;
    double phi; /* used in spherical approx to intersection points */
    double R0 = SPHERE_RADIUS;
    int i, k, n1, flag; //, online
    double locPtDistToEndPts;
    double projDist = 1.0;//nm
    LLPoint newGeoStart, newGeoEnd;
    double newStartDist, newEndDist;
    double temp;
    double stepSize = 9.0e99;
    double az12, az21, tempDist;

    *n = 0;
    incr = 1.001;

    if (err |= invCrs(loc.locusStart, loc.locusEnd, &fcrs, &bcrs, eps))
        return err;

    //    /* Find intersection of line defined by locus start and end points
    //     * with arc as a starting point for finding intersections of locus
    //     * and arc. */
    //    err |= geoArcIntx(loc.locusStart, fcrs, center, radius,
    //                                 intx, &n1, tol, eps);
    //    /* This is not quite correct, since it is possible for the arc to not intersect
    //     * the geodesic approximation, and yet still intersect the locus */

    /* Improved starting solution:
     * Project center of arc to locus.  If projection distance is > radius, then no intersection
     * If = radius, then tangent; if < radius, then two intersection candidates. */

    err |= projectToLocus(loc, center, &locpt, &crsFromPoint,
            &distFromPoint, tol, eps);

    if (err) {
        return err;
    }

    //printf("%f %f\n", distFromPoint, radius);

    if (distFromPoint > radius + tol)
    {
        /* No intersections */
        *n = 0;
        return err;
    }
    else if (fabs(distFromPoint - radius) < tol)
    {
        /* Locus is tangent to arc */
        if (ptIsOnLocus(orgloc, locpt, &perpint, &err, tol, eps))
        {
            *n = 1;
            intx[0] = locpt;
        }
        else
        {
            *n = 0;
        }
        return err;
    }
    else
    {
        /* two possible intersections */
        n1 = 2;
        err |= invDist(center, loc.locusStart, &distcent, eps);
        if (fabs(distcent - radius) < tol) {
        	intx[*n] = loc.locusStart;
        	*n = *n + 1;
        	//n1--;
        }

        err |= invDist(center, loc.locusEnd, &distcent, eps);
		if (fabs(distcent - radius) < tol) {
			intx[*n] = loc.locusEnd;
			*n = *n + 1;
			//n1--;
		}
		if (*n == n1) {
			return err;
		}
        //Check if locpt and center point are the same, if so then zero out crsFromPoint
        //This should handle the case where the arc is centered on the locus itself

		if (*n == 0) {
			if (distFromPoint <= tol)
			{
				/* approximate intersection points by moving a distance of radius along the locus from the center point */
				//err |= invCrs(loc.locusStart, center, &az12, &az21, eps);
					  az12 = locusCrsAtPt(loc,loc.locusStart,&geoPt,&temp,&err,tol,eps);
				/* one approx intersection to right of projected point (viewed from arc center) */
				err |= direct(center, az12, radius, &(intx[0]), eps);
				/*one to left */
				err |= direct(center, fmod(az12 + M_PI, M_2PI), radius, &(intx[1]), eps);
				crsFromPoint = 0;
			} else {
				/* Find approximate intersection points using spherical trig */
				/* Angle between perp projection geodesic and radial geodesic to approx points */
				phi = acos(tan(distFromPoint / R0) / tan(radius / R0));

				/* one approx intersection to right of projected point (viewed from arc center) */
				err |= direct(center, crsFromPoint + phi, radius, &(intx[0]), eps);
				/* one to left */
				err |= direct(center, crsFromPoint - phi, radius, &(intx[1]), eps);
			}
		} else {
			if (distFromPoint <= tol) {
				err |= invCrs(center, intx[0], &fcrs, &bcrs, eps);
				err |= direct(center, fcrs + M_PI, radius, &(intx[1]), eps);
			} else {
				err |= inverse(locpt, intx[0], &fcrs, &bcrs, &distFromPoint, eps);
				err |= direct(locpt, fcrs + M_PI, distFromPoint, &(intx[1]), eps);
			}
		}
    }
    /* Check if distance from locPt to locStart or locEnd is less than radius.
     * If so, then extend locus in one or both directions so that intersection point
     * is not almost abeam the geodesic's start/end point.  This will ensure that
     * distbase is large enough for accurate inv/dest calculations
     */

    //check distance between locPt and locus start point
    err |= invDist(locpt, loc.locusStart, &locPtDistToEndPts, eps);
    //checking with respect to 2*radius instead of radius to give adequate buffer
    if (fabs(locPtDistToEndPts) <= 2 * radius)
    {

        //handle "small" arcs so that we aren't extending by some small distance which may not help address the problem
        if (2 * radius <= 1)
            projDist = 1;
        else
            projDist = 2 * radius;

        //extend the locus geo start point out to ensure intersection is not "close" to locus start point
        err |= direct(loc.geoStart, loc.geoAz + M_PI, projDist,
                &newGeoStart, eps);

        //find the new locus start distance in the event that the locus is splay
        newStartDist = distToLocusFromGeoPt(orgloc, newGeoStart, &temp, &err, tol,
                eps);

    }
    else
    {
        newGeoStart = orgloc.geoStart;
        newStartDist = orgloc.startDist;
    }

    //check distance between locPt and locus end point
    err |= invDist(locpt, loc.locusEnd, &locPtDistToEndPts, eps);
    //checking with respect to 2*radius instead of radius to give adequate buffer
    if (fabs(locPtDistToEndPts) <= 2 * radius)
    {

        //handle "small" arcs so that we aren't extending by some small distance which may not help address the problem
        if (2 * radius <= 1)
            projDist = 1;
        else
            projDist = 2 * radius;

        //extend the locus geo end point out to ensure intersection is not "close" to locus end point
        err |= direct(loc.geoEnd, loc.geoRevAz + M_PI, projDist, &newGeoEnd,
                eps);

        //find the new locus end distance in the event that the locus is splay
        newEndDist
                = distToLocusFromGeoPt(orgloc, newGeoEnd, &temp, &err, tol, eps);
    }
    else
    {
        newGeoEnd = orgloc.geoEnd;
        newEndDist = orgloc.endDist;
    }

    //Update the locus in case it was extended
    err |= createLocus(&loc, newGeoStart, newGeoEnd, newStartDist,
            newEndDist, orgloc.lineType, tol, eps);

    /* find forward & reverse azimuths of locus's geodesic */
    err |= invCrs(loc.geoStart, loc.geoEnd, &fcrs1, &bcrs, eps);

    for (i = *n; i < n1; i++)
    {
        tmpErr = 0;
        k = 0;

        /* Find perp intercept of line-arc intersection point and the
         * geodesic of the locus. */
        tmpErr |= projectToGeo(loc.geoStart, fcrs1, intx[i], &perpint,
                &crsFromPoint, &distFromPoint, tol, eps);
        err |= invDist(loc.locusStart, loc.locusEnd, &tempDist, eps);
        if (tmpErr)
            continue;
        /* Compute initial error in approx point by moving it to locus */

        /* Find the distance along the locus geodesic of the perp intercept point. */
        err |= invDist(perpint, loc.geoStart, &distbase, eps);
        /* Find the point on the locus corresponding to the perp intercept point. */
        err |= ptOnLocusFromGeoPt(loc, perpint, &locpt, &crsFromPoint, tol, eps);
        //printf("locpt lat = %f lon = %f\n",locpt.latitude*180.0/M_PI, locpt.longitude*180.0/M_PI);
        err |= invDist(locpt, center, &distcent, eps);

        error = distcent - radius;
        if (fabs(error) < tol)
        {
            if ((orgloc.lineType == INFINITE) || (ptIsOnLocus(orgloc, locpt, &perpint, &tmpErr, tol, eps)))
            {
              intx[*n] = locpt;
              *n = *n + 1;
            }
            continue;
        }

        errarray[1] = error;
        geodarray[1] = distbase;
        newdistbase = incr * distbase;

        err |= invDist(intx[0], intx[1], &tempDist, eps);
        //printf("%f %f\n", newdistbase, tempDist);

        // Find the course from the locus geo start point to the perp intercept point
        err |= invCrs(loc.geoStart, perpint, &iterationCrs,
                &iterationInvCrs, eps);
        // Determine which course to iterate along

        if (abs((int) fcrs1 - (int) iterationCrs))
        {
            //perp intercept point is in the opposite direction of the locus geo start course
            iterationCrs = modcrs(fcrs1 + M_PI);
        }
        else
        {
            //perp intercept point is along the course from the locus geo start to locus geo end
            iterationCrs = fcrs1;
        }

        /* We iterate by varying a point (pt1) along the geodesic of the
         * locus until the corresponding point on the locus (locpt1) also
         * lies on the arc. */

        while ((k++ < MAX_ITERATIONS) && ((fabs(error) > tol) || (stepSize > tol)))
        {
            err |= direct(loc.geoStart, iterationCrs, newdistbase, &pt1, eps);
            err |= ptOnLocusFromGeoPt(loc, pt1, &locpt1, &crsFromPoint, tol, eps);
            err |= invDist(locpt1, center, &dist1, eps);
            error = dist1 - radius;

            geodarray[0] = geodarray[1];
            geodarray[1] = newdistbase;
            errarray[0] = errarray[1];
            errarray[1] = error;

            newdistbase = findRootSecantMethod(geodarray, errarray, &err);

            stepSize = fabs(geodarray[1] - geodarray[0]);
        } /* end while */

        if (fabs(error) >= MAX_DISTANCE_ERROR)
        {
            err |= ERROR_MAX_REACHED_ERR;
        }

        if (k >= MAX_ITERATIONS)
        {
            err |= SEC_NOT_CONVERGED_ERR;
            //return err;
        }
        else
        {
            /* Check if found intersection lies on original(not extended) locus */
            if ((orgloc.lineType == INFINITE) || (ptIsOnLocus(orgloc,
                    locpt1, &perpint, &tmpErr, tol, eps)))
            {
                intx[*n] = locpt1;
                *n = *n + 1;
            }
        }
    } /* end for i */

    return err;
}

/*******************************************************************************
 * Find the intersection of a locus and a geodesic
 */
ErrorSet locusGeoIntx(LLPoint geost, LLPoint geoend, Locus loc,
                            LLPoint* pint, double tol, double eps)
{

    ErrorSet err = 0;

    LLPoint p1;
    Locus newLoc;

    double geodarray[2], errarray[2];
    double crs12, fcrs, tcrs, crsFrompt, distFrompt;
    double crsBase, theta, faz, costheta;
    double distloc, error, distbase, newdistbase;
    double endCrsToApproxPt, endCrsBase, endDistBase;
    double angleToApproxPt, crsToApproxPt;
    int k = 0, online;
    LLPoint npPint;

    /* Assign local storage if optional pointers are not provided */
    if (NULL == pint) pint = &npPint;


//    double incr = 0.01;

    /* Find intersection of line defined by locus start and end points
     * with geodesic as a starting point for finding intersection of
     * locus and geodesic. */
//    printf("locusStart: %.15f, %.15f\n",loc.locusStart.latitude*180.0/M_PI,
//           loc.locusStart.longitude*180.0/M_PI);
//    printf("locusEnd: %.15f, %.15f\n",loc.locusEnd.latitude*180.0/M_PI,
//           loc.locusEnd.longitude*180.0/M_PI);
//    printf("geoStart: %.15f, %.15f\n",geost.latitude*180.0/M_PI,
//           geost.longitude*180.0/M_PI);
//    printf("  geoEnd: %.15f, %.15f\n",geoend.latitude*180.0/M_PI,
//           geoend.longitude*180.0/M_PI);

    /* Convert given locus to new locus with respect to given geodesic */
    inverse(geost,geoend,&crs12,NULL,&distbase,eps);
    if (distbase < tol)
    {
        /* Given geodesic is too short to compute intersection */
        err |= INVALID_SHAPE_ERR;
        return err;
    }
    err |= locusFromGeoAndPts(&newLoc,geost,crs12,loc.locusStart,loc.locusEnd,tol,eps);
    {
        double maxDist = SPHERE_RADIUS_NMI*M_PI_2;
        if (fabs(newLoc.slope) < fabs(newLoc.startDist/maxDist))
        {
            /* Locus is practically parallel to geodesic, so they won't intersect in hemisphere */
            err |= NO_INTERSECTION_ERR;
            return err;
        }
        distbase = -newLoc.startDist/newLoc.slope;
        if (fabs(distbase) < tol)
        {
            p1 = newLoc.geoStart;
        } else
        {
            /* Negative distabase handled by direct (M_PI will be added to course) */
            err |= direct(newLoc.geoStart, newLoc.geoAz, distbase, &p1, eps);
        }

        /* Place p1 precisely on geodesic */
        err |= inverse(geost,p1,&fcrs,NULL,&distbase,eps);
        if (fabs(modlon(fcrs-crs12)) < M_PI_2)
            /* Direction to p1 is same as direction of geodesic */
            fcrs = crs12;
        else
            /* Direction to p1 is opposite of geodesic direction */
            fcrs = modcrs(crs12+M_PI);

        err |= direct(geost, fcrs, distbase, &p1, eps);

    }


//    err |= geoIntx(loc.locusStart, loc.locusEnd, INFINITE, &crs31,
//            &dist13, geost, geoend, INFINITE, &crs32, &dist23, &p1, tol, eps);
    if (err)
    {
        return err;
    }

//    printf("geo intx: %.15f, %.15f\n",p1.latitude*180.0/M_PI,p1.longitude*180.0/M_PI);
    err |= inverse(p1,geost, &crsBase, NULL, NULL, eps);
    if (err)
        return err;
    // to capture crsBase = course from approx point to geoStart
    err |= inverse(geoend, p1, &endCrsToApproxPt, &endCrsBase,
            &endDistBase, eps);
    if (err)
        return err;

    /* Setup iteration to work from farthest reference point on geodesic.  This
     * makes the iteration less sensitive to azimuth noise. */
    if (distbase < endDistBase)
    {
        geost = geoend;
        distbase = endDistBase;
        fcrs = endCrsToApproxPt;
        crsBase = endCrsBase;
    }

    /* Find course and length of locus's geodesic */
    tcrs = loc.geoAz;
    /* Find perp intercept of line-geodesic intersection point and the
     * geodesic of the locus. */
    err |= projectToGeo(loc.geoStart, tcrs, p1, pint, &crsFrompt,
            &distFrompt, tol, eps);
    if (err)
        return err;

    distloc = distToLocusFromGeoPt(loc, *pint, &faz, &err, tol, eps);
    if (err & RADIUS_OUT_OF_RANGE_ERR)
    {
        /* In this case, approx intersection is so far from defining geodesic that distance
         * computations and other operations will not be reliable.  We will return a flag
         * that indicates no intersection exists in this hemisphere. */
        /* Mask off RADIUS err */
        err &= ~RADIUS_OUT_OF_RANGE_ERR;
        /* Add error flag indicating no intersection */
        err |= NO_INTERSECTION_ERR;
        /* return */
        return err;
    }
    if (distloc < 0.0)
        distloc = -distloc;

    error = distFrompt - distloc;
    errarray[1] = error;
    geodarray[1] = distbase;

    /* Initial guess is updated as described in algorithm document */
    theta = fabs(modlon(crsFrompt - crsBase));
    costheta = cos(theta);
    if (fabs(costheta) > 0.001)
      newdistbase = distbase - errarray[1] / costheta;
    else
      newdistbase = distbase - 0.5 * error; 

    /* We iterate by varying a point (pt1) along the geodesic until that
     * point also lies on the locus. */
    /* NOTE: 3.0*MAX_DISTANCE_ERROR is used because for some parallel configurations,
     * the initial guess is on the wrong end of the locus.  The iteration has been found to
     * converge after a large jump in these cases. */
    //TODO: Solve this problem by detecting when the initial solution is on the wrong side
    while ((k++ < MAX_ITERATIONS) && (fabs(error) > tol))
    {
        err |= direct(geost, fcrs, newdistbase, &p1, eps);
        if (fabs(newdistbase) > 10740.0) //represents 179 degree angle
          err |= NO_INTERSECTION_ERR;
        if (err)
            return err;

        /* Determine which side of loc.geodesic approx pt lies */
        err |= inverse(loc.geoStart,p1,&crsToApproxPt,NULL,NULL,eps);
        angleToApproxPt = modlon(crsToApproxPt - loc.geoAz);

        err |= projectToGeo(loc.geoStart, tcrs, p1, pint, &crsFrompt,
                &distFrompt, tol, eps);
        if (err)
        {
            err |= UNEXPECTED_ERR;
            return err;
        }

        /* Change sign of distance to match side of loc.geodesic point is on */
        if (angleToApproxPt < 0) distFrompt = -distFrompt;

        /* Compute expected distance to locus */
        distloc = distToLocusFromGeoPt(loc, *pint, &faz, &err, tol, eps);
        if (err & RADIUS_OUT_OF_RANGE_ERR)
        {
            /* Mask off RADIUS err */
            err &= ~RADIUS_OUT_OF_RANGE_ERR;
            /* Add error flag indicating no intersection */
            err |= NO_INTERSECTION_ERR;
            /* return */
            return err;
        }
        if (err)
            return err;

        error = distFrompt - distloc;
        geodarray[0] = geodarray[1];
        geodarray[1] = newdistbase;
        errarray[0] = errarray[1];
        errarray[1] = error;
//        printf("%d\t%.15e\t%.15e\t%.15e\n",k,newdistbase,error,distFrompt);
        distbase = newdistbase;
        newdistbase = findRootSecantMethod(geodarray, errarray, &err);

//        if (fabs(newdistbase-distbase) > fabs(distbase))
//        {
//            newdistbase = distbase + sgn(newdistbase-distbase)*distbase;
//        }
    } /* end while */
    if (k >= MAX_ITERATIONS)
    {
        err |= SEC_NOT_CONVERGED_ERR;
    }
    if (fabs(error) >= MAX_DISTANCE_ERROR)
    {
        err |= ERROR_MAX_REACHED_ERR;
    }

    /* MJM: I don't think we should check this until the exact point is found */
    /* Check that intersection point lies on geodesic.  Return null if not */
    online = ptIsOnGeo(loc.geoStart, loc.geoEnd, *pint, loc.lineType,
            &err, tol, eps);
    if (!online)
        err |= NO_INTERSECTION_ERR;

    *pint = p1;

    return err;

}

/*******************************************************************************
 * Find the intersection of two loci
 */
ErrorSet locusIntx(Locus loc1, Locus loc2, LLPoint* intx, double tol,
                         double eps)
{

    ErrorSet err = 0;

    Locus orgLoc1 = loc1, orgLoc2 = loc2;
    LLPoint p1;
    LLPoint pint1;
    LLPoint pint2;
    LLPoint ploc1, ploc2, ploc1array[2];
    //    LLPoint pt1;
    double geodarray[2], errarray[2];
    double crs31, dist13, crs32, dist23, crsFrompt, distFrompt;
    double error, incr, distbase;
    double testCrs, testRecipCrs, loc2Dist, angle;

    int k = 0, online;

    Shape common;
    LLPoint commonPt;

    LLPoint npIntx;

    /* Assign local storage if optional pointers are not provided */
    if (NULL == intx) intx = &npIntx;


    /*
     * check if locus 1 and locus 2 coincide.
     *
     * This section of code will address the cases when the
     * geodesi are collinear for the given loci and there
     * exists an intersection point or a common locus.
     *
     * The case when no intersection exists and the geodesi
     * are collinear is handled by the original locus intersect
     * algorithm.  See bug 18896 for additional details.
     */
    if (lociCoincide(loc1, loc2, &common, &err, tol, eps))
    {
        /*
         * check the shape type of the common shape.
         * If it's a point then is the intersection of the loci.
         * If it's a locus then error out since locus intersect only returns one intersection point.
         */
        switch (common.type) {
        case LLPOINT:
            commonPt = *((LLPoint*) common.this_shape);
            err |= createPt(intx, commonPt.latitude, commonPt.longitude);
            return err;
        case LOCUS:
            /* Should this return some sort of LOCI_COINCIDE_ERR??
             * Don't forget to take into account error code masking in the JNI.
             * It's special for the locus intersect function.  - jamezcua
             */
            err |= SHAPE_NOT_DEFINED_ERR;//TODO  Need to update error code so it's more informative.
            return err;
        default:
            //unexpected shape returned, error out
            err |= INVALID_TYPE_ERR;
            return err;
        }
    }

    incr = 1.0001;
    /* Find intersection of geodesics defined by locus start and end points
     * (p1) as a starting point for finding intersection of two loci  */
    err |= geoIntx(loc1.locusStart, loc1.locusEnd, INFINITE, &crs31,
            &dist13, loc2.locusStart, loc2.locusEnd, INFINITE, &crs32, &dist23,
            &p1, tol, eps);

    if (err)
        return err;

    if (dist13 < 10.0 / 1852.0) /* Ten meters in NM */
    {
        /* Intersection point is within ten meters of locus1's start point.
         * This will cause problems for iteration scheme because azimuths lose
         * precision over short distances.
         * Move start point of Locus1 back 1 NM */
        {
            LLPoint newGeoStart;
            double newStartDist;
            err |= direct(loc1.geoStart, loc1.geoAz + M_PI, 1.0,
                    &newGeoStart, eps);
            newStartDist = distToLocusFromGeoDist(loc1, -1.0);
            err |= createLocus(&loc1, newGeoStart, loc1.geoEnd, newStartDist,
                    loc1.endDist, loc1.lineType, tol, eps);
        }
    }

    if (dist23 < 10.0 / 1852.0) /* Ten meters in NM */
    {
        /* Intersection point is within ten meters of locus2's start point.
         * This will cause problems for iteration scheme because azimuths lose
         * precision over short distances.
         * Move start point of Locus2 back 1 NM */
        {
            LLPoint newGeoStart;
            double newStartDist;
            err |= direct(loc2.geoStart, loc2.geoAz + M_PI, 1.0,
                    &newGeoStart, eps);
            newStartDist = distToLocusFromGeoDist(loc2, -1.0);
            err |= createLocus(&loc2, newGeoStart, loc2.geoEnd, newStartDist,
                    loc2.endDist, loc2.lineType, tol, eps);
        }
    }

    //    printf("intx lat, lon = %.15lf, %.15lf\n", p1.latitude * 180.0 / M_PI,
    //            p1.longitude * 180.0 / M_PI);

    /* Find perp intercept (pint1) of intersection of lines and the
     * geodesic of locus1. */
    err |= projectToGeo(loc1.geoStart, loc1.geoAz, p1, &pint1,
            &crsFrompt, &distFrompt, tol, eps);

    if (err)
        return err;

    online = ptIsOnGeo(loc1.geoStart, loc1.geoEnd, pint1, loc1.lineType,
            &err, tol, eps);
    /* With Mike's approval, the following test was taken out on 7-22-11 */
    /*if (!online)
    {
        //TODO MJM believes this check is premature--may lead to missed intersections
        err |= NO_INTERSECTION_ERR;
        return err;
    }*/
    err |= invDist(loc1.geoStart, pint1, &distbase, eps);
    if (err)
        return err;

    /* We iterate by varying a point (pint1) along the geodesic of locus1
     * until the corresponding point on locus1 (ploc1) is also on
     * locus2. */
    /* Initialize the error function arrays */
    /* These values will never be used      */
    geodarray[1] = 0.0;
    errarray[1] = 1.0e99;
    ploc1array[1].latitude = 1.0e99;
    ploc1array[1].longitude = 1.0e99;
    /* Don't stop until intersection point is within tol of loc2, AND
     * intersection point changes by less than tol
     */
    while ((k == 0) || ((k < MAX_ITERATIONS) &&
             ((fabs(error) > tol) || (fabs(geodarray[1] - geodarray[0]) > tol))))
    {

        //        printf("++ Loc Int: k = %d\n", k);

        if (k > 0)
        {
            /* Updated position on loc1's geodesic */
            err |= direct(loc1.geoStart, loc1.geoAz, distbase, &pint1, eps);
        }
        /* Point on loc1 corresponding to pint1 */
        err |= ptOnLocusFromGeoPt(loc1, pint1, &ploc1, &crsFrompt, tol, eps);

        err |= projectToGeo(loc2.geoStart, loc2.geoAz, ploc1, &pint2,
                &crsFrompt, &distFrompt, tol, eps);
        loc2Dist = distToLocusFromGeoPt(loc2, pint2, NULL, &err, tol, eps);
        /* On which side of Locus2's geodesic is ploc1? */
        err |= invCrs(loc2.geoStart, ploc1, &testCrs, &testRecipCrs, eps);
        angle = modlon(loc2.geoAz - testCrs);
        if (angle > 0)
        {
            /* ploc1 (estimated intersection) is to left of Locus2's geodesic */
            distFrompt = -distFrompt;
        }

        //        err |= ptOnLocusFromGeoPt(loc2, pint2, &ploc2, &crsFrompt, tol, eps);
        //        err |= invDist(ploc1, ploc2, &error, eps);
        error = distFrompt - loc2Dist;
        //        printf("k = %d, error = %.15e, ", k, error);
        ptIsOnLocus(loc2, ploc1, &p1, &err, tol, eps);

        //        printf("k = %d error = %20.16f\n",k,error);
        geodarray[0] = geodarray[1];
        geodarray[1] = distbase;
        errarray[0] = errarray[1];
        errarray[1] = error;
        if (k == 0)
        {
            /* For this first step, we project ploc2 back to loc1's
             * geodesic to get a new estimate of distbase */
            err |= ptOnLocusFromGeoPt(loc2, pint2, &ploc2, &crsFrompt, tol, eps);
            err |= projectToGeo(loc1.geoStart, loc1.geoAz, ploc2, &pint1,
                    &crsFrompt, &distFrompt, tol, eps);
            err |= invDist(loc1.geoStart, pint1, &distbase, eps);
        }
        else
        {
            distbase = findRootSecantMethod(geodarray, errarray, &err);
        }

        /* If ploc1 stops changing shift distbase just a little
         * to get the algorithm moving again.
         */
        ploc1array[0] = ploc1array[1];
        ploc1array[1] = ploc1;
        if (ploc1array[0].latitude == ploc1array[1].latitude &&
        		ploc1array[0].longitude == ploc1array[1].longitude){
        	distbase += tol;
        }

        k++;

    } /* end while */

    if (fabs(error) >= MAX_DISTANCE_ERROR)
    {
        err |= ERROR_MAX_REACHED_ERR;
    }

    if (k >= MAX_ITERATIONS)
    {
        err |= SEC_NOT_CONVERGED_ERR;
        return err;
    }
    //    printf("Inf intx found: ");
    //    _display(ploc1);

    if (ptIsOnLocus(orgLoc1, ploc1, &p1, &err, tol, eps)
            && ptIsOnLocus(orgLoc2, ploc1, &p1, &err, tol, eps))
    {
        *intx = ploc1;
    }
    else
    {
        err |= NO_INTERSECTION_ERR;
    }

    return err;
}

static ErrorSet initSpiralGeoIntx(Spiral sp, Geodesic line, LLPointSet* pts, double tol, double eps) {

	double perpAz, perpDist, rad, lineDist, lineCrs, az12, az21, dist, crs31, dist31, crs32, dist32, testAz;
	double crs, angle;
	LLPoint perpPt, midChord, linePt, intPt, spPt;
	ErrorSet err = 0;
	int i = 0, index = 0;

	double distarray[2] = { 0.0, 0.0 };
	double errarray[2] = { 0.0, 0.0 };

	err |= spiralMidChord(sp, line, &midChord, tol, eps);

	err = 0;

	//Tangent Intersection
	if (ptIsOnGeo(line.startPoint, line.endPoint, midChord, line.lineType, &err, tol, eps)) {
		err |= addPtToPtSet(pts, &midChord);
		return err;
	}

	err |= projectToGeo(line.startPoint, line.startAz, sp.centerPoint, &perpPt, &perpAz, &perpDist, tol, eps);
	err |= createGeo(&line, perpPt, line.startPoint, INFINITE, eps);

	testAz = perpAz + M_PI / 2;

	if (calculateSubtendedAngle(sp.startAz, testAz, sp.dir) > sp.subtendedAngle) {
		testAz = sp.startAz;
	}
	err |= ptOnSpiral(sp, testAz, &spPt, eps);
	err |= projectToGeo(line.startPoint, line.startAz, spPt, &linePt, &az12, &dist, tol, eps);
	err |= inverse(line.startPoint, linePt, &lineCrs, NULL, &lineDist, eps);

	if (fabs(lineDist) < .0001) {
		err |= direct(line.startPoint, line.startAz + M_PI, 10, &linePt, eps);
		err |= inverse(linePt, line.startPoint, &az12, &az21, &lineDist, eps);
		line.startPoint = linePt;
		line.startAz = az12;
		lineCrs = az12;
	}

	distarray[0] = lineDist;

	err |= inverse(sp.centerPoint, linePt, &az12, &az21, &dist, eps);
	err |= spiralRadius(sp, az12, &rad);
	errarray[0] = rad - dist;

	lineDist = lineDist * 1.01;
	err |= direct(line.startPoint, lineCrs, lineDist,  &linePt, eps);
	distarray[1] = lineDist;
	err |= inverse(sp.centerPoint, linePt, &az12, &az21, &dist, eps);
	err |= spiralRadius(sp, az12, &rad);
	errarray[1] = rad - dist;

	while (1) {

		lineDist = findRootSecantMethod(distarray, errarray, &err);

		if (fabs(lineDist) < .0001) {
			err |= direct(line.startPoint, line.startAz + M_PI, 10, &linePt, eps);
			err |= inverse(linePt, line.startPoint, &az12, &az21, &lineDist, eps);
			line.startPoint = linePt;
			line.startAz = az12;
			lineCrs = az12;
		}

		err |= direct(line.startPoint, lineCrs, lineDist, &linePt, eps);
		err |= inverse(sp.centerPoint, linePt, &az12, &az21, &dist, eps);
		spiralRadius(sp, az12, &rad);

		distarray[0] = distarray[1];
		errarray[0] = errarray[1];

		distarray[1] = lineDist;
		errarray[1] = rad - dist;

		if ((fabs(rad - dist) < tol) && (fabs(distarray[1] - distarray[0]) < tol)) {
			err |= addPtToPtSet(pts, &linePt);
			i++;

			if (i > 1) {
				return err;
			}

			testAz = perpAz - M_PI / 2;
			if (calculateSubtendedAngle(sp.startAz, testAz, sp.dir) > sp.subtendedAngle) {
				testAz = sp.startAz;
			}

			err |= ptOnSpiral(sp, testAz, &spPt, eps);
			err |= spiralRadius(sp, testAz, &rad);

			err |= projectToGeo(line.startPoint, line.startAz, spPt, &linePt, &az12, &dist, tol, eps);
			err |= inverse(line.startPoint, linePt, &lineCrs, NULL, &lineDist, eps);

			if (fabs(lineDist) < .0001) {
				err |= direct(line.startPoint, line.startAz, -50, &linePt, eps);
				err |= inverse(linePt, line.startPoint, &az12, &az21, &lineDist, eps);
				line.startPoint = linePt;
				line.startAz = az12;
				lineCrs = az12;
			}

			distarray[0] = lineDist;

			err |= inverse(sp.centerPoint, linePt, &az12, &az21, &dist, eps);
			err |= spiralRadius(sp, az12, &rad);
			errarray[0] = rad - dist;

			lineDist = lineDist * 1.01;
			err |= direct(line.startPoint, lineCrs, lineDist,  &linePt, eps);
			distarray[1] = lineDist;
			err |= inverse(sp.centerPoint, linePt, &az12, &az21, &dist, eps);
			err |= spiralRadius(sp, az12, &rad);
			errarray[1] = rad - dist;

			index = 0;
		}
		if (index > 100) {
			err |= ITERATION_MAX_REACHED_ERR;
			return err;
		}
		index += 1;
	}


//	err |= spiralRadius(sp, line.startAz, &rad);
//
//	err |= invCrs(line.startPoint, perpPt, &az12, &lineAz, eps);
//
//	//This value is based on a circular approximation.  If the distances are too similar, the value is not useful
//	//The approximation becomes less accurate the closer to the center of the spiral
//	if (perpDist < rad * .9) {
//		dtheta = asin(perpDist / rad);
//	} else {
//		dtheta = 1;
//	}
//	dr = sp.growthRate * dtheta;
//
//	if (azDifference(lineAz, perpAz) > 0) {
//		lineDir = -1;
//	} else {
//		lineDir = 1;
//	}
//
//	printf("\n\n\n");
//
//	lineDist = cos(dtheta) * rad - lineDir * sp.dir * dr;
//	distarray[0] = lineDist;
//	err |= direct(perpPt, lineAz, lineDist, &linePt, eps);
//	err |= inverse(sp.centerPoint, linePt, &az12, &az21, &dist, eps);
//	if (calculateSubtendedAngle(sp.startAz, az12, sp.dir) > M_PI * 3 / 2) {
//		printf("Point on wrong section of spiral (start rad = 0)\n");
//		//Initial test point on the wrong section of the spiral (Currently occurs for spirals with 0 start radius)
//		err |= crsIntx(sp.centerPoint, sp.startAz, &crs31, &dist31, line.startPoint, line.startAz, &crs32, &lineDist, &linePt, tol, eps);
//		err |= inverse(sp.centerPoint, linePt, &az12, &az21, &dist, eps);
//		distarray[0] = lineDist;
//	} else if (calculateSubtendedAngle(sp.startAz, az12, sp.dir) < M_PI / 2){
//		printf("Point on wrong section of spiral (end rad = 0)\n");
//		//Initial test point on the wrong section of the spiral (Might occur for spirals with 0 end radius, but not seen yet)
//		err |= crsIntx(sp.centerPoint, sp.endAz, &crs31, &dist31, line.startPoint, line.startAz, &crs32, &lineDist, &linePt, tol, eps);
//		err |= inverse(sp.centerPoint, linePt, &az12, &az21, &dist, eps);
//		distarray[0] = lineDist;
//	}
//	err |= spiralRadius(sp, az12, &rad);
//	errArray[0] = (rad - dist);
//
//	lineDist = lineDist * 1.01;
//	distarray[1] = lineDist;
//	err |= direct(perpPt, lineAz, lineDist, &linePt, eps);
//	err |= inverse(sp.centerPoint, linePt, &az12, &az21, &dist, eps);
//	err |= spiralRadius(sp, az12, &rad);
//	errArray[1] = rad - dist;
//
////	printf("Initial Az:  %f\n", az12);
////	printf("Subtended Angle: %f\n", calculateSubtendedAngle(sp.startAz, az12, sp.dir));
////
////	printf("\nBeginning Loop\n\n");
//
//	while (1) {
//
//		lineDist = findRootSecantMethod(distarray, errArray, &err);
//
//		err |= direct(perpPt, lineAz, lineDist, &linePt, eps);
//		err |= inverse(sp.centerPoint, linePt, &az12, &az21, &dist, eps);
//		spiralRadius(sp, az12, &rad);
//
////		printf("Azimuth:  %f\n", az12);
////
////		printf("Subtended Angle: %f\n", calculateSubtendedAngle(sp.startAz, az12, sp.dir));
//
//		distarray[0] = distarray[1];
//		errArray[0] = errArray[1];
//
//		distarray[1] = lineDist;
//		errArray[1] = rad - dist;
//
//		if ((fabs(rad - dist) < tol) && (fabs(distarray[1] - distarray[0]) < tol)) {
//			err |= addPtToPtSet(pts, &linePt);
//			i++;
//
//			if (i > 1) {
//				return err;
//			}
//
////			printf("\nFound Point\n\n");
//
//			lineAz = fmod(lineAz + M_PI, 2*M_PI);
//			err |=spiralRadius(sp, fmod(line.startAz + M_PI, 2*M_PI), &rad);
//			lineDir = -1*lineDir;
//			lineDist = cos(dtheta) * rad - lineDir * sp.dir * dr;
//			distarray[0] = lineDist;
//			err |= direct(perpPt, lineAz, lineDist, &linePt, eps);
//			err |= inverse(sp.centerPoint, linePt, &az12, &az21, &dist, eps);
//			if (calculateSubtendedAngle(sp.startAz, az12, sp.dir) > M_PI * 3 / 2) {
//				//Initial test point on the wrong section of the spiral (Currently occurs for spirals with 0 start radius)
//				err |= crsIntx(sp.centerPoint, sp.startAz, &crs31, &dist31, line.startPoint, line.startAz, &crs32, &lineDist, &linePt, tol, eps);
//				err |= inverse(sp.centerPoint, linePt, &az12, &az21, &dist, eps);
//				distarray[0] = lineDist;
//			} else if (calculateSubtendedAngle(sp.startAz, az12, sp.dir) < M_PI / 2){
//				//Initial test point on the wrong section of the spiral (Might occur for spirals with 0 end radius, but not seen yet)
//				err |= crsIntx(sp.centerPoint, sp.endAz, &crs31, &dist31, line.startPoint, line.startAz, &crs32, &lineDist, &linePt, tol, eps);
//				err |= inverse(sp.centerPoint, linePt, &az12, &az21, &dist, eps);
//				distarray[0] = lineDist;
//			}
//			err |= spiralRadius(sp, az12, &rad);
//			errArray[0] = (rad - dist);
//
////			printf("Initial Az:  %f\n", az12);
////			printf("Subtended Angle: %f\n", calculateSubtendedAngle(sp.startAz, az12, sp.dir));
//
//			lineDist = lineDist * 1.01;
//			distarray[1] = lineDist;
//			err |= direct(perpPt, lineAz, lineDist, &linePt, eps);
//			err |= inverse(sp.centerPoint, linePt, &az12, &az21, &dist, eps);
//			err |= spiralRadius(sp, az12, &rad);
//			errArray[1] = rad - dist;
//
//			index = 0;
//		}
//		if (index > 100) {
//			err |= ITERATION_MAX_REACHED_ERR;
//			return err;
//		}
//		index += 1;
//	}
}

ErrorSet spiralGeoIntx(Spiral sp, Geodesic line, LLPointSet* pts, double tol, double eps) {

LLPoint tempPt;
Spiral tempSp;
LLPointSet testPts = createPtSet();
double initAz, dist, rad, az21;
int i;
LineType lineType;
ErrorSet err = 0;
ErrorSet temperr1 = 0;
ErrorSet temperr2 = 0;

if (ptIsOnSpiral(sp, line.startPoint, tol, eps)) {
	err |= addPtToPtSet(&testPts, &line.startPoint);
}
if (ptIsOnSpiral(sp, line.endPoint, tol, eps)) {
	err |= addPtToPtSet(&testPts, &line.endPoint);
}
if (ptIsOnGeo(line.startPoint, line.endPoint, sp.startPoint, LineType::INFINITE, &err, tol, eps)) {
	err |= addPtToPtSet(&testPts, &sp.startPoint);
}
if (ptIsOnGeo(line.startPoint, line.endPoint, sp.endPoint, LineType::INFINITE, &err, tol, eps)) {
	err |= addPtToPtSet(&testPts, &sp.endPoint);
}

err |= projectToGeo(line.startPoint, line.startAz, sp.centerPoint, &tempPt, &initAz, &dist, tol, eps);

if (dist == 0) {

	lineType = line.lineType;
	line.lineType = INFINITE;
	initAz = geoCrs(line, sp.centerPoint, NULL, NULL, NULL, &err, tol, eps);
	ptOnSpiral(sp, initAz, &tempPt, eps);
	addPtToPtSet(&testPts, &tempPt);

	line.lineType = lineType;

	initAz += M_PI;
	ptOnSpiral(sp, initAz, &tempPt, eps);
	addPtToPtSet(&testPts, &tempPt);

	for (i=0; i<testPts.length; i++) {
		tempPt.latitude = testPts.elements[i]->latitude;
		tempPt.longitude = testPts.elements[i]->longitude;
		if (ptIsOnSpiral(sp, tempPt, tol, eps)) {
			if (ptIsInSet(tempPt, *pts, tol) == 0) {
				err |= addPtToPtSet(pts, &tempPt);
			}
		}
	}

	return err;
}

err = 0;

err |= spiralRadius(sp, initAz, &rad);
if (fabs(rad - sp.startRadius) > sp.growthRate * M_PI) {
	//initAz is on the wrong spiral section (behind the start point)
	//recalculate rad
   ArcDirection oppositeDir(ArcDirection::CLOCKWISE);
   if (sp.dir == ArcDirection::CLOCKWISE) oppositeDir = ArcDirection::COUNTERCLOCKWISE;
	err |= createSpiral(&sp, sp.centerPoint, sp.endRadius, sp.startRadius, sp.endAz, sp.startAz, oppositeDir, eps);
	err |= spiralRadius(sp, initAz, &rad);
}
temperr1 |= createSpiralSection(sp, initAz, rad, &tempSp, eps);
//displayMatlabSpiral(tempSp,"tempSp1",0);

if (temperr1 == 0) {
	temperr1 |= initSpiralGeoIntx(tempSp, line, &testPts, tol, eps);
}

if ((tempSp.startRadius - sp.startRadius) < 0) {
	rad = rad + fabs(sp.growthRate*2*M_PI);
} else {
	rad = rad - fabs(sp.growthRate*2*M_PI);
}

temperr2 |= createSpiralSection(sp, initAz, rad, &tempSp, eps);
//displayMatlabSpiral(tempSp,"tempSp2",0);

if (temperr2 == 0) {
	temperr2 |= initSpiralGeoIntx(tempSp, line, &testPts, tol, eps);
}

if (((temperr1) && (temperr2)) && (testPts.length == 0)) {
	err |= NO_INTERSECTION_ERR;
	return err;
}

for (i=0; i<testPts.length; i++) {
	tempPt.latitude = testPts.elements[i]->latitude;
	tempPt.longitude = testPts.elements[i]->longitude;
//	displayMatlabPt(tempPt,"tmpt",0);
	if (ptIsOnSpiral(sp, tempPt, tol, eps)) {
		if (ptIsInSet(tempPt, *pts, tol) == 0) {
			err |= addPtToPtSet(pts, &tempPt);
		}
	}
}

return err;
}

static ErrorSet initSpiralLocusIntx(Spiral sp, Locus locus, LLPointSet* pts, double tol, double eps) {

double  crs31, dist31, crs32, dist32, az12, az21, dist, rad, geoAz, lineDist, perpCrs;
int index = 0, count = 0;
Geodesic geo;
LLPoint midChord, testPt1, testPt2, geoPt1, geoPt2, intPt, locPt, perpPt;
LLPointSet testPts = createPtSet();
ErrorSet err = 0;

double distArray[2] = { 0.0, 0.0 };
double errArray[2] = { 0.0, 0.0 };

if (ptIsOnLocus(locus, sp.centerPoint, &geoPt1, &err, tol, eps)) {
	az12 = locusCrsAtPt(locus, sp.centerPoint, &geoPt1, &perpCrs, &err, tol, eps);
	err |= ptOnSpiral(sp, az12, &testPt1, eps);
	err |= ptOnSpiral(sp, az12 + M_PI, &testPt2, eps);
} else {
	err |= createGeo(&geo, locus.locusStart, locus.locusEnd, INFINITE, eps);

	err |= spiralMidChord(sp, geo, &midChord, tol, eps);

	if (ptIsOnLocus(locus, midChord, &geoPt1, &err, tol, eps)) {
		err |= addPtToPtSet(pts, &midChord);
		return err;
	}

	if (geoIntx(geo.startPoint, geo.endPoint, INFINITE, &crs31, &dist31, sp.centerPoint, midChord, SEGMENT, &crs32, &dist32, &intPt, tol, eps)) {
		return err;
	}

	err |= initSpiralGeoIntx(sp, geo, &testPts, tol, eps);

	if(testPts.length == 2){
		testPt1.latitude = testPts.elements[0]->latitude;
		testPt1.longitude = testPts.elements[0]->longitude;

		testPt2.latitude = testPts.elements[1]->latitude;
		testPt2.longitude = testPts.elements[1]->longitude;
	} else {
		return err;
	}
}

if (ptIsOnLocus(locus, testPt1, &geoPt1, &err, tol, eps) && ptIsOnLocus(locus, testPt2, &geoPt2, &err, tol, eps)) {
	err |= addPtToPtSet(pts, &testPt1);
	err |= addPtToPtSet(pts, &testPt2);
	return err;
}
err |= projectToLocus(locus, testPt2, &locPt, &az12, &dist, tol, eps);
err |= projectToGeo(locus.geoStart, locus.geoAz, locPt, &geoPt2, &az12, &dist, tol, eps);

err |= projectToLocus(locus, testPt1, &locPt, &az12, &dist, tol, eps);
err |= projectToGeo(locus.geoStart, locus.geoAz, locPt, &geoPt1, &az12, &dist, tol, eps);

err |= projectToGeo(locus.geoStart, locus.geoAz, sp.centerPoint, &perpPt, &az12, &dist, tol, eps);

err |= inverse(perpPt, geoPt1, &geoAz, &az21, &lineDist, eps);
err |= inverse(sp.centerPoint, locPt, &az12, &az21, &dist, eps);
err |= spiralRadius(sp, az12, &rad);

distArray[0] = lineDist;
errArray[0] = dist - rad;

lineDist = lineDist * 1.01;
err |= direct(perpPt, geoAz, lineDist, &geoPt1, eps);
err |= ptOnLocusFromGeoPt(locus, geoPt1, &locPt, &perpCrs, tol, eps);
err |= inverse(sp.centerPoint, locPt, &az12, &az21, &dist, eps);
err |= spiralRadius(sp, az12, &rad);

distArray[1] = lineDist;
errArray[1] = dist - rad;

lineDist = findRootSecantMethod(distArray, errArray, &err);

while (1){

	lineDist = findRootSecantMethod(distArray, errArray, &err);
	err |= direct(perpPt, geoAz, lineDist, &geoPt1, eps);
	err |= ptOnLocusFromGeoPt(locus, geoPt1, &locPt, &perpCrs, tol, eps);
	err |= inverse(sp.centerPoint, locPt, &az12, &az21, &dist, eps);
	err |= spiralRadius(sp, az12, &rad);

	distArray[0] = distArray[1];
	errArray[0] = errArray[1];

	distArray[1] = lineDist;
	errArray[1] = dist - rad;

	if ((fabs(dist - rad) < tol) && (fabs(distArray[1] - distArray[0]) < tol)) {
		if (count == 1) {
			err |= addPtToPtSet(pts, &locPt);
			return err;
		}
		err |= addPtToPtSet(pts, &locPt);
		count++;
		err |= inverse(perpPt, geoPt2, &geoAz, &az21, &lineDist, eps);
		err |= inverse(sp.centerPoint, testPt2, &az12, &az21, &dist, eps);
		err |= spiralRadius(sp, az12, &rad);

		distArray[0] = lineDist;
		errArray[0] = dist - rad;

		lineDist = lineDist * 1.01;
		err |= direct(perpPt, geoAz, lineDist, &geoPt1, eps);
		err |= ptOnLocusFromGeoPt(locus, geoPt1, &locPt, &perpCrs, tol, eps);
		err |= inverse(sp.centerPoint, locPt, &az12, &az21, &dist, eps);
		err |= spiralRadius(sp, az12, &rad);

		distArray[1] = lineDist;
		errArray[1] = dist - rad;
	}

	index++;
	if (index == 100) {
		err |= ITERATION_MAX_REACHED_ERR;
		return err;
	}
}

}

//Finds the points where the input spiral and geodesic intersect.  This function calls findGeoIntx to determine
//all candidate points, and then checks to determine which of them lie on the input spiral.
ErrorSet spiralLocusIntx(Spiral sp, Locus locus, LLPointSet* pts, double tol, double eps)
{
LLPoint tempPt;
Spiral tempSp;
LLPointSet testPts = createPtSet();
double initAz, tempDist, rad;
int i;
ErrorSet err = 0, temperr1 = 0, temperr2 = 0;
locus.lineType = INFINITE;

if (ptIsOnSpiral(sp, locus.locusStart, tol, eps)) {
	err |= addPtToPtSet(&testPts, &locus.locusEnd);
}
if (ptIsOnSpiral(sp, locus.locusEnd, tol, eps)) {
	err |= addPtToPtSet(&testPts, &locus.locusEnd);
}
/*if (ptIsOnLocus(locus, sp.startPoint)) {
	err |= addPtToPtSet(pts, sp.startPoint);
}
if (ptIsOnLocus(locus, sp.endPoint)) {
	err |= addPtToPtSet(pts, sp.endPoint);
}
*/
err |= projectToLocus(locus, sp.centerPoint, &tempPt, &initAz, &tempDist, tol, eps);
err |= spiralRadius(sp, initAz, &rad);
err |= createSpiralSection(sp, initAz, rad, &tempSp, eps);
temperr1 |= initSpiralLocusIntx(tempSp, locus, &testPts, tol, eps);
if ((tempSp.startRadius - sp.startRadius) < 0) {
	rad = rad + fabs(sp.growthRate*2*M_PI);
} else {
	rad = rad - fabs(sp.growthRate*2*M_PI);
}

err |= createSpiralSection(sp, initAz, rad, &tempSp, eps);
temperr2 |= initSpiralLocusIntx(tempSp, locus, &testPts, tol, eps);

if ((temperr1) && (temperr2)) {
	err |= NO_INTERSECTION_ERR;
	return err;
}

for (i=0; i<testPts.length; i++) {
	tempPt.latitude = testPts.elements[i]->latitude;
	tempPt.longitude = testPts.elements[i]->longitude;
	if (ptIsOnSpiral(sp, tempPt, tol, eps)) {
		if (ptIsInSet(tempPt, *pts, tol) == 0) {
			err |= addPtToPtSet(pts, &tempPt);
		}
	}
}
return err;
}

ErrorSet initSpiralArcIntx(Spiral sp, Arc arc, double initAz, double span, LLPointSet* pts, double tol, double eps) {

	double az12, az21, dist, rad, testAz, diff, azDiff, spCentCrs, arcCentCrs;
	LLPoint spPt;
	int index = 0, count = 0;
	ErrorSet err = 0;

	double azArray[2] = { 0.0, 0.0 };
	double errArray[2] = { 0.0, 0.0 };

	err |= inverse(sp.centerPoint, arc.centerPoint, &spCentCrs, &arcCentCrs, &dist, eps);

	//Arc is far from spiral, no intersection occurs
	if ((dist > arc.radius + sp.startRadius + tol) && (dist > arc.radius + sp.endRadius + tol)) {
		return err;
	}
	//Arc is fully inside spiral
	if ((dist + arc.radius + tol < sp.startRadius) && (dist + arc.radius + tol < sp.endRadius)) {
		return err;
	}
	//Spiral fully inside arc
	if ((dist + sp.startRadius + tol < arc.radius) && (dist + sp.endRadius + tol < arc.radius)) {
		return err;
	}

	testAz = initAz + (span / 2);
	azArray[0] = testAz;

	err |= ptOnSpiral(sp, testAz, &spPt, eps);
	err |= inverse(spPt, arc.centerPoint, &az12, &az21, &dist, eps);
	diff = dist - arc.radius;
	errArray[0] = diff;

	testAz = testAz + fabs(diff) / diff * M_PI / 180;
	azArray[1] = testAz;

	err |= ptOnSpiral(sp, testAz, &spPt, eps);
	err |= inverse(spPt, arc.centerPoint, &az12, &az21, &dist, eps);
	diff = dist - arc.radius;
	errArray[1] = diff;

	while (1) {

		testAz = fmod(findRootSecantMethod(azArray, errArray, &err)+2*M_PI, 2*M_PI);

		err |= ptOnSpiral(sp, testAz, &spPt, eps);
		err |= inverse(spPt, arc.centerPoint, &az12, &az21, &dist, eps);
		diff = dist - arc.radius;

		azArray[0] = azArray[1];
		errArray[0] = errArray[1];

		azArray[1] = testAz;
		errArray[1] = diff;

		err |= spiralRadius(sp, testAz, &rad);

		if (fabs(diff) < tol) {
			if (count) {
				err |= addPtToPtSet(pts, &spPt);
				return err;
			}

			err |= addPtToPtSet(pts, &spPt);
			count++;

			err |= invCrs(sp.centerPoint, arc.centerPoint, &az12, &az21, eps);
			azDiff = azDifference(az12, testAz);

			testAz = fmod(az12 - azDiff, 2*M_PI);

			azArray[0] = testAz;

			err |= ptOnSpiral(sp, testAz, &spPt, eps);
			err |= inverse(spPt, arc.centerPoint, &az12, &az21, &dist, eps);
			diff = dist - arc.radius;
			errArray[0] = diff;

			testAz = testAz + fabs(diff) / diff * azDiff / 10;
			azArray[1] = testAz;

			err |= ptOnSpiral(sp, testAz, &spPt, eps);
			err |= inverse(spPt, arc.centerPoint, &az12, &az21, &dist, eps);
			diff = dist - arc.radius;
			errArray[1] = diff;

			index = 0;
		}

		index++;

		if (index > 50) {
			printf("Iteration Max Reached:  %e\n", diff);
			err |= NO_INTERSECTION_ERR;
			return err;
		}
	}
}

ErrorSet spiralArcIntx(Spiral sp, Arc arc, LLPointSet* pts, double tol, double eps) {

	double az12, az21, startAz, endAz, midAz, startRad, endRad, diff, split, dist, rad, crs1, crs2, initAz, span;
	LLPoint tempPt;
	LLPointSet testPts = createPtSet();
	Spiral tempSp;
	ErrorSet err = 0, tempErr = 0;
	int i, j = 0, count;

	LLPointPair pair;
	int n = 0;
	double bfR;

	//TODO:  Check if arc start/end points lie on spiral.

	err |= inverse(sp.centerPoint, arc.centerPoint, &az12, &az21, &dist, eps);

	//Arc is far from spiral, no intersection occurs
	if ((dist > arc.radius + sp.startRadius + tol) && (dist > arc.radius + sp.endRadius + tol)) {
		return err;
	}

	//Arc is fully inside spiral
	if ((dist + arc.radius + tol < sp.startRadius) && (dist + arc.radius + tol < sp.endRadius)) {
		return err;
	}
	//Spiral fully inside arc
	if ((dist + sp.startRadius + tol < arc.radius) && (dist + sp.endRadius + tol < arc.radius)) {
		return err;
	}

	startAz = sp.startAz;
	if (sp.subtendedAngle > M_PI) {
		//Split spiral into thirds
		split = 3;
	} else {
		split = 2;
	}
	diff = sp.subtendedAngle / split;
	for (j=0;j<split;j++) {
		endAz = startAz + sp.dir * diff;
		err |= spiralRadius(sp, startAz, &startRad);
		err |= spiralRadius(sp, endAz, &endRad);

		//TODO:  Possibly add logic to decide which arc (interior or exterior) to check first for intersections.
		err |= initArcIntx(arc.centerPoint, arc.radius, sp.centerPoint, startRad, pair, &n, &bfR, eps);

		midAz = (startAz + endAz) / 2;

		for (i=0;i<2;i++) {
			count = 0;
			if (n == 2) {
				err |= invCrs(sp.centerPoint, pair[0], &crs1, &az21, eps);
				err |= invCrs(sp.centerPoint, pair[1], &crs2, &az21, eps);
				initAz = (crs1 + crs2) / 2;
				span = fabs(azDifference(crs1, crs2));
				count = 1;
				break;
			} else if (n == 1) {
				count = 1;
				err |= invCrs(sp.centerPoint, pair[0], &initAz, &az21, eps);
				span = 0;
				break;
			} else {
				tempErr |= initArcIntx(arc.centerPoint, arc.radius, sp.centerPoint, endRad, pair, &n, &bfR, eps);
				if (i) {
					printf("Circle Intx Error\n");
					break;
				}
				i = 1;
			}
		}

		if (count > 0) {
			err |= spiralRadius(sp, midAz, &rad);
			err |= createSpiralSection(sp, midAz, rad, &tempSp, eps);
			tempErr |= initSpiralArcIntx(tempSp, arc, midAz, span, &testPts, tol, eps);
		}

		startAz = endAz;
	}

	if (testPts.length == 0) {
		err = tempErr;
	}

	for (i=0;i<testPts.length;i++) {
		tempPt.latitude = testPts.elements[i]->latitude;
		tempPt.longitude = testPts.elements[i]->longitude;
		if ((ptIsOnSpiral(sp, tempPt, tol, eps)) && (ptIsOnArc(arc.centerPoint, arc.radius, arc.startAz, arc.endAz, arc.dir, tempPt, &err, tol, eps))) {
			if (ptIsInSet(tempPt, *pts, tol) == 0) {
				err |= addPtToPtSet(pts, &tempPt);
			}
		}
	}
	return err;

}

ErrorSet initSpiralIntx(Spiral sp1, Spiral sp2, double initAz, double span, LLPointSet* pts, double tol, double eps) {

	double az12, az21, dist, rad, testAz, diff, azDiff, spCentCrs, revCentCrs;
	LLPoint spPt;
	int index = 0, count = 0;
	ErrorSet err = 0;

	double azArray[2] = { 0.0, 0.0 };
	double errArray[2] = { 0.0, 0.0 };

	err |= inverse(sp1.centerPoint, sp2.centerPoint, &spCentCrs, &revCentCrs, &dist, eps);

	//Spirals far from each other, no intersection occurs
	if ((dist > sp1.startRadius + sp2.startRadius + tol) && (dist > sp1.startRadius + sp2.endRadius + tol) && (dist > sp1.endRadius + sp2.startRadius + tol) && (dist > sp1.endRadius + sp2.endRadius + tol)) {
		printf("Arc far from spiral\n");
		return err;
	}

	//Spiral1 is fully inside spiral2
	if ((dist + sp1.startRadius + tol < sp2.startRadius) && (dist + sp1.startRadius + tol < sp2.endRadius) && (dist + sp1.endRadius + tol < sp2.startRadius) && (dist + sp1.endRadius + tol < sp2.endRadius)) {
		printf("Arc inside spiral\n");
		return err;
	}
	//Spiral2 is fully inside spiral1
	if ((dist + sp2.startRadius + tol < sp1.startRadius) && (dist + sp2.startRadius + tol < sp1.endRadius) && (dist + sp2.endRadius + tol < sp1.startRadius) && (dist + sp2.endRadius + tol < sp1.endRadius)) {
		printf("Arc inside spiral\n");
		return err;
	}

	/*if (ptIsOnSpiral(sp, arc.startPoint, tol, eps)) {
		err |= addPtToPtSet(pts, &arc.startPoint);
		count++;
	}
	if (ptIsOnSpiral(sp, arc.endPoint, tol, eps)) {
		err |= addPtToPtSet(pts, &arc.endPoint);
		count++;
	}
	if (ptIsOnArc(arc.centerPoint, arc.radius, arc.startAz, arc.endAz, arc.dir, sp.startPoint, &err, tol, eps)) {
		err |= addPtToPtSet(pts, &sp.startPoint);
		count++;
	}
	if (ptIsOnArc(arc.centerPoint, arc.radius, arc.startAz, arc.endAz, arc.dir, sp.endPoint, &err, tol, eps)) {
		err |= addPtToPtSet(pts, &sp.endPoint);
		count++;
	}

	count = 0;

	if (count >= 2) {
		return err;
	}
	*/

/*
	err |= spiralRadius(sp, spCentCrs, &rad);
	if (rad > dist) {
		inside = 1;
	}

	err |= projectToSpiral(sp, arc.centerPoint, &perpPt, tol, eps);
	if (err) {
		printf("Values:  %f %f %f %f\n", arc.centerPoint.latitude, arc.centerPoint.longitude, sp.centerPoint.latitude, sp.centerPoint.longitude);
	}
	err |= invDist(arc.centerPoint, perpPt, &dist, eps);

	//printf("Values:  %f %f\n", dist, arc.radius);
	//printf("LL:  %f %f %f %f\n", perpPt.latitude / (M_PI / 180), perpPt.longitude / (M_PI / 180), arc.centerPoint.latitude / (M_PI / 180), arc.centerPoint.longitude / (M_PI / 180));
	if (dist - arc.radius > tol) {
		//printf("No intersection\n");
		return err;
	} else if (fabs(dist - arc.radius) < tol){
		//printf("Tangent Intersection\n");
		//Tangent Intersection
		err |= addPtToPtSet(pts, &perpPt);
		return err;
	}
*/
	//printf("InitAz/Span:  %f %f\n", initAz, span);

	testAz = initAz + (span / 2);
	azArray[0] = testAz;

	err |= ptOnSpiral(sp1, testAz, &spPt, eps);
	err |= inverse(spPt, sp2.centerPoint, &az12, &az21, &dist, eps);
	err |= spiralRadius(sp2, az21, &rad);
	diff = dist - rad;
	errArray[0] = diff;

	testAz = testAz + fabs(diff) / diff * M_PI / 180;
	azArray[1] = testAz;

	err |= ptOnSpiral(sp1, testAz, &spPt, eps);
	err |= inverse(spPt, sp2.centerPoint, &az12, &az21, &dist, eps);
	err |= spiralRadius(sp2, az21, &rad);
	diff = dist - rad;
	errArray[1] = diff;

	//printf("TestAz: %f\n", testAz);

	while (1) {

		testAz = fmod(findRootSecantMethod(azArray, errArray, &err)+2*M_PI, 2*M_PI);

		err |= ptOnSpiral(sp1, testAz, &spPt, eps);
		err |= inverse(spPt, sp2.centerPoint, &az12, &az21, &dist, eps);
		err |= spiralRadius(sp2, az21, &rad);
		diff = dist - rad;

		azArray[0] = azArray[1];
		errArray[0] = errArray[1];

		azArray[1] = testAz;
		errArray[1] = diff;

		//err |= spiralRadius(sp, testAz, &rad);

		//printf("Test Az/Diff:  %f %f %f\n", testAz, diff, dist, rad);

		if (fabs(diff) < tol) {

			//printf("Dist/Rad/diff:  %f %f %e\n", dist, rad, diff);

			//printf("\n");
			if (count) {
				err |= addPtToPtSet(pts, &spPt);
				return err;
			}

			err |= addPtToPtSet(pts, &spPt);
			count++;

			err |= invCrs(sp1.centerPoint, sp2.centerPoint, &az12, &az21, eps);
			azDiff = azDifference(az12, testAz);

			//printf("Azs:  %f %f %f %f\n", testAz, fmod(az12 - azDiff, 2*M_PI), az12, azDiff);

			testAz = fmod(az12 - azDiff, 2*M_PI);
			//printf("TestAz: %f\n", testAz);

			azArray[0] = testAz;

			err |= ptOnSpiral(sp1, testAz, &spPt, eps);
			err |= inverse(spPt, sp2.centerPoint, &az12, &az21, &dist, eps);
			err |= spiralRadius(sp2, az21, &rad);
			diff = dist - rad;
			errArray[0] = diff;

			testAz = testAz + fabs(diff) / diff * azDiff / 10;
			azArray[1] = testAz;
			//printf("TestAz: %f\n", testAz);

			err |= ptOnSpiral(sp1, testAz, &spPt, eps);
			err |= inverse(spPt, sp2.centerPoint, &az12, &az21, &dist, eps);
			err |= spiralRadius(sp2, az21, &rad);
			diff = dist - rad;
			errArray[1] = diff;

			index = 0;
		}

		index++;

		if (index > 100) {
			err |= NO_INTERSECTION_ERR;
			return err;
		}
	}
}

ErrorSet spiralIntx(Spiral sp1, Spiral sp2, LLPointSet* pts, double tol, double eps) {

	double az12, az21, startAz, endAz, diff, dist, rad, crs1, crs2, initAz, span, count, split, midAz, rad1, rad2;
	LLPoint tempPt;
	LLPointSet testPts = createPtSet();
	Spiral tempSp;
	ErrorSet err = 0;
	int i, j, k;

	LLPointPair pair;
	int n = 0;
	double bfR;

	//TODO:  Check if arc start/end points lie on spiral.

	err |= inverse(sp1.centerPoint, sp2.centerPoint, &az12, &az21, &dist, eps);

	if (ptIsOnSpiral(sp1, sp2.startPoint, tol, eps)) {
		err |= addPtToPtSet(&testPts, &sp2.startPoint);
	}
	if (ptIsOnSpiral(sp1, sp2.endPoint, tol, eps)) {
		err |= addPtToPtSet(&testPts, &sp2.endPoint);
	}
	if (ptIsOnSpiral(sp2, sp1.startPoint, tol, eps)) {
		err |= addPtToPtSet(&testPts, &sp1.startPoint);
	}
	if (ptIsOnSpiral(sp2, sp1.endPoint, tol, eps)) {
		err |= addPtToPtSet(&testPts, &sp1.endPoint);
	}

	//Spirals far from each other, no intersection occurs
	if ((dist > sp1.startRadius + sp2.startRadius + tol) && (dist > sp1.startRadius + sp2.endRadius + tol) && (dist > sp1.endRadius + sp2.startRadius + tol) && (dist > sp1.endRadius + sp2.endRadius + tol)) {
		printf("Arc far from spiral\n");
		return err;
	}

	//Spiral1 is fully inside spiral2
	if ((dist + sp1.startRadius + tol < sp2.startRadius) && (dist + sp1.startRadius + tol < sp2.endRadius) && (dist + sp1.endRadius + tol < sp2.startRadius) && (dist + sp1.endRadius + tol < sp2.endRadius)) {
		printf("Arc inside spiral\n");
		return err;
	}
	//Spiral2 is fully inside spiral1
	if ((dist + sp2.startRadius + tol < sp1.startRadius) && (dist + sp2.startRadius + tol < sp1.endRadius) && (dist + sp2.endRadius + tol < sp1.startRadius) && (dist + sp2.endRadius + tol < sp1.endRadius)) {
		printf("Arc inside spiral\n");
		return err;
	}

	//TODO:  Possibly add logic to decide which arc (interior or exterior) to check first for intersections.
	err |= initArcIntx(sp1.centerPoint, sp1.startRadius, sp2.centerPoint, sp2.startRadius, pair, &n, &bfR, eps);

	//printf("N:  %i\n", n);

	if (sp1.subtendedAngle > M_PI) {
		//Split spiral into thirds
		split = 3;
	} else {
		split = 2;
	}

	diff = sp1.subtendedAngle / split;
	startAz = sp1.startAz;

	for (i=0;i<split;i++) {
		endAz = startAz + sp1.dir * diff;
		midAz = (startAz + endAz) / 2;
		for (j=0;j<2;j++) {
			err |= spiralRadius(sp1, startAz, &rad1);
			for (k=0;k<2;k++) {
				count = 0;
				rad2 = sp2.startRadius;
				err |= initArcIntx(sp1.centerPoint, rad1, sp2.centerPoint, rad2, pair, &n, &bfR, eps);
				if (n == 2) {
					err |= invCrs(sp1.centerPoint, pair[0], &crs1, &az21, eps);
					err |= invCrs(sp1.centerPoint, pair[1], &crs2, &az21, eps);
					if (fabs(azDifference(crs1, midAz)) < diff / 2 + M_PI / 18) {
						initAz = crs1;
						count++;
					}
					if (fabs(azDifference(crs2, midAz)) < diff / 2 + M_PI / 18) {
						initAz = crs2;
						count++;
					}
					if (count == 2) {
						initAz = (crs1 + crs2) / 2;
						span = fabs(azDifference(crs1, crs2));
					} else if (count == 1) {
						span = 0;
					}
					break;
				} else if (n == 1) {
					if (fabs(azDifference(crs2, midAz)) < diff / 2 + M_PI / 18) {
						count = 1;
						err |= invCrs(sp1.centerPoint, pair[0], &initAz, &az21, eps);
						span = 0;
					}
					break;
				}
				rad2 = sp2.endRadius;
			}
			if (count > 0) {
				break;
			}
			err |= spiralRadius(sp1, endAz, &rad2);
		}
		if (count > 0) {
			err |= spiralRadius(sp1, initAz, &rad);
			err |= createSpiralSection(sp1, initAz, rad, &tempSp, eps);
			err |= initSpiralIntx(tempSp, sp2, initAz, span, &testPts, tol, eps);
		}
		startAz = endAz;
	}
	for (i=0;i<testPts.length;i++) {
		tempPt.latitude = testPts.elements[i]->latitude;
		tempPt.longitude = testPts.elements[i]->longitude;
		/*int x, y;
		x = ptIsOnSpiral(sp, tempPt, tol, eps);
		y = ptIsOnArc(arc.centerPoint, arc.radius, arc.startAz, arc.endAz, arc.dir, tempPt, &err, tol, eps);
		//if (x + y < 2) {
			err |= invDist(arc.centerPoint, tempPt, &dist, eps);
			//printf("Spiral/Arc:  %i %i %e\n", x, y, arc.radius - dist);
		//}
		 */
		if ((ptIsOnSpiral(sp1, tempPt, tol, eps)) && (ptIsOnSpiral(sp2, tempPt, tol, eps))) {
			if (ptIsInSet(tempPt, *pts, tol) == 0) {
				err |= addPtToPtSet(pts, &tempPt);
			}
		}
	}
	return err;

}

} //namespace
