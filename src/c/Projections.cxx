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
 * Perpendicular Intercept
 *
 * Find location on geodesic at which course to given point is 90 degrees
 * different from geodesic course.
 *
 */
ErrorSet projectToGeo(LLPoint pt1, double geoStartAz, LLPoint pt3,
		LLPoint* pt2, double *crsFromPoint,
		double* distFromPoint, double tol, double eps)
{

	ErrorSet err = 0;

	// Spherical solution is first approximation
	LLPoint newPt1 = { 0.0, 0.0 };
	//    LLPoint testPt3 = { 0.0, 0.0 };
	double crs13, dist13, crs23, crs32, tmpCrs12;
	double crs21, dist12, crs31;
	double angle, error;
	double newDist;
	double dist23 = 0;
	double a, b, A, B, c; /* Angles of spherical triangle */
	double errarray[2];
	double distarray[2];
	double approxDist23;
	double npCrsFromPoint, npDistFromPoint;
	LLPoint npPt2;
	double startNbhdRadius = 1.0 / 1852.0; /* one meter in NM */
	double delta = 9.e99;
	double perpDistUpperBound = 1.0 / 1852.0; /* one meter in NM *///TODO Should this be part of the config struct???


	int pt1IsAtPole = 0;

	int k = 0;

#ifdef USE_BEST_FIT_ROC

	double sphereRad = lookUpROC(geocentricLat(pt3.latitude));
#else

	double sphereRad = SPHERE_RADIUS;
#endif

	/* Assign local storage if optional pointers are not provided */
	if (NULL == crsFromPoint) crsFromPoint = &npCrsFromPoint;
	if (NULL == distFromPoint) distFromPoint = &npDistFromPoint;
	if (NULL == pt2) pt2 = &npPt2;

	*crsFromPoint = 0;
	*distFromPoint = 0;

	/*determine if pt1 is at a pole */
	pt1IsAtPole = ptIsAtPole(pt1, &err, tol, eps);

	if (err |= inverse(pt1, pt3, &crs13, &crs31, &dist13, eps))
		return err;

	/* Check for perp intercept "behind" pt1 */
	angle = fabs(modlon(geoStartAz - crs13));

	/* Do approximate check for pt3 on geodesic */
	/* If it's close, then we check more carefully */
	approxDist23 = sphereRad * fabs(asin(sin(dist13 / sphereRad) * sin(angle)));

	if (dist13 <= tol)
	{
		/* pt2 is same as pt1 */
		*pt2 = pt1;
		return err;
	}
	else if (approxDist23 < (300.0 / 6076.0) && pt1IsAtPole == 0)
	{
		/* pt3 is near geodesic.  Move start point back 10 nm and check again */
		err |= direct(pt1, geoStartAz + M_PI, 10.0, &newPt1, eps);
		if (ptIsOnGeo(pt1, newPt1, pt3, INFINITE, &err, tol, eps))
		{
			/* point to be projected already lies on geodesic, so return it */
			/* NOTE: crsFromPoint undefined, distFromPoint == 0 in this case */
			*pt2 = pt3;
			return err;
		}
		else
		{
			/* point is near geodesic, but not within tol
			 * Use special approximation for small angle
			 */
		}
	}

	/* Check for orientation of start point/test point geometry
	 * Approximate spherical solution relies on having correct supplement
	 * of angle
	 */
	if (angle > M_PI_2)
	{
		B = M_PI - angle;
	}
	else
	{
		B = angle;
	}

	//determine the distance d12 of the first guess of pt2 from the pt1
	/* Check for situation where perp projected point is near p1
	 * This must be handled as special case to avoid numerical instabilities
	 */
	//    if (fabs(B - M_PI_2) < startNbhdRadius/sphereRad)

	if ((dist13 < startNbhdRadius) || (B > acos(
			tan(startNbhdRadius / sphereRad) / tan(dist13 / sphereRad))))
	{
		//B must not be allowed to equal A or the calculation of dist12 will
		//blow up (divide by 0). Need to handle this numerical boundary condition.
		//B = M_PI_2 - 1e-10; //shift slightly
		/* Approximate projected point will be within startNbhdRadius of pt1
		 * In this situation, spherical solution may fail, so set dist12 to trigger
		 * re-location of pt1 backwards along geodesic
		 */
		dist12 = 0.0;
	}
	else
	{
		/* Calculate spherical approximation to distance from pt1 to perp projection */
		a = dist13 / sphereRad;
		//        A = M_PI_2;
		//        b = asin(sin(B) * sin(a)); /* sin(A) = 1, so omitted from denominator */
		//        c = 2.0 * atan(tan(0.5 * (a - b)) * sin(0.5 * (A + B)) / sin(0.5 * (A
		//                - B))); //Napier's analogies (identities) from Spherical Trig
		c = atan(cos(B) * tan(a));
		if (c < 0.0)
			c = c + M_PI;
		dist12 = c * sphereRad;
	}

	if (angle > M_PI_2)
	{
		/* pt3 was behind pt1.  Need to move pt1 1.0 NM behind point that would
		 * be abeam pt3 */
		err |= direct(pt1, geoStartAz + M_PI, 5.0 + dist12, &newPt1, eps);
		dist12 = 5.0;
		err |= invCrs(newPt1, pt1, &geoStartAz, &crs21, eps);
		pt1 = newPt1;
		if (err)
			return err;
	}
	else if (fabs(dist12) < 5.0)
	{
		/* pt3 is within 5.0 nmi of being abeam pt1
		 * move pt1 backward 5 nmi to give the algorithms room to work */
		if (err |= direct(pt1, geoStartAz + M_PI, 5.0, &newPt1, eps))
			return err;
		dist12 = 5.0 + dist12;
		if (err |= invCrs(newPt1, pt1, &geoStartAz, &crs21, eps))
			return err;
		pt1 = newPt1;
	}

	//    //TODO Following case is not valid; uses spherical approximation as ellipsoidal solution
	//    else if (fabs(dist12) < tol)
	//    {
	//        /* pt3 is abeam pt1, so pt1 is projected point */
	//        *crsFromPoint = crs31;
	//        *distFromPoint = dist13;
	//        *pt2 = pt1;
	//        return err;
	//    }

	/*
	 * check if pt3 is very close to geodesic but not quite on the geodesic
	 * if so, then find the spherical approximation of pt2 and return this value
	 * A spherical approximation is sufficient for cases that fall within this case
	 * The perpDistUpperBound value was verified using the analysis described in the following bug
	 * See Geolib Bug 19576 for additional details
	 */
	if (approxDist23 > tol && approxDist23 < perpDistUpperBound)
	{

		//recalculate distances, courses, angles in case pt1 moved
		err |= inverse(pt1, pt3, &crs13, &crs31, &dist13, eps);
		angle = fabs(modlon(geoStartAz - crs13));
		if (angle > M_PI_2)
		{
			B = M_PI - angle;
		}
		else
		{
			B = angle;
		}

		/* Calculate spherical approximation of distance from pt1 to perp projection */
		a = dist13 / sphereRad;
		//        A = M_PI_2;
		//        b = asin(sin(B) * sin(a)); /* sin(A) = 1, so omitted from denominator */
		//        c = 2.0 * atan(tan(0.5 * (a - b)) * sin(0.5 * (A + B)) / sin(0.5 * (A
		//                - B))); //Napier's analogies (identities) from Spherical Trig
		c = atan(cos(B) * tan(a));
		if (c < 0.0)
			c = c + M_PI;
		dist12 = c * sphereRad;

		//find the projection point of pt3 on the geodesic using the spherical distance approx
		err |= direct(pt1, geoStartAz, dist12, pt2, eps);

		//determine the course and distance info with respect to pt2
		err
		|= inverse(*pt2, pt3, &crs23, crsFromPoint, distFromPoint,
				eps);
		return err;
	}
	else if (err |= direct(pt1, geoStartAz, dist12, pt2, eps))
		return err;

	/* Calculate angle between radial and approximate perpendicular */
	err |= inverse(*pt2, pt1, &crs21, &tmpCrs12, &dist12, eps);
	err |= inverse(*pt2, pt3, &crs23, &crs32, &dist23, eps);
	if (err)
		return err;
	/* Cast angle between main course and perpendicular into range [-Pi,Pi] */
	angle = fabs(modlon(crs21 - crs23));
	errarray[0] = angle - M_PI_2;
	distarray[0] = dist12;

	distarray[1] = distarray[0] + errarray[0] * dist23;

	if (err |= direct(pt1, geoStartAz, distarray[1], pt2, eps))
		return err;

	// Calculate angle between radial and approximate perpendicular
	err |= inverse(*pt2, pt1, &crs21, &tmpCrs12, &dist12, eps);
	err |= inverse(*pt2, pt3, &crs23, &crs32, &dist23, eps);
	if (err)
		return err;

	/* Cast angle between main course and perpendicular into range [-Pi,Pi] */
	angle = modlon(crs21 - crs23);
	errarray[1] = dist23 * (fabs(angle) - M_PI_2);
	error = errarray[1];

	while ((fabs(error) > tol || (delta > tol)) &&
			(k < MAX_ITERATIONS)
	)
	{
		newDist = findRootSecantMethod(distarray, errarray, &err);

		if(delta == 0){
			// If the iteration stops progressing but error > tol then move a little bit to get restarted.
			newDist += tol;
		}

		if (err |= direct(pt1, geoStartAz, newDist, pt2, eps))
			return err;
		/* Calculate angle between given line and approximate perpendicular */
		if (err |= inverse(*pt2, pt1, &crs21, &tmpCrs12, &dist12, eps))
			return err;
		if (err |= inverse(*pt2, pt3, &crs23, &crs32, &dist23, eps))
			return err;
		/* Cast angle between main course and perpendicular into range [-Pi,Pi] */
		angle = modlon(crs21 - crs23);

		//        if (angle < 0.0)
		//        {
		//            crs23 = crs21 + M_PI_2;
		//        }
		//        else
		//        {
		//            crs23 = crs21 - M_PI_2;
		//        }
		/* Error is distance between given test point and point projected out
		 * at 90 degree angle from geodesic */
		//        err |= direct(*pt2, crs23, dist23, &testPt3, eps);
		//        err |= inverse(pt3, testPt3, NULL, NULL, &error, eps);

		error = dist23 * (fabs(angle) - M_PI_2);

		errarray[0] = errarray[1];
		distarray[0] = distarray[1];

		/* error function has same shape as absolute value function
		 * Need to make it smooth for convergence */
		//        if (fabs(angle) < M_PI_2)
		//        {
		//            errarray[1] = -error;
		//        }
		//        else
		//        {
		errarray[1] = error;
		//        }
		//        errarray[1] = dist23 * (fabs(angle) - M_PI_2);
		distarray[1] = newDist;

		delta = fabs(distarray[0] - distarray[1]);
		//        if ((fabs(error) <= tol) && (fabs(delta) <= tol))
		//        {
		//            break;
		//        }
		//        error = fabs(distarray[1] - distarray[0]);

		//        if (k>0)
		//        printf("%d, %.15e, %.15e, %.15e, %.15e\n",k,newDist,(fabs(angle)-M_PI_2)*180/M_PI,dist23,error);

		k++;

	}

	if (k >= MAX_ITERATIONS)
	{
		err |= ITERATION_MAX_REACHED_ERR;
		//        printf("Error: ITERATION_MAX_REACHED in %s\n",__FUNCTION__);
	}

	if (fabs(error) >= MAX_DISTANCE_ERROR)
	{
		err |= ERROR_MAX_REACHED_ERR;
	}

	*crsFromPoint = crs32;
	*distFromPoint = dist23;

	return err;

}

/*******************************************************************************
 *
 * */

ErrorSet projectToGeoAtAngle(LLPoint pt1, double crs12, LLPoint pt3,
		double intAngle, LLPoint* pt2, double* crsFromPoint,
		double* distFromPoint, double tol, double eps)
{

	LLPoint newStart;
	LLPoint tmpPt3, tmpEndPt;
	double actAngle, angleError;
	double tmpCrs12, tmpCrs21, tmpDist12;
	double c, crs21, crs23;
	double crs13, crs31, dist13, angle13;
	double distarray[2], errarray[2];
	double error;
	double stepSize = 9.0e99;
	ErrorSet errCode = 0;
	int k;

	double npCrsFromPoint, npDistFromPoint;
	LLPoint npPt2;

	/* Assign local storage if optional pointers are not provided */
	if (NULL == crsFromPoint) crsFromPoint = &npCrsFromPoint;
	if (NULL == distFromPoint) distFromPoint = &npDistFromPoint;
	if (NULL == pt2) pt2 = &npPt2;

	if (fabs(intAngle) < tol)
	{
		return INVALID_CRS_ERR;
	}

	/* Ignore sign of intAngle for now */
	intAngle = fabs(intAngle);

	/* Find perpendicular projection (use smaller tolerance for better precision */
	errCode |= projectToGeo(pt1, crs12, pt3, pt2, crsFromPoint,
			distFromPoint, tol, eps);

	/* If intercept angle is 90 degrees, then projectToGeo
	 * gives us the answer */
	if ((errCode) || (fabs(intAngle - M_PI_2) < tol))
	{
		return errCode;
	}

	/* Calculate approx distance from perp projection to angled projection *
	 * This approximation uses spherical triangles and Napier's Circular Parts */
	/* NOTE: If intAngle > pi/2, then c < 0 */
	c = SPHERE_RADIUS * asin(tan(*distFromPoint / SPHERE_RADIUS) / tan(
			intAngle));

	/* Move startPt to perpendicular projected point */
	errCode |= inverse(pt1, *pt2, &tmpCrs12, &tmpCrs21, &tmpDist12, eps);
	if (errCode)
		return errCode;

	if (fabs(modlon(crs12 - tmpCrs12)) > M_PI_2)
	{
		/* Treat distance as negative if perp. point lies behind pt1. */
		tmpDist12 = c - tmpDist12;
	}
	else
	{
		tmpDist12 = c + tmpDist12;
	}

	/* Case where either perpendicular projected point is behind pt1, or approx. angled
	 * projected point is behind pt1, or or is too close to
	 * pt1, causes problems.  Move pt1 back sufficient and recalculate inputs */
	if (tmpDist12 < 1.0) // 2.0 nm threshold was chosen arbitrarily
	{
		/* Move start point backwards, away from angled projection, to give the
		 * algorithm more "room to work."  Makes azimuth calculations more accurate */
		//        printf("c = %20.15f, tmpDist12 = %20.15f, c + tmpDist12 = %20.15f\n",c,tmpDist12,c+tmpDist12);
		/* Step 1: Set temporary end point sufficiently far away to give precise
		 *         course from new location of pt1. */
		errCode |= direct(pt1, crs12, 5.0, &tmpEndPt, eps);
		/* Step 2: Relocate pt1 correct distance "behind" given location */
		errCode
		|= direct(pt1, crs12 + M_PI, 2.0 - tmpDist12, &newStart, eps);
		/* Step 3: Recalculate crs12 (azimuth of geodesic at pt1) for new location */
		errCode |= inverse(newStart, tmpEndPt, &tmpCrs12, &tmpCrs21,
				&tmpDist12, eps);
		if (errCode)
			return errCode;
		pt1 = newStart;
		crs12 = tmpCrs12;
		/* pt1 has been moved to 2.0 nm behind that approx projected point */
		tmpDist12 = 2.0;
	}

	/* Set correct sign on input angle */
	errCode |= inverse(pt1, pt3, &crs13, &crs31, &dist13, eps);
	angle13 = modlon(crs12 - crs13);
	intAngle = -sgn(angle13) * fabs(intAngle);

	k = 0;

	//    printf("[ ");

	while ((k == 0) || (((fabs(error) > tol) || (stepSize > tol)) && (k < MAX_ITERATIONS)
	))
	{

		/* If tmpDist12 < 0, then direct will use reciprocal course. */
		errCode |= direct(pt1, crs12, tmpDist12, pt2, eps);
		if (errCode)
			return errCode;

		//		printf("k = %d: ",k);
		//		_display(*pt2);

		/* Compute azimuths from approx point to pt1 and pt3 */
		errCode |= invCrs(*pt2, pt1, &crs21, &tmpCrs12, eps);
		errCode |= inverse(*pt2, pt3, &crs23, crsFromPoint,
				distFromPoint, eps);
		if (errCode)
			return errCode;

		actAngle = modlon(crs21 - crs23);
		//		printf("actAngle = %20.15f; ",actAngle*180.0/M_PI);

		/* Place point from pt2 along desired angle.  This will match pt3 if pt2 is correct */
		errCode |= direct(*pt2, crs21 - intAngle, *distFromPoint, &tmpPt3,
				eps);
		errCode |= invDist(pt3, tmpPt3, &error, eps);
		if (errCode)
			return errCode;

		/* Give sign to error */
		/* If intAngle < actAngle, then pt2 is too close to p1 and error < 0 */
		/* Otherwise, pt2 is too far and error > 0. */
		angleError = fabs(intAngle) - fabs(actAngle);
		//        error = sgn(angleError) * error;

		//        error = sgn(angleError)*SPHERE_RADIUS*asin(sin(angleError)*sin(*distFromPoint/SPHERE_RADIUS)/sin(fabs(intAngle)));

		//        error = *distFromPoint * sin(angleError);

		if (k == 0)
		{
			distarray[1] = tmpDist12;
			errarray[1] = error;

			/* Estimate improved distance */
			/* Uses spherical triangle approximation */
			error = SPHERE_RADIUS * asin(sin(angleError) * sin(*distFromPoint
					/ SPHERE_RADIUS) / sin(intAngle));
			/* error < 0 if pt2 is too close to pt1 */

			tmpDist12 = tmpDist12 + error;
		}
		else
		{
			distarray[0] = distarray[1];
			distarray[1] = tmpDist12;
			errarray[0] = errarray[1];
			errarray[1] = error;

			tmpDist12 = findRootSecantMethod(distarray, errarray, &errCode);
			if (errCode)
				return errCode;
		}
		stepSize = fabs(distarray[1] - distarray[0]);
		//        printf("%20.15f %e;\n",distarray[1],errarray[1]);

		k++;

	}

	if (k >= MAX_ITERATIONS)
	{
		errCode |= ITERATION_MAX_REACHED_ERR;
		//        printf("Error: ITERATION_MAX_REACHED in %s\n",__FUNCTION__);

	}

	if (fabs(error) >= MAX_DISTANCE_ERROR)
	{
		errCode |= ERROR_MAX_REACHED_ERR;
	}

	//    printf("];\n");

	return errCode;

}

/******************************************************************************
 * Locus Perpendicular Intercept
 *
 * Find location on Locus at which course to given point is 90 degrees
 * different from Locus course.
 *
 */

ErrorSet projectToLocus(Locus loc, LLPoint pt2, LLPoint* locPt,
		double *crsFromPoint, double* distFromPoint,
		double tol, double eps)
{

	ErrorSet err = 0;

	// Spherical solution is first approximation
	LLPoint perpPt;
	LLPoint geoPt;

	double gcrs, bgcrs, gdist;
	double lcrs, blcrs;
	double locAngle;

	double crsToPoint, angle;
	double crsFromGeoStart;
	double newDist = 0.0;

	double errarray[2];
	double distarray[2];

	int k = 0;
	int online = 0;

	LLPoint npLocPt;
	double npCrsFromPoint;
	double npDistFromPoint;

	/* Assign local storage if optional pointers are not provided */
	if (NULL == locPt) locPt = &npLocPt;
	if (NULL == crsFromPoint) crsFromPoint = &npCrsFromPoint;
	if (NULL == distFromPoint) distFromPoint = &npDistFromPoint;

	*crsFromPoint = 0;
	*distFromPoint = 0;

	/* Find course & length of defining geodesic.  Must do this in any case */
	//Un-necessary call because course & length are part of the locus struct.
	//err |= inverse(loc.geoStart, loc.geoEnd, &gcrs, &bgcrs, &gdist, eps);
	gcrs = loc.geoAz;
	bgcrs = loc.geoRevAz;
	gdist = loc.geoLength;
	if (fabs(loc.startDist - loc.endDist) < tol)
		/* Locus is "parallel" to geodesic; perp point is easily found from
		 * corresponding point on geodesic */
	{

		/* Project from given point to geodesic      */
		err |= projectToGeo(loc.geoStart, loc.geoAz, pt2, &perpPt,
				crsFromPoint, distFromPoint, tol, eps);
		if (err)
		{
			err |= NO_PROJECTED_POINT_ERR;
			return err;
		}

		online = ptIsOnGeo(loc.geoStart, loc.geoEnd, perpPt,
				loc.lineType, &err, tol, eps);
		if (err)
		{
			err |= NO_PROJECTED_POINT_ERR;
			return err;
		}
		/* Projected point on locus is on geodesic joining perpPt to pt2*/
		/* Find position of point on locus */
		err |= ptOnLocusFromGeoPt(loc, perpPt, locPt, &crsToPoint, tol, eps);
		if (err)
		{
			err |= NO_PROJECTED_POINT_ERR;
			return err;
		}

		/* Calculate distance and course from given point to projected point.
		 * This value is returned by reference */
		err |= inverse(*locPt, pt2, &crsToPoint, crsFromPoint,
				distFromPoint, eps);

	}
	else /* Locus is NOT "parallel" to geodesic. Must iterate */
	{
		/* Find approximate point based on geodesic approx. to locus */
		err |= invCrs(loc.locusStart, loc.locusEnd, &lcrs, &blcrs, eps);
		//printf("Geo approx Lcrs = %f\n", lcrs*180.0/M_PI);
		/* first approximation: project pt2 onto geodesic approximation *
		 * NOTE: this point is not on locus yet */
		if (err |= projectToGeo(loc.locusStart, lcrs, pt2, locPt,
				crsFromPoint, distFromPoint, tol, eps))
			return err;

		/* Locus inclination angle, relative to defining geodesic.
		 * This is used for first improvement step */
		locAngle = atan((loc.startDist - loc.endDist) / gdist);

		/* point projected on geodesic from first approximation */
		if (err |= projectToGeo(loc.geoStart, gcrs, *locPt, &geoPt,
				crsFromPoint, distFromPoint, tol, eps))
			return err;

		err |= inverse(loc.geoStart, geoPt, &crsFromGeoStart, NULL, &distarray[1], eps);
		if (fabs(modlon(loc.geoAz-crsFromGeoStart)) > M_PI_2)
		{
			/* Projected point is behind start of locus, so negate distance */
			distarray[1] = -distarray[1];
		}
		newDist = distarray[1];

		k = 0;

		//        printf("\n%% Begin Loop -----------------\n");


		while ((k == 0) ||   ( (fabs(errarray[1]) > tol) &&
				(k < MAX_ITERATIONS) ))
		{

			//            if (k > 0) /* project distarray along geodesic */
			//            {
			//                err |= direct(loc.geoStart, gcrs, newDist, &geoPt, eps);
			//                if (err)
			//                    return err;
			//            }
			//
			//            /* Find point on locus */
			//            if (err
			//                    |= ptOnLocusFromGeoPt(loc, geoPt, locPt, &crsToPoint, tol, eps))
			//                return err;
			//
			//            /* Compute course of locus at current locPt.
			//             * Also returns current correspoding point on geodesic */
			//            lcrs = locusCrsAtPt(loc, *locPt, &geoPt, &crsToPoint, &err,
			//                    tol, eps);

			/* This function also returns a point on the locus that is a candidate for the
			 * perp projection of the given point.  */
			lcrs = locusCrsAtGeoDist(loc,newDist,locPt,NULL,&err,tol,eps);

			/* Return value = -1 from locusCrsAtPt is redundant signal that computation failed */
			if ((err > 0) || (lcrs < 0.0))
			{
				/* lcrs set to -1 only if err != 0 is returned */
				/* Value of err gives reason for failure       */
				return err;
			}

			if (err |= inverse(*locPt, pt2, &crsToPoint, crsFromPoint,
					distFromPoint, eps))
				return err;

			/* Independent variable in iteration is along-geodesic distance to current approx. point */
			/* error metric is difference from 90 deg angle */
			angle = modlon(lcrs - crsToPoint);

			errarray[1] = -(*distFromPoint) * cos(angle);

			//            if (k>0)
			//            {
			//                printf("%d\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\n",k,newDist,errarray[1],distToLocusFromGeoDist(loc,newDist),lcrs,locPt->latitude,locPt->longitude,*distFromPoint);
			//                fflush(stdout);
			//            }

			/* Find distance along geodesic for improved approximation to perp locus point */
			if (k == 0)
			{
				/* Use direct trig calculation for first update */
				newDist = distarray[1] + errarray[1] * cos(locAngle);
			}
			else
			{
				newDist = findRootSecantMethod(distarray, errarray, &err);
			}
			distarray[0] = distarray[1];
			errarray[0] = errarray[1];
			distarray[1] = newDist;

			k++;
		}
		//        printf("\n%% End Loop -------------------\n");

		if (k >= MAX_ITERATIONS)
		{
			err |= ITERATION_MAX_REACHED_ERR;
			//            printf("Error: ITERATION_MAX_REACHED in %s\n",__FUNCTION__);

		}

		if (fabs(errarray[1]) >= MAX_DISTANCE_ERROR)
		{
			err |= ERROR_MAX_REACHED_ERR;
		}

	}

	return err;

}

/*******************************************************************************
 * Given an arc and a geodesic, find tangent points on arc such that geodesics
 * through those points cross given geodesic at right angle.  Also return the
 * geodesic intersection points.
 *
 * INPUT:
 * lineStart: start point of geodesic
 *       crs: azimuth of geodesic at lineStart
 *    center: center point of arc
 *    radius: radius of arc
 *   linePts: Two element array to store and return geodesic intersection points
 *    tanPts: Two element array to store and return tangent points
 *       tol: Required accuracy tolerance
 *       eps: forward/direct algorithm tolerance
 *            (soon to be deprecated and eventually removed)
 *
 * OUTPUT:
 *  none
 *
 */

ErrorSet projectArcTanPtsToGeo(LLPoint lineStart, double crs, LLPoint center,
		double radius, LLPointPair linePts,
		LLPointPair tanPts, double tol, double eps)
{

	ErrorSet err = 0;

	LLPoint perpPt;
	LLPoint tempPt;

	double perpCrs, perpDist;
	double crs12, crs21, angle1, dist12;
	double crsStartToCenter, distStartToCenter, temp;
	double crsCenterToStart;
	double error, delta;
	double radCrs, radDist, strCrs;
	int i, k;
	LLPointPair npTanPts;
	LLPointPair npLinePts;

	/* Assign local storage if optional pointers are not provided */
	if (NULL == tanPts) tanPts = npTanPts;
	if (NULL == linePts) linePts = npLinePts;

	//TODO: With hard-coded eps, this check should not be necessary
	if (tol < eps)
	{
		err |= TOL_TOO_SMALL_ERR;
		return err;
	}

	/* Handle the case for radius < tol, in which the arc is simply a point
	 * and is therefore not a valid shape. */

	if (radius < tol)
	{
		err |= INVALID_SHAPE_ERR;
		return err;
	}

	/* Find point on line perpendicular from arc center */
	if (err |= projectToGeo(lineStart, crs, center, &perpPt, &perpCrs,
			&perpDist, tol, eps))
	{
		err |= NO_PROJECTED_POINT_ERR;
		return err;
	}

	/* Find point on line perpendicular from arc center */
	if (err |= projectToGeo(lineStart, crs, center, &perpPt, &perpCrs,
			&perpDist, tol, eps))
	{
		err |= NO_PROJECTED_POINT_ERR;
		return err;
	}

	if (err |= inverse(lineStart, center, &crsStartToCenter,
			&crsCenterToStart, &distStartToCenter, eps))
		return err;

	/* angle1 tells us on which side of the line lie the points we seek */
	angle1 = modlon(crs - crsStartToCenter);

	if (fabs(crsStartToCenter - crs) < INTERNAL_ZERO)
	{
		/* given line is diameter of circle */
		err |= direct(lineStart, crs, distStartToCenter - radius,
				&tanPts[0], eps);
		err |= direct(lineStart, crs, distStartToCenter + radius,
				&tanPts[1], eps);
		linePts[0] = tanPts[0];
		linePts[1] = tanPts[1];
		return err;
	}

	err |= inverse(lineStart, perpPt, &crs12, &crs21, &dist12, eps);

	if (dist12 < 1.0 / 1852.0)
	{
		/* Projection of arc center onto line is too close to
		 * line start point to allow precise course computations.
		 * Find new point tmpStart on line 1 NM from lineStart and compute
		 * courses using this point.
		 */
		LLPoint tmpStart; /* This declaration is inside block, so OK in Visual C */
		err |= direct(lineStart, crs, 1.0, &tmpStart, eps);
		err |= inverse(tmpStart, perpPt, &crs12, &crs21, &dist12, eps);

	}

	/* Now crs21 is course of given geodesic at perp projection of centerpoint */

	if (err)
	{
		return err;
	}

	delta = radius;

	/* Loop through this twice to find two perp/tangent point pairs */
	for (i = 0; i < 2; i++)
	{

		k = 0;

		/* Iterate to improve approximation of i-th pair */
		while ((k == 0) || ((fabs(error) > tol) && (k < MAX_ITERATIONS)
		))
		{

			if (i == 0)
			{
				if (err |= direct(perpPt, crs21 + M_PI, delta, &linePts[i],
						eps))
					return err;
			}
			else
			{
				if (err |= direct(perpPt, crs21, delta, &linePts[i], eps))
					return err;
			}

			if (err |= invCrs(linePts[i], perpPt, &strCrs, &temp, eps))
				return err;

			if (angle1 > 0)
			{
				/* circle center is to left of line */
				perpCrs = strCrs - M_PI_2;
			}
			else
			{
				/* center is to right of line */
				perpCrs = strCrs + M_PI_2;
			}

			/* find perp intercept from arc center to tangent line */
			if (err |= projectToGeo(linePts[i], perpCrs, center, &tempPt,
					&radCrs, &radDist, tol, eps))
				return err;

			tanPts[i] = tempPt;

			/* too far from center == positive error */
			error = radDist - radius;

			/* positive error => smaller delta => move closer to center */
			delta = delta - error;

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

	}

	return err;

}

//Finds the point on the spiral such that the course from the spiral point to the input point is perpendicular
//to the spiral course at the spiral point.  This function is called by the higher level projectToSpiral.
//The input spiral is created by projectToSpiral such that no more than 1 point exists.
static ErrorSet initProjectToSpiral(Spiral sp, LLPoint pt, LLPoint* point, double testAz, double tol, double eps) {

	ErrorSet err = 0;
	double rad, az12, az21, ptDist, spCrs, azDiff, dPhi, theta, ptCentDist;
	int inSpiral = 1;
	LLPoint spPt;

	int counter = 0;

	err |= spiralTanCrs(sp, testAz, &spCrs);
	err |= ptOnSpiral(sp, testAz, &spPt, eps);
	err |= spiralRadius(sp, testAz, &rad);

	err |= inverse(sp.centerPoint, pt, &az12, &az21, &ptCentDist, eps);

	if (ptCentDist < rad) {
		inSpiral = -1;
	}

	err |= inverse(spPt, pt, &az12, &az21, &ptDist, eps);

	azDiff = azDifference(az12, spCrs);

	while (1) {

		err |= spiralTanCrs(sp, testAz, &spCrs);
		err |= ptOnSpiral(sp, testAz, &spPt, eps);
		err |= spiralRadius(sp, testAz, &rad);
		err |= inverse(spPt, pt, &az12, &az21, &ptDist, eps);

		if (inSpiral < 0) {
			az12 = fmod(az12 + M_PI, 2*M_PI);
		}



		azDiff = azDifference(az12, spCrs);

		theta = M_PI / 2 - fabs(azDiff);

		if (inSpiral < 0) {
			dPhi = theta / (1 + (rad / ptDist));
		} else {
			dPhi = theta * (1 - inSpiral * (rad / (ptDist + rad)));
		}

		testAz = fmod(testAz + inSpiral * fabs(azDiff) / azDiff * dPhi, 2*M_PI);

		if ((fabs(ptDist * theta) * 2 < tol) && (fabs(dPhi) * ptDist * 2 < tol)) {
			*point = spPt;
			return err;
		}
		if (counter == 1000) {
			err |= NO_INTERSECTION_ERR;
			return err;
		}
		counter++;
	}
}

//Finds the point on the spiral such that the course from the spiral point to the input point is perpendicular
//to the spiral course at the spiral point.  This function calls initProjectToSpiral to determine all candidate points
//and then checks to determine which of them lie on the input spiral.
ErrorSet projectToSpiral(Spiral sp, LLPoint pt, LLPoint* perpPt, double tol, double eps) {

	ErrorSet err = 0;
	double testAz, az21, tempDist, rad;
	LLPoint pt1, pt2;
	Spiral tempSp;

	err |= inverse(sp.centerPoint, pt, &testAz, &az21, &tempDist, eps);
	err |= spiralRadius(sp, testAz, &rad);
	err |= createSpiralSection(sp, testAz, rad, &tempSp, eps);
	err |= initProjectToSpiral(tempSp, pt, &pt1, testAz, tol, eps);

	if (ptIsOnSpiral(sp, pt1, tol, eps)) {
		*perpPt = pt1;
		return err;
	}

	if ((tempSp.startRadius - sp.startRadius) < 0) {
		rad = rad - fabs(sp.growthRate*2*M_PI);
	} else {
		rad = rad + fabs(sp.growthRate*2*M_PI);
	}
	err |= createSpiralSection(sp, testAz, rad, &tempSp, eps);
	err |= initProjectToSpiral(tempSp, pt, &pt2, testAz, tol, eps);

	if (ptIsOnSpiral(sp, pt1, tol, eps)) {
		*perpPt = pt2;
		return err;
	}
	err |= NO_INTERSECTION_ERR;
	return err;
}

ErrorSet distBetweenGeos(Geodesic line1, Geodesic line2, double *dist, double tol, double eps){

	ErrorSet err = 0;
	double closestDistance = 99e99;
	double crs31, dist13;
	double crs32, dist23;
	LLPoint intersectionPoint;

	double line1pt1dist;
	double line1pt2dist;
	double line2pt1dist;
	double line2pt2dist;
	LLPoint line1pt1proj;
	LLPoint line1pt2proj;
	LLPoint line2pt1proj;
	LLPoint line2pt2proj;

	double line1pt1line2pt1dist;
	double line1pt2line2pt1dist;
	double line1pt1line2pt2dist;
	double line1pt2line2pt2dist;

	// Check for intersection
	err |= geoIntx(line1.startPoint, line1.endPoint, SEGMENT, &crs31, &dist13,
			line2.startPoint, line2.endPoint, SEGMENT, &crs32, &dist23, &intersectionPoint, tol, eps);

	if(err != NO_INTERSECTION_ERR){
		closestDistance = 0.0;
	} else {
		err = 0;

		// Check end point projections
		err |= projectToGeo(line2.startPoint, line2.startAz, line1.startPoint, &line1pt1proj, NULL, &line1pt1dist, tol, eps);
		err |= projectToGeo(line2.startPoint, line2.startAz, line1.endPoint, &line1pt2proj, NULL, &line1pt2dist, tol, eps);
		err |= projectToGeo(line1.startPoint, line1.startAz, line2.startPoint, &line2pt1proj, NULL, &line2pt1dist, tol, eps);
		err |= projectToGeo(line1.startPoint, line1.startAz, line2.endPoint, &line2pt2proj, NULL, &line2pt2dist, tol, eps);

		if(ptIsOnGeo(line2.startPoint, line2.endPoint, line1pt1proj, SEGMENT, &err, tol, eps))
			closestDistance = fmin(closestDistance, line1pt1dist);
		if(ptIsOnGeo(line2.startPoint, line2.endPoint, line1pt2proj, SEGMENT, &err, tol, eps))
			closestDistance = fmin(closestDistance, line1pt2dist);
		if(ptIsOnGeo(line1.startPoint, line1.endPoint, line2pt1proj, SEGMENT, &err, tol, eps))
			closestDistance = fmin(closestDistance, line2pt1dist);
		if(ptIsOnGeo(line1.startPoint, line1.endPoint, line2pt2proj, SEGMENT, &err, tol, eps))
			closestDistance = fmin(closestDistance, line2pt2dist);

		if(closestDistance == 99e99){
			// Check end points to end points
			err |= invDist(line1.startPoint, line2.startPoint, &line1pt1line2pt1dist, eps);
			err |= invDist(line1.endPoint, line2.startPoint, &line1pt2line2pt1dist, eps);
			err |= invDist(line1.startPoint, line2.endPoint, &line1pt1line2pt2dist, eps);
			err |= invDist(line1.endPoint, line2.endPoint, &line1pt2line2pt2dist, eps);

			closestDistance = fmin(closestDistance, line1pt1line2pt1dist);
			closestDistance = fmin(closestDistance, line1pt2line2pt1dist);
			closestDistance = fmin(closestDistance, line1pt1line2pt2dist);
			closestDistance = fmin(closestDistance, line1pt2line2pt2dist);
		}
	}

	*dist = closestDistance;

	return err;
}

ErrorSet distBetweenArcGeo(Arc arc, Geodesic line, double *dist, double tol, double eps){

	ErrorSet err = 0;
	double closestDistance = 99e99;
	LLPointPair intx;
	int n;
	int i = 0;
	LLPoint centerProj;
	double crsToLine, arcCenterToLineDist;
	double arcCenterToLineStartCrs, arcCenterToLineStartDist;
	double arcCenterToLineEndCrs, arcCenterToLineEndDist;
	double arcStartProjDist, arcEndProjDist;
	LLPoint arcStartProj, arcEndProj;
	double lineStartArcStartDist, lineStartArcEndDist, lineEndArcStartDist, lineEndArcEndDist;

	// Check for intersection
	err = geoArcIntx(line.startPoint, line.startAz, arc.centerPoint, arc.radius, intx, &n, tol, eps);
	for(i = 0; i < n; i++){
		if(ptIsOnGeo(line.startPoint, line.endPoint, intx[i], SEGMENT, &err, tol, eps)
				&& ptIsOnArc(arc.centerPoint, arc.radius, arc.startAz, arc.endAz, arc.dir, intx[i], &err, tol, eps)){
			*dist = 0;
			return err;
		}
	}

	// Check projections
	err |= projectToGeo(line.startPoint, line.startAz, arc.centerPoint, &centerProj, &crsToLine, &arcCenterToLineDist, tol, eps);
	if(azIsInArcExtent(arc, crsToLine)
			&& ptIsOnGeo(line.startPoint, line.endPoint, centerProj, SEGMENT, &err, tol, eps)){
		closestDistance = fmin(closestDistance,
				fabs(arcCenterToLineDist - arc.radius));
	}

	err |= inverse(arc.centerPoint, line.startPoint, &arcCenterToLineStartCrs, NULL, &arcCenterToLineStartDist, eps);
	if(azIsInArcExtent(arc, arcCenterToLineStartCrs) || arcCenterToLineStartDist == 0){
		closestDistance = fmin(closestDistance,
				fabs(arcCenterToLineStartDist - arc.radius));
	}

	err |= inverse(arc.centerPoint, line.endPoint, &arcCenterToLineEndCrs, NULL, &arcCenterToLineEndDist, eps);
	if(azIsInArcExtent(arc, arcCenterToLineEndCrs) || arcCenterToLineEndDist == 0){
		closestDistance = fmin(closestDistance,
				fabs(arcCenterToLineEndDist - arc.radius));
	}

	err |= projectToGeo(line.startPoint, line.startAz, arc.startPoint, &arcStartProj, NULL, &arcStartProjDist, tol, eps);
	if(ptIsOnGeo(line.startPoint, line.endPoint, arcStartProj, SEGMENT, &err, tol, eps)){
		closestDistance = fmin(closestDistance, arcStartProjDist);
	}

	err |= projectToGeo(line.startPoint, line.startAz, arc.endPoint, &arcEndProj, NULL, &arcEndProjDist, tol, eps);
	if(ptIsOnGeo(line.startPoint, line.endPoint, arcEndProj, SEGMENT, &err, tol, eps)){
		closestDistance = fmin(closestDistance, arcEndProjDist);
	}

	if(closestDistance == 99e99){
		// Check end points to end points
		err |= invDist(line.startPoint, arc.startPoint, &lineStartArcStartDist, eps);
		err |= invDist(line.startPoint, arc.endPoint, &lineStartArcEndDist, eps);
		err |= invDist(line.endPoint, arc.startPoint, &lineEndArcStartDist, eps);
		err |= invDist(line.endPoint, arc.endPoint, &lineEndArcEndDist, eps);

		closestDistance = fmin(closestDistance, lineStartArcStartDist);
		closestDistance = fmin(closestDistance, lineStartArcEndDist);
		closestDistance = fmin(closestDistance, lineEndArcStartDist);
		closestDistance = fmin(closestDistance, lineEndArcEndDist);
	}

	*dist = closestDistance;
	return err;
}
} //namespace
