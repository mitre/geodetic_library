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
 * Carries out spherical inverse computations, returning course and distance
 */

void sphereInverse(LLPoint org, LLPoint dest, double* crs, double* dist,
                   double* userROC, double approxZero)
{

    if (crs != NULL)
        *crs = sphereInvCrs(org, dest, approxZero);
    if (dist != NULL)
        *dist = sphereInvDist(org, dest, userROC);

}

/*******************************************************************************
 * Carries out inverse computation but returns only the distance.  Slightly
 * faster than sphereInverse if you don't need to know the course.
 */

double sphereInvDist(LLPoint org, LLPoint dest, double* userROC)
{

    double dist = 0.0;
    double arg1, arg2;
#ifdef USE_BEST_FIT_ROC

    double sphereRad =
    lookUpROC(0.5*fabs(org.latitude+dest.latitude));
#else

    double sphereRad = SPHERE_RADIUS;
#endif

    /* If ROC is supplied, then use it to override best fit and default */
    if (userROC != NULL)
    {
        sphereRad = *userROC;
    }

    arg1 = sin(0.5 * (org.latitude - dest.latitude));
    arg1 = arg1 * arg1;
    arg2 = sin(0.5 * (org.longitude - dest.longitude));
    arg2 = arg2 * arg2;
    arg2 = cos(org.latitude) * cos(dest.latitude) * arg2;

    dist = sqrt(arg1 + arg2);
    dist = asin(dist);
    dist = 2.0 * sphereRad * dist;

    //    printf("InvDist sphereRad = %20.10f\n",sphereRad);

    return dist;

}

/*******************************************************************************/
/*
 * Carries out inverse computation but returns only the course.  Slightly
 * faster than sphereInverse if you don't need to know the distance.
 *
 */

double sphereInvCrs(LLPoint org, LLPoint dest, double approxZero)
{

    double course, yarg, xarg;

	yarg = sin(org.longitude - dest.longitude) * cos(dest.latitude);
	xarg = cos(org.latitude) * sin(dest.latitude) - (sin(org.latitude)
			* cos(dest.latitude) * cos(org.longitude - dest.longitude));
	//	  course = modpos(atan2(yarg,xarg),2*M_PI);
	course = M_2PI - modpos(atan2(yarg, xarg), M_2PI);
    course = modpos(course, M_2PI);

    return course;

}

/*******************************************************************************/
/* This function returns a destination point given a starting point, course, and
 * distance.  It is called by the initArcIntx function.
 */
ErrorSet sphereDirect(LLPoint org, double course, double dist,
                           LLPoint* dest)
{

    double beta = dist / SPHERE_RADIUS;
    double newLat, newLon, dLon;
    ErrorSet err = 0;
    int smallBetaForm = 0;
    double betaMax;
    double sinBeta;
    double angleToPole, angleToGo;
    double complementA;
    double sphereHalfCircum = M_PI * SPHERE_RADIUS;
    double sphereCircum = M_2PI * SPHERE_RADIUS;
    double dlat;
    LLPoint npDest;

    /* Assign local storage if optional pointers are not provided */
    if (NULL == dest) dest = &npDest;

    if(dist < 0.0){
    	/* Make dist positive */
    	dist = fabs(dist);
    	course = course + M_PI;
    }
    dist = fmod(dist, sphereCircum); /* Do not wrap around earth */
    if (dist > sphereHalfCircum)
    {
        // Need to go farther than half way around globe, so turn around and
        // go shorter distance along reciprocal course
        course = course + M_PI;
        dist = sphereCircum - dist;
    }
    beta = dist/SPHERE_RADIUS;

    course = modpos(course, M_2PI);

	/* Check for start at pole */
	if(fabs(org.latitude) == M_PI_2){
		dlat = (org.latitude > 0) ? -beta : beta;
		newLat = org.latitude + dlat;
		newLon = (org.latitude > 0) ? modlon((org.longitude + M_PI) - course): modlon((org.longitude) + course);
	}
	else if (course == 0.0)
    {
        /* Due north path.  May go over pole, so check how far to project beyond */
        angleToPole = M_PI_2 - org.latitude;
        angleToGo = beta;
        if (angleToGo <= angleToPole)
        {
            newLat = org.latitude + angleToGo;
            newLon = org.longitude;
        }
        else
        {
            newLat = M_PI_2 - (angleToGo - angleToPole);
            newLon = org.longitude + M_PI;
        }
    }
    else if (course == M_PI)
    {
        /* Due south path.  May go over pole, so check how far to project beyond */
        angleToPole = fabs(-M_PI_2 - org.latitude);
        angleToGo = beta;
        if (angleToGo <= angleToPole)
        {
            newLat = org.latitude - angleToGo;
            newLon = org.longitude;
        }
        else
        {
            newLat = -M_PI_2 + (angleToGo - angleToPole);
            newLon = org.longitude + M_PI;
        }
    }
    else
    {
        sinBeta = sin(beta);
        newLat = asin(cos(beta) * sin(org.latitude) + sinBeta * cos(
                org.latitude) * cos(course));

        /* betaMax is distance equivalent to 90 degree change in latitude.
         * If the given distance is less than this, then we can use the formula
         * for dLon directly.
         * If the given distance is greater than this, then we have to find the supplement
         * of the computed angle.
         * This is because asin is not singularly valued, and will return [-Pi/2,Pi/2] for
         * arguments in the range [-1,1].  Thus, changes in longitude greater than Pi/2
         * cannot be solved directly using this formula.
         */
        complementA = asin(sin(org.latitude) * sin(course));
        betaMax = acos(tan(complementA) / tan(course));
        /* The following formula is accurate for small angles, but will always return
         * -Pi/2 <= dLon <= Pi/2. Checking the required distance against the 90 degree
         * distance (betaMax) will tell us when to use the supplement to the
         * returned angle. */
        dLon = asin(sin(course) * sinBeta / cos(newLat));
        if (beta > fabs(betaMax))
        {
            if (dLon > 0.0)
                dLon = M_PI - dLon;
            else
                dLon = -M_PI - dLon;
        }
        newLon = org.longitude + dLon;
        smallBetaForm = 1;
    }

    newLon = modlon(newLon);

    dest->latitude = newLat;
    dest->longitude = newLon;

    return err;

}

/*******************************************************************************/
/*  Computes intersection points of two small circles on sphere
 *  Input points must be in geocentric coordinates
 *  Outputs will also be in geocentric coordinates
 *
 *  Tol represents the neighborhood around the intersection points. If the circles
 *  are separated by more than tol nmi, then no intersection will be found.  If the circles
 *  overlap by less than tol, then one tangent intersection will be returned.
 *
 */
ErrorSet initArcIntx(LLPoint center1, double r1, LLPoint center2, double r2,
                        LLPointPair intx, int* n, double* bestFitROC,
                        double tol)
{

    ErrorSet err = 0;

    Vector temp;
    Vector norm1;
    Vector norm2;
    Vector xLineDir;
    Vector x0 = { 0.0, 0.0, 0.0 };

    double n1_dot_n2;
    double D;

    double t, lx0, delta, crs, dist, midptErr, insideErr;
    double nn;
    double conditionNumber = 0.0;

    /* For debugging MATLAB calls */

#ifdef USE_BEST_FIT_ROC
    double sphereRad =
    lookUpROC(0.5*fabs(center1.latitude+center2.latitude));
#else
    double sphereRad = SPHERE_RADIUS;
#endif

    double halfCircumference = M_PI * sphereRad;

    double RHS[2];
    double C[2];

#ifdef USE_BEST_FIT_ROC

    if (bestFitROC != NULL)
    {
        *bestFitROC = sphereRad;
    }
#endif

	int npN;
	LLPointPair npIntx;

    /* Assign local storage if optional pointers are not provided */
    if (NULL == n) n = &npN;
    if (NULL == intx) intx = npIntx;

    if ((r1 >= halfCircumference) || (r2 >= halfCircumference))
    {
        err |= RADIUS_OUT_OF_RANGE_ERR;
        return err;
    }

    /* If radius is large enough, then the center point is farther
     * from the circle than its antipode.  In this case,
     * the checks for non-intersecting arcs are invalid.
     * Here we replace the center point with its antipode and recalculate the
     * radius to give the same circle. */
    if (r1 > halfCircumference / 2.0)
    {
        center1.latitude = -center1.latitude;
        center1.longitude = modlon(center1.longitude + M_PI);
        r1 = halfCircumference - r1;
    }
    if (r2 > halfCircumference / 2.0)
    {
        center2.latitude = -center2.latitude;
        center2.longitude = modlon(center2.longitude + M_PI);
        r2 = halfCircumference - r2;
    }

    norm1 = mapToUnitSphere(center1); // unit vector
    norm2 = mapToUnitSphere(center2); // unit vector

    xLineDir = cross(norm1, norm2); // v in documentation
    //    printf("xLineDir: (%20.15f, %20.15f, %20.15f)\n",xLineDir.x,xLineDir.y,xLineDir.z);
    n1_dot_n2 = dot(norm1, norm2); // n1.n2 in the documentation

    sphereInverse(center1, center2, &crs, &dist, &sphereRad, tol);
    //    printf("SphereInverse:\n");
    //    printf("  crs = %20.15f\n",crs);
    //    printf(" dist = %20.15f\n",crs);

    conditionNumber = (1 + n1_dot_n2) / (1 - n1_dot_n2);
    if (conditionNumber > 1.0e5)
    {
        /* Matrix is ill-conditioned and may lead to incorrect results,
         * so use spherical trig solution */
        double numerator = cos(r2 / sphereRad) - cos(dist / sphereRad) * cos(r1
                / sphereRad);
        double denominator = sin(dist / sphereRad) * sin(r1 / sphereRad);
        double arg = numerator / denominator;

        if (arg <= -1.0)
        { /* Protects against one case where roundoff leads to non-physical solution */
            *n = 1;
            err |= sphereDirect(center1, crs + M_PI, r1, &(intx[0]));
        }
        else if (arg >= 1.0)
        { /* Protects against other case of non-physical solution */
            *n = 1;
            err |= sphereDirect(center1, crs, r1, &(intx[0]));
        }
        else
        {
            double rho = acos(arg);
            *n = 2;
            err |= sphereDirect(center1, crs + rho, r1, &(intx[0]));
            err |= sphereDirect(center1, crs - rho, r1, &(intx[1]));
        }

        return err;

    }

    nn = pow(n1_dot_n2, 2.0);
    nn = 1.0 - nn;
    D = ((double) 1.0) / nn;

    //    D = 1.0 / (1.0 - n1_dot_n2 * n1_dot_n2); // Inverse of determinant of matrix to solve for C1, C2
    //    printf("    D: %20.15f\n",D);
    RHS[0] = sphereRad * cos(r1 / sphereRad);
    RHS[1] = sphereRad * cos(r2 / sphereRad);
    /* Doing these calcs in stages gives better agreement among compilers */
    //    C[0] = D * (RHS[0] - n1_dot_n2 * RHS[1]); // c1 in documentation
    nn = n1_dot_n2 * RHS[1];
    nn = RHS[0] - nn;
    C[0] = D * nn; // c1 in documentation
    //    C[1] = D * (RHS[1] - n1_dot_n2 * RHS[0]); // c2 in documentation
    nn = n1_dot_n2 * RHS[0];
    nn = RHS[1] - nn;
    C[1] = D * nn; // c2 in documentation

    /* Distance between circles along line joining centers */
    /* This is the dist between circles if neither center is inside
     * the other circle.   */
    midptErr = dist - (r1 + r2);

    //    printf(" midptErr = %20.15f\n",midptErr);

    /* if one circle inside other, distance between them
     * depends on which is bigger */
    if (r2 >= r1)
    {
        insideErr = r2 - r1 - dist;
    }
    else
    {
        insideErr = r1 - r2 - dist;
    }
    //    printf(" insideErr = %20.15f\n",insideErr);

    if ((midptErr > tol) || (insideErr > tol)) //TODO Check for tol in second clause

    {
        /* Circles don't intersect (first clause => widely separated;
         * second clause => smaller circle contained within larger circle) */
        *n = 0;
        err |= NO_INTERSECTION_ERR;
        return err;
    }
    else if ((fabs(midptErr) <= tol) || (fabs(insideErr) <= tol))
    {
        /* Circles just barely touch at one tangent intersection */
        *n = 1;
    }
    else
    {
        *n = 2;
    }

    scalarMultiply(&norm1, C[0]);
    scalarMultiply(&norm2, C[1]);
    //    printf("scaled norm1: (%20.15f, %20.15f, %20.15f)\n",norm1.x,norm1.y,norm1.z);
    //    printf("scaled norm2: (%20.15f, %20.15f, %20.15f)\n",norm2.x,norm2.y,norm2.z);

    x0 = vectorAdd(norm1, norm2); /* Vector to center of chord common to both circles */
    //    printf("x0: (%20.15f, %20.15f, %20.15f)\n",x0.x,x0.y,x0.z);
    normalize(&xLineDir); /* make intersection vector have unit length */

    lx0 = sqrt(dot(x0, x0)); //magnitude (length) of x0
    //    printf("lx0: %20.15f\n",lx0);
    delta = sphereRad - lx0;

    if (*n == 2) //(delta > tol)

    {
        if (delta < 0.0)
        {
            /* this should never happen if *n == 2 */
            /* Try spherical triangles solution */

            err |= UNEXPECTED_ERR;
            //            printf("initArcIntx: delta < 0.0\n");
            //            printf(" ---  midptErr = %e\n", midptErr);
            //            printf(" --- insideErr = %e\n", insideErr);
            return err;
        }
        t = sqrt(delta) * sqrt(sphereRad + lx0);
        if (t < tol)
        {
            *n = 1;
        }
        temp.x = x0.x + t * xLineDir.x;
        temp.y = x0.y + t * xLineDir.y;
        temp.z = x0.z + t * xLineDir.z;
        /* Will return geocentric coordinates */
        intx[0] = mapVectorToSphere(temp);
        temp.x = x0.x - t * xLineDir.x;
        temp.y = x0.y - t * xLineDir.y;
        temp.z = x0.z - t * xLineDir.z;
        /* Will return geocentric coordinates */
        intx[1] = mapVectorToSphere(temp);
    }
    else if (*n == 1)
    {
        /* Tangent case */
        if ((dist > r1) && (dist > r2))
        { /* intersection lies on line between circle centers */
            err |= sphereDirect(center1, crs, r1 + midptErr / 2.0, &(intx[0]));
        }
        else
        {
            /* One circle is inside another, but not concentric */
            if (r1 > r2) /* implies r1 > dist > r2, or circle2 is inside circle 1*/
            {
                err |= sphereDirect(center1, crs, r1 - insideErr / 2.0,
                        &(intx[0]));
            }
            else /* implies r2 < dist < r1, or circle 1 is inside circle 2 */
            {
                err |= sphereDirect(center1, crs + M_PI, r1 + insideErr / 2.0,
                        &(intx[0]));
            }
        }
    }

    return err;

}
} //namespace
