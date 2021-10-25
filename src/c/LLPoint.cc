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

/*  \file LLPoint.c
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

using namespace geolib_idealab;
/*******************************************************************************
 * Allocate and initialize a new LLPoint structure from given input parameters
 */

ErrorSet geolib_idealab::createPt(LLPoint* llpoint, double lat, double lon)
{

    ErrorSet err = 0;
    LLPoint npLlpoint;
    if(NULL == llpoint) llpoint = &npLlpoint;

    llpoint->latitude = lat;
    llpoint->longitude = lon;

    return err;

}

int geolib_idealab::ptsAreSame(LLPoint p1, LLPoint p2, double tol)
{
    double approxDist = 0.0;
//    smallDistInverse(p1, p2, NULL, &approxDist);
    inverse(p1,p2,NULL,NULL,&approxDist,1e-20);
    return (int) (approxDist <= tol);
}

ErrorSet geolib_idealab::ptsAreAntipodal(LLPoint p1, LLPoint p2)
{
    ErrorSet err = 0;
    double antipodalTest = sin(p1.latitude) * sin(p2.latitude) + cos(
            p1.latitude) * cos(p2.latitude) * cos(p1.longitude - p2.longitude);

    if (antipodalTest < cos(ANTIPODAL_TOL))

    {
        //        printf("antipodalTest = %.15f; ANTIPODAL_TOL = %.15f\n",antipodalTest,cos(ANTIPODAL_TOL));
        //        printf("P1 = (%15f, %15f)\n",p1.latitude*180.0/M_PI, p1.longitude*180.0/M_PI);
        //        printf("P2 = (%15f, %15f)\n",p2.latitude*180.0/M_PI, p2.longitude*180.0/M_PI);

        err = ANTIPODAL_POINTS_ERR;
    }

    return err;
}

/*
 * Returns an integer that indicates whether a testPt is at a pole and which one.
 * Valid returns:
 *      -1 -> South Pole
 *       0 -> NOT at a pole
 *       1 -> North Pole
 * */
int geolib_idealab::ptIsAtPole(LLPoint testPt, ErrorSet* err, double tol, double eps)
{
    double distFromNPole, distFromSPole;
    ErrorSet newErr = 0;
    LLPoint northPolePt, southPolePt;
    northPolePt.latitude = M_PI_2;
    northPolePt.longitude = 0;
    southPolePt.latitude = -M_PI_2;
    southPolePt.longitude = 0;

    /**\section Algorithm Algorithm Description
     *  <ol>
     *  <li> Get distance from North Pole using invDist, distFromNPole
     */
    newErr |= invDist(northPolePt, testPt, &distFromNPole, eps);
    if (newErr != SUCCESS && newErr != ANTIPODAL_POINTS_ERR) //filter out the ANTIPODAL_POINTS_ERR status code
    {
        *err |= newErr;
        return 0;
    }
    /** <ul><li>if |distFromNPole| < tol, return 1 </ul> */
    if (fabs(distFromNPole) < tol)
    {
        //Yes, at north pole.
        return 1;
    }

    /** <li> Get distance from South Pole using invDist, distFromSPole */
    newErr |= invDist(southPolePt, testPt, &distFromSPole, eps);
    if (newErr != SUCCESS && newErr != ANTIPODAL_POINTS_ERR) //filter out the ANTIPODAL_POINTS_ERR status code
    {
        *err |= newErr;
        return 0;
    }
    /** <ul><li>if |distFromSPole| < tol, return -1 </ul> */
    if (fabs(distFromSPole) < tol)
    {
        //Yes, at south pole.
        return -1;
    }

    //No, not at pole.
    /** <li>otherwise, return 0 </ol> */
    return 0;

}

LLPointSet geolib_idealab::createPtSet()
{
    LLPointSet set = { 0, NULL, LLPOINTSET_ARRAY_INCREMENT};
    set.elements = (LLPoint**) realloc(set.elements, LLPOINTSET_ARRAY_INCREMENT * sizeof(LLPoint*));
    return set;
}

static ErrorSet growPtSetByN(LLPointSet* set, int n)
{
    ErrorSet err = 0;
    int d = 0;
    int temp;

    if (NULL == set)
    {
        err |= NO_MEMORY_ALLOCATED_ERR;
    }
    else if (n < 0)
    {
        err |= MALLOC_ERR; //TODO: Should be different error code
    }
    else
    {
        d = (set->length + n) - set->maxLength; /* min number of additional array elements needed */
        if (d > 0)
        {
            /* Array will become longer than maxLength.
             * Reallocate memory to store additional LLPoints */
            set->maxLength += (d / LLPOINTSET_ARRAY_INCREMENT + 1)*LLPOINTSET_ARRAY_INCREMENT;
            set->elements = (LLPoint**) realloc(set->elements, set->maxLength
                    * sizeof(LLPoint*));
        }

        if (NULL == set->elements) {
            err |= MALLOC_ERR;
        }
        else
            set->length += n;
    }
    return err;
}
/** Increase size of array by 1 element */
static ErrorSet growPtSet(LLPointSet* set)
{
    return growPtSetByN(set, 1);
}

/** Frees all of the dynamically allocated memory associated with a
 * LLPointSet structure.  Resulting set has length == 0 and elements pointer
 * set to NULL.
 * @param set Pointer to an LLPointSet structure
 * @returns Nothing
 */
void geolib_idealab::clearPtSet(LLPointSet* set)
{
    int i = 0;

    if ((NULL == set) || (0 == set->length) || (NULL == set->elements))
    {
        return;
    }

    for (i = 0; i < set->length; i++)
    {
    	if (NULL != set->elements[i])
        {
            // printf("Freeing i=%d, mem = %#x, type = %d\n", i,
            /* Deallocate memory for element */
    		free(set->elements[i]);
    		set->elements[i] = NULL;
        }
    }
    if (NULL != set->elements)
    {
        /* Deallocate memory for array */
        free(set->elements);
        set->elements = NULL;
    }
    set->length = 0;

}

ErrorSet geolib_idealab::addPtToPtSet(LLPointSet* set, LLPoint* point)
{
    int n, i;
    ErrorSet err = 0;
    LLPoint* newPoint;
    for (i=0;i<set->length;i++) {
    	if ((point->latitude == set->elements[i]->latitude) && (point->longitude == set->elements[i]->longitude)) {
    		return err;
    	}
    }
    err |= growPtSet(set);
    if (err == 0)
    {
        n = set->length;
        if (NULL != point)
        {
            newPoint = static_cast<LLPoint *>(malloc(sizeof(LLPoint)));
            err |= geolib_idealab::createPt(newPoint, point->latitude, point->longitude);
        	set->elements[n - 1] = newPoint;
        }
        else
        {
            printf("\nNULL element supplied\n");
            err |= MALLOC_ERR;
        }
    }
    return err;

}

