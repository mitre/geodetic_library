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


/* Exit codes to trigger special action */
enum {
    NORMAL_EXIT = 0,
    POINT_IS_ON_BOUNDARY = 1,
    NEW_LINE_NEEDED = 2,
    TANGENT_ARC = 3,
    COMMON_ARC = 4
};


static ErrorSet arcShapeIntx(Arc a, Shape s, LLPointPair intx, int* n,
                                     int* exitCode, double tol, double eps);

static int ptIsOnShape(Shape s, LLPoint p, ErrorSet* err, double tol,
                          double eps);

/** Increase size of array by n elements
 *
 */
static ErrorSet growBndryArrayByN(Boundary* b, int n)
{
    ErrorSet err = 0;
    int d = 0;

    if (NULL == b)
    {
        err |= NO_MEMORY_ALLOCATED_ERR;
    }
    else if (n < 0)
    {
        err |= MALLOC_ERR; //TODO: Should be different error code
    }
    else
    {
        d = (b->length + n) - b->maxLength; /* min number of additional array elements needed */
        if (d > 0)
        {
            /* Array will become longer than maxLength.
             * Reallocate memory to store additional shapes */
            b->maxLength += (d / BOUNDARY_ARRAY_INCREMENT + 1)*BOUNDARY_ARRAY_INCREMENT;
            b->elements = (Shape*) realloc(b->elements, b->maxLength
                    * sizeof(Shape));
        }

        if (NULL == b->elements)
            err |= MALLOC_ERR;
        else
            b->length += n;
    }
    return err;

}
/** Increase size of array by 1 element */
static ErrorSet growBndryArray(Boundary* b)
{
    return growBndryArrayByN(b, 1);
}

/** Frees all of the dynamically allocated memory associated with a
 * Boundary structure.  Resulting boundary has length == 0 and elements pointer
 * set to NULL.
 * @param b Pointer to a Boundary structure
 * @returns Nothing
 */
void clearBndry(Boundary* b)
{
    int i = 0;

//    return;

    if ((NULL == b) || (0 == b->length) || (NULL == b->elements))
    {
        return;
    }

    for (i = 0; i < b->length; i++)
    {
        if (NULL != b->elements[i].this_shape)
        {
            // printf("Freeing i=%d, mem = %#x, type = %d\n", i,
            //        (int) b->elements[i].this, (int) b->elements[i].type);
            /* Deallocate memory for boundary element */
            free(b->elements[i].this_shape);
            b->elements[i].this_shape = NULL;
        }
    }
    if (NULL != b->elements)
    {
        /* Deallocate memory for Shape* array */
        free(b->elements);
        b->elements = NULL;
    }
    b->length = 0;

}

static ErrorSet appendPtToArray(LLPoint** parray, int* n, LLPoint p)
{
    LLPoint* newArray = NULL;
    LLPoint* oldArray = NULL;
    int i = 0;
    ErrorSet err = 0;

    /* make sure element count is an int */
    if (n == NULL)
    {
        return NO_MEMORY_ALLOCATED_ERR;
    }

    /* Allocate space for array with new element */
    if ((newArray = calloc(sizeof(LLPoint), *n + 1)) == NULL)
    {
        return MALLOC_ERR;
    }

    /* Add given point to end of new array */
    newArray[*n] = p;

    /* Copy over existing points */
    for (i = 0; i < *n; i++)
    {
        newArray[i] = (*parray)[i];
    }

    /* Reassign array address and free old memory */
    oldArray = *parray;
    *parray = newArray;
    if (oldArray != NULL)
    {
        free(oldArray);
    }

    /* Increment count (will update count in calling function) */
    *n = *n + 1;

    return err;

}

static void* dupBndryElement(const void* el, ShapeType type)
{

    void *newShape = NULL;

    if (NULL != el)
    {
        switch (type) {
        case GEODESIC:
            newShape = malloc(sizeof(Geodesic));
            *((Geodesic*) newShape) = *((Geodesic*) el);
            break;
        case ARC:
            newShape = malloc(sizeof(Arc));
            *((Arc*) newShape) = *((Arc*) el);
            break;
        case LOCUS:
            newShape = malloc(sizeof(Locus));
            *((Locus*) newShape) = *((Locus*) el);
            break;
        case SPIRAL:
            newShape = malloc(sizeof(Spiral));
            *((Spiral*) newShape) = *((Spiral*) el);
            break;
        case LLPOINT:
            newShape = malloc(sizeof(LLPoint));
            *((LLPoint*) newShape) = *((LLPoint*) el);
            break;
        case LOCUS2NDORDER:
            newShape = malloc(sizeof(Locus2ndOrder));
            *((Locus2ndOrder*) newShape) = *((Locus2ndOrder*) el);
            break;
        default:
        	break;
        }

    }
    return newShape;

}

/**
 * Name: sortPtsByAz
 *
 * Use this to sort array of points on an arc.  Array is
 * sorted in order of the arc direction with clockwise
 * starting at azimuth 0 radians
 * and counterclockwise stopping at azimuth 0 radians
 *
 * Currently this function does not account for modular arithmetic
 * so two points at azimuths x and (x + 2pi) will not be
 * treated as the same point.
 *
 * This function assumes that the azimuths have been normalized modulo 2pi,
 * i.e. azimuths in the range [0,2pi),
 * prior to calling the function otherwise unexpected results may occur.
 *
 * Input:
 *      points - array points to be sorted
 *      pointsLength - int - the number of points in the array
 *      azVals - array of azimuths [radians] associated with each point
 *      respectively in the range [0, 2pi)
 *      dir - the arc direction to sort the points in -1 = clockwise, 1 = counterclockwise
 * Output:
 * return - ErrorSet - error code 0x200 if and only if arguments; points or azVals is NULL.
 *              Otherwise return error code 0.
 */
static ErrorSet sortPtsByAz(LLPoint* points, int pointsLength, double* azVals, int dir)
{
    ErrorSet err = 0;
    LLPoint ip;
    double iaz;
    int i = 0, j = 0;

    if (points == NULL || azVals == NULL)
    {
        err |= NO_MEMORY_ALLOCATED_ERR;
        return err;
    }

    for (i = 1; i < pointsLength; i++)
    {
        if (&points[i] != NULL)
        {
            ip = points[i];
            iaz = azVals[i];
            j = i;
            while ((j > 0) && (dir == CLOCKWISE ? (azVals[j - 1] > iaz) : (azVals[j - 1]
                    < iaz)))
            {
                points[j] = points[j - 1];
                azVals[j] = azVals[j - 1];
                j = j - 1;
            }
            points[j] = ip;
            azVals[j] = iaz;
        }
    }

    return err;
}

ErrorSet addElementToBndry(Boundary* b, const void* element, ShapeType type)
{
    int n;
    ErrorSet err = 0;

    err |= growBndryArray(b);

    if (err == 0)
    {
        n = b->length;

        if (NULL != element)
        {
            //printf("Adding element of type ");
            switch (type) {
            case ARC:
                b->elements[n - 1].this_shape = (Arc*) dupBndryElement((Arc*) element,
                        type);
                //printf("ARC");
                break;
            case GEODESIC:
                b->elements[n - 1].this_shape = (Geodesic*) dupBndryElement(
                        (Geodesic*) element, type);
                //printf("GEODESIC");
                break;
            case LOCUS:
                b->elements[n - 1].this_shape = (Locus*) dupBndryElement((Locus*) element,
                        type);
                //printf("LOCUS");
                break;
            case SPIRAL:
                b->elements[n - 1].this_shape = (Spiral*) dupBndryElement(
                        (Spiral*) element, type);
                //printf("SPIRAL");
                break;
            case LLPOINT:
                b->elements[n - 1].this_shape = (LLPoint*) dupBndryElement(
                        (LLPoint*) element, type);
                //printf("LLPOINT");
                break;
            case LOCUS2NDORDER:
                b->elements[n - 1].this_shape = (Locus2ndOrder*) dupBndryElement(
                        (Locus2ndOrder*) element, type);
                break;
            default:
                return INVALID_TYPE_ERR; /* error */
            }
            //printf(": Success!\n");
            b->elements[n - 1].type = type;
        }
        else
        {
            printf("NULL element supplied\n");
        }
        /* Allow return of error code, not yet implemented */
    }

    return err;

}

/** Initialize a Boundary struct and return it. */
Boundary createBndry()
{
    Boundary b = { 0, NULL, BOUNDARY_ARRAY_INCREMENT};
    b.elements = (Shape*) realloc(b.elements, BOUNDARY_ARRAY_INCREMENT * sizeof(Shape));
    return b;
}

ErrorSet addGeoToBndry(Boundary* b, const Geodesic* g)
{
    return addElementToBndry(b, g, GEODESIC);
}

ErrorSet addArcToBndry(Boundary* b, const Arc* a)
{
    return addElementToBndry(b, a, ARC);
}

ErrorSet addLocusToBndry(Boundary* b, const Locus* l)
{
    return addElementToBndry(b, l, LOCUS);
}

ErrorSet addLocus2ndOrderToBndry(Boundary* b, const Locus2ndOrder* l)
{
    return addElementToBndry(b, l, LOCUS2NDORDER);
}

ErrorSet addSpiralToBndry(Boundary* b, const Spiral* s)
{
    return addElementToBndry(b, s, SPIRAL);
}


/* Return list of non-tangent intersections
 * NOTE: Intersections are only found on forward direction of geodesic,
 * starting at point p following azimuth az.  This is sufficient for the
 * inside/out algorithm.
 * b = Boundary structure of interest
 * p = Start point of geodesic that will be intersected
 * az = start azimuth of geodesic
 * intx = array of pointers to intersection points
 * intxCount = integer number of intersections found
 * exitCode = integer code used to return intersection status and alert
 *            inside/out algorithm if vertex is hit or if given point
 *            lies on boundary.
 * tol = the usual convergence tolerance
 * eps = the usual forward/inverse convergence tolerance */
ErrorSet bndryGeoIntx(Boundary b, LLPoint p, double az,
                                 LLPoint** intx, int* intxCount, int* exitCode,
                                 double tol, double eps)
{
    Arc* thisArc = NULL;
    Geodesic* thisGeo = NULL;
    Locus* thisLocus = NULL;
    Spiral*   thisSpiral = NULL;
    Shape* thisShape = NULL;

    ShapeType thisType;

    LLPoint tmpIntx;
    LLPoint geoPt;
    LLPoint pEnd; /* End point of test line for use with ptIsOnGeo */
    LLPointPair tmpIntx2;
    Geodesic geo;
	LLPointSet pts = createPtSet();
    double crs31, dist13, crs32, dist23;

    ErrorSet err = 0;
    ErrorSet locErr = 0;

    int tangentIntersection = 0;
    int i = 0, j = 0, n = 0;

    if ((intx == NULL) || (intxCount == NULL) || (exitCode == NULL))
    {
        return NO_MEMORY_ALLOCATED_ERR;
    }

    *intxCount = 0;

    /* Exit code is used to signal special cases to calling function
     *  exitCode = POINT_IS_ON_BOUNDARY => Start point of geodesic is on boundary
     *  exitCode = NEW_LINE_NEEDED      => Intersection point is a tangent point or a
     *                                     shape end point
     * */
    *exitCode = NORMAL_EXIT;

    err |= direct(p, az, 10.0, &pEnd, eps);

    for (i = 0; i < b.length; i++)
    {
        thisShape = b.elements[i].this_shape;
        if (thisShape == NULL)
        {
            return SHAPE_NOT_DEFINED_ERR;
        }
        thisType = b.elements[i].type;
        locErr = 0;
        switch (thisType) {
        case ARC:
            thisArc = (Arc*) thisShape;
            /* Check for intersection with line */
            locErr |= geoArcIntx(p, az, thisArc->centerPoint,
                    thisArc->radius, tmpIntx2, &n, tol, eps);

            if (locErr > 0)
            {
                /* Following errors should not trigger a failure */
                if (locErr & LINE_TOO_FAR_FROM_ARC_ERR)
                    locErr = locErr - LINE_TOO_FAR_FROM_ARC_ERR;
                if (locErr & CONCENTRIC_CIRCLE_ERR)
                    locErr = locErr - CONCENTRIC_CIRCLE_ERR;
                err = locErr;
            }

            if (n == 1)
            {
                /* May be tangent intersection, set return flag */
                tangentIntersection = 1;
            }
            for (j = 0; j < n; j++)
            {
                /* geoArcIntx extends line in unpredictable ways to
                 * ensure accurate intersection is found.  Must check that intx
                 * is actually on both arc and semi-infinite line extending from p */
                if ((ptIsOnGeo(p, pEnd, tmpIntx2[j], SEMIINFINITE, &err,
                        tol, eps)) && (ptIsOnArc(thisArc->centerPoint,
                        thisArc->radius, thisArc->startAz, thisArc->endAz,
                        thisArc->dir, tmpIntx2[j], &err, tol, eps)))
                {
                    if (ptsAreSame(tmpIntx2[j], p, tol))
                    {
                        /* test point is on boundary */
                        /* Should also handle case where p is start/end point */
                        *exitCode = POINT_IS_ON_BOUNDARY;
                    }
                    else if (ptsAreSame(tmpIntx2[j],
                            thisArc->startPoint, tol) || ptsAreSame(
                            tmpIntx2[j], thisArc->endPoint, tol))
                    {
                        /* Line hits start/end point, set return flag to pick different line */
                        *exitCode = NEW_LINE_NEEDED;
                    }
                    else if (tangentIntersection > 0)
                    {
                        *exitCode = NEW_LINE_NEEDED;
                    }
                    /* intxPoint is on both line and arc */
                    appendPtToArray(intx, intxCount, tmpIntx2[j]);
                }
            }
            break;

        case GEODESIC:
            thisGeo = (Geodesic*) thisShape;
            locErr |= geoIntx(p, pEnd, SEMIINFINITE, &crs31,
                    &dist13, thisGeo->startPoint, thisGeo->endPoint, SEGMENT,
                    &crs32, &dist23, &tmpIntx, tol, eps);
            if (locErr == 0)
            {
                /* Intersection was found, store it in list */
                appendPtToArray(intx, intxCount, tmpIntx);

                if (dist13 <= tol)
                {
                    /* p is on geodesic */
                    *exitCode = POINT_IS_ON_BOUNDARY;
                }
                else if (ptsAreSame(tmpIntx, thisGeo->startPoint, tol)
                        || ptsAreSame(tmpIntx, thisGeo->endPoint, tol))
                {
                    /* test line passes through vertex */
                    *exitCode = NEW_LINE_NEEDED;
                }
            }
            break;

        case LOCUS:
            thisLocus = (Locus*) thisShape;
            /* KNOWN BUG: GeoLocusIntersect will miss all double intersections */
            if (ptIsOnLocus(*thisLocus, p, &geoPt, &err, tol, eps))
            {
              tmpIntx = p;
              locErr = err;
            }
            else
              locErr |= locusGeoIntx(p, pEnd, *thisLocus, &tmpIntx, tol, eps);
            if ((locErr == 0) && (ptIsOnGeo(p, pEnd, tmpIntx,
                    SEMIINFINITE, &err, tol, eps)))
            {
                /* Intersection was found, store it in list */
                appendPtToArray(intx, intxCount, tmpIntx);

                if (ptIsOnLocus(*thisLocus, p, &geoPt, &err, tol, eps))
                {
                    /* p is on locus */
                    *exitCode = POINT_IS_ON_BOUNDARY;
                }
                else if (ptsAreSame(tmpIntx, thisLocus->locusStart, tol)
                        || ptsAreSame(tmpIntx, thisLocus->locusEnd, tol))
                {
                    /* test line passes through vertex */
                    *exitCode = NEW_LINE_NEEDED;
                }

            }
            break;
        case SPIRAL:
            thisSpiral = (Spiral*) thisShape;
            locErr |= createGeo(&geo, p, pEnd, SEMIINFINITE, eps);
            /* Check for intersection with line */
            locErr |= spiralGeoIntx(*thisSpiral, geo, &pts, tol, eps);
            if (locErr > 0)
            {
                /* Following errors should not trigger a failure */
                if (locErr & LINE_TOO_FAR_FROM_ARC_ERR)
                    locErr = locErr - LINE_TOO_FAR_FROM_ARC_ERR;
                if (locErr & CONCENTRIC_CIRCLE_ERR)
                    locErr = locErr - CONCENTRIC_CIRCLE_ERR;
                if (locErr & NO_INTERSECTION_ERR)
                    locErr = locErr - NO_INTERSECTION_ERR;
                err = locErr;

            }

//            printf("Pointset length:  %i\n", pts.length);
            for (j = 0; j < pts.length; j++)
            {
            	LLPoint pt;
            	pt.latitude = pts.elements[j]->latitude;
            	pt.longitude = pts.elements[j]->longitude;

                /* geoArcIntx extends line in unpredictable ways to
                 * ensure accurate intersection is found.  Must check that intx
                 * is actually on both arc and semi-infinite line extending from p */
//            	double c1, c2, d;
//            	locErr |= inverse(p, pt, &c1, &c2, &d, eps);
//            	ptIsOnGeo(p, pEnd, pt, SEMIINFINITE, &err, tol, eps);
//            	ptIsOnSpiral(*thisSpiral, pt, tol, eps);
                if ((ptIsOnGeo(p, pEnd, pt, SEMIINFINITE, &locErr, tol, eps)) && (ptIsOnSpiral(*thisSpiral, pt, tol, eps)))
                {
                    if (ptsAreSame(pt, p, tol))
                    {
                        /* test point is on boundary */
                        /* Should also handle case where p is start/end point */
                        *exitCode = POINT_IS_ON_BOUNDARY;
                    }
                    else if (ptsAreSame(pt,
                            thisSpiral->startPoint, tol) || ptsAreSame(
                            pt, thisSpiral->endPoint, tol))
                    {
                        /* Line hits start/end point, set return flag to pick different line */
                        *exitCode = NEW_LINE_NEEDED;
                    }
                    else if (tangentIntersection > 0)
                    {
                        *exitCode = NEW_LINE_NEEDED;
                    }
                    /* intxPoint is on both line and spiral */
                    appendPtToArray(intx, intxCount, pt);
                }
            }
            break;
        case LLPOINT:
            //TODO do something
            break;
        default:
        	break;
        }
    }
    return err;

}

/* Return list of intersections
 * b = Boundary structure of interest
 * locus that will be intersected
 * intx = array of pointers to intersection points
 * intxCount = integer number of intersections found
 * tol = the usual convergence tolerance
 * eps = the usual forward/inverse convergence tolerance */
ErrorSet bndryLocusIntx(Boundary b, Locus loc, 
                                 LLPoint** intx, int* intxCount, 
                                 double tol, double eps)
{
    Arc* thisArc = NULL;
    Geodesic* thisGeo = NULL;
    Locus* thisLocus = NULL;
    Spiral* thisSpiral = NULL;
    //    Spiral*   thisSpiral = NULL;  /* Needed in future */
    Shape* thisShape = NULL;
    ShapeType thisType;

    LLPoint tmpIntx, geoPt;
    LLPointPair tmpIntx2;
    LLPointSet pts;

    ErrorSet err = 0;
    ErrorSet locErr = 0;

    int i = 0, j = 0, n = 0, noIntersection = 0;

    if ((intx == NULL) || (intxCount == NULL))
    {
        return NO_MEMORY_ALLOCATED_ERR;
    }

    *intxCount = 0;

    for (i = 0; i < b.length; i++)
    {
        thisShape = b.elements[i].this_shape;
        if (thisShape == NULL)
        {
            return SHAPE_NOT_DEFINED_ERR;
        }
        thisType = b.elements[i].type;
        switch (thisType) {
        case ARC:
            thisArc = (Arc*) thisShape;
            err |= locusArcIntx(loc, thisArc->centerPoint, thisArc->radius, tmpIntx2, &n, tol, eps);
            if (err == 0)
            for (j = 0; j < n; j++)
            {
                    /* Intersection was found, store it in list */
                    if (ptIsOnArc(thisArc->centerPoint, thisArc->radius, thisArc->startAz, thisArc->endAz, thisArc->dir, tmpIntx2[j], &err, tol, eps))
                      appendPtToArray(intx, intxCount, tmpIntx2[j]);
            }
            break;

        case GEODESIC:
            thisGeo = (Geodesic*) thisShape;
            locErr |= locusGeoIntx(thisGeo->startPoint, thisGeo->endPoint, loc, &tmpIntx, tol, eps);
            noIntersection = 0;
            if (locErr > 0)
            {
              /* Following error should not trigger a failure */
              if (locErr & NO_INTERSECTION_ERR)
              {
                locErr = locErr - NO_INTERSECTION_ERR;
                noIntersection = 1;
              }
              err = locErr;
            }
            if (err == 0 && noIntersection == 0)
            {
                /* Intersection was found, store it in list */
                if (ptIsOnGeo(thisGeo->startPoint, thisGeo->endPoint, tmpIntx, SEGMENT, &err, tol, eps) && ptIsOnLocus(loc, tmpIntx, &geoPt, &err, tol, eps) )
                  appendPtToArray(intx, intxCount, tmpIntx);

            }
            break;

        case LOCUS:
            thisLocus = (Locus*) thisShape;
            locErr |= locusIntx(loc, *thisLocus, &tmpIntx,
                    tol, eps);
            noIntersection = 0;
            if (locErr > 0)
            {
              /* Following error should not trigger a failure */
              if (locErr & NO_INTERSECTION_ERR)
              {
                locErr = locErr - NO_INTERSECTION_ERR;
                noIntersection = 1;
              }
              err = locErr; 
            }
            if (err == 0 && noIntersection == 0)
            {
                /* Intersection was found, store it in list */
                appendPtToArray(intx, intxCount, tmpIntx);
            }
            break;

        case SPIRAL:
            thisSpiral = (Spiral*) thisShape;
            pts = createPtSet();
            /* Check for intersection with line */
            locErr |= spiralLocusIntx(*thisSpiral, loc, &pts, tol, eps);

            if (locErr > 0)
            {
                /* Following errors should not trigger a failure */
                if (locErr & LINE_TOO_FAR_FROM_ARC_ERR)
                    locErr = locErr - LINE_TOO_FAR_FROM_ARC_ERR;
                if (locErr & CONCENTRIC_CIRCLE_ERR)
                    locErr = locErr - CONCENTRIC_CIRCLE_ERR;
                if (locErr & NO_INTERSECTION_ERR)
                    locErr = locErr - NO_INTERSECTION_ERR;
                err = locErr;
            }

            for (j = 0; j < pts.length; j++)
            {
            	LLPoint pt;
            	pt.latitude = pts.elements[j]->latitude;
            	pt.longitude = pts.elements[j]->longitude;
                /* geoArcIntx extends line in unpredictable ways to
                 * ensure accurate intersection is found.  Must check that intx
                 * is actually on both arc and semi-infinite line extending from p */

                if ((ptIsOnLocus(loc, pt, &geoPt, &locErr, tol, eps)) && (ptIsOnSpiral(*thisSpiral, pt, tol, eps)))
                {
                    appendPtToArray(intx, intxCount, pt);
                }
            }
            break;
        case LLPOINT:
            //TODO do something
            break;
        default:
        	break;
        }
    }

    return err;

}

/* Intersect an arc (or circle) with Boundary structure
 * b = Boundary structure of interest
 * a = Arc to be intersected with b
 * intx = array of pointers to intersection LLPoints
 * intxCount = count of pointers in intx
 * exitCode = integer code used to return intersection status (purpose TBD)
 * tol = convergence tolerance (in nmi)
 * eps = forward/inverse algorithm convergence tolerance */
ErrorSet bndryArcIntx(Boundary b, Arc a, LLPoint** intx,
                                 int* intxCount, double tol, double eps)
{

    Arc* thisArc = NULL;
    Geodesic* thisGeo = NULL;
    Locus* thisLocus = NULL;
    Spiral* thisSpiral = NULL;
    Shape* thisShape = NULL;
    ShapeType thisType;

    LLPointSet pts;
    LLPoint pt;

    LLPointPair tmpIntx;

    int i = 0;
    int j = 0;
    int n = 0;
    int exitCode;
    int noIntersection = 0;
    int concCircErr = 0;
    ErrorSet err = 0;
    ErrorSet locErr = 0;
    *intxCount = 0;

    if ((NULL == intx) || (NULL == intxCount))
    {
        /* CANNOT CONTINUE */
        err |= NO_MEMORY_ALLOCATED_ERR;
        return err;
    }

    for (i = 0; i < b.length; i++)
    {
        /* Find intersections of arc with each boundary element */
        thisShape = b.elements[i].this_shape;
        if (thisShape == NULL)
        {
          return SHAPE_NOT_DEFINED_ERR;
        }
        thisType = b.elements[i].type;
        n = 0;
        locErr = 0;
        switch (thisType) {
          case ARC:
              thisArc = (Arc*) thisShape; 
              err |= arcIntx(a.centerPoint, a.radius, thisArc->centerPoint, thisArc->radius, tmpIntx, &n, tol, eps);
              if (err == 0)
              for (j = 0; j < n; j++)
              {
                    /* Intersection was found, store it in list */
                    if (ptIsOnArc(thisArc->centerPoint, thisArc->radius, thisArc->startAz, thisArc->endAz, thisArc->dir, tmpIntx[j], &err, tol, eps) &&
                    ptIsOnArc(a.centerPoint, a.radius, a.startAz, a.endAz, a.dir, tmpIntx[j], &err, tol, eps))
                       appendPtToArray(intx, intxCount, tmpIntx[j]);
              }
              break;
          case GEODESIC:
              thisGeo = (Geodesic*) thisShape;
              err |= geoArcIntx(thisGeo->startPoint, thisGeo->startAz, a.centerPoint, a.radius, tmpIntx, &n, tol, eps);
              if (err == 0)
              for (j = 0; j < n; j++)
              {
                    /* Intersection was found, store it in list */
                    if (ptIsOnGeo(thisGeo->startPoint, thisGeo->endPoint, tmpIntx[j], thisGeo->lineType, &err, tol, eps) && 
                    ptIsOnArc(a.centerPoint, a.radius, a.startAz, a.endAz, a.dir, tmpIntx[j], &err, tol, eps))
                       appendPtToArray(intx, intxCount, tmpIntx[j]);
              }
              break;
          case LOCUS:
              thisLocus = (Locus*) thisShape;
              /* This routine checks that intersection pt is on Locus */
              err |= locusArcIntx(*thisLocus, a.centerPoint, a.radius, tmpIntx, &n, tol, eps);
              if (err == 0)
              for (j = 0; j < n; j++)
              {
                  /* Intersection was found, store it in list */
                  if (ptIsOnArc(a.centerPoint, a.radius, a.startAz, a.endAz, a.dir, tmpIntx[j], &err, tol, eps))
                    appendPtToArray(intx, intxCount, tmpIntx[j]);
              }
              break;
          case SPIRAL:
       	      thisSpiral = (Spiral*) thisShape;
              pts = createPtSet();
              locErr |= spiralArcIntx(*thisSpiral, a, &pts, tol, eps);
              noIntersection = 0;
              if (locErr > 0)
              {
                /* Following error should not trigger a failure */
                if (locErr & NO_INTERSECTION_ERR)
                {
                  locErr = locErr - NO_INTERSECTION_ERR;
                  if (pts.length == 0)
                    noIntersection = 1;
                }
                err = locErr; 
              }
              if (err == 0 && noIntersection == 0)
              {
                for (j = 0; j < pts.length; j++) {
        	  pt.latitude = pts.elements[j]->latitude;
        	  pt.longitude = pts.elements[j]->longitude;
        	  appendPtToArray(intx, intxCount, pt);
                }
              }
              break;
          case LLPOINT:
              //TODO something
              break;
          default:
          	break;
        }

        /* Set exitCode to handle special cases if concentric arc error
         * is raised, or if arc is tangent to boundary */

        if (err & CONCENTRIC_CIRCLE_ERR)
        {
            /* handle special case */
            /* Add end points of common arc section to array of return points */
            /* The intersection of concentric arcs can be:  one arc, two arcs, an arc and a point, or two points.  Calling arcsCoincide can determine which case it is.  The code below was designed for the case of the intersection being one arc or two points.  What should be returned in the other cases needs to be determined.  A flag may be needed to indicate what case it is. The code below is commented out pending further review.  RES 9/30/11. */
            /*if (ptIsOnArc(thisArc->centerPoint, thisArc->radius, thisArc->startAz, thisArc->endAz, thisArc->dir, a.startPoint, &err, tol, eps))
        	  appendPtToArray(intx, intxCount, a.startPoint);
            else if (ptIsOnArc(a.centerPoint, a.radius, a.startAz, a.endAz, a.dir, thisArc->startPoint, &err, tol, eps))
        	  appendPtToArray(intx, intxCount, thisArc->startPoint);

            if (ptIsOnArc(thisArc->centerPoint, thisArc->radius, thisArc->startAz, thisArc->endAz, thisArc->dir, a.endPoint, &err, tol, eps))
        	  appendPtToArray(intx, intxCount, a.endPoint);
            else if (ptIsOnArc(a.centerPoint, a.radius, a.startAz, a.endAz, a.dir, thisArc->endPoint, &err, tol, eps))
        	  appendPtToArray(intx, intxCount, thisArc->endPoint);*/
            err = err - CONCENTRIC_CIRCLE_ERR;
            concCircErr = 1;
        }

    }//for i

    if (concCircErr)
      err = err + CONCENTRIC_CIRCLE_ERR;

    if (*intxCount == 1)
    {
        /* Arc is either tangent to boundary or an error has occurred */
        /* Test for tangency.  Set error>0 if this test fails */
        /* Algorithm TBD */
        /* Set exitCode to indicate this case */
    }
    else if (*intxCount > 1)
    {
        /* Sort list of points according to azimuth from arc center
         * NOTE: Spherical azimuth is sufficient, as we're only using it for comparative
         * purposes */
        /* NOTE ALSO: This is an insertion sort algorithm */
        LLPoint ip;
        double* azVals;
        double iaz;
        j = 0;
        azVals = calloc(sizeof(double), *intxCount);
        for (i = 0; i < *intxCount; i++)
        {
            /* Store azimuth values */
            azVals[i] = sphereInvCrs(a.centerPoint, ip, eps);
        }
        for (i = 1; i < *intxCount; i++)
        {
            ip = (*intx)[i]; /* intx is a pointer to an array */
            iaz = azVals[i];
            j = i;
            while ((j > 0) && (azVals[j - 1] > iaz))
            {
                (*intx)[j] = (*intx)[j - 1];
                azVals[j] = azVals[j - 1];
                j = j - 1;
            }
            (*intx)[j] = ip;
            azVals[j] = iaz;

        }

        if (azVals)
            free(azVals);

    }

    return err;

}

//ErrorSet bndrySpiralIntx(Boundary b, Spiral s, LLPoint** intx,
//                                 int* intxCount, double tol, double eps)
//{
//
//    Arc* thisArc = NULL;
//    Geodesic* thisGeo = NULL;
//    Locus* thisLocus = NULL;
//    Spiral* thisSpiral = NULL;
//    Shape* thisShape = NULL;
//    ShapeType thisType;
//
//    LLPointSet pts;
//    LLPoint pt;
//
//    LLPointPair tmpIntx;
//
//    int i = 0;
//    int j = 0;
//    int n = 0;
//    int exitCode;
//    int noIntersection = 0;
//    int concCircErr = 0;
//    ErrorSet err = 0;
//    ErrorSet locErr = 0;
//    *intxCount = 0;
//
//    if ((NULL == intx) || (NULL == intxCount))
//    {
//        /* CANNOT CONTINUE */
//        err |= NO_MEMORY_ALLOCATED_ERR;
//        return err;
//    }
//
//    for (i = 0; i < b.length; i++)
//    {
//        /* Find intersections of arc with each boundary element */
//        thisShape = b.elements[i].this;
//        if (thisShape == NULL)
//        {
//          return SHAPE_NOT_DEFINED_ERR;
//        }
//        thisType = b.elements[i].type;
//        n = 0;
//        locErr = 0;
//        switch (thisType) {
//          case ARC:
//        	thisArc = (Arc*) thisShape;
//			pts = createPtSet();
//			locErr |= spiralArcIntx(s, *thisArc, &pts, tol, eps);
//			noIntersection = 0;
//			if (locErr > 0)
//			{
//			  /* Following error should not trigger a failure */
//			  if (locErr & NO_INTERSECTION_ERR)
//			  {
//				locErr = locErr - NO_INTERSECTION_ERR;
//				if (pts.length == 0)
//				  noIntersection = 1;
//			  }
//			  err = locErr;
//			}
//			if (err == 0 && noIntersection == 0)
//			{
//			  for (j = 0; j < pts.length; j++) {
//				pt.latitude = pts.elements[j]->latitude;
//				pt.longitude = pts.elements[j]->longitude;
//				appendPtToArray(intx, intxCount, pt);
//			  }
//			}
//			break;
////          case GEODESIC:
////              thisGeo = (Geodesic*) thisShape;
////              err |= spiralGeoIntx(s, *thisGeo, &pts, tol, eps);
////              if (locErr > 0) {
////				  /* Following errors should not trigger a failure */
////				  if (locErr & LINE_TOO_FAR_FROM_ARC_ERR)
////					  locErr = locErr - LINE_TOO_FAR_FROM_ARC_ERR;
////				  if (locErr & CONCENTRIC_CIRCLE_ERR)
////					  locErr = locErr - CONCENTRIC_CIRCLE_ERR;
////				  if (locErr & NO_INTERSECTION_ERR)
////					  locErr = locErr - NO_INTERSECTION_ERR;
////				  err = locErr;
////			  }
////
////			  for (j = 0; j < pts.length; j++)
////			  {
////				LLPoint pt;
////				pt.latitude = pts.elements[j]->latitude;
////				pt.longitude = pts.elements[j]->longitude;
////				  /* geoArcIntx extends line in unpredictable ways to
////				   * ensure accurate intersection is found.  Must check that intx
////				   * is actually on both arc and semi-infinite line extending from p */
////				  if ((ptIsOnGeo(p, pEnd, pt, SEMIINFINITE, &err, tol, eps)) && (ptIsOnSpiral(*thisSpiral, pt, tol, eps)))
////				  {
////					  if (ptsAreSame(pt, p, tol))
////					  {
////						  /* test point is on boundary */
////						  /* Should also handle case where p is start/end point */
////						  *exitCode = POINT_IS_ON_BOUNDARY;
////					  }
////					  else if (ptsAreSame(pt,
////							  thisSpiral->startPoint, tol) || ptsAreSame(
////							  pt, thisSpiral->endPoint, tol))
////					  {
////						  /* Line hits start/end point, set return flag to pick different line */
////						  *exitCode = NEW_LINE_NEEDED;
////					  }
////					  else if (tangentIntersection > 0)
////					  {
////						  *exitCode = NEW_LINE_NEEDED;
////					  }
////					  /* intxPoint is on both line and spiral */
////					  appendPtToArray(intx, intxCount, pt);
////				  }
////			  }
////			  break;
////          case LOCUS:
////              thisLocus = (Locus*) thisShape;
////              /* This routine checks that intersection pt is on Locus */
////              err |= locusArcIntx(*thisLocus, a.centerPoint, a.radius, tmpIntx, &n, tol, eps);
////              if (err == 0)
////              for (j = 0; j < n; j++)
////              {
////                  /* Intersection was found, store it in list */
////                  if (ptIsOnArc(a.centerPoint, a.radius, a.startAz, a.endAz, a.dir, tmpIntx[j], &err, tol, eps))
////                    appendPtToArray(intx, intxCount, tmpIntx[j]);
////              }
////              break;
////          case SPIRAL:
////       	      thisSpiral = (Spiral*) thisShape;
////              pts = createPtSet();
////              locErr |= spiralArcIntx(*thisSpiral, a, &pts, tol, eps);
////              noIntersection = 0;
////              if (locErr > 0)
////              {
////                /* Following error should not trigger a failure */
////                if (locErr & NO_INTERSECTION_ERR)
////                {
////                  locErr = locErr - NO_INTERSECTION_ERR;
////                  if (pts.length == 0)
////                    noIntersection = 1;
////                }
////                err = locErr;
////              }
////              if (err == 0 && noIntersection == 0)
////              {
////                for (j = 0; j < pts.length; j++) {
////        	  pt.latitude = pts.elements[j]->latitude;
////        	  pt.longitude = pts.elements[j]->longitude;
////        	  appendPtToArray(intx, intxCount, pt);
////                }
////              }
////              break;
//          case LLPOINT:
//              //TODO something
//              break;
//        }
//
//        /* Set exitCode to handle special cases if concentric arc error
//         * is raised, or if arc is tangent to boundary */
//
//        if (err & CONCENTRIC_CIRCLE_ERR)
//        {
//            /* handle special case */
//            /* Add end points of common arc section to array of return points */
//            /* The intersection of concentric arcs can be:  one arc, two arcs, an arc and a point, or two points.  Calling arcsCoincide can determine which case it is.  The code below was designed for the case of the intersection being one arc or two points.  What should be returned in the other cases needs to be determined.  A flag may be needed to indicate what case it is. The code below is commented out pending further review.  RES 9/30/11. */
//            /*if (ptIsOnArc(thisArc->centerPoint, thisArc->radius, thisArc->startAz, thisArc->endAz, thisArc->dir, a.startPoint, &err, tol, eps))
//        	  appendPtToArray(intx, intxCount, a.startPoint);
//            else if (ptIsOnArc(a.centerPoint, a.radius, a.startAz, a.endAz, a.dir, thisArc->startPoint, &err, tol, eps))
//        	  appendPtToArray(intx, intxCount, thisArc->startPoint);
//
//            if (ptIsOnArc(thisArc->centerPoint, thisArc->radius, thisArc->startAz, thisArc->endAz, thisArc->dir, a.endPoint, &err, tol, eps))
//        	  appendPtToArray(intx, intxCount, a.endPoint);
//            else if (ptIsOnArc(a.centerPoint, a.radius, a.startAz, a.endAz, a.dir, thisArc->endPoint, &err, tol, eps))
//        	  appendPtToArray(intx, intxCount, thisArc->endPoint);*/
//            err = err - CONCENTRIC_CIRCLE_ERR;
//            concCircErr = 1;
//        }
//
//    }//for i
//
//    if (concCircErr)
//      err = err + CONCENTRIC_CIRCLE_ERR;
//
//    if (*intxCount == 1)
//    {
//        /* Arc is either tangent to boundary or an error has occurred */
//        /* Test for tangency.  Set error>0 if this test fails */
//        /* Algorithm TBD */
//        /* Set exitCode to indicate this case */
//    }
//    else if (*intxCount > 1)
//    {
//        /* Sort list of points according to azimuth from arc center
//         * NOTE: Spherical azimuth is sufficient, as we're only using it for comparative
//         * purposes */
//        /* NOTE ALSO: This is an insertion sort algorithm */
//        LLPoint ip;
//        double* azVals;
//        double iaz;
//        j = 0;
//        azVals = calloc(sizeof(double), *intxCount);
//        for (i = 0; i < *intxCount; i++)
//        {
//            /* Store azimuth values */
//            azVals[i] = sphereInvCrs(a.centerPoint, ip, eps);
//        }
//        for (i = 1; i < *intxCount; i++)
//        {
//            ip = (*intx)[i]; /* intx is a pointer to an array */
//            iaz = azVals[i];
//            j = i;
//            while ((j > 0) && (azVals[j - 1] > iaz))
//            {
//                (*intx)[j] = (*intx)[j - 1];
//                azVals[j] = azVals[j - 1];
//                j = j - 1;
//            }
//            (*intx)[j] = ip;
//            azVals[j] = iaz;
//
//        }
//
//        if (azVals)
//            free(azVals);
//
//    }
//
//    return err;
//
//}

/*
 * Intersect a circle with a simple connected Boundary structure
 *
 * Return a new boundary that defines the area that is common to both the circle
 * and to the input boundary.
 *
 * A simple boundary implies one that contains no holes in the boundary area.  This also
 * means that the boundary area is not partitioned (contains sub areas) in any way.
 *
 *                    not partitioned           partioned
 *                       ---------              ---------
 *                       |       |              |   |   |
 *                       |       |              |   |   |
 *                       ---------              ---------
 *
 * A connected boundary implies that all the start/end points are "connected" to one
 * another.  Two shapes are connected if they share an identical point or if the
 * Vincenty distance between the two points is less than tol in nautical miles.
 *
 * The speed of this algorithm is directly proportional to the number of shapes that
 * comprise the boundary.  More boundary shapes implies more time to process.
 *
 * inputs:
 * 		boundary = Boundary structure of interest
 * 		circle = Circle to be intersected with boundary
 * 		createBndry = Pointer to the new boundary.
 * 		tol = distance criteria (in nmi)
 * 		eps = Vincenty forward/inverse algorithm convergence criteria
 */
ErrorSet bndryCircleIntx(Boundary boundary, Arc circle,
                                    Boundary* newBoundary, double tol,
                                    double eps)
{
    ErrorSet err = 0; //The error code, 0 if successful
    int i = 0, j = 0; //for-loop indices
    double circleExtent; //The signed arc extent of Arc circle
    Shape* segment = NULL; //The shape currently being processed
    int segmentCount; // Number of shapes in given boundary
    ShapeType segClass; //The shape type of the current shape being processed
    Arc* segmentArc = NULL; //Temp object to hold segment recast as an Arc
    Geodesic* segmentGeo = NULL; //Temp object to hold segment recast as an Geodesic
    Locus* segmentLocus = NULL; //Temp object to hold segment recast as an Locus
    LLPoint* segmentLLPoint = NULL; //Temp object to hold segment recast an an LLPoint
    LLPointPair intPoints; //Holds the intersection points for a particluar segment-cirlce intersection
    LLPoint tempIntPoint; //Holds the intersection point temporarily
    int n = -1; //The number of intersection points returned by a particular segment-circle intersection
    int nFnd = -1; //The number of intersection points found by an intersection method.  Not necessarily same as n.
    int coincide = 0; //1 = arc and circle coincide
    int nCoincide = 0; //The number of shapes that coincide between a circle and arc
    Shape coincideCommonShapes[2]; //Common shapes if arc and circle coincide
    int tangent = 0;
    LLPoint* discretePoints = NULL; //Array of pointers for discrete intersection points as we find them
    int nDiscretePoints = 0; //The number of discrete points found
    LLPoint* arcPoints; //Array of pointers for the intersection points that lie on the circle
    int nArcPoints = 0; //Number of intersection points that lie on the circle
    double* arcPointsAngle; //Array of pointers to the start true course from the arcPoint to the circle center point
    Arc testArc; //The test arc constructed from the arc points
    LLPoint testPoint; //The point on the test arc that lies along the angle bisector of the test arc
    Shape* newSegment = NULL; //A temporary memory location that holds the new shape to be created for the new boundary
    LLPoint newStartPoint; //The start point of the new segment (This is the geo start point for a new segment of type Locus)
    LLPoint newEndPoint; //The end point of the new segment  (This is the geo end point for a new segment of type Locus)
    int startPointIsInside = 0; //1 if start point is inside, 0 if point is outside
    int endPointIsInside = 0; //1 if end point is inside, 0 if point is outside
    int pointsOnly;
    LLPoint tempStartPoint;
    LLPoint tempEndPoint;
    double tempStartDistance;
    double tempEndDistance;
    double angleBtwn;
    double halfAngleBtwn;
    double angleBisectorCourse;
    LineType tempLineType;
    LLPoint tempPt; //Generic temporary point
    LLPoint dummyPtIsOnLocusPt;
    double dummyCourse;
    double arcExtent; //The signed arc extent of a shape of type Arc
    Geodesic testLine;
	LLPoint testArcMidPt;
	double crsSegCtrToNewStart;
	double crsSegCtrToNewEnd;
	double crsSegCtrToMidPt;
	double subtendedAngle;
	double tempCrs;

    // Check for nulls
    if ((NULL == &boundary) || (NULL == &circle) || (NULL == newBoundary))
    {
        /* CANNOT CONTINUE */
        err |= NO_MEMORY_ALLOCATED_ERR;
        return err;
    }

    //intPoints[0] = NULL;
    //intPoints[1] = NULL;

    // Get the angle extent of the circle
    err |= getArcExtent(circle.centerPoint, circle.radius,
            circle.startAz, circle.endAz, circle.dir, &circleExtent, tol, eps);
    if (err)
        return err;

    // Check if full circle
    if (fabs(circleExtent) < M_2PI)
    {
        /* CANNOT CONTINUE */
        err |= SUBTENDED_ANGLE_OUT_OF_RANGE_ERR;
        return err;
    }

    // Set the number of elements in given boundary
    segmentCount = boundary.length;

    /* Allocate memory
     * Each non-llpoint segment can have at most 2 intersections with circle.
     * Add an extra memory location for creating test arcs later.
     */
    arcPoints = calloc(2 * segmentCount + 1, sizeof(LLPoint));
    arcPointsAngle = calloc(2 * segmentCount + 1, sizeof(double));
    discretePoints = calloc(segmentCount + 1, sizeof(LLPoint));

    for (i = 0; i < segmentCount; i++)
    {
        //Get the shape from the boundary
        segment = boundary.elements[i].this_shape;

        //Get the shape type
        segClass = boundary.elements[i].type;

        //Set the value of n that won't trigger further analysis of current segment
        n = -1;

        newSegment = NULL;
        startPointIsInside = 0;
        endPointIsInside = 0;
        coincide = 0;
        tangent = 0;

        switch (segClass) {
        case LOCUS:

            //Cast shape as a Locus
            segmentLocus = (Locus*) segment;

            //Find intersections of segment with circle

            //First store the original lineType
            tempLineType = segmentLocus->lineType;

            //segmentLocus must be of lineType INFINITE to implement LocusArcIntersect properly
            segmentLocus->lineType = INFINITE;

            //Get the intersection points of the circle and locus
            err |= locusArcIntx(*segmentLocus, circle.centerPoint,
                    circle.radius, intPoints, &nFnd, tol, eps);
            //Set lineType back to original
            segmentLocus->lineType = tempLineType;

            //Check if found intersection points are on the shape
            n = nFnd;

            if ((n > 1) && !ptIsOnLocus(*segmentLocus, intPoints[1],
                    &tempPt, &err, tol, eps))
            {
                //If point is not on shape then do not count
                intPoints[1].latitude = 0;
                intPoints[1].longitude = 0;
                n--;
            }

            if ((n > 0) && !ptIsOnLocus(*segmentLocus, intPoints[0],
                    &tempPt, &err, tol, eps))
            {
                //If point is not on shape then do not count
                intPoints[0].latitude = 0;
                intPoints[0].longitude = 0;
                n--;
            }

            //Determine location of segment start/end points with respect to circle
            startPointIsInside = ptIsInsideArc(circle.centerPoint,
                    circle.radius, circle.startAz, circle.endAz, circle.dir,
                    segmentLocus->locusStart, &err, tol, eps);
            endPointIsInside = ptIsInsideArc(circle.centerPoint,
                    circle.radius, circle.startAz, circle.endAz, circle.dir,
                    segmentLocus->locusEnd, &err, tol, eps);

            break;

        case GEODESIC:

            //Cast shape as a Geodesic
            segmentGeo = (Geodesic*) segment;

            //Find intersections of segment with circle

            // Mask out error codes CONCENTRIC_CIRCLE_ERR.
            // In this case, CONCENTRIC_CIRCLE_ERR is an exit indicator only and not an actual error.
            // See bug 16362
            err |= getMaskedError(geoArcIntx(segmentGeo->startPoint,
                    segmentGeo->startAz, circle.centerPoint, circle.radius,
                    intPoints, &nFnd, tol, eps), getMask(0, 1, 0, 0, 0, 0, 0));

            //Check if found intersection points are on the shape
            n = nFnd;

            if (n > 1)
            {
                if (!ptIsOnGeo(segmentGeo->startPoint,
                        segmentGeo->endPoint, intPoints[1],
                        segmentGeo->lineType, &err, tol, eps))
                {
                    //If point is not on shape then do not count
                    intPoints[1].latitude = 0;
                    intPoints[1].longitude = 0;
                    n--;

                    //check if intersection point is a start/end point
                }
                else
                {
                    if (ptsAreSame(segmentGeo->startPoint,
                            intPoints[1], tol))
                    {
                        intPoints[1] = segmentGeo->startPoint;
                    }
                    else if (ptsAreSame(segmentGeo->endPoint,
                            intPoints[1], tol))
                    {
                        intPoints[1] = segmentGeo->endPoint;
                    }
                }
            }

            if (n > 0)
            {
                if (!ptIsOnGeo(segmentGeo->startPoint,
                        segmentGeo->endPoint, intPoints[0],
                        segmentGeo->lineType, &err, tol, eps))
                {
                    //If point is not on shape then do not count
                    intPoints[0].latitude = 0;
                    intPoints[0].longitude = 0;
                    n--;
                    //check if intersection point is a start/end point
                }
                else
                {
                    if (ptsAreSame(segmentGeo->startPoint,
                            intPoints[0], tol))
                    {
                        intPoints[0] = segmentGeo->startPoint;
                    }
                    else if (ptsAreSame(segmentGeo->endPoint,
                            intPoints[0], tol))
                    {
                        intPoints[0] = segmentGeo->endPoint;
                    }
                }
            }

            //Determine location of segment start/end points with respect to circle
            startPointIsInside = ptIsInsideArc(circle.centerPoint,
                    circle.radius, circle.startAz, circle.endAz, circle.dir,
                    segmentGeo->startPoint, &err, tol, eps);
            endPointIsInside = ptIsInsideArc(circle.centerPoint,
                    circle.radius, circle.startAz, circle.endAz, circle.dir,
                    segmentGeo->endPoint, &err, tol, eps);

            break;

        case ARC:

            //Cast shape as an Arc
            segmentArc = (Arc*) segment;

            //Check if arcs coincide
            coincide = arcsCoincide(circle, *segmentArc, coincideCommonShapes,
                    &nCoincide, &err, tol, eps);

            if (coincide)
            {
                /*
                 * If the arc is a "circle" and it coincides with the boundary then it cannot be connected to
                 * any other shapes and still comprise a simple boundary.  Thus the boundary and circle must
                 * be the same.  Continue processing if this is the case.
                 *
                 * If the arc is not a "circle" and it coincides with the boundary then any shape connected
                 * to the arc's start/end points will generate an intersection point to be processed later on,
                 * assuming the boundary is a simple and connected boundary.
                 * In this case, ignore the arc and continue processing the next shape in the boundary.
                 *
                 * Note - arcIntx should return zero intersections for concentric arcs which
                 *        implies zero intersections will be found for arcs that coincide.
                 */

                //Get the signed arc extent with respect to tol
                err |= getArcExtent(segmentArc->centerPoint,
                        segmentArc->radius, segmentArc->startAz,
                        segmentArc->endAz, segmentArc->dir, &arcExtent, tol,
                        eps);

                //if arc is not a "circle" then move onto next shape
                if (fabs(arcExtent) < M_2PI)
                {
                    continue;
                }

            }

            //Get the intersection points of the circle and arc

            // Mask out error codes CIRCLE_INSIDE_CIRCLE_ERR.
            // In this case, CIRCLE_INSIDE_CIRCLE_ERR is an exit indicator only and not an actual error.
            // See bug 16362

            err |= getMaskedError(arcIntx(segmentArc->centerPoint,
                    segmentArc->radius, circle.centerPoint, circle.radius,
                    intPoints, &nFnd, tol, eps), getMask(0, 0, 0, 0, 0, 0, 1));

            //Check if found intersection points are on the shape
            n = nFnd;

            if (n > 1)
            {
                if (!ptIsOnArc(segmentArc->centerPoint,
                        segmentArc->radius, segmentArc->startAz,
                        segmentArc->endAz, segmentArc->dir, intPoints[1], &err,
                        tol, eps))
                {
                    //If point is not on shape then do not count
                    intPoints[1].latitude = 0;
                    intPoints[1].longitude = 0;
                    n--;

                    //check if intersection point is a start/end point
                }
                else
                {
                    if (ptsAreSame(segmentArc->startPoint,
                            intPoints[1], tol))
                    {
                        intPoints[1] = segmentArc->startPoint;
                    }
                    else if (ptsAreSame(segmentArc->endPoint,
                            intPoints[1], tol))
                    {
                        intPoints[1] = segmentArc->endPoint;
                    }
                }
            }

            if (n > 0)
            {
                if (!ptIsOnArc(segmentArc->centerPoint,
                        segmentArc->radius, segmentArc->startAz,
                        segmentArc->endAz, segmentArc->dir, intPoints[0], &err,
                        tol, eps))
                {
                    //If point is not on shape then do not count
                    intPoints[0].latitude = 0;
                    intPoints[0].longitude = 0;
                    n--;
                    //check if intersection point is a start/end point
                }
                else
                {
                    if (ptsAreSame(segmentArc->startPoint,
                            intPoints[0], tol))
                    {
                        intPoints[0] = segmentArc->startPoint;
                    }
                    else if (ptsAreSame(segmentArc->endPoint,
                            intPoints[0], tol))
                    {
                        intPoints[0] = segmentArc->endPoint;
                    }
                }
            }

            //Determine location of segment start/end points with respect to circle
            startPointIsInside = ptIsInsideArc(circle.centerPoint,
                    circle.radius, circle.startAz, circle.endAz, circle.dir,
                    segmentArc->startPoint, &err, tol, eps);
            endPointIsInside = ptIsInsideArc(circle.centerPoint,
                    circle.radius, circle.startAz, circle.endAz, circle.dir,
                    segmentArc->endPoint, &err, tol, eps);

            break;

        case LLPOINT:

            //Cast shape as an LLPoint
            segmentLLPoint = (LLPoint*) segment;

            //Determine location of segment point with respect to circle
            startPointIsInside = ptIsInsideArc(circle.centerPoint,
                    circle.radius, circle.startAz, circle.endAz, circle.dir,
                    *segmentLLPoint, &err, tol, eps);

            break;

        default:

            // CANNOT CONTINUE //
            err |= SHAPE_NOT_DEFINED_ERR;
            return err;
        }
        //Process intersection points
        switch (n) {

        case 0://No intersection points found

            //Determine location of segment start/end points with respect to circle
            switch (segClass) {
            case LOCUS:

                //Check if segment is inside of circle
                if (startPointIsInside && endPointIsInside)
                {

                    newStartPoint = segmentLocus->locusStart;
                    newEndPoint = segmentLocus->locusEnd;

                }
                else
                {
                    //No intersections found and segment is outside circle.  Move on to next shape.
                    continue;
                }

                break;

            case GEODESIC:

                //Check if segment is inside of circle
                if (startPointIsInside && endPointIsInside)
                {

                    newStartPoint = segmentGeo->startPoint;
                    newEndPoint = segmentGeo->endPoint;

                }
                else
                {
                    //No intersections found and segment is outside circle.  Move on to next shape.
                    continue;
                }

                break;

            case ARC:

                //Check if segment is inside of circle
                if (startPointIsInside && endPointIsInside)
                {

                    newStartPoint = segmentArc->startPoint;
                    newEndPoint = segmentArc->endPoint;

                }
                else
                {
                    //No intersections found and segment is outside circle.  Move on to next shape.
                    continue;
                }

                break;

            case LLPOINT:

                if (startPointIsInside)
                {

                    intPoints[0] = *segmentLLPoint;
                    newStartPoint = *segmentLLPoint;
                }
                else
                {
                    //No intersections found and segment is outside circle.  Move on to next shape.
                    continue;
                }

                break;

            default:

                // CANNOT CONTINUE //
                err |= SHAPE_NOT_DEFINED_ERR;
                return err;
            }

            break;

        case 1: //One intersection point found

            //Set the new segment start/end points with respect to the intersection point found
            switch (segClass) {
            case LOCUS:
                //Determine which intersection point is on locus
            	if ((intPoints[0].latitude != 0) && (intPoints[0].longitude != 0)) {
            		if (ptIsOnLocus(*segmentLocus, intPoints[0],
            			&dummyPtIsOnLocusPt, &err, tol, eps))
            		{
            			tempIntPoint = intPoints[0];
            		}
            	} else if ((intPoints[1].latitude != 0) && (intPoints[1].longitude != 0)) {
                	if (ptIsOnLocus(*segmentLocus, intPoints[1],
                        &dummyPtIsOnLocusPt, &err, tol, eps))
                	{
                		tempIntPoint = intPoints[1];
					}
            	} else
                {
                    // CANNOT CONTINUE //
                    err |= UNEXPECTED_ERR;
                    return err;
                }
                if (startPointIsInside && endPointIsInside)
                {
                    //entire segment is inside circle
                    newStartPoint = segmentLocus->locusStart;
                    newEndPoint = segmentLocus->locusEnd;

                }
                else if (startPointIsInside)
                {
                    //check if intersection point and start point are the same
                    if (ptsAreSame(tempIntPoint,
                            segmentLocus->locusStart, tol))
                    {
                        err |= createPt(&newStartPoint,
                                segmentLocus->locusStart.latitude,
                                segmentLocus->locusStart.longitude);
                        segClass = LLPOINT;
                    }
                    else
                    {
                        newStartPoint = segmentLocus->locusStart;
                        newEndPoint = tempIntPoint;
                    }

                }
                else if (endPointIsInside)
                {
                    //check if intersection point and end point are the same
                    if (ptsAreSame(tempIntPoint,
                            segmentLocus->locusEnd, tol))
                    {
                        err |= createPt(&newStartPoint,
                                segmentLocus->locusEnd.latitude,
                                segmentLocus->locusEnd.longitude);
                        segClass = LLPOINT;
                    }
                    else
                    {
                        newStartPoint = tempIntPoint;
                        newEndPoint = segmentLocus->locusEnd;
                    }
                }
                else
                {
                    //start/end points are outside of circle and an intersection exists
                    //Therefore must be tangent
                    tangent = 1;
                    break;
                }
                break;

            case GEODESIC:

                //Determine which intersection point is on geodesic
                if (ptIsOnGeo(segmentGeo->startPoint,
                        segmentGeo->endPoint, intPoints[0],
                        segmentGeo->lineType, &err, tol, eps))
                {
                    tempIntPoint = intPoints[0];
                }
                else if (ptIsOnGeo(segmentGeo->startPoint,
                        segmentGeo->endPoint, intPoints[1],
                        segmentGeo->lineType, &err, tol, eps))
                {
                    tempIntPoint = intPoints[1];
                }
                else
                {
                    // CANNOT CONTINUE //
                    err |= UNEXPECTED_ERR;
                    return err;
                }

                if (startPointIsInside && endPointIsInside)
                {
                    //entire segment is inside circle
                    newStartPoint = segmentGeo->startPoint;
                    newEndPoint = segmentGeo->endPoint;
                }
                else if (startPointIsInside)
                {
                    //check if intersection point and start point are the same
                    if (ptsAreSame(tempIntPoint,
                            segmentGeo->startPoint, tol))
                    {
                        err |= createPt(&newStartPoint,
                                segmentGeo->startPoint.latitude,
                                segmentGeo->startPoint.longitude);
                        segClass = LLPOINT;
                    }
                    else
                    {
                        newStartPoint = segmentGeo->startPoint;
                        newEndPoint = tempIntPoint;
                    }
                }
                else if (endPointIsInside)
                {
                    //check if intersection point and end point are the same
                    if (ptsAreSame(tempIntPoint, segmentGeo->endPoint,
                            tol))
                    {
                        err |= createPt(&newStartPoint,
                                segmentGeo->endPoint.latitude,
                                segmentGeo->endPoint.longitude);
                        segClass = LLPOINT;
                    }
                    else
                    {
                        newStartPoint = tempIntPoint;
                        newEndPoint = segmentGeo->endPoint;
                    }
                }
                else
                {
                    //start/end points are outside of circle and an intersection exists
                    //Therefore must be tangent
                    tangent = 1;
                }

                break;

            case ARC:
                //Determine which intersection point is on geodesic
                if (ptIsOnArc(segmentArc->centerPoint,
                        segmentArc->radius, segmentArc->startAz,
                        segmentArc->endAz, segmentArc->dir, intPoints[0], &err,
                        tol, eps))
                {
                    tempIntPoint = intPoints[0];
                }
                else if (ptIsOnArc(segmentArc->centerPoint,
                        segmentArc->radius, segmentArc->startAz,
                        segmentArc->endAz, segmentArc->dir, intPoints[1], &err,
                        tol, eps))
                {
                    tempIntPoint = intPoints[1];
                }
                else
                {
                    // CANNOT CONTINUE //
                    err |= UNEXPECTED_ERR;
                    return err;
                }

                if (startPointIsInside && endPointIsInside)
                {
                    //entire segment is inside circle
                    newStartPoint = segmentArc->startPoint;
                    newEndPoint = segmentArc->endPoint;
                }
                else if (startPointIsInside)
                {
                    //check if intersection point and start point are the same
                    if (ptsAreSame(tempIntPoint,
                            segmentArc->startPoint, tol))
                    {
                        err |= createPt(&newStartPoint,
                                segmentArc->startPoint.latitude,
                                segmentArc->startPoint.longitude);
                        segClass = LLPOINT;
                    }
                    else
                    {
                        newStartPoint = segmentArc->startPoint;
                        newEndPoint = tempIntPoint;
                    }
                }
                else if (endPointIsInside)
                {
                    //check if intersection point and end point are the same
                    if (ptsAreSame(tempIntPoint, segmentArc->endPoint,
                            tol))
                    {
                        err |= createPt(&newStartPoint,
                                segmentArc->endPoint.latitude,
                                segmentArc->endPoint.longitude);
                        segClass = LLPOINT;
                    }
                    else
                    {
                        newStartPoint = tempIntPoint;
                        newEndPoint = segmentArc->endPoint;
                    }
                }
                else
                {
                    //start/end points are outside of circle and an intersection exists
                    //Therefore must be tangent
                    tangent = 1;
                }

                break;

            case LLPOINT:
                break;

            default:
                // CANNOT CONTINUE //
                err |= SHAPE_NOT_DEFINED_ERR;
                return err;
            }

            //Save the intersection point on the circle
            arcPoints[nArcPoints] = tempIntPoint;

            //Save the start true course from the intersection point to the circle center point
            err |= invCrs(tempIntPoint, circle.centerPoint,
                    &arcPointsAngle[nArcPoints], &dummyCourse, eps);

            nArcPoints++;

            if (tangent)
                continue;

            break;

        case 2: //Two intersection points found

            newStartPoint = intPoints[0];
            newEndPoint = intPoints[1];

            //check if correct orientation is set correctly for arcs and loci
            if (segClass == ARC)
            {

            	//find the mid point between newStartPoint and newEndPoint
            	err |= invCrs(segmentArc->centerPoint, newStartPoint, &crsSegCtrToNewStart, &tempCrs, eps);
            	err |= invCrs(segmentArc->centerPoint, newEndPoint, &crsSegCtrToNewEnd, &tempCrs, eps);
            	subtendedAngle = computeSubtendedAngle(crsSegCtrToNewStart, crsSegCtrToNewEnd, segmentArc->dir);
            	crsSegCtrToMidPt = crsSegCtrToNewStart + 0.5*subtendedAngle;
            	err |= direct(segmentArc->centerPoint, crsSegCtrToMidPt, segmentArc->radius, &testArcMidPt, eps);

            	//check if midpoint is on segment arc, if not then reverse new start/end points
                if (!(ptIsOnArc(segmentArc->centerPoint, segmentArc->radius, segmentArc->startAz, segmentArc->endAz, segmentArc->dir,
                        testArcMidPt, &err, tol, eps)))
                {
                    newStartPoint = intPoints[1];
                    newEndPoint = intPoints[0];
                }

            }
            else if (segClass == LOCUS)
            {
                //Find the point on the geodesic that corresponds to intersection point intPoints[0]
                ptIsOnLocus(*segmentLocus, newStartPoint, &tempPt, &err,
                        tol, eps);

                //Construct test locus
                err |= createGeo(&testLine, (*segmentLocus).geoStart,
                        tempPt, (*segmentLocus).lineType, eps);

                //Find the point on the geodesic that corresponds to intersection point intPoints[1]
                ptIsOnLocus(*segmentLocus, newEndPoint, &tempPt, &err,
                        tol, eps);

                //check if new EndPoint is between the locus start point and intPoints[0].  If it is then reverse new start/end points.
                if (ptIsOnGeo(testLine.startPoint, testLine.endPoint,
                        tempPt, testLine.lineType, &err, tol, eps))
                {
                    newStartPoint = intPoints[1];
                    newEndPoint = intPoints[0];
                }
            }

            //If segment lies ON circle, then save the common arc as
            //part of the boundary segment.  It would be redundant to save
            //the circle portion
            if (!coincide)
            {
                for (j = 0; j < n; j++)
                {
                    //Save the intersection point on the arc
                    arcPoints[nArcPoints] = intPoints[j];

                    //Save the start true course from the intersection point to the circle center point
                    err |= invCrs(intPoints[j], circle.centerPoint,
                            &arcPointsAngle[nArcPoints], &dummyCourse, eps);

                    nArcPoints++;
                }
            }

            break;

        default:
            // CANNOT CONTINUE //
            err |= UNEXPECTED_ERR;
            return err;
        }
        //Create the new segment
        switch (segClass) {
        case LOCUS:

            //Set the temporary start point equal to point on geodesic from locus point
            ptIsOnLocus(*segmentLocus, newStartPoint, &tempStartPoint,
                    &err, tol, eps);

            //Set the temporary end point equal to point on geodesic from locus point
            ptIsOnLocus(*segmentLocus, newEndPoint, &tempEndPoint, &err,
                    tol, eps);

            //Set the distance
            tempStartDistance = distToLocusFromGeoPt(*segmentLocus,
                    tempStartPoint, NULL, &err, tol, eps);

            tempEndDistance = distToLocusFromGeoPt(*segmentLocus, tempEndPoint,
                    NULL, &err, tol, eps);

            //Allocate new memory for new segment
            newSegment = malloc(sizeof(Locus));

            //Create the new segment locus
            err |= createLocus((Locus*) newSegment, tempStartPoint,
                    tempEndPoint, tempStartDistance, tempEndDistance,
                    segmentLocus->lineType, tol, eps);

            err |= addLocusToBndry(newBoundary, (Locus*) newSegment);

            break;

        case GEODESIC:

            //Allocate new memory for new segment
            newSegment = malloc(sizeof(Geodesic));

            //Create the new segment geodesic
            err |= createGeo((Geodesic*) newSegment, newStartPoint,
                    newEndPoint, segmentGeo->lineType, eps);

            err |= addGeoToBndry(newBoundary, (Geodesic*) newSegment);

            break;

        case ARC:

            //Allocate new memory for new segment
            newSegment = malloc(sizeof(Arc));

            err |= createArc((Arc*) newSegment, segmentArc->centerPoint,
                    newStartPoint, newEndPoint, segmentArc->dir, tol, eps);

            err |= addArcToBndry(newBoundary, (Arc*) newSegment);

            break;

        case LLPOINT:
            discretePoints[nDiscretePoints] = newStartPoint;
            nDiscretePoints++;

            break;

        default:
            // CANNOT CONTINUE //
            err |= SHAPE_NOT_DEFINED_ERR;
            return err;
        }

    } //End for i = 1:segmentCount
    // Identify which subarcs need to be returned
    if (nArcPoints > 1)
    {

        /*
         * Put the intersection points in order based on azimuth from center of
         * arc (see below for sort algorithm). This is okay because
         * circles always have positive orientation.
         */
        err |= sortPtsByAz(arcPoints, nArcPoints, arcPointsAngle, CLOCKWISE);

        /*
         * Now neighboring pairs of points define subarcs that may need to be
         * included.  Test each arc; if it is inside boundary then include it
         * in createBndry.
         */

        // copy first point to end of list for complete subarc construction
        arcPoints[nArcPoints] = arcPoints[0];

        // copy first point azimuth to end of list for complete subarc construction
        arcPointsAngle[nArcPoints] = arcPointsAngle[0];

        //Create test arcs from the list of ordered arc points.
        //Test each test arc to determine whether it is inside boundary.
        //If the test arc is inside the boundary then add the test arc to the new boundary.
        for (i = 0; i < nArcPoints; i++)
        {
            //Create the arc to be tested
            tempStartPoint = arcPoints[i];
            tempEndPoint = arcPoints[i + 1];

            //Check if points are identical
            if (ptsAreSame(tempStartPoint, tempEndPoint, tol))
            {
                //if so then move onto next pair of points
                continue;
            }

            err |= createArc(&testArc, circle.centerPoint, tempStartPoint,
                    tempEndPoint, CLOCKWISE, tol, eps);

            //Get the midpoint on the test arc

            //Determine the magnitude of the angle between the start and end true course
            err |= minSubtendedAngle(testArc.startAz, testArc.endAz,
                    &angleBtwn);

            //Determine the half angle
            halfAngleBtwn = angleBtwn / 2.0;

            //Determine the angle bisector course
            if (testArc.dir == CLOCKWISE ) //arc direction is clockwise
            {
                angleBisectorCourse = fmod(testArc.startAz + halfAngleBtwn,
                        M_2PI); //modulo 2Pi (i.e. 360 degrees)
            }
            else
            {
                angleBisectorCourse = fmod(testArc.startAz - halfAngleBtwn,
                        M_2PI); //modulo 2Pi (i.e. 360 degrees)
            }

            //Determine the point that is the arc radius distance along the angle bisector course from the arc center point
            err |= direct(testArc.centerPoint, angleBisectorCourse,
                    testArc.radius, &testPoint, eps);

            //check if test arc is inside boundary
            if (ptIsInsideBndry(boundary, testPoint, &err, tol, eps))
            {
                err |= addArcToBndry(newBoundary, &testArc);
            }
        }

    }
    /*
     * No intersections (other than discrete points) were found.  If center
     * of circle is inside boundary, then entire circle must be inside
     * boundary, so return all of circle as createBndry.
     */
    else if (newBoundary->length == 0)
    {
        if (ptIsInsideBndry(boundary, circle.centerPoint, &err, tol, eps))
        {
            err |= addArcToBndry(newBoundary, &circle);
        }
        else if (nArcPoints == 1)
        {
            //If an intersection exists and circle is outside of boundary then
            //circle must just graze the boundary.  Treat this as a discrete point.
            discretePoints[nDiscretePoints] = arcPoints[0];
            nDiscretePoints++;
        }
    }

    /*
     * Now check discrete points against Arcs and LineSegments in createBndry.
     * Keep only the discrete points that do not lie on an Arc or LineSegment of
     * createBndry.
     */
    if (newBoundary->length != 0)
    {
        for (i = 0; i < nDiscretePoints; i++)
        {
            testPoint = discretePoints[i];
            if (!ptIsInsideBndry(*newBoundary, testPoint, &err, tol, eps))
            {
                err |= addElementToBndry(newBoundary, &testPoint, LLPOINT);
            }
        }
    }
    else if (nDiscretePoints)
    {

        for (i = 0; i < nDiscretePoints; i++)
        {
            //Check if points are identical (Note first point is always included)
            if (i > 0 && ptsAreSame(discretePoints[i], discretePoints[i
                    - 1], tol))
            {
                //if so then move onto next point
                continue;
            }
            err |= addElementToBndry(newBoundary, &discretePoints[i],
                    LLPOINT);
        }
    }

    //check if circle is inside boundary -- fix for Bug 25861
    if (newBoundary->length == 0)
    {
    	if(ptIsInsideBndry(boundary, circle.centerPoint, &err, tol, eps))
          addElementToBndry(newBoundary, &circle, ARC);
    }
    else
    {
      pointsOnly = 1;
      for (i = 0; i < newBoundary->length; i++)
      {
        if (newBoundary->elements[i].type != LLPOINT)
        {
          pointsOnly = 0;
          break;
        } 
      }

      if (pointsOnly)
      {
    	if(ptIsInsideBndry(boundary, circle.centerPoint, &err, tol, eps))
        {
          newBoundary->length = 0;
          addElementToBndry(newBoundary, &circle, ARC);
        }
      }
    }

    return err;
}

/*
 * Determine whether a circle "intersects" a simple Boundary structure.
 *
 * Case 1: checkSurface is true (i.e. non-zero)
 * In this case, the methods treats the boundary as a surface
 * and if the circle lies on the boundary or intersects the boundary then true is returned.
 * Otherwise false is returned.
 *
 * Case 2: checkSurface is false. (i.e. zero)
 * In this case, the method only evaluates against the shapes that comprise the boundary
 * itself.  If any intersection occurs with any shape that comprises the boundary, then
 * true is returned.  Otherwise false is returned.
 *
 *
 * A simple boundary implies one that contains no holes in the boundary area.  This also
 * means that the boundary area is not partitioned (contains sub areas) in any way.
 *
 *                    not partitioned          partitioned
 *                       ---------              ---------
 *                       |       |              |   |   |
 *                       |       |              |   |   |
 *                       ---------              ---------
 *
 * A connected boundary implies that all the start/end points are "connected" to one
 * another.  Two shapes are connected if they share an identical point or if the
 * Vincenty distance between the two points is less than tol in nautical miles.
 *
 * The speed of this algorithm is directly proportional to the number of shapes that
 * comprise the boundary.  More boundary shapes implies more time to process.
 *
 * inputs:
 * 		boundary = Boundary structure of interest
 * 		circle = Circle to be evaluated against boundary
 * 		checkSurface = true: evaluate boundary as a surface, false: evaluate boundary shapes only
 * 		tol = distance criteria (in nmi)
 * 		eps = Vincenty forward/inverse algorithm convergence criteria
 */
ErrorSet bndryCircleIntxExists(Boundary boundary, Arc circle, int checkSurface,
                                    int* intersectionExists, double tol, double eps)
{
    ErrorSet err = 0; //The error code, 0 if successful
    int i = 0;
    double circleExtent; //The signed arc extent of Arc circle
    Shape* segment = NULL; //The shape currently being processed
    ShapeType segClass; //The shape type of the current shape being processed
    Arc* segmentArc = NULL; //Temp object to hold segment recast as an Arc
    Geodesic* segmentGeo = NULL; //Temp object to hold segment recast as an Geodesic
    Locus* segmentLocus = NULL; //Temp object to hold segment recast as an Locus
    LLPoint* segmentLLPoint = NULL; //Temp object to hold segment recast an an LLPoint
    LLPointPair intPoints; //Holds the intersection points for a particluar segment-cirlce intersection
    int nFnd = -1; //The number of intersection points found by an intersection method.  Not necessarily same as n.
    int coincide = 0; //1 = arc and circle coincide
    int nCoincide = 0; //The number of shapes that coincide between a circle and arc
    Shape coincideCommonShapes[2]; //Common shapes if arc and circle coincide
    LineType tempLineType;
    LLPoint tempPt; //Generic temporary point
    double arcExtent; //The signed arc extent of a shape of type Arc
    int intersectFound = 0;//0 = false, 1 = true

    // Check for nulls
    if ((NULL == &boundary) || (NULL == &circle) || (NULL == createBndry))
    {
        /* CANNOT CONTINUE */
        err |= NO_MEMORY_ALLOCATED_ERR;
        return err;
    }

    // Get the angle extent of the circle
    err |= getArcExtent(circle.centerPoint, circle.radius,
            circle.startAz, circle.endAz, circle.dir, &circleExtent, tol, eps);
    if (err)
        return err;

    // Check if full circle
    if (fabs(circleExtent) < M_2PI)
    {
        /* CANNOT CONTINUE */
        err |= SUBTENDED_ANGLE_OUT_OF_RANGE_ERR;
        return err;
    }

    //check if circle center point lies on the surface of boundary
    if (checkSurface)
    {
    	//use ptIsInsideBndry method to complete this check for the moment, may need to change depending on performance
    	intersectFound = ptIsInsideBndry(boundary, circle.centerPoint, &err, tol, eps);
    }


    //Check if circle intersects any element of the boundary
    while (!intersectFound && i < boundary.length)
    {
        //Get the shape from the boundary
        segment = boundary.elements[i].this;

        //Get the shape type
        segClass = boundary.elements[i].type;

        coincide = 0;

        switch (segClass) {
        case LOCUS:

            //Cast shape as a Locus
            segmentLocus = (Locus*) segment;

            //Find intersections of segment with circle

            //First store the original lineType
            tempLineType = segmentLocus->lineType;

            //segmentLocus must be of lineType INFINITE to implement LocusArcIntersect properly
            segmentLocus->lineType = INFINITE;

            //Get the intersection points of the circle and locus
            err |= locusArcIntx(*segmentLocus, circle.centerPoint,
                    circle.radius, intPoints, &nFnd, tol, eps);

            //Set lineType back to original
            segmentLocus->lineType = tempLineType;

            //Check if found intersection points are on the shape
            intersectFound = nFnd;
            if ((intersectFound > 1) && !ptIsOnLocus(*segmentLocus, intPoints[1],
                    &tempPt, &err, tol, eps))
            {
                //If point is not on shape then do not count
                intersectFound--;
            }

            if ((intersectFound > 0) && !ptIsOnLocus(*segmentLocus, intPoints[0],
                    &tempPt, &err, tol, eps))
            {
                //If point is not on shape then do not count
                intersectFound--;
            }

            break;

        case GEODESIC:

            //Cast shape as a Geodesic
            segmentGeo = (Geodesic*) segment;

            //Find intersections of segment with circle

            // Mask out error codes CONCENTRIC_CIRCLE_ERR.
            // In this case, CONCENTRIC_CIRCLE_ERR is an exit indicator only and not an actual error.
            // See bug 16362
            err |= getMaskedError(geoArcIntx(segmentGeo->startPoint,
                    segmentGeo->startAz, circle.centerPoint, circle.radius,
                    intPoints, &nFnd, tol, eps), getMask(0, 1, 0, 0, 0, 0, 0));

            //Check if found intersection points are on the shape
            intersectFound = nFnd;
            if (intersectFound > 1)
            {
                if (!ptIsOnGeo(segmentGeo->startPoint,
                        segmentGeo->endPoint, intPoints[1],
                        segmentGeo->lineType, &err, tol, eps))
                {
                    //If point is not on shape then do not count
                    intersectFound--;
                }
            }

            if (intersectFound > 0)
            {
                if (!ptIsOnGeo(segmentGeo->startPoint,
                        segmentGeo->endPoint, intPoints[0],
                        segmentGeo->lineType, &err, tol, eps))
                {
                    //If point is not on shape then do not count
                    intersectFound--;
                }
            }

            break;

        case ARC:

            //Cast shape as an Arc
            segmentArc = (Arc*) segment;

            //Check if arcs coincide
            coincide = arcsCoincide(circle, *segmentArc, coincideCommonShapes,
                    &nCoincide, &err, tol, eps);

            if (coincide)
            {
                /*
                 * If the arc is a "circle" and it coincides with the boundary then it cannot be connected to
                 * any other shapes and still comprise a simple boundary.  Thus the boundary and circle must
                 * be the same.  Continue processing if this is the case.
                 *
                 * If the arc is not a "circle" and it coincides with the boundary then any shape connected
                 * to the arc's start/end points will generate an intersection point to be processed later on,
                 * assuming the boundary is a simple and connected boundary.
                 * In this case, ignore the arc and continue processing the next shape in the boundary.
                 *
                 * Note - arcIntx should return zero intersections for concentric arcs which
                 *        implies zero intersections will be found for arcs that coincide.
                 */

                //Get the signed arc extent with respect to tol
                err |= getArcExtent(segmentArc->centerPoint,
                        segmentArc->radius, segmentArc->startAz,
                        segmentArc->endAz, segmentArc->dir, &arcExtent, tol,
                        eps);

                //if arc is not a "circle" then move onto next shape
                if (fabs(arcExtent) < M_2PI)
                {
                    continue;
                }

            }

            //Get the intersection points of the circle and arc

            // Mask out error codes CIRCLE_INSIDE_CIRCLE_ERR.
            // In this case, CIRCLE_INSIDE_CIRCLE_ERR is an exit indicator only and not an actual error.
            // See bug 16362

            err |= getMaskedError(arcIntx(segmentArc->centerPoint,
                    segmentArc->radius, circle.centerPoint, circle.radius,
                    intPoints, &nFnd, tol, eps), getMask(0, 0, 0, 0, 0, 0, 1));

            //Check if found intersection points are on the shape
            intersectFound = nFnd;

            if (intersectFound > 1)
            {
                if (!ptIsOnArc(segmentArc->centerPoint,
                        segmentArc->radius, segmentArc->startAz,
                        segmentArc->endAz, segmentArc->dir, intPoints[1], &err,
                        tol, eps))
                {
                    //If point is not on shape then do not count
                    intersectFound--;
                }
            }

            if (intersectFound > 0)
            {
                if (!ptIsOnArc(segmentArc->centerPoint,
                        segmentArc->radius, segmentArc->startAz,
                        segmentArc->endAz, segmentArc->dir, intPoints[0], &err,
                        tol, eps))
                {
                    //If point is not on shape then do not count
                    intersectFound--;
                }
            }

            break;

        case LLPOINT:

            //Cast shape as an LLPoint
            segmentLLPoint = (LLPoint*) segment;

            //Determine location of segment point with respect to circle
            intersectFound = ptIsInsideArc(circle.centerPoint,
                    circle.radius, circle.startAz, circle.endAz, circle.dir,
                    *segmentLLPoint, &err, tol, eps);

            break;

        default:

            // CANNOT CONTINUE //
            err |= SHAPE_NOT_DEFINED_ERR;
            return err;
        }

        i++;
    }


    *intersectionExists = intersectFound;
    return err;

}

/*
 * Order a simple connected boundary by the sequence in which the shapes connect
 * at their respective start/end points, starting with the shape at index zero
 * of the original boundary's shape list.
 *
 * A simple boundary implies one that contains no holes in the boundary area.  This also
 * means that the boundary area is not partitioned (contains sub areas) in any way.
 *
 *                    not partitioned          partitioned
 *                       ---------              ---------
 *                       |       |              |   |   |
 *                       |       |              |   |   |
 *                       |       |              |   |   |
 *                       ---------              ---------
 *
 * Let unordered and ordered be Boundary structs.
 * Let the numbers suffix indicate the array index at which the given boundary shape is located.
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
 * in the shapes array of the unordered boundary, i.e. S_0 = unordered.elements[0], S_1 = unordered.elements[1], ...
 * In this sequence, there exists pairs of consecutive shapes in the array such that these consecutive shapes
 * are not always connected.  For example the shapes located at indexes 0 and 1 (S_0 and S_1) do not share a
 * common start/end point.
 *
 *
 * inputs:
 * 		boundary = Boundary to be ordered by connectivity
 * 		orderedBoundary = Boundary ordered by the connecting shapes.  Each shape in the boundary has a start/end point connection.
 * 		tol = distance criteria (in nmi)
 * 		eps = Vincenty forward/inverse algorithm convergence criteria
 */
ErrorSet orderBndry(Boundary boundary, Boundary* orderedBoundary, double tol, double eps)
{

	ErrorSet err = 0;
    Shape* segment = NULL; //The shape currently being processed
    Arc* segmentArc = NULL; //Temp object to hold segment recast as an Arc
    Geodesic* segmentGeo = NULL; //Temp object to hold segment recast as an Geodesic
    Locus* segmentLocus = NULL; //Temp object to hold segment recast as an Locus
    Spiral* segmentSpiral = NULL; //Temp object to hold segment recast as an Spiral
    Arc tempArc;
    Geodesic tempGeo;
    Locus tempLocus;
    Spiral tempSpiral;
    int i = 0;
    int j = 0;
    Shape* shapes = NULL;
    LLPoint firstPoint;
    LLPoint terminationPoint;
    LLPoint connectionPoint;
    LLPoint connectedPoint;
    LLPoint startPt;
    LLPoint endPt;
    int connection;//0 = no connection found, non-zero = connection found

    if ((NULL == orderedBoundary) )
    {
        err |= NO_MEMORY_ALLOCATED_ERR;
        return err;
    }

    //get the list of shapes
    shapes = boundary.elements;

    //initialize the termination point to be the start point of the first shape in the shape list
    //initialize the connection point to be the end point of the first shape in the shape list
    segment = shapes[0].this;
    switch (shapes[0].type) {
		case LOCUS:
			segmentLocus = (Locus*) segment;
			firstPoint = segmentLocus->locusEnd;
			terminationPoint = segmentLocus->locusStart;
			connectionPoint = segmentLocus->locusEnd;
			connectedPoint = terminationPoint;
			break;

		case GEODESIC:
			segmentGeo = (Geodesic*) segment;
			firstPoint = segmentGeo->endPoint;
			terminationPoint = segmentGeo->startPoint;
			connectionPoint = segmentGeo->endPoint;
			connectedPoint = terminationPoint;
			break;

		case ARC:
			segmentArc = (Arc*) segment;
			firstPoint = segmentArc->endPoint;
			terminationPoint = segmentArc->startPoint;
			connectionPoint = segmentArc->endPoint;
			connectedPoint = terminationPoint;
			break;

		case LLPOINT:
			//TODO do something
			break;

		case SPIRAL:
			segmentSpiral = (Spiral*) segment;
			firstPoint = segmentSpiral->endPoint;
			terminationPoint = segmentSpiral->startPoint;
			connectionPoint = segmentSpiral->endPoint;
			connectedPoint = terminationPoint;
			break;
        default:
        	break;
    }

    //add the first shape in the shape list to the ordered boundary
    err |= addElementToBndry(orderedBoundary, segment, shapes[0].type);

	//start looking for a connecting shape by searching through the entire shape list
    for (i = 1; i < (boundary.length); i++)
    {

    	connection = 0;
    	j = 0;
        while ( (!connection) && (j < boundary.length) )
        {
        	//get the shape evaluation points
            segment = shapes[j].this;
			switch (shapes[j].type)
			{
				case LOCUS:
					segmentLocus = (Locus*) segment;
					startPt = segmentLocus->locusStart;
					endPt = segmentLocus->locusEnd;
					break;

				case GEODESIC:
					segmentGeo = (Geodesic*) segment;
					startPt = segmentGeo->startPoint;
					endPt = segmentGeo->endPoint;
					break;

				case ARC:
					segmentArc = (Arc*) segment;
					startPt = segmentArc->startPoint;
					endPt = segmentArc->endPoint;
					break;

				case LLPOINT:
					//TODO do something
					break;

				case SPIRAL:
					segmentSpiral = (Spiral*) segment;
					startPt = segmentSpiral->startPoint;
					endPt = segmentSpiral->endPoint;
					break;
		        default:
		        	break;
			}

			//check the evaluation points against the current connection point
			if (ptsAreSame(connectionPoint, startPt, tol))
			{
				//possibly connected at evalPoint1 to the current shape

				//check the other evaluation point against most recent connection point
				if (!ptsAreSame(connectedPoint, endPt, tol) || (boundary.length == 2 && j == 1))
				{
					//valid connection, store this shape as the next shape in ordered boundary
					err |= addElementToBndry(orderedBoundary, segment, shapes[j].type);

					//update the connection points
					connectedPoint = startPt;
					connectionPoint = endPt;

					//set the connection flag
					connection = 1;
				}

			}
			else if (ptsAreSame(connectionPoint, endPt, tol))
			{
				//possibly connected at evalPoint2 to the current shape

				//check the other evaluation point against most recent connection point
				if (!ptsAreSame(connectedPoint, startPt, tol) || (boundary.length == 2 && j == 1))
				{
					//valid connection, store this shape as the next shape in ordered boundary

					//reorder start/end points
					switch(shapes[j].type){
						case GEODESIC:
							err |= createGeo(&tempGeo, segmentGeo->endPoint, segmentGeo->startPoint, segmentGeo->lineType, eps);
							err |= addGeoToBndry(orderedBoundary, &tempGeo);
							break;
						case LOCUS:
							err |= createLocus(&tempLocus, segmentLocus->geoEnd, segmentLocus->geoStart, -1 * segmentLocus->endDist, -1 * segmentLocus->startDist, segmentLocus->lineType, tol, eps);
							err |= addLocusToBndry(orderedBoundary, &tempLocus);
							break;
						case ARC:
							err |= createArc(&tempArc, segmentArc->centerPoint, segmentArc->endPoint, segmentArc->startPoint, -1 * (segmentArc->dir), tol, eps);
							err |= addArcToBndry(orderedBoundary, &tempArc);
							break;
						case LLPOINT:
							//TODO do something
							break;
						case SPIRAL:
							err |= createSpiral(&tempSpiral, segmentSpiral->centerPoint, segmentSpiral->endRadius, segmentSpiral->startRadius, segmentSpiral->endAz, segmentSpiral->startAz, -1 * (segmentSpiral->dir), eps);
							err |= addSpiralToBndry(orderedBoundary, &tempSpiral);
							break;
				        default:
				        	break;
					}

					//update the connection points
					connectedPoint = endPt;
					connectionPoint = startPt;

					//set the connection flag
					connection = 1;
				}
			}

			j++;
        }

    }

    //verify that the ordered boundary is of the same size as the original
    if (orderedBoundary->length != boundary.length)
    {
    	err |= UNEXPECTED_ERR;
    }

    return err;
}

/*
 * Separate a boundary into simple ordered boundaries
 *
 *
 * inputs:
 * 		boundary = Boundary to be ordered by connectivity
 * 		orderedBoundary[] = Array of separated boundaries
 * 		tol = distance criteria (in nmi)
 * 		eps = Vincenty forward/inverse algorithm convergence criteria
 */
ErrorSet separateBndry(Boundary boundary, Boundary separatedBoundaries[], int* numberOfBoundaries, double tol, double eps)
{

	ErrorSet err = 0;
    Arc* segmentArc = NULL; //Temp object to hold segment recast as an Arc
    Geodesic* segmentGeo = NULL; //Temp object to hold segment recast as an Geodesic
    Locus* segmentLocus = NULL; //Temp object to hold segment recast as an Locus
    Spiral* segmentSpiral = NULL; //Temp object to hold segment recast as an Spiral
    Arc tempArc;
    Geodesic tempGeo;
    Locus tempLocus;
    Spiral tempSpiral;

    int i = 0;
    Shape* shapes = NULL;
    LLPoint firstPoint;
    LLPoint connectionPoint;

    LLPoint startPoints[boundary.length];
    LLPoint endPoints[boundary.length];
    int used[boundary.length];
    int unused = boundary.length;
    int connectingSegmentFound = 0;

    if ((NULL == separatedBoundaries) )
    {
        err |= NO_MEMORY_ALLOCATED_ERR;
        return err;
    }

    //get the list of shapes
    shapes = boundary.elements;

    for(i = 0; i < boundary.length; i++){
		switch (shapes[i].type) {
			case LOCUS:
				segmentLocus = (Locus*) shapes[i].this;
				startPoints[i] = segmentLocus->locusStart;
				endPoints[i] = segmentLocus->locusEnd;
				break;

			case GEODESIC:
				segmentGeo = (Geodesic*) shapes[i].this;
				startPoints[i] = segmentGeo->startPoint;
				endPoints[i] = segmentGeo->endPoint;
				break;

			case ARC:
				segmentArc = (Arc*) shapes[i].this;
				startPoints[i] = segmentArc->startPoint;
				endPoints[i] = segmentArc->endPoint;
				break;

			case LLPOINT:
				//TODO do something
				break;

			case SPIRAL:
				segmentSpiral = (Spiral*) shapes[i].this;
				startPoints[i] = segmentSpiral->startPoint;
				endPoints[i] = segmentSpiral->endPoint;
				break;
	        default:
	        	break;
		}
		used[i] = 0;
    }

    int k = 0;
    while(unused > 0){
    	separatedBoundaries[k] = createBndry();
    	i = 0;
    	// Find the first unused element and start there
    	for(i = 0; i < boundary.length; i++){
    		if(used[i] == 0){
    			// Add first unused element to separatedBoundaries
    			err |= addElementToBndry(&separatedBoundaries[k], shapes[i].this, shapes[i].type);

    			// Assign first point and connection Point
				firstPoint = startPoints[i];
				connectionPoint = endPoints[i];

				used[i] = 1;
    			unused--;
    			break;
    		}
    	}

    	// Search for all the connected segments
    	while(!ptsAreSame(firstPoint,connectionPoint,tol)){

    		connectingSegmentFound = 0;
    		for(i = 0; i < boundary.length; i++){
    			if(used[i] == 0){
    				if(ptsAreSame(connectionPoint,startPoints[i],tol)){

    	    			err |= addElementToBndry(&separatedBoundaries[k], shapes[i].this, shapes[i].type);

    					connectionPoint = endPoints[i];

    					used[i] = 1;
    	    			unused--;
    	    			connectingSegmentFound = 1;
    	    			break;

    				} else if(ptsAreSame(connectionPoint,endPoints[i],tol)){
    					//reorder start/end points
    					switch(shapes[i].type){
    						case GEODESIC:
    							segmentGeo = (Geodesic*) shapes[i].this;
    							err |= createGeo(&tempGeo, segmentGeo->endPoint, segmentGeo->startPoint, segmentGeo->lineType, eps);
    							err |= addGeoToBndry(&separatedBoundaries[k], &tempGeo);
    							break;
    						case LOCUS:
    							segmentLocus = (Locus*) shapes[i].this;
    							err |= createLocus(&tempLocus, segmentLocus->geoEnd, segmentLocus->geoStart, -1 * segmentLocus->endDist, -1 * segmentLocus->startDist, segmentLocus->lineType, tol, eps);
    							err |= addLocusToBndry(&separatedBoundaries[k], &tempLocus);
    							break;
    						case ARC:
    							segmentArc = (Arc*) shapes[i].this;
    							err |= createArc(&tempArc, segmentArc->centerPoint, segmentArc->endPoint, segmentArc->startPoint, -1 * (segmentArc->dir), tol, eps);
    							err |= addArcToBndry(&separatedBoundaries[k], &tempArc);
    							break;
    						case LLPOINT:
    							//TODO do something
    							break;
    						case SPIRAL:
    							segmentSpiral = (Spiral*) shapes[i].this;
    							err |= createSpiral(&tempSpiral, segmentSpiral->centerPoint, segmentSpiral->endRadius, segmentSpiral->startRadius, segmentSpiral->endAz, segmentSpiral->startAz, -1 * (segmentSpiral->dir), eps);
    							err |= addSpiralToBndry(&separatedBoundaries[k], &tempSpiral);
    							break;
    				        default:
    				        	break;
    					}

    					connectionPoint = startPoints[i];

    					used[i] = 1;
    	    			unused--;
    	    			connectingSegmentFound = 1;
    	    			break;

    				}

    			}

    		}

    		if(connectingSegmentFound == 0){
    			err |= SHAPE_NOT_DEFINED_ERR;
    			return err;
    		}

    	}
    	k++;
    }

    *numberOfBoundaries = k;

    return err;
}


int ptIsInsideBndry(Boundary b, LLPoint p, ErrorSet* err, double tol, double eps)
{

    double az = 0.0, dist, fcrs, temp;
    ErrorSet newErr = 0;

    LLPoint* intxList = NULL;
    Arc* thisArc = NULL;
    Geodesic* thisGeo = NULL;
    Locus* thisLocus = NULL;
    Spiral* thisSpiral = NULL;
    Shape* thisShape = NULL;
    ShapeType thisType;
    LLPoint tempLLPoint;

    int intxCount = 0;
    int intxCode = 0;

    int done = 0;
    int isInside = 0;
    int i, azChanged = 0;
    
    for (i = 0; i < b.length; i++)
    {
      thisShape = b.elements[i].this;
      if (thisShape == NULL)
      {
        return SHAPE_NOT_DEFINED_ERR;
      }
      thisType = b.elements[i].type;
      switch (thisType) {
      case ARC:
        thisArc = (Arc*) thisShape;
        if (ptIsOnArc(thisArc->centerPoint, thisArc->radius, thisArc->startAz, thisArc->endAz, thisArc->dir, p, &newErr, tol, eps))
        {
          isInside = 1; 
          done = 1;
          break;
        }
        newErr |= invDist(p, thisArc->startPoint, &dist, eps);
        if (dist > tol && dist < 5.0 * tol)
        {
          newErr |= invCrs(p, thisArc->startPoint, &fcrs, &temp, eps);
          az = modcrs(fcrs + M_PI);
          break;
        }
        newErr |= invDist(p, thisArc->endPoint, &dist, eps);
        if (dist > tol && dist < 5.0 * tol)
        {
          newErr |= invCrs(p, thisArc->endPoint, &fcrs, &temp, eps);
          az = modcrs(fcrs + M_PI);
          break;
        }
        break;
      case GEODESIC:
        thisGeo = (Geodesic*) thisShape;    
        if (ptIsOnGeo(thisGeo->startPoint, thisGeo->endPoint, p, thisGeo->lineType, &newErr, tol, eps))
        {
          isInside = 1; 
          done = 1;
          break;
        }
        if (!azChanged)
        {
        newErr |= invDist(p, thisGeo->startPoint, &dist, eps);
        if (dist > tol && dist < 5.0 * tol)
        {
          newErr |= invCrs(p, thisGeo->startPoint, &fcrs, &temp, eps);
          az = modcrs(fcrs + M_PI);
          azChanged = 1;
          break;
        }
        }
        if (!azChanged)
        { 
        newErr |= invDist(p, thisGeo->endPoint, &dist, eps);
        if (dist > tol && dist < 5.0 * tol)
        {
          newErr |= invCrs(p, thisGeo->endPoint, &fcrs, &temp, eps);
          az = modcrs(fcrs + M_PI);
          azChanged = 1;
          break;
        }
        }
        if (!azChanged)
        {
          newErr |= projectToGeo(thisGeo->startPoint, thisGeo->startAz, p, &tempLLPoint, &temp, &dist, tol, eps);
          if (dist > tol && dist < 5.0 * tol)
          {
            az = modcrs(temp + M_PI);
            azChanged = 1;
          }
        }
        break;
      case LOCUS:
        thisLocus = (Locus*) thisShape;
        if (ptIsOnLocus(*thisLocus, p, &tempLLPoint, &newErr, tol, eps))
        {
          isInside = 1; 
          done = 1;
          break;
        }
        newErr |= invDist(p, thisLocus->locusStart, &dist, eps);
        if (dist > tol && dist < 5.0 * tol)
        {
          newErr |= invCrs(p, thisLocus->locusStart, &fcrs, &temp, eps);
          az = modcrs(fcrs + M_PI);
          break;
        }
        newErr |= invDist(p, thisLocus->locusEnd, &dist, eps);
        if (dist > tol && dist < 5.0 * tol)
        {
          newErr |= invCrs(p, thisLocus->locusEnd, &fcrs, &temp, eps);
          az = modcrs(fcrs + M_PI);
          break;
        }
        break;
      case SPIRAL:
    	thisSpiral = (Spiral*) thisShape;
		if (ptIsOnSpiral(*thisSpiral, p, tol, eps))
		  {
			isInside = 1;
			done = 1;
			break;
		  }
		  newErr |= invDist(p, thisSpiral->startPoint, &dist, eps);
		  if (dist > tol && dist < 5.0 * tol)
		  {
			newErr |= invCrs(p, thisSpiral->startPoint, &fcrs, &temp, eps);
			az = modcrs(fcrs + M_PI);
			break;
		  }
		  newErr |= invDist(p, thisSpiral->endPoint, &dist, eps);
		  if (dist > tol && dist < 5.0 * tol)
		  {
			newErr |= invCrs(p, thisSpiral->endPoint, &fcrs, &temp, eps);
			az = modcrs(fcrs + M_PI);
			break;
		  }
        break;
      case LLPOINT:
        break;
      default:
      	break;
      } //switch
        
    } //for i
    newErr = 0;  //We don't want any errors in the prescreening phase to prevent further processing.

    while (!done)
    {

        newErr |= bndryGeoIntx(b, p, az, &intxList, &intxCount,
                &intxCode, tol, eps);

        if (newErr)
        {
            *err |= newErr;
            return 0;
        }

        switch (intxCode) {
        case POINT_IS_ON_BOUNDARY:
            isInside = 1;
            done = 1;
            break;
        case NEW_LINE_NEEDED:
            /* Rotate test line by 1 degree */
            //az = az + 1.0 * M_PI / 180.0;
            az = az + 37.0 * M_PI / 180.0;
            done = 0;
            break;
        default:
            /* point is inside if odd number of intersections */
            isInside = (intxCount & 1);
            done = 1;
        }

    }

    return isInside;

}

/* Return intersection points only if they lie on both arcs */
static ErrorSet arcShapeIntx(Arc a, Shape s, LLPointPair intx, int* n,
                                     int* exitCode, double tol, double eps)
{

    Geodesic* thisGeo = NULL;
    Arc* thisArc = NULL;
    Locus* thisLocus = NULL;

    LLPointPair tmpIntx;

    int m = 0;
    int i = 0;
    ErrorSet err = 0;
    ErrorSet newErr = 0;

    if ((NULL == intx) || (NULL == n))
    {
        err |= NO_MEMORY_ALLOCATED_ERR;
        return err;
    }

    switch (s.type) {
    case ARC:
        thisArc = (Arc*) s.this;
        err |= arcIntx(a.centerPoint, a.radius, thisArc->centerPoint,
                thisArc->radius, tmpIntx, &m, tol, eps);
        break;
    case GEODESIC:
        thisGeo = (Geodesic*) s.this;
        err |= geoArcIntx(thisGeo->startPoint, thisGeo->startAz,
                a.centerPoint, a.radius, tmpIntx, &m, tol, eps);
        break;
    case LOCUS:
        thisLocus = (Locus*) s.this;
        err |= locusArcIntx(*thisLocus, a.centerPoint, a.radius,
                tmpIntx, &m, tol, eps);
        break;
    case SPIRAL:
        break;
    case LLPOINT:
        //TODO do something
        break;
    }

    /* Keep only the points that are on the given arc */
    *n = 0;
    if (!err)
    {
        for (i = 0; i < m; i++)
        {
            newErr = 0;
            if (ptIsOnArc(a.centerPoint, a.radius, a.startAz, a.endAz,
                    a.dir, tmpIntx[i], &newErr, tol, eps))
            {
                /* point is on arc */
                if (ptIsOnShape(s, tmpIntx[i], &newErr, tol, eps))
                {
                    /* point is also on boundary shape */
                    if (!newErr)
                    {
                        intx[*n] = tmpIntx[i];
                        *n = *n + 1;
                    }
                }
            }
            err |= newErr;
        }
    }

    return err;

}

static int ptIsOnShape(Shape s, LLPoint p, ErrorSet* err, double tol,
                          double eps)
{

    Geodesic* thisGeo = NULL;
    Arc* thisArc = NULL;
    Locus* thisLocus = NULL;
    //	Spiral* thisSpiral = NULL;

    LLPoint geoPt;

    int result = 0;
    ErrorSet newErr = 0;

    if (NULL == s.this)
    {
        *err |= NO_MEMORY_ALLOCATED_ERR;
        return 0;
    }

    switch (s.type) {
    case GEODESIC:
        thisGeo = (Geodesic*) s.this;
        result = ptIsOnGeo(thisGeo->startPoint, thisGeo->endPoint, p,
                thisGeo->lineType, &newErr, tol, eps);
        break;
    case ARC:
        thisArc = (Arc*) s.this;
        result = ptIsOnArc(thisArc->centerPoint, thisArc->radius,
                thisArc->startAz, thisArc->endAz, thisArc->dir, p, &newErr,
                tol, eps);
        break;
    case LOCUS:
        thisLocus = (Locus*) s.this;
        result = ptIsOnLocus(*thisLocus, p, &geoPt, &newErr, tol, eps);
        break;
    case SPIRAL:
        break;
    case LLPOINT:
        //TODO do something
        break;
    }

    if (newErr)
    {
        *err |= newErr;
        return 0;
    }

    return result;

}


/*******************************************************************************
 * Test whether given arc lies on other given arc
 * Input: 	arc1 - The first arc
 * 			arc2 - The second arc
 * 			commonShapes - If arc1 and arc2 coincide then add the common arc or
 * 				common point to this list of pointers.
 * 			shapeCount - The number of common shapes (arcs and points) found
 * 			err - The pointer to the error variable
 * 			return 0 if arcs do not coincide or if error resulted.  return 1 otherwise.
 *
 */
int arcsCoincide(Arc arc1, Arc arc2, Shape commonShapes[], int* shapeCount,
                 ErrorSet* err, double tol, double eps)
{

    Arc newArc1 = arc1;
    int arc1IsCircle = 0;
    int arc2IsCircle = 0;
    LLPoint arcPoints[5];
    double arcPointAzs[5];
    LLPoint commonStart;
    LLPoint commonEnd;
    LLPoint commonMid;
    Arc testArc;
    Arc* arcPtr;
    int i = 0, j, skipPoint;
    int midPtOnArc1 = 0;
    int midPtOnArc2 = 0;
    ErrorSet newErr = 0;
    LLPoint* commonPoint = NULL;
    LLPoint testPoint;
    double angleBtwn;
    double halfAngleBtwn;
    double angleBisectorCourse;

    if (commonShapes == NULL || shapeCount == NULL || err == NULL)
    {
        if (err == NULL)
        {
            err = malloc(sizeof(ErrorSet));
        }
        *err |= NO_MEMORY_ALLOCATED_ERR;
        return 0;
    }

    //initialize shape count to zero
    *shapeCount = 0;

    //check if arcs are concentric.  Not possible to coincide if not concentric.
    if (!ptsAreSame(arc1.centerPoint, arc2.centerPoint, tol))
    {
        return 0;
    }

    //check if arc radii are same.  Not possible to coincide if not of same radius.
    if (fabs(arc1.radius - arc2.radius) > tol)
    {
        return 0;
    }

    //Make arcs same orientation if not already
    if (arc1.dir * arc2.dir < 0)
    {
        newErr |= createArc(&newArc1, arc1.centerPoint, arc1.endPoint,
                arc1.startPoint, -arc1.dir, tol, eps);

        //check if new arc was created without errors
        if (newErr)
        {
            *err |= newErr;
            return 0;
        }
    }

    //Determine if arcs are circles
    arc1IsCircle = (fabs(arc1.subtendedAngle) >= M_2PI);
    arc2IsCircle = (fabs(arc2.subtendedAngle) >= M_2PI);

    //Check if both arcs are circles
    if (arc1IsCircle && arc2IsCircle)
    {
        //Return true since at this point we know that they are concentric and have the same radius and are both circles
        commonShapes[0].this = (Arc*) malloc(sizeof(Arc));
        newErr |= createArc((Arc*) commonShapes[0].this, arc1.centerPoint,
                arc1.startPoint, arc1.endPoint, arc1.dir, tol, eps);
        commonShapes[0].type = ARC;

        //check if new arc was created without errors
        if (newErr)
        {
            *err |= newErr;
            return 0;
        }
        else
        {
            *shapeCount = 1;
            return 1;
        }
    }
    //If only one arc is a circle, return the other arc
    else if (arc1IsCircle)
    {
        commonShapes[0].this = (Arc*) malloc(sizeof(Arc));
        newErr |= createArc((Arc*) commonShapes[0].this, arc2.centerPoint,
                arc2.startPoint, arc2.endPoint, arc2.dir, tol, eps);
        commonShapes[0].type = ARC;

        //check if new arc was created without errors
        if (newErr)
        {
            *err |= newErr;
            return 0;
        }
        else
        {
            *shapeCount = 1;
            return 1;
        }
    }
    else if (arc2IsCircle)
    {
        commonShapes[0].this = (Arc*) malloc(sizeof(Arc));
        newErr |= createArc((Arc*) commonShapes[0].this, arc1.centerPoint,
                arc1.startPoint, arc1.endPoint, arc1.dir, tol, eps);
        commonShapes[0].type = ARC;

        //check if new arc was created without errors
        if (newErr)
        {
            *err |= newErr;
            return 0;
        }
        else
        {
            *shapeCount = 1;
            return 1;
        }
    }
    else
    {
        //Add the arc points to the array for sorting
        createPt(&arcPoints[0], newArc1.startPoint.latitude,
                newArc1.startPoint.longitude);
        createPt(&arcPoints[1], newArc1.endPoint.latitude,
                newArc1.endPoint.longitude);
        createPt(&arcPoints[2], arc2.startPoint.latitude,
                arc2.startPoint.longitude);
        createPt(&arcPoints[3], arc2.endPoint.latitude,
                arc2.endPoint.longitude);

        arcPointAzs[0] = newArc1.startAz;
        arcPointAzs[1] = newArc1.endAz;
        arcPointAzs[2] = arc2.startAz;
        arcPointAzs[3] = arc2.endAz;

        //Sort the points with respect to the arc direction (should be same for both arcs)
        newErr |= sortPtsByAz(arcPoints, 4, arcPointAzs, arc2.dir);

        if (newErr)
        {
            *err |= newErr;
            return 0;
        }

        //copy first sorted point to end
        createPt(&arcPoints[4], arcPoints[0].latitude,
                arcPoints[0].longitude);

        //copy first sorted point azimuth to end
        arcPointAzs[4] = arcPointAzs[0];

        //construct test arcs and determine if common
        for (i = 0; i < 4; i++)
        {

            commonStart = arcPoints[i];
            commonEnd = arcPoints[i + 1];

            //Check if points are the same
            if (ptsAreSame(commonStart, commonEnd, tol))
            {
                if (commonPoint == NULL)
                {
                    commonPoint = malloc(sizeof(LLPoint));
                    //Save the common end point
                    *commonPoint = commonEnd;
                    continue;//start next iteration
                }
                else
                {
                    //At least three points are very close to one another.  This implies at least one arc has identical end points.
                    *err |= SHAPE_NOT_DEFINED_ERR;
                    return 0;
                }
            }

            //Create the arc to be tested
            newErr |= createArc(&testArc, newArc1.centerPoint, commonStart,
                    commonEnd, arc2.dir, tol, eps);

            if (newErr)
            {
                *err |= newErr;
                return 0;
            }

            //Get the midpoint on the test arc
            //Determine the magnitude of the angle between the start and end true course
            newErr |= minSubtendedAngle(testArc.startAz, testArc.endAz,
                    &angleBtwn);

            if (newErr)
            {
                *err |= newErr;
                return 0;
            }

            //Determine the half angle
            halfAngleBtwn = angleBtwn / 2.0;

            //Determine the angle bisector course
            if (testArc.dir == CLOCKWISE) //arc direction is clockwise
            {
                angleBisectorCourse = fmod(testArc.startAz + halfAngleBtwn,
                        M_2PI); //modulo 2Pi (i.e. 360 degrees)
            }
            else
            {
                angleBisectorCourse = fmod(testArc.startAz - halfAngleBtwn,
                        M_2PI); //modulo 2Pi (i.e. 360 degrees)
            }

            //Determine the point an arc radius distance along the angle bisector course from the arc center point
            newErr |= direct(testArc.centerPoint, angleBisectorCourse,
                    testArc.radius, &testPoint, eps);
            commonMid = testPoint;

            if (newErr)
            {
                *err |= newErr;
                return 0;
            }

            //check if midpoint of test arc is on both arc1 and arc2
            midPtOnArc1 = ptIsOnArc(newArc1.centerPoint, newArc1.radius,
                    newArc1.startAz, newArc1.endAz, newArc1.dir, commonMid,
                    &newErr, tol, eps);
            midPtOnArc2 = ptIsOnArc(arc2.centerPoint, arc2.radius,
                    arc2.startAz, arc2.endAz, arc2.dir, commonMid, &newErr,
                    tol, eps);

            if (newErr)
            {
                *err |= newErr;
                return 0;
            }

            if (midPtOnArc1 && midPtOnArc2)
            {
                //test arc is common to both
                commonShapes[*shapeCount].this = (Arc*) malloc(sizeof(Arc));
                newErr |= createArc((Arc*) commonShapes[*shapeCount].this,
                        testArc.centerPoint, testArc.startPoint,
                        testArc.endPoint, testArc.dir, tol, eps);

                if (newErr)
                {
                    *err |= newErr;
                    return 0;
                }
                arcPtr = (Arc*)commonShapes[*shapeCount].this;
                if ( (commonPoint != NULL) && ptIsOnArc(arcPtr->centerPoint, arcPtr->radius, arcPtr->startAz, arcPtr->endAz, arcPtr->dir, *commonPoint, &newErr, tol, eps))
                  commonPoint = NULL;

                commonShapes[*shapeCount].type = ARC;
                (*shapeCount)++;
            }
            else if (commonPoint != NULL)
            {
                skipPoint = 0; 
                for (j = 0; j < *shapeCount; j++)
                {
                  if (commonShapes[j].type == ARC)
                  {
                    arcPtr = (Arc*)commonShapes[j].this;
                    if (ptIsOnArc(arcPtr->centerPoint, arcPtr->radius, arcPtr->startAz, arcPtr->endAz, arcPtr->dir, *commonPoint, &newErr, tol, eps))
                    {
                      commonPoint = NULL;
                      skipPoint = 1;
                      break;
                    }
                  }
                }
                if (skipPoint)
                  continue;
                //discrete point was found and is common to both arcs
                commonShapes[*shapeCount].this = (LLPoint*) malloc(
                        sizeof(LLPoint));
                createPt((LLPoint*) commonShapes[*shapeCount].this,
                        commonPoint->latitude, commonPoint->longitude);
                commonShapes[*shapeCount].type = LLPOINT;
                (*shapeCount)++;

                //clear the common point
                commonPoint = NULL;
            }

        } // for i

        if (*shapeCount > 1)
        {
          if ((commonShapes[0].type == LLPOINT) && (commonShapes[1].type == ARC))
          {
            commonPoint = (LLPoint*)commonShapes[0].this;
            arcPtr = (Arc*)commonShapes[1].this;
            if (ptIsOnArc(arcPtr->centerPoint, arcPtr->radius, arcPtr->startAz, arcPtr->endAz, arcPtr->dir, *commonPoint, &newErr, tol, eps))
            {
              for (j = 0; j < *shapeCount - 1; j++)
              {
                commonShapes[j].type = commonShapes[j+1].type;
                commonShapes[j].this = commonShapes[j+1].this;
              }
              (*shapeCount)--;
            }
          }
        } 

        return 1;
    }
}

ErrorSet bndryIntxExists(Boundary B1, Boundary B2, int* intersectionExists, double tol, double eps){

	int i,j,k;
	LLPointPair tempIntx;
	LLPointSet tempIntxSet = createPtSet();
	int numberOfIntx;
	Geodesic* geo1 =  NULL;
	Geodesic* geo2 =  NULL;
	Arc* arc1 =  NULL;
	Arc* arc2 =  NULL;
	Locus* locus1 =  NULL;
	Locus* locus2 =  NULL;
	Spiral* spiral1 =  NULL;
	Spiral* spiral2 =  NULL;
	ErrorSet err = 0;
	Shape* tempShape;

	// Assume no intersection
	*intersectionExists = 0;

	// Test all segment combinations for intersection

	for(i = 0; i < B1.length; i++){
		tempShape = B1.elements[i].this;
		switch(B1.elements[i].type){
			case GEODESIC:
			geo1 = (Geodesic*) tempShape;
			for(j = 0; j < B2.length; j++){
				tempShape = B2.elements[j].this;
				switch(B2.elements[j].type){
					case GEODESIC:
						geo2 = (Geodesic*) tempShape;
						err = geoIntx(geo1 -> startPoint, geo1 -> endPoint, SEGMENT, NULL, NULL,
									  geo2 -> startPoint, geo2 -> endPoint, SEGMENT, NULL, NULL,
									  &tempIntx[0], tol, eps); // Applies segment bounds
						if(err == SUCCESS){
							*intersectionExists = 1;
							return err;
						} else if(err == NO_INTERSECTION_ERR) {
							err = 0;
						} else {
							*intersectionExists = -1;
							return err;
						}
						break;

					case LOCUS:
						locus2 = (Locus*) tempShape;
						err = locusGeoIntx(geo1->startPoint, geo1->endPoint,
										   *locus2,
										   &tempIntx[0], tol, eps); // Geo bounds must be checked
						if(err == SUCCESS){
							if(ptIsOnGeo(geo1->startPoint,geo1->endPoint,tempIntx[0],SEGMENT,&err,tol,eps)){
								*intersectionExists = 1;
								return err;
							}
						} else if(err == NO_INTERSECTION_ERR) {
							err = 0;
						} else {
							*intersectionExists = -1;
							return err;
						}
						break;

					case ARC:
						arc2 = (Arc*) tempShape;
						err = geoArcIntx(geo1 -> startPoint, geo1 -> startAz,
								         arc2 -> centerPoint, arc2 -> radius,
								         tempIntx, &numberOfIntx, tol, eps); // Arc and geo bounds must be checked
						if(err == SUCCESS){
							for(k = 0; k < numberOfIntx; k++){
								if(ptIsOnArc(arc2->centerPoint,arc2->radius,arc2->startAz,arc2->endAz,arc2->dir,tempIntx[k],&err,tol,eps)
										&&
								   ptIsOnGeo(geo1->startPoint,geo1->endPoint,tempIntx[k],SEGMENT,&err,tol,eps)){
										if(numberOfIntx > 0){
											*intersectionExists = 1;
											return err;
										}
								}
							}
						} else if(err == NO_INTERSECTION_ERR) {
							err = 0;
						} else {
							*intersectionExists = -1;
							return err;
						}
						break;

					case SPIRAL:
						spiral2 = (Spiral*) tempShape;
						err = spiralGeoIntx(*spiral2,
											*geo1,
											&tempIntxSet, tol, eps);
						if(err == SUCCESS){
							if(tempIntxSet.length > 0){
								*intersectionExists = 1;
								return err;
							}
						} else if(err == NO_INTERSECTION_ERR) {
							err = 0;
						} else {
							*intersectionExists = -1;
							return err;
						}
						clearPtSet(&tempIntxSet);
						break;
			        default:
			        	break;
				}

			}
			break;

			case LOCUS:
				locus1 = (Locus*) tempShape;
				for(j = 0; j < B2.length; j++){

					tempShape = B2.elements[j].this;
					switch(B2.elements[j].type){
						case GEODESIC:
							geo2 = (Geodesic*) tempShape;
							err = locusGeoIntx(geo2->startPoint, geo2->endPoint,
											   *locus1,
											   &tempIntx[0], tol, eps); // Geo bounds must be checked
							if(err == SUCCESS){
								if(ptIsOnGeo(geo2->startPoint,geo2->endPoint,tempIntx[0],SEGMENT,&err,tol,eps)){
									*intersectionExists = 1;
									return err;
								}
							} else if(err == NO_INTERSECTION_ERR) {
								err = 0;
							} else {
								*intersectionExists = -1;
								return err;
							}
							break;

						case LOCUS:
							locus2 = (Locus*) tempShape;
							err = locusIntx(*locus1,
											*locus2,
											&tempIntx[0], tol, eps); // Applies segment bounds
							if(err == SUCCESS){
								*intersectionExists = 1;
								return err;
							} else if(err == NO_INTERSECTION_ERR) {
								err = 0;
							} else {
								*intersectionExists = -1;
								return err;
							}
							break;

						case ARC:
							arc2 = (Arc*) tempShape;
							err = locusArcIntx(*locus1,
												arc2->centerPoint, arc2->radius,
												tempIntx, &numberOfIntx, tol, eps); // Locus bounds are applied but Arc bounds must be checked
							if(err == SUCCESS){
								if(numberOfIntx > 0){
									for(k = 0; k < numberOfIntx; k++){
										if(ptIsOnArc(arc2->centerPoint,arc2->radius,arc2->startAz,arc2->endAz,arc2->dir,tempIntx[k],&err,tol,eps)){
											*intersectionExists = 1;
											return err;
										}
									}
								}
							} else if(err == NO_INTERSECTION_ERR) {
								err = 0;
							} else {
								*intersectionExists = -1;
								return err;
							}
							break;

						case SPIRAL:
							spiral2 = (Spiral*) tempShape;
							err = spiralLocusIntx(*spiral2,
												*locus1,
												&tempIntxSet, tol, eps);
							if(err == SUCCESS){
								if(tempIntxSet.length > 0){
									*intersectionExists = 1;
									return err;
								}
							} else if(err == NO_INTERSECTION_ERR) {
								err = 0;
							} else {
								*intersectionExists = -1;
								return err;
							}
							clearPtSet(&tempIntxSet);
							break;
				        default:
				        	break;
					}

				}
				break;

			case ARC:
				arc1 = (Arc*) tempShape;
				for(j = 0; j < B2.length; j++){

					tempShape = B2.elements[j].this;
					switch(B2.elements[j].type){
						case GEODESIC:
							geo2 = (Geodesic*) tempShape;
							err = geoArcIntx(geo2 -> startPoint, geo2 -> startAz,
									         arc1 -> centerPoint, arc1 -> radius,
									         tempIntx, &numberOfIntx, tol, eps); // Arc and geo bounds must be checked
							if(err == SUCCESS){
								if(numberOfIntx > 0){
									for(k = 0; k < numberOfIntx; k++){
										if(ptIsOnArc(arc1->centerPoint,arc1->radius,arc1->startAz,arc1->endAz,arc1->dir,tempIntx[k],&err,tol,eps)
												&&
										   ptIsOnGeo(geo2->startPoint,geo2->endPoint,tempIntx[k],SEGMENT,&err,tol,eps)){
												*intersectionExists = 1;
												return err;
										}
									}
								}
							} else if(err == NO_INTERSECTION_ERR) {
								err = 0;
							} else {
								*intersectionExists = -1;
								return err;
							}
							break;

						case LOCUS:
							locus2 = (Locus*) tempShape;
							err = locusArcIntx(*locus2,
												arc1->centerPoint, arc1->radius,
												tempIntx, &numberOfIntx, tol, eps); // Locus bounds are applied but Arc bounds must be checked
							if(err == SUCCESS){
								if(numberOfIntx > 0){
									for(k = 0; k < numberOfIntx; k++){
										if(ptIsOnArc(arc1->centerPoint,arc1->radius,arc1->startAz,arc1->endAz,arc1->dir,tempIntx[k],&err,tol,eps)){
											*intersectionExists = 1;
											return err;
										}
									}
								}
							} else if(err == NO_INTERSECTION_ERR) {
								err = 0;
							} else {
								*intersectionExists = -1;
								return err;
							}
							break;

						case ARC:
							arc2 = (Arc*) tempShape;
							err = arcIntx(arc1 -> centerPoint, arc1 -> radius,
									      arc2 -> centerPoint, arc2 -> radius,
									      tempIntx, &numberOfIntx, tol, eps); // Both arc bounds must be checked
							if(err == SUCCESS){
								if(numberOfIntx > 0){
									for(k = 0; k < numberOfIntx; k++){
										if(ptIsOnArc(arc1->centerPoint,arc1->radius,arc1->startAz,arc1->endAz,arc1->dir,tempIntx[k],&err,tol,eps)
												&&
										   ptIsOnArc(arc2->centerPoint,arc2->radius,arc2->startAz,arc2->endAz,arc2->dir,tempIntx[k],&err,tol,eps)){
												*intersectionExists = 1;
												return err;
										}
									}
								}
							} else if(err == NO_INTERSECTION_ERR) {
								err = 0;
							} else {
								*intersectionExists = -1;
								return err;
							}
							break;

						case SPIRAL:
							spiral2 = (Spiral*) tempShape;
							err = spiralArcIntx(*spiral2,
												*arc1,
												&tempIntxSet, tol, eps);
							if(err == SUCCESS){
								if(tempIntxSet.length > 0){
									*intersectionExists = 1;
									return err;
								}
							} else if(err == NO_INTERSECTION_ERR) {
								err = 0;
							} else {
								*intersectionExists = -1;
								return err;
							}
							clearPtSet(&tempIntxSet);
							break;
				        default:
				        	break;

					}

				}
				break;

			case SPIRAL:
				spiral1 = (Spiral*) tempShape;
				for(j = 0; j < B2.length; j++){
					tempShape = B2.elements[j].this;
					switch(B2.elements[j].type){
						case GEODESIC:
							geo2 = (Geodesic*) tempShape;
							err = spiralGeoIntx(*spiral1,
												*geo2,
												&tempIntxSet, tol, eps);
							if(err == SUCCESS){
								if(tempIntxSet.length > 0){
									*intersectionExists = 1;
									return err;
								}
							} else if(err == NO_INTERSECTION_ERR) {
								err = 0;
							} else {
								*intersectionExists = -1;
								return err;
							}
							clearPtSet(&tempIntxSet);
							break;

						case LOCUS:
							locus2 = (Locus*) tempShape;
							err = spiralLocusIntx(*spiral1,
												  *locus2,
												  &tempIntxSet, tol, eps);
							if(err == SUCCESS){
								if(tempIntxSet.length > 0){
									*intersectionExists = 1;
									return err;
								}
							} else if(err == NO_INTERSECTION_ERR) {
								err = 0;
							} else {
								*intersectionExists = -1;
								return err;
							}
							clearPtSet(&tempIntxSet);
							break;

						case ARC:
							arc2 = (Arc*) tempShape;
							err = spiralArcIntx(*spiral1,
												*arc2,
												&tempIntxSet, tol, eps);
							if(err == SUCCESS){
								if(tempIntxSet.length > 0){
									*intersectionExists = 1;
									return err;
								}
							} else if(err == NO_INTERSECTION_ERR) {
								err = 0;
							} else {
								*intersectionExists = -1;
								return err;
							}
							clearPtSet(&tempIntxSet);
							break;

						case SPIRAL:
							spiral2 = (Spiral*) tempShape;
							err = spiralIntx(*spiral1,
											 *spiral2,
										     &tempIntxSet, tol, eps);
							if(err == SUCCESS){
								if(tempIntxSet.length > 0){
									*intersectionExists = 1;
									return err;
								}
							} else if(err == NO_INTERSECTION_ERR) {
								err = 0;
							} else {
								*intersectionExists = -1;
								return err;
							}
							clearPtSet(&tempIntxSet);
							break;
				        default:
				        	break;
					}

				}
				break;
	        default:
	        	break;
		}
	}

	tempShape = B1.elements[0].this;
	switch(B1.elements[0].type){
		case GEODESIC:
			geo1 = (Geodesic*) tempShape;
			if(ptIsInsideBndry(B2, geo1->startPoint, &err, tol, eps)){
				*intersectionExists = 2;
				return err;
			}
			break;

		case LOCUS:
			locus1 = (Locus*) tempShape;
			if(ptIsInsideBndry(B2, locus1->locusStart, &err, tol, eps)){
				*intersectionExists = 2;
				return err;
			}
			break;

		case ARC:
			arc1 = (Arc*) tempShape;
			if(ptIsInsideBndry(B2, arc1->startPoint, &err, tol, eps)){
				*intersectionExists = 2;
				return err;
			}
			break;

		case SPIRAL:
			spiral1 = (Spiral*) tempShape;
			if(ptIsInsideBndry(B2, spiral1->startPoint, &err, tol, eps)){
				*intersectionExists = 2;
				return err;
			}
			break;
        default:
        	break;
	}

	tempShape = B2.elements[0].this;
	switch(B2.elements[0].type){
		case GEODESIC:
			geo2 = (Geodesic*) tempShape;
			if(ptIsInsideBndry(B1, geo2->startPoint, &err, tol, eps)){
				*intersectionExists = 3;
				return err;
			}
			break;

		case LOCUS:
			locus2 = (Locus*) tempShape;
			if(ptIsInsideBndry(B1, locus2->locusStart, &err, tol, eps)){
				*intersectionExists = 3;
				return err;
			}
			break;

		case ARC:
			arc2 = (Arc*) tempShape;
			if(ptIsInsideBndry(B1, arc2->startPoint, &err, tol, eps)){
				*intersectionExists = 3;
				return err;
			}
			break;

		case SPIRAL:
			spiral2 = (Spiral*) tempShape;
			if(ptIsInsideBndry(B1, spiral2->startPoint, &err, tol, eps)){
				*intersectionExists = 3;
				return err;
			}
			break;
        default:
        	break;

	}

	return err;
}

ErrorSet projectToBndry(Boundary b, LLPoint pt, LLPoint* perpPt, double* perpDist, double tol, double eps) {

	ErrorSet err = 0;
	int i, flag = 0;
	double az12, az21, dist, minDist = 1000;
	LLPoint tempPt, p;
	Shape* thisShape;
	Geodesic* thisGeo = NULL;
	Locus* thisLocus = NULL;
	Arc* thisArc = NULL;
	Spiral* thisSpiral = NULL;

    for (i = 0; i < b.length; i++)
    {
		thisShape = b.elements[i].this;
		if (thisShape == NULL)
		{
		return SHAPE_NOT_DEFINED_ERR;
		}
		ShapeType thisType = b.elements[i].type;
		switch (thisType) {
		case ARC:
			thisArc = (Arc*) thisShape;
			err |= inverse(thisArc->centerPoint, pt, &az12, &az21, &dist, eps);
			err |= direct(thisArc->centerPoint, az12, thisArc->radius, &tempPt, eps);
			dist = dist - thisArc->radius;
			//TODO:  Add logic to see if point is within arc extent
			if ((dist < minDist) && (ptIsOnArc(thisArc->centerPoint, thisArc->radius, thisArc->startAz, thisArc->endAz, thisArc->dir, tempPt, &err, tol, eps))) {
				minDist = dist;
				*perpPt = tempPt;
			}
			break;
        case GEODESIC:
            thisGeo = (Geodesic*) thisShape;
            err |= projectToGeo(thisGeo->startPoint, thisGeo->startAz, pt, &tempPt, &az12, &dist, tol, eps);
            if (err == 0) {
    			if ((dist < minDist) && (ptIsOnGeo(thisGeo->startPoint, thisGeo->endPoint, tempPt, SEGMENT, &err, tol, eps))) {
    				minDist = dist;
    				*perpPt = tempPt;
    			}
            } else {
            	err = 0;
            }
            break;
        case LOCUS:
            thisLocus = (Locus*) thisShape;
            err |= projectToLocus(*thisLocus, pt, &tempPt, &az12, &dist, tol, eps);
            LLPoint geoPt;
            if (err == 0) {
                if ((dist < minDist) && ptIsOnLocus(*thisLocus, tempPt, &geoPt, &err, tol, eps)) {
                	minDist = dist;
                	*perpPt = tempPt;
                }
            } else {
              	err = 0;
            }
            break;
        case SPIRAL:
        	thisSpiral = (Spiral*) thisShape;
        	err |= projectToSpiral(*thisSpiral, pt, &tempPt, tol, eps);
        	err |= invDist(pt, tempPt, &dist, eps);
            if (err == 0) {
                if ((dist < minDist) && ptIsOnSpiral(*thisSpiral,pt, tol, eps)) {
                	minDist = dist;
                	*perpPt = tempPt;
                }
            } else {
               	err = 0;
            }
            break;
        default:
        	break;
        } //switch thisType
    } //for i

    if (minDist == 1000) {
    	err = NO_PROJECTED_POINT_ERR;
    }

    return err;
}
