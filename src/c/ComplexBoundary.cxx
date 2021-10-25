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
#include "ComplexBoundary.h"



/** Initialize a ComplexBoundary struct and return it. */
ComplexBoundary newComplexBoundary()
{
    ComplexBoundary b = { 0, NULL, COMPLEX_BOUNDARY_ARRAY_INCREMENT };
    b.elements = (Boundary**) realloc(b.elements, COMPLEX_BOUNDARY_ARRAY_INCREMENT * sizeof(Boundary*));
    return b;
}

/** Frees all of the dynamically allocated memory associated with a
 * ComplexBoundary structure.  Resulting complex boundary has length == 0 and elements pointer
 * set to NULL.
 * @param c Pointer to a ComplexBoundary structure
 * @returns Nothing
 */
void clearComplexBoundary(ComplexBoundary* c)
{
    int i = 0;
    Boundary* b;

    if ((NULL == c) || (0 == c->length) || (NULL == c->elements))
    {
        return;
    }

    for (i = 0; i < c->length; i++)
    {
    	b = (Boundary*) (c->elements[i]);
		clearBndry(b);
    }

    if (NULL != c->elements)
    {
        /* Deallocate memory for Boundary array */
        free(c->elements);
        c->elements = NULL;
    }

    c->length = 0;

}

/** Increase size of array by n elements
 *
 */
ErrorSet growComplexBoundaryArrayByN(ComplexBoundary* b, int n)
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
            b->maxLength += (d / COMPLEX_BOUNDARY_ARRAY_INCREMENT + 1)*COMPLEX_BOUNDARY_ARRAY_INCREMENT;
            b->elements = (Boundary**) realloc(b->elements, b->maxLength * sizeof(Boundary*));
        }

        if (NULL == b->elements)
            err |= MALLOC_ERR;
        else
            b->length += n;
    }
    return err;

}

/** Increase size of array by 1 element */
ErrorSet growComplexBoundaryArray(ComplexBoundary* b)
{
    return growComplexBoundaryArrayByN(b, 1);
}


ErrorSet addElementToComplexBoundary(ComplexBoundary* c, Boundary* element)
{
    int n;
    ErrorSet err = 0;

    if (NULL == c || NULL == element)
    {
    	err |= NO_MEMORY_ALLOCATED_ERR;
        return err;
    }

    err |= growComplexBoundaryArray(c);

    if (err == 0)
    {
        n = c->length;

		c->elements[n - 1] = (element);

    }

    return err;

}


/*
 * Return a list of common areas for a list of circles and a simple connected Boundary structure
 *
 * Find the boundary that defines the area that is common to both the circle
 * and to the input boundary.  This boundary is added as an element to the
 * commonAreas output.
 *
 * Note - If any circle in the list does not intersect the boundary, an entry is still
 * 		  created in the commonAreas list.  This entry will be a new instance of a
 * 		  boundary struct with length zero.
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
 * comprise the boundary and the number of cirlces  to be evaluated.
 * More boundary shapes implies more time to process as does more circles.
 *
 *
 * inputs:
 * 		boundary = Boundary structure to be evaluated against
 * 		circles = Circle array to be evaluated for intersections with boundary
 *      circleSize = the dimension of the circle array
 * 		commonAreas = Pointer to the complex boundary that holds a list of common areas.
 * 		tol = distance criteria (in nmi)
 * 		eps = Vincenty forward/inverse algorithm convergence criteria
 */
ErrorSet boundaryCircleListIntersections(Boundary boundary, Arc circles[], int circlesSize,
                                    ComplexBoundary* commonAreas, double tol, double eps)
{
	ErrorSet err = 0;
	Boundary* common = NULL;
	int i = 0;

	for (i = 0; i < circlesSize; i++)
	{
		common = (Boundary*) malloc(sizeof(Boundary));
		*common = createBndry();
		err |= bndryCircleIntx(boundary, circles[i], common, tol, eps);
		err |= addElementToComplexBoundary(commonAreas, common);
	}

	return err;
}

/*
 * Partition a simple connected ordered boundary using the geodesics that lie on the meridian
 * segments of a longitude integer value.  Please see the Boundary class for a description of
 * simple boundaries and for a description of ordered boundaries.
 *
 * A simple boundary implies one that contains no holes in the boundary area.  This also
 * means that the boundary area is not originally partitioned (contains sub areas) in any way.
 *
 *                    not partitioned          partitioned
 *
 *                          S_0                     X
 *                       ---------             -----x-----
 *                       |       |             |    x    |
 *                   S_3 |       | S_1         | A  x  B |
 *                       |       |             |    x    |
 *                       ---------             -----x-----
 *                          S_2                     X
 *
 * In the above partition diagram, the sequence of x's (capital and lower case) represent the meridian for a
 * calculated integer longitude value.  The geodesic that lies on this meridian path partitions
 * the simple boundary into two partitions (A and B).  The sequence of lower case x's represent
 * the geodesic segment that partitions the boundary such that this segment's start/end points lie on the
 * original simple boundary.
 *
 * Let the integer longitude domain be equal to [-180, 180] degrees.
 *
 * This algorithm should also handle whether the start/end shape connections are traversed on the ordered boundary in
 * either a clockwise or counter-clockwise direction.
 *
 * The algorithm for this function is as follows:
 * 		Step 1 - Find the minimum and maximum integer longitude from the set of longitude values used
 *               to define the start/end points for each shape in the simple boundary.  Note that these values
 *               are the least integer upper bound and greatest integer lower bound and may not necessarily be
 *               actual longitudinal values found in any shape's start/end point.
 *               For example suppose the actual start/end points for the shapes in the above unpartitioned diagram
 *               have the following set L of unique longitude values.
 *                           L = {-73.235, -74.299, -73.211}
 *               Note that this set has no integer values and that the min = -74.299 and the max = -73.211.
 *               These are lower and upper bounds for the set but they are not the least integer upper bound and
 *               greatest integer lower bound.
 *               Since -75 is an integer, -75 < -74.299, and -74 >= x where x is any integer less than or equal to
 *               -74.299 (not a typo) then this algorithm determines that -75 is the greatest integer lower bound.
 *               Similarly, since -73 is an integer, -73.211 < -73, and -73 <= y where y is any integer greater than
 *               or equal to -73.211 (not a typo) then this algorithm determines that -73 is the least integer upper bound.
 *
 *      Step 2 - Create a set of infinite geodesics that lie on the meridians that correspond to an integer
 *               set of longitude values bounded by the upper/lower bounds found in step 1.  The set of integer longitude
 *               values are uniformly separated by one degree.
 *               Continuing with the above example, using the lower/upper bound value of -75 and -73 (respectively), the
 *               integer set {-75, -74, -73} of longitude values is identified.  At which point three infinite geodesics
 *               that lie on the meridians that correspond to each integer longitude value from this set are created.
 *      Step 3 - Use the set of infinite geodesics found in step 2, to begin partitioning an ordered simple boundary.
 *      		 Let the infinite geodesic that corresponds to the greatest integer lower bound be called the
 *               least infinite geodesic.
 *               Let the infinite geodesic that corresponds to the least integer upper bound be called the greatest
 *               infinite geodesic.
 *               Compute the number of partitions to be equal to one minus the order of the integer longitude set.
 *               Identify the two infinite geodesics that bound each partition.
 *               Process the shapes that comprise each partition one at a time starting with the partition bounded by the
 *               least infinite geodesic and the next infinite geodesic from that set.  Do this finding the intersections
 *               of the ordered simple boundary with the two bounding infinite geodesics.  Using properties that result
 *               from the fact that the input boundary is simple and is well-ordered, begin truncating input boundary
 *               shapes and creating segments from the bounding infinite geodesics such that the truncated shapes
 *               lie within the bounding infinite longitudinal geodesics.  Furthermore, since the segments from the
 *               infinite bounding geodesics lie on the infinite geodesic then they must also be part of the partition.
 *
 *               Continuing with the above example, the order of the integer longitude set is 3.  Therefore,
 *               the number of partitions is 2 = 3 - 1.
 *
 *               The least infinite geodesic in this example is the infinite geodesic that corresponds to -75 degrees.
 *
 *               The greatest infinite geodesic in this example is the infinite geodesic that corresponds to -73 degrees.
 *
 *               Let partition A be bounded by the infinite geodesics that correspond to -75 and -74 degrees.
 *
 *               Let partition B be bounded by the infinite geodesics that correspond to -74 and -73 degrees.
 *
 *               Note that since partition partition A is bounded by the infinite geodesic that corresponds to -75 degrees
 *               then it must be that partition A is bounded by the least infinite geodesic.  Therefore partition A must be
 *               the first partition to be processed.
 *
 *               Note that since partition  partition B is bounded by the infinite geodesic that corresponds to -73 degrees
 *               then it must be that partition B is bounded by the greatest infinite geodesic.  Therefore partition B must be
 *               the last partition to be processed.
 *
 *               Observe that in this case since number of partitions is two, we have that the first partition created is partition
 *               A and the second (and last) partition created is partition B.  If more than two partitions were to be created
 *               then they would continue to be bounded in a similar logical order as described above with some number of parititions
 *               created before partition B was created.
 *       Step 4 - Once a partition has been created, add the partition to the returned complex boundary.  Continue creating
 *       		partitions in similar manner until all the computed partitions have been created.
 *
 * inputs:
 * 		boundary = The simple ordered boundary to be partitioned
 * 		tol = distance criteria (in nmi)
 * 		eps = Vincenty forward/inverse algorithm convergence criteria
 * outputs:
 *      complexBoundary = A pointer to a collection of boundaries.  Each element from this collection represents the boundary that
 *      contains the area from the input boundary which has been partitioned by the meridians that correspond to integer
 *      longitude values.
 */
ErrorSet longitudinallyPartitionBoundary(Boundary orderedBoundary, ComplexBoundary* complexBoundary, double tol, double eps)
{
	const double RAD2DEG = 180.0 / M_PI;
	const double DEG2RAD = M_PI / 180.0;
	ErrorSet err = 0;
	int longitudeSetSize;//the dimension of the array that holds the set of longitudes from the shapes that comprise the ordered boundary
	double* orderedBoundaryLongitudes = NULL;//the array that holds the set of longitudes from the shapes that comprise the ordered boundary
	double longitudeMax;//the maximum longitude from the array that holds the set of longitudes from the shapes that comprise the ordered boundary
	double longitudeMin;//the minimum longitude from the array that holds the set of longitudes from the shapes that comprise the ordered boundary
	int longitudeLowerBound;//the integer lower bound from the array that holds the set of longitudes from the shapes that comprise the ordered boundary
	int longitudeUpperBound;//the integer upper bound from the array that holds the set of longitudes from the shapes that comprise the ordered boundary
	int i = 0;
	int j = 0;
    Shape* segment = NULL; //The shape currently being processed
    Arc* segmentArc = NULL; //Temp object to hold segment recast as an Arc
    Geodesic* segmentGeo = NULL; //Temp object to hold segment recast as an Geodesic
    Locus* segmentLocus = NULL; //Temp object to hold segment recast as an Locus
    LLPoint startPt, endPt;
    Geodesic** infLongitudeSet = NULL;//pointer array of infinite geodesics that lie on the meridians that correspond to integer longitudes bounded by the upper/lower integer longitude bounds
    int infLongitudeSetSize;//the size of the pointer array that holds the infinite geodesics on the meridians
    LLPoint* infLongStartPt = NULL;//start point for an infinite geodesic on a meridian
    LLPoint* infLongEndPt = NULL;//end point for an infinite geodesic on a meridian
    Geodesic* infLongGeo = NULL;//infinite geodesic on a meridian
    double tempLongitude;
    double tempLatitude;
    Boundary* partition = NULL;//a partitioned boundary created from infinite longitudinal geodesics and the shapes from the ordered boundary
    double longitudeLowerBoundForPartition;//integer lower bound longitude value for the partition being processed
    double longitudeUpperBoundForPartition;//integer upper bound longitude value for the partition being processed
    int startPtInside = 1;//1 = true if the start point of the ordered boundary shape lies between the longitudes for the infinite geodesics that bound the partition being process, 0 = false otherwise
    int endPtInside = 1;//1 = true if the end point of the ordered boundary shape lies between the longitudes for the infinite geodesics that bound the partition being process, 0 = false otherwise
    Geodesic lowerInfLongGeo;//an infinite geodesic from the set of all infinite geodesics used to act as a "lower" geometric bound for the partition being processed.
    Geodesic upperInfLongGeo;//an infinite geodesic from the set of all infinite geodesics used to act as a "upper" geometric bound for the partition being processed.
    Geodesic* lowerLongSegment = NULL;//the "lower" infinite geodesic segment that lies "in" the area defined by the ordered boundary, which comprises a segment of the partition being processed
    Geodesic* upperLongSegment = NULL;//the "upper" infinite geodesic segment that lies "in" the area defined by the ordered boundary, which comprises a segment of the partition being processed
    int lowerLongSegmentFound = 0;//1 = true if an intersection has been found with the "lower" infinite geodesic as the algorithm traverses the ordered boundary, 0 = false otherwise
    int upperLongSegmentFound = 0;//1 = true if an intersection has been found with the "upper" infinite geodesic as the algorithm traverses the ordered boundary, 0 = false otherwise
    LLPoint lowerLongitudinalIntx;//the intersection at the "lower" infinite geodesic with a shape from the ordered boundary
    LLPoint upperLongitudinalIntx;//the intersection at the "upper" infinite geodesic with a shape from the ordered boundary
    ErrorSet lowerLongError, upperLongError;//local error code use for processing intersections with lower/upper infinite geodesics while processing a partition
    Geodesic tempGeo;
    Locus tempLocus;
    Arc tempArc;
    LLPoint tempPt;
    LLPoint tempPt2;
    LLPoint* tempLowerLongSegStartPt;//holds the first intersection with the lower infinite geodesic until a second intersection with the lower geodesic is found
    LLPoint* tempUpperLongSegStartPt;//holds and second intersection with the upper infinite geodesic until a second intersection with the upper geodesic is found
    double tempCrs1, tempCrs2, tempDist1, tempDist2;



	//find longitudinal max/min
    longitudeSetSize = 2 * orderedBoundary.length;
	orderedBoundaryLongitudes = (double*) realloc(orderedBoundaryLongitudes, longitudeSetSize * sizeof(double));
	for(i = 0; i < orderedBoundary.length; i++)
	{
		segment = orderedBoundary.elements[i].this_shape;
		switch(orderedBoundary.elements[i].type)
		{
			case GEODESIC:
	            segmentGeo = (Geodesic*) segment;
				startPt = segmentGeo->startPoint;
				endPt = segmentGeo->endPoint;
				break;

			case LOCUS:
	            segmentLocus = (Locus*) segment;
				startPt = segmentLocus->locusStart;
				endPt = segmentLocus->locusEnd;
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
				//TODO do something
				break;
			default:
				err |= SHAPE_NOT_DEFINED_ERR;
		}
		orderedBoundaryLongitudes[2*i] = startPt.longitude;
		orderedBoundaryLongitudes[2*i + 1] = endPt.longitude;
		tempLatitude = startPt.latitude;//get any latitude, we just need one
	}
	err |= findSetMaxAndMin(orderedBoundaryLongitudes, longitudeSetSize, &longitudeMax, &longitudeMin);
	if (err)
	{
		free(orderedBoundaryLongitudes);
		orderedBoundaryLongitudes = NULL;
		return err;
	}
	//find the integer longitudinal extrema
	longitudeLowerBound = (int)floor(longitudeMin*RAD2DEG);//using degrees here
	longitudeUpperBound = (int)ceil(longitudeMax*RAD2DEG);//using degrees here

	//create set of "infinite" longitudinal lines with respect to integer longitudinal extrema
	infLongitudeSetSize = longitudeUpperBound - longitudeLowerBound + 1;
	infLongitudeSet = (Geodesic**) malloc(infLongitudeSetSize * sizeof(Geodesic*));
	for(i = 0; i < infLongitudeSetSize; i++)
	{
		tempLongitude = (longitudeLowerBound + i)*DEG2RAD;//degrees converted to radians

		infLongStartPt = (LLPoint*) malloc(sizeof(LLPoint));
		err |= createPt(infLongStartPt, tempLatitude, tempLongitude);

		infLongEndPt = (LLPoint*) malloc(sizeof(LLPoint));
		err |= createPt(infLongEndPt, tempLatitude + (1*DEG2RAD), tempLongitude);//degrees converted to radians

		infLongGeo = (Geodesic*) malloc(sizeof(Geodesic));
		err |= createGeo(infLongGeo, *infLongStartPt, *infLongEndPt, INFINITE, eps);

		infLongitudeSet[i] = infLongGeo;
	}

	//start partitioning ordered boundary
	for(i = 0; i < (infLongitudeSetSize - 1); i++)
	{
		//set the longitudinal bounds for the current partition
		longitudeLowerBoundForPartition = (longitudeLowerBound + i)*DEG2RAD;//degrees converted to radians
		longitudeUpperBoundForPartition = (longitudeLowerBound + i + 1)*DEG2RAD;//degrees converted to radians

		//get the infinite longitudinal geodesics to use for partitioning the ordered boundary
		lowerInfLongGeo = *(infLongitudeSet[i]);
		upperInfLongGeo = *(infLongitudeSet[i + 1]);

		//create instance of boundary partition
		partition = (Boundary*) malloc(sizeof(Boundary));
		(*partition) = createBndry();

		//populate boundary partition with respect to the current infinite longitudinal geodesics
		for(j = 0; j < orderedBoundary.length; j++)
		{
			segment = orderedBoundary.elements[j].this_shape;
			switch(orderedBoundary.elements[j].type)
			{
				case GEODESIC:
		            segmentGeo = (Geodesic*) segment;
					startPt = segmentGeo->startPoint;
					endPt = segmentGeo->endPoint;

					//check if start/end points lie between the longitudinal extrema for this boundary partition
					startPtInside = (longitudeLowerBoundForPartition <= startPt.longitude && startPt.longitude <= longitudeUpperBoundForPartition);
					endPtInside = (longitudeLowerBoundForPartition <= endPt.longitude && endPt.longitude <= longitudeUpperBoundForPartition);

					//find infinite longitudinal intersections
					lowerLongError = geoIntx(segmentGeo->startPoint, segmentGeo->endPoint, segmentGeo->lineType, &tempCrs1, &tempDist1,
							lowerInfLongGeo.startPoint, lowerInfLongGeo.endPoint, lowerInfLongGeo.lineType, &tempCrs2, &tempDist2, &lowerLongitudinalIntx, tol, eps);
					upperLongError = geoIntx(segmentGeo->startPoint, segmentGeo->endPoint, segmentGeo->lineType, &tempCrs1, &tempDist1,
							upperInfLongGeo.startPoint, upperInfLongGeo.endPoint, upperInfLongGeo.lineType, &tempCrs2, &tempDist2, &upperLongitudinalIntx, tol, eps);

					//add masked intersection errors to master error log for this function
					err |= getMaskedError(lowerLongError, getMask(0, 0, 1, 0, 0, 0, 0));
					err |= getMaskedError(upperLongError, getMask(0, 0, 1, 0, 0, 0, 0));

					//process the shape
					if (startPtInside && endPtInside)
					{
						//add shape to partition, no truncation needed
						err |= addElementToBndry(partition, segment, orderedBoundary.elements[j].type);

					}
					else if (startPtInside && !endPtInside)
					{
						//shape should intersect only one infinite longitudinal geodesic

						//process shape
						if ( !(lowerLongError & NO_INTERSECTION_ERR) )//using bitwise logical AND operation
						{
							//truncate shape at intersection
							err |= createGeo(&tempGeo, segmentGeo->startPoint, lowerLongitudinalIntx, SEGMENT, eps);

							//add shape to boundary
							err |= addElementToBndry(partition, &tempGeo, orderedBoundary.elements[j].type);

							//process longitudinal segment
							if (!lowerLongSegmentFound)
							{
								//save the first intersection point for the longitudinal segment that intersects at the current shape
								tempLowerLongSegStartPt = (LLPoint*) malloc(sizeof(LLPoint));
								err |= createPt(tempLowerLongSegStartPt, lowerLongitudinalIntx.latitude, lowerLongitudinalIntx.longitude);

								lowerLongSegmentFound = 1;
							} else {

								//create instance of longitudinal segment
								lowerLongSegment = (Geodesic*) malloc(sizeof(Geodesic));
								err |= createGeo(lowerLongSegment, *tempLowerLongSegStartPt, lowerLongitudinalIntx, SEGMENT, eps);

								//add longitudinal segment to boundary
								err |= addElementToBndry(partition, lowerLongSegment, GEODESIC);

								//reset longitudinal segment switch
								lowerLongSegment = 0;

								//free resources since the addElementToBndry method generates copies in memory
								free(tempLowerLongSegStartPt);
								tempLowerLongSegStartPt = NULL;
								free(lowerLongSegment);
								lowerLongSegment = NULL;

							}

						}
						else if ( !(upperLongError & NO_INTERSECTION_ERR) )//using bitwise logical AND operation
						{
							//truncate shape at intersection
							err |= createGeo(&tempGeo, segmentGeo->startPoint, upperLongitudinalIntx, SEGMENT, eps);

							//add shape to boundary
							err |= addElementToBndry(partition, &tempGeo, orderedBoundary.elements[j].type);

							//process longitudinal segment
							if (!upperLongSegmentFound)
							{
								//save the first intersection point for the longitudinal segment that intersects at the current shape
								tempUpperLongSegStartPt = (LLPoint*) malloc(sizeof(LLPoint));
								err |= createPt(tempUpperLongSegStartPt, upperLongitudinalIntx.latitude, upperLongitudinalIntx.longitude);

								upperLongSegmentFound = 1;
							} else {

								//create instance of longitudinal segment
								upperLongSegment = (Geodesic*) malloc(sizeof(Geodesic));
								err |= createGeo(upperLongSegment, *tempUpperLongSegStartPt, upperLongitudinalIntx, SEGMENT, eps);

								//add longitudinal segment to boundary
								err |= addElementToBndry(partition, upperLongSegment, GEODESIC);

								//reset longitudinal segment switch
								upperLongSegment = 0;

								//free resources since the addElementToBndry method generates copies in memory
								free(tempUpperLongSegStartPt);
								tempUpperLongSegStartPt = NULL;
								free(upperLongSegment);
								upperLongSegment = NULL;

							}
						}
						else
						{
							//intersection exepected but not found
							err |= UNEXPECTED_ERR;
							return err;
						}

					}
					else if (!startPtInside && endPtInside)
					{

						//shape should intersect only one infinite longitudinal geodesic

						//process shape
						if ( !(lowerLongError & NO_INTERSECTION_ERR) )//using bitwise logical AND operation
						{
							//truncate shape at intersection
							err |= createGeo(&tempGeo, lowerLongitudinalIntx, segmentGeo->endPoint, SEGMENT, eps);

							//add shape to boundary
							err |= addElementToBndry(partition, &tempGeo, orderedBoundary.elements[j].type);

							//process longitudinal segment
							if (!lowerLongSegmentFound)
							{
								//save the first intersection point for the longitudinal segment that intersects at the current shape
								tempLowerLongSegStartPt = (LLPoint*) malloc(sizeof(LLPoint));
								err |= createPt(tempLowerLongSegStartPt, lowerLongitudinalIntx.latitude, lowerLongitudinalIntx.longitude);

								lowerLongSegmentFound = 1;
							} else {

								//create instance of longitudinal segment
								lowerLongSegment = (Geodesic*) malloc(sizeof(Geodesic));
								err |= createGeo(lowerLongSegment, *tempLowerLongSegStartPt, lowerLongitudinalIntx, SEGMENT, eps);

								//add longitudinal segment to boundary
								err |= addElementToBndry(partition, lowerLongSegment, GEODESIC);

								//reset longitudinal segment switch
								lowerLongSegment = 0;

								//free resources since the addElementToBndry method generates copies in memory
								free(tempLowerLongSegStartPt);
								tempLowerLongSegStartPt = NULL;
								free(lowerLongSegment);
								lowerLongSegment = NULL;

							}

						}
						else if ( !(upperLongError & NO_INTERSECTION_ERR) )//using bitwise logical AND operation
						{
							//truncate shape at intersection
							err |= createGeo(&tempGeo, upperLongitudinalIntx, segmentGeo->endPoint, SEGMENT, eps);

							//add shape to boundary
							err |= addElementToBndry(partition, &tempGeo, orderedBoundary.elements[j].type);

							//process longitudinal segment
							if (!upperLongSegmentFound)
							{
								//save the first intersection point for the longitudinal segment that intersects at the current shape
								tempUpperLongSegStartPt = (LLPoint*) malloc(sizeof(LLPoint));
								err |= createPt(tempUpperLongSegStartPt, upperLongitudinalIntx.latitude, upperLongitudinalIntx.longitude);

								upperLongSegmentFound = 1;
							} else {

								//create instance of longitudinal segment
								upperLongSegment = (Geodesic*) malloc(sizeof(Geodesic));
								err |= createGeo(upperLongSegment, *tempUpperLongSegStartPt, upperLongitudinalIntx, SEGMENT, eps);

								//add longitudinal segment to boundary
								err |= addElementToBndry(partition, upperLongSegment, GEODESIC);

								//reset longitudinal segment switch
								upperLongSegment = 0;

								//free resources since the addElementToBndry method generates copies in memory
								free(tempUpperLongSegStartPt);
								tempUpperLongSegStartPt = NULL;
								free(upperLongSegment);
								upperLongSegment = NULL;

							}
						}
						else
						{
							//intersection exepected but not found
							err |= UNEXPECTED_ERR;
							return err;
						}

					}
					else if (!startPtInside && !endPtInside)
					{
						//check if shape straddles the longitudinal boundary, i.e. two intersections expected
						if ( !(lowerLongError & NO_INTERSECTION_ERR) && !(upperLongError & NO_INTERSECTION_ERR) )
						{

							//TODO handle the two cases: 1. shape.startPt closer to lowerLongitudinal bounding geo; shape.startPt closer to upperLongitudinal bounding geo
							//truncate shape at intersection
							err |= createGeo(&tempGeo, lowerLongitudinalIntx, upperLongitudinalIntx, SEGMENT, eps);

							//add shape to boundary
							err |= addElementToBndry(partition, &tempGeo, orderedBoundary.elements[j].type);

							//process longitudinal segments
							if (!lowerLongSegmentFound)
							{
								//save the first intersection point for the longitudinal segment that intersects at the current shape
								tempLowerLongSegStartPt = (LLPoint*) malloc(sizeof(LLPoint));
								err |= createPt(tempLowerLongSegStartPt, lowerLongitudinalIntx.latitude, lowerLongitudinalIntx.longitude);

								lowerLongSegmentFound = 1;
							} else {

								//create instance of longitudinal segment
								lowerLongSegment = (Geodesic*) malloc(sizeof(Geodesic));
								err |= createGeo(lowerLongSegment, *tempLowerLongSegStartPt, lowerLongitudinalIntx, SEGMENT, eps);

								//add longitudinal segment to boundary
								err |= addElementToBndry(partition, lowerLongSegment, GEODESIC);

								//reset longitudinal segment switch
								lowerLongSegment = 0;

								//free resources since the addElementToBndry method generates copies in memory
								free(tempLowerLongSegStartPt);
								tempLowerLongSegStartPt = NULL;
								free(lowerLongSegment);
								lowerLongSegment = NULL;

							}

							if (!upperLongSegmentFound)
							{
								//save the first intersection point for the longitudinal segment that intersects at the current shape
								tempUpperLongSegStartPt = (LLPoint*) malloc(sizeof(LLPoint));
								err |= createPt(tempUpperLongSegStartPt, upperLongitudinalIntx.latitude, upperLongitudinalIntx.longitude);

								upperLongSegmentFound = 1;
							} else {

								//create instance of longitudinal segment
								upperLongSegment = (Geodesic*) malloc(sizeof(Geodesic));
								err |= createGeo(upperLongSegment, *tempUpperLongSegStartPt, upperLongitudinalIntx, SEGMENT, eps);

								//add longitudinal segment to boundary
								err |= addElementToBndry(partition, upperLongSegment, GEODESIC);

								//reset longitudinal segment switch
								upperLongSegment = 0;

								//free resources since the addElementToBndry method generates copies in memory
								free(tempUpperLongSegStartPt);
								tempUpperLongSegStartPt = NULL;
								free(upperLongSegment);
								upperLongSegment = NULL;

							}

						}

					}
					else
					{
						//something really bad happened, this code should never be reached
						err |= UNEXPECTED_ERR;
						return err;
					}

					break;

				case LOCUS:
		            segmentLocus = (Locus*) segment;
					startPt = segmentLocus->locusStart;
					endPt = segmentLocus->locusEnd;

					//check if start/end points lie between the longitudinal extrema for this boundary partition
					startPtInside = (longitudeLowerBoundForPartition <= startPt.longitude && startPt.longitude <= longitudeUpperBoundForPartition);
					endPtInside = (longitudeLowerBoundForPartition <= endPt.longitude && endPt.longitude <= longitudeUpperBoundForPartition);

					//find infinite longitudinal intersections
					lowerLongError = locusGeoIntx(lowerInfLongGeo.startPoint, lowerInfLongGeo.endPoint, *segmentLocus, &lowerLongitudinalIntx, tol, eps);
					upperLongError = locusGeoIntx(upperInfLongGeo.startPoint, upperInfLongGeo.endPoint, *segmentLocus, &upperLongitudinalIntx, tol, eps);

					//add masked intersection errors to master error log for this function
					err |= getMaskedError(lowerLongError, getMask(0, 0, 1, 0, 0, 0, 0));
					err |= getMaskedError(upperLongError, getMask(0, 0, 1, 0, 0, 0, 0));

					//process the shape
					if (startPtInside && endPtInside)
					{
						//add shape to partition, no truncation needed
						err |= addElementToBndry(partition, segment, orderedBoundary.elements[j].type);

					}
					else if (startPtInside && !endPtInside)
					{
						//shape should intersect only one infinite longitudinal geodesic

						//process shape
						if ( !(lowerLongError & NO_INTERSECTION_ERR) )//using bitwise logical AND operation
						{
							//truncate shape at intersection
							ptIsOnLocus(*segmentLocus, lowerLongitudinalIntx, &tempPt, &err, tol, eps);
							tempDist1 = distToLocusFromGeoPt(*segmentLocus, tempPt, NULL, &err, tol, eps);
							err |= createLocus(&tempLocus, segmentLocus->geoStart, tempPt, segmentLocus->startDist, tempDist1, SEGMENT, tol, eps);

							//add shape to boundary
							err |= addElementToBndry(partition, &tempLocus, orderedBoundary.elements[j].type);

							//process longitudinal segment
							if (!lowerLongSegmentFound)
							{
								//save the first intersection point for the longitudinal segment that intersects at the current shape
								tempLowerLongSegStartPt = (LLPoint*) malloc(sizeof(LLPoint));
								err |= createPt(tempLowerLongSegStartPt, lowerLongitudinalIntx.latitude, lowerLongitudinalIntx.longitude);

								lowerLongSegmentFound = 1;
							} else {

								//create instance of longitudinal segment
								lowerLongSegment = (Geodesic*) malloc(sizeof(Geodesic));
								err |= createGeo(lowerLongSegment, *tempLowerLongSegStartPt, lowerLongitudinalIntx, SEGMENT, eps);

								//add longitudinal segment to boundary
								err |= addElementToBndry(partition, lowerLongSegment, GEODESIC);

								//reset longitudinal segment switch
								lowerLongSegment = 0;

								//free resources since the addElementToBndry method generates copies in memory
								free(tempLowerLongSegStartPt);
								tempLowerLongSegStartPt = NULL;
								free(lowerLongSegment);
								lowerLongSegment = NULL;

							}

						}
						else if ( !(upperLongError & NO_INTERSECTION_ERR) )//using bitwise logical AND operation
						{
							//truncate shape at intersection
							ptIsOnLocus(*segmentLocus, upperLongitudinalIntx, &tempPt, &err, tol, eps);
							tempDist1 = distToLocusFromGeoPt(*segmentLocus, tempPt, NULL, &err, tol, eps);
							err |= createLocus(&tempLocus, segmentLocus->geoStart, tempPt, segmentLocus->startDist, tempDist1, SEGMENT, tol, eps);

							//add shape to boundary
							err |= addElementToBndry(partition, &tempLocus, orderedBoundary.elements[j].type);

							//process longitudinal segment
							if (!upperLongSegmentFound)
							{
								//save the first intersection point for the longitudinal segment that intersects at the current shape
								tempUpperLongSegStartPt = (LLPoint*) malloc(sizeof(LLPoint));
								err |= createPt(tempUpperLongSegStartPt, upperLongitudinalIntx.latitude, upperLongitudinalIntx.longitude);

								upperLongSegmentFound = 1;
							} else {

								//create instance of longitudinal segment
								upperLongSegment = (Geodesic*) malloc(sizeof(Geodesic));
								err |= createGeo(upperLongSegment, *tempUpperLongSegStartPt, upperLongitudinalIntx, SEGMENT, eps);

								//add longitudinal segment to boundary
								err |= addElementToBndry(partition, upperLongSegment, GEODESIC);

								//reset longitudinal segment switch
								upperLongSegment = 0;

								//free resources since the addElementToBndry method generates copies in memory
								free(tempUpperLongSegStartPt);
								tempUpperLongSegStartPt = NULL;
								free(upperLongSegment);
								upperLongSegment = NULL;

							}
						}
						else
						{
							//intersection exepected but not found
							err |= UNEXPECTED_ERR;
							return err;
						}

					}
					else if (!startPtInside && endPtInside)
					{

						//shape should intersect only one infinite longitudinal geodesic

						//process shape
						if ( !(lowerLongError & NO_INTERSECTION_ERR) )//using bitwise logical AND operation
						{
							//truncate shape at intersection
							ptIsOnLocus(*segmentLocus, lowerLongitudinalIntx, &tempPt, &err, tol, eps);
							tempDist1 = distToLocusFromGeoPt(*segmentLocus, tempPt, NULL, &err, tol, eps);
							err |= createLocus(&tempLocus, tempPt, segmentLocus->geoEnd, tempDist1, segmentLocus->endDist, SEGMENT, tol, eps);

							//add shape to boundary
							err |= addElementToBndry(partition, &tempLocus, orderedBoundary.elements[j].type);

							//process longitudinal segment
							if (!lowerLongSegmentFound)
							{
								//save the first intersection point for the longitudinal segment that intersects at the current shape
								tempLowerLongSegStartPt = (LLPoint*) malloc(sizeof(LLPoint));
								err |= createPt(tempLowerLongSegStartPt, lowerLongitudinalIntx.latitude, lowerLongitudinalIntx.longitude);

								lowerLongSegmentFound = 1;
							} else {

								//create instance of longitudinal segment
								lowerLongSegment = (Geodesic*) malloc(sizeof(Geodesic));
								err |= createGeo(lowerLongSegment, *tempLowerLongSegStartPt, lowerLongitudinalIntx, SEGMENT, eps);

								//add longitudinal segment to boundary
								err |= addElementToBndry(partition, lowerLongSegment, GEODESIC);

								//reset longitudinal segment switch
								lowerLongSegment = 0;

								//free resources since the addElementToBndry method generates copies in memory
								free(tempLowerLongSegStartPt);
								tempLowerLongSegStartPt = NULL;
								free(lowerLongSegment);
								lowerLongSegment = NULL;

							}

						}
						else if ( !(upperLongError & NO_INTERSECTION_ERR) )//using bitwise logical AND operation
						{
							//truncate shape at intersection
							ptIsOnLocus(*segmentLocus, upperLongitudinalIntx, &tempPt, &err, tol, eps);
							tempDist1 = distToLocusFromGeoPt(*segmentLocus, tempPt, NULL, &err, tol, eps);
							err |= createLocus(&tempLocus, tempPt, segmentLocus->geoEnd, tempDist1, segmentLocus->endDist, SEGMENT, tol, eps);

							//add shape to boundary
							err |= addElementToBndry(partition, &tempLocus, orderedBoundary.elements[j].type);

							//process longitudinal segment
							if (!upperLongSegmentFound)
							{
								//save the first intersection point for the longitudinal segment that intersects at the current shape
								tempUpperLongSegStartPt = (LLPoint*) malloc(sizeof(LLPoint));
								err |= createPt(tempUpperLongSegStartPt, upperLongitudinalIntx.latitude, upperLongitudinalIntx.longitude);

								upperLongSegmentFound = 1;
							} else {

								//create instance of longitudinal segment
								upperLongSegment = (Geodesic*) malloc(sizeof(Geodesic));
								err |= createGeo(upperLongSegment, *tempUpperLongSegStartPt, upperLongitudinalIntx, SEGMENT, eps);

								//add longitudinal segment to boundary
								err |= addElementToBndry(partition, upperLongSegment, GEODESIC);

								//reset longitudinal segment switch
								upperLongSegment = 0;

								//free resources since the addElementToBndry method generates copies in memory
								free(tempUpperLongSegStartPt);
								tempUpperLongSegStartPt = NULL;
								free(upperLongSegment);
								upperLongSegment = NULL;

							}
						}
						else
						{
							//intersection exepected but not found
							err |= UNEXPECTED_ERR;
							return err;
						}

					}
					else if (!startPtInside && !endPtInside)
					{
						//check if shape straddles the longitudinal boundary, i.e. two intersections expected
						if ( !(lowerLongError & NO_INTERSECTION_ERR) && !(upperLongError & NO_INTERSECTION_ERR) )
						{

							//TODO handle the two cases: 1. shape.startPt closer to lowerLongitudinal bounding geo; shape.startPt closer to upperLongitudinal bounding geo
							//truncate shape at intersection
							ptIsOnLocus(*segmentLocus, lowerLongitudinalIntx, &tempPt, &err, tol, eps);
							tempDist1 = distToLocusFromGeoPt(*segmentLocus, tempPt, NULL, &err, tol, eps);
							ptIsOnLocus(*segmentLocus, upperLongitudinalIntx, &tempPt2, &err, tol, eps);
							tempDist2 = distToLocusFromGeoPt(*segmentLocus, tempPt2, NULL, &err, tol, eps);
							err |= createLocus(&tempLocus, tempPt, tempPt2, tempDist1, tempDist2, SEGMENT, tol, eps);

							//add shape to boundary
							err |= addElementToBndry(partition, &tempLocus, orderedBoundary.elements[j].type);

							//process longitudinal segments
							if (!lowerLongSegmentFound)
							{
								//save the first intersection point for the longitudinal segment that intersects at the current shape
								tempLowerLongSegStartPt = (LLPoint*) malloc(sizeof(LLPoint));
								err |= createPt(tempLowerLongSegStartPt, lowerLongitudinalIntx.latitude, lowerLongitudinalIntx.longitude);

								lowerLongSegmentFound = 1;
							} else {

								//create instance of longitudinal segment
								lowerLongSegment = (Geodesic*) malloc(sizeof(Geodesic));
								err |= createGeo(lowerLongSegment, *tempLowerLongSegStartPt, lowerLongitudinalIntx, SEGMENT, eps);

								//add longitudinal segment to boundary
								err |= addElementToBndry(partition, lowerLongSegment, GEODESIC);

								//reset longitudinal segment switch
								lowerLongSegment = 0;

								//free resources since the addElementToBndry method generates copies in memory
								free(tempLowerLongSegStartPt);
								tempLowerLongSegStartPt = NULL;
								free(lowerLongSegment);
								lowerLongSegment = NULL;

							}

							if (!upperLongSegmentFound)
							{
								//save the first intersection point for the longitudinal segment that intersects at the current shape
								tempUpperLongSegStartPt = (LLPoint*) malloc(sizeof(LLPoint));
								err |= createPt(tempUpperLongSegStartPt, upperLongitudinalIntx.latitude, upperLongitudinalIntx.longitude);

								upperLongSegmentFound = 1;
							} else {

								//create instance of longitudinal segment
								upperLongSegment = (Geodesic*) malloc(sizeof(Geodesic));
								err |= createGeo(upperLongSegment, *tempUpperLongSegStartPt, upperLongitudinalIntx, SEGMENT, eps);

								//add longitudinal segment to boundary
								err |= addElementToBndry(partition, upperLongSegment, GEODESIC);

								//reset longitudinal segment switch
								upperLongSegment = 0;

								//free resources since the addElementToBndry method generates copies in memory
								free(tempUpperLongSegStartPt);
								tempUpperLongSegStartPt = NULL;
								free(upperLongSegment);
								upperLongSegment = NULL;

							}

						}

					}
					else
					{
						//something really bad happened, this code should never be reached
						err |= UNEXPECTED_ERR;
						return err;
					}

					break;

				case ARC:
//		            segmentArc = (Arc*) segment;
//					startPt = segmentArc->startPoint;
//					endPt = segmentArc->endPoint;
					//TODO incomplete - finish this code
					break;

				case LLPOINT:
					//TODO do something
					break;

				case SPIRAL:
					//TODO do something
					break;
				default:
					err |= SHAPE_NOT_DEFINED_ERR;
			}
		}

		//add boundary partition to complex boundary
		err |= addElementToComplexBoundary(complexBoundary, partition);
	}

	return err;
}


//TODO check whether this algorithm holds for cases along the prime meridian where modular arithmetic problems may occur (i.e. -180 degrees = 180 degrees, or close to these regions)
ErrorSet complexBoundaryCircleIntersectionExists(Boundary boundary, Arc circles[], int circlesSize, int checkSurface,
                                    Boundary* intersectionList, Boundary* noIntersectionList, double tol, double eps)
{
	ErrorSet err = 0;
	Boundary orderedBoundary;
	const double RAD2DEG = 180.0 / M_PI;
	int longitudeSetSize;//the dimension of the array that holds the set of longitudes from the shapes that comprise the ordered boundary
	double* orderedBoundaryLongitudes = NULL;//the array that holds the set of longitudes from the shapes that comprise the ordered boundary
	double longitudeMax;//the maximum longitude from the array that holds the set of longitudes from the shapes that comprise the ordered boundary
	double longitudeMin;//the minimum longitude from the array that holds the set of longitudes from the shapes that comprise the ordered boundary
	int longitudeLowerBound;//the integer lower bound from the array that holds the set of longitudes from the shapes that comprise the ordered boundary
	int longitudeUpperBound;//the integer upper bound from the array that holds the set of longitudes from the shapes that comprise the ordered boundary
	int i = 0;
	int j = 0;
    Shape* segment = NULL; //The shape currently being processed
    Arc* segmentArc = NULL; //Temp object to hold segment recast as an Arc
    Geodesic* segmentGeo = NULL; //Temp object to hold segment recast as an Geodesic
    Locus* segmentLocus = NULL; //Temp object to hold segment recast as an Locus
    LLPoint startPt, endPt;
    ComplexBoundary partitionedBoundarySet;
    Arc tempCircle;
    LLPoint tempPt;
    int lessThanLongMin = 0;
    int greaterThanLongMax = 0;
    int circleBoundaryIntersection = 0;
    int circleLongMin;
    int circleLongMax;
    int evaluationPartitionIndex;
    int partitionSpanOfCircle;


	orderedBoundary = createBndry();

	//order boundary so shapes are sequentially connected by start/end points
	err |= orderBndry(boundary, &orderedBoundary, tol, eps);

	//find longitudinal max/min
    longitudeSetSize = 2 * orderedBoundary.length;
	orderedBoundaryLongitudes = (double*) realloc(orderedBoundaryLongitudes, longitudeSetSize * sizeof(double));
	for(i = 0; i < orderedBoundary.length; i++)
	{
		segment = orderedBoundary.elements[i].this_shape;
		switch(orderedBoundary.elements[i].type)
		{
			case GEODESIC:
	            segmentGeo = (Geodesic*) segment;
				startPt = segmentGeo->startPoint;
				endPt = segmentGeo->endPoint;
				break;

			case LOCUS:
	            segmentLocus = (Locus*) segment;
				startPt = segmentLocus->locusStart;
				endPt = segmentLocus->locusEnd;
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
				//TODO do something
				break;
			default:
				err |= SHAPE_NOT_DEFINED_ERR;
		}
		orderedBoundaryLongitudes[2*i] = startPt.longitude;
		orderedBoundaryLongitudes[2*i + 1] = endPt.longitude;
	}
	err |= findSetMaxAndMin(orderedBoundaryLongitudes, longitudeSetSize, &longitudeMax, &longitudeMin);
	if (err)
	{
		free(orderedBoundaryLongitudes);
		orderedBoundaryLongitudes = NULL;
		return err;
	}
	//find the integer longitudinal bounds
	longitudeLowerBound = (int)floor(longitudeMin*RAD2DEG);//using degrees here
	longitudeUpperBound = (int)ceil(longitudeMax*RAD2DEG);//using degrees here

	//partition ordered boundary longitudinally
	partitionedBoundarySet = newComplexBoundary();
	err |= longitudinallyPartitionBoundary(orderedBoundary, &partitionedBoundarySet, tol, eps);

	//start processing circles using the partitions from an ordered boundary
	for (i = 0; i < circlesSize; i++)
	{
		//get the circle being evaluated
		err |= createArc(&tempCircle, circles[i].centerPoint, circles[i].startPoint, circles[i].endPoint, circles[i].dir, tol, eps);

		//exclude circles that are outside longitudinal extrema
		//check against longitude min/max instead of lower/upper bounds to get a "tighter" evaluation
		lessThanLongMin = tempCircle.centerPoint.longitude < longitudeMin;
		greaterThanLongMax = tempCircle.centerPoint.longitude > longitudeMax;

		if (lessThanLongMin)
		{
			//check circle point closest to the longitudinal min of the ordered boundary
			err |= direct(tempCircle.centerPoint, M_PI_2, tempCircle.radius, &tempPt, eps);
			if (tempPt.longitude < longitudeMin)
			{
				//circle is completely outside the ordered boundary
				circleBoundaryIntersection = 0;
			}
			else
			{
				//circle may intersect ordered boundary, evaluate circle against the least partition
				evaluationPartitionIndex = 0;
				err |= bndryCircleIntxExists(*(partitionedBoundarySet.elements[evaluationPartitionIndex]), tempCircle, checkSurface, &circleBoundaryIntersection, tol, eps);
			}
		}
		else if (greaterThanLongMax)
		{
			//check circle point closest to the longitudinal max of the ordered boundary
			err |= direct(tempCircle.centerPoint, M_PI + M_PI_2, tempCircle.radius, &tempPt, eps);
			if (tempPt.longitude > longitudeMax)
			{
				//circle is completely outside the ordered boundary
				circleBoundaryIntersection = 0;
			}
			else
			{
				//circle may intersect ordered boundary, evaluate circle against the greatest partition
				evaluationPartitionIndex = partitionedBoundarySet.length - 1;
				err |= bndryCircleIntxExists(*(partitionedBoundarySet.elements[evaluationPartitionIndex]), tempCircle, checkSurface, &circleBoundaryIntersection, tol, eps);
			}
		}
		else
		{
			//determine the span of the circle across the partition set
			err |= direct(tempCircle.centerPoint, M_PI_2, tempCircle.radius, &tempPt, eps);
			circleLongMin = (int) floor(tempPt.longitude * RAD2DEG);//using degrees here
			err |= direct(tempCircle.centerPoint, M_PI + M_PI_2, tempCircle.radius, &tempPt, eps);
			circleLongMax = (int) ceil(tempPt.longitude * RAD2DEG);//using degrees here
			partitionSpanOfCircle = circleLongMax - circleLongMin;
			//maximum function handles cases where circle long min is less than the lower bound
			//this may happen for circles that whose center is close to the lower long bound but have large radii that spans multiple meridians
			evaluationPartitionIndex = (int) maximum(circleLongMin - longitudeLowerBound, 0);

			//handle special case if span >= size of partition set
			if (partitionSpanOfCircle >= partitionedBoundarySet.length)
			{
				//circle spans across entire longitudinal domain thus across the entire partition set
				//in this case just evaluate the circle against the original ordered boundary w/o using partitions
				err |= bndryCircleIntxExists(orderedBoundary, tempCircle, checkSurface, &circleBoundaryIntersection, tol, eps);
			}
			else
			{
				circleBoundaryIntersection = 0;
				//start evaluating circle with respect to longitudinal partition
				//note, the circle may need to be evaluated against multiple partitions depending on whether the ordered boundary is concave and the radius of the circle
				//stop evaluating against the partition set after the first intersection is encountered
				while(!circleBoundaryIntersection && j < partitionSpanOfCircle)
				{
					err |= bndryCircleIntxExists(*(partitionedBoundarySet.elements[evaluationPartitionIndex]), tempCircle, checkSurface, &circleBoundaryIntersection, tol, eps);
					evaluationPartitionIndex++;
					j++;
				}
			}
		}

		//place the evaluated circle in the appropriate list
		if (circleBoundaryIntersection)
		{
			//circle intersects ordered boundary
			err |= addElementToBndry(intersectionList, &tempCircle, ARC);
		}
		else
		{
			//circle does not intersect ordered boundary
			err |= addElementToBndry(noIntersectionList, &tempCircle, ARC);
		}


	}

	return err;
}
