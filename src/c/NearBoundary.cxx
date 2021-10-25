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

ErrorSet ptOnLocus2FromGeoDist(Locus loc, double secondDist, double geoDist, LLPoint* loc2Pt, double tol, double eps){

	ErrorSet err = 0;
	LLPoint locPt;
	double locCrs;
	double crsFromLocPtToLoc2Pt;

	locCrs = locusCrsAtGeoDist(loc, geoDist, &locPt, NULL, &err, tol, eps);
	crsFromLocPtToLoc2Pt = locCrs + M_PI_2;
	err |= direct(locPt, crsFromLocPtToLoc2Pt, secondDist, loc2Pt, eps);

	return err;
}

ErrorSet locus2ndOrderMeridianIntx(Locus2ndOrder locus2ndOrder, double longitude, LLPoint *intx, double tol, double eps){

	Locus locus = locus2ndOrder.locus;
	double secondDist = locus2ndOrder.secondDist;

	ErrorSet err = 0;
	double dL;
	double R = SPHERE_RADIUS_NMI;
	double error;
	double errorArray[2], distArray[2];
	double angle, crs;

	int k = 0;
	double delta;
	double geodist, locdist;
	LLPoint tempPt, geoTempPt, locTempPt;

	double alpha;
	double originalLon;
	double lonTolStart, lonTolEnd;

	LLPoint loc2Start, loc2End;

	// Find the end points of loc2
	ptOnLocus2FromGeoDist(locus, secondDist, 0, &loc2Start, tol, eps);
	ptOnLocus2FromGeoDist(locus, secondDist, locus.geoLength, &loc2End, tol, eps);
	originalLon = loc2Start.longitude;

	// Rotate co-ordinates to make the loc2Start.longitude the prime meridian
	loc2Start.longitude = 0.0;
	loc2End.longitude = modlon(loc2End.longitude - originalLon);
	locus.locusStart.longitude = modlon(locus.locusStart.longitude - originalLon);
	locus.locusEnd.longitude = modlon(locus.locusEnd.longitude - originalLon);
	locus.geoStart.longitude = modlon(locus.geoStart.longitude - originalLon);
	locus.geoEnd.longitude = modlon(locus.geoEnd.longitude - originalLon);
	longitude = modlon(longitude - originalLon);


	// Test for intx existence
	lonTolStart = tol/findN(loc2Start)*cos(loc2Start.latitude); // Extend the locus by tol for the edge case
	lonTolEnd = tol/findN(loc2End)*cos(loc2End.latitude); // Extend the locus by tol for the edge case
	err |= invCrs(loc2Start, loc2End, &alpha, NULL, eps);
	if (0 < alpha && alpha < M_PI){
		if (longitude > (loc2End.longitude + lonTolEnd) || longitude < (0.0 - lonTolStart)){
			err = NO_INTERSECTION_ERR;
			return err;
		}
	} else {
		if (longitude < (loc2End.longitude - lonTolEnd) || longitude > (0.0 + lonTolStart)){
			err = NO_INTERSECTION_ERR;
			return err;
		}
	}

	// 1st Approximation
	dL = findN(loc2Start)*cos(loc2Start.latitude)*((longitude));
	angle = fabs(M_PI_2 - alpha);
	locdist = R*atan(tan(dL/R)/cos(angle));
	err |= direct(locus.locusStart, alpha, locdist, &tempPt, eps);
	err |= projectToLocus(locus, tempPt, &locTempPt, NULL, NULL, tol, eps);
	err |= projectToGeo(locus.geoStart, locus.geoAz, locTempPt, &geoTempPt, NULL, NULL, tol, eps);
	err |= inverse(locus.geoStart, geoTempPt, &crs, NULL, &geodist, eps);
	if(fabs(crs - locus.geoAz) > M_PI_2)
		geodist = -geodist;
	err |= ptOnLocus2FromGeoDist(locus, secondDist, geodist, intx, tol, eps);
	error = findN(*intx)*cos(intx->latitude)*((longitude - intx->longitude));

	errorArray[0] = error;
	distArray[0] = geodist;

	// 2nd Approximation
	err |= invCrs(loc2Start, *intx, NULL, &alpha, eps);
	alpha = alpha - M_PI;
	angle = fabs(M_PI_2 - alpha);
	locdist = R*atan(tan(error/R)/cos(angle));
	err |= direct(*intx, alpha, locdist, &tempPt, eps);
	err |= projectToLocus(locus, tempPt, &locTempPt, NULL, NULL, tol, eps);
	err |= projectToGeo(locus.geoStart, locus.geoAz, locTempPt, &geoTempPt, NULL, NULL, tol, eps);
	err |= inverse(locus.geoStart, geoTempPt, &crs, NULL, &geodist, eps);
	if(fabs(crs - locus.geoAz) > M_PI_2)
		geodist = -geodist;
	err |= ptOnLocus2FromGeoDist(locus, secondDist, geodist, intx, tol, eps);
	error = findN(*intx)*cos(intx->latitude)*((longitude - intx->longitude));

	errorArray[1] = error;
	distArray[1] = geodist;

    // Calculate Step Size
    delta = distArray[1] - distArray[0];

    /* Iteration */
    while ( ((k == 0) || (fabs(error) > tol) || (fabs(delta) > tol))
    		&& (k < MAX_ITERATIONS)
           )
    {

        // Solve for the root of the error function
        geodist = findRootSecantMethod(distArray, errorArray, &err);

    	// Place intx
    	err |= ptOnLocus2FromGeoDist(locus, secondDist, geodist, intx, tol, eps);
//        err |= direct(locus.geoStart, locus.geoAz, geodist, &tempPt, eps);
//        err |= ptOnLocusFromGeoPt(locus, tempPt, intx, NULL, tol, eps);

        // Calculate Error
        error = findN(*intx)*cos(intx->latitude)*((longitude - intx->longitude));

        // Update the error function arrays
        distArray[0] = distArray[1];
        errorArray[0] = errorArray[1];
        distArray[1] = geodist;
        errorArray[1] = error;

        // Calculate Step Size
        delta = distArray[1] - distArray[0];

        k++;
	}


    // Check the edge case
    //todo improve
    if (geodist < -tol || geodist > locus.geoLength + tol){
    	err = NO_INTERSECTION_ERR;
    }

    // Undo co-ordinate change
    intx->longitude = modlon(longitude + originalLon);

    if (k >= MAX_ITERATIONS)
    {
        err |= ITERATION_MAX_REACHED_ERR;
    }

    if (fabs(error) >= MAX_DISTANCE_ERROR)
    {
        err |= ERROR_MAX_REACHED_ERR;
    }

	return err;
}

ErrorSet areTerrainPtsNearPt(LLPoint pt, double nearDist, double* latList, double* lonList, int* near, int numberOfPts, double tol, double eps){

	ErrorSet err = 0;
	Boundary b = createBndry();
	Arc arc;
	LLPoint startEndPt;

	err |= direct(pt, 0.0, nearDist, &startEndPt, eps);
	err |= createArc(&arc, pt, startEndPt, startEndPt, CLOCKWISE, tol, eps);

	err |= addArcToBndry(&b, &arc);

	err |= ptsAreInsideBndry(b, latList, lonList, near, numberOfPts, tol, eps);
	clearBndry(&b);

	return err;
}

ErrorSet areTerrainPtsNearGeoAbeam(Geodesic geo, double nearDist, double* latList, double* lonList, int* near, int numberOfPts, double tol, double eps){

	ErrorSet err = 0;
	Boundary b = createBndry();
	Locus loc1, loc3;
	Geodesic geo2, geo4;

	err |= createLocus(&loc1, geo.startPoint, geo.endPoint, nearDist, nearDist, SEGMENT, tol, eps);
	err |= createLocus(&loc3, geo.endPoint, geo.startPoint, nearDist, nearDist, SEGMENT, tol, eps);

	err |= createGeo(&geo2, loc1.locusEnd, loc3.locusStart, SEGMENT, eps);
	err |= createGeo(&geo4, loc3.locusEnd, loc1.locusStart, SEGMENT, eps);

	err |= addLocusToBndry(&b, &loc1);
	err |= addGeoToBndry(&b, &geo2);
	err |= addLocusToBndry(&b, &loc3);
	err |= addGeoToBndry(&b, &geo4);

	err |= ptsAreInsideBndry(b, latList, lonList, near, numberOfPts, tol, eps);
	clearBndry(&b);

	return err;
}

ErrorSet areTerrainPtsNearLocusAbeam(Locus locus, double nearDist, double* latList, double* lonList, int* near, int numberOfPts, double tol, double eps){

	ErrorSet err = 0;
	Boundary b = createBndry();
	Locus2ndOrder loc2nd1, loc2nd3;
	LLPoint loc2nd1Start, loc2nd1End;
	LLPoint loc2nd3Start, loc2nd3End;
	Geodesic geo2, geo4;

	loc2nd1.locus = locus;
	loc2nd1.secondDist  = nearDist;
	loc2nd3.locus = locus;
	loc2nd3.secondDist  = -nearDist;

	err |= ptOnLocus2FromGeoDist(locus, nearDist, 0.0, &loc2nd1Start, tol, eps);
	err |= ptOnLocus2FromGeoDist(locus, nearDist, locus.geoLength, &loc2nd1End, tol, eps);
	err |= ptOnLocus2FromGeoDist(locus, -nearDist, locus.geoLength, &loc2nd3Start, tol, eps);
	err |= ptOnLocus2FromGeoDist(locus, -nearDist, 0.0, &loc2nd3End, tol, eps);

	err |= createGeo(&geo2, loc2nd1End, loc2nd3Start, SEGMENT, eps);
	err |= createGeo(&geo4, loc2nd3End, loc2nd1Start, SEGMENT, eps);

	err |= addLocus2ndOrderToBndry(&b, &loc2nd1);
	err |= addGeoToBndry(&b, &geo2);
	err |= addLocus2ndOrderToBndry(&b, &loc2nd3);
	err |= addGeoToBndry(&b, &geo4);

	err |= ptsAreInsideBndry(b, latList, lonList, near, numberOfPts, tol, eps);
	clearBndry(&b);

	return err;
}

ErrorSet areTerrainPtsNearArcAbeam(Arc arc, double nearDist, double* latList, double* lonList, int* near, int numberOfPts, double tol, double eps){

	ErrorSet err = 0;
	Boundary b = createBndry();
	Arc arc1, arc3;
	LLPoint arc1Start, arc1End;
	LLPoint arc3Start, arc3End;
	Geodesic geo2, geo4;

	// Create outside arc
	err |= direct(arc.centerPoint, arc.endAz, arc.radius + nearDist, &arc3Start, eps);
	err |= direct(arc.centerPoint, arc.startAz, arc.radius + nearDist, &arc3End, eps);

	err |= createArc(&arc3, arc.centerPoint, arc3Start, arc3End, -arc.dir, tol, eps);

	if (arc.radius - nearDist <= tol) {
		arc1Start = arc.centerPoint;
		arc1End = arc.centerPoint;
	} else {
		err |= direct(arc.centerPoint, arc.startAz, arc.radius - nearDist, &arc1Start, eps);
		err |= direct(arc.centerPoint, arc.endAz, arc.radius - nearDist, &arc1End, eps);

		err |= createArc(&arc1, arc.centerPoint, arc1Start, arc1End, arc.dir, tol, eps);

		err |= addArcToBndry(&b, &arc1);
	}

	err |= createGeo(&geo2, arc1End, arc3.startPoint, SEGMENT, eps);
	err |= createGeo(&geo4, arc3.endPoint, arc1Start, SEGMENT, eps);
	
	err |= addGeoToBndry(&b, &geo2);
	err |= addArcToBndry(&b, &arc3);
	err |= addGeoToBndry(&b, &geo4);

	err |= ptsAreInsideBndry(b, latList, lonList, near, numberOfPts, tol, eps);
	clearBndry(&b);

	return err;
}

ErrorSet areTerrainPtsNearSpiralAbeam(Spiral spiral, double nearDist, double* latList, double* lonList, int* near, int numberOfPts, double tol, double eps){

	ErrorSet err = 0;
	Boundary b = createBndry();
	Spiral spiral1, spiral3;
	Geodesic geo2, geo4;

	err |= createSpiral(&spiral1, spiral.centerPoint, spiral.startRadius - nearDist, spiral.endRadius - nearDist, spiral.startAz, spiral.endAz, spiral.dir, eps);
	err |= createSpiral(&spiral3, spiral.centerPoint, spiral.endRadius + nearDist, spiral.startRadius + nearDist, spiral.endAz, spiral.startAz, -spiral.dir, eps);

	err |= createGeo(&geo2, spiral1.endPoint, spiral3.startPoint, SEGMENT, eps);
	err |= createGeo(&geo4, spiral3.endPoint, spiral1.startPoint, SEGMENT, eps);

	err |= addSpiralToBndry(&b, &spiral1);
	err |= addGeoToBndry(&b, &geo2);
	err |= addSpiralToBndry(&b, &spiral3);
	err |= addGeoToBndry(&b, &geo4);

	err |= ptsAreInsideBndry(b, latList, lonList, near, numberOfPts, tol, eps);
	clearBndry(&b);

	return err;
}

ErrorSet areTerrainPtsNearBoundary(Boundary b, double nearDist, double* latList, double * lonList, int* near, int* inside, int numberOfPts, double tol, double eps){

	ErrorSet err = 0;
	int* nearElement;
	nearElement = malloc(numberOfPts*sizeof(int));

	Shape* tempShape;
	Geodesic* geo;
	Arc* arc;
	Locus* locus;
	Spiral* spiral;

	// Find which points are inside
	err |= ptsAreInsideBndry(b, latList, lonList, inside, numberOfPts, tol, eps);

	// Find which points are near each element (only tests abeam)
	int i, j;
	for(i = 0; i < b.length; i++){

		tempShape = b.elements[i].this_shape;
		switch(b.elements[i].type){
		case GEODESIC:
			geo = (Geodesic*) tempShape;
			err = areTerrainPtsNearGeoAbeam(*geo, nearDist, latList, lonList, nearElement, numberOfPts, tol, eps);
			break;

		case LOCUS:
			locus = (Locus*) tempShape;
			err = areTerrainPtsNearLocusAbeam(*locus, nearDist, latList, lonList, nearElement, numberOfPts, tol, eps);
			break;

		case ARC:
			arc = (Arc*) tempShape;
			err = areTerrainPtsNearArcAbeam(*arc, nearDist, latList, lonList, nearElement, numberOfPts, tol, eps);
			break;

		case SPIRAL:
			spiral = (Spiral*) tempShape;
			err = areTerrainPtsNearSpiralAbeam(*spiral, nearDist, latList, lonList, nearElement, numberOfPts, tol, eps);
			break;

		default:
			break;
		}

		for(j = 0; j < numberOfPts; j++){
			near[j] = nearElement[j] || near[j];
			nearElement[j] = 0;
		}
	}

	// Find which points are near vertices
	Boundary orderedBndry = createBndry();
	orderBndry(b, &orderedBndry, tol, eps);
	for(i = 0; i < orderedBndry.length; i++){

		tempShape = orderedBndry.elements[i].this_shape;
		switch(orderedBndry.elements[i].type){
		case GEODESIC:
			geo = (Geodesic*) tempShape;
			err = areTerrainPtsNearPt(geo->endPoint, nearDist, latList, lonList, nearElement, numberOfPts, tol, eps);
			break;

		case LOCUS:
			locus = (Locus*) tempShape;
			err = areTerrainPtsNearPt(locus->locusEnd, nearDist, latList, lonList, nearElement, numberOfPts, tol, eps);
			break;

		case ARC:
			arc = (Arc*) tempShape;
			err = areTerrainPtsNearPt(arc->endPoint, nearDist, latList, lonList, nearElement, numberOfPts, tol, eps);
			break;

		case SPIRAL:
			spiral = (Spiral*) tempShape;
			err = areTerrainPtsNearPt(spiral->endPoint, nearDist, latList, lonList, nearElement, numberOfPts, tol, eps);
			break;

		default:
			break;
		}

		for(j = 0; j < numberOfPts; j++){
			near[j] = nearElement[j] || near[j];
			nearElement[j] = 0;
		}
	}

	free(nearElement);

	return err;
}














