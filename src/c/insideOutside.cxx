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

ErrorSet geoMeridianIntx(Geodesic geo, double longitude, LLPoint *intx, double tol, double eps){

	ErrorSet err = 0;
	double dL;
	double R = SPHERE_RADIUS_NMI;
	double error;
	double errorArray[2], distArray[2];
	double angle;

	int k = 0;
	double delta;
	double dist, distCorrection;

	double alpha;
	double originalLon = geo.startPoint.longitude;
	double lontol;

	// Rotate co-ordinates to make the geo.startPoint.longitude the prime meridian
	geo.startPoint.longitude = 0.0;
	geo.endPoint.longitude = modlon(geo.endPoint.longitude - originalLon);
	longitude = modlon(longitude - originalLon);


	// Test for intx existence
	lontol = tol/findN(geo.endPoint)*cos(geo.endPoint.latitude); // Extend the line by tol to catch the edge case
	alpha = modcrs(geo.startAz);
	if (0 < alpha && alpha < M_PI){
		if (longitude > (geo.endPoint.longitude + lontol) || longitude < (0.0 - lontol)){
			err = NO_INTERSECTION_ERR;
			return err;
		}
	} else {
		if (longitude < (geo.endPoint.longitude - lontol) || longitude > (0.0 + lontol)){
			err = NO_INTERSECTION_ERR;
			return err;
		}
	}

	// 1st Approximation
	dL = findN(geo.startPoint)*cos(geo.startPoint.latitude)*((longitude));
	angle = fabs(M_PI_2 - geo.startAz);
	dist = R*atan(tan(dL/R)/cos(angle));
	err |= direct(geo.startPoint, geo.startAz, dist, intx, eps);
	error = findN(*intx)*cos(intx->latitude)*((longitude - intx->longitude));

	errorArray[0] = error;
	distArray[0] = dist;

	// 2nd Approximation
	err |= invCrs(geo.startPoint, *intx, NULL, &angle, eps);
	angle = angle - M_PI;
	angle = fabs(M_PI_2 - angle);
	distCorrection = R*atan(tan(error/R)/cos(angle));
	dist = dist + distCorrection;
	err |= direct(geo.startPoint, geo.startAz, dist, intx, eps);
	error = findN(*intx)*cos(intx->latitude)*((longitude - intx->longitude));

	errorArray[1] = error;
	distArray[1] = dist;

    // Calculate Step Size
    delta = distArray[1] - distArray[0];

    /* Iteration */
    while ( ((k == 0) || (fabs(error) > tol) || (fabs(delta) > tol))
    		&& (k < MAX_ITERATIONS)
           )
    {

        // Solve for the root of the error function
        dist = findRootSecantMethod(distArray, errorArray, &err);

    	// Place intx
        err |= direct(geo.startPoint, geo.startAz, dist, intx, eps);

        // Calculate Error
        error = findN(*intx)*cos(intx->latitude)*((longitude - intx->longitude));

        // Update the error function arrays
        distArray[0] = distArray[1];
        errorArray[0] = errorArray[1];
        distArray[1] = dist;
        errorArray[1] = error;

        // Calculate Step Size
        delta = distArray[1] - distArray[0];

        k++;
	}


    // Make sure the edge case is actually on the line
	if (0 < alpha && alpha < M_PI){
		if (longitude > geo.endPoint.longitude || longitude < 0){
		    if(!(ptIsOnGeo(geo.startPoint,geo.endPoint, *intx, SEGMENT, &err, tol, eps)))
				err = NO_INTERSECTION_ERR;
		}
	} else {
		if (longitude < geo.endPoint.longitude || longitude > 0){
		    if(!(ptIsOnGeo(geo.startPoint,geo.endPoint, *intx, SEGMENT, &err, tol, eps)))
				err = NO_INTERSECTION_ERR;
		}
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

ErrorSet locusMeridianIntx(Locus locus, double longitude, LLPoint *intx, double tol, double eps){

	ErrorSet err = 0;
	double dL;
	double R = SPHERE_RADIUS_NMI;
	double error;
	double errorArray[2], distArray[2];
	double angle, crs;

	int k = 0;
	double delta;
	double geodist, locdist;
	LLPoint tempPt, geoTempPt;

	double alpha;
	double originalLon = locus.locusStart.longitude;
	double lontol;

	// Rotate co-ordinates to make the geo.startPoint.longitude the prime meridian
	locus.locusStart.longitude = 0.0;
	locus.locusEnd.longitude = modlon(locus.locusEnd.longitude - originalLon);
	locus.geoStart.longitude = modlon(locus.geoStart.longitude - originalLon);
	locus.geoEnd.longitude = modlon(locus.geoEnd.longitude - originalLon);
	longitude = modlon(longitude - originalLon);


	// Test for intx existence
	lontol = tol/findN(locus.locusEnd)*cos(locus.locusEnd.latitude); // Extend the locus by tol for the edge case
	err |= invCrs(locus.locusStart, locus.locusEnd, &alpha, NULL, eps);
	if (0 < alpha && alpha < M_PI){
		if (longitude > (locus.locusEnd.longitude + lontol) || longitude < (0.0 - lontol)){
			err = NO_INTERSECTION_ERR;
			return err;
		}
	} else {
		if (longitude < (locus.locusEnd.longitude - lontol) || longitude > (0.0 + lontol)){
			err = NO_INTERSECTION_ERR;
			return err;
		}
	}

	// 1st Approximation
	dL = findN(locus.locusStart)*cos(locus.locusStart.latitude)*((longitude));
	angle = fabs(M_PI_2 - alpha);
	locdist = R*atan(tan(dL/R)/cos(angle));
	err |= direct(locus.locusStart, alpha, locdist, &tempPt, eps);
	err |= projectToLocus(locus, tempPt, intx, NULL, NULL, tol, eps);
	err |= projectToGeo(locus.geoStart, locus.geoAz, *intx, &geoTempPt, NULL, NULL, tol, eps);
//	err |= invDist(locus.geoStart, geoTempPt, &geodist, eps);
	err |= inverse(locus.geoStart, geoTempPt, &crs, NULL, &geodist, eps);
	if(fabs(crs - locus.geoAz) > M_PI_2)
		geodist = -geodist;
	error = findN(*intx)*cos(intx->latitude)*((longitude - intx->longitude));

	errorArray[0] = error;
	distArray[0] = geodist;

	// 2nd Approximation
	err |= invCrs(locus.locusStart, *intx, NULL, &alpha, eps);
	alpha = alpha - M_PI;
	angle = fabs(M_PI_2 - alpha);
	locdist = R*atan(tan(error/R)/cos(angle));
	err |= direct(*intx, alpha, locdist, &tempPt, eps);
	err |= projectToLocus(locus, tempPt, intx, NULL, NULL, tol, eps);
	err |= projectToGeo(locus.geoStart, locus.geoAz, *intx, &geoTempPt, NULL, NULL, tol, eps);
//	err |= invDist(locus.geoStart, geoTempPt, &geodist, eps);
	err |= inverse(locus.geoStart, geoTempPt, &crs, NULL, &geodist, eps);
	if(fabs(crs - locus.geoAz) > M_PI_2)
		geodist = -geodist;
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
        err |= direct(locus.geoStart, locus.geoAz, geodist, &tempPt, eps);
        err |= ptOnLocusFromGeoPt(locus, tempPt, intx, NULL, tol, eps);

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


    // Make sure the edge case is really on the locus
	if (0 < alpha && alpha < M_PI){
		if (longitude > locus.locusEnd.longitude || longitude < 0.0){
			if(!ptIsOnLocus(locus, *intx, NULL, &err, tol, eps))
				err = NO_INTERSECTION_ERR;
		}
	} else {
		if (longitude < locus.locusEnd.longitude || longitude > 0.0){
			if(!ptIsOnLocus(locus, *intx, NULL, &err, tol, eps))
				err = NO_INTERSECTION_ERR;
		}
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


ErrorSet arcMeridianIntx(Arc arc, double longitude, LLPointPair intx, int *numberOfIntx, double tol, double eps){

	ErrorSet err = 0;
	LLPoint meridianStart = {0.0, longitude};
	double meridianCrs = 0.0;
	double crsToMeridian;
	double distToMeridian;
	double a, da;
	double crs;
	double error = 1.0;
	double errorArray[2];
	double crsArray[2];
	double R = SPHERE_RADIUS_NMI;

	double delta;
	int k = 0;
	int i = 0;

	int intx0OnArc, intx1OnArc;

	// Test for intx existence
	err |= projectToGeo(meridianStart, meridianCrs, arc.centerPoint, NULL, &crsToMeridian, &distToMeridian, tol, eps);
	if(distToMeridian - arc.radius > tol/2){
		// No Intx
		*numberOfIntx = 0;
		return NO_INTERSECTION_ERR;
	}
	if(fabs(distToMeridian - arc.radius) < tol/2){
		// Tangent
		*numberOfIntx = 1;
	}
	if(distToMeridian - arc.radius < -tol/2){
		// 2 Intx
		*numberOfIntx = 2;
	}


	for(i = 0; i < *numberOfIntx; i++){
		// 1st Approximation
		a = acos(tan(distToMeridian/R)/tan(arc.radius/R));
		crs = i == 0 ? crsToMeridian - a : crsToMeridian + a;
		err |= direct(arc.centerPoint, crs, arc.radius, &intx[i], eps);

		// Provision for when the arc encompasses a pole
		if(fabs(intx[i].longitude - longitude) > M_PI_2)
			longitude = modlon(longitude + M_PI);

		error = findN(intx[i])*cos(intx[i].latitude)*((longitude - intx[i].longitude));

		errorArray[0] = error;
		crsArray[0] = crs;

		// 2nd Approximation
//		da = acos((cos(error/R) - pow(cos(arc.radius),2))/pow(sin(arc.radius),2));
		da = error/(R*sin(arc.radius/R));
//		crs = i == 0 ? crs + sgn(error)*da : crs - sgn(error)*da;
		crs = i == 0 ? crs + da : crs - da;
		err |= direct(arc.centerPoint, crs, arc.radius, &intx[i], eps);
		error = findN(intx[i])*cos(intx[i].latitude)*((longitude - intx[i].longitude));

		errorArray[1] = error;
		crsArray[1] = crs;

		// Calculate Step Size
		delta = (crsArray[1] - crsArray[0])*arc.radius;

		/* Iteration */
		while ( ((k == 0) || (fabs(error) > tol) || (fabs(delta) > tol))
				&& (k < MAX_ITERATIONS)
			   )
		{

			// Solve for the root of the error function
			crs = findRootSecantMethod(crsArray, errorArray, &err);

			// Place intx
			err |= direct(arc.centerPoint, crs, arc.radius, &intx[i], eps);

			// Calculate Error
			error = findN(intx[0])*cos(intx[i].latitude)*((longitude - intx[i].longitude));

			// Update the error function arrays
			crsArray[0] = crsArray[1];
			errorArray[0] = errorArray[1];
			crsArray[1] = crs;
			errorArray[1] = error;

			// Calculate Step Size
			delta = (crsArray[1] - crsArray[0])*arc.radius;

			k++;
		}
	}

	// Check to see if intx are on arc
	if(*numberOfIntx == 2){
		intx0OnArc = ptIsOnArc(arc.centerPoint, arc.radius, arc.startAz, arc.endAz, arc.dir, intx[0], &err, tol, eps);
		intx1OnArc = ptIsOnArc(arc.centerPoint, arc.radius, arc.startAz, arc.endAz, arc.dir, intx[1], &err, tol, eps);
		if(!intx0OnArc && !intx1OnArc){
			*numberOfIntx = 0;
		} else if (intx0OnArc && !intx1OnArc){
			*numberOfIntx = 1;
		} else if (!intx0OnArc && intx1OnArc){
			*numberOfIntx = 1;
			intx[0] = intx[1];
		}
	} else if(*numberOfIntx == 1) {
		intx0OnArc = ptIsOnArc(arc.centerPoint, arc.radius, arc.startAz, arc.endAz, arc.dir, intx[0], &err, tol, eps);
		if(!intx0OnArc){
			*numberOfIntx = 0;
		}
	}

	if (k >= MAX_ITERATIONS)
	{
		err |= ITERATION_MAX_REACHED_ERR;
	}

	if (fabs(error) >= MAX_DISTANCE_ERROR)
	{
		err |= ERROR_MAX_REACHED_ERR;
	}

	if(*numberOfIntx == 0)
		err |= NO_INTERSECTION_ERR;

	return err;
}

ErrorSet spiralMeridianIntx(Spiral spiral, double longitude, LLPointPair intx, int *numberOfIntx, double tol, double eps){

	ErrorSet err = 0;
	LLPoint meridianStart = {0.0, longitude};
	double meridianCrs = 0.0;
	double crsToMeridian;
	double distToMeridian;
	double a, da;
	double crs;
	double error = 1.0;
	double errorArray[2];
	double crsArray[2];
	double R = SPHERE_RADIUS_NMI;

	double delta;
	int k = 0;
	int i = 0;

	int intx0OnArc, intx1OnArc;

	double radiusTowardsMeridian, radius;

	// Test for intx existence
	err |= projectToGeo(meridianStart, meridianCrs, spiral.centerPoint, NULL, &crsToMeridian, &distToMeridian, tol, eps);
	spiralRadius(spiral, crsToMeridian, &radiusTowardsMeridian);
	if(distToMeridian - radiusTowardsMeridian > tol/2){
		// No Intx
		*numberOfIntx = 0;
		return err;
	}
	if(fabs(distToMeridian - radiusTowardsMeridian) < tol/2){
		// Tangent
		*numberOfIntx = 1;
	}
	if(distToMeridian - radiusTowardsMeridian < -tol/2){
		// 2 Intx
		*numberOfIntx = 2;
	}


	for(i = 0; i < *numberOfIntx; i++){
		// 1st Approximation
		a = acos(tan(distToMeridian/R)/tan(radiusTowardsMeridian/R));
		crs = i == 0 ? crsToMeridian - a : crsToMeridian + a;
		spiralRadius(spiral, crs, &radius);
		err |= direct(spiral.centerPoint, crs, radius, &intx[i], eps);

		// Provision for when the spiral encompasses a pole
		if(fabs(intx[i].longitude - longitude) > M_PI_2)
			longitude = modlon(longitude + M_PI);

		error = findN(intx[i])*cos(intx[i].latitude)*((longitude - intx[i].longitude));

		errorArray[0] = error;
		crsArray[0] = crs;

		// 2nd Approximation
//		da = acos((cos(error/R) - pow(cos(spiral.radius),2))/pow(sin(spiral.radius),2));
		da = error/(R*sin(radius/R));
//		crs = i == 0 ? crs + sgn(error)*da : crs - sgn(error)*da;
		crs = i == 0 ? crs + da : crs - da;
		spiralRadius(spiral, crs, &radius);
		err |= direct(spiral.centerPoint, crs, radius, &intx[i], eps);
		error = findN(intx[i])*cos(intx[i].latitude)*((longitude - intx[i].longitude));

		errorArray[1] = error;
		crsArray[1] = crs;

		// Calculate Step Size
		delta = (crsArray[1] - crsArray[0])*radius;

		/* Iteration */
		while ( ((k == 0) || (fabs(error) > tol) || (fabs(delta) > tol))
				&& (k < MAX_ITERATIONS)
			   )
		{

			// Solve for the root of the error function
			crs = findRootSecantMethod(crsArray, errorArray, &err);
			spiralRadius(spiral, crs, &radius);
			// Place intx
			err |= direct(spiral.centerPoint, crs, radius, &intx[i], eps);

			// Calculate Error
			error = findN(intx[0])*cos(intx[i].latitude)*((longitude - intx[i].longitude));

			// Update the error function arrays
			crsArray[0] = crsArray[1];
			errorArray[0] = errorArray[1];
			crsArray[1] = crs;
			errorArray[1] = error;

			// Calculate Step Size
			delta = (crsArray[1] - crsArray[0])*radius;

			k++;
		}
	}

	// Check to see if intx are on spiral
	if(*numberOfIntx == 2){
		intx0OnArc = ptIsOnSpiral(spiral, intx[0], tol, eps);
		intx1OnArc = ptIsOnSpiral(spiral, intx[1], tol, eps);
		if(!intx0OnArc && !intx1OnArc){
			*numberOfIntx = 0;
		} else if (intx0OnArc && !intx1OnArc){
			*numberOfIntx = 1;
		} else if (!intx0OnArc && intx1OnArc){
			*numberOfIntx = 1;
			intx[0] = intx[1];
		}
	} else if(*numberOfIntx == 1) {
		intx0OnArc = ptIsOnSpiral(spiral, intx[0], tol, eps);
		if(!intx0OnArc){
			*numberOfIntx = 0;
		}
	}

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

// This may be useful to deal with the polar case
void removeAntipodalLonPts(LLPoint intx[], int* intxNumber, double longitude, double tol){

	int i = 0;
	int j;

	//Remove AntopodalLonPts
	while(i < *intxNumber - 1){
		if(fabs(intx[i].longitude - longitude) > M_PI_2){
			for(j = i + 1; j < *intxNumber - 1; j++){
				intx[j] = intx[j + 1];
			}
			*intxNumber = *intxNumber - 1;
		} else {
			i++;
		}
	}

}


void orderMeridianPts(LLPoint intx[], int* intxNumber, double tol){

	LLPoint tempPt;
	int i, j, k;

	switch(*intxNumber){
		case 0:{
			break;
		}
		case 1:{
			break;
		}
		default:{
			// Selection Sort
			for(i = 0; i < *intxNumber; i++){
				k = i;
				for(j = i + 1; j < *intxNumber; j++){
					if( intx[k].latitude > intx[j].latitude){
						k = j;
					}
				}
				if(k != i){
					tempPt = intx[i];
					intx[i] = intx[k];
					intx[k] = tempPt;
				}
			}

			break;
		}
	}

	//Remove Duplicates
	i = 0;
	while(i < *intxNumber - 1){
		if(ptsAreSame(intx[i], intx[i + 1], tol)){
			for(j = i + 1; j < *intxNumber - 1; j++){
				intx[j] = intx[j + 1];
			}
			*intxNumber = *intxNumber - 1;
		} else {
			i++;
		}
	}

}

ErrorSet lineBoundarySet(Boundary B, double longitude, LLPoint intxSet[], int* numberOfIntx, int elementsIntx[], int* numElementsIntx,  double tol, double eps){

	ErrorSet err = 0;
	int i = 0, k = 0, n = 0;

	Geodesic* geo;
	Arc* arc;
	Locus* locus;
	Spiral* spiral;
	Locus2ndOrder* locus2ndOrder;

	Geodesic meridian;
	LLPoint geoStartPoint = {0, longitude};
	LLPoint geoEnd = {M_PI_2, longitude};
	createGeo(&meridian, geoStartPoint, geoEnd, INFINITE, eps);
	LLPoint intx;
	int intxCount = 0;
	int elementCount = 0;
	LLPointSet tempIntxSet = createPtSet();
	LLPoint tempIntx[2];
	Shape* tempShape;

	for(i = 0; i < B.length; i++){
			tempShape = B.elements[i].this_shape;
			switch(B.elements[i].type){
				case GEODESIC:
					geo = (Geodesic*) tempShape;
					err = geoMeridianIntx(*geo, longitude, &intx, tol, eps); // Applies segment bounds
					if(err == SUCCESS){
						intxSet[intxCount] = intx;
						intxCount++;
						elementsIntx[elementCount] = i;
						elementCount++;
					} else if(err == NO_INTERSECTION_ERR){
						err = 0;
					}
					break;

				case LOCUS:
					locus = (Locus*) tempShape;
					err = locusMeridianIntx(*locus, longitude, &intx, tol, eps); // Applies segment bounds
					if(err == SUCCESS){
						intxSet[intxCount] = intx;
						intxCount++;
						elementsIntx[elementCount] = i;
						elementCount++;
					} else if(err == NO_INTERSECTION_ERR){
						err = 0;
					}
					break;

				case ARC:
					arc = (Arc*) tempShape;
					err = arcMeridianIntx(*arc, longitude, tempIntx, &n, tol, eps); // Applies segment bounds
					if(err == SUCCESS){
						for(k = 0; k < n; k++){
							intxSet[intxCount] = tempIntx[k];
							intxCount++;
						}
						elementsIntx[elementCount] = i;
						elementCount++;
					} else if(err == NO_INTERSECTION_ERR){
						err = 0;
					}
					break;

				case SPIRAL:
					spiral = (Spiral*) tempShape;
					err = spiralGeoIntx(*spiral,
										meridian,
										&tempIntxSet, tol, eps);  // Applies segment bounds
					if(err == SUCCESS){
						for(k = 0; k < tempIntxSet.length; k++){
							if(ptIsOnSpiral(*spiral, *tempIntxSet.elements[k], tol, eps)){
									intxSet[intxCount] = *tempIntxSet.elements[k];
									intxCount++;
							}
							elementsIntx[elementCount] = i;
							elementCount++;
						}
					} else if(err == NO_INTERSECTION_ERR){
						err = 0;
					}
					clearPtSet(&tempIntxSet);
					tempIntxSet = createPtSet();
					break;

				case LOCUS2NDORDER:
					locus2ndOrder = (Locus2ndOrder*) tempShape;
					err = locus2ndOrderMeridianIntx(*locus2ndOrder, longitude, &intx, tol, eps); // Applies segment bounds
					if(err == SUCCESS){
						intxSet[intxCount] = intx;
						intxCount++;
						elementsIntx[elementCount] = i;
						elementCount++;
					} else if(err == NO_INTERSECTION_ERR){
						err = 0;
					}
					break;

				default:
					break;
			}

	}
	orderMeridianPts(intxSet, &intxCount, tol);
//	removeAntipodalLonPts(intxSet, &intxCount, longitude, tol);
	*numberOfIntx = intxCount;
	*numElementsIntx = elementCount;
	return err;
}

int checkEdgeCases(Boundary b, LLPoint point, LLPoint intxSet[], int numberOfIntx, int elementsIntx[], int numElementsIntx, double tol, double eps){

	ErrorSet err = 0;
	int inside = 0;
	int i;
	Shape* tempShape;

	Geodesic* geo;
	Arc* arc;
	Locus* locus;
	Spiral* spiral;


	for(i = 0; i < numElementsIntx; i++){
		if(inside == 1)
			break;

		tempShape = b.elements[elementsIntx[i]].this_shape;
		switch(b.elements[elementsIntx[i]].type){
			case GEODESIC:
				geo = (Geodesic*) tempShape;
				inside = ptIsOnGeo(geo->startPoint, geo->endPoint, point, geo->lineType, &err, tol, eps);
				break;
			case LOCUS:
				locus = (Locus*) tempShape;
				inside = ptIsOnLocus(*locus, point, NULL, &err, tol, eps);
				break;
			case ARC:
				arc = (Arc*) tempShape;
				inside = ptIsOnArc(arc->centerPoint, arc->radius, arc->startAz, arc->endAz, arc->dir, point, &err, tol, eps);
				break;
			case SPIRAL:
				spiral = (Spiral*) tempShape;
				inside = ptIsOnSpiral(*spiral, point, tol, eps);
				break;
			default:
				break;
		}
	}

	return inside;
}


ErrorSet ptsAreInsideBndry(Boundary b, double latList[], double lonList[], int inside[], int numberOfPts, double tol, double eps){

	ErrorSet err = 0;
	int i, j;
	double longitude = lonList[0];
	int numberOfIntx;
	LLPoint intxSet[10];
	LLPoint pt;
	int elementsIntx[b.length];
	int numElementsIntx;
	double lattol = tol/SPHERE_RADIUS_NMI;

	// Find intersections for the first meridian
	err |= lineBoundarySet(b, lonList[0], intxSet, &numberOfIntx, elementsIntx, &numElementsIntx, tol, eps);

	i = 0;
	while(i < numberOfPts){

		// If point is on the same meridian check for inside/outside
		if(lonList[i] == longitude){

			// Check for edge cases
			inside[i] = 0;
			for(j = 0; j < numberOfIntx; j++){
				if(fabs(latList[i] - intxSet[j].latitude) < 2000*lattol){
					err = createPt(&pt, latList[i], lonList[i]);
					inside[i] = checkEdgeCases(b, pt, intxSet, numberOfIntx, elementsIntx, numElementsIntx, tol, eps);
					break;
				} else {
					inside[i] = 0;
				}
			}

			// Check for inside/outside
			if(inside[i] == 0){
				switch(numberOfIntx){
					case 0:
						inside[i] = 0;
						break;
					case 1:
						inside[i] = 0;
						break;
					default:
						if(latList[i] < intxSet[0].latitude){
							inside[i] = 0;
							break;
						}
						for(j = 1; j < numberOfIntx; j++){
							if((latList[i] > intxSet[j - 1].latitude) && (latList[i] < intxSet[j].latitude)){
								inside[i] = j % 2; // Between even and odd inside
								break;
							}
						}
						if(latList[i] > intxSet[numberOfIntx - 1].latitude){
							inside[i] = 0;
							break;
						}
						break;
					}
			}

			i++;

		} else {

			// If point is not on the same meridian find intersections for the new meridian
			err |= lineBoundarySet(b, lonList[i], intxSet, &numberOfIntx, elementsIntx, &numElementsIntx, tol, eps);
			longitude = lonList[i];
		}
	}

	return err;
}

//void sortPoints(double latList[], double lonList[], int idx[], int len) {
//	int numLons = 0, i = 1, pos = 0, lonCount = 0;
//	double tempLat[len], tempLon[len];
//	int tempIdx[len];
//	double lat;
//	lat = latList[0];
//	while (numLons == 0) {
//		if (latList[i] != lat){
//			numLons = i;
//		}
//		i++;
//	}
//	printf("Number of lons:  %i\n", numLons);
//
//	for (i=0;i<len;i++) {
//		tempLat[i] = latList[pos+lonCount];
//		tempLon[i] = lonList[pos+lonCount];
//		tempIdx[i] = idx[pos+lonCount];
//
//		pos = pos + numLons;
//		if (pos >= len) {
//			pos = 0;
//			lonCount++;
//		}
//	}
//
//	latList = tempLat;
//	lonList = tempLon;
//	idx = tempIdx;
//	return;
//}

///* The following algorithm uses a bubble sort method to sort the array.
// * If this proves to be too slow a better sorting algorithm could
// * be implemented
// */
//
void sortPoints(double latList[], double lonList[], int idx[], int len) {
	int flag = 1, i = 0;
	double temp;

	for (i=0;i<len;i++) {
		idx[i] = i;
	}
	while (flag) {
		flag = 0;
		for (i=0;i<len-1;i++) {
			if (lonList[i+1] > lonList[i]) {
				temp = lonList[i];
				lonList[i] = lonList[i+1];
				lonList[i+1] = temp;

				temp = latList[i];
				latList[i] = latList[i+1];
				latList[i+1] = temp;

				temp = idx[i];
				idx[i] = idx[i+1];
				idx[i+1] = temp;
				flag = 1;
			}
		}
	}
	return;
}

void sortArrayByIndex(int array[], int idx[], int len) {
	int temp[len];
	int i;

	for (i=0;i<len;i++) {
		temp[idx[i]] = array[i];
	}

	array = temp;
}

