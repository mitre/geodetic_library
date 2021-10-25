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
//#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#if REPLACE_WITH_AMDLIBM
#include "amdlibm.h"
#endif
#include "Geolib.h"

namespace geolib_idealab {
/* This equivalence is also made in libWGS84.c.  Should figure out a better way
 * to make it so that it's not repeated.
 */

/*******************************************************************************
 * Allocate and initialize a new Spiral structure from given input parameters
 */

int ptIsInSet(LLPoint pt, LLPointSet pts, double tol) {
	int i;
	LLPoint temp;

	for (i=0; i < pts.length; i++) {
		temp.latitude = pts.elements[i]->latitude;
		temp.longitude = pts.elements[i]->longitude;
		if (ptsAreSame(pt, temp, 100*tol)) {
			return 1;
		}
	}
	return 0;
}

double principalAngle( double angle )
{
  /* Determine/return the principal value (-PI < pangle <= PI) of an angle
     Input and return value in radians */

  double pangle = angle;

  while (pangle > M_PI)
    pangle -=  M_2PI;
  while (pangle < -M_PI)
    pangle +=  M_2PI;

  return( pangle);
}

//TODO:  libWGS84 has computeSubtendedAngle, could replace this function.
double calculateSubtendedAngle(double startAngle, double endAngle, int dir)
{
	double angle;

	if (dir == COUNTERCLOCKWISE) {
		angle = endAngle;
		endAngle = startAngle;
		startAngle = angle;
	}

	if (startAngle > endAngle) {
		startAngle = startAngle - M_2PI;
	}

	angle = endAngle - startAngle;

	return angle;
}

ErrorSet createSpiral(Spiral* sp, LLPoint center, double startRad, double endRad,
                    double startAz, double endAz, ArcDirection dir, double eps)
{

    LLPoint endPoint, startPoint;

    ErrorSet err = 0;

    if (sp == NULL)
    {
        err |= NO_MEMORY_ALLOCATED_ERR;
        return err;
    }

    sp->centerPoint = center;

    err = direct(center, startAz, startRad, &startPoint, eps);
    err |= direct(center, endAz, endRad, &endPoint, eps);

    sp->startPoint = startPoint;
    sp->endPoint = endPoint;

    sp->startRadius = startRad;
    sp->endRadius = endRad;

    startAz = fmod(startAz, 2*M_PI);
    endAz = fmod(endAz, 2*M_PI);
    if (fabs(endAz - startAz) <= 5e-6) {
    	endAz = startAz - dir * 5e-6;
    }
    if (startAz < 0) {
    	startAz = startAz + 2*M_PI;
    }
    if (endAz < 0) {
        endAz = endAz + 2*M_PI;
    }

    sp->startAz = startAz;
    sp->endAz = endAz;

    sp->dir = dir;
    sp->subtendedAngle = calculateSubtendedAngle(startAz, endAz, dir);

    sp->growthRate = (endRad - startRad) / sp->subtendedAngle;

    return err;

}

double azDifference(double az1, double az2)

{
	double difference;

	difference = az2 - az1;

	if (fabs(difference) - M_PI < .00001)

	{
		return difference;
	}

	az1 = principalAngle(az1);
	az2 = principalAngle(az2);

	difference = az2 - az1;

	return difference;

}

ErrorSet spiralRadius(Spiral sp, double az, double* radius)

{
	ErrorSet err = 0;
	double subAngle;

	if (radius == NULL)
	{
		return NO_MEMORY_ALLOCATED_ERR;
	}

	if (az == sp.startAz)
	{
		*radius = sp.startRadius;
		return err;
	}
	if (az == sp.endAz)
	{
		*radius = sp.endRadius;
		return err;
	}

	az = fmod(az, 2*M_PI);

	if (az < 0)
	{
		az = 2*M_PI + az;
	}

	subAngle = calculateSubtendedAngle(sp.startAz, az, sp.dir);
	*radius = sp.startRadius + subAngle * sp.growthRate;

	return err;

}

int azIsInSpiralExtent(Spiral sp, double az)

{
	double angle;
	int val = 0;

	angle = calculateSubtendedAngle(sp.startAz, az, sp.dir);
	if (fabs(sp.subtendedAngle) >= fabs(angle))
	{
		val = 1;
	}
	return val;
}

ErrorSet ptOnSpiral(Spiral sp, double az, LLPoint* pt, double eps)

{
	ErrorSet err = 0;
	double rad;

	err |= spiralRadius(sp, az, &rad);
    err |= direct(sp.centerPoint, az, rad, pt, eps);
    return err;
}

int ptIsOnSpiral(Spiral sp, LLPoint point, double tol, double eps)

{
	ErrorSet err = 0;
	int val = 0;
	double az12, az21, dist, rad;
	if (err |= inverse(sp.centerPoint, point, &az12, &az21,
		     &dist, eps)) {
		return 0;
	}

	if (ptsAreSame(point, sp.startPoint, tol)) {
		val = 1;
		return val;
	}

	if (ptsAreSame(point, sp.endPoint, tol)) {
			val = 1;
			return val;
	}

	if (azIsInSpiralExtent(sp, az12)) {
		err |= spiralRadius(sp, az12, &rad);
		if (fabs(dist - rad) < tol)
		{
	        val = 1;
		}
	}
	return val;
}

ErrorSet moveSpiralEndToPt(Spiral sp, LLPoint pt, Spiral* newSp, double tol, double eps)

{
	ErrorSet err = 0;
	double az12, az21, dist;
	LLPoint spiralPt;
	Spiral tempSp;

	if (newSp == NULL) {
		return NO_MEMORY_ALLOCATED_ERR;
	}

	err |= inverse(sp.centerPoint, pt, &az12, &az21,
		     &dist, eps);

	err |= ptOnSpiral(sp, az12, &spiralPt, eps);
	err |= invDist(pt, spiralPt, &dist, eps);

	if (dist > tol)
	{
		return POINT_NOT_ON_ARC_ERR ;
	}

	err |= spiralRadius(sp, az12, &dist);
	err |= createSpiral(&tempSp, sp.centerPoint, sp.startRadius, dist, sp.startAz, az12, sp.dir, eps);

	*newSp = tempSp;

	return err;

}

ErrorSet moveSpiralStartToPt(Spiral sp, LLPoint pt, Spiral* newSp, double tol, double eps)
{
	ErrorSet err = 0;
	double az12, az21, dist;
	LLPoint spiralPt;
	Spiral tempSp;

	if (newSp == NULL) {
		return NO_MEMORY_ALLOCATED_ERR;
	}

	err = inverse(sp.centerPoint, pt, &az12, &az21,
		     &dist, eps);

	err |= ptOnSpiral(sp, az12, &spiralPt, eps);
	err |= invDist(pt, spiralPt, &dist, eps);

	if (dist > tol)
	{
		return POINT_NOT_ON_ARC_ERR ;
	}
	err |= spiralRadius(sp, az12, &dist);
	err |= createSpiral(&tempSp, sp.centerPoint, dist, sp.endRadius, az12, sp.endAz, sp.dir, eps);
	*newSp = tempSp;

	return err;

}

ErrorSet spiralTanCrs(Spiral sp, double az, double* tanCrs)
{
	ErrorSet err = 0;
	double az12, az21;
	double rad;
	double perpAz, tanAngle;
	LLPoint pt;
	double eps = 1e-20;

	if (tanCrs == NULL) {
		return NO_MEMORY_ALLOCATED_ERR;
	}

	if (err |= spiralRadius(sp, az, &rad)) {
		return err;
	}

	if (rad == 0) {
		tanCrs = &az;
		return err;
	}

	err |= ptOnSpiral(sp, az, &pt, eps);
	err |= invCrs(sp.centerPoint, pt, &az12, &az21, eps);
	perpAz = fmod(az21 + M_PI, 2*M_PI) + sp.dir * M_PI / 2;
	tanAngle = atan(fabs(sp.growthRate)/ rad);
	if(sp.growthRate == 0){
		*tanCrs = fmod(perpAz, 2*M_PI);
	} else {
		*tanCrs = fmod((perpAz - (fabs(sp.growthRate)/sp.growthRate) * sp.dir * tanAngle), 2*M_PI);
	}

	return err;
}

ErrorSet createSpiralSection(Spiral inputsp, double az, double rad, Spiral* sp, double eps) {

	double srad, erad, saz, eaz, daz, temp;
	Spiral tempSp;
	ErrorSet err = 0;

	saz = fmod(az - inputsp.dir*3*M_PI/4, 2*M_PI);
	eaz = fmod(az + inputsp.dir*3*M_PI/4, 2*M_PI);

	srad = rad - inputsp.growthRate*3*M_PI/4;
	erad = rad + inputsp.growthRate*3*M_PI/4;

	//printf("Radii:  %f %f %f\n", rad, srad, erad);

	if (srad < 0) {
		daz = fabs(srad / inputsp.growthRate);
		saz = fmod(saz + inputsp.dir*daz, 2*M_PI);
		srad = 0;
	}
	if (erad < 0) {
		daz = fabs(erad / inputsp.growthRate);
		eaz = fmod(eaz - inputsp.dir*daz, 2*M_PI);
		erad = 0;
		if (inputsp.dir == ArcDirection::CLOCKWISE) {
		   inputsp.dir = ArcDirection::COUNTERCLOCKWISE;
		}
		else {
         inputsp.dir = ArcDirection::CLOCKWISE;
		}
		temp = srad;
		srad = erad;
		erad = temp;
		temp = saz;
		saz = eaz;
		eaz = temp;
	}

	if ((srad == 0) && (erad == 0)) {
		err |= RADIUS_OUT_OF_RANGE_ERR;
		//printf("Error!!!!\n");
		return err;
	}

	//printf("Azs: %f %f\n", saz, eaz);

	err |= createSpiral(&tempSp, inputsp.centerPoint, srad, erad, saz, eaz, inputsp.dir, eps);

	*sp = tempSp;

	return err;
}

ErrorSet spiralMidChord(Spiral sp, Geodesic line, LLPoint* midChord, double tol, double eps) {

	double azDiff, testAz, lineAz;
		double spCrs, dist13, az23, dist23, az12, az21, dist, crs2sp, crs2line, lineCrs, revLineCrs, lineDiff, spDiff;

		LLPoint spPt, linePt, perpPt;
		ErrorSet err = 0;
		int index = 0;

		double distarray[2] = { 0.0, 0.0 };
		double errArray[2] = { 0.0, 0.0 };

		azDiff = 100;

		err |= projectToGeo(line.startPoint, line.startAz, sp.centerPoint, &perpPt, &testAz, &dist, tol, eps);

		if (dist < tol) {
			err |= UNEXPECTED_ERR;
			return err;
		}

		err |= inverse(line.startPoint, perpPt, &az12, &az21, &dist, eps);

		testAz = testAz + sp.dir * sp.growthRate * .0138;

		if (fabs(azDifference(az21, line.startAz) < M_PI / 2)) {
			lineAz  = az21;
		} else {
			lineAz = fmod(az21 + M_PI / 2, 2*M_PI);
		}

		err |= crsIntx(sp.centerPoint, testAz, &az12, &dist13, line.startPoint, line.startAz, &az23, &dist23, &linePt, tol, eps);

		err |= inverse(perpPt, linePt, &lineCrs, &revLineCrs, &dist, eps);
		err |= invCrs(sp.centerPoint, linePt, &testAz, &crs2line, eps);

		err |= ptOnSpiral(sp, testAz, &spPt, eps);

		err |= spiralTanCrs(sp, testAz, &spCrs);
		err |= invCrs(sp.centerPoint, spPt, &az12, &crs2sp, eps);

		if (fabs(azDifference(revLineCrs, spCrs)) > M_PI / 2) {
			revLineCrs = fmod(revLineCrs + M_PI, 2*M_PI);
		}

		spDiff = azDifference(spCrs, crs2sp);
		lineDiff = azDifference(revLineCrs, crs2line);

		distarray[0] = dist;
		errArray[0] = spDiff - lineDiff;

		dist = dist * 1.01;
		err |= direct(perpPt, lineCrs, dist, &linePt, eps);
		err |= inverse(perpPt, linePt, &lineCrs, &revLineCrs, &dist, eps);
		err |= invCrs(sp.centerPoint, linePt, &testAz, &crs2line, eps);

		err |= ptOnSpiral(sp, testAz, &spPt, eps);
		err |= spiralTanCrs(sp, testAz, &spCrs);
		err |= invCrs(sp.centerPoint, spPt, &az12, &crs2sp, eps);

		if (fabs(azDifference(revLineCrs, spCrs)) > M_PI / 2) {
			revLineCrs = fmod(revLineCrs + M_PI, 2*M_PI);
		}

		spDiff = azDifference(spCrs, crs2sp);
		lineDiff = azDifference(revLineCrs, crs2line);

		distarray[1] = dist;
		errArray[1] = spDiff - lineDiff;

		index = 0;

		while (1) {

			if (fabs(spDiff - lineDiff) < 1e-10) {
				break;
			}

			dist = findRootSecantMethod(distarray, errArray, &err);
			err |= direct(perpPt, lineCrs, dist, &linePt, eps);
			err |= inverse(perpPt, linePt, &lineCrs, &revLineCrs, &dist, eps);
			err |= invCrs(sp.centerPoint, linePt, &testAz, &crs2line, eps);

			err |= ptOnSpiral(sp, testAz, &spPt, eps);
			err |= spiralTanCrs(sp, testAz, &spCrs);
			err |= invCrs(sp.centerPoint, spPt, &az12, &crs2sp, eps);

			if (fabs(azDifference(revLineCrs, spCrs)) > M_PI / 2) {
				revLineCrs = fmod(revLineCrs + M_PI, 2*M_PI);
			}

			spDiff = azDifference(spCrs, crs2sp);
			lineDiff = azDifference(revLineCrs, crs2line);

	        distarray[0] = distarray[1];
	        errArray[0] = errArray[1];

	        distarray[1] = dist;
	        errArray[1] = spDiff - lineDiff;

	        index += 1;

			if (index > 100) {
				err |= ITERATION_MAX_REACHED_ERR;
				break;
			}
		}

		if (azIsInSpiralExtent(sp, testAz)) {
			*midChord = spPt;
		} else {
			err |= POINT_NOT_ON_ARC_ERR;
		}

		return err;
	}

ErrorSet spiralTanPts(Spiral sp, Geodesic line, LLPointSet* pts, double tol, double eps) {

Spiral tempSp;
double perpAz, rad;
LLPoint tempPt;
ErrorSet err = 0;

perpAz = fmod(line.startAz + M_PI / 2, 2*M_PI);

err |= spiralRadius(sp, perpAz, &rad);

err |= createSpiral(&tempSp, sp.centerPoint, rad - sp.growthRate*3*M_PI/4, rad + sp.growthRate*3*M_PI/4,
		fmod(perpAz - sp.dir*3*M_PI/4, 2*M_PI), fmod(perpAz + sp.dir*3*M_PI/4, 2*M_PI), sp.dir, eps);

err |= spiralMidChord(tempSp, line, &tempPt, tol, eps);
if (err) {
	err = 0;
}
err |= addPtToPtSet(pts, &tempPt);

perpAz = fmod(line.startAz - M_PI / 2, 2*M_PI);

err |= spiralRadius(sp, perpAz, &rad);

err |= createSpiral(&tempSp, sp.centerPoint, rad - sp.growthRate*3*M_PI/4, rad + sp.growthRate*3*M_PI/4,
		fmod(perpAz - sp.dir*3*M_PI/4, 2*M_PI), fmod(perpAz + sp.dir*3*M_PI/4, 2*M_PI), sp.dir, eps);

err |= spiralMidChord(tempSp, line, &tempPt, tol, eps);
if (err) {
	err = 0;
}
err |= addPtToPtSet(pts, &tempPt);

return err;

}


} //namespace
