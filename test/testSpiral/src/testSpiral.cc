/*
 * Copyright 2021 The MITRE Corporation.  All Rights reserved.
 *
 * This is the copyright work of The MITRE Corporation, and was produced
 * for the U.S. Government under Contract Number DTFA01-93-C-00001,
 * and is subject to the Federal Acquisition Regulation Clause 52.227-14.
 * Rights in Data--General, Alt. III(JUN 1987) and Alt. IV (JUN 1987).
 * No use other than that granted to the U.S. Government, or to those
 * acting on behalf of the U.S. Government, under that Clause is authorized
 * without the express written permission of The MITRE Corporation.
 * For further information, please contact The MITRE Corporation, Contracts
 * Office, 7515 Colshire Drive, McLean VA 22102, (703) 883-6000.
 */


#include <sys/types.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#if REPLACE_WITH_AMDLIBM
#include "amdlibm.h"
#endif
#include <stdarg.h>
#include <float.h>
#include <limits.h>
#include "Geolib.h"
#include "testGeolib.h"

namespace geolib_idealab {

#ifndef EPS
#define EPS 0.5e-13
#endif

#ifdef TOL
#undef EPS
#endif
#define TOL 1.37e-9

#ifndef DEG2RAD
#define DEG2RAD M_PI / 180
#endif
//double DEG2RAD = M_PI / 180;

double PrinAngle( double angle )
{
  // Determine/return the principal value (-PI < pangle <= PI) of an angle
  //   Input and return value in radians

  double pangle = angle;

  while (pangle > M_PI)
    pangle -=  M_2PI;
  while (pangle < -M_PI)
    pangle +=  M_2PI;

  return( pangle);
}

double AzDiff(double az1, double az2)

{
	double difference;

	//az1 = fmod(az1 + 2*M_PI, 2*M_PI);
	//az2 = fmod(az2 + 2*M_PI, 2*M_PI);

	difference = az2 - az1;

	if (fabs(difference) < M_PI)

	{
		return difference;
	}

	az1 = PrinAngle(az1);
	az2 = PrinAngle(az2);

	difference = az2 - az1;

	return difference;

}

TestSet testSpiralMidChord_Set1(){

		TestSet set = newTestSet("testSpiralMidChord_Set1");
		LLPoint geoStartPoint, geoEndPoint, centerPoint, midChord, testPt, spPt, perpPt;
		double latxExp, lonxExp, azDiff;
		double startRad, endRad, startAz, endAz, testAz, spCrs, perpCrs;
		ArcDirection direction;
		int val;
		double eps = 1.0e-20;
		double tol = 1.37e-9;
		long err = 0;
		Geodesic geo;
		Spiral spiral;
		double az12, az21, dist;
		double mcAz, mcDist;

		int testNum = 0, successCount = 0, errCount = 0, failCount = 0, setupCount = 0, unverifiedCount = 0;
		int i;

		srand(04012011);

		double radius[10] = {0.1,0.2,0.3,0.5,1.0,2.0,5.0,10.0,15.0,20.0};
        //double radius[10] = {0.01,0.1,0.2,0.5,1.0,5.0,10.0,20.0,35.0,50.0};

		for (i=0;i<100;i++) {

	        latxExp = randLat();
	        lonxExp = randLon();

			err |= createPt(&centerPoint, latxExp*DEG2RAD, lonxExp*DEG2RAD);
			if ((rand() % 100) > 49) {
				direction = ArcDirection::CLOCKWISE;
			} else {
				direction = ArcDirection::COUNTERCLOCKWISE;
			}
			startRad = 0.01 * randDist(); //0-54 nm
			endRad = 0.01 * randDist(); //0-54 nm
			startAz = randAzimuth();
			endAz = fmod(startAz + direction * 270 + 360, 360);
			testAz = fmod(startAz + direction * 135 + 360, 360);
			//printf("Azs:  %f %f %f %i\n", startAz, endAz, testAz, direction);

			err |= createSpiral(&spiral, centerPoint, startRad, endRad, startAz * DEG2RAD, endAz * DEG2RAD, direction, eps);

			//printf("Start/End Az:  %f %f\n", spiral.startAz, spiral.endAz);

			err |= spiralTanCrs(spiral, testAz * DEG2RAD, &spCrs);
			//printf("Tangent Course:  %f\n", spCrs);
			err |= ptOnSpiral(spiral, testAz * DEG2RAD, &spPt, eps);

			err |= inverse(spiral.centerPoint, spPt, &mcAz, &az21, &mcDist, eps);

			azDiff = AzDiff(spCrs, az21);

			testNum = 0;

			//printf("Test Az:  %f %f %f\n", mcAz, testAz * DEG2RAD, az21);

			while (testNum < 10) {

				midChord.latitude = 0;
				midChord.longitude = 0;
				//printf("Courses:  %f %f\n", spCrs, az21);
				err |= direct(centerPoint, mcAz, mcDist * radius[testNum], &geoStartPoint, eps);
				err |= invCrs(centerPoint, geoStartPoint, &az12, &az21, eps);
				err |= direct(geoStartPoint, az21 - direction * fabs(azDiff), 1, &geoEndPoint, eps);
				//printf("Az21:  %f %f\n", az21, azDiff);
				err |= createGeo(&geo, geoStartPoint, geoEndPoint, SEGMENT, eps);
				azDiff = AzDiff(az21 - direction * azDiff, az21);
				err |= inverse(centerPoint, spPt, &az12, &az21, &dist, eps);
				err |= inverse(centerPoint, geoStartPoint, &az12, &az21, &dist, eps);
				//printf("Az/Rad/LineAz:  %f %f %f\n", mcAz, mcDist, geo.startAz);
				if (err) {
					setupCount += 1;
				}
				err = 0;

				err |= spiralMidChord(spiral, geo, &midChord, tol, eps);
				//printf("Geo Start Az:  %f\n", geo.startAz);

				if (err) {
					errCount += 1;
				}

				err |= inverse(centerPoint, midChord, &az12, &az21, &dist, eps);

				//printf("Azs:  %f %f\n", az12, testAz * DEG2RAD);

				err |= ptOnSpiral(spiral, mcAz, &testPt, eps);

				val = ptsAreSame(testPt, midChord, 10*tol);
				err |= invDist(testPt, midChord, &dist, eps);
				//printf("Distance:  %e\n", dist);
				err |= inverse(centerPoint, testPt, &az12, &az21, &dist, eps);
				err |= inverse(spPt, testPt, &az12, &az21, &dist, eps);

				successCount = successCount + val;

				if (val == 0) {
					failCount += 1;
					err |= invDist(testPt, midChord, &dist, eps);
					//printf("Az/Rad/LineAz:  %f %f %f\n", mcAz, mcDist, geo.startAz);
					//printf("Direction:  %i\n", direction);
					printf("Distance between points:  %e %i %i\n", dist, i, testNum);
					printf("Geo Start Point:  %f %f\n", geoStartPoint.latitude / DEG2RAD, geoStartPoint.longitude / DEG2RAD);
					printf("Error:  %s", formatErrorMessage(err));
					printf("Spiral radii:  %f %f\n", spiral.startRadius, spiral.endRadius);
					printf("Spiral az's:  %f %f %f\n", spiral.startAz, spiral.endAz, testAz * DEG2RAD);
					//printf("Spiral center:  %f %f\n", spiral.centerPoint.latitude / DEG2RAD, spiral.centerPoint.longitude / DEG2RAD);
					err |= inverse(centerPoint, midChord, &az12, &az21, &dist, eps);
					printf("Test Az/Found Az:  %f %f\n", testAz * DEG2RAD, az12);
					if (ptIsAtPole(testPt, &err, tol, eps)) {
						printf("Point is at pole\n");
					}
					printf("Mid Chord:  %f %f\n", midChord.latitude / DEG2RAD, midChord.longitude / DEG2RAD);
					printf("Test Point:  %f %f\n\n", testPt.latitude / DEG2RAD, testPt.longitude / DEG2RAD);

					//printf("\n");
				} else {
					err |= projectToGeo(geo.startPoint, geo.startAz, centerPoint, &perpPt, &perpCrs, &dist, tol, eps);
					err |= invCrs(centerPoint, midChord, &az12, &az21, eps);
					//if (spiral.growthRate > 0) {
						//printf("Spiral Course: %f\n", spCrs);
					//	printf("Courses/Dist:  %f %f %f %f %i %i\n", perpCrs, az12, AzDiff(perpCrs, az12), spiral.growthRate, testNum, i);
					//}
				}
				//printf("\n");
				//printf("Error: %s\n", WGS84FormatErrorMessage(err));

				testNum++;
			}
		}
		set.errors = errCount;
		set.fail = failCount;
		set.pass = successCount;
		set.testCases = 1000;
		set.setupFailures = setupCount;
		set.unverified = unverifiedCount;

		return set;
}

TestSuite testSpiralMidChord_AllSets() {

	TestSuite suite;
	TestSet set1;

    printf("\nStart testSpiralMidChord_AllSets\n");

    suite = newTestSuite("testSpiralMidChord_AllSets");

    set1 = testSpiralMidChord_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("Finish testSpiralMidChord_AllSets\n\n\n");

    return suite;
}

TestSet testSpiralGeoIntx_Set1(){

		TestSet set = newTestSet("testSpiralLineIntersect_Set1");
		LLPoint geoStartPoint, geoEndPoint, centerPoint;
		double startRad, endRad, startAz, endAz, dist, totDist, spAngle, epAngle;
		ArcDirection direction;
		double eps = 1.0e-20;
		double tol = 1.37e-9;
		long err = 0;
		int i, j, spcounter = 0, epcounter = 0;
		Geodesic geo;
		Spiral spiral;

		int testNum = 0, successCount = 0, errCount = 0, failCount = 0, setupCount = 0, unverifiedCount = 0;

		double az12, az21;

		double total;

		double latxExp, lonxExp;

		testNum = 1000;

		srand(04012011);

		for (i=0;i<testNum;i++) {

			LLPointSet pointSet = createPtSet();

			latxExp = randLat();
			lonxExp = randLon();

			err |= createPt(&centerPoint, latxExp*DEG2RAD, lonxExp*DEG2RAD);
			if ((rand() % 100) > 49) {
				direction = CLOCKWISE;
			} else {
				direction = COUNTERCLOCKWISE;
			}
			startRad = 30 + 0.01 * randDist() / 2; //0-54 nm
			endRad = 30 + 0.01 * randDist() / 2; //0-54 nm
			startAz = randAzimuth();
			endAz = startAz + direction * 270;

			spAngle = rand() % 270;
			epAngle = rand() % 270;

			//printf("Azs:  %f %f %f %f\n", startAz, endAz, spAngle, epAngle);

			//printf("Current Test:  %i\n", i);

			//These values should be generated randomly
			//err |= createPt(&geoStartPoint, 25.694592395967856*DEG2RAD, -80.10724447560116*DEG2RAD);
			//err |= createPt(&geoEndPoint, 25.701834182442678*DEG2RAD, -80.09061424611792*DEG2RAD);

			err |= createSpiral(&spiral, centerPoint, startRad, endRad, startAz * DEG2RAD, endAz * DEG2RAD, direction, eps);

			//printf("Spiral Growthrate:  %f\n", spiral.growthRate);

			err |= ptOnSpiral(spiral, (startAz + direction * spAngle) * DEG2RAD, &geoStartPoint, eps);
			err |= ptOnSpiral(spiral, (startAz + direction * epAngle) * DEG2RAD, &geoEndPoint, eps);

			//printf("Azs:  %f %f\n", (startAz - direction * spAngle) * DEG2RAD, (startAz - direction * epAngle) * DEG2RAD);

			err |= createGeo(&geo, geoStartPoint, geoEndPoint, INFINITE, eps);

			err |= inverse(centerPoint, geoStartPoint, &az12, &az21, &dist, eps);
			//printf("Point 1:  %f %f\n", az12, dist);
			err |= inverse(centerPoint, geoEndPoint, &az12, &az21, &dist, eps);
			//printf("Point 2:  %f %f\n", az12, dist);

			//printf("Start Point:  %f %f\n", geo.startPoint.latitude / DEG2RAD, geo.startPoint.longitude / DEG2RAD);
			//printf("End Point:  %f %f\n", geo.endPoint.latitude / DEG2RAD, geo.endPoint.longitude / DEG2RAD);
			if (err) {
				setupCount++;
				err = 0;
			}

			err |= spiralGeoIntx(spiral, geo, &pointSet, tol, eps);

			if (err) {
				errCount++;
			}

			spcounter = 0;
			epcounter = 0;

			totDist = 0;

			double d, h, d1, d2, x, R, rad, perpCrs;
			LLPoint perpPt;
			R = 3440;

			err |= projectToGeo(geoStartPoint, geo.startAz, spiral.centerPoint, &perpPt, &perpCrs, &h, tol, eps);
			err |= invDist(spiral.centerPoint, geoStartPoint, &rad, eps);

			d2 = R * acos(cos(rad / R) / cos(h / R));
			d1 = R * acos(cos((rad - tol) / R) / cos((h + tol) / R));
			d = d2 - d1;
			x = R * acos(cos(d / R) * cos(tol / R));
			//if (x != 0) {
				//printf("Values:  %f %f %e %f %f %e %f %e\n", d1, d2, cos(d / R) * cos(tol / R) - 1.0, h, rad, d, R, x);
			//}

			//printf("Start Point:  %f %f\n", spiral.startPoint.latitude / DEG2RAD, spiral.startPoint.longitude / DEG2RAD);
			//printf("End Point:  %f %f\n", spiral.endPoint.latitude / DEG2RAD, spiral.endPoint.longitude / DEG2RAD);
			for (j=0;j<pointSet.length;j++) {
				LLPoint point;
				point.latitude = pointSet.elements[j]->latitude;
				point.longitude = pointSet.elements[j]->longitude;
				err |= invDist(point, geoStartPoint, &dist, eps);
				totDist = totDist + dist;
				err |= invDist(point, geoEndPoint, &dist, eps);
				totDist = totDist + dist;
				if (ptsAreSame(point, geoStartPoint, tol)) {
					spcounter = 1;
				}
				if (ptsAreSame(point, geoEndPoint, tol)) {
					epcounter = 1;
				}
			}

			total += x;

			//printf("Pointset Length:  %i\n", pointSet.length);

			if (pointSet.length > 2) {
				printf("Pointset:  %i %i\n", pointSet.length, i);
				spcounter = 0;
				epcounter = 0;
				for (j=0;j<pointSet.length;j++) {
					LLPoint point;
					point.latitude = pointSet.elements[j]->latitude;
					point.longitude = pointSet.elements[j]->longitude;
					err |= invDist(point, geoStartPoint, &dist, eps);
					totDist = totDist + dist;
					err |= invDist(point, geoEndPoint, &dist, eps);
					totDist = totDist + dist;
					if (ptsAreSame(point, geoStartPoint, tol)) {
						spcounter++;
					}
					if (ptsAreSame(point, geoEndPoint, tol)) {
						epcounter++;
					}
				}
				printf("Counters:  %i %i\n", spcounter, epcounter);
			}

			if (spcounter + epcounter < 2) {
				failCount++;
				printf("Pointset:  %i %i\n", pointSet.length, i);
				//printf("Azs/Dir:  %f %f %f %f %i\n", startAz, endAz, spAngle, epAngle, direction);
				//printf("Error!  %s", formatErrorMessage(err));
				//printf("Direction:  %i\n", direction);
				//printf("Failure:  %i\n", i);
				//printf("Counters:  %i %i\n", spcounter, epcounter);
				//printf("Spiral:  %f %f %f %f\n", spiral.startRadius, spiral.endRadius, spiral.growthRate, spiral.subtendedAngle);
				for (j=0;j<pointSet.length;j++) {
					LLPoint point;
					point.latitude = pointSet.elements[j]->latitude;
					point.longitude = pointSet.elements[j]->longitude;
					err |= invDist(point, geoStartPoint, &dist, eps);
					printf("Start Point:  %e\n", dist);
					err |= invDist(point, geoEndPoint, &dist, eps);
					printf("End Point:  %e\n", dist);
				}
				if (pointSet.length > 0) {
					//printf("Spiral:  %f %f %f %f\n", spiral.startRadius, spiral.endRadius, spiral.startAz, spiral.endAz);
					//printf("Start Point:  %f %f\n", spiral.startPoint.latitude / DEG2RAD, spiral.startPoint.longitude / DEG2RAD);
					//printf("End Point:  %f %f\n", spiral.endPoint.latitude / DEG2RAD, spiral.endPoint.longitude / DEG2RAD);
					//printf("Start/End:  %i %i\n", spcounter, epcounter);
					//printf("Distances:  %f %f\n", totDist, geo.length * pointSet.length);
					//printf("Point Set:  %i\n", pointSet.length);
					//printf("Failure\n\n");
				} else {
					//printf("No Points Returned\n\n");
				}
			} else {
				successCount++;
			}

			//printf("\n");

			clearPtSet(&pointSet);

			//val = ptsAreSame(testPt, midChord, tol);
			//printf("%i\n", val);

		}

		set.errors = errCount;
		set.fail = failCount;
		set.pass = successCount;
		set.testCases = testNum;
		set.setupFailures = setupCount;
		set.unverified = unverifiedCount;

		return set;

		//printf("Spiral Line Intersect:  %i/1000\n", successCount);
		//printf("Average ROA:  %e\n", total / 1000);
}

TestSuite testSpiralGeoIntx_AllSets() {

	TestSuite suite;
	TestSet set1;

    printf("\nStart testSpiralGeoIntx_AllSets\n");

    suite = newTestSuite("testSpiralGeoIntx_AllSets");

    set1 = testSpiralGeoIntx_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("Finish testSpiralGeoIntx_AllSets\n\n\n");

    return suite;
}

TestSet testSpiralLocusIntx_Set1(){

		TestSet set = newTestSet("testSpiralLocusIntersect_Set1");
		LLPoint geoStartPoint, geoEndPoint, centerPoint;
		double startRad, endRad, locCrs, startDist, endDist, startAz, endAz, dist, totDist;
		ArcDirection direction;
		double eps = 1.0e-20;
		double tol = 1.37e-9;
		long err = 0;
		int testNum = 0, successCount = 0, errCount = 0, failCount = 0, setupCount = 0, unverifiedCount = 0;
		int i, j, spcounter = 0, epcounter = 0;
		Locus loc;
		Spiral spiral;

		double az21;

		double latxExp, lonxExp;

		testNum = 1000;

		srand(04012011);

		for (i=0;i<testNum;i++) {

			err = 0;
			//printf("%i\n", i);

			LLPointSet pointSet = createPtSet();

			latxExp = randLat();
			lonxExp = randLon();

			err |= createPt(&centerPoint, latxExp*DEG2RAD, lonxExp*DEG2RAD);
			if ((rand() % 100) > 49) {
				direction = CLOCKWISE;
			} else {
				direction = COUNTERCLOCKWISE;
			}
			startRad = 30 + 0.01 * randDist() / 2; //0-54 nm
			endRad = 30 + 0.01 * randDist() / 2; //0-54 nm
			startDist = 0.01 * randDist(); //0-54 nm
			endDist = 0.01 * randDist(); //0-54 nm
			startAz = randAzimuth();
			endAz = startAz + direction * 270 ;

			//These values should be generated randomly
			//err |= createPt(&geoStartPoint, 25.694592395967856*DEG2RAD, -80.10724447560116*DEG2RAD);
			//err |= createPt(&geoEndPoint, 25.701834182442678*DEG2RAD, -80.09061424611792*DEG2RAD);

			err |= createSpiral(&spiral, centerPoint, startRad, endRad, startAz * DEG2RAD, endAz * DEG2RAD, direction, eps);

			err |= invCrs(spiral.startPoint, spiral.endPoint, &locCrs, &az21, eps);

			err |= ptOnSpiral(spiral, (startAz + direction * 60) * DEG2RAD, &geoStartPoint, eps);
			err |= ptOnSpiral(spiral, (startAz + direction * 120) * DEG2RAD, &geoEndPoint, eps);

			err |= createLocus(&loc, geoStartPoint, geoEndPoint, 0, 0, SEGMENT, tol, eps);

			//printf("Start Point:  %f %f\n", geo.startPoint.latitude / DEG2RAD, geo.startPoint.longitude / DEG2RAD);
			//printf("End Point:  %f %f\n", geo.endPoint.latitude / DEG2RAD, geo.endPoint.longitude / DEG2RAD);

			if (err) {
				printf("Setup Error:\n");
				printf("%s\n", formatErrorMessage(err));
				setupCount++;
				err = 0;
			}

			err |= spiralLocusIntx(spiral, loc, &pointSet, tol, eps);

			if (err) {
				printf("%s\n", formatErrorMessage(err));
				errCount++;
			}

			spcounter = 0;
			epcounter = 0;

			totDist = 0;

			//printf("Start Point:  %f %f\n", spiral.startPoint.latitude / DEG2RAD, spiral.startPoint.longitude / DEG2RAD);
			//printf("End Point:  %f %f\n", spiral.endPoint.latitude / DEG2RAD, spiral.endPoint.longitude / DEG2RAD);
			for (j=0;j<pointSet.length;j++) {
				LLPoint point;
				point.latitude = pointSet.elements[j]->latitude;
				point.longitude = pointSet.elements[j]->longitude;
				err |= invDist(point, geoStartPoint, &dist, eps);
				totDist = totDist + dist;
				err |= invDist(point, geoEndPoint, &dist, eps);
				totDist = totDist + dist;
				if (ptsAreSame(point, geoStartPoint, tol)) {
					spcounter = 1;
				}
				if (ptsAreSame(point, geoEndPoint, tol)) {
					epcounter = 1;
				}
			}

			if (spcounter + epcounter < 2) {
				failCount++;
				printf("Pointset Length:  %i\n", pointSet.length);
				printf("Spiral:  %f %f %f %f\n", spiral.startRadius, spiral.endRadius, spiral.growthRate, spiral.subtendedAngle);
				printf("Start Point:  %f %f\n", geoStartPoint.latitude / DEG2RAD, geoStartPoint.longitude / DEG2RAD);
				printf("End Point:  %f %f\n", geoEndPoint.latitude / DEG2RAD, geoEndPoint.longitude / DEG2RAD);
				for (j=0;j<pointSet.length;j++) {
					LLPoint point;
					point.latitude = pointSet.elements[j]->latitude;
					point.longitude = pointSet.elements[j]->longitude;
					err |= invDist(point, geoStartPoint, &dist, eps);
					printf("Start Point:  %e\n", dist);
					err |= invDist(point, geoEndPoint, &dist, eps);
					printf("End Point:  %e\n", dist);
				}
				if (pointSet.length > 0) {
					printf("Spiral:  %f %f %f %f\n", spiral.startRadius, spiral.endRadius, spiral.growthRate, spiral.subtendedAngle);
					printf("Start Point:  %f %f\n", spiral.startPoint.latitude / DEG2RAD, spiral.startPoint.longitude / DEG2RAD);
					printf("End Point:  %f %f\n", spiral.endPoint.latitude / DEG2RAD, spiral.endPoint.longitude / DEG2RAD);
					printf("Start/End:  %i %i\n", spcounter, epcounter);
					printf("Distances:  %f %i\n", totDist, pointSet.length);
					//printf("Point Set:  %i\n", pointSet.length);
					printf("Failure\n\n");
				} else {
					printf("No Points Returned\n\n");
				}
			} else {
				successCount++;
			}

			//printf("\n");

			clearPtSet(&pointSet);

			//val = ptsAreSame(testPt, midChord, tol);
			//printf("%i\n", val);

		}

		set.errors = errCount;
		set.fail = failCount;
		set.pass = successCount;
		set.testCases = testNum;
		set.setupFailures = setupCount;
		set.unverified = unverifiedCount;

		return set;

		//printf("Spiral Locus Intersect:  %i/1000\n", successCount);

}

TestSuite testSpiralLocusIntx_AllSets() {

	TestSuite suite;
	TestSet set1;

    printf("\nStart testSpiralLocusIntx_AllSets\n");

    suite = newTestSuite("testSpiralLocusIntx_AllSets");

    set1 = testSpiralLocusIntx_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("Finish testSpiralLocusIntx_AllSets\n\n\n");

    return suite;
}

TestSet testGeoTanToSpiralAtAngleToGeo_Set1(){

	TestSet set = newTestSet("testGeoTanToSpiralAtAngleToGeo_Set1");
	LLPoint geoStartPoint, tanStartPoint, endPoint, centerPoint;
	double startRad, endRad, geoLength, spCrs, startAz, endAz, testAz, bCrs, dist;
	ArcDirection direction;
	double eps = 1.0e-20;
	double tol = 1.37e-9;
	long err = 0;
	int testNum = 0, successCount = 0, errCount = 0, failCount = 0, setupCount = 0, unverifiedCount = 0;
	int i, j;
	Geodesic geo, tanGeo;
	Spiral spiral;

	double angle[5] = {10,30,45,60,85};

	double az12, az21;

	double latxExp, lonxExp;

	testNum = 200;

	srand(04012011);

	for (i=0;i<testNum;i++) {

		latxExp = randLat();
		lonxExp = randLon();

		err |= createPt(&centerPoint, latxExp*DEG2RAD, lonxExp*DEG2RAD);
		if ((rand() % 100) > 49) {
			direction = CLOCKWISE;
		} else {
			direction = COUNTERCLOCKWISE;
		}
		startRad = 20 + 0.01 * randDist(); //0-54 nm
		endRad = 20 + 0.01 * randDist(); //0-54 nm
		geoLength = randDist() * .01 + 20;
		startAz = randAzimuth();
		endAz = fmod(startAz + direction * 270, 360);

		err |= createSpiral(&spiral, centerPoint, startRad, endRad, startAz * DEG2RAD, endAz * DEG2RAD, direction, eps);

		testAz = fmod(startAz + direction * 135 + 360, 360) * DEG2RAD;

		err |= ptOnSpiral(spiral, testAz, &tanStartPoint, eps);
		err |= spiralTanCrs(spiral, testAz, &spCrs);

		err |= direct(tanStartPoint, spCrs, geoLength, &endPoint, eps);
		err |= inverse(tanStartPoint, endPoint, &spCrs, &bCrs, &dist, eps);
		//printf("Back Course:  %f\n", bCrs);

		for (j=0;j<5;j++) {
			err |= direct(endPoint, bCrs - angle[j] * DEG2RAD, geoLength, &geoStartPoint, eps);
			err |= createGeo(&geo, geoStartPoint, endPoint, INFINITE, eps);
			err |= inverse(geoStartPoint, endPoint, &az12, &az21, &dist, eps);
			//printf("From Geo Start:  %f %f\n", az12, dist);
			//printf("Line Back Course:  %f\n", az21);
			err |= inverse(tanStartPoint, endPoint, &az12, &az21, &dist, eps);
			//printf("Courses:  %f %f %f %f\n", spCrs / DEG2RAD, bCrs / DEG2RAD, geo.startAz / DEG2RAD, (bCrs - angle[j] * DEG2RAD) / DEG2RAD);
			//printf("Spiral Course:  %f %f %f\n", spCrs, az12 / DEG2RAD, az21 / DEG2RAD);
			//printf("Angle:  %f\n", angle[j] * DEG2RAD);
			if (err) {
				setupCount++;
				err = 0;
			}
			err |= geoTanToSpiralAtAngleToGeo(spiral, geo, angle[j] * DEG2RAD, &tanGeo, tol, eps);

			double rad, R = 3440.0, d, x;
			err |= invCrs(spiral.centerPoint, tanGeo.startPoint, &az12, &az21, eps);
			err |= spiralRadius(spiral, az12, &rad);
			d = R * acos(cos((rad + tol) / R) / cos((rad - tol) / R));
			x = R * acos(cos(d / R) * cos(tol / R));
			//printf("X:  %e\n", x);
			if (err) {
				errCount++;
			}
			if (ptsAreSame(tanStartPoint, tanGeo.startPoint, tol)) {
				//printf("Success!");
				successCount++;
			} else {
				failCount++;
				err |= inverse(tanStartPoint, tanGeo.startPoint, &az12, &az21, &dist, eps);
				printf("Distance/Growth/Az/lineCrs/dir:  %e %f %f %f %i\n", dist, spiral.growthRate * M_2PI, az12, spCrs, spiral.dir);
			}
			//printf("\n");
			//printf("%f %f %f %f\n", tanGeo.startPoint.latitude, tanGeo.startPoint.longitude, tanStartPoint.latitude, tanStartPoint.longitude);
		}
	}
	set.errors = errCount;
	set.fail = failCount;
	set.pass = successCount;
	set.testCases = 1000;
	set.setupFailures = setupCount;
	set.unverified = unverifiedCount;

	return set;

//	printf("Line Tangent To Spiral Angle From Line:  %i/1000\n", successCount);
}

TestSuite testGeoTanToSpiralAtAngleToGeo_AllSets() {

	TestSuite suite;
	TestSet set1;

    printf("\nStart testGeoTanToSpiralAtAngleToGeo_AllSets\n");

    suite = newTestSuite("testGeoTanToSpiralAtAngleToGeo_AllSets");

    set1 = testGeoTanToSpiralAtAngleToGeo_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("Finish testGeoTanToSpiralAtAngleToGeo_AllSets\n\n\n");

    return suite;
}

TestSet testPtsOnSpiralOnTanThruPt_Set1(){

	TestSet set = newTestSet("testPointToSpiralTangents_Set1");
	double latxExp, lonxExp, startRad, endRad, startAz, endAz, az1, az2, crs1, crs2, crs13, dist13, crs23, dist23, dist, rad;
	ArcDirection direction;
	int i, j, k;
	LLPoint centerPoint, pt1, pt2, testPt;
	Spiral spiral;
	double eps = 1.0e-20;
	double tol = 1.37e-9;
	double angle[8] = {30,45,60,75,90,105,120,135};
	long err = 0;

	double az12, az21;

	double total;

	int count1 = 0, count2 = 0, counter = 0;
	int testNum = 0, successCount = 0, errCount = 0, failCount = 0, setupCount = 0, unverifiedCount = 0;

	testNum = 125;

	srand(04012011);

	for (i=0;i<testNum;i++) {

		latxExp = randLat();
		lonxExp = randLon();

		err |= createPt(&centerPoint, latxExp*DEG2RAD, lonxExp*DEG2RAD);
		if ((rand() % 100) > 49) {
			direction = CLOCKWISE;
		} else {
			direction = COUNTERCLOCKWISE;
		}
		startRad = 20 + 0.01 * randDist(); //0-54 nm
		endRad = 20 + 0.01 * randDist(); //0-54 nm
		startAz = randAzimuth();
		endAz = fmod(startAz + direction * 270, 360);

		az1 = fmod(startAz + direction * 45, 360);

		err |= createSpiral(&spiral, centerPoint, startRad, endRad, startAz * DEG2RAD, endAz * DEG2RAD, direction, eps);

		err |= ptOnSpiral(spiral, az1 * DEG2RAD, &pt1, eps);
		err |= spiralTanCrs(spiral, az1 * DEG2RAD, &crs1);

		for (j=0;j<8;j++) {
			LLPointSet pointSet = createPtSet();

			az2 = fmod(az1 + direction * angle[j], 360);
			err |= ptOnSpiral(spiral, az2 * DEG2RAD, &pt2, eps);
			err |= spiralTanCrs(spiral, az2 * DEG2RAD, &crs2);
			crs2 = fmod(crs2 + M_PI, 2 * M_PI);

			err |= spiralRadius(spiral, az1 * DEG2RAD, &rad);

			//printf("Rad 1:  %f\n", rad);

			err |= spiralRadius(spiral, az2 * DEG2RAD, &rad);

			//printf("Rad 2:  %f\n", rad);

			err |= crsIntx(pt1, crs1, &crs13, &dist13, pt2, crs2, &crs23, &dist23, &testPt, tol, eps);

			//printf("Course Intersect:  %f %f\n", crs13, dist13);
			//printf("Course Intersect:  %f %f\n", crs23, dist23);

			err |= inverse(pt1, testPt, &az12, &az21, &dist, eps);
			//printf("Inverse:  %f %f\n", az21, dist);
			err |= inverse(pt2, testPt, &az12, &az21, &dist, eps);
			//printf("Inverse:  %f %f\n", az21, dist);

			err |= invDist(spiral.centerPoint, testPt, &rad, eps);

			//printf("Test Rad:  %f\n", rad);

			if (err) {
				setupCount++;
				err = 0;
			}

			err |= ptsOnSpiralOnTanThruPt(spiral, testPt, &pointSet, tol, eps);

			if (err) {
				//printf("Error:  %s\n", WGS84FormatErrorMessage(err));
				//printf("i, Angle:  %i %f\n", i, angle[j]);
				//printf("Input Spiral:  %f %f %f %f %f %i\n", spiral.startAz, spiral.endAz, spiral.startRadius, spiral.endRadius, spiral.growthRate, i);
				errCount++;
			}

			double d, x, R;

			R = 3440;

			for (k=0;k<pointSet.length;k++) {
				LLPoint point;
				point.latitude = pointSet.elements[k]->latitude;
				point.longitude = pointSet.elements[k]->longitude;
				err |= invDist(spiral.centerPoint, pt1, &rad, eps);
				//printf("Distance to Point 1:  %e\n", dist);

				d = R * acos(cos((rad + tol) / R) / cos((rad - tol) / R));
				x = R * acos(cos(d / R) * cos(tol / R));

				err |= invDist(point, pt2, &dist, eps);
				//printf("Distance to Point 2:  %e\n", dist);

				if (ptsAreSame(point, pt1, tol)) {
					count1++;
				} else {
					if (ptsAreSame(point, pt2, tol)) {
						count2++;
					}
				}
			}
			total += x;
			//printf("X:  %e\n", x);
			//printf("Counts:  %i %i\n", count1, count2);
			if (count1 + count2 >= 2) {
				successCount++;
			} else {
				if ((count1 + count2) > 2) {

				}
				//printf("i, j:  %i %i\n", i, j);
				//printf("Spiral:  %f\n", spiral.growthRate);
				//printf("Counts:  %i %i\n", count1, count2);
				//printf("Point Set:  %i\n", pointSet.length);
				failCount++;
				for (k=0;k<pointSet.length;k++) {
					LLPoint point;
					point.latitude = pointSet.elements[k]->latitude;
					point.longitude = pointSet.elements[k]->longitude;
					//printf("LL:  %f %f\n", point.latitude / DEG2RAD, point.longitude / DEG2RAD);
					err |= invDist(point, pt1, &dist, eps);
					err |= invDist(point, pt2, &dist23, eps);
					printf("Distances :  %e %e\n", dist, dist23);

					//printf("Distance to Point 2:  %e\n", dist);
				}

			}

			counter++;

			count1 = 0;
			count2 = 0;

			/*
			printf("Courses:  %f %f\n", crs3, crs4);

			err |= inverse(pt1, pt3, &az12, &az21, &dist, eps);
			printf("Point:  %f %e\n", az12, dist);
			err |= GetRadius(spiral, az12, &dist);
			printf("Radius:  %f\n", dist);

			err |= inverse(pt1, pt4, &az12, &az21, &dist, eps);
			printf("Point:  %f %e\n", az12, dist);
			err |= GetRadius(spiral, az12, &dist);
			printf("Radius:  %f\n", dist);

			err |= inverse(pt2, pt3, &az12, &az21, &dist, eps);
			printf("Point:  %f %e\n", az12, dist);
			err |= GetRadius(spiral, az12, &dist);
			printf("Radius:  %f\n", dist);

			err |= inverse(pt2, pt4, &az12, &az21, &dist, eps);
			printf("Point:  %f %e\n", az12, dist);
			err |= GetRadius(spiral, az12, &dist);
			printf("Radius:  %f\n", dist);

			printf("Pointset Length:  %i\n", pointSet.length);
			printf("Error: %s\n", WGS84FormatErrorMessage(err));
			*/
		}
	}

	set.errors = errCount;
	set.fail = failCount;
	set.pass = successCount;
	set.testCases = 1000;
	set.setupFailures = setupCount;
	set.unverified = unverifiedCount;

	return set;
	//printf("Total Tests Run:  %i\n", counter);
	//printf("Point To Spiral Tangents:  %i/1000\n", successCount);
	//printf("Average ROA:  %e\n", total / 1000);
	//printf("Errors:  %i\n", errCount);
}

TestSuite testPtsOnSpiralOnTanThruPt_AllSets() {

	TestSuite suite;
	TestSet set1;

    printf("\nStart testPtsOnSpiralOnTanThruPt_AllSets\n");

    suite = newTestSuite("testPtsOnSpiralOnTanThruPt_AllSets");

    set1 = testPtsOnSpiralOnTanThruPt_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("Finish testPtsOnSpiralOnTanThruPt_AllSets\n\n\n");

    return suite;
}

TestSet testGeoTanToTwoSpirals_Set1(){

	TestSet set = newTestSet("testLineTangentToTwoSpirals_Set1");
	double latxExp, lonxExp, startRad, endRad, startAz, endAz;
	double testAz, spCrs, az12, az21, dist, rad, azDiff;
	ArcDirection direction;
	int i;
	LLPoint centerPoint, pt1, pt2, center2;
	Spiral spiral, spiral2;
	double eps = 1.0e-20;
	double tol = 1.37e-9;
	ErrorSet err = 0;

	int count1 = 0, count2 = 0, k;

	int testNum = 0, successCount = 0, errCount = 0, failCount = 0, setupCount = 0, unverifiedCount = 0;

	double az1, az2;

	testNum = 1000;

	srand(04012011);

	for (i=0;i<testNum;i++) {

		LLPointSet pointSet = createPtSet();

		latxExp = randLat();
		lonxExp = randLon();

		err |= createPt(&centerPoint, latxExp*DEG2RAD, lonxExp*DEG2RAD);
		if ((rand() % 100) > 49) {
			direction = CLOCKWISE;
		} else {
			direction = COUNTERCLOCKWISE;
		}
		startRad = 20 + 0.01 * randDist(); //0-54 nm
		endRad = 20 + 0.01 * randDist(); //0-54 nm
		startAz = randAzimuth();
		endAz = fmod(startAz + direction * 270, 360);
		testAz = fmod(startAz + direction * 135, 360) * DEG2RAD;

		err |= createSpiral(&spiral, centerPoint, startRad, endRad, startAz * DEG2RAD, endAz * DEG2RAD, direction, eps);

		err |= ptOnSpiral(spiral, testAz, &pt1, eps);
		err |= spiralTanCrs(spiral, testAz, &spCrs);
		az1 = testAz;
		err |= inverse(centerPoint, pt1, &az12, &az21, &rad, eps);
		azDiff = AzDiff(spCrs, az21);

		err |= direct(pt1, spCrs, randDist() * .1 + 30, &pt2, eps);
		err |= inverse(pt1, pt2, &az12, &az21, &dist, eps);
		//printf("Distance between spirals:  %f\n", dist);
		err |= direct(pt2, az21 + M_PI + azDiff, rad, &center2, eps);

		err |= inverse(center2, pt2, &az12, &az21, &dist, eps);
		//printf("Point 2 Az:  %f\n", az12);
		az2 = az12;
		err |= createSpiral(&spiral2, center2, startRad, endRad, az12 - direction * 135 * DEG2RAD, az12 + direction * 135 * DEG2RAD, direction, eps);

		err |= spiralTanCrs(spiral2, az12, &spCrs);
		//printf("Tangent:  %f\n", spCrs);
		azDiff = AzDiff(spCrs, az21);
		//printf("Az Difference:  %f\n", azDiff);
		//printf("Az/Tangent/Rad:  %f %f %f\n", az12, spCrs, rad);

		err |= invCrs(pt1, pt2, &az12, &az21, eps);
		//printf("Tangents:  %f %f\n", az12, az21);

		if (err) {
			setupCount++;
			err = 0;
		}

		err |= geoTanToTwoSpirals(spiral, spiral2, &pointSet, tol, eps);

		if (err) {
			errCount++;
		}

		LLPoint tempStart, tempEnd, mc1, mc2;
		Geodesic tempGeo;

		err |= inverse(centerPoint, center2, &az12, &az21, &dist, eps);
		err |= direct(centerPoint, az12 - M_PI / 2, 5, &tempStart, eps);
		err |= direct(center2, az21 + M_PI / 2, 5, &tempEnd, eps);
		err |= createGeo(&tempGeo, tempStart, tempEnd, INFINITE, eps);

		err |= spiralMidChord(spiral, tempGeo, &mc1, tol, eps);
		err |= spiralMidChord(spiral2, tempGeo, &mc2, tol, eps);

		err |= invDist(mc1, pt1, &dist, eps);
		err |= invDist(mc1, pt1, &az12, eps);

		//printf("Distances:  %f %f\n", dist, az12);
		err |= invCrs(mc1, mc2, &az12, &az21, eps);
		//printf("Course\Back Course:  %f %f\n", az12, az21);
		//printf("Pointset:  %i\n", pointSet.length);

		for (k=0;k<pointSet.length;k++) {
			LLPoint point;
			err |= inverse(point, mc1, &az12, &az21, &dist, eps);
			err |= inverse(point, mc2, &az12, &az21, &dist, eps);
			point.latitude = pointSet.elements[k]->latitude;
			point.longitude = pointSet.elements[k]->longitude;

			if (ptsAreSame(point, pt1, tol)) {
				count1++;
			} else {
				if (ptsAreSame(point, pt2, tol)) {
					count2++;
				}
			}
		}
		if (count1 + count2 == 2) {
			successCount++;
		} else {
			failCount++;
			err |= invCrs(spiral.centerPoint, spiral2.centerPoint, &az12, &az21, eps);
			printf("Azs:  %f %f %f\n", az1, az2, az12);
			for (k=0;k<pointSet.length;k++) {
				LLPoint point;
				point.latitude = pointSet.elements[k]->latitude;
				point.longitude = pointSet.elements[k]->longitude;
				if (k == 0) {
					err |= inverse(spiral.centerPoint, point, &az12, &az21, &dist, eps);
					err |= invDist(point, pt1, &dist, eps);
					printf("Distance to Point 1:  %e\n", dist);
				} else {
					err |= inverse(spiral2.centerPoint, point, &az12, &az21, &dist, eps);
					err |= invDist(point, pt2, &dist, eps);
					printf("Distance to Point 2:  %e\n", dist);
				}
			}
		}
		err = 0;
		count1 = 0;
		count2 = 0;
	}

	set.errors = errCount;
	set.fail = failCount;
	set.pass = successCount;
	set.testCases = 1000;
	set.setupFailures = setupCount;
	set.unverified = unverifiedCount;

	return set;

	//printf("Line Tangent To Two Spirals:  %i/1000\n", successCount);
}

TestSuite testGeoTanToTwoSpirals_AllSets() {

	TestSuite suite;
	TestSet set1;

    printf("\nStart testGeoTanToTwoSpirals_AllSets\n");

    suite = newTestSuite("testGeoTanToTwoSpirals_AllSets");

    set1 = testGeoTanToTwoSpirals_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("Finish testGeoTanToTwoSpirals_AllSets\n\n\n");

    return suite;
}

TestSet testProjectToSpiral_Set1(){

	TestSet set = newTestSet("testProjectToSpiral_Set1");
	double latxExp, lonxExp, startRad, endRad, startAz, endAz, testAz, spCrs, ptCrs, rad, ptDist;
	ArcDirection direction;
	int i;
	LLPoint centerPoint, spPt, testPt, perpPt;
	Spiral spiral;
	double eps = 1.0e-20;
	double tol = 1.37e-9;
	long err = 0;

	int testNum = 0, successCount = 0, errCount = 0, failCount = 0, setupCount = 0, unverifiedCount = 0;

	double az12, az21, dist;

	testNum = 1000;

	srand(04012011);

	for (i=0;i<testNum;i++) {

		//printf("\n%i\n", i);

		//printf("Test %i\n\n", i);

		latxExp = randLat();
		lonxExp = randLon();

		err |= createPt(&centerPoint, latxExp*DEG2RAD, lonxExp*DEG2RAD);
		if ((rand() % 100) > 49) {
			direction = CLOCKWISE;
		} else {
			direction = COUNTERCLOCKWISE;
		}
		startRad = 20 + 0.01 * randDist(); //0-54 nm
		endRad = 20 + 0.01 * randDist(); //0-54 nm
		startAz = randAzimuth();
		endAz = fmod(startAz + direction * 270, 360);
		testAz = fmod(startAz + direction * 135, 360) * DEG2RAD;

		err |= createSpiral(&spiral, centerPoint, startRad, endRad, startAz * DEG2RAD, endAz * DEG2RAD, direction, eps);

		err |= ptOnSpiral(spiral, testAz, &spPt, eps);
		err |= spiralTanCrs(spiral, testAz, &spCrs);
		err |= spiralRadius(spiral, testAz, &rad);
		//printf("Point Radius/Az:  %f %f\n", rad, testAz);
		ptCrs = spCrs - spiral.dir * M_PI / 2;
		ptDist = randDist() - (rad * .95);
		//ptDist = rad;
		err |= direct(spPt, ptCrs, ptDist, &testPt, eps);

		if (err) {
			setupCount++;
			err = 0;
		}

		err |= projectToSpiral(spiral, testPt, &perpPt, tol, eps);

		err |= inverse(spiral.centerPoint, perpPt, &az12, &az21, &dist, eps);
		if (err) {
			//printf("Input Spiral:  %f %f %f %f %f %i\n", spiral.startAz, spiral.endAz, spiral.startRadius, spiral.endRadius, spiral.growthRate, i);
			//printf("Point LL/Distance:  %f %f %f\n", testPt.latitude / DEG2RAD, testPt.longitude / DEG2RAD, dist);
			errCount++;
		}
		if (ptsAreSame(spPt, perpPt, tol)) {
		//if (ptsAreSame(spPt, perpPt, 1e-7)) {
			successCount++;
		} else {
			failCount++;
			printf("Input Spiral:  %f %f %f %f %f %i\n", spiral.startAz, spiral.endAz, spiral.startRadius, spiral.endRadius, spiral.growthRate, i);
			err |= inverse(spPt, perpPt, &az12, &az21, &dist, eps);
			printf("Distance from Point:  %e\n", dist);
			err |= inverse(spPt, testPt, &az12, &az21, &dist, eps);
			printf("Distance from Spiral:  %f\n", dist);
			printf("LLs:  %f %f %f %f\n", spPt.latitude / DEG2RAD, spPt.longitude / DEG2RAD, perpPt.latitude / DEG2RAD, perpPt.longitude / DEG2RAD);
			printf("Point Distance:  %f\n", ptDist);
		}
	}


	set.errors = errCount;
	set.fail = failCount;
	set.pass = successCount;
	set.testCases = 1000;
	set.setupFailures = setupCount;
	set.unverified = unverifiedCount;

	return set;

	//printf("Project To Spiral:  %i/1000\n", successCount);
	//printf("Errors:  %i\n", errCount);
}

TestSuite testProjectToSpiral_AllSets() {

	TestSuite suite;
	TestSet set1;

    printf("\nStart testProjectToSpiral_AllSets\n");

    suite = newTestSuite("testProjectToSpiral_AllSets");

    set1 = testProjectToSpiral_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("Finish testProjectToSpiral_AllSets\n\n\n");

    return suite;
}


/*void testSpiralBug(){

	//double DEG2RAD = M_PI / 180.0;
	long err = 0;
	double eps = 1.0e-20;
	double tol = 1.37e-9;
	double startCrs, endCrs, startRadius, endRadius, angle;
	LLPoint spiralCenter, lineStart, lineEnd;
	Spiral spiral;
	Geodesic inputGeo;
	Geodesic tangentGeo;

	err |= createPt(&lineStart, 25.808469346293656*DEG2RAD,
-80.16772673371962*DEG2RAD);
	err |= createPt(&lineEnd, 25.709659058763766*DEG2RAD,
-80.11523923849838*DEG2RAD);

	err |= createPt(&spiralCenter, 25.77513248597836*DEG2RAD,
-80.1385420066993*DEG2RAD);

	startCrs = 357.40843153496894*DEG2RAD;
	endCrs = 357.40843153496894*DEG2RAD;
	startRadius = 4.671160165218935;
	endRadius = 11.368944903201141;
	angle = 30.0*DEG2RAD;

	err |= createSpiral(&spiral, spiralCenter, startRadius, endRadius, startCrs,
endRadius, -1, eps);
	err |= createGeo(&inputGeo,  lineStart, lineEnd, SEGMENT, eps);

	if (err) printf("\nError(s) occurred during setup: 0x%lx", err);

	err |= geoTanToSpiralAtAngleToGeo(spiral, inputGeo, angle,
&tangentGeo, tol, eps);

	if (err) printf("Error: %s\n", formatErrorMessage(err));
	else printf("Finished\n");

}
*/

void testSpiralArcIntx_Set1(){

	LLPoint arcStartPoint, arcEndPoint, centerPoint;
	double startRad, endRad, startAz, endAz, dist, totDist, spAngle, epAngle;
	ArcDirection direction;
	double eps = 1.0e-20;
	double tol = 1.37e-9;
	long err = 0;
	int i, j, testNum, successCount = 0, spcounter = 0, epcounter = 0, similarCount = 0, simFailCount = 0;
	Arc arc;
	Spiral spiral;

	double az12, az21, geoAz, rad;

	double latxExp, lonxExp;

	testNum = 1000;

	srand(04012011);

	double minRad, maxRad;

	for (i=0;i<testNum;i++) {

		err = 0;

		LLPointSet pointSet = createPtSet();

		latxExp = randLat();
		lonxExp = randLon();

		err |= createPt(&centerPoint, latxExp*DEG2RAD, lonxExp*DEG2RAD);
		if ((rand() % 100) > 49) {
			direction = CLOCKWISE;
		} else {
			direction = COUNTERCLOCKWISE;
		}
		startRad = 30 + 0.01 * randDist() / 2; //0-54 nm
		endRad = 30 + 0.01 * randDist() / 2; //0-54 nm
		startAz = randAzimuth();
		endAz = startAz + direction * 270;

		spAngle = rand() % 270;
		epAngle = rand() % 270;

		//printf("Azs:  %f %f %f %f\n", startAz, endAz, spAngle, epAngle);

		//printf("Current Test:  %i\n", i);

		//These values should be generated randomly
		//err |= createPt(&geoStartPoint, 25.694592395967856*DEG2RAD, -80.10724447560116*DEG2RAD);
		//err |= createPt(&geoEndPoint, 25.701834182442678*DEG2RAD, -80.09061424611792*DEG2RAD);

		err |= createSpiral(&spiral, centerPoint, startRad, endRad, startAz * DEG2RAD, endAz * DEG2RAD, direction, eps);

		//printf("Spiral Growthrate:  %f\n", spiral.growthRate);

		//printf("Angles:  %f %f\n", (startAz + direction * spAngle) * DEG2RAD, (startAz + direction * epAngle) * DEG2RAD);

		err |= ptOnSpiral(spiral, (startAz + direction * spAngle) * DEG2RAD, &arcStartPoint, eps);
		err |= ptOnSpiral(spiral, (startAz + direction * epAngle) * DEG2RAD, &arcEndPoint, eps);

		//printf("LL:  %f %f\n", arcStartPoint.latitude, arcStartPoint.longitude);
		//printf("LL:  %f %f\n", arcEndPoint.latitude, arcEndPoint.longitude);

		if (err) {
			printf("Error:  %s\n", formatErrorMessage(err));
		}

		//printf("Azs:  %f %f\n", (startAz - direction * spAngle) * DEG2RAD, (startAz - direction * epAngle) * DEG2RAD);

		geoAz = randAzimuth();
		err |= arcFromStartAndEnd(arcStartPoint, geoAz * DEG2RAD, arcEndPoint, &arc, tol, eps);
		err |= createArc(&arc, arc.centerPoint, arc.startPoint, arc.startPoint, arc.dir, tol, eps);
		//(startAz + direction * spAngle) * DEG2RAD
		err |= inverse(centerPoint, arcStartPoint, &az12, &az21, &dist, eps);
		//printf("Point 1:  %f %f\n", az12, dist);
		err |= inverse(centerPoint, arcEndPoint, &az12, &az21, &dist, eps);
		//printf("Point 2:  %f %f\n", az12, dist);

		//printf("Start Point:  %f %f\n", geo.startPoint.latitude / DEG2RAD, geo.startPoint.longitude / DEG2RAD);
		//printf("End Point:  %f %f\n", geo.endPoint.latitude / DEG2RAD, geo.endPoint.longitude / DEG2RAD);

		err |= spiralArcIntx(spiral, arc, &pointSet, tol, eps);

		spcounter = 0;
		epcounter = 0;

		totDist = 0;
		minRad = spiral.startRadius;
		maxRad = spiral.endRadius;
		if (spiral.startRadius > spiral.endRadius) {
			minRad = spiral.endRadius;
			maxRad = spiral.startRadius;
		}
		err |= invDist(spiral.centerPoint, arc.centerPoint, &dist, eps);
		if ((dist + arc.radius > minRad) && (dist + arc.radius < maxRad)  && (arc.radius > maxRad - minRad)) {
			similarCount++;
		}

		//printf("Start Point:  %f %f\n", spiral.startPoint.latitude / DEG2RAD, spiral.startPoint.longitude / DEG2RAD);
		//printf("End Point:  %f %f\n", spiral.endPoint.latitude / DEG2RAD, spiral.endPoint.longitude / DEG2RAD);
		for (j=0;j<pointSet.length;j++) {
			LLPoint point;
			point.latitude = pointSet.elements[j]->latitude;
			point.longitude = pointSet.elements[j]->longitude;
			err |= invDist(point, arcStartPoint, &dist, eps);
			totDist = totDist + dist;
			err |= invDist(point, arcEndPoint, &dist, eps);
			totDist = totDist + dist;
			if (ptsAreSame(point, arcStartPoint, 100 * tol)) {
				spcounter = 1;
			}
			if (ptsAreSame(point, arcEndPoint, 100 * tol)) {
				epcounter = 1;
			}
		}

		//printf("Pointset Length:  %i\n", pointSet.length);

		if ((spcounter + epcounter < 2)) {
			err |= inverse(spiral.centerPoint, arc.centerPoint, &az12, &az21, &dist, eps);
			minRad = spiral.startRadius;
			maxRad = spiral.endRadius;
			if (spiral.startRadius > spiral.endRadius) {
				minRad = spiral.endRadius;
				maxRad = spiral.startRadius;
			}
			if ((dist + arc.radius > minRad) && (dist + arc.radius < maxRad)  && (arc.radius > maxRad - minRad)) {
				printf("Similar Shapes:  %f %f %f %f\n", dist, arc.radius, minRad, maxRad);
				simFailCount++;
			} else {
				printf("#####################################################################\n");
				printf("Pointset Length:  %i\n", pointSet.length);
				/*printf("cp = CoordinatePoint([%f %f 0]);\n", spiral.centerPoint.latitude / (M_PI / 180), spiral.centerPoint.longitude / (M_PI / 180));
				err |= inverse(spiral.centerPoint, spiral.startPoint, &az12, &az21, &dist, eps);
				//printf("sp = movePointByCourse(cp, %f, 0, %f);\n", az12, dist);
				printf("sp = CoordinatePoint([%f %f 0]);\n", spiral.startPoint.latitude / (M_PI / 180), spiral.startPoint.longitude / (M_PI / 180));
				err |= inverse(spiral.centerPoint, spiral.endPoint, &az12, &az21, &dist, eps);
				//printf("ep = movePointByCourse(cp, %f, 0, %f);\n", az12, dist);
				printf("ep = CoordinatePoint([%f %f 0]);\n", spiral.endPoint.latitude / (M_PI / 180), spiral.endPoint.longitude / (M_PI / 180));
				err |= inverse(spiral.centerPoint, arc.centerPoint, &az12, &az21, &dist, eps);
				//printf("ac = movePointByCourse(cp, %f, 0, %f);\n", az12, dist);
				printf("ac = CoordinatePoint([%f %f 0]);\n", arc.centerPoint.latitude / (M_PI / 180), arc.centerPoint.longitude / (M_PI / 180));
				err |= inverse(spiral.centerPoint, arcStartPoint, &az12, &az21, &dist, eps);
				//printf("as = movePointByCourse(cp, %f, 0, %f);\n", az12, dist);
				printf("as = CoordinatePoint([%f %f 0]);\n", arc.startPoint.latitude / (M_PI / 180), arc.startPoint.longitude / (M_PI / 180));
				err |= inverse(spiral.centerPoint, arcEndPoint, &az12, &az21, &dist, eps);
				//printf("ae = movePointByCourse(cp, %f, 0, %f);\n", az12, dist);
				printf("sp1 = Spiral(cp, sp, ep, %i);\n", spiral.dir * -1);
				printf("dist = xyDistance(ac, as);\n");
				printf("arc = Arc(ac, dist);\n");
				printf("drawFlat(sp1);\ndrawFlat(arc);\n");
				*/
				err |= spiralRadius(spiral, az12, &rad);
				//printf("Radius:  %f\n", rad);
			}
			err |= spiralRadius(spiral, az12, &rad);
			//printf("Inside:  %i\n", rad - dist > 0);
			//printf("Failure!!!\n");
			printf("Azimuths:  %f %f\n", spAngle, epAngle);
			printf("Azs/Dir:  %f %f\n", fmod((startAz + direction * spAngle) * DEG2RAD + 2*M_PI, 2*M_PI), fmod((startAz + direction * epAngle) * DEG2RAD + 2*M_PI, 2*M_PI));
			//printf("Radius:  %f\n", arc.radius);
			//printf("Error!  %s", formatErrorMessage(err));
			//printf("Direction:  %i\n", direction);
			printf("Failure:  %i\n", i);
			//printf("Counters:  %i %i\n", spcounter, epcounter);

			for (j=0;j<pointSet.length;j++) {
				LLPoint point;
				point.latitude = pointSet.elements[j]->latitude;
				point.longitude = pointSet.elements[j]->longitude;
				err |= invDist(point, arcStartPoint, &dist, eps);
				printf("Start Point:  %e\n", dist);
				err |= invDist(point, arcEndPoint, &dist, eps);
				printf("End Point:  %e\n", dist);
				err |= inverse(spiral.centerPoint, point, &az12, &az21, &dist, eps);
				printf("Azimuth/Distance:  %f %f\n", az12, dist);
				err |= spiralRadius(spiral, az12, &rad);
				//printf("Radius:  %f\n", rad);
				err |= invDist(point, arc.centerPoint, &dist, eps);
				//printf("Arc Dist:  %f\n", dist);
			}
			printf("Spiral:  %f %f %f %f\n", spiral.startRadius, spiral.endRadius, spiral.startAz, spiral.endAz);
			err |= inverse(spiral.centerPoint, arc.centerPoint, &az12, &az21, &dist, eps);
			printf("Arc Rad/Dist/Az:  %f %f %f\n", arc.radius, dist, az12);
			if (pointSet.length > 0) {
				//printf("Spiral:  %f %f %f %f\n", spiral.startRadius, spiral.endRadius, spiral.startAz, spiral.endAz);
				//printf("Start Point:  %f %f\n", spiral.startPoint.latitude / DEG2RAD, spiral.startPoint.longitude / DEG2RAD);
				//printf("End Point:  %f %f\n", spiral.endPoint.latitude / DEG2RAD, spiral.endPoint.longitude / DEG2RAD);
				//printf("Start/End:  %i %i\n", spcounter, epcounter);
				//printf("Distances:  %f %f\n", totDist, geo.length * pointSet.length);
				//printf("Point Set:  %i\n", pointSet.length);
				//printf("Failure\n\n");
			} else {
				//printf("No Points Returned\n\n");
			}
			//printf("\n");
		} else {
			successCount++;
		}

		//printf("\n");

		clearPtSet(&pointSet);

		//val = ptsAreSame(testPt, midChord, tol);
		//printf("%i\n", val);

	}

	printf("Spiral Arc Intersect:  %i/1000\n", successCount);
	printf("Similar Shapes:  %i %i\n", similarCount, simFailCount);

}

void testSpiralIntx_Set1(){

	LLPoint arcStartPoint, arcEndPoint, centerPoint;
	double startRad, endRad, endAz, startAz, dist, totDist, spAngle, epAngle, dtheta;
	ArcDirection direction;
	double eps = 1.0e-20;
	double tol = 1.37e-9;
	long err = 0;
	int i, j, testNum, successCount = 0, spcounter = 0, epcounter = 0;
	Arc arc;
	Spiral spiral, spiral2;

	double az12, az21, geoAz, rad;

	double latxExp, lonxExp;

	testNum = 1000;

	srand(04012011);

	for (i=0;i<testNum;i++) {

		err = 0;

		LLPointSet pointSet = createPtSet();

		latxExp = randLat();
		lonxExp = randLon();

		err |= createPt(&centerPoint, latxExp*DEG2RAD, lonxExp*DEG2RAD);
		if ((rand() % 100) > 49) {
			direction = CLOCKWISE;
		} else {
			direction = COUNTERCLOCKWISE;
		}
		startRad = 30 + 0.01 * randDist() / 2; //0-54 nm
		endRad = 30 + 0.01 * randDist() / 2; //0-54 nm
		startAz = randAzimuth();
		endAz = startAz + direction * 270;

		spAngle = rand() % 270;
		epAngle = rand() % 270;

		//printf("Azs:  %f %f %f %f\n", startAz, endAz, spAngle, epAngle);

		//printf("Current Test:  %i\n", i);

		//These values should be generated randomly
		//err |= createPt(&geoStartPoint, 25.694592395967856*DEG2RAD, -80.10724447560116*DEG2RAD);
		//err |= createPt(&geoEndPoint, 25.701834182442678*DEG2RAD, -80.09061424611792*DEG2RAD);

		err |= createSpiral(&spiral, centerPoint, startRad, endRad, startAz * DEG2RAD, endAz * DEG2RAD, direction, eps);

		//printf("Spiral Growthrate:  %f\n", spiral.growthRate);

		//printf("Angles:  %f %f\n", (startAz + direction * spAngle) * DEG2RAD, (startAz + direction * epAngle) * DEG2RAD);

		err |= ptOnSpiral(spiral, (startAz + direction * spAngle) * DEG2RAD, &arcStartPoint, eps);
		err |= ptOnSpiral(spiral, (startAz + direction * epAngle) * DEG2RAD, &arcEndPoint, eps);

		//printf("LL:  %f %f\n", arcStartPoint.latitude, arcStartPoint.longitude);
		//printf("LL:  %f %f\n", arcEndPoint.latitude, arcEndPoint.longitude);

		if (err) {
			printf("Error:  %s\n", formatErrorMessage(err));
		}

		//printf("Azs:  %f %f\n", (startAz - direction * spAngle) * DEG2RAD, (startAz - direction * epAngle) * DEG2RAD);

		geoAz = randAzimuth();
		err |= arcFromStartAndEnd(arcStartPoint, geoAz * DEG2RAD, arcEndPoint, &arc, tol, eps);

		dtheta = rand() % 45;

		err |= ptOnSpiral(spiral, (startAz + direction * spAngle + dtheta) * DEG2RAD, &arcStartPoint, eps);
		err |= inverse(arc.centerPoint, arcStartPoint, &startAz, &az21, &startRad, eps);
		err |= inverse(arc.centerPoint, arcEndPoint, &endAz, &az21, &endRad, eps);

		//printf("Azs:  %f %f\n", startAz, endAz);

		if ((rand() % 100) > 49) {
			direction = CLOCKWISE;
		} else {
			direction = COUNTERCLOCKWISE;
		}

		err |= createSpiral(&spiral2, arc.centerPoint, startRad, endRad, startAz, endAz, direction, eps);

		//err |= createSpiralSection(tempSp);


		//(startAz + direction * spAngle) * DEG2RAD
		err |= inverse(centerPoint, arcStartPoint, &az12, &az21, &dist, eps);
		//printf("Point 1:  %f %f\n", az12, dist);
		err |= inverse(centerPoint, arcEndPoint, &az12, &az21, &dist, eps);
		//printf("Point 2:  %f %f\n", az12, dist);

		//printf("Start Point:  %f %f\n", geo.startPoint.latitude / DEG2RAD, geo.startPoint.longitude / DEG2RAD);
		//printf("End Point:  %f %f\n", geo.endPoint.latitude / DEG2RAD, geo.endPoint.longitude / DEG2RAD);

		err |= spiralIntx(spiral, spiral2, &pointSet, tol, eps);

		spcounter = 0;
		epcounter = 0;

		totDist = 0;

		//printf("Spiral 1:  %f %f %f %f %f %f\n", spiral.startRadius, spiral.endRadius, spiral.startAz, spiral.endAz, spiral.growthRate, spiral.subtendedAngle);
		//printf("Spiral 2:  %f %f %f %f %f %f\n", spiral2.startRadius, spiral2.endRadius, spiral2.startAz, spiral2.endAz, spiral2.growthRate, spiral2.subtendedAngle);

		//printf("Start Point:  %f %f\n", spiral.startPoint.latitude / DEG2RAD, spiral.startPoint.longitude / DEG2RAD);
		//printf("End Point:  %f %f\n", spiral.endPoint.latitude / DEG2RAD, spiral.endPoint.longitude / DEG2RAD);
		for (j=0;j<pointSet.length;j++) {
			LLPoint point;
			point.latitude = pointSet.elements[j]->latitude;
			point.longitude = pointSet.elements[j]->longitude;
			err |= invDist(point, arcStartPoint, &dist, eps);
			totDist = totDist + dist;
			err |= invDist(point, arcEndPoint, &dist, eps);
			totDist = totDist + dist;
			if (ptsAreSame(point, arcStartPoint, 100 * tol)) {
				spcounter = 1;
			}
			if (ptsAreSame(point, arcEndPoint, 100 * tol)) {
				epcounter = 1;
			}
		}

		//printf("Pointset Length:  %i\n", pointSet.length);

		if ((spcounter + epcounter < 2)) {
			err |= inverse(spiral.centerPoint, arc.centerPoint, &az12, &az21, &dist, eps);
			err |= spiralRadius(spiral, az12, &rad);
			//printf("Inside:  %i\n", rad - dist > 0);
			printf("Pointset Length:  %i\n", pointSet.length);
			//printf("Failure!!!\n");
			//printf("Azimuths:  %f %f\n", spAngle, epAngle);
			//printf("Azs/Dir:  %f %f %f %f %f %i\n", startAz * DEG2RAD, endAz * DEG2RAD, fmod((startAz + direction * spAngle) * DEG2RAD + 2*M_PI, 2*M_PI), fmod((startAz + direction * epAngle) * DEG2RAD + 2*M_PI, 2*M_PI), geoAz, direction);
			//printf("Radius:  %f\n", arc.radius);
			//printf("Error!  %s", formatErrorMessage(err));
			//printf("Direction:  %i\n", direction);
			//printf("Failure:  %i\n", i);
			//printf("Counters:  %i %i\n", spcounter, epcounter);
			/*printf("cp = CoordinatePoint([%f %f 0]);\n", spiral.centerPoint.latitude / (M_PI / 180), spiral.centerPoint.longitude / (M_PI / 180));
			err |= inverse(spiral.centerPoint, spiral.startPoint, &az12, &az21, &dist, eps);
			printf("sp = movePointByCourse(cp, %f, 0, %f);\n", az12, dist);
			err |= inverse(spiral.centerPoint, spiral.endPoint, &az12, &az21, &dist, eps);
			printf("ep = movePointByCourse(cp, %f, 0, %f);\n", az12, dist);
			err |= inverse(spiral.centerPoint, spiral2.centerPoint, &az12, &az21, &dist, eps);
			printf("ac = movePointByCourse(cp, %f, 0, %f);\n", az12, dist);
			err |= inverse(spiral.centerPoint, arcStartPoint, &az12, &az21, &dist, eps);
			printf("as = movePointByCourse(cp, %f, 0, %f);\n", az12, dist);
			err |= inverse(spiral.centerPoint, arcEndPoint, &az12, &az21, &dist, eps);
			printf("ae = movePointByCourse(cp, %f, 0, %f);\n", az12, dist);
			err |= spiralRadius(spiral, az12, &rad);
			printf("sp1 = Spiral(cp, sp, ep, %i);\n", spiral.dir * -1);
			printf("sp2 = Spiral(ac, as, ae, %i);\n", spiral2.dir * -1);
			printf("drawFlat(sp1);\n");
			printf("drawFlat(sp2);\n");
			*/
			//printf("Radius:  %f\n", rad);

			for (j=0;j<pointSet.length;j++) {
				LLPoint point;
				point.latitude = pointSet.elements[j]->latitude;
				point.longitude = pointSet.elements[j]->longitude;
				err |= invDist(point, arcStartPoint, &dist, eps);
				printf("Start Point:  %e\n", dist);
				err |= invDist(point, arcEndPoint, &dist, eps);
				printf("End Point:  %e\n", dist);
				err |= inverse(spiral.centerPoint, point, &az12, &az21, &dist, eps);
				//printf("pt%i = movePointByCourse(cp, %f, 0, %f);\n", j, az12, dist);
				//printf("drawFlat(pt%i);\n", j);
				//printf("Azimuth/Distance:  %f %f\n", az12, dist);
				err |= spiralRadius(spiral, az12, &rad);
				//printf("Radius:  %f\n", rad);
				err |= invDist(point, arc.centerPoint, &dist, eps);
				//printf("Arc Dist:  %f\n", dist);
			}

			if (pointSet.length > 0) {
				//printf("Spiral:  %f %f %f %f\n", spiral.startRadius, spiral.endRadius, spiral.startAz, spiral.endAz);
				//printf("Start Point:  %f %f\n", spiral.startPoint.latitude / DEG2RAD, spiral.startPoint.longitude / DEG2RAD);
				//printf("End Point:  %f %f\n", spiral.endPoint.latitude / DEG2RAD, spiral.endPoint.longitude / DEG2RAD);
				//printf("Start/End:  %i %i\n", spcounter, epcounter);
				//printf("Distances:  %f %f\n", totDist, geo.length * pointSet.length);
				//printf("Point Set:  %i\n", pointSet.length);
				//printf("Failure\n\n");
			} else {
				//printf("No Points Returned\n\n");
			}
			//printf("\n");
		} else {
			successCount++;
		}

		//printf("\n");

		clearPtSet(&pointSet);

		//val = ptsAreSame(testPt, midChord, tol);
		//printf("%i\n", val);

	}

	printf("Spiral Spiral Intersect:  %i/1000\n", successCount);

}

void SpiralBoundary() {
	double latxExp, lonxExp, len, lineAz, startRad, endRad, az, rad, subAngle;
	LLPoint sp, ep, mc1, mc2, tempPt1, tempPt2, testPt;
	Geodesic line1, line2;
	Spiral spiral1, spiral2, mcSpiral, tempSp1, tempSp2;
	double az12, az21, dist;
	int i, j, val;
	long err = 0;

	Boundary b;

	LLPoint perpPt1, perpPt2;
	double perpCrs1, perpDist1, perpCrs2, perpDist2;

	Locus outLoc1, outLoc2;

	double eps = 1.0e-20;
	double tol = 1.37e-9;

	double distArray[7] = {0, .5, .99, .999999, 1.000001, 1.01, 1.5};

	srand(04012011);

	for (i=0;i<1;i++) {

		b = createBndry();

		latxExp = randLat();
		lonxExp = randLon();

		err |= createPt(&sp, latxExp*DEG2RAD, lonxExp*DEG2RAD);
		startRad = 30 + 0.01 * randDist() / 2; //0-54 nm
		endRad = 30 + 0.01 * randDist() / 2; //0-54 nm
		lineAz = randAzimuth();
		len = 120 + .01 * randDist();

		err |= direct(sp, lineAz * DEG2RAD, len, &ep, eps);
		err |= invCrs(sp, ep, &az12, &az21, eps);

		err |= createSpiral(&spiral1, sp, startRad, endRad, az12, az12, CLOCKWISE, eps);
		err |= createSpiral(&spiral2, ep, startRad, endRad, az21, az21, COUNTERCLOCKWISE, eps);

		err |= ptOnSpiral(spiral1, az12 + spiral1.dir * M_PI/4, &tempPt1, eps);
		err |= ptOnSpiral(spiral2, az21 + spiral2.dir * M_PI/4, &tempPt2, eps);
		err |= createGeo(&line1, tempPt1, tempPt2, INFINITE, eps);
		err |= ptOnSpiral(spiral1, az12 - spiral1.dir * M_PI/4, &tempPt1, eps);
		err |= ptOnSpiral(spiral2, az21 - spiral2.dir * M_PI/4, &tempPt2, eps);
		err |= createGeo(&line2, tempPt1, tempPt2, INFINITE, eps);

		err |= spiralRadius(spiral1, az12 + spiral1.dir * M_PI / 2, &rad);
		err |= createSpiralSection(spiral1, az12 + spiral1.dir * M_PI / 2, rad, &mcSpiral, eps);
		err |= spiralMidChord(mcSpiral, line1, &mc1, tol, eps);

		err |= spiralRadius(spiral2, az21 + spiral2.dir * M_PI / 2, &rad);
		err |= createSpiralSection(spiral2, az21 + spiral2.dir * M_PI / 2, rad, &mcSpiral, eps);
		err |= spiralMidChord(mcSpiral, line1, &mc2, tol, eps);

		err |= projectToGeo(sp, lineAz * DEG2RAD, mc1, &perpPt1, &perpCrs1, &perpDist1, tol, eps);
		err |= projectToGeo(sp, lineAz * DEG2RAD, mc2, &perpPt2, &perpCrs2, &perpDist2, tol, eps);

		err |= createLocus(&outLoc1, perpPt1, perpPt2, perpDist1, perpDist2, SEGMENT, tol, eps);

		err |= addLocusToBndry(&b, &outLoc1);

		err |= moveSpiralStartToPt(spiral1, mc1, &tempSp1, tol, eps);
		err |= moveSpiralStartToPt(spiral2, mc2, &tempSp2, tol, eps);

		err |= spiralRadius(spiral1, az12 - spiral1.dir * M_PI / 2, &rad);
		err |= createSpiralSection(spiral1, az12 - spiral1.dir * M_PI / 2, rad, &mcSpiral, eps);
		err |= spiralMidChord(mcSpiral, line2, &mc1, tol, eps);

		err |= spiralRadius(spiral2, az21 - spiral2.dir * M_PI / 2, &rad);
		err |= createSpiralSection(spiral2, az21 - spiral2.dir * M_PI / 2, rad, &mcSpiral, eps);
		err |= spiralMidChord(mcSpiral, line2, &mc2, tol, eps);

		err |= projectToGeo(sp, lineAz * DEG2RAD, mc1, &perpPt1, &perpCrs1, &perpDist1, tol, eps);
		err |= projectToGeo(sp, lineAz * DEG2RAD, mc2, &perpPt2, &perpCrs2, &perpDist2, tol, eps);
		err |= createLocus(&outLoc2, perpPt1, perpPt2, -perpDist1, -perpDist2, SEGMENT, tol, eps);

		err |= addLocusToBndry(&b, &outLoc2);

		err |= moveSpiralEndToPt(tempSp1, mc1, &spiral1, tol, eps);
		err |= moveSpiralEndToPt(tempSp2, mc2, &spiral2, tol, eps);

		err |= addSpiralToBndry(&b, &spiral1);
		err |= addSpiralToBndry(&b, &spiral2);

		//Create Points that are inside/outside
		err |= createGeo(&line1, outLoc1.geoStart, outLoc1.geoEnd, SEGMENT, eps);
		err |= createGeo(&line2, outLoc2.geoStart, outLoc2.geoEnd, SEGMENT, eps);

		for (j=0;j<7;j++) {
			dist = rand() % (int) floor(line1.length);
			err |= direct(line1.startPoint, line1.startAz, dist, &tempPt1, eps);
			err |= ptOnLocusFromGeoPt(outLoc1, tempPt1, &testPt, &az12, tol, eps);
			err |= direct(tempPt1, az12, distArray[j] * outLoc1.startDist, &testPt, eps);
			err |= projectToLocus(outLoc1, testPt, &tempPt1, &az12, &dist, tol, eps);
			val = ptIsInsideBndry(b, testPt, &err, tol, eps);
			printf("Val:  %i\n", val);
		}
		printf("\n");
		for (j=0;j<7;j++) {
			dist = rand() % (int) floor(line2.length);
			err |= direct(line2.startPoint, line2.startAz, dist, &tempPt1, eps);
			err |= ptOnLocusFromGeoPt(outLoc2, tempPt1, &testPt, &az12, tol, eps);
			err |= direct(tempPt1, az12, distArray[j] * fabs(outLoc2.startDist), &testPt, eps);
			err |= projectToLocus(outLoc2, testPt, &tempPt1, &az12, &dist, tol, eps);
			val = ptIsInsideBndry(b, testPt, &err, tol, eps);
			printf("Val:  %i\n", val);
		}
		printf("\n");
		for (j=0;j<7;j++) {
			subAngle = spiral1.subtendedAngle / DEG2RAD;
			az = spiral1.startAz + (spiral1.dir * (rand() % (int) floor(subAngle)) * DEG2RAD);
			err |= spiralRadius(spiral1, az, &rad);
			err |= direct(spiral1.centerPoint, az, distArray[j] * rad, &testPt, eps);
			err |= inverse(spiral1.centerPoint, testPt, &az12, &az21, &dist, eps);
			//printf("Az/Dist/Rad:  %f %f %f\n", az12, dist, rad);
			val = ptIsInsideBndry(b, testPt, &err, tol, eps);
			printf("Val:  %i\n", val);
		}
		printf("\n");
		for (j=0;j<7;j++) {
			subAngle = spiral2.subtendedAngle / DEG2RAD;
			az = spiral2.startAz + (spiral2.dir * (rand() % (int) floor(subAngle))) * DEG2RAD;
			err |= spiralRadius(spiral2, az, &rad);
			err |= direct(spiral2.centerPoint, az, distArray[j] * rad, &testPt, eps);
			err |= inverse(spiral2.centerPoint, testPt, &az12, &az21, &dist, eps);
			//printf("Az/Dist/Rad:  %f %f %f\n", az12, dist, rad);
			val = ptIsInsideBndry(b, testPt, &err, tol, eps);
			printf("Val:  %i\n", val);

		}
	}
}
} //namespace
