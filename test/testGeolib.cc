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

/*
 * Author: Balakrishna Babu
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

#ifndef FILEROOT
#define FILEROOT "."
#endif

using namespace geolib_idealab;

const double DEG2RAD = M_PI / 180.0;
const double RAD2DEG = 180.0 / M_PI;
const double NMTOL = 0.03 / 100.0 / 1852.0; //0.03 cm or ~1.62e-7 nm
const double LLDEGTOL = 5.0e-8; //5e-8 degrees or 0.00018 sec
const double AZDEGTOL = 5.0e-7; //5e-7 deg or 8.72665e-9 rad or 0.162 cm at 100 nm
const double NMSHORTDIST = SPHERICAL_AZIMUTH_CUTOFF_DISTANCE; //5e-4 nm or 0.926 m
const double NMLARGEDIST = 10000.0; //10000 nm




/*
 * Note - Cases 0 - 99 are reserved specifically for testing suites that are to be
 * 		included in the master testing suite for geolib.
 *
 * 		Any other testing suites/sets should be placed under the "individual testing"
 * 		section, starting with case number 100;
 *
 */

/*
 * Note - Cases 0 - 99 are reserved specifically for testing suites that are to be
 * 		included in the master testing suite for geolib.
 *
 * 		Any other testing suites/sets should be placed under the "individual testing"
 * 		section, starting with case number 100;
 *
 */

//The method below tests the maximum distance error between a locus and a
//geodesic with the same start and end points.

void testLocusDifferences() {

	FILE *fp;
	double locDist;
	double locLength; //10 feet
	int i,j,k, failCount = 0, testCount = 0;
	double testLat, testLon, testCourse, dist, dist2, projDist, geoTestCourse;
	LLPoint geoStart, geoEnd, testPt, projPt;
	Locus testLoc;

	LLPoint geoPt, testLocPt;

	double maxDist = 0, maxLocDist = 0, maxCourse, maxCourseDiff, courseDiff;

	ErrorSet err = 0;
	double tol = 1.37e-9, eps = 1e-20;

	srand(04012012);

	fp = fopen("C:/Users/jheidrich/Documents/locusErrorTestOutput.csv", "w");

	locLength = 200;
	projDist = 200 / 6076.12;

	locDist = 5;
	printf("Locus Length: %e\n", locLength);

	for (k=0;k<11;k=k+2) {
		locDist = (double)k;
		for (i=0;i<40;i++) {
			maxDist = 0;
			maxLocDist = 0;
			maxCourseDiff = 0;
			failCount = 0;
			for (j=0;j<1000;j++) {
				testLat = randLat() * M_PI / 180;
				testLon = randLon() * M_PI / 180;
				testCourse = randAzimuth() * M_PI / 180;
				err = createPt(&geoStart, testLat, testLon);
				if (!ptIsAtPole(geoStart, &err, tol, eps)) {
					testCount = testCount + 1;

					err = direct(geoStart, testCourse, locLength, &geoEnd, eps);
					err = createLocus(&testLoc, geoStart, geoEnd, locDist, locDist, LineType::SEGMENT, tol, eps);

					err = direct(geoStart, testCourse, projDist, &geoPt, eps);

					err = ptOnLocusFromGeoPt(testLoc, geoPt, &projPt, NULL, tol, eps);
					err = projectToGeo(geoStart, testCourse, projPt, &testPt, NULL, NULL, tol, eps);

					err = inverse(geoPt, testPt, &geoTestCourse, NULL, &dist, eps);

					err = minSubtendedAngle(geoTestCourse, testCourse, &courseDiff);
					courseDiff = fabs(courseDiff);
					if (courseDiff > M_PI / 2) {
						courseDiff = fabs(courseDiff - M_PI);
					}

					if ((projDist > 0) && (courseDiff > maxCourseDiff)) {
						maxCourseDiff = courseDiff;
						maxCourse = geoTestCourse;
					}

					if (dist > maxDist) {
						maxDist = dist;
					}

					err = ptOnLocusFromGeoPt(testLoc, testPt, &testLocPt, NULL, tol, eps);

					err = inverse(testLocPt, projPt, NULL, NULL, &dist2, eps);

//					fprintf(fp,"%f,%f,%f,%f,%f,%f,%f,%f\n",locDist,projDist,dist/tol,dist2/tol,testLat*180/M_PI,testLon*180/M_PI,testCourse, maxCourse);
					if (dist2 > tol) {
						failCount++;
						fprintf(fp,"%f,%f,%f,%f,%f,%f,%f,%f\n",locDist,projDist,dist/tol,dist2/tol,testLat*180/M_PI,testLon*180/M_PI,testCourse, maxCourseDiff);
//						fprintf(fp,"%f,%f,%f,%f,%f,%f,%f,%f\n",locDist,projDist,dist/tol,dist2/tol,testLat*180/M_PI,testLon*180/M_PI,testCourse, maxCourse);
					}

					if (dist > maxLocDist) {
						maxLocDist = dist;
					}
				}
			}
//			fprintf(fp,"%f,%f,%f,%f,%f,%f,%f,%f\n",locDist,projDist,dist/tol,dist2/tol,testLat*180/M_PI,testLon*180/M_PI,testCourse, maxCourseDiff);
			if (maxLocDist > 1.37e-9) {
//				fprintf(fp,"%f,%f,%f,%f,%f,%f,%f,%f\n",locDist,projDist,dist/tol,dist2/tol,testLat*180/M_PI,testLon*180/M_PI,testCourse, maxCourseDiff);
//				printf("Distance from geo start:  %e\n", projDist);
//				printf("Maximum Geo Distance Error: %e\n", maxDist);
//				printf("Maximum Locus Distance Error: %e\n", maxLocDist);
//				printf("Maximum Course Difference: %e\n", maxCourse);
			}
			printf("Fail Count:  %i\n", failCount);
			projDist = projDist - 5 / 6076.12;
			printf("Proj Dist: %f\n", projDist * 6076.12);
		}
		fprintf(fp,"\n");
		printf("\n\n");
		//locLength = locLength / 2;
		projDist = 200 / 6076.12;
		printf("Locus Length: %e\n", locDist);
	}
	fclose(fp);

	printf("Failed: %i\n", failCount);

//	fp = fopen("C:/Users/jheidrich/Documents/shortLocusTestOutput.csv", "w");
//
//	locLength = 0.00164579; //10 Feet
//
//	minLocDist = 0.0;
//	maxLocDist = 50.0;
//
//	for (k=2;k<11;k=k+2) {
//		locLength = 20.0;
//		locDist = (double)k;
//		for (j=1;j<20;j++) {
//			for (i=0;i<numTests;i++) {
//				testLat = randLat() * M_PI / 180;
//				testLon = randLon() * M_PI / 180;
//				testCourse = randAzimuth();
//
//				err = createPt(&geoStart, testLat, testLon);
//
//				if (!ptIsAtPole(geoStart, &err, tol, eps)) {
//					testCount = testCount + 1;
//
//					err = direct(geoStart, testCourse, locLength, &geoEnd, eps);
//
//					err = createLocus(&testLoc, geoStart, geoEnd, locDist, locDist, 0, tol, eps);
//
//					if (err) {
//						printf(formatErrorMessage(err));
//					}
//
//					err = createPt(&testGeoStart, testLoc.locusStart.latitude, testLoc.locusStart.longitude);
//					err = createPt(&testGeoEnd, testLoc.locusEnd.latitude, testLoc.locusEnd.longitude);
//
//					err = inverse(testGeoStart, testGeoEnd, &testGeoCourse, NULL, &dist, eps);
//
//					err = direct(testGeoStart, testGeoCourse, dist / 2, &testPt, eps);
//
//					err = projectToLocus(testLoc, testPt, &projPt, NULL, &projDist, tol, eps);
//
//					if ((projDist > 0.000164579)) {
//						failCount = failCount + 1;
//					}
//
//					if ((projDist > maxDist)) {
//						maxDist = projDist;
//						maxLat = geoStart.latitude;
//						maxLon = geoStart.longitude;
//						maxCourse = testCourse;
//					}
//				}
//			}
//			fprintf(fp, "%f,%f,%f,%f,%f,%f\n",locDist,maxDist*6076.12*30,locLength, maxLat,maxLon,maxCourse);
//			printf("Locus Length:  %f\n", locLength);
//			printf("Max Dist: %e\n", maxDist);
//			printf("Failures:  %i / %i\n", failCount, testCount);
//			maxDist = 0;
//			failCount = 0;
//			testCount = 0;
//			locLength = locLength + 3;
//		}
//	}
//	fclose(fp);
}

//void tempTest() {
//	LLPoint pts[5];
//	double lonList[5], latList[5];
//	int idx[5];
//	int i = 0;
//	ErrorSet err = 0;
//
//	err = createPt(&pts[4], .5, .2);
//	err = createPt(&pts[3], .5, .4);
//	err = createPt(&pts[2], .5, .2);
//	err = createPt(&pts[1], .5, .4);
//	err = createPt(&pts[0], .5, .8);
//
//	for (i=0;i<5;i++) {
//		lonList[i] = pts[i].longitude;
//		idx[i] = i;
//	}
//
//	for (i=0;i<5;i++) {
//		printf("Longitude: %f\n", lonList[i]);
//	}
//
//	sortPoints(latList,lonList,idx,5);
//
//	for (i=0;i<5;i++) {
//		printf("Longitude: %f\n", lonList[i]);
//	}

//	sortPoints(pts, idx, 5);
//
//	for (i=0;i<5;i++) {
//		printf("Longitude:  %f\n", pts[i].longitude);
//		printf("Index: %i\n", idx[i]);
//	}
//}

void set_fpu (unsigned int mode)
{
  asm ("fldcw %0" : : "m" (*&mode));
}

int main(int argc, char* argv[])
{

	int opt = -1;
    int i = 0;
    int saveSummary = 0;
    char* outputFileName = NULL;
    time_t start, stop;
    double diff;

    TestSuite masterSuite;
    TestSuite suite;
    TestSet set;

#ifdef DOUBLE
  set_fpu (0x27F);  /* use double-precision rounding */
#endif


    masterSuite = newTestSuite("testGeolib - master testing suite");

    printf("FileRoot is %s.\n", FILEROOT);

    if (argc > 1) // At least on option specified
    {
        /* Parse command line options */
        i = 1;
        while (i < argc)
        {
            if (!strcasecmp(argv[i], "-t") || !strcasecmp(argv[i], "--test"))
            {
                /* Test option given on command line */
                if ((i + 1 < argc) && (argv[i + 1][0] != '-'))
                {
                    i++;
                    opt = atoi(argv[i]);
                }
                else
                {
                    /* No option given, use default of zero */
                    opt = 0;
                }
            }
            else if (!strcasecmp(argv[i], "-f") || !strcasecmp(argv[i],
                    "--file"))
            {
                saveSummary = 1;
                /* Check for output file name */
                if ((i + 1 < argc) && (argv[i + 1][0] != '-'))
                {
                    /* -f followed by non-dash string
                     * string is assumed to be file name */
                    i++;
                    outputFileName = argv[i];
                }
                /* NB: If no file name specified, then output will be saved in default file */
            }
            i++;
        }
    }

    /* If option not set on command line, prompt user */
    if (opt < 0)
    {
        printf("Input Test option (Enter 0 to run all) => ");
        fflush(NULL);
        scanf("%d", &opt);
    }

    printf("option = %d\n", opt);
//	time(&start);
    	double time_c;
        clock_t clock_start=clock();

    switch (opt) {
    case 0: //Master Test Suite Case
        //*****KEEP THIS CASE EMPTY********//

        /***************************************************************************************
         ********************************** Testing Suites *************************************
         ***************************************************************************************/

        //TODO get new sphere test data
//    case 1:
//        suite = testSphereInverse_AllSets();
//        addTestSuite(suite, &masterSuite);
//        if (opt)
//            break;
//    case 2:
//        suite = testsphereInvDist_AllSets();
//        addTestSuite(suite, &masterSuite);
//        if (opt)
//            break;
//    case 3:
//        suite = testsphereInvCrs_AllSets();
//        addTestSuite(suite, &masterSuite);
//        if (opt)
//            break;
    case 4:
        suite = testInverse_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
    case 5:
        suite = testInvCrs_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
    case 6:
        suite = testInvDist_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
    case 7:
        suite = testDirect_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
    case 8:
        suite = testDirectLat_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
    case 9:
        suite = testDirectLon_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
    case 10:
        suite = testPtsAreSame_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
    case 11:
        suite = testArcTanToArcAndGeo_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
    case 12:
        suite = testArcFromStartAndEnd_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
    case 13:
        suite = testArcTanToCrs_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
    case 14:
        suite = testProjectToGeo_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
    case 15:
        suite = testArcTanToTwoGeos_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
    case 16:
        suite = testPtIsOnGeo_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
    case 17:
        suite = testGeoCrs_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
    case 18: //formerly case 175
        //TODO Create new test set
        suite = testProjectToGeoAtAngle_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
    case 21: //formerly case 18
        //TODO Recombine test data creation logic
        suite = testInitArcIntx_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
    case 22: //formerly case 19
        suite = testArcIntx_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
    case 23: 
        suite = testArcEndFromStartAndRadius_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
    case 24: 
        suite = testArcEndFromStartAndCenter_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
    case 25:
        suite = testGeoArcIntx_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
    case 28:
        suite = testPtsOnArcOnTanThruPt_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
    case 29:
        suite = testGetArcExtent_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
    case 30:
        suite = testArcLength_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
    case 31:
        suite = testProjectArcTanPtsToGeo_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
    case 32:
        suite = testPtIsOnArc_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
    case 33:
        suite = testArcFromLength_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
    case 36:
        suite = testPtIsInsideArc_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
    case 39:
        suite = testGeoTanToTwoCircles_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
    case 41: //formerly case 201
        suite = testCreateLocus_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
    case 42: //formerly case 202
        suite = testDistToLocusFromGeoDist_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
    case 43: //formerly case 203
        suite = testDistToLocusFromGeoPt_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
    case 44: //formerly case 204
        suite = testPtOnLocusFromGeoPt_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
    case 45: //formerly case 205
        suite = testPtIsOnLocus_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
    case 46: //formerly case 206
        suite = testLocusArcIntx_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
    case 47: //formerly case 207
        suite = testLocusGeoIntx_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
    case 48: //formerly case 208
        suite = testLocusIntx_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
    case 49: //formerly case 210
        suite = testLocusCrsAtPt_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
    case 50: //formerly case 211
        suite = testProjectToLocus_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
    case 51: //formerly case 214
        suite = testArcTanToTwoLoci_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
    case 52:
        suite = testPtIsOnCrs_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
    case 53:
        suite = testComputeSubtendedAngle_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
            /*
    case 60: 
        suite = testBndryCircleIntx_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
    case 61:
        suite = testBndryIntxExists_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
    case 63:
        suite = testArcsCoincide_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
    case 64:
    	suite = testGeoTanToArcAtAngleToGeo_AllSets();
    	addTestSuite(suite, &masterSuite);
    	if (opt)
    		break;
    case 65:
    	suite = testLociCoincide_AllSets();
    	addTestSuite(suite, &masterSuite);
    	if (opt)
    		break;
    case 66:
    	suite = testPtIsInsideBndry_AllSets();
    	addTestSuite(suite, &masterSuite);
    	if (opt)
    		break;
    case 67:
    	suite = testBndryGeoIntx_AllSets();
    	addTestSuite(suite, &masterSuite);
    	if (opt)
    		break;
    case 68:
    	suite = testBndryArcIntx_AllSets();
    	addTestSuite(suite, &masterSuite);
    	if (opt)
    		break;
    case 69:
    	suite = testBndryLocusIntx_AllSets();
    	addTestSuite(suite, &masterSuite);
    	if (opt)
    		break;
    case 70: 
        suite = testBndryCircleIntxExists_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
    case 71: 
        suite = testOrderBndry_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
    case 72:
        suite = testSeparateBndry_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
    case 73:
        suite = testSpiralGeoIntx_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
	case 74:
        suite = testSpiralLocusIntx_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
	case 75:
        suite = testGeoTanToSpiralAtAngleToGeo_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
	case 76:
        suite = testPtsOnSpiralOnTanThruPt_AllSets();
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
	case 77:
		suite = testGeoTanToTwoSpirals_AllSets();
		addTestSuite(suite, &masterSuite);
		if (opt)
			break;
	case 78:
		suite = testProjectToSpiral_AllSets();
		addTestSuite(suite, &masterSuite);
		if (opt)
			break;
            */
    case 96:
        suite = testMinSubtendedAngle_AllSets();//TODO determine why this suite is causing discrepancies in the master test suite if not at end of master suite - jamezcua
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
    case 97:
        suite = testCrsIntx_AllSets();//TODO determine why this suite is causing discrepancies in the master test suite if not at end of master suite - jamezcua
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;
    case 98:
        suite = testGeoIntx_AllSets();//TODO determine why this suite is causing discrepancies in the master test suite if not at end of master suite - jamezcua
        addTestSuite(suite, &masterSuite);
        if (opt)
            break;

    case 99: //Master Test Suite Case
        //****Do not update this case****//
        break;
        /***************************************************************************************
         ********************* Individual Test Sets/Suites *************************************
         ***************************************************************************************/
        /*
    case 105:
    	set = testPtIsInsideBndry_Set2();
    	break;
    case 106:
    	set = testPtsAreInsideBndry_Set3();
    	break;
    case 107:
    	suite = testPtsAreInsideBndry_AllSets();
    	addTestSuite(suite, &masterSuite);
    	break;
    case 108: //formerly case 209
        set = testLocusIntx_Set2();
        break;
    case 116:
    	set = testGeoArcIntx_Set1();
    	break;
    case 117:
    	set = testGeoArcIntx_Set2();
    	break;
//    case 119:
//    	set = testInsertionSort_Set1();
//    	break;
    case 120:
    	set = testComputeSubtendedAngle_Set1();
    	break;
//    case 121:
//    	set = testSortPtsByAz_Set1();
//    	break;
    case 122:
    	set = testArcsCoincide_Set1();
    	break;
    case 124:
    	set = testGeoTanToArcAtAngleToGeo_Set1();
    	break;
    case 127:
    	set = testPtIsInsideBndry_Set1();
    	break;
    case 128:
    	set = testOrderBndry_Set1();
    	break;
    case 130:
    	set = testFindSetMaxAndMin_Set1();
    	break;
    case 134:
    	set = testAddPointToLLPointSet_Set1();
    	break;
    case 135:
    	set = testSpiralMidChord_Set1();
    	break;
    case 136:
    	set = testSpiralGeoIntx_Set1();
    	break;
    case 137:
    	set = testSpiralLocusIntx_Set1();
    	break;
    case 138:
    	set = testGeoTanToSpiralAtAngleToGeo_Set1();
    	break;
    case 139:
    	set = testPtsOnSpiralOnTanThruPt_Set1();
    	break;
    case 140:
    	set = testGeoTanToTwoSpirals_Set1();
    	break;
    case 141:
    	set = testProjectToSpiral_Set1();
    	break;
    case 142:
    	testSpiralArcIntx_Set1();
    	break;
    case 143:
        testSpiralIntx_Set1();
        break;
        */
    case 144:
    	testLocusDifferences();
    	break;
    case 145:
    	testDirectInverseConsistency();
    	break;
    case 146:
    	testDirectInverseMathematicaData();
    	break;
    default:
        printf("*** Invalid Test option, Exiting ***\n");
        break;
    }

//    time(&stop);
//    diff = difftime(stop,start);

    time_c=(double)(clock()- clock_start)/CLOCKS_PER_SEC;

    printf("\n");

    if (saveSummary)
    {
        outputTestSuiteMetrics(masterSuite, outputFileName);
    }
    else if (!opt)
    {
        reportConstants(outputFileName);  /* NULL arg will default to STDOUT */
        displayTestSuite(masterSuite);
        //Pause the console output. No need to set a breakpoint to see test results
        // This breaks autobuilds.  --RKIRKMAN
        //newflush(stdin);
    }

    printf("Time Elapsed: %lf sec\n",time_c);

    return 0;

}

