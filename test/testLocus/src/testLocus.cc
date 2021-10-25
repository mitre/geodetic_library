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


/*
 * NAME: testCreateLocus_Set1
 *
 * DESCRIPTION:
 * 		This function is used to test the createLocus function.
 *
 * 		This function runs the test data created by Rich Snow
 *
 * 		Approach for testing createLocus is as follows:
 *
 * 		1. Create input data. Input data consists of start and end
 *         points of the defining geodesic, start and end distances that define
 * 		   the locus, and the lineType.
 * 		2. For verification, compute the distance from start of geodesic to start
 * 	       of locus, compute the distance from end of geodesic to end of locus,
 * 		   compute the angle from start of locus to end of geodesic at the
 *  	   start of geodesic, and compute the angle from end of locus to start of
 *    	   geodesic at the end of geodesic.
 * 		3. The distances should match the input distances and the angles should be
 *   	   90 degrees.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testCreateLocus_Set1(TestSet) - A test set with the folling metrics:
 *  	testCases(int): The number of test cases executed.
 * 		pass(int): The number of test cases that passed.
 * 		fail(int): The number of test cases that failed.
 * 		setupFailures(int): The number of tests cases that generated an error prior
 * 					to actually executing the given test case.
 * 		errors(int): The number of test cases that generated an error.  Test cases that
 * 					generated an error prior to executing the given test are included in
 * 					this count.  Test cases that generated an error while executing the
 * 					given test are included as well.
 * 		unverified(int): The number of invalid records encountered while processing
 * 					the set of test data.
 * 		testingError(int): Return 1 if unable to execute testing function.  Return 0 otherwise.
 *
 */
TestSet testCreateLocus_Set1()
{
    double DEG2RAD = M_PI / 180.0;
    int jj;
    double geoLatS, geoLonS, startDist, endDist;

    LLPoint geoStart, geoEnd, locStartExp, locEndExp;
    LineType lineType;  //0=SEGMENT,1=SEMIINFINITE,2=INFINITE
    Locus locus;
    Locus* locusPtr;
    double geoSLocSDist, geoELocEDist;
    double crsGSLS, crsLSGS, crsGELE, crsLEGE;
    double geoSCrsDiff, geoECrsDiff;
    double geolen, geoAz, dir1, dir2, tempAz, tempCrs;
    int passedCount=0, failedCount=0;
    int errorCount = 0;
    int unverifiedCount = 0;
    int setupFailureCount = 0;
    int testCaseCount = 0;
    ErrorSet err=0;

    TestSet set;

    printf("Start testCreateLocus_Set1\n");

    set = newTestSet("testCreateLocus_Set1");

    srand(04132011);

    for (jj = 0; jj < 100; jj++)
    {
    	err = 0;
        geoLatS = randLat();
        geoLonS = randLon();
        geoStart.latitude = geoLatS * DEG2RAD;
        geoStart.longitude = geoLonS * DEG2RAD;
        geoAz = DEG2RAD * randAzimuth();
        geolen = 100.0 + 0.01 * randDist();
        err |= direct(geoStart, geoAz, geolen, &geoEnd, EPS);
        err |= invCrs(geoEnd, geoStart, &tempAz, &tempCrs, EPS);
        if (jj == 0)
          geoEnd = geoStart;
        dir1 = 1.0;
        if ((rand() % 2) == 0)
          dir1 = -1.0;
        dir2= 1.0;
        if ((rand() % 2) == 0)
          dir2 = -1.0;
        startDist = dir1 * 0.02 * randDist();
        endDist = dir2 * 0.02 * randDist();
        locusPtr = &locus;
        lineType = SEGMENT;
        if (jj == 98)
         lineType = SEMIINFINITE;
        else if (jj == 99)
         lineType = INFINITE;
        if (startDist >= 0.0)
          direct(geoStart,geoAz + M_PI_2, startDist, &locStartExp, EPS); 
        else
          direct(geoStart,geoAz - M_PI_2, -startDist, &locStartExp, EPS); 
        if (endDist >= 0.0)
          direct(geoEnd,tempAz - M_PI_2, endDist, &locEndExp, EPS);
        else
          direct(geoEnd, tempAz + M_PI_2, -endDist, &locEndExp, EPS);
        //printf("jj = %d geoAz %4.8f geolen %4.5f sdist %4.5f edist %4.5f\n",jj,geoAz,geolen,startDist,endDist);

        err |= createLocus(locusPtr, geoStart, geoEnd, startDist, endDist, lineType, TOL, EPS);

        if ( getMaskedError(err, getMaskAll()) )
        {
                if ((jj == 0) && (err == INVALID_SHAPE_ERR))
                {
                  passedCount++;
        	  testCaseCount++;
        	  continue;
                }
                else
                {
        	  printf("Error occurred in createLocus err=0x%lx\n",err);
        	  failedCount++;
        	  errorCount++;
        	  testCaseCount++;
        	  continue;
                }
        }

        err |= inverse(geoStart, locus.locusStart, &crsGSLS, &crsLSGS, &geoSLocSDist, EPS);
        err |= inverse(geoEnd, locus.locusEnd, &crsGELE, &crsLEGE, &geoELocEDist, EPS);

        if ( getMaskedError(err, getMaskAll()) )
        {
        	printf("Error occurred in post-createLocus err=0x%lx\n",err);
        	failedCount++;
        	errorCount++;
        	testCaseCount++;
        	continue;
        }

        geoSCrsDiff = locus.geoAz - crsGSLS;
        if (fabs(geoSCrsDiff) < M_PI)
        {
          if (geoSCrsDiff > 0.0)
            geoSLocSDist = -geoSLocSDist;
        } 
        else
        {
          tempAz = locus.geoAz - M_PI;
          tempCrs = crsGSLS - M_PI;
          if (tempAz < 0.0)
            tempAz = tempAz + M_2PI;
          if (tempCrs < 0.0)
            tempCrs = tempCrs + M_2PI;
          if (tempAz - tempCrs > 0.0)
            geoSLocSDist = -geoSLocSDist;
        }
        geoECrsDiff = locus.geoRevAz - crsGELE;
        if (fabs(geoECrsDiff) < M_PI)
        { 
          if (geoECrsDiff < 0.0)
            geoELocEDist = -geoELocEDist;
        }
        else
        {
          tempAz = locus.geoRevAz - M_PI;
          tempCrs = crsGELE - M_PI;
          if (tempAz < 0.0)
            tempAz = tempAz + M_2PI;
          if (tempCrs < 0.0)
            tempCrs = tempCrs + M_2PI;
          if (tempAz - tempCrs < 0.0)
            geoELocEDist = -geoELocEDist;
        }

        if ( (fabs(startDist-geoSLocSDist)<TESTTOL) &&
        	 (fabs(endDist-geoELocEDist)<TESTTOL) &&
               ptsAreSame(geoStart, locus.geoStart, TESTTOL) &&
               ptsAreSame(geoEnd, locus.geoEnd, TESTTOL) &&
               (fabs(locus.geoLength - geolen) < TESTTOL) && 
                 (locus.lineType == lineType) &&
                  ptsAreSame(locus.locusStart, locStartExp, TESTTOL) &&
                  ptsAreSame(locus.locusEnd, locEndExp, TESTTOL)) {

        	passedCount++;
        }
        else {

                printf("test jj = %d failed\n",jj);
        	failedCount++;
        }

        testCaseCount++;
    }

    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("Finish testCreateLocus_Set1\n\n\n");

    return set;

}

/*
 * NAME: testCreateLocus_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the createLocus function.
 *
 * 		This function runs all the test cases from set1 and returns
 * 		the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testCreateLocus_AllSets(TestSuite) - A test suite with the folling metrics:
 *  	testCases(int): The cumulative number of test cases executed.
 * 		pass(int): The cumulative number of test cases that passed.
 * 		fail(int): The cumulative number of test cases that failed.
 * 		setupFailures(int): The cumulative number of tests cases that generated an error prior
 * 					to actually executing the given test case.
 * 		errors(int): The cumulative number of test cases that generated an error.  Test cases that
 * 					generated an error prior to executing the given test are included in
 * 					this count.  Test cases that generated an error while executing the
 * 					given test are included as well.
 * 		unverified(int): The cumulative number of invalid records encountered while processing
 * 					the set of test data.
 * 		testingSuiteError(int): Return 1 if unable to execute one or more testing functions.
 * 					Return 0 otherwise.
 *
 */
TestSuite testCreateLocus_AllSets()
{

	TestSuite suite;
	TestSet set1;

    printf("\nStart testCreateLocus_AllSets\n");

    suite = newTestSuite("testCreateLocus_AllSets");

    set1 = testCreateLocus_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("Finish testCreateLocus_AllSets\n\n\n");

    return suite;
}

/*
 * NAME: testDistToLocusFromGeoDist_Set1
 *
 * DESCRIPTION:
 * 		This function is used to test the distToLocusFromGeoDist function.
 *
 * 		This function runs the test data created by Rich Snow
 *
 *      Approach for testing distToLocusFromGeoDist is as follows:
 *
 *      1. Create input data. Input data consists of locus (start and
 *         end points of the defining geodesic, start and end distances that define
 *         the locus, and the lineType) and distance.
 *      2. Create NewLocus using input data.
 *      3. Verify the result by comparing to the linear interpolation.
 *
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testDistToLocusFromGeoDist_Set1(TestSet) - A test set with the folling metrics:
 *  	testCases(int): The number of test cases executed.
 * 		pass(int): The number of test cases that passed.
 * 		fail(int): The number of test cases that failed.
 * 		setupFailures(int): The number of tests cases that generated an error prior
 * 					to actually executing the given test case.
 * 		errors(int): The number of test cases that generated an error.  Test cases that
 * 					generated an error prior to executing the given test are included in
 * 					this count.  Test cases that generated an error while executing the
 * 					given test are included as well.
 * 		unverified(int): The number of invalid records encountered while processing
 * 					the set of test data.
 * 		testingError(int): Return 1 if unable to execute testing function.  Return 0 otherwise.
 *
 */
TestSet testDistToLocusFromGeoDist_Set1()
{
	//Locus loc, double distance);
    double DEG2RAD = M_PI / 180.0;
    int jj;
    double geoLatS, geoLonS, startDist, endDist;
    double distance;
    double geolen, geoAz, dir1, dir2;
    LineType lineType;  //0=SEGMENT,1=SEMIINFINITE,2=INFINITE

    LLPoint geoStart, geoEnd;
    Locus locus;
    Locus* locusPtr;
    double distToLoc, distToLocExp;
    int passedCount=0, failedCount=0;
    int testCaseCount = 0, setupFailureCount = 0, errorCount = 0, unverifiedCount = 0;
    ErrorSet err=0;

    TestSet set;

    printf("Start testDistToLocusFromGeoDist_Set1\n");

    set = newTestSet("testDistToLocusFromGeoDist_Set1");

    srand(05152011);

    for (jj = 0; jj < 100; jj++)
    {
    	err = 0;
        geoLatS = randLat();
        geoLonS = randLon();
        geoStart.latitude = geoLatS * DEG2RAD;
        geoStart.longitude = geoLonS * DEG2RAD;
        geoAz = DEG2RAD * randAzimuth();
        geolen = 100.0 + 0.01 * randDist();
        err |= direct(geoStart, geoAz, geolen, &geoEnd, EPS);
        dir1 = 1.0;
        if ((rand() % 2) == 0)
          dir1 = -1.0;
        dir2= 1.0;
        if ((rand() % 2) == 0)
          dir2 = -1.0;
        startDist = dir1 * 0.02 * randDist();
        endDist = dir2 * 0.02 * randDist();
        distance = (double)(rand() % (int)geolen) + 0.1;
        if (jj == 1)
          distance = geolen + 10.0;
        locusPtr = &locus;
        lineType = SEGMENT;
        if (jj == 98)
         lineType = SEMIINFINITE;
        else if (jj == 99)
         lineType = INFINITE;
        //printf("jj = %d geolen %4.5f dist %4.5f sdist %4.5f edist %4.5f\n",jj,geolen,distance,startDist,endDist);
        

        err |= createLocus(locusPtr, geoStart, geoEnd, startDist, endDist, lineType, TOL, EPS);

        if ( getMaskedError(err, getMaskAll()) )
        {
        	printf("Error occurred in createLocus prior to testing distToLocusFromGeoDist err=0x%lx\n",err);
        	setupFailureCount++;
            testCaseCount++;
            errorCount++;
        	continue;
        }

        if (jj == 0)
          locus.geoLength = 0.0;
        if (locus.geoLength > 0.0) {
            distToLocExp = startDist + distance * (endDist - startDist) / geolen;
        }
        else {
            distToLocExp = 0.0; 
        }

        distToLoc = distToLocusFromGeoDist(*locusPtr,distance);

        if (fabs(distToLocExp-distToLoc) < TESTTOL) {
        	passedCount++;
        }
        else {
        	failedCount++;
        }
        testCaseCount++;
    }

    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("Finish testDistToLocusFromGeoDist_Set1\n");

    return set;

}

/*
 * NAME: testDistToLocusFromGeoDist_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the distToLocusFromGeoDist function.
 *
 * 		This function runs all the test cases from all the test sets and returns
 * 		the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testDistToLocusFromGeoDist_AllSets(TestSuite) - A test suite with the folling metrics:
 *  	testCases(int): The cumulative number of test cases executed.
 * 		pass(int): The cumulative number of test cases that passed.
 * 		fail(int): The cumulative number of test cases that failed.
 * 		setupFailures(int): The cumulative number of tests cases that generated an error prior
 * 					to actually executing the given test case.
 * 		errors(int): The cumulative number of test cases that generated an error.  Test cases that
 * 					generated an error prior to executing the given test are included in
 * 					this count.  Test cases that generated an error while executing the
 * 					given test are included as well.
 * 		unverified(int): The cumulative number of invalid records encountered while processing
 * 					the set of test data.
 * 		testingSuiteError(int): Return 1 if unable to execute one or more testing functions.
 * 					Return 0 otherwise.
 *
 */
TestSuite testDistToLocusFromGeoDist_AllSets()
{

	TestSuite suite;
	TestSet set1;

    printf("\nStart testDistToLocusFromGeoDist_AllSets\n");

    suite = newTestSuite("testDistToLocusFromGeoDist_AllSets");

    set1 = testDistToLocusFromGeoDist_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("Finish testDistToLocusFromGeoDist_AllSets\n\n\n");

    return suite;
}

/*
 * NAME: testDistToLocusFromGeoPt_Set1
 *
 * DESCRIPTION:
 * 		This function is used to test the distToLocusFromGeoPt function.
 *
 * 		This function runs the test data created by Rich Snow
 *
 *      Approach for testing distToLocusFromGeoPt is as follows:
 *
 *      1. Create input data. Input data consists of locus (start and
 *         end points of the defining geodesic, start and end distances that define
 *         the locus, and the lineType) and a point (latitude, longitude).
 *      2. Create NewLocus using input data.
 *      3. Verify the result by comparing to the linear interpolation.
 *
 *
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testDistToLocusFromGeoPt(TestSet) - A test set with the folling metrics:
 *  	testCases(int): The number of test cases executed.
 * 		pass(int): The number of test cases that passed.
 * 		fail(int): The number of test cases that failed.
 * 		setupFailures(int): The number of tests cases that generated an error prior
 * 					to actually executing the given test case.
 * 		errors(int): The number of test cases that generated an error.  Test cases that
 * 					generated an error prior to executing the given test are included in
 * 					this count.  Test cases that generated an error while executing the
 * 					given test are included as well.
 * 		unverified(int): The number of invalid records encountered while processing
 * 					the set of test data.
 * 		testingError(int): Return 1 if unable to execute testing function.  Return 0 otherwise.
 *
 */
TestSet testDistToLocusFromGeoPt_Set1()
{

    double DEG2RAD = M_PI / 180.0;
    int jj;
    double geoLatS, geoLonS, startDist, endDist;
    double distance;
    double geolen, geoAz, dir1, dir2, distTogeoPt, crsTogeoPt;
    LineType lineType;  //0=SEGMENT,1=SEMIINFINITE,2=INFINITE

    LLPoint geoStart, geoEnd, geoPt, geoPtProj;
    Locus locus;
    Locus* locusPtr;
    double faz, distToLoc, distToLocExp;
    double distToStart, distToEnd;
    double crsFromPoint, distFromPoint;
    int passedCount=0, failedCount=0;
    int testCaseCount = 0, setupFailureCount = 0, errorCount = 0, unverifiedCount = 0;
    ErrorSet err=0;

    TestSet set;

    printf("Start testDistToLocusFromGeoPt_Set1\n");

    set = newTestSet("testDistToLocusFromGeoPt_Set1");

    srand(05202011);

    for (jj = 0; jj < 100; jj++)
    {
    	err = 0;
        geoLatS = randLat();
        geoLonS = randLon();
        geoStart.latitude = geoLatS * DEG2RAD;
        geoStart.longitude = geoLonS * DEG2RAD;
        geoAz = DEG2RAD * randAzimuth();
        geolen = 100.0 + 0.01 * randDist();
        err |= direct(geoStart, geoAz, geolen, &geoEnd, EPS);
        dir1 = 1.0;
        if ((rand() % 2) == 0)
          dir1 = -1.0;
        dir2= 1.0;
        if ((rand() % 2) == 0)
          dir2 = -1.0;
        startDist = dir1 * 0.02 * randDist();
        endDist = dir2 * 0.02 * randDist();
        locusPtr = &locus;
        distTogeoPt = (double)(rand() % (int)geolen) + 0.1;
        if (jj == 0)
          distTogeoPt = 6000.0; 
        else if (jj == 1)
          distTogeoPt = 0.5 * TOL;
        else if (jj == 2)
          distTogeoPt = 1.1 * TOL;
        else if (jj == 3)
        {
          startDist = 5500.0;
          endDist = 5500.0;
        }
        if (jj <= 25)
          crsTogeoPt = geoAz;
        else if (jj > 25 && jj <= 50)
          crsTogeoPt = modcrs(geoAz + M_PI);
        else if (jj > 50)
          crsTogeoPt = DEG2RAD * randAzimuth();
        err |= direct(geoStart, crsTogeoPt, distTogeoPt, &geoPt, EPS);
        lineType = SEGMENT;
        if (jj == 98)
         lineType = SEMIINFINITE;
        else if (jj == 99)
         lineType = INFINITE;
        //printf("jj = %d distTogeoPt %4.12f crsTogeoPt = %4.6f startDist %4.5f endDist %4.5f\n",jj,distTogeoPt,crsTogeoPt,startDist,endDist);

        err |= createLocus(locusPtr, geoStart, geoEnd, startDist, endDist, lineType, TOL, EPS);

        err |= projectToGeo(geoStart, geoAz, geoPt, &geoPtProj,
                                  &crsFromPoint, &distFromPoint, TOL, EPS);
        err |= invDist(geoStart, geoPtProj, &distToStart, EPS);
        err |= invDist(geoEnd, geoPtProj, &distToEnd, EPS);
        if ( (distToStart<distToEnd) &&
        	 (fabs(distToEnd-distToStart-geolen)<TOL) ) {
        	distance = -distToStart;
        }
        else {
        	distance = distToStart;
        }

        distToLocExp = startDist + distance * (endDist - startDist) / geolen;

        if ( getMaskedError(err, getMaskAll()) )
        {
        	printf("Error occurred in pre-distToLocusFromGeoPt err=0x%lx\n",err);
        	testCaseCount++;
            errorCount++;
            setupFailureCount++;
        	continue;
        }

        distToLoc = distToLocusFromGeoPt(locus,geoPt,&faz,&err,TOL,EPS);

        if ( getMaskedError(err, getMaskAll()) )
        {
          if (((jj == 0) || (jj == 3)) && (err == RADIUS_OUT_OF_RANGE_ERR))
          {
            passedCount++; 
            testCaseCount++;
            continue;
          }
          else
          {
        	printf("Error occurred in distToLocusFromGeoPt err=0x%lx\n",err);
        	failedCount++;
            errorCount++;
            testCaseCount++;
        	continue;
          }
        }

        if (fabs(distToLocExp -distToLoc) < TESTTOL) {
        	passedCount++;
        }
        else {
        	failedCount++;
        }
        testCaseCount++;
    }

    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("Finish testDistToLocusFromGeoPt_Set1\n");

    return set;

}

/*
 * NAME: testDistToLocusFromGeoPt_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the distToLocusFromGeoPt function.
 *
 * 		This function runs all the test cases from all the test sets and returns
 * 		the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testDistToLocusFromGeoPt_AllSets(TestSuite) - A test suite with the folling metrics:
 *  	testCases(int): The cumulative number of test cases executed.
 * 		pass(int): The cumulative number of test cases that passed.
 * 		fail(int): The cumulative number of test cases that failed.
 * 		setupFailures(int): The cumulative number of tests cases that generated an error prior
 * 					to actually executing the given test case.
 * 		errors(int): The cumulative number of test cases that generated an error.  Test cases that
 * 					generated an error prior to executing the given test are included in
 * 					this count.  Test cases that generated an error while executing the
 * 					given test are included as well.
 * 		unverified(int): The cumulative number of invalid records encountered while processing
 * 					the set of test data.
 * 		testingSuiteError(int): Return 1 if unable to execute one or more testing functions.
 * 					Return 0 otherwise.
 *
 */
TestSuite testDistToLocusFromGeoPt_AllSets()
{

	TestSuite suite;
	TestSet set1;

    suite = newTestSuite("testDistToLocusFromGeoPt_AllSets");

    printf("\nStart %s\n", suite.name.c_str());

    set1 = testDistToLocusFromGeoPt_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("Finish %s\n\n\n", suite.name.c_str());

    return suite;
}

/*
 * NAME: testPtOnLocusFromGeoPt_Set1
 *
 * DESCRIPTION:
 * 		This function is used to test the ptOnLocusFromGeoPt function.
 *
 * 		This function runs the test data created by Rich Snow
 *
 * Approach for testing ptOnLocusFromGeoPt is as follows:
 *
 *      1. Create input data. Input data consists of locus (start and
 *         end points of the defining geodesic, start and end distances that define
 *         the locus, and the lineType) and a point (latitude, longitude).
 *      2. Create NewLocus using input data.
 *      3. Calculate the expected point on the Locus.
 *      4. Call the function under test using the created locus and other input data to compute the point.
 *      5. Verify that the expected and computed points are the same.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testPtOnLocusFromGeoPt_Set1(TestSet) - A test set with the folling metrics:
 *  	testCases(int): The number of test cases executed.
 * 		pass(int): The number of test cases that passed.
 * 		fail(int): The number of test cases that failed.
 * 		setupFailures(int): The number of tests cases that generated an error prior
 * 					to actually executing the given test case.
 * 		errors(int): The number of test cases that generated an error.  Test cases that
 * 					generated an error prior to executing the given test are included in
 * 					this count.  Test cases that generated an error while executing the
 * 					given test are included as well.
 * 		unverified(int): The number of invalid records encountered while processing
 * 					the set of test data.
 * 		testingError(int): Return 1 if unable to execute testing function.  Return 0 otherwise.
 *
 */
TestSet testPtOnLocusFromGeoPt_Set1()
{

    double DEG2RAD = M_PI / 180.0;
    int jj;
    double geoLatS, geoLonS, startDist, endDist;
    double geolen, geoAz, dir1, dir2, distTogeoPt, fcrs, bcrs;
    LineType lineType;  //0=SEGMENT,1=SEMIINFINITE,2=INFINITE

    LLPoint geoStart, geoEnd, geoPt, newGeoStart;
    LLPoint ptonloc, ptonlocExp;
    Locus locus;
    Locus* locusPtr;
    double perpCrs;
    double crsToPtonloc, distToPtonloc;
    int passedCount=0, failedCount=0;
    int testCaseCount = 0, errorCount = 0, setupFailureCount = 0, unverifiedCount = 0;
    ErrorSet err=0;

    TestSet set;

    printf("\nStart testPtOnLocusFromGeoPt_Set1\n");

    set = newTestSet("testPtOnLocusFromGeoPt_Set1");

    srand(05232011);

    for(jj = 0; jj < 100; jj++)
    {
    	err = 0;
        geoLatS = randLat();
        geoLonS = randLon();
        geoStart.latitude = geoLatS * DEG2RAD;
        geoStart.longitude = geoLonS * DEG2RAD;
        geoAz = DEG2RAD * randAzimuth();
        geolen = 100.0 + 0.01 * randDist();
        err |= direct(geoStart, geoAz, geolen, &geoEnd, EPS);
        dir1 = 1.0;
        if ((rand() % 2) == 0)
          dir1 = -1.0;
        dir2= 1.0;
        if ((rand() % 2) == 0)
          dir2 = -1.0;
        startDist = dir1 * 0.02 * randDist();
        endDist = dir2 * 0.02 * randDist();
        locusPtr = &locus;
        distTogeoPt = 0.01 * randDist();
        lineType = SEGMENT;
        if (jj == 98)
         lineType = SEMIINFINITE;
        else if (jj == 99)
         lineType = INFINITE;
        if ((rand() % 2) == 0)
          distTogeoPt = -distTogeoPt;
        if (jj == 0)
          distTogeoPt = 6000.0;
        else if (jj == 1)
          distTogeoPt = 0.0001; // SMALL_DIST_THRESHOLD is 0.00054
        else if (jj == 2)
          distTogeoPt = -0.0001;
        else if (jj == 3)
        {
          startDist = 5500.0;
          endDist = 5500.0;
        }
        //printf("jj = %d geolen %4.5f distTogeoPt %4.5f startDist %4.5f endDist %4.5f\n",jj,geolen,distTogeoPt,startDist,endDist);

        err |= createLocus(locusPtr, geoStart, geoEnd, startDist, endDist, lineType, TOL, EPS);

        err |= direct(geoStart, geoAz, distTogeoPt, &geoPt, EPS);

        distToPtonloc = startDist + distTogeoPt * (endDist - startDist) / geolen;
        err |= invCrs(geoPt, geoStart, &fcrs, &bcrs, EPS);

        if (fabs(distTogeoPt) <= SMALL_DIST_THRESHOLD)
        {
          err |= direct(geoEnd, locusPtr->geoRevAz, geolen + 100.0 / 1852.0, &newGeoStart, EPS); 
          err |= invCrs(geoPt, newGeoStart, &fcrs, NULL, EPS);
          distTogeoPt = fabs(distTogeoPt);
        }

        if (distTogeoPt > 0.0)
        {
        if (distToPtonloc > 0.0)
          crsToPtonloc = modcrs(fcrs + M_PI + M_PI_2);
        else
          crsToPtonloc = modcrs(fcrs + M_PI - M_PI_2);
        }

        if (distTogeoPt < 0.0)
        {
        if (distToPtonloc > 0.0)
          crsToPtonloc = modcrs(fcrs + M_PI_2);
        else
          crsToPtonloc = modcrs(fcrs - M_PI_2);
        }

        err |= direct(geoPt,crsToPtonloc,fabs(distToPtonloc),&ptonlocExp, EPS);

        if ( getMaskedError(err, getMaskAll()) )
        {
        	printf("Error occurred in pre-ptOnLocusFromGeoPt err=0x%lx\n",err);
        	testCaseCount++;
            setupFailureCount++;
            errorCount++;
        	continue;
        }

        err |= ptOnLocusFromGeoPt(locus,geoPt,&ptonloc,&perpCrs,TOL,EPS);

        if ( getMaskedError(err, getMaskAll()) )
        {
           if (((jj == 0) || (jj == 3)) && (err == RADIUS_OUT_OF_RANGE_ERR))
           {
             passedCount++;
             testCaseCount++;
             continue;
           }
           else
           {
             printf("Error occurred in ptOnLocusFromGeoPt err=0x%lx\n",err);
             failedCount++;
             errorCount++;
             testCaseCount++;
             continue;
           }
        }

        if (ptsAreSame(ptonlocExp, ptonloc, TESTTOL)) {
        	passedCount++;
        }
        else {
        	err |= invDist(ptonloc, ptonlocExp, &distToPtonloc, EPS);
        		printf("Distance between points: %e\n", distToPtonloc);
                printf("jj = %d failed\n",jj);
        	failedCount++;
        }
        testCaseCount++;
    }

    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("\nFinish testPtOnLocusFromGeoPt_Set1\n");

    return set;

}

/*
 * NAME: testPtOnLocusFromGeoPt_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the ptOnLocusFromGeoPt function.
 *
 * 		This function runs all the test cases from all the test sets and returns
 * 		the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testDistToLocusFromGeoPt_AllSets(TestSuite) - A test suite with the folling metrics:
 *  	testCases(int): The cumulative number of test cases executed.
 * 		pass(int): The cumulative number of test cases that passed.
 * 		fail(int): The cumulative number of test cases that failed.
 * 		setupFailures(int): The cumulative number of tests cases that generated an error prior
 * 					to actually executing the given test case.
 * 		errors(int): The cumulative number of test cases that generated an error.  Test cases that
 * 					generated an error prior to executing the given test are included in
 * 					this count.  Test cases that generated an error while executing the
 * 					given test are included as well.
 * 		unverified(int): The cumulative number of invalid records encountered while processing
 * 					the set of test data.
 * 		testingSuiteError(int): Return 1 if unable to execute one or more testing functions.
 * 					Return 0 otherwise.
 *
 */
TestSuite testPtOnLocusFromGeoPt_AllSets()
{

	TestSuite suite;
	TestSet set1;

    printf("\nStart testPtOnLocusFromGeoPt_AllSets\n");

    suite = newTestSuite("testPtOnLocusFromGeoPt_AllSets");

    set1 = testPtOnLocusFromGeoPt_Set1();
    addTestSet(set1,&suite);


    displayTestSuite(suite);

    printf("Finish testPtOnLocusFromGeoPt_AllSets\n\n\n");

    return suite;
}

/*
 * NAME: testPtIsOnLocus_Set1
 *
 * DESCRIPTION:
 * 		This function is used to test the ptIsOnLocus function.
 *
 * 		This function runs the test data created by Balakrishna Babu and uses
 * 		his original test plan to execute the actual tests.
 *
 *	    Approach for testing ptIsOnLocus is as follows:
 *
 *	    1)	Use random number to generate input data for locus.  Input data consists of
 *	        locus (start and end points of the defining geodesic, start and end distances
 *	        that define the locus, and the lineType).  Select latitude and longitude for
 *	        geoStart, and select azimuth and distance to compute geoEnd.
 *	    2)	Create NewLocus using input data.
 *	    3)	Select test points using the following three steps (4,5,6).
 *	    4)	Select a distance (dAlong, nm values of -1000, -5, 0, 5, 1000, 0.5*geoLength,
 *	        geoLength) and compute geoBase along the defining geodesic.
 *	    5)	Using the geoBase, compute the point on the locus (ptOnLocusFromGeoPt).
 *	        Save ptonloc and perpCrs.
 *	    6)	Compute test points perpendicular to geodesic (geoBase, perpCrs) with distances of
 *	        -1000, -10, 0, 10, 1000, distToLocus nm.  The last test point is the same as ptonloc.
 *	    7)	Call the function under test with testPoints.
 *	    8)	All testPoints should return false, except for ptonloc which should return true
 *	        depending on the lineType.  onlineExp array is setup with expected result of 0/1.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testPtIsOnLocus_Set1(TestSet) - A test set with the folling metrics:
 *  	testCases(int): The number of test cases executed.
 * 		pass(int): The number of test cases that passed.
 * 		fail(int): The number of test cases that failed.
 * 		setupFailures(int): The number of tests cases that generated an error prior
 * 					to actually executing the given test case.
 * 		errors(int): The number of test cases that generated an error.  Test cases that
 * 					generated an error prior to executing the given test are included in
 * 					this count.  Test cases that generated an error while executing the
 * 					given test are included as well.
 * 		unverified(int): The number of invalid records encountered while processing
 * 					the set of test data.
 * 		testingError(int): Return 1 if unable to execute testing function.  Return 0 otherwise.
 *
 */
TestSet testPtIsOnLocus_Set1()
{
    double DEG2RAD = M_PI / 180.0;
    int testNum=0;
    char testname[256];
    double geoLatS, geoLonS, geoLatE, geoLonE, startDist, endDist;
    double geoAz, geoLength;
    double dAlong[7] = {-1000.0, -5.0, 0.0, 5.0, 1000.0, 0.0, 0.0};
    double dCross[6] = {-1000.0, -5.0, 0.0, 5.0, 1000.0, 0.0};
    LineType lineType;  //0=SEGMENT,1=SEMIINFINITE,2=INFINITE

    LLPoint geoStart, geoEnd, geoBase, testPoint;
    LLPoint ptonloc, pt2;
    Locus locus;
    double perpCrs;
    //[lineType][dAlong][dCross]
    int onlineExp[3][7][6] =
								{	{	{0,0,0,0,0,0},  {0,0,0,0,0,0},   {0,0,0,0,0,1},
										{0,0,0,0,0,1},  {0,0,0,0,0,1},   {0,0,0,0,0,1},   {0,0,0,0,0,1}	},
									{	{0,0,0,0,0,0},  {0,0,0,0,0,0},   {0,0,0,0,0,1},
										{0,0,0,0,0,1},  {0,0,0,0,0,1},   {0,0,0,0,0,1},   {0,0,0,0,0,1}	},
									{	{0,0,0,0,0,1},  {0,0,0,0,0,1},   {0,0,0,0,0,1},
										{0,0,0,0,0,1},  {0,0,0,0,0,1},   {0,0,0,0,0,1},   {0,0,0,0,0,1}	}	};
    int i, j, online;
    int passedCount=0, failedCount=0;
    int testCaseCount = 0, setupFailureCount = 0, errorCount = 0, unverifiedCount = 0;
	int onlineExp035, onlineExp045;
    ErrorSet err=0;

    TestSet set;

    long newSeed=20080322;

    printf("Start testPtIsOnLocus_Set1\n");

    set = newTestSet("testPtIsOnLocus_Set1");

	srand(newSeed);  //Initialize the random number generator
	
	onlineExp035 = onlineExp[0][3][5];
	onlineExp045 = onlineExp[0][4][5];

	while (testNum<100)
    {
		//Initialize or reset to defaults
    	err = 0;
    	perpCrs = 0.0;
	    onlineExp[0][3][5] = onlineExp035; 
	    onlineExp[0][4][5] = onlineExp045; 

		//Select start point, azimuth and distance to end point for the defining geodesic
    	//Select start and end distances for defining locus start and end points
		geoLatS = randLat();
		geoLonS = randLon();
		geoAz = DEG2RAD * randAzimuth();
    	geoLength = 0.1 * randDist();  //about 540 nm
    	dAlong[5] = 0.5 * geoLength;
    	dAlong[6] = geoLength;
    	startDist = 0.001 * randDist();  //about 5.40 nm
    	endDist = 0.001 * randDist();  //about 5.40 nm
    	testNum++;

        sprintf(testname,"TEST%-d",testNum);

        geoStart.latitude = geoLatS * DEG2RAD;
        geoStart.longitude = geoLonS * DEG2RAD;

        //Compute end point
        err |= direct(geoStart, geoAz, geoLength, &geoEnd, EPS);

        geoLatE = geoEnd.latitude / DEG2RAD;
        geoLonE = geoEnd.longitude / DEG2RAD;
        lineType = SEGMENT;

        if (dAlong[3] > geoLength) { //lineType=SEGMENT
        	onlineExp[0][3][5] = 0;
        }
        if (dAlong[4] > geoLength) { //lineType=SEGMENT
        	onlineExp[0][4][5] = 0;
        }

        for (int lineTypeI=0; lineTypeI < 3; lineTypeI++) {
            lineType = static_cast<LineType>(lineTypeI);
            err |= createLocus(&locus, geoStart, geoEnd, startDist, endDist, lineType, TOL, EPS);

	        //For each point along geodesic
	        for (i=0; i<7; i++) {

	        	dCross[5] = distToLocusFromGeoDist (locus, dAlong[i]);
	        	err |= direct(geoStart,geoAz,dAlong[i],&geoBase,EPS);
	        	err |= ptOnLocusFromGeoPt(locus,geoBase,&ptonloc,&perpCrs,TOL,EPS);

	            //perpCrs is always from base geodesic to locus
	            //All other dCross values have corresponding value with opposite sign
	            if (dCross[5] < 0.0) {
	            	perpCrs += M_PI;
	            }

		        //For each point along perpendicular course
		        for (j=0; j<6; j++) {

		        	err |= direct(geoBase,perpCrs,dCross[j],&testPoint,EPS);

		            if ( getMaskedError(err, getMaskAll()) )
		            {
			        	printf("Error occurred in pre-ptIsOnLocus err=0x%lx lineTyep=%d i=%d j=%d\n",
			        			err,lineType,i,j);
			        	testCaseCount++;
                        setupFailureCount++;
                        errorCount++;
			        	continue;
			        }

		            online = ptIsOnLocus(locus,testPoint,&pt2,&err,TOL*5,EPS);

		            if ( getMaskedError(err, getMaskAll()) )
		            {
			        	printf("Error occurred in ptIsOnLocus err=0x%lx lineTyep=%d i=%d j=%d\n",
			        			err,lineType,i,j);
			        	failedCount++;
                        errorCount++;
                        testCaseCount++;
			        	continue;
			        }

			        if ( online == onlineExp[lineType][i][j] ) {
			        	if (PRINT_PASSED_CASES) printLLPoint("   geoStart:  ",geoStart);
			        	if (PRINT_PASSED_CASES) printLLPoint("   geoEnd:    ",geoEnd);
			        	if (PRINT_PASSED_CASES) printf("%s geoAz=%14.8f geoLength=%14.8f startDist=%14.8f endDist=%14.8f\n",
			            		testname, geoAz/DEG2RAD, geoLength, startDist, endDist);
			        	if (PRINT_PASSED_CASES) printLLPoint("   locusStart:",locus.locusStart);
			        	if (PRINT_PASSED_CASES) printLLPoint("   locusEnd:  ",locus.locusEnd);
			        	if (PRINT_PASSED_CASES) printLLPoint("   geoBase:   ",geoBase);
			        	if (PRINT_PASSED_CASES) printLLPoint("   ptonloc:   ",ptonloc);
			        	if (PRINT_PASSED_CASES) printf("%s perpCrs=%14.8f\n", testname, perpCrs/DEG2RAD);
			        	if (PRINT_PASSED_CASES) printLLPoint("   testPoint: ",testPoint);

			        	if (PRINT_PASSED_CASES) printf("%s passed lineType=%d  dAlong=%14.8f dCross=%14.8f  onlineExp=%d online=%d\n",
			        			testname,lineType,dAlong[i],dCross[j],onlineExp[lineType][i][j],online);
			        	passedCount++;
			        }
			        else {
			            printLLPoint("   geoStart:  ",geoStart);
			            printLLPoint("   geoEnd:    ",geoEnd);
			            printf("%s geoAz=%14.8f geoLength=%14.8f startDist=%14.8f endDist=%14.8f\n",
			            		testname, geoAz/DEG2RAD, geoLength, startDist, endDist);
			            printLLPoint("   locusStart:",locus.locusStart);
			            printLLPoint("   locusEnd:  ",locus.locusEnd);
			            printLLPoint("   geoBase:   ",geoBase);
			            printLLPoint("   ptonloc:   ",ptonloc);
			            printf("%s perpCrs=%14.8f\n", testname, perpCrs/DEG2RAD);
			            printLLPoint("   testPoint: ",testPoint);

			        	printf("%s failed lineType=%d  dAlong=%14.8f dCross=%14.8f  onlineExp=%d online=%d\n",
			        			testname,lineType,dAlong[i],dCross[j],onlineExp[lineType][i][j],online);
			        	failedCount++;
			        }
			        testCaseCount++;
		        }
	        }
        }
    }
    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("Finish testPtIsOnLocus_Set1\n");

    return set;

}

/*
 * NAME: testPtIsOnLocus_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the ptIsOnLocus function.
 *
 * 		This function runs all the test cases from all the test sets and returns
 * 		the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testPtIsOnLocus_AllSets(TestSuite) - A test suite with the folling metrics:
 *  	testCases(int): The cumulative number of test cases executed.
 * 		pass(int): The cumulative number of test cases that passed.
 * 		fail(int): The cumulative number of test cases that failed.
 * 		setupFailures(int): The cumulative number of tests cases that generated an error prior
 * 					to actually executing the given test case.
 * 		errors(int): The cumulative number of test cases that generated an error.  Test cases that
 * 					generated an error prior to executing the given test are included in
 * 					this count.  Test cases that generated an error while executing the
 * 					given test are included as well.
 * 		unverified(int): The cumulative number of invalid records encountered while processing
 * 					the set of test data.
 * 		testingSuiteError(int): Return 1 if unable to execute one or more testing functions.
 * 					Return 0 otherwise.
 *
 */
TestSuite testPtIsOnLocus_AllSets()
{

	TestSuite suite;
	TestSet set1;

    printf("\nStart testPtIsOnLocus_AllSets\n");

    suite = newTestSuite("testPtIsOnLocus_AllSets");

    set1 = testPtIsOnLocus_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("Finish testPtIsOnLocus_AllSets\n\n\n");

    return suite;
}

/*
 * NAME: testLocusArcIntx_Set1
 *
 * DESCRIPTION:
 * 		This function is used to test the locusArcIntx function.
 *
 * 		This function runs the test data created by Balakrishna Babu and uses
 * 		his original test plan to execute the actual tests.
 *
 * 		The approach for testing geoArcIntx is as follows:
 *
 *      1)	Use random number to generate input data for Arc (center C, radius r1).
 *      2)	Select two azimuths and find those points on the Arc.
 *          These will be the expected intersection points (X1, X2).
 *      3)	Build the locus that goes through X1, X2
 * 	    	a.	Select two distances to extend CX1 and CX2 to define the start and
 * 	    		end points of the defining geodesic (GS, GE).
 * 	    	b.	Project X1 and X2 on to the geodesic (GS, GE).
 * 	    		Save the points XP1, XP2 and the distances h1, h2.
 * 	    	c.	Compute the distances along the geodesic of points XP1 and XP2 from GS (d1, d2).
 * 	    	d.	Compute slope for locus (h2-h1)/(d2-d1).
 * 	    	e.	Compute start and end distances for the locus (hs, he) using h2, d2, h1, d1, slope.
 * 	    	f.	Create NewLocus using GS, GE, hs, he.
 *      4)	Call the function under test with Arc and Locus.
 *      5)	The result should contain two intersection points which should match with X1 and X2.
 *      6)	Also check that the points are on Arc and Locus.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testLocusArcIntx_Set1(TestSet) - A test set with the folling metrics:
 *  	testCases(int): The number of test cases executed.
 * 		pass(int): The number of test cases that passed.
 * 		fail(int): The number of test cases that failed.
 * 		setupFailures(int): The number of tests cases that generated an error prior
 * 					to actually executing the given test case.
 * 		errors(int): The number of test cases that generated an error.  Test cases that
 * 					generated an error prior to executing the given test are included in
 * 					this count.  Test cases that generated an error while executing the
 * 					given test are included as well.
 * 		unverified(int): The number of invalid records encountered while processing
 * 					the set of test data.
 * 		testingError(int): Return 1 if unable to execute testing function.  Return 0 otherwise.
 *
 */
TestSet testLocusArcIntx_Set1()
{

    int testNum=0;
    char testname[256];
    double DEG2RAD = M_PI / 180.0;

    double latC, lonC;
    double r1, az1, az2;
    double crsCX1, crsCX2;
    LLPoint center, X1, X2;

    double distCgeoStart1, distCgeoEnd1;
    double crsX1XP1, crsX2XP2, crsGeo12, crsGeo21, crsGeo1X1, crsGeo2X2;
    double geoLength, d1, d2, h1, h2, hs, he, slope1;

    LLPoint XP1, XP2;
    LLPoint geoStart1, geoEnd1;
    LineType lineType1=INFINITE;  //0=SEGMENT,1=SEMIINFINITE,2=INFINITE
    Locus locus1;
    Locus* locusPtr1 = &locus1;

    LLPointPair intx;
    int nX;

    double fcrs1, bcrs1, fcrs2, bcrs2;
    int passedCount=0, failedCount=0;
    int testCaseCount = 0, errorCount = 0, setupFailureCount = 0;
    ErrorSet err=0;

    TestSet set;

    long newSeed=20080523;

    printf("Start testLocusArcIntx_Set1\n");

    set = newTestSet("testLocusArcIntx_Set1");


    srand(newSeed);  //Initialize the random number generator

    while (testNum<500)
    {
    	err = 0;

    	//generate input data for Arc (center C, radius r1)
        latC = randLat();
        lonC = randLon();
        r1 = 0.1 * randDist();  //about 540 nm
        az1 = randAzimuth();
        az2 = randAzimuth();
        testNum++;

        sprintf(testname,"TEST%-d",testNum);

        //Select two azimuths and find expected intersection points (X1, X2) on the Arc
        center.latitude = latC * DEG2RAD;
        center.longitude = lonC * DEG2RAD;
        crsCX1 = az1 * DEG2RAD;
        crsCX2 = az2 * DEG2RAD;

        err |= direct(center, crsCX1, r1, &X1, EPS);
        err |= direct(center, crsCX2, r1, &X2, EPS);

        //Select two distances to extend CX1 and CX2 to define the start and end
        //points of the defining geodesic (GS, GE)
        distCgeoStart1 = 0.1 * randDist();  //about 540 nm
        distCgeoEnd1 = 0.1 * randDist();  //about 540 nm

        err |= direct(center, crsCX1, distCgeoStart1, &geoStart1, EPS);
        err |= direct(center, crsCX2, distCgeoEnd1, &geoEnd1, EPS);
        err |= inverse(geoStart1, geoEnd1, &crsGeo12, &crsGeo21, &geoLength, EPS);

        //Project X1 and X2 on to the geodesic (GS, GE).
        //Save the points XP1, XP2 and the distances h1, h2.
        err |= projectToGeo(geoStart1, crsGeo12, X1,
                                &XP1, &crsX1XP1, &h1, TOL, EPS);
        err |= projectToGeo(geoStart1, crsGeo12, X2,
                                &XP2, &crsX2XP2, &h2, TOL, EPS);

        //Compute the distances along the geodesic of points XP1 and XP2 from GS (d1, d2).
        err |= inverse(geoStart1, XP1, &fcrs1, &bcrs1, &d1, EPS);
        err |= inverse(geoStart1, XP2, &fcrs2, &bcrs2, &d2, EPS);

        //Correct the sign of h1, h2, d1 and d2 based on the courses crsGeo12, fcrs1, fcrs2...
          
        err |= inverse(geoStart1, X1, &crsGeo1X1, &bcrs1, NULL, EPS);
        if (fabs(crsGeo12 - crsGeo1X1) > M_PI)
        {
            if (crsGeo1X1 > crsGeo12)
              h1 = -h1;
        }
        else if (crsGeo12 > crsGeo1X1)
          h1 = -h1;

        err |= inverse(geoEnd1, X2, &crsGeo2X2, &bcrs1, NULL, EPS);
        if (fabs(crsGeo2X2 - crsGeo21) > M_PI)
        {
          if (crsGeo21 > crsGeo2X2)
            h2 = -h2;
        }  
        else if ( crsGeo2X2 > crsGeo21) 
          h2 = -h2;

        if (fabs(crsGeo12 - fcrs1) > 0.1)
          d1 = -d1;
        if (fabs(crsGeo12 - fcrs2) > 0.1)
          d2 = -d2;

        //Compute slope for locus (h2-h1)/(d2-d1)
        slope1 = (h2-h1)/(d2-d1);

        //Compute start and end distances for the locus
        //h = h1 + (d-d1)*slope
        hs = h1 + (-d1)*slope1;
        he = h1 + (geoLength-d1)*slope1;
        if (fabs(hs) > 1000.0 || fabs(he) > 1000.0)
        {
          testNum--;
          continue;
        } 

        //Create locus
        err |= createLocus(locusPtr1,geoStart1,geoEnd1,hs,he,lineType1, TOL, EPS);

        if ( getMaskedError(err, getMaskAll()) )
        {
            printf("\n--------------------------------\n");
            printf("\n%s\n", testname);
            printf("Error occurred in pre-locusArcIntx err=0x%lx\n",err);
            testCaseCount++;
            errorCount++;
            setupFailureCount++;
            printf("\n--------------------------------\n");
            continue;
        }

        //call the function under test with locus and Arc
        nX = 0;
        intx[0].latitude = 0.0;
        intx[0].longitude = 0.0;
        intx[1].latitude = 0.0;
        intx[1].longitude = 0.0;

        err |= locusArcIntx(locus1, center, r1, intx, &nX, TOL, EPS);

        if ( getMaskedError(err, getMaskAll()) )
        {
            printf("\n--------------------------------\n");
            printf("\n%s\n", testname);
            printf("Error occurred in locusArcIntx err=0x%lx\n",err);
            printf("\nLocus:\n");
            displayBoundaryLocus(locus1, 1);
            printf("\nArc Center:\n");
            displayBoundaryLLPoint(center,1);
            printf("\nArc Radius:%.20lf [nm]", r1);
            printf("\ntol=%le",TOL);
            printf("\neps=%le",EPS);
            failedCount++;
            errorCount++;
            testCaseCount++;
            printf("\n--------------------------------\n");
            continue;
        }

        //Verify results against expected intersection points X1 and X2
        if ( (nX == 2) &&
        	 ( ptsAreSame(intx[0], X1, TESTTOL) || ptsAreSame(intx[1], X1, TESTTOL) ) &&
        	 ( ptsAreSame(intx[0], X2, TESTTOL) || ptsAreSame(intx[1], X2, TESTTOL) ) )
        {
        	if (PRINT_PASSED_CASES) printf("\n--------------------------------\n");
        	if (PRINT_PASSED_CASES) printf("%s center: %14.8f %14.8f  r1=%14.8f  az1=%14.8f  az2=%14.8f\n",
                   testname, latC, lonC, r1, az1, az2);
        	if (PRINT_PASSED_CASES) printLocus("  locus1",locus1);
        	if (PRINT_PASSED_CASES) printLLPoint("    X1",X1);
        	if (PRINT_PASSED_CASES) printLLPoint("    X2",X2);
        	if (PRINT_PASSED_CASES) printLLPoint("    intx[0]",intx[0]);
        	if (PRINT_PASSED_CASES) printLLPoint("    intx[1]",intx[1]);

        	if (PRINT_PASSED_CASES) printf("%s passed nX: %d\n",
                   testname,nX);
        	if (PRINT_PASSED_CASES) printf("--------------------------------\n");
            passedCount++;
        }
        else if ( (nX == 1) &&
           	 ( ptsAreSame(intx[0], X1, TESTTOL) || ptsAreSame(intx[0], X2, TESTTOL) ) )
        {
        	if (PRINT_PASSED_CASES) printf("\n--------------------------------\n");
        	if (PRINT_PASSED_CASES) printf("%s center: %14.8f %14.8f  r1=%14.8f  az1=%14.8f  az2=%14.8f\n",
                   testname, latC, lonC, r1, az1, az2);
        	if (PRINT_PASSED_CASES) printLocus("  locus1",locus1);
        	if (PRINT_PASSED_CASES) printLLPoint("    X1",X1);
        	if (PRINT_PASSED_CASES) printLLPoint("    X2",X2);
        	if (PRINT_PASSED_CASES) printLLPoint("    intx[0]",intx[0]);
        	if (PRINT_PASSED_CASES) printLLPoint("    intx[1]",intx[1]);

        	if (PRINT_PASSED_CASES) printf("%s passed nX: %d\n",
                   testname,nX);
        	if (PRINT_PASSED_CASES) printf("--------------------------------\n");
            passedCount++;
        }
        else
        {
        	printf("\n--------------------------------\n");
            printf("%s center: %14.8f %14.8f  r1=%14.8f  az1=%14.8f  az2=%14.8f\n",
                   testname, latC, lonC, r1, az1, az2);
            printf("h1 = %14.8f h2 = %14.8f d1 = %14.8f d2 = %14.8f crsGeo12 = %14.8f crsGeo1X1 = %14.8f\n",h1,h2,d1,d2,crsGeo12,crsGeo1X1);
            printLocus("  locus1",locus1);
            printLLPoint("    X1",X1);
            printLLPoint("    X2",X2);
            printLLPoint("    intx[0]",intx[0]);
            printLLPoint("    intx[1]",intx[1]);
            //displayMatlabLocus(locus1,"loc",0);
            //displayMatlabPt(X1,"X1",0);
            //displayMatlabPt(X2,"X2",0);
            //displayMatlabPt(intx[0],"intx0",0);
            //displayMatlabPt(intx[1],"intx1",0);

            printf("%s failed nX: %d\n",
                   testname,nX);
            printf("--------------------------------\n");
            failedCount++;
        }
        testCaseCount++;

//        printf("--------------------------------\n");
    }

    //Bug 27723
    center.latitude = 34.986577628206916*DEG2RAD;
    center.longitude = -80.85502297678627*DEG2RAD;
    XP1.latitude = 34.97895881015197*DEG2RAD;
    XP1.longitude = -80.98694475787796*DEG2RAD;
    XP2.latitude = 34.88742878906488*DEG2RAD;
    XP2.longitude = -80.8004936507796*DEG2RAD; 
    err |= invDist(center, XP1, &r1, EPS);
    X1.latitude = 35.11308239677041*DEG2RAD;
    X1.longitude = -80.95777365783344*DEG2RAD;
    X2.latitude = 35.05093608760776*DEG2RAD;
    X2.longitude = -80.95244301459037*DEG2RAD;
    err |= createLocus(&locus1,X1,X2,2.0,2.0,SEMIINFINITE,TOL,EPS);
    err |= locusArcIntx(locus1,center,r1,intx,&nX,TOL,EPS);
    err |= invCrs(center, XP1, &fcrs1, &bcrs1, EPS);
    err |= invCrs(center, XP2, &fcrs2, &bcrs2, EPS);
    if((nX == 1) && ptIsOnArc(center,r1,fcrs1,fcrs2,CLOCKWISE,intx[0],&err, TESTTOL,EPS))
      passedCount++;
    else
      failedCount++;
    testCaseCount++;

    //Bug 24707
    center.latitude = 36.0149712221257 * DEG2RAD;
    center.longitude = -79.0539309886121 * DEG2RAD;
    XP1.latitude = 36.07955064184488*DEG2RAD;
    XP1.longitude = -78.97379261207271*DEG2RAD;
    XP2.latitude = 36.0676030438308*DEG2RAD;
    XP2.longitude = -78.96149558916392*DEG2RAD;
    err |= invDist(center, XP1, &r1, EPS);
    X1.latitude = 35.94661733190422*DEG2RAD;
    X1.longitude = -78.86958692402968*DEG2RAD;
    X2.latitude = 35.990936634934776*DEG2RAD;
    X2.longitude = -78.92354877580625*DEG2RAD;
    err |= createLocus(&locus1,X1,X2,2.0,2.0,SEMIINFINITE,TOL,EPS);
    err |= locusArcIntx(locus1,center,r1,intx,&nX,TOL,EPS);
    Arc ob;
    err |= createArc(&ob, center, XP1, XP2, CLOCKWISE, TOL, EPS);    
    err |= invCrs(center, XP1, &fcrs1, &bcrs1, EPS);
    err |= invCrs(center, XP2, &fcrs2, &bcrs2, EPS);
    if((nX == 2) && ptIsOnArc(center,r1,fcrs1,fcrs2,CLOCKWISE,intx[0],&err, TESTTOL,EPS) && ptIsOnArc(center,r1,fcrs1,fcrs2,CLOCKWISE,intx[1],&err, TESTTOL,EPS))
      passedCount++;
    else
      failedCount++;
    testCaseCount++;

    //Bug 24825
    X1.latitude = -64.77375985979448*DEG2RAD;
    X1.longitude = 64.23598361181546*DEG2RAD;
    X2.latitude = -64.84446933617082*DEG2RAD;
    X2.longitude = 64.26809854346578*DEG2RAD;
    center.latitude = -64.83856852756524*DEG2RAD;
    center.longitude = 64.34164744230075*DEG2RAD;
    err |= createLocus(&locus1,X1,X2,-0.8382411059745324,-2.0,SEGMENT,TOL,EPS);
    r1 = 0.08228941684665227;
    err |= locusArcIntx(locus1,center,r1,intx,&nX,TOL,EPS);
    err |= direct(center, 0.0, r1, &XP1, EPS);
    err |= createArc(&ob, center, XP1, XP1, CLOCKWISE, TOL, EPS);    
    if (nX == 0)
      passedCount++;
    else
      failedCount++;
    testCaseCount++;
    //displayMatlabLocus(locus1, "locus1", 0);
    //displayMatlabArc(ob, "ob", 0);

    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = 0;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("Finish testLocusArcIntx_Set1\n");

    return set;

}

/*
 * NAME: testLocusArcIntx_Set2
 *
 * DESCRIPTION:
 * 		This function is used to test the locusArcIntx function.
 *
 * 		This function runs the test cases created by Juan Amezcua.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testLocusArcIntx_Set2(TestSet) - A test set with the following metrics:
 *  	testCases(int): The number of test cases executed.
 * 		pass(int): The number of test cases that passed.
 * 		fail(int): The number of test cases that failed.
 * 		setupFailures(int): The number of tests cases that generated an error prior
 * 					to actually executing the given test case.
 * 		errors(int): The number of test cases that generated an error.  Test cases that
 * 					generated an error prior to executing the given test are included in
 * 					this count.  Test cases that generated an error while executing the
 * 					given test are included as well.
 * 		unverified(int): The number of invalid records encountered while processing
 * 					the set of test data.
 * 		testingError(int): Return 1 if unable to execute testing function.  Return 0 otherwise.
 *
 */
TestSet testLocusArcIntx_Set2()
{

    char testname[256];
    int testCase = 0;
    double DEG2RAD = M_PI / 180.0;

    LLPoint* center;
    LLPoint* arcStart;
    LLPoint* arcEnd;
    LLPoint* tempPt1;
    LLPoint* tempPt2;
    LLPoint* geoMidPt;
    LLPoint* originCenter;
    LLPoint* originArcStart;
    LLPoint* originArcEnd;
    Arc* arc;
    Arc* originArc;
    double radius;
    double crsGStLSt, crsLStGst, crsGEnLEn, crsLEnGen;
    double crsGeoMidLocMid, distGeoMidLocMid, delta;
    double d1, d2, d3, crs31, dist13, crs32, dist23;
    double tempCrs1, tempDist, tempCrs2, crsTangent;
    double slopeLocAngle, thetaRad;
    int maxGeoLength = 50; //nm
    int locusRightOrient; //0 = left, 1 = right
    ArcDirection dir;
    int numType0 = 0, maxType0 = 50, thetaDeg, caseType;
    int maxLocusStartEndDist = 50; // nm
    int maxHorizTranslationDist = 10; //nm

    LineType lineType=INFINITE;
    double geoStartLat, geoStartLon;
    double geoLen, geoCrs;
    double startDist, endDist;
    LLPoint* geoStart;
    LLPoint* geoEnd;
    Locus* locus;

    LLPointPair intx;
    int nX,expNX;
    LLPoint* X1;
    LLPoint* X2;


    int passedCount=0, failedCount=0;
    int testCaseCount = 0, errorCount = 0, setupFailureCount = 0, unverifiedCount = 0;
    ErrorSet err=0;
    long seed = 20080523;
    int passed = 0;

    TestSet set;
    srand(seed);

    printf("Start testLocusArcIntx_Set2\n");

    set = newTestSet("testLocusArcIntx_Set2");

    while(testCase < 1000)
    {
	err = 0;
	passed = 0;

    	//generate random locus geodesic
        geoStartLat = randLat()*DEG2RAD; //rad
        geoStartLon = randLon()*DEG2RAD; //rad
        geoLen = (rand() % maxGeoLength) + 1;  //nm  (the +1 is keep it positive due to modular arithmetic)
        geoCrs = randAzimuth()*DEG2RAD; //rad

        geoStart = NULL;
        geoStart = new LLPoint();
        err |= createPt(geoStart, geoStartLat, geoStartLon);

        geoEnd = NULL;
        geoEnd = new LLPoint();
        err |= direct(*geoStart, geoCrs, geoLen, geoEnd, EPS);

        //generate random locus

        //pick locus orientation (left = 0/right = 1)
        locusRightOrient = rand() % 2;

        dir = ArcDirection::CLOCKWISE;
        if (locusRightOrient)
          dir = ArcDirection::COUNTERCLOCKWISE;
        //pick random start/end distance
        startDist = (rand() % maxLocusStartEndDist) + 1; //nm  (the +1 is keep it positive due to modular arithmetic)
        if ( !locusRightOrient ) startDist *= -1; //locus is left orientation so sign is negative
        endDist = startDist;  //nm

        locus = new Locus();
        err |= createLocus(locus,*geoStart,*geoEnd,startDist,endDist,lineType, TOL, EPS);

        /*
         * choose a random azimuth theta in the range (0,90] degrees
         * (use of degrees in this comment and throughout this function is for explanatory purposes only)
         * we work strictly in radians except when using modular arithmetic as necessary
         *
         * Theta directly determines the properties of the test case arc
         * where theta = 90 degrees will describe an arc centered on the locus and
         * as theta goes to zero the arc center will move away from the locus in the
         * direction of the locus' geodesic
         */
        //in degrees when using modular arithmetic and then converting to radians
        thetaDeg = ((rand() % 90) + 1);//keep positive in range (0,90] degrees
        thetaRad = thetaDeg*DEG2RAD;

        if (thetaDeg != 90){
  	    //find temp courses for use in finding infinite geodesic that bisects the locus
	    err |= invCrs(locus->geoStart, locus->locusStart, &crsGStLSt, &crsLStGst, EPS);
	    err |= invCrs(locus->geoEnd, locus->locusEnd, &crsGEnLEn, &crsLEnGen, EPS);

	    if (locusRightOrient)
	    {
		crsLStGst = modcrs(crsLStGst + thetaRad);
		crsLEnGen = modcrs(crsLEnGen - thetaRad);
	    } else {
		crsLStGst = modcrs(crsLStGst - thetaRad);
		crsLEnGen = modcrs(crsLEnGen + thetaRad);
	    }

		//find temp points for use in defining the infinite geodesics for use in finding the arc center point
		tempPt1 = NULL;
		tempPt1 = new LLPoint();
		tempPt2 = NULL;
		tempPt2 = new LLPoint();
		err |= direct(locus->locusStart, crsLStGst, 2*maxGeoLength, tempPt1, EPS);
		err |= direct(locus->locusEnd, crsLEnGen, 2*maxGeoLength, tempPt2, EPS);

		//find arc center that intersects locus start/end pts, this will be the intersection of the infinite geodesics
		originCenter = NULL;
		originCenter = new LLPoint();
		err |= geoIntx(locus->locusStart, *tempPt1, SEGMENT, &crs31, &dist13, locus->locusEnd, *tempPt2, SEGMENT, &crs32, &dist23, originCenter, TOL, EPS);
        } else {
          //tangency so find locus midpoint

          //find geodesic mid point (tempPt1)
          delta = locus->geoLength / 2;

          geoMidPt = NULL;
          geoMidPt = new LLPoint();
          err |= direct(locus->geoStart, locus->geoAz, delta, geoMidPt, EPS);

          //project geodesic mid point onto locus to find locus midpoint
          originCenter = NULL;
          originCenter = new LLPoint();
          err |= projectToLocus(*locus, *geoMidPt, originCenter, &crsGeoMidLocMid, &distGeoMidLocMid, TOL, EPS);
        }

        //find arc radius
        err |= invDist(*originCenter, locus->locusStart, &radius, EPS);//should be able to use either locusStart or locusEnd due to the special geometry here

        //find arc start/end points
        originArcStart = NULL;
        originArcStart = new LLPoint();
        originArcEnd = NULL;
        originArcEnd = new LLPoint();
        err |= direct(*originCenter, 0,radius, originArcStart, EPS);
        err |= direct(*originCenter, 0,radius, originArcEnd, EPS);

        //create the arc
        originArc = NULL;
        originArc = new Arc();
        err |= createArc(originArc, *originCenter, *originArcStart, *originArcEnd, dir, TOL, EPS);

        /*
         * Determine the case type we are trying to create:
         * 		type 0: zero intersections by a horizontal translation away from the locus
         * 		type 1: one intersection by a horizontal translation away from the locus
         * 		type 2: two intersections that are not locus start/end points by a horizontal translation away from the locus
         * 		type 3: two intersections that are the locus start/end points
         */
        caseType = rand() % 4;

        /* We are only considering types 0 and 1 in this test set */
        if (caseType > 1)
        {
          testCase++;
          continue;
        }
        else if (caseType == 0)
        {
          numType0++;
          if (numType0 > maxType0)
          {
            testCase++;
            continue;
          }
        }

        switch(caseType)
        {
	  case 0://type 0: zero intersections

          if (thetaDeg != 90)
	  {
	    //project the arc center onto the locus, this should be the locus mid point
	    tempPt1 = NULL;
	    tempPt1 = new LLPoint();
	    //tempDist is the distance from the arc center to the locus
	    err |= projectToLocus(*locus, originArc->centerPoint, tempPt1, &tempCrs1, &tempDist, TOL, EPS);

	    //Calculate the difference so we know the lower bound of the horizontal translation away from the locus
	    delta = originArc->radius - fabs(tempDist);  //nm, fabs in case the distance is signed

	    //Determine the forward course (tempCrs1) from the arc center to its projection onto the locus
	    err |= invCrs(originArc->centerPoint, *tempPt1, &tempCrs1, &tempCrs2, EPS);

	    //We want to move in the opposite direction of tempCrs1 a minimum distance of delta + 1 [nm]
	    //Find the translated arc center point
	    tempDist = delta + (rand() % maxHorizTranslationDist + 1);
	    center = NULL;
	    center = new LLPoint();
	    err |= direct(originArc->centerPoint, tempCrs1 + M_PI, tempDist, center, EPS);
	  }
	  else
	  {
	    center = NULL;
	    center = new LLPoint();
	    err |= direct(originArc->centerPoint, crsGeoMidLocMid + M_PI, radius+1.0, center, EPS);
 	  }

  	  //Create the test case arc
	  arcStart = NULL;
	  arcStart = new LLPoint();
	  arcEnd = NULL;
	  arcEnd = new LLPoint();
	  err |= direct(*center, 0,originArc->radius, arcStart, EPS);
	  err |= direct(*center, 0,originArc->radius, arcEnd, EPS);

	  arc = NULL;
	  arc = new Arc();
	  err |= createArc(arc, *center, *arcStart, *arcEnd, dir, TOL, EPS);

	  //set the expected results
	  expNX = 0;
          X1 = new LLPoint();
          X2 = new LLPoint();
	  err |= createPt(X1,0,0);//zero out X1
	  err |= createPt(X2,0,0);//zero out X2

	  break;
	  case 1: //type 1: one intersection by creating a locus tangent to the arc

	  tempPt1 = NULL;
	  tempPt1 = new LLPoint();
	  err |= direct(originArc->centerPoint, geoCrs, radius, tempPt1, EPS);
	  err |= invCrs(*tempPt1, originArc->centerPoint, &tempCrs1, &tempCrs2, EPS);
	  crsTangent = crsZeroToTwoPI(tempCrs1 + originArc->dir*M_PI_2);
	  //create locus
	  slopeLocAngle = 10.0 * DEG2RAD;
	  d1 = 0.01 * randDist();
	  d2 = 0.01 * randDist();
	  d3 = 0.01 * randDist();
	  createLocusThroughPoint(locus, *tempPt1, crsTangent, slopeLocAngle, d1, d2, d3);
	  center = NULL;
	  center = new LLPoint();
	  center->latitude = originArc->centerPoint.latitude;
	  center->longitude = originArc->centerPoint.longitude;

	  //Create the test case arc
	  arcStart = NULL;
	  arcStart = new LLPoint();
	  arcEnd = NULL;
	  arcEnd = new LLPoint();
	  err |= direct(*center, 0,originArc->radius, arcStart, EPS);
	  err |= direct(*center, 0,originArc->radius, arcEnd, EPS);

	  arc = NULL;
	  arc = new Arc();
	  err |= createArc(arc, *center, *arcStart, *arcEnd, dir, TOL, EPS);

	  //set the expected results
	  expNX = 1;
	  X1 = tempPt1; //locus midpoint
          X2 = new LLPoint();
	  err |= createPt(X2,0,0);//zero out X2

	  break;
        }//switch


	sprintf(testname,"TEST%i-%i",testCase,caseType);

	if ( getMaskedError(err, getMaskAll()) )
	{
		printf("\n--------------------------------\n");
		printf("\n%s\n", testname);
		printf("Error occurred in pre-locusArcIntx err=0x%lx\n",err);
		printf("\nLocus:\n");
		displayBoundaryLocus(*locus, 0);
		printf("\nArc Center:\n");
		displayBoundaryLLPoint(*center,1);
		printf("\nArc Radius:%.20lf [nm]", radius);
		printf("\ntol=%le",TOL);
		printf("\neps=%le",EPS);
		printf("\n\nExpected Results:");
		if (expNX > 0)
		{
			printf("\nX1:");
			displayBoundaryLLPoint(*X1,0);
		}
		if (expNX > 1)
		{
			printf("X2:");
			displayBoundaryLLPoint(*X2,0);
		}
		printf("\n--------------------------------\n");
		testCaseCount++;
		errorCount++;
		setupFailureCount++;
		continue;
	}

	//call the function under test with locus and Arc
	nX = 0;
	intx[0].latitude = 0.0;
	intx[0].longitude = 0.0;
	intx[1].latitude = 0.0;
	intx[1].longitude = 0.0;

	err |= locusArcIntx(*locus, *center, radius, intx, &nX, TOL, EPS);

	if ( getMaskedError(err, getMaskAll()) )
	{
		printf("\n--------------------------------\n");
		printf("\n%s\n", testname);
		printf("Error occurred in locusArcIntx err=0x%lx\n",err);
		printf("\nLocus:\n");
		displayBoundaryLocus(*locus, 0);
		printf("\nArc Center:\n");
		displayBoundaryLLPoint(*center,1);
		printf("\nArc Radius:%.20lf [nm]", radius);
		printf("\ntol=%le",TOL);
		printf("\neps=%le",EPS);
		printf("\n\nExpected Results:");
		if (expNX > 0)
		{
			printf("\nX1:");
			displayBoundaryLLPoint(*X1,0);
		}
		if (expNX > 1)
		{
			printf("X2:");
			displayBoundaryLLPoint(*X2,0);
		}
		printf("Actual Results:");
		printf("\nintx[0]:");
		displayBoundaryLLPoint(intx[0],0);
		printf("intx[1]:");
		displayBoundaryLLPoint(intx[1],0);
		printf("\n--------------------------------\n");
		failedCount++;
		errorCount++;
		testCaseCount++;
		continue;
	}

	//Verify results against expected intersection points X1 and X2
	if (expNX == nX)
	{
		if (expNX > 1)
		{
			if
			(
				( ptsAreSame(intx[0],*X2, TESTTOL) || ptsAreSame(intx[1],*X2, TESTTOL) ) &&
				( ptsAreSame(intx[0],*X1, TESTTOL) || ptsAreSame(intx[1],*X1, TESTTOL) )
			)
				passed = 1;
		}
		else if (expNX > 0)
		{
			if ( ptsAreSame(intx[0],*X1, TESTTOL) || ptsAreSame(intx[1],*X1, TESTTOL) )
				passed = 1;
		}
		else
			passed = 1;
	}

	if (!passed)
	{
		printf("\n--------------------------------\n");
		printf("\n%s FAILED\n", testname);
		printf("\ntol=%le",TOL);
		printf("\neps=%le",EPS);
		printf("\n\nLocus:\n");
                printLocus(" ",*locus);
		//displayBoundaryLocus(*locus, 0);
		printf("\nArc Center:\n");
		//displayBoundaryLLPoint(*center,0);
                printLLPoint(" ",*center);
		printf("Arc Radius[nm]: %.20lf", radius);
		printf("\n\nExpected Results:");
		printf("\nexpNX:%i", expNX);
		if (expNX > 0)
		{
			printf("\nX1:");
			//displayBoundaryLLPoint(*X1,0);
                        printLLPoint(" ",*X1);
		}
		if (expNX > 1)
		{
			printf("X2:");
			//displayBoundaryLLPoint(*X2,0);
                        printLLPoint(" ",*X2);
		}
		printf("Actual Results:");
		printf("\nnX:%i",nX);
		printf("\nintx[0]:");
		//displayBoundaryLLPoint(intx[0],0);
                printLLPoint(" ",intx[0]);
		printf("intx[1]:");
		//displayBoundaryLLPoint(intx[1],0);
                printLLPoint(" ",intx[1]);
		printf("\n--------------------------------\n");
		failedCount++;
	} else
		passedCount++;

	testCaseCount++;
        testCase++;

    }//while

    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("Finish testLocusArcIntx_Set2\n");

    return set;

}

/*
 * NAME: testLocusArcIntx_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the locusArcIntx function.
 *
 * 		This function runs all the test cases from all the test sets and returns
 * 		the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testLocusArcIntx_AllSets(TestSuite) - A test suite with the folling metrics:
 *  	testCases(int): The cumulative number of test cases executed.
 * 		pass(int): The cumulative number of test cases that passed.
 * 		fail(int): The cumulative number of test cases that failed.
 * 		setupFailures(int): The cumulative number of tests cases that generated an error prior
 * 					to actually executing the given test case.
 * 		errors(int): The cumulative number of test cases that generated an error.  Test cases that
 * 					generated an error prior to executing the given test are included in
 * 					this count.  Test cases that generated an error while executing the
 * 					given test are included as well.
 * 		unverified(int): The cumulative number of invalid records encountered while processing
 * 					the set of test data.
 * 		testingSuiteError(int): Return 1 if unable to execute one or more testing functions.
 * 					Return 0 otherwise.
 *
 */
TestSuite testLocusArcIntx_AllSets()
{

	TestSuite suite;
	TestSet set1;
	TestSet set2;

    printf("\nStart testLocusArcIntx_AllSets\n");

    suite = newTestSuite("testLocusArcIntx_AllSets");

    set1 = testLocusArcIntx_Set1();
    addTestSet(set1,&suite);

    set2 = testLocusArcIntx_Set2();
    addTestSet(set2,&suite);


    displayTestSuite(suite);

    printf("Finish testLocusArcIntx_AllSets\n\n\n");

    return suite;
}

/*
 * NAME: testLocusIntx_Set1
 *
 * DESCRIPTION:
 * 		This function is used to test the locusIntx function.
 *
 * 		This function runs the test data created by Balakrishna Babu and uses
 * 		his original test plan to execute the actual tests.
 *
 *      1. Select an arbitrary center of interest
 *      2. Select 4 azimuths and 4 distances to locate the start and end
 *         points of the defining geodesics for two loci.
 *      3. Select start and end distances for the two loci.
 *      4. Create the two loci.
 *      5. Find locus intersect using INFINITE lineType for both loci.
 *      6. Find locus intersect using all lineTypes for both loci.
 *      7. Verify if NO_INTERSECTION_ERR is expected.
 *      8. Verify if the intesection point lies on both loci.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testLocusIntx_Set1(TestSet) - A test set with the folling metrics:
 *  	testCases(int): The number of test cases executed.
 * 		pass(int): The number of test cases that passed.
 * 		fail(int): The number of test cases that failed.
 * 		setupFailures(int): The number of tests cases that generated an error prior
 * 					to actually executing the given test case.
 * 		errors(int): The number of test cases that generated an error.  Test cases that
 * 					generated an error prior to executing the given test are included in
 * 					this count.  Test cases that generated an error while executing the
 * 					given test are included as well.
 * 		unverified(int): The number of invalid records encountered while processing
 * 					the set of test data.
 * 		testingError(int): Return 1 if unable to execute testing function.  Return 0 otherwise.
 *
 */
TestSet testLocusIntx_Set1()
{

    int testNum=0;
    char testname[256];
    double DEG2RAD = M_PI / 180.0;

    double latC, lonC;
    double az1, d1, az2, d2, az3, d3, az4, d4;
    double startDist1, endDist1;
    double startDist2, endDist2;

    LLPoint center, proj;
    LLPoint geoStart1, geoEnd1;
    LLPoint geoStart2, geoEnd2;
    LineType lineType1, lineType2;  //0=SEGMENT,1=SEMIINFINITE,2=INFINITE
    Locus locus1, locus2;
    Locus* locusPtr1;
    Locus* locusPtr2;

    LLPoint intx, intxINF;
    double crsC1, crsC2, crsC3, crsC4;
    int expectValidIntx1, expectValidIntx2;
    int passedCount=0, failedCount=0;
    int testCaseCount = 0, setupFailureCount = 0, errorCount = 0, unverifiedCount = 0;
    ErrorSet err=0;

    TestSet set;

    long newSeed=20080424;

    printf("Start testLocusIntx_Set1\n");

    set = newTestSet("testLocusIntx_Set1");

	srand(newSeed);  //Initialize the random number generator

	locusPtr1 = &locus1;
	locusPtr2 = &locus2;

	while (testNum<=1000)
    {
    	err = 0;

    	//Select an arbitrary center of interest
    	latC = randLat();
    	lonC = randLon();

    	//Select 4 azimuths and 4 distances to locate the start and end points
    	//of the defining geodesics for two loci.
    	az1 = randAzimuth();
    	d1 = 0.01*randDist();
    	az2 = randAzimuth();
    	d2 = 0.01*randDist();
    	az3 = randAzimuth();
    	d3 = 0.01*randDist();
    	az4 = randAzimuth();
    	d4 = 0.01*randDist();

    	//Select start and end distances for the two loci.
    	startDist1 = 0.001*randDist();
    	endDist1 = 0.001*randDist();
    	startDist2 = 0.001*randDist();
    	endDist2 = 0.001*randDist();
    	lineType1 = INFINITE;
    	lineType2 = INFINITE;
    	testNum++;

        sprintf(testname,"TEST%-d",testNum);

        center.latitude = latC * DEG2RAD;
        center.longitude = lonC * DEG2RAD;
        crsC1 = az1 * DEG2RAD;
        crsC2 = az2 * DEG2RAD;
        crsC3 = az3 * DEG2RAD;
        crsC4 = az4 * DEG2RAD;

        err |= direct(center, crsC1, d1, &geoStart1, EPS);
        err |= direct(center, crsC2, d2, &geoEnd1, EPS);
        err |= direct(center, crsC3, d3, &geoStart2, EPS);
        err |= direct(center, crsC4, d4, &geoEnd2, EPS);

        //Create the two loci.
        err |= createLocus(locusPtr1,geoStart1,geoEnd1,startDist1,endDist1,lineType1,TOL,EPS);
        err |= createLocus(locusPtr2,geoStart2,geoEnd2,startDist2,endDist2,lineType2,TOL,EPS);

        //Find locus intersect using INFINITE lineType for both loci.
        intx.latitude = 0.0;
        intx.longitude = 0.0;
        err |= locusIntx(locus1,locus2,&intxINF,TOL,EPS);
        if ((err > 0) && (err != NO_INTERSECTION_ERR))
          continue;

        //Find locus intersect using all lineTypes for both loci.
        for (int locus1lineType=0; locus1lineType<3; locus1lineType++) {
         locus1.lineType = static_cast<LineType>(locus1lineType);
        	if (ptIsOnLocus(locus1,intxINF,&proj,&err,TOL,EPS)) {
        		expectValidIntx1 = 1;  //true
        	}
        	else {
        		expectValidIntx1 = 0;  //false
        	}
            for (int locus2lineType=0; locus2lineType<3; locus2lineType++) {
               locus2.lineType = static_cast<LineType>(locus2lineType);
            	if (ptIsOnLocus(locus2,intxINF,&proj,&err,TOL,EPS)) {
            		expectValidIntx2 = 1;  //true
            	}
            	else {
            		expectValidIntx2 = 0;  //false
            	}

                if ( getMaskedError(err, getMaskAll()) )
                {
                	printf("Error occurred in pre-locusIntx testNum %d err=0x%lx\n",testNum,err);
                    setupFailureCount++;
                    errorCount++;
                    testCaseCount++;
                	continue;
                }

            	intx.latitude = 0.0;
                intx.longitude = 0.0;
                err |= locusIntx(locus1,locus2,&intx,TOL,EPS);

                if ( getMaskedError(err, getMaskAll()) )
                {
                	printf("Error occurred in locusIntx err=0x%lx\n",err);
                	failedCount++;
                    errorCount++;
                    testCaseCount++;
                	continue;
                }

                //Verify if NO_INTERSECTION_ERR is expected.
                if ( (expectValidIntx1 == 0) || (expectValidIntx2 == 0) ) {
                	if (err == NO_INTERSECTION_ERR) {
                		if (PRINT_PASSED_CASES) printf("%s C: %14.8f %14.8f  az1=%14.8f d1=%14.8f  az2=%14.8f d2=%14.8f  startDist1=%14.8f endDist1=%14.8f\n",
                          		testname,latC,lonC,az1,d1,az2,d2,startDist1,endDist1);
                		if (PRINT_PASSED_CASES) printf("%s C: %14.8f %14.8f  az3=%14.8f d3=%14.8f  az4=%14.8f d4=%14.8f  startDist2=%14.8f endDist2=%14.8f\n",
                          		testname,latC,lonC,az3,d3,az4,d4,startDist2,endDist2);
                		if (PRINT_PASSED_CASES) printLLPoint("geoStart1",geoStart1);
                		if (PRINT_PASSED_CASES) printLLPoint("geoEnd1",geoEnd1);
                		if (PRINT_PASSED_CASES) printLLPoint("geoStart2",geoStart2);
                		if (PRINT_PASSED_CASES) printLLPoint("geoEnd2",geoEnd2);

                		if (PRINT_PASSED_CASES) printf("%s passed lineType1=%i lineType2=%i expIntx1=%i expIntx2=%i\n",
	                			testname,locus1.lineType,locus2.lineType,expectValidIntx1,expectValidIntx2);
	                	passedCount++;
                        testCaseCount++;
                        err = 0;
                	}
                	else {
                        printf("%s C: %14.8f %14.8f  az1=%14.8f d1=%14.8f  az2=%14.8f d2=%14.8f  startDist1=%14.8f endDist1=%14.8f\n",
                          		testname,latC,lonC,az1,d1,az2,d2,startDist1,endDist1);
                        printf("%s C: %14.8f %14.8f  az3=%14.8f d3=%14.8f  az4=%14.8f d4=%14.8f  startDist2=%14.8f endDist2=%14.8f\n",
                          		testname,latC,lonC,az3,d3,az4,d4,startDist2,endDist2);
                        printLLPoint("geoStart1",geoStart1);
                        printLLPoint("geoEnd1",geoEnd1);
                        printLLPoint("geoStart2",geoStart2);
                        printLLPoint("geoEnd2",geoEnd2);

	                	printf("%s failed lineType1=%i lineType2=%i expIntx1=%i expIntx2=%i\n",
	                			testname,locus1.lineType,locus2.lineType,expectValidIntx1,expectValidIntx2);
	                	failedCount++;
                        testCaseCount++;
                	}
                }
                //Verify if the intersection point lies on both loci.
                else {
	                if ( ptIsOnLocus(locus1,intx,&proj,&err,TESTTOL,EPS) &&
	                     ptIsOnLocus(locus2,intx,&proj,&err,TESTTOL,EPS) ) {

	                	if ( getMaskedError(err, getMaskAll()) )
	                	{
	                		printf("Error occurred in post-locusIntx err=0x%lx\n",err);
	                		failedCount++;
                            errorCount++;
                            testCaseCount++;
	                		continue;
	                	}
                		if (PRINT_PASSED_CASES) printf("%s C: %14.8f %14.8f  az1=%14.8f d1=%14.8f  az2=%14.8f d2=%14.8f  startDist1=%14.8f endDist1=%14.8f\n",
                          		testname,latC,lonC,az1,d1,az2,d2,startDist1,endDist1);
                		if (PRINT_PASSED_CASES) printf("%s C: %14.8f %14.8f  az3=%14.8f d3=%14.8f  az4=%14.8f d4=%14.8f  startDist2=%14.8f endDist2=%14.8f\n",
                          		testname,latC,lonC,az3,d3,az4,d4,startDist2,endDist2);
                		if (PRINT_PASSED_CASES) printLLPoint("geoStart1",geoStart1);
                		if (PRINT_PASSED_CASES) printLLPoint("geoEnd1",geoEnd1);
                		if (PRINT_PASSED_CASES) printLLPoint("geoStart2",geoStart2);
                		if (PRINT_PASSED_CASES) printLLPoint("geoEnd2",geoEnd2);

	                	if (PRINT_PASSED_CASES) printf("%s passed %14.8f %14.8f\n",
	                			testname,intx.latitude/DEG2RAD,intx.longitude/DEG2RAD);
	                	passedCount++;
                        testCaseCount++;
	                }
	                else {
                        printf("%s C: %14.8f %14.8f  az1=%14.8f d1=%14.8f  az2=%14.8f d2=%14.8f  startDist1=%14.8f endDist1=%14.8f\n",
                          		testname,latC,lonC,az1,d1,az2,d2,startDist1,endDist1);
                        printf("%s C: %14.8f %14.8f  az3=%14.8f d3=%14.8f  az4=%14.8f d4=%14.8f  startDist2=%14.8f endDist2=%14.8f\n",
                          		testname,latC,lonC,az3,d3,az4,d4,startDist2,endDist2);
                        printLLPoint("geoStart1",geoStart1);
                        printLLPoint("geoEnd1",geoEnd1);
                        printLLPoint("geoStart2",geoStart2);
                        printLLPoint("geoEnd2",geoEnd2);

	                	printf("%s failed %14.8f %14.8f err=%ld\n",
	                			testname,intx.latitude/DEG2RAD,intx.longitude/DEG2RAD,err);
	                	failedCount++;
                        testCaseCount++;
                        err = 0;
	                }
                }
            }
        }
    }
    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("Finish testLocusIntx_Set1\n");
    return set;
}

/*
 * NAME: testLocusIntx_Set2
 *
 * DESCRIPTION:
 * 		This function is used to test the locusIntx function.
 *
 * 		This function runs the test data created by Balakrishna Babu and uses
 * 		his original test plan to execute the actual tests.
 *
 * 		This function was originally called testLocusIntxB but was
 * 		renamed as part of the test rearchitecturing of geolib.
 *
 *      1. Select an intersection point using random.
 *      2. Construct a locus that passes through the selected intersection point.
 *            2.1. Azimuth and a distance to find the point on the defining geodesic
 *            2.2. Compute azimuth to the point on locus from point on geodesic
 *            2.3. Compute two azimuths 90 degress on either side of azimuth to point
 *                 and select two distances to locate the start and end points of the
 *                 defining geodesic
 *            2.4. Select startDist and compute endDist
 *            2.5. Create locus that passed through the selected intersection point
 *      3. Construct a second locus that passes through the same intersection point.
 *      4. Find locus intersect using all lineTypes for both loci.
 *      5. Verify if NO_INTERSECTION_ERR is expected.
 *      6. Verify if the intescetion point lies on both loci.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testLocusIntx_Set2(TestSet) - A test set with the folling metrics:
 *  	testCases(int): The number of test cases executed.
 * 		pass(int): The number of test cases that passed.
 * 		fail(int): The number of test cases that failed.
 * 		setupFailures(int): The number of tests cases that generated an error prior
 * 					to actually executing the given test case.
 * 		errors(int): The number of test cases that generated an error.  Test cases that
 * 					generated an error prior to executing the given test are included in
 * 					this count.  Test cases that generated an error while executing the
 * 					given test are included as well.
 * 		unverified(int): The number of invalid records encountered while processing
 * 					the set of test data.
 * 		testingError(int): Return 1 if unable to execute testing function.  Return 0 otherwise.
 *
 */
//TODO seperate the test data creation logic from the actual testing logic
TestSet testLocusIntx_Set2()
{

    int testNum=0;
    char testname[256];
    double DEG2RAD = M_PI / 180.0;

    double latX, lonX;
    double az1, d1, d2, d3;
    double startDist1, endDist1, slopeAngle1;
    double startDist2, endDist2, slopeAngle2;

    LLPoint proj, proj1, proj2;
    LLPoint geoStart1, geoEnd1;
    LLPoint geoStart2, geoEnd2;
    LineType lineType1, lineType2;  //0=SEGMENT,1=SEMIINFINITE,2=INFINITE
    Locus locus1, locus2;
    Locus* locusPtr1;
    Locus* locusPtr2;

    double fcrs,bcrs;
    LLPoint intxExp, intx, intxINF;
    double crsC1, crsC2, crsC3;
    int expectValidIntx1, expectValidIntx2;
    int passedCount=0, failedCount=0;
    int testCaseCount = 0, errorCount = 0, setupFailureCount = 0, unverifiedCount = 0;
    ErrorSet err=0;

    TestSet set;

    long newSeed=20080501;

    printf("Start testLocusIntx_Set2\n");

    set = newTestSet("testLocusIntx_Set2");

	srand(newSeed);  //Initialize the random number generator

	locusPtr1 = &locus1;
	locusPtr2 = &locus2;

	while (testNum<=1000)
    {
    	err = 0;

    	testNum++;
        sprintf(testname,"TEST%-d",testNum);

        //Select an intersection point
    	latX = randLat();
    	lonX = randLon();
        intxExp.latitude = latX * DEG2RAD;
        intxExp.longitude = lonX * DEG2RAD;

    	//Construct a locus that passes through the selected intersection point
        //1. Azimuth and a distance to find the point on the defining geodesic
    	az1 = randAzimuth();
    	d1 = 0.001*randDist();
        crsC1 = az1 * DEG2RAD;
        err |= direct(intxExp, crsC1, d1, &proj1, EPS);

        //2. Compute azimuth to the point on locus from point on geodesic
        err |= invCrs(proj1,intxExp,&fcrs,&bcrs,EPS);

        //3. Compute two azimuths 90 degress on either side of azimuth to point
        //   and select two distances to locate the start and end points of the
        //   defining geodesic
        crsC2 = fcrs + M_PI_2;
        if (crsC2 > M_2PI) {
        	crsC2 -= M_2PI;
        }
        crsC3 = fcrs - M_PI_2;
        if (crsC3 < 0.0) {
        	crsC3 += M_2PI;
        }
    	d2 = 0.01*randDist();
    	d3 = 0.01*randDist();
        err |= direct(proj1, crsC2, d2, &geoStart1, EPS);
        err |= direct(proj1, crsC3, d3, &geoEnd1, EPS);

        //4. Select startDist and compute endDist
    	startDist1 = 0.001*randDist();
    	endDist1 = d1 + d3 * (d1-startDist1)/d2;
    	slopeAngle1 = atan((d1-startDist1)/d2);
    	lineType1 = INFINITE;

    	//5. Create locus that passed through the selected intersection point
        err |= createLocus(locusPtr1,geoStart1,geoEnd1,startDist1,endDist1,lineType1,TOL,EPS);

    	//Construct a locus that passes through the selected intersection point
        //1. Azimuth and a distance to find the point on the defining geodesic
    	az1 = randAzimuth();
    	d1 = 0.001*randDist();
        crsC1 = az1 * DEG2RAD;
        err |= direct(intxExp, crsC1, d1, &proj2, EPS);

        //2. Compute azimuth to the point on locus from point on geodesic
        err |= invCrs(proj2,intxExp,&fcrs,&bcrs,EPS);

        //3. Compute two azimuths 90 degress on either side of azimuth to point
        //   and select two distances to locate the start and end points of the
        //   defining geodesic
        crsC2 = fcrs + M_PI_2;
        if (crsC2 > M_2PI) {
        	crsC2 -= M_2PI;
        }
        crsC3 = fcrs - M_PI_2;
        if (crsC3 < 0.0) {
        	crsC3 += M_2PI;
        }
    	d2 = 0.01*randDist();
    	d3 = 0.01*randDist();
        err |= direct(proj2, crsC2, d2, &geoStart2, EPS);
        err |= direct(proj2, crsC3, d3, &geoEnd2, EPS);

        //4. Select startDist and compute endDist
    	startDist2 = 0.001*randDist();
    	endDist2 = d1 + d3 * (d1-startDist2)/d2;
    	slopeAngle2 = atan((d1-startDist2)/d2);
    	lineType2 = INFINITE;

    	//5. Create locus that passed through the selected intersection point
        err |= createLocus(locusPtr2,geoStart2,geoEnd2,startDist2,endDist2,lineType2,TOL,EPS);

        intxINF.latitude = 0.0;
        intxINF.longitude = 0.0;
        err |= locusIntx(locus1,locus2,&intxINF,TOL,EPS);

        //Find locus intersect using all lineTypes for both loci.
        for (int locus1lineType=0; locus1lineType<3; locus1lineType++) {
           locus1.lineType = static_cast<LineType>(locus1lineType);
        	if (ptIsOnLocus(locus1,intxINF,&proj,&err,TOL,EPS)) {
        		expectValidIntx1 = 1;  //true
        	}
        	else {
        		expectValidIntx1 = 0;  //false
        	}
            for (int locus2lineType=0; locus2lineType<3; locus2lineType++) {
               locus2.lineType = static_cast<LineType>(locus2lineType);
            	if (ptIsOnLocus(locus2,intxINF,&proj,&err,TOL,EPS)) {
            		expectValidIntx2 = 1;  //true
            	}
            	else {
            		expectValidIntx2 = 0;  //false
            	}

                if ( getMaskedError(err, getMaskAll()) )
                {
                	printf("Error occurred in pre-locusIntx err=0x%lx\n",err);
                    errorCount++;
                    setupFailureCount++;
                    testCaseCount++;
                	continue;
                }

            	intx.latitude = 0.0;
                intx.longitude = 0.0;
                err |= locusIntx(locus1,locus2,&intx,TOL,EPS);

                if ( getMaskedError(err, getMaskAll()) )
                {
                	printf("Error occurred in locusIntx err=0x%lx\n",err);
                	failedCount++;
                    errorCount++;
                    testCaseCount++;
                	continue;
                }

                //Verify if NO_INTERSECTION_ERR is expected.
                if ( (expectValidIntx1 == 0) || (expectValidIntx2 == 0) ) {
                	if (err == NO_INTERSECTION_ERR) {
                		if (PRINT_PASSED_CASES) printf("%s latX=%14.8f lonX=%14.8f\n",
                          		testname,latX,lonX);
                		if (PRINT_PASSED_CASES) printLLPoint("Expected Intersection point",intxExp);
                		if (PRINT_PASSED_CASES) printLLPoint("point on defining geodesic",proj1);
                		if (PRINT_PASSED_CASES) printLLPoint("geoStart1",geoStart1);
                		if (PRINT_PASSED_CASES) printLLPoint("geoEnd1",geoEnd1);
                		if (PRINT_PASSED_CASES) printf("%s X: %14.8f %14.8f  az1=%14.8f d1=%14.8f d2=%14.8f d3=%14.8f\n",
                          		testname,latX,lonX,az1,d1,d2,d3);
                		if (PRINT_PASSED_CASES) printf("%s startDist1=%14.8f endDist1=%14.8f slopeAngle1=%14.8f\n",
                          		testname,startDist1,endDist1,slopeAngle1/DEG2RAD);
                		if (PRINT_PASSED_CASES) printLLPoint("locusStart1",locus1.locusStart);
                		if (PRINT_PASSED_CASES) printLLPoint("locusEnd1",locus1.locusEnd);
                		if (PRINT_PASSED_CASES)  printLLPoint("point on defining geodesic 2",proj2);
                		if (PRINT_PASSED_CASES)  printLLPoint("geoStart2",geoStart2);
                		if (PRINT_PASSED_CASES) printLLPoint("geoEnd2",geoEnd2);
                		if (PRINT_PASSED_CASES)  printf("%s X: %14.8f %14.8f  az1=%14.8f d1=%14.8f d2=%14.8f d3=%14.8f\n",
                          		testname,latX,lonX,az1,d1,d2,d3);
                		if (PRINT_PASSED_CASES) printf("%s startDist2=%14.8f endDist2=%14.8f slopeAngle2=%14.8f\n",
                          		testname,startDist2,endDist2,slopeAngle2/DEG2RAD);
                		if (PRINT_PASSED_CASES) printLLPoint("locusStart2",locus2.locusStart);
                		if (PRINT_PASSED_CASES)  printLLPoint("locusEnd2",locus2.locusEnd);
                		if (PRINT_PASSED_CASES) printLLPoint("intxINF",intxINF);

                		if (PRINT_PASSED_CASES) printf("%s passed lineType1=%i lineType2=%i expIntx1=%i expIntx2=%i slAng1=%8.2f slAng2=%8.2f\n",
	                			testname,locus1.lineType,locus2.lineType,expectValidIntx1,expectValidIntx2,
	                			slopeAngle1/DEG2RAD,slopeAngle2/DEG2RAD);
	                	passedCount++;
                        testCaseCount++;
                	}
                	else {
                        printf("%s latX=%14.8f lonX=%14.8f\n",
                          		testname,latX,lonX);
                        printLLPoint("Expected Intersection point",intxExp);
                        printLLPoint("point on defining geodesic",proj1);
                        printLLPoint("geoStart1",geoStart1);
                        printLLPoint("geoEnd1",geoEnd1);
                        printf("%s X: %14.8f %14.8f  az1=%14.8f d1=%14.8f d2=%14.8f d3=%14.8f\n",
                          		testname,latX,lonX,az1,d1,d2,d3);
                        printf("%s startDist1=%14.8f endDist1=%14.8f slopeAngle1=%14.8f\n",
                          		testname,startDist1,endDist1,slopeAngle1/DEG2RAD);
                        printLLPoint("locusStart1",locus1.locusStart);
                        printLLPoint("locusEnd1",locus1.locusEnd);
                        printLLPoint("point on defining geodesic 2",proj2);
                        printLLPoint("geoStart2",geoStart2);
                        printLLPoint("geoEnd2",geoEnd2);
                        printf("%s X: %14.8f %14.8f  az1=%14.8f d1=%14.8f d2=%14.8f d3=%14.8f\n",
                          		testname,latX,lonX,az1,d1,d2,d3);
                        printf("%s startDist2=%14.8f endDist2=%14.8f slopeAngle2=%14.8f\n",
                          		testname,startDist2,endDist2,slopeAngle2/DEG2RAD);
                        printLLPoint("locusStart2",locus2.locusStart);
                        printLLPoint("locusEnd2",locus2.locusEnd);
                        printLLPoint("intxINF",intxINF);

	                	printf("%s failed lineType1=%i lineType2=%i expIntx1=%i expIntx2=%i slAng1=%8.2f slAng2=%8.2f\n",
	                			testname,locus1.lineType,locus2.lineType,expectValidIntx1,expectValidIntx2,
            					slopeAngle1/DEG2RAD,slopeAngle2/DEG2RAD);
	                	failedCount++;
                        testCaseCount++;
                	}
                }
                //Verify if the intescetion point lies on both loci.
                else {
	                if ( ptIsOnLocus(locus1,intx,&proj,&err,TESTTOL,EPS) &&
	                     ptIsOnLocus(locus2,intx,&proj,&err,TESTTOL,EPS) ) {

	                	if ( getMaskedError(err, getMaskAll()) )
	                	{
	                		printf("Error occurred in locusIntx err=0x%lx\n",err);
	                		failedCount++;
                            errorCount++;
                            testCaseCount++;
	                		continue;
	                	}

                		if (PRINT_PASSED_CASES) printf("%s latX=%14.8f lonX=%14.8f\n",
                          		testname,latX,lonX);
                		if (PRINT_PASSED_CASES) printLLPoint("Expected Intersection point",intxExp);
                		if (PRINT_PASSED_CASES) printLLPoint("point on defining geodesic",proj1);
                		if (PRINT_PASSED_CASES) printLLPoint("geoStart1",geoStart1);
                		if (PRINT_PASSED_CASES) printLLPoint("geoEnd1",geoEnd1);
                		if (PRINT_PASSED_CASES) printf("%s X: %14.8f %14.8f  az1=%14.8f d1=%14.8f d2=%14.8f d3=%14.8f\n",
                          		testname,latX,lonX,az1,d1,d2,d3);
                		if (PRINT_PASSED_CASES) printf("%s startDist1=%14.8f endDist1=%14.8f slopeAngle1=%14.8f\n",
                          		testname,startDist1,endDist1,slopeAngle1/DEG2RAD);
                		if (PRINT_PASSED_CASES) printLLPoint("locusStart1",locus1.locusStart);
                		if (PRINT_PASSED_CASES) printLLPoint("locusEnd1",locus1.locusEnd);
                		if (PRINT_PASSED_CASES)  printLLPoint("point on defining geodesic 2",proj2);
                		if (PRINT_PASSED_CASES)  printLLPoint("geoStart2",geoStart2);
                		if (PRINT_PASSED_CASES) printLLPoint("geoEnd2",geoEnd2);
                		if (PRINT_PASSED_CASES)  printf("%s X: %14.8f %14.8f  az1=%14.8f d1=%14.8f d2=%14.8f d3=%14.8f\n",
                          		testname,latX,lonX,az1,d1,d2,d3);
                		if (PRINT_PASSED_CASES) printf("%s startDist2=%14.8f endDist2=%14.8f slopeAngle2=%14.8f\n",
                          		testname,startDist2,endDist2,slopeAngle2/DEG2RAD);
                		if (PRINT_PASSED_CASES) printLLPoint("locusStart2",locus2.locusStart);
                		if (PRINT_PASSED_CASES) printLLPoint("locusEnd2",locus2.locusEnd);
                		if (PRINT_PASSED_CASES) printLLPoint("intxINF",intxINF);

	                	if (PRINT_PASSED_CASES) printf("%s passed %14.8f %14.8f\n",
	                			testname,intx.latitude/DEG2RAD,intx.longitude/DEG2RAD);
	                	passedCount++;
                        testCaseCount++;
	                }
	                else {
                        printf("%s latX=%14.8f lonX=%14.8f\n",
                          		testname,latX,lonX);
                        printLLPoint("Expected Intersection point",intxExp);
                        printLLPoint("point on defining geodesic",proj1);
                        printLLPoint("geoStart1",geoStart1);
                        printLLPoint("geoEnd1",geoEnd1);
                        printf("%s X: %14.8f %14.8f  az1=%14.8f d1=%14.8f d2=%14.8f d3=%14.8f\n",
                          		testname,latX,lonX,az1,d1,d2,d3);
                        printf("%s startDist1=%14.8f endDist1=%14.8f slopeAngle1=%14.8f\n",
                          		testname,startDist1,endDist1,slopeAngle1/DEG2RAD);
                        printLLPoint("locusStart1",locus1.locusStart);
                        printLLPoint("locusEnd1",locus1.locusEnd);
                        printLLPoint("point on defining geodesic 2",proj2);
                        printLLPoint("geoStart2",geoStart2);
                        printLLPoint("geoEnd2",geoEnd2);
                        printf("%s X: %14.8f %14.8f  az1=%14.8f d1=%14.8f d2=%14.8f d3=%14.8f\n",
                          		testname,latX,lonX,az1,d1,d2,d3);
                        printf("%s startDist2=%14.8f endDist2=%14.8f slopeAngle2=%14.8f\n",
                          		testname,startDist2,endDist2,slopeAngle2/DEG2RAD);
                        printLLPoint("locusStart2",locus2.locusStart);
                        printLLPoint("locusEnd2",locus2.locusEnd);
                        printLLPoint("intxINF",intxINF);

	                	printf("%s failed %14.8f %14.8f\n",
	                			testname,intx.latitude/DEG2RAD,intx.longitude/DEG2RAD);
	                	failedCount++;
                        testCaseCount++;
	                }
                }
            }
        }
    }
    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("Finish testLocusIntx_Set2\n");
    return set;
}

/*
 * NAME: testLocusIntx_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the locusIntx function.
 *
 * 		This function runs all the test cases from all the test sets and returns
 * 		the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testLocusIntx_AllSets(TestSuite) - A test suite with the folling metrics:
 *  	testCases(int): The cumulative number of test cases executed.
 * 		pass(int): The cumulative number of test cases that passed.
 * 		fail(int): The cumulative number of test cases that failed.
 * 		setupFailures(int): The cumulative number of tests cases that generated an error prior
 * 					to actually executing the given test case.
 * 		errors(int): The cumulative number of test cases that generated an error.  Test cases that
 * 					generated an error prior to executing the given test are included in
 * 					this count.  Test cases that generated an error while executing the
 * 					given test are included as well.
 * 		unverified(int): The cumulative number of invalid records encountered while processing
 * 					the set of test data.
 * 		testingSuiteError(int): Return 1 if unable to execute one or more testing functions.
 * 					Return 0 otherwise.
 *
 */
TestSuite testLocusIntx_AllSets()
{

	TestSuite suite;
	TestSet set1;
	TestSet set2;

    printf("\nStart testLocusIntx_AllSets\n");

    suite = newTestSuite("testLocusIntx_AllSets");

    set1 = testLocusIntx_Set1();
    addTestSet(set1,&suite);

    set2 = testLocusIntx_Set2();
    addTestSet(set2,&suite);

    displayTestSuite(suite);

    printf("Finish testLocusIntx_AllSets\n\n\n");

    return suite;
}

/*
 * NAME: testArcTanToTwoLoci_Set1
 *
 * DESCRIPTION:
 * 		This function is used to test the arcTanToTwoLoci function.
 *
 * 		This function runs the test data created by Balakrishna Babu and uses
 * 		his original test plan to execute the actual tests.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testArcTanToTwoLoci_Set1(TestSet) - A test set with the folling metrics:
 *  	testCases(int): The number of test cases executed.
 * 		pass(int): The number of test cases that passed.
 * 		fail(int): The number of test cases that failed.
 * 		setupFailures(int): The number of tests cases that generated an error prior
 * 					to actually executing the given test case.
 * 		errors(int): The number of test cases that generated an error.  Test cases that
 * 					generated an error prior to executing the given test are included in
 * 					this count.  Test cases that generated an error while executing the
 * 					given test are included as well.
 * 		unverified(int): The number of invalid records encountered while processing
 * 					the set of test data.
 * 		testingError(int): Return 1 if unable to execute testing function.  Return 0 otherwise.
 *
 */
TestSet testArcTanToTwoLoci_Set1()
{

    double DEG2RAD = M_PI / 180.0;
    int testNum=0;
    char testname[256];
    double latxExp, lonxExp;
    double az1, r1;
    Locus locus1, locus2;
    LLPoint pp[3], qq[3];
    double radius1;
//    LineType lineType1=INFINITE;

    int ir=0;
    int irMax=10;
    double radius[10] = {0.01,0.1,0.2,0.5,1.0,5.0,10.0,20.0,50.0,100.0};

    int iz=0;
    int izMax=20;

    LLPoint startExp, centerExp, endExp;
    int dirExp;  //CLOCKWISE = 1, COUNTERCLOCKWISE = -1
    LLPoint start, center, end;
    ArcDirection dir;

    double slopeLocAngle1, slopeLocAngle2;  //slope angles for both loci
    double d1, d2, d3, e1, e2, e3;  //distances to define locus 1 and 2

    double crsCS, crsSC, distCS;
    double crsCE, crsEC, distCE;
    double crsSTangent, crsETangent;
    double centerDistErr, startDistErr, endDistErr;
    int passedCount=0, failedCount=0;
    int testCaseCount = 0, setupFailureCount = 0, errorCount = 0, unverifiedCount = 0;
    ErrorSet err=0;

    double azmin, azmax, azr;

    TestSet set;

    long newSeed=20080703;

    printf("Start testArcTanToTwoLoci_Set1\n");

    set = newTestSet("testArcTanToTwoLoci_Set1");

    err = 0;
    pp[0].latitude = 45.0*DEG2RAD;
    pp[0].longitude = -77.0*DEG2RAD;
    pp[1].latitude = 45.0*DEG2RAD;
    pp[1].longitude = -75.0*DEG2RAD;
    pp[2].latitude = 45.0*DEG2RAD;
    pp[2].longitude = -73.0*DEG2RAD;
    qq[0].latitude = 43.0*DEG2RAD;
    qq[0].longitude = -75.0*DEG2RAD;
    qq[1].latitude = 45.0*DEG2RAD;
    qq[1].longitude = -75.0*DEG2RAD;
    qq[2].latitude = 47.0*DEG2RAD;
    qq[2].longitude = -75.0*DEG2RAD;
    radius1 = 50.0;

    srand(newSeed);  //Initialize the random number generator

    while (testNum<10)
    {
        latxExp = randLat();
        lonxExp = randLon();
        az1 = randAzimuth();
        testNum++;

        sprintf(testname,"TEST%-d",testNum);

        centerExp.latitude = latxExp * DEG2RAD;
        centerExp.longitude = lonxExp * DEG2RAD;
        crsCS = az1 * DEG2RAD;

        //Varying radius for the Arc
        for (ir=0; ir<irMax; ir++) {
        	r1 = radius[ir];

	        //same center and start points, end point varies
	        for (iz=0; iz<izMax; iz++) {

	        	err = 0;
	        	azmin = 0.01*DEG2RAD; //set by the accuracy limit of CrsIntersect
	        	azmax = 179.6*DEG2RAD; //Experimental accuracy limit of TangentFixedRadiusArc
	        	azr = azmin + (azmax - azmin) * (double) rand() / (double) RAND_MAX;
	        	if(iz < irMax/2) azr = -azr;
	        	//Set the direction of the turn or Arc
	            if (azr < 0.0)
	            {
	            	dirExp = COUNTERCLOCKWISE;  //COUNTERCLOCKWISE = -1
	            }
	            else
	            {
	                dirExp = CLOCKWISE;  //CLOCKWISE = 1,
	            }

	            //Set the azimuth from center to end of Arc
	            crsCE = crsZeroToTwoPI(crsCS + azr);

		        //Compute the start and end points of Arc
		        err |= direct(centerExp, crsCS, r1, &startExp, EPS);
		        err |= direct(centerExp, crsCE, r1, &endExp, EPS);

		        sprintf(testname,"TEST%-d-%-d-%-d",testNum,ir+1,iz+1);

		        //Compute course at start of Arc
		        err |= inverse(startExp, centerExp, &crsSC, &crsCS, &distCS, EPS);
		        crsSTangent = crsZeroToTwoPI(crsSC - dirExp*M_PI_2);

		        //Compute course at end of Arc
		        err |= inverse(endExp, centerExp, &crsEC, &crsCE, &distCE, EPS);
		        crsETangent = crsZeroToTwoPI(crsEC - dirExp*M_PI_2);

		        //Create loci
		        slopeLocAngle1 = 10.0 * DEG2RAD;
		        d1 = 0.01 * randDist(); //0-54 nm
		        d2 = 0.01 * randDist(); //0-54 nm
		        d3 = 0.01 * randDist(); //0-54 nm

		        createLocusThroughPoint(&locus1, startExp, crsSTangent,
		        						slopeLocAngle1,d1,d2,d3);

		        slopeLocAngle2 = -5.0 * DEG2RAD;
		        e1 = 0.01 * randDist(); //0-54 nm
		        e2 = 0.01 * randDist(); //0-54 nm
		        e3 = 0.01 * randDist(); //0-54 nm

		        createLocusThroughPoint(&locus2, endExp, crsETangent,
		        						slopeLocAngle2, e1,e2,e3);

		        if ( getMaskedError(err, getMaskAll()) )
		        {
		            printf("%s - Error occurred in pre-arcTanToTwoLoci err=0x%lx\n", testname, err);
                    setupFailureCount++;
                    errorCount++;
                    testCaseCount++;
		            continue;
		        }

		        center.latitude = 0.0;
		        center.longitude = 0.0;
		        start.latitude = 0.0;
		        start.longitude = 0.0;
		        end.latitude = 0.0;
		        end.longitude = 0.0;
		        //dir = 0;
		        err |= arcTanToTwoLoci(locus1,locus2,r1,
		                                          &center,&start,&end,&dir,TOL,EPS);

		        if ( getMaskedError(err, getMaskAll()) )
		        {
		            printf("%s - Error occurred in arcTanToTwoLoci err=0x%lx\n", testname, err);
		            failedCount++;
                    errorCount++;
                    testCaseCount++;
		            continue;
		        }

		        err |= invDist(center,centerExp,&centerDistErr,EPS);
		        err |= invDist(start,startExp,&startDistErr,EPS);
		        err |= invDist(end,endExp,&endDistErr,EPS);

		        if ( getMaskedError(err, getMaskAll()) )
		        {
		            printf("%s - Error occurred in post-arcTanToTwoLoci err=0x%lx\n", testname, err);
		            failedCount++;
                    errorCount++;
                    testCaseCount++;
		            continue;
		        }

		        if ( ptsAreSame(center,centerExp,TESTTOL) &&
		                ptsAreSame(start,startExp,TESTTOL) &&
		                ptsAreSame(end,endExp,TESTTOL) &&
		                (dirExp == dir) )
		        {
		        	if (PRINT_PASSED_CASES) printf("%s passed  dirExp= %i dir= %i  centerDistErr=%16.8e  startDistErr=%16.8e  endDistErr=%16.8e err=0x%lx\n",
			               testname,dirExp,dir,centerDistErr,startDistErr,endDistErr,err);
		        	if (PRINT_PASSED_CASES) printf("%s Center: %14.8f %14.8f  az1=%14.8f\n",
		                   testname,latxExp,lonxExp,az1);
		        	if (PRINT_PASSED_CASES) printf("%s r1=%14.8f  crsCS=%14.8f  crsCE=%14.8f\n",
			               testname,r1,crsCS/DEG2RAD,crsCE/DEG2RAD);
		        	if (PRINT_PASSED_CASES) printLLPoint("   startExp  ",startExp);
		        	if (PRINT_PASSED_CASES) printLLPoint("   centerExp ",centerExp);
		        	if (PRINT_PASSED_CASES) printLLPoint("   endExp    ",endExp);
		        	if (PRINT_PASSED_CASES) printLocus("   Locus1    ", locus1);
		        	if (PRINT_PASSED_CASES) printLocus("   Locus2    ", locus2);
		        	if (PRINT_PASSED_CASES) printLLPoint("   start    ",start);
		        	if (PRINT_PASSED_CASES) printLLPoint("   center   ",center);
		        	if (PRINT_PASSED_CASES) printLLPoint("   end      ",end);

		            passedCount++;
		        }
		        else
		        {
			        printf("%s failed  dirExp= %i dir= %i  centerDistErr=%16.8e  startDistErr=%16.8e  endDistErr=%16.8e err=0x%lx\n",
			               testname,dirExp,dir,centerDistErr,startDistErr,endDistErr,err);
		            printf("%s Center: %14.8f %14.8f  az1=%14.8f azr=%14.8f\n",
		                   testname,latxExp,lonxExp,az1,azr*180.0/M_PI);
			        printf("%s r1=%14.8f  crsCS=%14.8f  crsCE=%14.8f\n",
			               testname,r1,crsCS/DEG2RAD,crsCE/DEG2RAD);
			        printLLPoint("   startExp  ",startExp);
			        printLLPoint("   centerExp ",centerExp);
			        printLLPoint("   endExp    ",endExp);
			        printLocus("   Locus1    ", locus1);
			        printLocus("   Locus2    ", locus2);
			        printLLPoint("   start    ",start);
			        printLLPoint("   center   ",center);
			        printLLPoint("   end      ",end);
			        failedCount++;
		        }
		        testCaseCount++;

		    }
        }
    }

    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("\nFinish testArcTanToTwoLoci_Set1\n");

    return set;
}

/*
 * NAME: testArcTanToTwoLoci_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the arcTanToTwoLoci function.
 *
 * 		This function runs all the test cases from all the test sets and returns
 * 		the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testArcTanToTwoLoci_AllSets(TestSuite) - A test suite with the folling metrics:
 *  	testCases(int): The cumulative number of test cases executed.
 * 		pass(int): The cumulative number of test cases that passed.
 * 		fail(int): The cumulative number of test cases that failed.
 * 		setupFailures(int): The cumulative number of tests cases that generated an error prior
 * 					to actually executing the given test case.
 * 		errors(int): The cumulative number of test cases that generated an error.  Test cases that
 * 					generated an error prior to executing the given test are included in
 * 					this count.  Test cases that generated an error while executing the
 * 					given test are included as well.
 * 		unverified(int): The cumulative number of invalid records encountered while processing
 * 					the set of test data.
 * 		testingSuiteError(int): Return 1 if unable to execute one or more testing functions.
 * 					Return 0 otherwise.
 *
 */
TestSuite testArcTanToTwoLoci_AllSets()
{

	TestSuite suite;
	TestSet set1;

    printf("\nStart testArcTanToTwoLoci_AllSets\n");

    suite = newTestSuite("testArcTanToTwoLoci_AllSets");

    set1 = testArcTanToTwoLoci_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("Finish testArcTanToTwoLoci_AllSets\n\n\n");

    return suite;
}

/*
 * NAME: testLocusCrsAtPt_Set1
 *
 * DESCRIPTION:
 * 		This function is used to test the locusCrsAtPt function.
 *
 * 		This function runs the test data created by Balakrishna Babu and uses
 * 		his original test plan to execute the actual tests.
 *
 *		Approach for testing locusCrsAtPt is as follows:
 *
 *	1)	Use random number to generate input data for locus.  Input data consists of
 *	    locus (start and end points of the defining geodesic, start and end distances
 *	    that define the locus, and the lineType).  Select latitude and longitude for
 *	    geoStart, and select azimuth and distance to compute geoEnd.
 *	2)	Create NewLocus using input data. Compute the slope of locus
 *	    ((endDist-StartDist)/Length of geodesic).
 *	3)	Select test points using the following two steps (4,5).
 *	4)	Select a distance (dAlong, nm values of -1000, -5, 0, 5, 1000, 0.5*geoLength,
 *	    geoLength) and compute geoBase along the defining geodesic.
 *	5)	Using the geoBase, compute the point on the locus (ptOnLocusFromGeoPt).
 *	    Save ptonloc and perpCrs.  ptonloc will be the testPt.
 *	    Compute and save the course from testPt to geoBase.
 *	6)	Call the function under test with testPoints.
 *	7)	Verify that the geoPt lies on the geodesic.
 *	8)	Verify that the course of locus equals course from testPt to geoPt plus
 *		90 degrees plus locus slope.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testLocusCrsAtPt_Set1(TestSet) - A test set with the folling metrics:
 *  	testCases(int): The number of test cases executed.
 * 		pass(int): The number of test cases that passed.
 * 		fail(int): The number of test cases that failed.
 * 		setupFailures(int): The number of tests cases that generated an error prior
 * 					to actually executing the given test case.
 * 		errors(int): The number of test cases that generated an error.  Test cases that
 * 					generated an error prior to executing the given test are included in
 * 					this count.  Test cases that generated an error while executing the
 * 					given test are included as well.
 * 		unverified(int): The number of invalid records encountered while processing
 * 					the set of test data.
 * 		testingError(int): Return 1 if unable to execute testing function.  Return 0 otherwise.
 *
 */
//TODO separate the test data creation logic from the actual testing logic
TestSet testLocusCrsAtPt_Set1()
{

    int testNum=0;
    char testname[256];
    double geoLatS, geoLonS, startDist, endDist;
    double geoCrs, geoLength;
    double locSlope, locSlopeAngle, locCrs, locCrsExp;
    double dAlong[7] = {-1000.0, -5.0, 0.0, 5.0, 1000.0, 0.0, 0.0};
    double maxDistFromLoc = 1000.0;
    LineType lineType;  //0=SEGMENT,1=SEMIINFINITE,2=INFINITE
    double DEG2RAD = M_PI / 180.0;
    double AZDEGTOL = 5.0e-7; //5e-7 deg or 8.72665e-9 rad or 0.162 cm at 100 nm

    LLPoint geoStart, geoEnd, geoBase, geoPt, ptonloc;
    LLPoint testPt;
    Locus locus;
    double perpCrs, crsToBase, crsFromBase, distToLoc;
    int i;
    int passedCount=0, failedCount=0;
    int testCaseCount = 0, setupFailureCount = 0, errorCount = 0, unverifiedCount = 0;
    ErrorSet err=0;

    TestSet set;

    long newSeed=20080520;

    printf("Start testLocusCrsAtPt_Set1\n");

    set = newTestSet("testLocusCrsAtPt_Set1");

	srand(newSeed);  //Initialize the random number generator

	while (testNum<1000)
    {
		//Initialize or reset to defaults
    	err = 0;
    	perpCrs = 0.0;
        locCrs = 0.0;

		//Select start point, azimuth and distance to end point for the defining geodesic
    	//Select start and end distances for defining locus start and end points
		geoLatS = randLat();
		geoLonS = randLon();
		geoCrs = DEG2RAD * randAzimuth();
    	geoLength = 0.1 * randDist();  //about 540 nm
    	dAlong[5] = 0.5 * geoLength;
    	dAlong[6] = geoLength;
    	startDist = 0.1 * randSignedDist();  //about -540 to +540 nm
    	endDist = 0.1 * randSignedDist();  //about -540 to +540 nm
    	locSlope = (endDist-startDist)/geoLength;
    	locSlopeAngle = atan(locSlope);
        if (fabs(locSlopeAngle) > 89.0 * DEG2RAD)
          continue;
    	distToLoc = startDist - 1000.0 * locSlope;
    	if(fabs(distToLoc) > maxDistFromLoc){
    		dAlong[0] = (locSlope > 0) ? (-maxDistFromLoc - startDist)/locSlope : (maxDistFromLoc - startDist)/locSlope;
    	} else {
    		dAlong[0] = -1000.0;
    	}
    	distToLoc = startDist + 1000.0 * locSlope;
    	if(fabs(distToLoc) > maxDistFromLoc){
    		dAlong[4] = (locSlope > 0) ? (maxDistFromLoc - startDist)/locSlope : (-maxDistFromLoc - startDist)/locSlope;
		} else {
			dAlong[4] = 1000.0;
		}
    	testNum++;

        sprintf(testname,"TEST%-d",testNum);

        geoStart.latitude = geoLatS * DEG2RAD;
        geoStart.longitude = geoLonS * DEG2RAD;

        //Compute end point
        err |= direct(geoStart, geoCrs, geoLength, &geoEnd, EPS);


        for (int lineTypeI = 0; lineTypeI < 3; lineTypeI++) {
            lineType = static_cast<LineType>(lineTypeI);

            err |= createLocus(&locus, geoStart, geoEnd, startDist, endDist, lineType, TOL,EPS);

	        //For each point along geodesic
	        for (i=0; i<7; i++) {

	        	err |= direct(geoStart,geoCrs,dAlong[i],&geoBase,EPS);
	        	err |= ptOnLocusFromGeoPt(locus,geoBase,&ptonloc,&perpCrs,TOL,EPS); //ignore perpCrs
	        	testPt = ptonloc;

	            crsToBase = 0.0;
	            crsFromBase = 0.0;
	        	err |= invCrs(testPt,geoBase,&crsToBase,&crsFromBase,EPS);
	        	if(ptsAreSame(geoBase,testPt,TOL)){
	        		invCrs(geoStart,geoBase,NULL,&crsToBase,EPS);
                                if (i == 2)
                                  crsToBase = geoCrs;
	        		locCrsExp = crsToBase + locSlopeAngle;
				} else {
					locCrsExp = crsToBase + M_PI_2 + locSlopeAngle;
				}

	        	if (locCrsExp >= M_2PI) {
	        		locCrsExp -= M_2PI;
	        	}
	        	if (locCrsExp < 0.0) {
	        		locCrsExp += M_2PI;
	        	}

	            if ( ( getMaskedError(err, getMaskAll()) ) || (locCrsExp == -1.0) ) {
		        	printf("Error occurred in pre-locusCrsAtPt err=0x%lx lineTyep=%i i=%i locCrs=%14.8f\n",
		        			err,lineType,i,locCrs);
                    setupFailureCount++;
                    errorCount++;
                    testCaseCount++;
                    err = 0; // Clear the error
		        	continue;
		        }

	            locCrs = locusCrsAtPt(locus,testPt,&geoPt,&perpCrs,&err,TOL,EPS);

	            if ( ( getMaskedError(err, getMaskAll()) ) || (locCrs == -1.0) ) {
		        	printf("Error occurred in locusCrsAtPt err=0x%lx testNum %d lineType=%i i=%i locCrs=%14.8f\n",
		        			err,testNum,lineType,i,locCrs);
		        	failedCount++;
                    errorCount++;
                    testCaseCount++;
                    err = 0; // Clear the error

		        	continue;
		        }

		        if ( (fabs(locCrs-locCrsExp) <= AZDEGTOL*DEG2RAD) ||
		        	 (fabs(fabs(locCrs-locCrsExp)-M_PI) <= AZDEGTOL*DEG2RAD) ) {
		        	if (PRINT_PASSED_CASES) printLLPoint("   geoStart:  ",geoStart);
		        	if (PRINT_PASSED_CASES) printLLPoint("   geoEnd:    ",geoEnd);
		        	if (PRINT_PASSED_CASES) printLLPoint("   locusStart:",locus.locusStart);
		        	if (PRINT_PASSED_CASES) printLLPoint("   locusEnd:  ",locus.locusEnd);
		        	if (PRINT_PASSED_CASES) printLLPoint("   geoBase:   ",geoBase);
		        	if (PRINT_PASSED_CASES) printLLPoint("   testPt:   ",testPt);
		        	if (PRINT_PASSED_CASES) printf("%s geoCrs=%14.8f geoLength=%14.8f startDist=%14.8f endDist=%14.8f locSlopeAngle=%14.8f\n",
		            		testname, geoCrs/DEG2RAD, geoLength, startDist, endDist, locSlopeAngle/DEG2RAD);

		        	if (PRINT_PASSED_CASES) printf("%s passed lineType=%d  dAlong=%14.8f locCrsExp=%14.8f locCrs=%14.8f\n",
		        			testname,lineType,dAlong[i],locCrsExp/DEG2RAD,locCrs/DEG2RAD);
		        	passedCount++;
                    testCaseCount++;
		        }
		        else {
		            printLLPoint("   geoStart:  ",geoStart);
		            printLLPoint("   geoEnd:    ",geoEnd);
		            printLLPoint("   locusStart:",locus.locusStart);
		            printLLPoint("   locusEnd:  ",locus.locusEnd);
		            printLLPoint("   geoBase:   ",geoBase);
		            printLLPoint("   testPt:   ",testPt);
		            printLLPoint("   geoPt:    ",geoPt);
		            printf("crsToBase:%3.20lf, perpCrs:%3.20lf\n", crsToBase, perpCrs);
		            printf("%s geoCrs=%14.8f geoLength=%14.8f startDist=%14.8f endDist=%14.8f locSlopeAngle=%14.8f\n",
		            		testname, geoCrs/DEG2RAD, geoLength, startDist, endDist, locSlopeAngle/DEG2RAD);

		        	printf("%s failed lineType=%d  dAlong=%14.8f locCrsExp=%14.8f locCrs=%14.8f\n",
		        			testname,lineType,dAlong[i],locCrsExp/DEG2RAD,locCrs/DEG2RAD);
		        	failedCount++;
                    testCaseCount++;
		        }
	        }
        }
    }
    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("\nFinish testLocusCrsAtPt_Set1\n");
    return set;
}

/*
 * NAME: testLocusCrsAtPt_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the locusCrsAtPt function.
 *
 * 		This function runs all the test cases from all the test sets and returns
 * 		the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testLocusCrsAtPt_AllSets(TestSuite) - A test suite with the folling metrics:
 *  	testCases(int): The cumulative number of test cases executed.
 * 		pass(int): The cumulative number of test cases that passed.
 * 		fail(int): The cumulative number of test cases that failed.
 * 		setupFailures(int): The cumulative number of tests cases that generated an error prior
 * 					to actually executing the given test case.
 * 		errors(int): The cumulative number of test cases that generated an error.  Test cases that
 * 					generated an error prior to executing the given test are included in
 * 					this count.  Test cases that generated an error while executing the
 * 					given test are included as well.
 * 		unverified(int): The cumulative number of invalid records encountered while processing
 * 					the set of test data.
 * 		testingSuiteError(int): Return 1 if unable to execute one or more testing functions.
 * 					Return 0 otherwise.
 *
 */
TestSuite testLocusCrsAtPt_AllSets()
{

	TestSuite suite;
	TestSet set1;

    printf("\nStart testLocusCrsAtPt_AllSets\n");

    suite = newTestSuite("testLocusCrsAtPt_AllSets");

    set1 = testLocusCrsAtPt_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("Finish testLocusCrsAtPt_AllSets\n\n\n");

    return suite;
}

/*
 * NAME: testProjectToLocus_Set1
 *
 * DESCRIPTION:
 * 		This function is used to test the projectToLocus function.
 *
 * 		This function runs the test data created by Balakrishna Babu and uses
 * 		his original test plan to execute the actual tests.
 *
 *		Approach for testing projectToLocus is as follows:
 *
 *	1)  Use random number to select the expected projection point.
 *  2)	Construct a locus that passes through the selected projection point.
 * 		a.	Azimuth and a distance to find the point on the defining geodesic
 * 		b.	Compute azimuth to the point on locus from point on geodesic
 * 		c.	Compute two azimuths 90 degress on either side of azimuth to point
 *      	and select two distances to locate the start and end points of the
 *      	defining geodesic
 * 		d.	Select startDist and compute endDist and slope of locus
 *            ((endDist-StartDist)/Length of geodesic)
 * 		e.	Create locus using start and end points of defining geodesic and
 *          start and end distances for the locus
 *  3)	Compute course of locus at the expected projection point (azimuth from
 *      projection point to point on geodesic + 90 degrees + locus slope).
 * 4)	Compute the course to the testPt from the expected projection point by
 *      adding 90 degrees to the course computed in previous step, select a
 *      distance and compute the testPt.
 * 5)	Call the function under test with testPt.
 * 6)	Verify that the resulting projection point is the same as the expected
 *      projection point.  Also verify the course and distance to the projection point.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testProjectToLocus_Set1(TestSet) - A test set with the folling metrics:
 *  	testCases(int): The number of test cases executed.
 * 		pass(int): The number of test cases that passed.
 * 		fail(int): The number of test cases that failed.
 * 		setupFailures(int): The number of tests cases that generated an error prior
 * 					to actually executing the given test case.
 * 		errors(int): The number of test cases that generated an error.  Test cases that
 * 					generated an error prior to executing the given test are included in
 * 					this count.  Test cases that generated an error while executing the
 * 					given test are included as well.
 * 		unverified(int): The number of invalid records encountered while processing
 * 					the set of test data.
 * 		testingError(int): Return 1 if unable to execute testing function.  Return 0 otherwise.
 *
 */
//TODO separate the test data creation logic from the actual testing logic
TestSet testProjectToLocus_Set1() {

	double DEG2RAD = M_PI / 180.0;

    int testNum=0;
    char testname[256];

    double latP, lonP;
    double az1, d1, d2, d3, slope1, projDist;
    double startDist1, endDist1;
    double crsC1, crsC2, crsC3, crsC4, slopeAngle, locusCrs;

    LLPoint projExp, proj, proj1, proj2, testPt;
    LLPoint geoStart1, geoEnd1;//, locusStart1, locusEnd1;
    LineType lineType1;  //0=SEGMENT,1=SEMIINFINITE,2=INFINITE
    Locus locus1;
    Locus* locusPtr1;
    double crsFromPt, distFromPt;

    double fcrs,bcrs;
    int expectValidIntx1;

    int passedCount=0, failedCount=0;
    int errorCount = 0;
    int unverifiedCount = 0;
    int setupFailureCount = 0;
    int testCaseCount = 0;
    ErrorSet err=0;

    TestSet set;

    long newSeed=20080516;

    printf("Start testProjectToLocus_Set1\n");

    set = newTestSet("testProjectToLocus_Set1");

	srand(newSeed);  //Initialize the random number generator

	locusPtr1 = &locus1;

	while (testNum < 1000)
    {
    	err = 0;

    	testNum++;
        sprintf(testname,"TEST%-d",testNum);

        //Select an intersection point
    	latP = randLat();
    	lonP = randLon();
        projExp.latitude = latP * DEG2RAD;
        projExp.longitude = lonP * DEG2RAD;

    	//Construct a locus that passes through the selected projection point
        //1. Azimuth and a distance to find the point on the defining geodesic
    	az1 = randAzimuth();
    	d1 = 0.01*randDist();
        crsC1 = az1 * DEG2RAD;
        err |= direct(projExp, crsC1, d1, &proj1, EPS);

        //2. Compute azimuth to the point on locus from point on geodesic
        err |= invCrs(proj1,projExp,&fcrs,&bcrs,EPS);

        //3. Compute two azimuths and select two distances to locate the start and
        //   end points of the defining geodesic
        crsC2 = fcrs + M_PI_2;
        if (crsC2 > M_2PI) {
        	crsC2 -= M_2PI;
        }
        crsC3 = fcrs - M_PI_2;
        if (crsC3 < 0.0) {
        	crsC3 += M_2PI;
        }
    	d2 = 0.1*randDist();
    	d3 = 0.1*randDist();
        err |= direct(proj1, crsC2, d2, &geoStart1, EPS);
        err |= direct(proj1, crsC3, d3, &geoEnd1, EPS);

        //4. Select startDist and compute endDist
    	startDist1 = 0.01*randDist();
    	slope1 = (d1-startDist1)/d2;
        // slope angle (in radians)
        slopeAngle = atan(slope1);
    	endDist1 = d1 + d3 * (d1-startDist1)/d2;
    	lineType1 = INFINITE;

    	//5. Create locus that passed through the selected intersection point
        err |= createLocus(locusPtr1,geoStart1,geoEnd1,startDist1,endDist1,lineType1,TOL,EPS);

        //Compute locus course at the expected projection point
        locusCrs = crsC1 + M_PI_2 + slopeAngle;
        if (locusCrs >= M_2PI) {
        	locusCrs -= M_2PI;
        }
        if (locusCrs < 0.0) {
        	locusCrs += M_2PI;
        }

    	projDist = 0.1*randDist();
    	crsC4 = locusCrs + M_PI_2;
        err |= direct(projExp, crsC4, projDist, &testPt, EPS);

    	if (ptIsOnLocus(locus1,projExp,&proj2,&err,TOL,EPS)) {
    		expectValidIntx1 = 1;  //true
    	}
    	else {
    		expectValidIntx1 = 0;  //false
    	}

        if ( getMaskedError(err, getMaskAll()) )
        {
        	printf("Error occurred in pre-projectToLocus INFINITE err=0x%lx\n",err);
        	setupFailureCount++;
        	testCaseCount++;
        	errorCount++;
        	continue;
        }

    	proj.latitude = 0.0;
        proj.longitude = 0.0;
        crsFromPt = 0.0;
        distFromPt = 0.0;
        err |= projectToLocus(locus1,testPt,&proj,&crsFromPt,&distFromPt,TOL,EPS);

        if ( getMaskedError(err, getMaskAll()) )
        {
        	printf("Error occurred in projectToLocus INFINITE err=0x%lx\n",err);
        	errorCount++;
        	failedCount++;
        	testCaseCount++;
        	continue;
        }

    	if (ptIsOnLocus(locus1,proj,&proj2,&err,TOL,EPS)) {
    		expectValidIntx1 = 1;  //true
    	}
    	else {
    		expectValidIntx1 = 0;  //false
    	}

        if ( getMaskedError(err, getMaskAll()) )
        {
        	printf("Error occurred in post-projectToLocus INFINITE err=0x%lx\n",err);
        	failedCount++;
        	errorCount++;
        	testCaseCount++;
        	continue;
        }

        if ( ptsAreSame(projExp,proj,TESTTOL) ) {
        	if (PRINT_PASSED_CASES) printf("%s latP=%14.8f lonP=%14.8f\n",
              		testname,latP,lonP);
        	if (PRINT_PASSED_CASES) printLLPoint("Expected Projection point",projExp);
        	if (PRINT_PASSED_CASES) printLLPoint("point on defining geodesic",proj1);
        	if (PRINT_PASSED_CASES) printLLPoint("geoStart1",geoStart1);
        	if (PRINT_PASSED_CASES) printLLPoint("geoEnd1",geoEnd1);
        	if (PRINT_PASSED_CASES) printf("%s X: %14.8f %14.8f  az1=%14.8f d1=%14.8f d2=%14.8f d3=%14.8f startDist1=%14.8f endDist1=%14.8f slope1=%14.8f\n",
              		testname,latP,lonP,az1,d1,d2,d3,startDist1,endDist1,slope1);
        	if (PRINT_PASSED_CASES) printf("%s slopeAngle=%14.8f crsC1=%14.8f crsC2=%14.8f crsC3=%14.8f\n",
        			testname,slopeAngle/DEG2RAD,crsC1/DEG2RAD,crsC2/DEG2RAD,crsC3/DEG2RAD);
        	if (PRINT_PASSED_CASES) printLLPoint("locusStart1",locus1.locusStart);
        	if (PRINT_PASSED_CASES) printLLPoint("locusEnd1",locus1.locusEnd);
        	if (PRINT_PASSED_CASES) printf("%s projDist=%14.8f crsC4=%14.8f\n",
        			testname,projDist,crsC4/DEG2RAD);
        	if (PRINT_PASSED_CASES)  printLLPoint("testPt",testPt);
        	if (PRINT_PASSED_CASES) printf("%s lineType1=%i expIntx1=%i\n",
        			testname,locus1.lineType,expectValidIntx1);
        	if (PRINT_PASSED_CASES) printLLPoint("proj",proj);

        	if (PRINT_PASSED_CASES) printf("%s passed expIntx1=%i\n",
        			testname,expectValidIntx1);
        	passedCount++;
    	}
    	else {
            printf("%s latP=%14.8f lonP=%14.8f\n",
              		testname,latP,lonP);
            printLLPoint("Expected Projection point",projExp);
            printLLPoint("point on defining geodesic",proj1);
            printLLPoint("geoStart1",geoStart1);
            printLLPoint("geoEnd1",geoEnd1);
            printf("%s X: %14.8f %14.8f  az1=%14.8f d1=%14.8f d2=%14.8f d3=%14.8f startDist1=%14.8f endDist1=%14.8f slope1=%14.8f\n",
              		testname,latP,lonP,az1,d1,d2,d3,startDist1,endDist1,slope1);
        	printf("%s slopeAngle=%14.8f crsC1=%14.8f crsC2=%14.8f crsC3=%14.8f\n",
        			testname,slopeAngle/DEG2RAD,crsC1/DEG2RAD,crsC2/DEG2RAD,crsC3/DEG2RAD);
            printLLPoint("locusStart1",locus1.locusStart);
            printLLPoint("locusEnd1",locus1.locusEnd);
        	printf("%s projDist=%14.8f crsC4=%14.8f\n",
        			testname,projDist,crsC4/DEG2RAD);
            printLLPoint("testPt",testPt);
        	printf("%s lineType1=%i expIntx1=%i\n",
        			testname,locus1.lineType,expectValidIntx1);
            printLLPoint("proj",proj);

        	printf("%s failed expIntx1=%i\n",
        			testname,expectValidIntx1);
        	failedCount++;
    	}

    	testCaseCount++;
    }

    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("\nFinish testProjectToLocus_Set1\n\n\n");

    return set;
}

/*
 * NAME: testProjectToLocus_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the projectToLocus function.
 *
 * 		This function runs all the test cases from all the test sets for testing
 *      the above function and returns the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testProjectToLocus_AllSets(TestSuite) - A test suite with the folling metrics:
 *  	testCases(int): The cumulative number of test cases executed.
 * 		pass(int): The cumulative number of test cases that passed.
 * 		fail(int): The cumulative number of test cases that failed.
 * 		setupFailures(int): The cumulative number of tests cases that generated an error prior
 * 					to actually executing the given test case.
 * 		errors(int): The cumulative number of test cases that generated an error.  Test cases that
 * 					generated an error prior to executing the given test are included in
 * 					this count.  Test cases that generated an error while executing the
 * 					given test are included as well.
 * 		unverified(int): The cumulative number of invalid records encountered while processing
 * 					the set of test data.
 * 		testingSuiteError(int): Return 1 if unable to execute one or more testing functions.
 * 					Return 0 otherwise.
 *
 */
TestSuite testProjectToLocus_AllSets()
{
	TestSuite suite;
	TestSet set1;

    printf("\nStart testProjectToLocus_AllSets\n");

    suite = newTestSuite("testProjectToLocus_AllSets");

    set1 = testProjectToLocus_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("\nFinish testProjectToLocus_AllSets\n\n\n");

    return suite;
}

/*
 * NAME: testLociCoincide_Set1
 *
 * DESCRIPTION:
 * 		This function is used to test the lociCoincide function for the
 * 		cases created by Rich Snow.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testLociCoincide_Set1(TestSet) - A test set with the folling metrics:
 *  	testCases(int): The number of test cases executed.
 * 		pass(int): The number of test cases that passed.
 * 		fail(int): The number of test cases that failed.
 * 		setupFailures(int): The number of tests cases that generated an error prior
 * 					to actually executing the given test case.
 * 		errors(int): The number of test cases that generated an error.  Test cases that
 * 					generated an error prior to executing the given test are included in
 * 					this count.  Test cases that generated an error while executing the
 * 					given test are included as well.
 * 		unverified(int): The number of invalid records encountered while processing
 * 					the set of test data.
 * 		testingError(int): Return 1 if unable to execute testing function.  Return 0 otherwise.
 *
 */
TestSet testLociCoincide_Set1()
{

    double DEG2RAD = M_PI / 180.0;

    ShapeType expShapeType;
    double expStartDist, expEndDist;
    double geo1LatS, geo1LonS, geo1Az, geo1len, geo2len, dir1st, dir1end; 
    double start0Dist, end0Dist, start2Dist, end2Dist, dist, tempdbl;
    LineType expLineType, lineType;
    Locus expLocus;

    LLPoint geo1Start, geo1End, geo2Start, geo2End, geoPt;
    Locus locus1, locus2;
    Locus* locus1Ptr;
    Locus* locus2Ptr;


    LLPoint actIntX;

    LLPoint expPt1, expPt2;

    LLPoint dummyPt;

    int coincide = 0;
    Shape common;
    LLPoint commonPt;
    Locus commonLocus;



    int passedCount = 0, failedCount = 0, unverifiedCount = 0;
    int setupFailureCount = 0, errorCount = 0, testCaseCount = 0;
    ErrorSet err = 0;
    int ii, jj;
    int failed = 0;


    TestSet set;

    printf("Start testLociCoincide_Set1\n");

    set = newTestSet("testLociCoincide_Set1");

    srand(05252011);

    for (ii = 0; ii < 100; ii++)
    {
      err = 0;
      geo1LatS = randLat();
      geo1LonS = randLon();
      geo1Start.latitude = geo1LatS * DEG2RAD;
      geo1Start.longitude = geo1LonS * DEG2RAD;
      geo1Az = DEG2RAD * randAzimuth();
      geo1len = 100.0 + 0.01 * randDist();
      err |= direct(geo1Start, geo1Az, geo1len, &geo1End, EPS);
      dir1st = 1.0;
      if ((rand() % 2) == 0)
        dir1st = -1.0;
      dir1end = 1.0;
      if ((rand() % 2) == 0)
        dir1end = -1.0;
      start0Dist = dir1st * 0.01 * randDist();
      end0Dist = dir1end * 0.01 * randDist();
      locus1Ptr = &locus1;
      locus2Ptr = &locus2;
      lineType = SEGMENT;
      err |= createLocus(locus1Ptr, geo1Start, geo1End, start0Dist, end0Dist, lineType, TOL, EPS);
      //printf("ii = %d geo1LatS %4.5f geo1LonS %4.5f geo1Az %4.5f geo1len %4.5f start0Dist %4.5f end0Dist %4.5f\n",ii,geo1LatS,geo1LonS,geo1Az*180.0/M_PI,geo1len,start0Dist,end0Dist);
      for (jj = 0; jj < 9; jj++)
      {
        err = 0;
        if (jj == 0)
        {
          // identical loci
          err |= createLocus(locus2Ptr, geo1Start, geo1End, start0Dist, end0Dist, lineType, TOL, EPS);
          expPt1 = geo1Start;
          expPt2 = geo1End;
          expStartDist = start0Dist;
          expEndDist = end0Dist;
          expShapeType = LOCUS;
          expLineType = SEGMENT; 
        }
        else if (jj == 1)
        {
          //common start and locus2 lies on locus1
          geo2len = geo1len / 2.0;
          err |= direct(geo1Start, geo1Az, geo2len, &geo2End, EPS);
          end2Dist = distToLocusFromGeoDist(*locus1Ptr, geo2len);
          err |= createLocus(locus2Ptr, geo1Start, geo2End, start0Dist, end2Dist, lineType, TOL, EPS);
          expPt1 = geo1Start;
          expPt2 = geo2End;
          expStartDist = start0Dist;
          expEndDist = end2Dist;
          expShapeType = LOCUS;
          expLineType = SEGMENT; 
        }
        else if (jj == 2)
        {
          //no common start/end points but some part of loci is common
          err |= direct(geo1Start, geo1Az, 0.25*geo1len, &geo2Start, EPS);
          err |= direct(geo1Start, geo1Az, 0.75*geo1len, &geo2End, EPS);
          start2Dist = distToLocusFromGeoDist(*locus1Ptr, 0.25*geo1len);
          end2Dist = distToLocusFromGeoDist(*locus1Ptr, 0.75*geo1len);
          err |= createLocus(locus2Ptr, geo2Start, geo2End, start2Dist, end2Dist, lineType, TOL, EPS);
          expPt1 = geo2Start;
          expPt2 = geo2End;
          expStartDist = start2Dist;
          expEndDist = end2Dist;
          expShapeType = LOCUS;
          expLineType = SEGMENT; 
        }
        else if (jj == 3 || jj == 4)
        {
          //one start/end point is all that is common
          if (jj == 3) //geo directions the same
            err |= createLocus(locus2Ptr, geo1Start, geo1End, start0Dist, end0Dist / 2.0, lineType, TOL, EPS);
          else
          {
            //geo directions opposite 
            err |= direct(geo1Start, geo1Az + M_PI, geo1len, &geo2End, EPS);
            err |= createLocus(locus2Ptr, geo1Start, geo2End, -start0Dist, -end0Dist / 2.0, lineType, TOL, EPS);
          }
          expPt1 = locus1Ptr->locusStart;
          expShapeType = LLPOINT; 
        }
        else if (jj == 5 || jj == 6)
        {
          //no start/end point is common but there is an intersection
          if (jj == 5)
          {
            //geo directions the same
            err |= createLocus(locus2Ptr, geo1Start, geo1End, end0Dist, start0Dist, lineType, TOL, EPS);
            dist = (locus2Ptr->startDist - locus1Ptr->startDist) / (locus1Ptr->slope - locus2Ptr->slope);
          }
          else
          {
            //geo directions opposite 
            err |= createLocus(locus2Ptr, geo1End, geo1Start, -start0Dist, -end0Dist, lineType, TOL, EPS);
            dist = (-locus2Ptr->endDist - locus1Ptr->startDist) / (locus1Ptr->slope - locus2Ptr->slope);
          }
          err |= direct(geo1Start, geo1Az, dist, &geoPt, EPS);
          err |= ptOnLocusFromGeoPt(locus1, geoPt, &expPt1, &tempdbl, TOL, EPS);
          expShapeType = LLPOINT; 
        }
        else if (jj == 7)
        {
          //colinear geos but no intersection
          if ((start0Dist > 0.0 && end0Dist > 0.0) || (start0Dist < 0.0 && end0Dist < 0.0))
            err |= createLocus(locus2Ptr, geo1Start, geo1End, start0Dist/2.0, end0Dist/2.0, lineType, TOL, EPS);
          else if (start0Dist > 0.0 && end0Dist < 0.0)
            err |= createLocus(locus2Ptr, geo1Start, geo1End, start0Dist/2.0, 1.5*end0Dist, lineType, TOL, EPS);
          else if (start0Dist < 0.0 && end0Dist > 0.0)
            err |= createLocus(locus2Ptr, geo1Start, geo1End, 2.0*start0Dist, end0Dist/2.0, lineType, TOL, EPS);
          expShapeType = ARC;//to go through default expShapeType
        } 
        else if (jj == 8)
        {
          //geos not colinear
          err |= direct(geo1Start, geo1Az + M_PI/6.0, geo1len, &geo1End, EPS);
          err |= createLocus(locus2Ptr, geo1Start, geo1End, start0Dist, end0Dist, lineType, TOL, EPS);
          expShapeType = ARC;//to go through default expShapeType
        }

        if (expShapeType == LOCUS)
        	err |= createLocus(&expLocus, expPt1, expPt2, expStartDist, expEndDist, expLineType, TOL, EPS);


        if (getMaskedError(err, getMaskAll()))
        {
            printf("ii = %d jj = %d - Error occurred in setup err=0x%lx\n", ii, jj, err);
            setupFailureCount++;
            errorCount++;
            testCaseCount++;
            continue;
        }

        coincide = lociCoincide(locus1, locus2, &common, &err, TOL, EPS);

        failed = 0;

        if (getMaskedError(err, getMaskAll()))
        {
            printf("ii = %d jj = %d - Error occurred in lociCoincide err=0x%lx\n", ii, jj, err);
            failedCount++;
            errorCount++;
            testCaseCount++;
            continue;
        }

        switch(expShapeType)
        {
	  case LLPOINT:
 	    if (!coincide)
	      failed = 1;

	    if (common.type != LLPOINT)
	      failed = 1;
	    else
	    {
	      commonPt = *((LLPoint*) &common);
	      //fail if expected and actual intersection points are not the same, with respect to tol
	      if ( !ptsAreSame(expPt1, commonPt, TESTTOL) )
	        failed = 1;
	      //fail if the actual intersection does not lie on both loci
	      if ( (!ptIsOnLocus(locus1, commonPt, &dummyPt, &err, TESTTOL, EPS)) && (!ptIsOnLocus(locus2, actIntX, &dummyPt, &err, TESTTOL, EPS)) )
		failed = 1;
	      //fail if any error occurred because we are expecting a valid intersection in this case
	      if ( getMaskedError(err, getMaskAll()) )
		failed = 1;
	    }
	    break;
	  case LOCUS:
	    if (!coincide)
	      failed = 1;

	    if (common.type != LOCUS)
	      failed = 1;
	    else
	    {
	      commonLocus = *((Locus*) &common);

	      //orient the loci so that they are in the same direction
	      if (commonLocus.geoAz - expLocus.geoAz > M_PI_4) // locus 1 and locus 2 point in opposite directions
	      {
		err |= createLocus(&commonLocus, expPt2, expPt1, expEndDist, expStartDist, expLineType, TOL, EPS);
	      }

	      if ( !ptsAreSame(expLocus.geoStart, commonLocus.geoStart, TESTTOL))
		failed = 1;

	      if ( !ptsAreSame(expLocus.geoEnd, commonLocus.geoEnd, TESTTOL))
		failed = 1;

	      if ( fabs(expLocus.startDist - commonLocus.startDist) >= TESTTOL)
		failed = 1;

	      if ( fabs(expLocus.endDist - commonLocus.endDist) >= TESTTOL)
		failed = 1;

	      if (expLocus.lineType != commonLocus.lineType)
		failed = 1;

	    }
	    break;
	    default://no intersections at all
	      if (coincide)
	       failed = 1;
        }

        if (failed)
        {
	  printf("\nii = %d jj = %d failed\n", ii, jj);
	  failedCount++;
        } else {
          passedCount++;
        }

        testCaseCount++;
      }//for jj

    }//for ii

    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("Finish testLociCoincide_Set1\n");

    return set;

}

/*
 * NAME: testLociCoincide_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the lociCoincide function.
 *
 * 		This function runs all the test cases from all the test sets for testing
 *      the above function and returns the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testLociCoincide_AllSets(TestSuite) - A test suite with the following metrics:
 *  	testCases(int): The cumulative number of test cases executed.
 * 		pass(int): The cumulative number of test cases that passed.
 * 		fail(int): The cumulative number of test cases that failed.
 * 		setupFailures(int): The cumulative number of tests cases that generated an error prior
 * 					to actually executing the given test case.
 * 		errors(int): The cumulative number of test cases that generated an error.  Test cases that
 * 					generated an error prior to executing the given test are included in
 * 					this count.  Test cases that generated an error while executing the
 * 					given test are included as well.
 * 		unverified(int): The cumulative number of invalid records encountered while processing
 * 					the set of test data.
 * 		testingSuiteError(int): Return 1 if unable to execute one or more testing functions.
 * 					Return 0 otherwise.
 *
 */
TestSuite testLociCoincide_AllSets()
{
	TestSuite suite;
        TestSet set1;

    printf("\nStart testLociCoincide_AllSets\n");

    suite = newTestSuite("testLociCoincide_AllSets");

    set1 = testLociCoincide_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("\nFinish testLociCoincide_AllSets\n\n\n");

    return suite;
}

void testBug33911() {

	LLPoint geoS, geoE1, geoE2, intx, geoPt, testPt;
	Locus loc1, loc2;
	ErrorSet err = 0;
	double perpCrs, dist;

	double tol = 1.37e-9;
	double eps = 1e-20;

	err |= createPt(&geoS, 35.833645437300056 * M_PI / 180, -115.5713858704401 * M_PI / 180);
	err |= createPt(&geoE1, 36.0005555555556 * M_PI / 180, -115.562222222222 * M_PI / 180);
	err |= createPt(&geoE2, 35.89143551105784 * M_PI / 180, -115.56821748598247 * M_PI / 180);

	err |= createLocus(&loc1, geoS, geoE1, -6, -6, LineType::SEGMENT, tol, eps);
	err |= createLocus(&loc2, geoS, geoE2, -6, -3, LineType::SEGMENT, tol, eps);

	err |= locusIntx(loc1, loc2, &intx, tol, eps);

	printf("Start Point 1:  %f %f\n", loc1.locusStart.latitude * 180 / M_PI, loc1.locusStart.longitude * 180 / M_PI);

	printf("Start Point 2:  %f %f\n", loc2.locusStart.latitude * 180 / M_PI, loc2.locusStart.longitude * 180 / M_PI);

	printf("Intersection Point:  %f %f\n", intx.latitude * 180 / M_PI, intx.longitude * 180 / M_PI);

	ptIsOnLocus(loc1, intx, &geoPt, &err, tol, eps);

	err |= ptOnLocusFromGeoPt(loc1, geoPt, &testPt, &perpCrs, tol, eps);

	printf("Intersection Point:  %f %f\n", testPt.latitude * 180 / M_PI, testPt.longitude * 180 / M_PI);

	err |= invDist(testPt, intx, &dist, eps);

	printf("Distance:  %e\n", dist);

	int size = 1536;

	char t[size];

	aboutGeolib(t, size);

	printf("%s", t);

	//printf("%s, %u\n", t, sizeof(t));
}

void testBug33966() {

	LLPoint geoS1, geoS2, geoE1, geoE2, cp, sp, ep;
	ArcDirection dir;

	Locus pLoc1, pLoc2, sLoc1, sLoc2;

	double radius = 20112463720815415;

	double crs1, bcrs1, crs2, bcrs2;

	double tol = 1.37e-9;
	double eps = 1e-20;

	ErrorSet err = 0;

	err |= createPt(&geoS1, 36.065133948365904 * M_PI / 180, -114.59919492451887 * M_PI / 180);
	err |= createPt(&geoS2, 36.12778046412819 * M_PI / 180, -114.59241342956928 * M_PI / 180);

	err |= createPt(&geoE1, 36.12778046412819 * M_PI / 180, -114.59241342956928 * M_PI / 180);
	err |= createPt(&geoE2, 36.14329011568265 * M_PI / 180, -114.75321272084759 * M_PI / 180);

	err |= createLocus(&pLoc1, geoS1, geoE1, -2, -2, LineType::SEGMENT, tol, eps);
	err |= createLocus(&sLoc1, geoS1, geoE1, -3, -3, LineType::SEGMENT, tol, eps);
	err |= createLocus(&pLoc2, geoS2, geoE2, -2, -2, LineType::SEGMENT, tol, eps);
	err |= createLocus(&sLoc2, geoS2, geoE2, -3, -3, LineType::SEGMENT, tol, eps);

	err |= invCrs(geoS1, geoE1, &crs1, &bcrs1, eps);
	err |= invCrs(geoS2, geoE2, &crs2, &bcrs2, eps);

	printf("Course 1:  %f\n", crs1);
	printf("Course 2:  %f\n", crs2);

	double cca = modlon(crs2 - crs1);

	printf("Calculated:  %f\n", cca);

	err |= arcTanToTwoLoci(sLoc1, sLoc2, radius, &cp, &sp, &ep, &dir, tol, eps);


}
} //namespace
