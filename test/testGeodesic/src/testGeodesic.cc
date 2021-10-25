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
#include "Geolib.h"
#include "testGeolib.h"

namespace geolib_idealab {


#define nInputs 4
#define nInputs1 21
#define nInputs2 100


/*
 * NAME: testMinSubtendedAngle_Set1
 *
 * DESCRIPTION:
 * 		This function is used to test the minSubtendedAngle function.
 *
 * 		This function runs the test data created by Balakrishna Babu and uses
 * 		his original test plan to execute the actual tests.
 *
 *     	The approach for testing minSubtendedAngle is as follows:
 *
 *     	1. Select a start course and an end course.
 *     	2. Calculate the difference between start and end courses and normalize
 *     	   it to be between -180 to 180, and take the absolute value.
 *     	3. Using the selected data in step 1 as input, compute the Angle between
 *     	   courses.
 *     	4. The result from step 3 should match those calculated in steps 2.
 *
 *
 * INPUT(Type):
 *		None
 *
 * OUTPUT(Return Type):
 * 		testMinSubtendedAngle_Set1(TestSet) - A test set with the folling metrics:
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
TestSet testMinSubtendedAngle_Set1()
{

	double DEG2RAD = M_PI / 180.0;
	//double AZDEGTOL = 5.0e-7; //5e-7 deg or 8.72665e-9 rad or 0.162 cm at 100 nm

    int testNum=0;
    char testname[80];
    double az1, az2;

    double crsCS, crsCE, angleExp, angle;
    int passedCount=0, failedCount=0;
    int errorCount = 0;
    int unverifiedCount = 0;
    int setupFailureCount = 0;
    int testCaseCount = 0;
    ErrorSet err=0;
    TestSet set;

    long newSeed=20080116;

    printf("\n\nStart testMinSubtendedAngle_Set1\n");

    set = newTestSet("testMinSubtendedAngle_Set1");

    srand(newSeed);  //Initialize the random number generator

    while (testNum<1000)
    {
    	err = 0;
        //Select start and end azimuths
        az1 = randAzimuth();
        az2 = randAzimuth();
        testNum++;

        sprintf(testname,"TEST%-d",testNum);

        crsCS = az1 * DEG2RAD;
        crsCE = az2 * DEG2RAD;

        //Calculate angle between courses (always less than 360)
        if (crsCE<crsCS)
        {
            angleExp = (crsCS-crsCE);
        }
        else
        {
            angleExp = (crsCE-crsCS);
        }

        //Normalize to within 180
        if (angleExp>M_PI)
        {
            angleExp = M_2PI-angleExp;
        }

        err |= minSubtendedAngle(crsCS,crsCE,&angle);

        if ( getMaskedError(err, getMaskAll()) )
        {
        	printf("Error occurred in minSubtendedAngle err=0x%lx\n",err);
        	errorCount++;
        	failedCount++;
        	testCaseCount++;
        	continue;
        }

        if (angle == angleExp)
        /*if ( fabs(angle-angleExp) < AZDEGTOL*DEG2RAD )*/
        {
        	if (PRINT_PASSED_CASES) printf("%s azStart=%14.8f  azEnd=%14.8f\n",
                   testname, az1, az2);
        	if (PRINT_PASSED_CASES) printf("%s passed angleExp=%14.8f  angle=%14.8f\n",
                   testname,angleExp/DEG2RAD,angle/DEG2RAD);
            passedCount++;
        }
        else
        {
            printf("%s azStart=%14.8f  azEnd=%14.8f\n",
                   testname, az1, az2);
            printf("%s failed angleExp=%14.8f  angle=%14.8f\n",
                   testname,angleExp/DEG2RAD,angle/DEG2RAD);
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

    printf("\nFinish testMinSubtendedAngle_Set1\n\n\n");

    return set;
}

/*
 * NAME: testMinSubtendedAngle_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the minSubtendedAngle function.
 *
 * 		This function runs all the test cases from set1 and set2 and returns
 * 		the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testMinSubtendedAngle_AllSets(TestSuite) - A test suite with the folling metrics:
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
TestSuite testMinSubtendedAngle_AllSets()
{
    TestSet set1;
    TestSuite suite;

    printf("\nStart testMinSubtendedAngle_AllSets\n");

    suite = newTestSuite("testMinSubtendedAngle_AllSets");

    set1 = testMinSubtendedAngle_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("\nFinish testMinSubtendedAngle_AllSets\n\n\n");

    return suite;

}

/*
 * NAME: testCrsIntx_Set1
 *
 * DESCRIPTION:
 * 		This function is used to test the crsIntx function.
 *
 * 		This function runs the test data created by Balakrishna Babu and uses
 * 		his original test plan to execute the actual tests.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testCrsIntx_Set1(TestSet) - A test set with the folling metrics:
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
TestSet testCrsIntx_Set1()
{
	double DEG2RAD = M_PI / 180.0;

    int testNum=0;
    char testname[80];
    double latxExp, lonxExp;
    double az1, d1, az2, d2;
    LLPoint intxExp;

    double latx, lonx;
    LLPoint p1, p2, intx;
    double crs13, crs31, dist13;
    double crs23, crs32, dist23;
    int passedCount=0, failedCount=0;
    int errorCount = 0;
    int unverifiedCount = 0;
    int setupFailureCount = 0;
    int testCaseCount = 0;
    ErrorSet err=0;

    TestSet set;

    long newSeed=12262007;

    printf("\nStart testCrsIntx_Set1\n");

    set = newTestSet("testCrsIntx_Set1");

    srand(newSeed);  //Initialize the random number generator

    while (testNum < 10000)
    {
        err = 0;

        latxExp = randLat();
        lonxExp = randLon();
        az1 = randAzimuth();
        d1 = 0.8*randDist();
        az2 = randAzimuth();
        d2 = 0.8*randDist();
        testNum++;

        sprintf(testname,"TEST%-d",testNum);

        intxExp.latitude = latxExp * DEG2RAD;
        intxExp.longitude = lonxExp * DEG2RAD;
        crs31 = az1 * DEG2RAD;
        crs32 = az2 * DEG2RAD;

        err |= direct(intxExp, crs31, d1, &p1, EPS);
        err |= direct(intxExp, crs32, d2, &p2, EPS);

        err |= inverse(p1, intxExp, &crs13, &crs31, &dist13, EPS);
        err |= inverse(p2, intxExp, &crs23, &crs32, &dist23, EPS);

        if ( getMaskedError(err, getMaskAll()) )
        {
            printf("Error occurred in pre-crsIntx err=0x%lx\n",err);
            testCaseCount++;
            errorCount++;
            setupFailureCount++;
            continue;
        }

        err |= crsIntx(p1, crs13, &crs31, &dist13,
                                 p2, crs23, &crs32, &dist23,
                                 &intx, TOL, EPS);

        if (err == COLLINEAR_COURSE_ERR)
        {
          if ((fabs(az1 - az2) < TESTTOL) || (fabs(fabs(az1 - az2) - 180.0) < TESTTOL))
          {
            passedCount++;
            testCaseCount++;
            continue; 
          }
        }
 
        if ( getMaskedError(err, getMaskAll()) )
        {
            printf("Error occurred in crsIntx err=0x%lx\n",err);
            errorCount++;
            testCaseCount++;
            failedCount++;
            continue;
        }

        latx = intx.latitude / DEG2RAD;
        lonx = intx.longitude / DEG2RAD;

            if (ptsAreSame(intx, intxExp, TESTTOL))
            {
            	if (PRINT_PASSED_CASES) printf("%s Intx: %14.8f %14.8f  az1=%14.8f d1=%14.8f  az2=%14.8f d2=%14.8f\n",
                       testname, latxExp, lonxExp, az1, d1, az2, d2);
            	if (PRINT_PASSED_CASES) printf("%s p1: %14.8f =%14.8f  p2: %14.8f %14.8f\n",
                       testname,p1.latitude/DEG2RAD,p1.longitude/DEG2RAD,p2.latitude/DEG2RAD,p2.longitude/DEG2RAD);

            	if (PRINT_PASSED_CASES) printf("%s crs13=%14.8f crs31=%14.8f d13=%14.8f  crs23=%14.8f crs32=%14.8f d23=%14.8f\n",
                       testname,crs13/DEG2RAD,crs31/DEG2RAD,dist13,crs23/DEG2RAD,crs32/DEG2RAD,dist23);
            	if (PRINT_PASSED_CASES) printf("%s passed %14.8f %14.8e     %14.8f %14.8e\n",
                       testname,latx,fabs(latxExp-latx),lonx,fabs(lonxExp-lonx));
                passedCount++;
            }
            else
            {
                printf("%s Intx: %14.8f %14.8f  az1=%14.8f d1=%14.8f  az2=%14.8f d2=%14.8f\n",
                       testname, latxExp, lonxExp, az1, d1, az2, d2);
                printf("%s p1: %14.8f =%14.8f  p2: %14.8f %14.8f\n",
                       testname,p1.latitude/DEG2RAD,p1.longitude/DEG2RAD,p2.latitude/DEG2RAD,p2.longitude/DEG2RAD);

                printf("%s crs13=%14.8f crs31=%14.8f d13=%14.8f  crs23=%14.8f crs32=%14.8f d23=%14.8f\n",
                       testname,crs13/DEG2RAD,crs31/DEG2RAD,dist13,crs23/DEG2RAD,crs32/DEG2RAD,dist23);
                printf("%s crs13=%14.8f crs31=%14.8f d13=%14.8f  crs23=%14.8f crs32=%14.8f d23=%14.8f\n",
                       testname,crs13/DEG2RAD,crs31/DEG2RAD,dist13,crs23/DEG2RAD,crs32/DEG2RAD,dist23);
                printf("%s failed %14.8f %14.8f     %14.8f %14.8f\n",
                       testname,latxExp,latx,lonxExp,lonx);
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

    printf("\nFinish testCrsIntx_Set1\n\n\n");

    return set;

}

/*
 * NAME: testCrsIntx_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the crsIntx function.
 *
 * 		This function runs all the test cases from all the test sets for testing
 *      the above function and returns the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testCrsIntx_AllSets(TestSuite) - A test suite with the folling metrics:
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
TestSuite testCrsIntx_AllSets()
{
    TestSuite suite;
    TestSet set1;

    printf("\nStart testCrsIntx_AllSets\n");

    suite = newTestSuite("testCrsIntx_AllSets");

    set1 = testCrsIntx_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("\nFinish testCrsIntx_AllSets\n\n\n");

    return suite;
}

/*
 * NAME: testGeoIntx_Set1
 *
 * DESCRIPTION:
 * 		This function is used to test the geoIntx function.
 *
 * 		This function runs the test data created by Balakrishna Babu and uses
 * 		his original test plan to execute the actual tests.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testGeoIntx_Set1(TestSet) - A test set with the folling metrics:
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
TestSet testGeoIntx_Set1()
{
	double DEG2RAD = M_PI / 180.0;
	long seed = 20080108;
    int testNum=0;
    char testname[80];
    double latxExp, lonxExp;
    double az1, d1, az2, d2;
    LLPoint intxExp;


    double latx, lonx;
    LLPoint p1, p2, intx;
    LLPoint q1, q2;
    double crs13, crs31, dist13;
    double crs23, crs32, dist23;
    double d13, d23;
    LineType lineType1=INFINITE; //0=Segment,1=semiInfinite,2=Infinite
    LineType lineType2=INFINITE; //0=Segment,1=semiInfinite,2=Infinite
    int passedCount=0, failedCount=0;
    int errorCount = 0;
    int unverifiedCount = 0;
    int setupFailureCount = 0;
    int testCaseCount = 0;
    ErrorSet err = 0;

    TestSet set;

    printf("\nStart testGeoIntx_Set1\n");

    set = newTestSet("testGeoIntx_Set1");

    srand(seed);  //Initialize the random number generator

    while (testNum < 10000)
    {
        err = 0;
        latxExp = randLat();
        lonxExp = randLon();
        az1 = randAzimuth();
        d1 = 0.99 * randDist();
        az2 = randAzimuth();
        d2 = 0.99 * randDist();
        testNum++;

        sprintf(testname,"TEST%-d",testNum);

        intxExp.latitude = latxExp * DEG2RAD;
        intxExp.longitude = lonxExp * DEG2RAD;
        crs31 = az1 * DEG2RAD;
        crs32 = az2 * DEG2RAD;

        err |= direct(intxExp, crs31, d1, &p1, EPS);
        err |= direct(intxExp, crs32, d2, &p2, EPS);

        err |= inverse(p1, intxExp, &crs13, &crs31, &dist13, EPS);
        err |= inverse(p2, intxExp, &crs23, &crs32, &dist23, EPS);

        d13 = randDist();
        d23 = randDist();
        err |= direct(p1, crs13, d13, &q1, EPS);
        err |= direct(p2, crs23, d23, &q2, EPS);

        if ( getMaskedError(err, getMaskAll()) )
        {
            printf("Error occurred prior to calling geoIntx err=0x%lx\n",err);
            errorCount++;
            setupFailureCount++;
            testCaseCount++;
            continue;
        }

        err |= geoIntx(p1, q1, lineType1, &crs31, &dist13,
                                  p2, q2, lineType2, &crs32, &dist23,
                                  &intx, TOL, EPS);

        if (err == COLLINEAR_COURSE_ERR)
        {
          if (ptIsOnGeo(p1, q1, p2, INFINITE, &err, TESTTOL, EPS) && ptIsOnGeo (p2, q2, p1, INFINITE, &err, TESTTOL,EPS))
          {
            passedCount++;
            testCaseCount++;
            continue;
          }
        }

        if ( getMaskedError(err, getMaskAll()) )
        {
            printf("Error occurred in geoIntx err=0x%lx\n",err);
            failedCount++;
            errorCount++;
            testCaseCount++;
            continue;
        }

        latx = intx.latitude / DEG2RAD;
        lonx = intx.longitude / DEG2RAD;

        if (ptsAreSame(intx, intxExp, TESTTOL))
        {
        	if (PRINT_PASSED_CASES) printf("%s Intx: %14.8f %14.8f  az1=%14.8f d1=%14.8f  az2=%14.8f d2=%14.8f\n",
                   testname, latxExp, lonxExp, az1, d1, az2, d2);
        	if (PRINT_PASSED_CASES) printf("%s p1: %14.8f %14.8f  p2: %14.8f %14.8f\n",
                   testname,p1.latitude/DEG2RAD,p1.longitude/DEG2RAD,p2.latitude/DEG2RAD,p2.longitude/DEG2RAD);
        	if (PRINT_PASSED_CASES) printf("%s To intxExp:  crs13=%14.8f crs31=%14.8f d13=%14.8f  crs23=%14.8f crs32=%14.8f d23=%14.8f\n",
                   testname,crs13/DEG2RAD,crs31/DEG2RAD,dist13,crs23/DEG2RAD,crs32/DEG2RAD,dist23);
        	if (PRINT_PASSED_CASES) printf("%s q1: %14.8f %14.8f  q2: %14.8f %14.8f\n",
                   testname,q1.latitude/DEG2RAD,q1.longitude/DEG2RAD,q2.latitude/DEG2RAD,q2.longitude/DEG2RAD);
        	if (PRINT_PASSED_CASES) printf("%s d13=%14.8f  d23=%14.8f\n",
                   testname,d13,d23);

        	if (PRINT_PASSED_CASES) printf("%s passed %14.8f %14.8e     %14.8f %14.8e\n",
                   testname,latx,fabs(latxExp-latx),lonx,fabs(lonxExp-lonx));
            passedCount++;
        }
        else
        {
            printf("%s Intx: %14.8f %14.8f  az1=%14.8f d1=%14.8f  az2=%14.8f d2=%14.8f\n",
                   testname, latxExp, lonxExp, az1, d1, az2, d2);
            printf("%s p1: %14.8f %14.8f  p2: %14.8f %14.8f\n",
                   testname,p1.latitude/DEG2RAD,p1.longitude/DEG2RAD,p2.latitude/DEG2RAD,p2.longitude/DEG2RAD);
            printf("%s To intxExp:  crs13=%14.8f crs31=%14.8f d13=%14.8f  crs23=%14.8f crs32=%14.8f d23=%14.8f\n",
                   testname,crs13/DEG2RAD,crs31/DEG2RAD,dist13,crs23/DEG2RAD,crs32/DEG2RAD,dist23);
            printf("%s q1: %14.8f %14.8f  q2: %14.8f %14.8f\n",
                   testname,q1.latitude/DEG2RAD,q1.longitude/DEG2RAD,q2.latitude/DEG2RAD,q2.longitude/DEG2RAD);
            printf("%s d13=%14.8f  d23=%14.8f\n",
                   testname,d13,d23);
            printf("%s crs13=%14.8f crs31=%14.8f d13=%14.8f  crs23=%14.8f crs32=%14.8f d23=%14.8f\n",
                   testname,crs13/DEG2RAD,crs31/DEG2RAD,dist13,crs23/DEG2RAD,crs32/DEG2RAD,dist23);
            printf("%s failed %14.8f %14.8f     %14.8f %14.8f\n",
                   testname,latxExp,latx,lonxExp,lonx);
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

    printf("\nFinish testGeoIntx_Set1\n\n\n");

    return set;

}

/*
 * NAME: testGeoIntx_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the geoIntx function.
 *
 * 		This function runs all the test cases from all the test sets for testing
 *      the above function and returns the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testGeoIntx_AllSets(TestSuite) - A test suite with the folling metrics:
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
TestSuite testGeoIntx_AllSets()
{
	TestSuite suite;
	TestSet set1;

    printf("\nStart testGeoIntx_AllSets\n");

    suite = newTestSuite("testGeoIntx_AllSets");

    set1 = testGeoIntx_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("\nFinish testGeoIntx_AllSets\n\n\n");

    return suite;
}

/*
 * NAME: testProjectToGeo_Set1
 *
 * DESCRIPTION:
 * 		This function is used to test the projectToGeo function.
 *
 * 		This function runs the test data created by Balakrishna Babu and uses
 * 		his original test plan to execute the actual tests.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testProjectToGeo_Set1(TestSet) - A test set with the folling metrics:
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
TestSet testProjectToGeo_Set1()
{
	double DEG2RAD = M_PI / 180.0;

    int testNum=0;
    char testname[80];
    double latxExp, lonxExp;
    double az1, d1, az2, d2;
    LLPoint intxExp;

    double latx, lonx;
    LLPoint p1, p2, intx;
    double crs13, crs31, dist13;
    double crs23, crs32, dist23;
    double distErr=0.0;
    int passedCount=0, failedCount=0;
    int errorCount = 0;
    int unverifiedCount = 0;
    int setupFailureCount = 0;
    int testCaseCount = 0;
    ErrorSet err = 0;
    long newSeed = 20080108;
    TestSet set;

	double dmin = 0.0;
	double dmax = 5300; // There is more than one perpendicular intercept. The algorithm finds the closest.
						// At greater distances the algorithm begins to find other perpendicular intercepts

    printf("\nStart testProjectToGeo_Set1\n");

    set = newTestSet("testProjectToGeo_Set1");

    srand(newSeed);

    while (testNum < 1e4)
    {
        err = 0;

        //Create test data
        //  Start from expected intersection point
        //  azimuth and distance to find point1
        //  azimuth+90 and another distance to find point2
        latxExp = randLat();
        lonxExp = randLon();
        az1 = randAzimuth();
        d1 = dmin + (dmax - dmin) * (double) rand() / (double) RAND_MAX;
        az2 = az1+90;
        if (az2>360)
        {
            az2 -= 360;
        }
        d2 = dmin + (dmax - dmin) * (double) rand() / (double) RAND_MAX;
        testNum++;


        sprintf(testname,"TEST%-d",testNum);

        intxExp.latitude = latxExp * DEG2RAD;
        intxExp.longitude = lonxExp * DEG2RAD;
        crs31 = az1 * DEG2RAD;
        crs32 = az2 * DEG2RAD;

        err |= direct(intxExp, crs31, d1, &p1, EPS);
        err |= direct(intxExp, crs32, d2, &p2, EPS);

        if (d1 == 0.0)
		{
			dist13 = 0.0;
			crs13 = crs31 + M_PI;
		} else {
			err |= inverse(p1, intxExp, &crs13, &crs31, &dist13, EPS);
		}
        err |= inverse(p2, intxExp, &crs23, &crs32, &dist23, EPS);

        if ( getMaskedError(err, getMaskAll()) )
        {
            printf("Error occurred prior to calling projectToGeo err=0x%lx\n",err);
            setupFailureCount++;
            errorCount++;
            testCaseCount++;
            continue;
        }

        err |= projectToGeo(p1, crs13, p2, &intx, &crs23, &dist23, TOL, EPS);

        if ( getMaskedError(err, getMaskAll()) )
        {
            printf("Error occurred in projectToGeo err=0x%lx\n",err);
            failedCount++;
            errorCount++;
            testCaseCount++;
            continue;
        }

        latx = intx.latitude / DEG2RAD;
        lonx = intx.longitude / DEG2RAD;

        err |= invDist(intx, intxExp, &distErr, EPS);

        if ( getMaskedError(err, getMaskAll()) )
        {
            printf("Error occurred in invDist while testing projectToGeo err=0x%lx\n",err);
            failedCount++;
            errorCount++;
            testCaseCount++;
            continue;
        }

        if (ptsAreSame(intx, intxExp, TESTTOL))
        {
        	if (PRINT_PASSED_CASES) printf("%s,%14.8f,%14.8f,%14.8f,%14.8f,%14.8f,%14.8f,",
                   testname, latxExp, lonxExp, az1, d1, az2, d2);
        	if (PRINT_PASSED_CASES) printf("%14.8f,%14.8f,%14.8f,%14.8f,",
                   p1.latitude/DEG2RAD,p1.longitude/DEG2RAD,p2.latitude/DEG2RAD,p2.longitude/DEG2RAD);
        	if (PRINT_PASSED_CASES) printf("%14.8f,%14.8f,%14.8f,%14.8f,%14.8f,%14.8f,",
                   crs13/DEG2RAD,crs31/DEG2RAD,dist13,crs23/DEG2RAD,crs32/DEG2RAD,dist23);
        	if (PRINT_PASSED_CASES) printf("%s passed\n",testname);
        	if (PRINT_PASSED_CASES) printf("%14.8f,%14.8e,%14.8f,%14.8e,%16.8e\n",
                               latx,fabs(latxExp-latx),lonx,fabs(lonxExp-lonx),distErr);
            passedCount++;
        }
        else
        {
            printf("%s, %3.15f,  %3.15f,  %3.15f,  %3.15f,  %3.15e,  %3.15f,  %3.15f,  %3.15f,  %3.15f\n",
            		testname, latxExp,latx,lonxExp,lonx,distErr, az1, d1, az2, d2);

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

    printf("\nFinish testProjectToGeo_Set1\n\n\n");

    return set;
}

/*
 * NAME: testProjectToGeo_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the projectToGeo function.
 *
 * 		This function runs all the test cases from all the test sets for testing
 *      the above function and returns the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testProjectToGeo_AllSets(TestSuite) - A test suite with the folling metrics:
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
TestSuite testProjectToGeo_AllSets()
{
	TestSuite suite;
	TestSet set1;

    printf("\nStart testProjectToGeo_AllSets\n");

    suite = newTestSuite("testProjectToGeo_AllSets");

    set1 = testProjectToGeo_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("\nFinish testProjectToGeo_AllSets\n\n\n");

    return suite;
}

/*
 * NAME: testGeoArcIntx_Set1
 *
 * DESCRIPTION:
 * 		This function is used to test the geoArcIntx function.
 *
 * 		This function runs the test data created by Balakrishna Babu and uses
 * 		his original test plan to execute the actual tests.
 *
 * 		The approach for testing geoArcIntx is as follows:
 * 		1. Start with one expected intersection point.
 * 		2. Select an azimuth and distance to compute start point of Line using direct.
 * 		3. Compute the course from start point of Line to the expected intersection point.
 * 		4. Select an azimuth and radius to compute center point of Arc using direct.
 * 		5. Line and Arc form the input.
 * 		6. One of the resulting two intersection points should match the expected intersection point.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testGeoArcIntx_Set1(TestSet) - A test set with the folling metrics:
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
TestSet testGeoArcIntx_Set1()
{

	double DEG2RAD = M_PI / 180.0;

    int testNum=0;
    char testname[80];
    double latxExp, lonxExp;
    double az1, r1, az2, r2;
    LLPoint intxExp;

    LLPoint p1, c2;
    LLPointPair intx;
    int nX;
    double latx1, lonx1, latx2, lonx2;
    double crs12, crs21, dist12;
    double crs13, crs31, dist13;
    double crs23, crs32, dist23;
    double distExpx1,distExpx2;
    int passedCount=0, failedCount=0;
    int errorCount = 0;
    int unverifiedCount = 0;
    int setupFailureCount = 0;
    int testCaseCount = 0;
    ErrorSet err = 0;

    TestSet set;

    long newSeed=20080102;

    printf("\nStart testGeoArcIntx_Set1\n");

    set = newTestSet("testGeoArcIntx_Set1");

    srand(newSeed);  //Initialize the random number generator

    while (testNum < 10000)
    {
    	err = 0;
        latxExp = randLat();
        lonxExp = randLon();
        az1 = randAzimuth();
        r1 = 0.2 * randDist();  //about 1080 nm
        az2 = randAzimuth();
        r2 = 0.2 * randDist();
        testNum++;

        sprintf(testname,"TEST%-d",testNum);

        intxExp.latitude = latxExp * DEG2RAD;
        intxExp.longitude = lonxExp * DEG2RAD;
        crs31 = az1 * DEG2RAD;
        crs32 = az2 * DEG2RAD;

        err |= direct(intxExp, crs31, r1, &p1, EPS);
        err |= direct(intxExp, crs32, r2, &c2, EPS);

        err |= inverse(p1, intxExp, &crs13, &crs31, &dist13, EPS);
        err |= inverse(c2, intxExp, &crs23, &crs32, &dist23, EPS);
        err |= inverse(p1, c2, &crs12, &crs21, &dist12, EPS);

        if ( getMaskedError(err, getMaskAll()) )
        {
            printf("Error occurred in pre-geoArcIntx err=0x%lx\n",err);
            setupFailureCount++;
            errorCount++;
            testCaseCount++;
            continue;
        }

        intx[0].latitude = 0.0;
        intx[0].longitude = 0.0;
        intx[1].latitude = 0.0;
        intx[1].longitude = 0.0;
        err |= geoArcIntx(p1, crs13, c2, r2, intx, &nX, TOL, EPS);

        if ( getMaskedError(err, getMaskAll()) )
        {
            printf("Error occurred in geoArcIntx err=0x%lx\n",err);
            failedCount++;
            errorCount++;
            testCaseCount++;
            continue;
        }

        latx1 = intx[0].latitude / DEG2RAD;
        lonx1 = intx[0].longitude / DEG2RAD;
        latx2 = intx[1].latitude / DEG2RAD;
        lonx2 = intx[1].longitude / DEG2RAD;

        if ( (nX == 2) &&
             ( ptsAreSame(intx[0], intxExp, TESTTOL) ||
               ptsAreSame(intx[1], intxExp, TESTTOL) ) )
        {
        	if (PRINT_PASSED_CASES) printf("%s Intx: %14.8f %14.8f  az1=%14.8f r1=%14.8f  az2=%14.8f r2=%14.8f\n",
                   testname, latxExp, lonxExp, az1, r1, az2, r2);
        	if (PRINT_PASSED_CASES) printf("%s p1:   %14.8f %14.8f   r1: %14.8f\n",
                   testname,p1.latitude/DEG2RAD,p1.longitude/DEG2RAD,r1);

        	if (PRINT_PASSED_CASES) printf("%s c2:   %14.8f %14.8f   r2: %14.8f\n",
                   testname,c2.latitude/DEG2RAD,c2.longitude/DEG2RAD,r2);

        	if (PRINT_PASSED_CASES) printf("%s intx: %14.8f %14.8f   dist12: %14.8f\n",
                   testname,intxExp.latitude/DEG2RAD,intxExp.longitude/DEG2RAD,dist12);

        	if (PRINT_PASSED_CASES) printf("%s crs13=%14.8f crs31=%14.8f d13=%14.8f  crs23=%14.8f crs32=%14.8f d23=%14.8f\n",
                   testname,crs13/DEG2RAD,crs31/DEG2RAD,dist13,crs23/DEG2RAD,crs32/DEG2RAD,dist23);

        	if (PRINT_PASSED_CASES) printf("%s passed nX: %d  latx1: %14.8f %14.8e    lonx1: %14.8f %14.8e\n",
                   testname,nX,latx1,fabs(latxExp-latx1),lonx1,fabs(lonxExp-lonx1));
        	if (PRINT_PASSED_CASES) printf("%s passed         latx2: %14.8f %14.8e    lonx2: %14.8f %14.8e\n",
                   testname,latx2,fabs(latxExp-latx2),lonx2,fabs(lonxExp-lonx2));
            passedCount++;
        }
        else if ( (nX == 1) &&
        		  ptsAreSame(intx[0], intxExp, TESTTOL) )
        {
        	if (PRINT_PASSED_CASES) printf("%s Intx: %14.8f %14.8f  az1=%14.8f r1=%14.8f  az2=%14.8f r2=%14.8f\n",
                   testname, latxExp, lonxExp, az1, r1, az2, r2);
        	if (PRINT_PASSED_CASES) printf("%s p1:   %14.8f %14.8f   r1: %14.8f\n",
                   testname,p1.latitude/DEG2RAD,p1.longitude/DEG2RAD,r1);

        	if (PRINT_PASSED_CASES) printf("%s c2:   %14.8f %14.8f   r2: %14.8f\n",
                   testname,c2.latitude/DEG2RAD,c2.longitude/DEG2RAD,r2);

        	if (PRINT_PASSED_CASES) printf("%s intx: %14.8f %14.8f   dist12: %14.8f\n",
                   testname,intxExp.latitude/DEG2RAD,intxExp.longitude/DEG2RAD,dist12);

        	if (PRINT_PASSED_CASES) printf("%s crs13=%14.8f crs31=%14.8f d13=%14.8f  crs23=%14.8f crs32=%14.8f d23=%14.8f\n",
                   testname,crs13/DEG2RAD,crs31/DEG2RAD,dist13,crs23/DEG2RAD,crs32/DEG2RAD,dist23);

        	if (PRINT_PASSED_CASES) printf("%s passed nX: %d  latx1: %14.8f %14.8e    lonx1: %14.8f %14.8e\n",
                   testname,nX,latx1,fabs(latxExp-latx1),lonx1,fabs(lonxExp-lonx1));
            passedCount++;
        }
        else
        {
            printf("%s Intx: %.16lf %.16lf  az1=%.16lf r1=%.16lf  az2=%.16lf r2=%.16lf\n",
                   testname, latxExp, lonxExp, az1, r1, az2, r2);
            printf("%s p1:   %.16lf %.16lf   r1: %.16lf\n",
                   testname,p1.latitude/DEG2RAD,p1.longitude/DEG2RAD,r1);

            printf("%s c2:   %.16lf %.16lf   r2: %.16lf\n",
                   testname,c2.latitude/DEG2RAD,c2.longitude/DEG2RAD,r2);

            printf("%s intx: %.16lf %.16lf   dist12: %.16lf\n",
                   testname,intxExp.latitude/DEG2RAD,intxExp.longitude/DEG2RAD,dist12);

            printf("%s crs13=%.16lf crs31=%.16lf d13=%.16lf  crs23=%.16lf crs32=%.16lf d23=%.16lf\n",
                   testname,crs13/DEG2RAD,crs31/DEG2RAD,dist13,crs23/DEG2RAD,crs32/DEG2RAD,dist23);

            printf("%s failed nX: %d  Exp: %16.10f %16.10f   intx1: %16.10f %16.10f  intx2: %16.10f %16.10f\n",
                   testname,nX,latxExp,lonxExp,latx1,lonx1,latx2,lonx2);
            err |= invDist(intxExp,intx[0],&distExpx1,EPS);
            err |= invDist(intxExp,intx[1],&distExpx2,EPS);

            if ( getMaskedError(err, getMaskAll()) )
            {
                printf("Error occurred in post-geoArcIntx err=0x%lx\n",err);
                failedCount++;
                errorCount++;
                testCaseCount++;
                continue;
            }

            printf("%s failed         distExpx1=%18.8e  distExpx2=%18.8e\n",
                   testname,distExpx1,distExpx2);
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

    printf("\nFinish testGeoArcIntx_Set1\n\n\n");

    return set;

}

/*
 * NAME: testGeoArcIntx_Set2
 *
 * DESCRIPTION:
 * 		This function is used to test the geoArcIntx function.  Specifically,
 * 		this function is used to test for zero intersections between a line and arc.
 *
 * 		This function runs the test data created by Balakrishna Babu and uses
 * 		his original test plan to execute the actual tests.
 *
 *  	The approach for testing geoArcIntx is as follows:
 *     	1. Start with one expected intersection point.
 *     	2. Select an azimuth and distance to compute start point of Line using direct.
 *     	3. Compute the course from start point of Line to the expected intersection point.
 *     	4. Select an azimuth and radius to compute center point of Arc using direct.
 *     	5. Compute the perpendicular intercept from center of Arc on the Line using projectToGeo.
 *     	6. Set the radius of the Arc equal to the perpendicular intercept distance minus distance tolerance.
 *     	7. Line and Arc form the input.
 *     	8. Result should be zero intersection points.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testGeoArcIntx_Set2(TestSet) - A test set with the folling metrics:
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

TestSet testGeoArcIntx_Set2()
{

	double DEG2RAD = M_PI / 180.0;
	double NMTOL = 0.03 / 100.0 / 1852.0;  //0.03 cm or ~1.62e-7 nm

    int testNum=0;
    char testname[80];
    double latxExp, lonxExp;
    double az1, r1, az2, r2;
    LLPoint intxExp;

    LLPoint p1, c2;
    LLPointPair intx;
    int nX;
    double latx1, lonx1, latx2, lonx2;
    double crs12, crs21, dist12;
    double crs13, crs31, dist13;
    double crs23, crs32, dist23;
    int passedCount=0, failedCount=0;
    int errorCount = 0;
    int unverifiedCount = 0;
    int setupFailureCount = 0;
    int testCaseCount = 0;
    ErrorSet err = 0;

    TestSet set;

    long newSeed=20080102;

    printf("\nStart testGeoArcIntx_Set2\n");

    set = newTestSet("testGeoArcIntx_Set2");

    srand(newSeed);  //Initialize the random number generator

    while (testNum < 10000)
    {
    	err = 0;
        latxExp = randLat();
        lonxExp = randLon();
        az1 = randAzimuth();
        r1 = 0.2 * randDist();  //about 1080 nm
        az2 = randAzimuth();
        r2 = 0.2 * randDist();
        testNum++;

        sprintf(testname,"TEST%-d",testNum);

        intxExp.latitude = latxExp * DEG2RAD;
        intxExp.longitude = lonxExp * DEG2RAD;
        crs31 = az1 * DEG2RAD;
        crs32 = az2 * DEG2RAD;

        err |= direct(intxExp, crs31, r1, &p1, EPS);
        err |= direct(intxExp, crs32, r2, &c2, EPS);

        err |= inverse(p1, intxExp, &crs13, &crs31, &dist13, EPS);
        err |= inverse(c2, intxExp, &crs23, &crs32, &dist23, EPS);
        err |= inverse(p1, c2, &crs12, &crs21, &dist12, EPS);

        err |= projectToGeo(p1, crs13, c2, &intxExp, &crs23, &r2, TOL, EPS);
        latxExp = intxExp.latitude / DEG2RAD;
        lonxExp = intxExp.longitude / DEG2RAD;
        r2 -= NMTOL;

        if( getMaskedError(err, getMaskAll()) )
        {
            printf("Error occurred in pre-geoArcIntx err=0x%lx\n",err);
            setupFailureCount++;
            errorCount++;
            testCaseCount++;
            continue;
        }

        intx[0].latitude = 0.0;
        intx[0].longitude = 0.0;
        intx[1].latitude = 0.0;
        intx[1].longitude = 0.0;
        err |= geoArcIntx(p1, crs13, c2, r2, intx, &nX, TOL, EPS);

        if( getMaskedError(err, getMaskAll()) )
        {
            printf("Error occurred in geoArcIntx err=0x%lx\n",err);
            failedCount++;
            errorCount++;
            testCaseCount++;
            continue;
        }

        latx1 = intx[0].latitude / DEG2RAD;
        lonx1 = intx[0].longitude / DEG2RAD;
        latx2 = intx[1].latitude / DEG2RAD;
        lonx2 = intx[1].longitude / DEG2RAD;

        if ( nX == 0 )
        {
        	if (PRINT_PASSED_CASES) printf("%s Intx: %14.8f %14.8f  az1=%14.8f r1=%14.8f  az2=%14.8f r2=%14.8f\n",
                   testname, latxExp, lonxExp, az1, r1, az2, r2);
        	if (PRINT_PASSED_CASES) printf("%s p1:   %14.8f %14.8f   r1: %14.8f\n",
                   testname,p1.latitude/DEG2RAD,p1.longitude/DEG2RAD,r1);

        	if (PRINT_PASSED_CASES) printf("%s c2:   %14.8f %14.8f   r2: %14.8f\n",
                   testname,c2.latitude/DEG2RAD,c2.longitude/DEG2RAD,r2);

        	if (PRINT_PASSED_CASES) printf("%s intx: %14.8f %14.8f   dist12: %14.8f\n",
                   testname,intxExp.latitude/DEG2RAD,intxExp.longitude/DEG2RAD,dist12);

        	if (PRINT_PASSED_CASES) printf("%s crs13=%14.8f crs31=%14.8f d13=%14.8f  crs23=%14.8f crs32=%14.8f d23=%14.8f\n",
                   testname,crs13/DEG2RAD,crs31/DEG2RAD,dist13,crs23/DEG2RAD,crs32/DEG2RAD,dist23);

        	if (PRINT_PASSED_CASES) printf("%s passed nX: %d\n",
                   testname,nX);
            passedCount++;
        }
        else
        {
            printf("%s Intx: %14.8f %14.8f  az1=%14.8f r1=%14.8f  az2=%14.8f r2=%14.8f\n",
                   testname, latxExp, lonxExp, az1, r1, az2, r2);
            printf("%s p1:   %14.8f %14.8f   r1: %14.8f\n",
                   testname,p1.latitude/DEG2RAD,p1.longitude/DEG2RAD,r1);

            printf("%s c2:   %14.8f %14.8f   r2: %14.8f\n",
                   testname,c2.latitude/DEG2RAD,c2.longitude/DEG2RAD,r2);

            printf("%s intx: %14.8f %14.8f   dist12: %14.8f\n",
                   testname,intxExp.latitude/DEG2RAD,intxExp.longitude/DEG2RAD,dist12);

            printf("%s crs13=%14.8f crs31=%14.8f d13=%14.8f  crs23=%14.8f crs32=%14.8f d23=%14.8f\n",
                   testname,crs13/DEG2RAD,crs31/DEG2RAD,dist13,crs23/DEG2RAD,crs32/DEG2RAD,dist23);

            printf("%s failed nX: %d  Exp: %14.8f %14.8f   intx1: %14.8f %14.8f  intx2: %14.8f %14.8f\n",
                   testname,nX,latxExp,lonxExp,latx1,lonx1,latx2,lonx2);
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

    printf("\nFinish testGeoArcIntx_Set2\n\n\n");

    return set;

}

/*
 * NAME: testGeoArcIntx_Set3
 *
 * DESCRIPTION:
 * 		This function is used to test the geoArcIntx function.  Specifically,
 * 		this function is used to test for one intersection between a line and arc.
 *
 * 		This function runs the test data created by Balakrishna Babu and uses
 * 		his original test plan to execute the actual tests.
 *
 *    	The approach for testing geoArcIntx is as follows:
 *     	1. Start with one expected intersection point.
 *     	2. Select an azimuth and distance to compute start point of Line using direct.
 *     	3. Compute the course from start point of Line to the expected intersection point.
 *     	4. Select an azimuth and radius to compute center point of Arc using direct.
 *     	5. Compute the perpendicular intercept from center of Arc on the Line using projectToGeo.
 *     	6. Set the radius of the Arc equal to the perpendicular intercept distance.
 *     	7. Perpendicular intercept point is the touching point between the Line and the Arc.
 *     	8. Line and Arc form the input.
 *     	9. Resulting intersection point should match the touching point.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testGeoArcIntx_Set3(TestSet) - A test set with the folling metrics:
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

TestSet testGeoArcIntx_Set3()
{

	double DEG2RAD = M_PI / 180.0;

    int testNum=0;
    double latxExp, lonxExp;
    double az1, r1, az2, r2;
    LLPoint intxExp;

    LLPoint p1, c2;
    LLPointPair intx;
    int nX;
    double latx1, lonx1, latx2, lonx2;
    double crs12, crs21, dist12;
    double crs13, crs31, dist13;
    double crs23, crs32, dist23;
    int passedCount=0, failedCount=0;
    int errorCount = 0;
    int unverifiedCount = 0;
    int setupFailureCount = 0;
    int testCaseCount = 0;
    ErrorSet err = 0;

    TestSet set;

    char testname[80];

    long newSeed=20080102;

    printf("\nStart testGeoArcIntx_Set3\n");

    set = newTestSet("testGeoArcIntx_Set3");

    srand(newSeed);  //Initialize the random number generator

    while (testNum < 10000)
    {
    	err = 0;
        latxExp = randLat();
        lonxExp = randLon();
        az1 = randAzimuth();
        r1 = 0.2 * randDist();  //about 1080 nm
        az2 = randAzimuth();
        r2 = 0.2 * randDist();
        testNum++;

        sprintf(testname,"TEST%-d",testNum);

        intxExp.latitude = latxExp * DEG2RAD;
        intxExp.longitude = lonxExp * DEG2RAD;
        crs31 = az1 * DEG2RAD;
        crs32 = az2 * DEG2RAD;

        err |= direct(intxExp, crs31, r1, &p1, EPS);
        err |= direct(intxExp, crs32, r2, &c2, EPS);

        err |= inverse(p1, intxExp, &crs13, &crs31, &dist13, EPS);
        err |= inverse(c2, intxExp, &crs23, &crs32, &dist23, EPS);
        err |= inverse(p1, c2, &crs12, &crs21, &dist12, EPS);

        err |= projectToGeo(p1, crs13, c2, &intxExp, &crs23, &r2, TOL, EPS);
        latxExp = intxExp.latitude / DEG2RAD;
        lonxExp = intxExp.longitude / DEG2RAD;

        if( getMaskedError(err, getMaskAll()) )
        {
            printf("Error occurred in pre-geoArcIntx err=0x%lx\n",err);
            setupFailureCount++;
            errorCount++;
            testCaseCount++;
            continue;
        }

        intx[0].latitude = 0.0;
        intx[0].longitude = 0.0;
        intx[1].latitude = 0.0;
        intx[1].longitude = 0.0;
        err |= geoArcIntx(p1, crs13, c2, r2, intx, &nX, TOL, EPS);

        if( getMaskedError(err, getMaskAll()) )
        {
            printf("Error occurred in geoArcIntx err=0x%lx\n",err);
            failedCount++;
            errorCount++;
            testCaseCount++;
            continue;
        }

        latx1 = intx[0].latitude / DEG2RAD;
        lonx1 = intx[0].longitude / DEG2RAD;
        latx2 = intx[1].latitude / DEG2RAD;
        lonx2 = intx[1].longitude / DEG2RAD;

        if (ptsAreSame(intx[0], intxExp, TESTTOL) && (nX == 1))
        {
        	if (PRINT_PASSED_CASES) printf("%s passed nX: %d  latx1: %14.8f %14.8e    lonx1: %14.8f %14.8e\n",
                   testname,nX,latx1,fabs(latxExp-latx1),lonx1,fabs(lonxExp-lonx1));
        	if (PRINT_PASSED_CASES) printf("%d,%4.15f,%4.15f,%4.15f,%4.15f,%4.15f,%4.15f,%4.15f,%4.15f,%d,%d,%d\n", testNum, p1.latitude/DEG2RAD,
        	        			p1.longitude/DEG2RAD, crs13/DEG2RAD, c2.latitude/DEG2RAD, c2.longitude/DEG2RAD, r2, latx1,lonx1, -1, -1, nX);
            passedCount++;
        }
        else
        {
            printf("%s failed nX: %d  Exp: %14.8f %14.8f   intx1: %14.8f %14.8f  intx2: %14.8f %14.8f\n",
                   testname,nX,latxExp,lonxExp,latx1,lonx1,latx2,lonx2);
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

    printf("\nFinish testGeoArcIntx_Set3\n\n\n");

    return set;

}

/*
 * NAME: testGeoArcIntx_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the geoArcIntx function.
 *
 * 		This function runs all the test cases from all the test sets for testing
 *      the above function and returns the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testGeoArcIntx_AllSets(TestSuite) - A test suite with the folling metrics:
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
TestSuite testGeoArcIntx_AllSets()
{
	TestSuite suite;
	TestSet set1;
	TestSet set2;
	TestSet set3;

    printf("\nStart testGeoArcIntx_AllSets\n");

    suite = newTestSuite("testGeoArcIntx_AllSets");

    set1 = testGeoArcIntx_Set1();
    addTestSet(set1,&suite);

    set2 = testGeoArcIntx_Set2();
    addTestSet(set2,&suite);

    set3 = testGeoArcIntx_Set3();
    addTestSet(set3,&suite);

    displayTestSuite(suite);

    printf("\nFinish testGeoArcIntx_AllSets\n\n\n");

    return suite;
}

/*
 * NAME: testArcTanToTwoGeos_Set1
 *
 * DESCRIPTION:
 * 		This function is used to test the arcTanToTwoGeos function.
 *
 * 		This function runs the test data created by Balakrishna Babu and uses
 * 		his original test plan to execute the actual tests.
 *
 * 		This specific function tests strictly for the creation of the tangent
 * 		arc in the correct quadrant:
 *
 *                                2  |  1
 *                               ---------
 *                                3  |  4
 *
 *     	The approach for testing arcTanToTwoGeos is as follows:
 *     	1. Start with an expected center point.
 *     	2. Select an azimuth and radius to compute start point using direct.
 *     	3. Select an azimuth that differs from previous by 1,10,90,150,180,220,270,330,359
 *     	   degrees and same radius to compute end point using direct.
 *     	4. Compute azimuth from the start point to center and add 90 degrees for tangent's
 *     	   azimuth at start using inverse.
 *     	5. Select a distance from start point and compute point 1 along tangent using direct.
 *     	6. Repeat the previous two steps for the end point by subtracting 90 degrees
 *     	   instead of adding and compute point 3.
 *     	7. Compute azimuths from point 1 to start, and point 3 to end using inverse.
 *     	8. Points 1 and 3 and their respective azimuths, and the radius from step 2 form the input.
 *     	9. Verify that the resulting centerPoint, StartPoint, endPoint match the expected center,
 *     	   start and end points from steps 1,2,3.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testArcTanToTwoGeos_Set1(TestSet) - A test set with the folling metrics:
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
TestSet testArcTanToTwoGeos_Set1()
{
	double DEG2RAD = M_PI / 180.0;
    char testname[50];
    LLPoint start, center, end;
    ArcDirection dir;


    LLPoint pp[3], qq[3];
    double crsSTangent, crsETangent, crsS, crsT;
    int ip, iq;
    double radius1;

    int passedCount=0, failedCount=0;
    int errorCount = 0;
    int unverifiedCount = 0;
    int setupFailureCount = 0;
    int testCaseCount = 0;
    ErrorSet err = 0;

    TestSet set;

    printf("\nStart testArcTanToTwoGeos_Set1\n");

    set = newTestSet("testArcTanToTwoGeos_Set1");

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

    for (ip=0; ip<3; ip++) {

    	for (iq=0; iq<3; iq++) {

            for (crsS=90.0; crsS<360.0; crsS=crsS+180.0) {

	            crsSTangent = crsS * DEG2RAD;

	            for (crsT=0.0; crsT<360.0; crsT=crsT+180.0) {

	            	sprintf(testname, "TEST %i", testCaseCount+1);

		            crsETangent = crsT * DEG2RAD;

		            err = 0;
		            center.latitude = 0.0;
		            center.longitude = 0.0;
		            start.latitude = 0.0;
		            start.longitude = 0.0;
		            end.latitude = 0.0;
		            end.longitude = 0.0;
		            //dir = 0;
		            err |= arcTanToTwoGeos(pp[ip],crsSTangent,qq[iq],crsETangent,radius1,
		                                              &center,&start,&end,&dir,TOL,EPS);

		            if ( getMaskedError(err, getMaskAll()) )
		            {
		                printf("%s - Error occurred in arcTanToTwoGeos err=0x%lx\n", testname, err);
		                failedCount++;
		                errorCount++;
		                testCaseCount++;
		                continue;
		            }

		            if (getMaskedError(err, getMaskAll()) == SUCCESS)
		            {
		            	if ((crsS == 90.0) && (crsT == 0.0))
		            	{
			            	if ( (quadrant(center,pp[1])==2) && (dir==COUNTERCLOCKWISE) )
			            	{
			            		if (PRINT_PASSED_CASES) printf("%s ip=%d iq=%d crsS=%6.1f crsT=%6.1f  err=%ld\n",
				            			testname, ip,iq,crsS,crsT,err);
			            		if (PRINT_PASSED_CASES) printf("passed quadExp=2  dirExp=1  quad=%d  dir=%d\n",
			            				quadrant(center,pp[1]),dir);
			            		if (PRINT_PASSED_CASES) printLLPoint("  pt1  ",pp[ip]);
			            		if (PRINT_PASSED_CASES) printLLPoint("  pt3  ",qq[iq]);
			            		if (PRINT_PASSED_CASES) printf("   crs12=%16.10f  crs3=%16.10f\n",crsS,crsT);
			            		if (PRINT_PASSED_CASES) printLLPoint("   start     ",start);
			            		if (PRINT_PASSED_CASES) printLLPoint("   center    ",center);
			            		if (PRINT_PASSED_CASES) printLLPoint("   end       ",end);
			            		if (PRINT_PASSED_CASES) printf("     dir=%d\n",dir);



			            		passedCount++;
			            	}
			            	else
			            	{
				            	printf("%s ip=%d iq=%d crsS=%6.1f crsT=%6.1f  err=%ld\n",
				            			testname, ip,iq,crsS,crsT,err);
			            		printf("failed quadExp=2  dirExp=1  quad=%d  dir=%d\n",
			            				quadrant(center,pp[1]),dir);
					            printLLPoint("  pt1  ",pp[ip]);
					            printLLPoint("  pt3  ",qq[iq]);
					            printf("   crs12=%16.10f  crs3=%16.10f\n",crsS,crsT);
					            printLLPoint("   start     ",start);
					            printLLPoint("   center    ",center);
					            printLLPoint("   end       ",end);
					            printf("     dir=%d\n",dir);



			            		failedCount++;
			            	}
			            }
		            	else if ((crsS == 90.0) && (crsT == 180.0))
		            	{
			            	if ( (quadrant(center,pp[1])==3) && (dir==CLOCKWISE) )
			            	{
			            		if (PRINT_PASSED_CASES) printf("%s ip=%d iq=%d crsS=%6.1f crsT=%6.1f  err=%ld\n",
				            			testname, ip,iq,crsS,crsT,err);
			            		if (PRINT_PASSED_CASES) printf("passed quadExp=3  dirExp=-1  quad=%d  dir=%d\n",
			            				quadrant(center,pp[1]),dir);
			            		if (PRINT_PASSED_CASES) printLLPoint("  pt1  ",pp[ip]);
			            		if (PRINT_PASSED_CASES) printLLPoint("  pt3  ",qq[iq]);
			            		if (PRINT_PASSED_CASES) printf("   crs12=%16.10f  crs3=%16.10f\n",crsS,crsT);
			            		if (PRINT_PASSED_CASES) printLLPoint("   start     ",start);
			            		if (PRINT_PASSED_CASES) printLLPoint("   center    ",center);
			            		if (PRINT_PASSED_CASES) printLLPoint("   end       ",end);
			            		if (PRINT_PASSED_CASES) printf("     dir=%d\n",dir);



			            		passedCount++;
			            	}
			            	else
			            	{
				            	printf("%s ip=%d iq=%d crsS=%6.1f crsT=%6.1f  err=%ld\n",
				            			testname, ip,iq,crsS,crsT,err);
			            		printf("failed quadExp=3  dirExp=-1  quad=%d  dir=%d\n",
			            				quadrant(center,pp[1]),dir);
					            printLLPoint("  pt1  ",pp[ip]);
					            printLLPoint("  pt3  ",qq[iq]);
					            printf("   crs12=%16.10f  crs3=%16.10f\n",crsS,crsT);
					            printLLPoint("   start     ",start);
					            printLLPoint("   center    ",center);
					            printLLPoint("   end       ",end);
					            printf("     dir=%d\n",dir);



			            		failedCount++;
			            	}
		            	}
		            	else if ((crsS == 270.0) && (crsT == 0.0))
		            	{
			            	if ( (quadrant(center,pp[1])==1) && (dir==CLOCKWISE) )
			            	{
			            		if (PRINT_PASSED_CASES) printf("%s ip=%d iq=%d crsS=%6.1f crsT=%6.1f  err=%ld\n",
				            			testname, ip,iq,crsS,crsT,err);
			            		if (PRINT_PASSED_CASES) printf("passed quadExp=1  dirExp=-1  quad=%d  dir=%d\n",
			            				quadrant(center,pp[1]),dir);
			            		if (PRINT_PASSED_CASES) printLLPoint("  pt1  ",pp[ip]);
			            		if (PRINT_PASSED_CASES) printLLPoint("  pt3  ",qq[iq]);
			            		if (PRINT_PASSED_CASES) printf("   crs12=%16.10f  crs3=%16.10f\n",crsS,crsT);
			            		if (PRINT_PASSED_CASES) printLLPoint("   start     ",start);
			            		if (PRINT_PASSED_CASES) printLLPoint("   center    ",center);
			            		if (PRINT_PASSED_CASES) printLLPoint("   end       ",end);
			            		if (PRINT_PASSED_CASES) printf("     dir=%d\n",dir);



			            		passedCount++;
			            	}
			            	else
			            	{
				            	printf("%s ip=%d iq=%d crsS=%6.1f crsT=%6.1f  err=%ld\n",
				            			testname, ip,iq,crsS,crsT,err);
			            		printf("failed quadExp=1  dirExp=-1  quad=%d  dir=%d\n",
			            				quadrant(center,pp[1]),dir);
					            printLLPoint("  pt1  ",pp[ip]);
					            printLLPoint("  pt3  ",qq[iq]);
					            printf("   crs12=%16.10f  crs3=%16.10f\n",crsS,crsT);
					            printLLPoint("   start     ",start);
					            printLLPoint("   center    ",center);
					            printLLPoint("   end       ",end);
					            printf("     dir=%d\n",dir);



			            		failedCount++;
			            	}
		            	}
		            	else if ((crsS == 270.0) && (crsT == 180.0))
		            	{
			            	if ( (quadrant(center,pp[1])==4) && (dir==COUNTERCLOCKWISE) )
			            	{
			            		if (PRINT_PASSED_CASES) printf("%s ip=%d iq=%d crsS=%6.1f crsT=%6.1f  err=%ld\n",
				            			testname, ip,iq,crsS,crsT,err);
			            		if (PRINT_PASSED_CASES) printf("passed quadExp=4  dirExp=1  quad=%d  dir=%d\n",
			            				quadrant(center,pp[1]),dir);
			            		if (PRINT_PASSED_CASES) printLLPoint("  pt1  ",pp[ip]);
			            		if (PRINT_PASSED_CASES) printLLPoint("  pt3  ",qq[iq]);
			            		if (PRINT_PASSED_CASES) printf("   crs12=%16.10f  crs3=%16.10f\n",crsS,crsT);
			            		if (PRINT_PASSED_CASES) printLLPoint("   start     ",start);
			            		if (PRINT_PASSED_CASES) printLLPoint("   center    ",center);
			            		if (PRINT_PASSED_CASES) printLLPoint("   end       ",end);
			            		if (PRINT_PASSED_CASES) printf("     dir=%d\n",dir);



			            		passedCount++;
			            	}
			            	else
			            	{
				            	printf("%s ip=%d iq=%d crsS=%6.1f crsT=%6.1f  err=%ld\n",
				            			testname, ip,iq,crsS,crsT,err);
			            		printf("failed quadExp=4  dirExp=1  quad=%d  dir=%d\n",
			            				quadrant(center,pp[1]),dir);
					            printLLPoint("  pt1  ",pp[ip]);
					            printLLPoint("  pt3  ",qq[iq]);
					            printf("   crs12=%16.10f  crs3=%16.10f\n",crsS,crsT);
					            printLLPoint("   start     ",start);
					            printLLPoint("   center    ",center);
					            printLLPoint("   end       ",end);
					            printf("     dir=%d\n",dir);



			            		failedCount++;
			            	}
		            	}
		            }
		            else
		            {
		            	printf("failed\n");
		            	errorCount++;
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

    printf("\nFinish testArcTanToTwoGeos_Set1\n\n\n");

    return set;
}

/*
 * NAME: testArcTanToTwoGeos_Set2
 *
 * DESCRIPTION:
 * 		This function is used to test the arcTanToTwoGeos function.
 *
 * 		This function runs the test data created by Balakrishna Babu and uses
 * 		his original test plan to execute the actual tests.
 *
 *     	The approach for testing arcTanToTwoGeos is as follows:
 *     	1. Start with an expected center point.
 *     	2. Select an azimuth and radius to compute start point using direct.
 *     	3. Select an azimuth that differs from previous by 1,10,90,150,180,220,270,330,359
 *     	   degrees and same radius to compute end point using direct.
 *     	4. Compute azimuth from the start point to center and add 90 degrees for tangent's
 *     	   azimuth at start using inverse.
 *     	5. Select a distance from start point and compute point 1 along tangent using direct.
 *     	6. Repeat the previous two steps for the end point by subtracting 90 degrees
 *     	   instead of adding and compute point 3.
 *     	7. Compute azimuths from point 1 to start, and point 3 to end using inverse.
 *     	8. Points 1 and 3 and their respective azimuths, and the radius from step 2 form the input.
 *     	9. Verify that the resulting centerPoint, StartPoint, endPoint match the expected center,
 *     	   start and end points from steps 1,2,3.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testArcTanToTwoGeos_Set2(TestSet) - A test set with the folling metrics:
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
TestSet testArcTanToTwoGeos_Set2()
{
	double DEG2RAD = M_PI / 180.0;

    int testNum=0;
    char testname[80];
    double latxExp, lonxExp;
    double az1, r1;

    int ir=0;
    int irMax=10;
    double radius[10] = {0.01,0.1,0.2,0.5,1.0,5.0,10.0,20.0,50.0,100.0};

    int iz=0;
    int izMax=20;

    int ic=0;
    LLPoint pt1Arr[4],pt3Arr[4];
    double crs1SArr[4],crs3EArr[4];

    LLPoint start, center, end;
    ArcDirection dir;


    double crsSTangent, crsETangent;

    LLPoint pt1, pt3, intx;
    LLPoint startExp, centerExp, endExp;
    int dirExp;  //CLOCKWISE = 1, COUNTERCLOCKWISE = -1

    double crsX1, crsX3, dist1X, dist3X;
    double crs1S, crsS1, dist1S, distS1;
    double crs3E, crsE3, dist3E, distE3;
    double crsCS, crsSC, distCS;
    double crsCE, crsEC, distCE;

    double centerDistErr, startDistErr, endDistErr;

    int passedCount=0, failedCount=0;
    int errorCount = 0;
    int unverifiedCount = 0;
    int setupFailureCount = 0;
    int testCaseCount = 0;
    ErrorSet err = 0;
    double azmin, azmax, azr;
    TestSet set;

    long newSeed=20080108;

    printf("\nStart testArcTanToTwoGeos_Set2\n");

    set = newTestSet("testArcTanToTwoGeos_Set2");

    srand(newSeed);  //Initialize the random number generator

    while (testNum<10)
    {
        latxExp = randLat();
        lonxExp = randLon();
        az1 = randAzimuth();
        dist1S = 0.01 * randDist(); //0-54 nm
        dist3E = 0.01 * randDist(); //0-54 nm
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
	            crsCE = crsCS + azr;
		        if ( crsCE >= M_2PI )
		        {
		            crsCE -= M_2PI;
		        }
		        else if ( crsCE < 0.0 )
		        {
		            crsCE += M_2PI;
		        }

		        //Compute the start and end points of Arc
		        err |= direct(centerExp, crsCS, r1, &startExp, EPS);
		        err |= direct(centerExp, crsCE, r1, &endExp, EPS);

		        sprintf(testname,"TEST%-d-%-d-%-d",testNum,ir+1,iz+1);

		        //Compute point 1 and course from pt1 to start of Arc
		        err |= inverse(startExp, centerExp, &crsSC, &crsCS, &distCS, EPS);
		        crsSTangent = crsSC + dirExp*M_PI_2;
		        if ( crsSTangent >= M_2PI )
		        {
		        	crsSTangent -= M_2PI;
		        }
		        else if ( crsSTangent < 0.0 )
		        {
		        	crsSTangent += M_2PI;
		        }
		        err |= direct(startExp, crsSTangent, dist1S, &pt1, EPS);
		        err |= inverse(pt1, startExp, &crs1S, &crsS1, &distS1, EPS);

		        //Compute point 3 and course from pt3 to end of Arc
		        err |= inverse(endExp, centerExp, &crsEC, &crsCE, &distCE, EPS);
		        crsETangent = crsEC - dirExp*M_PI_2;
		        if ( crsETangent >= M_2PI )
		        {
		        	crsETangent -= M_2PI;
		        }
		        else if ( crsETangent < 0.0 )
		        {
		        	crsETangent += M_2PI;
		        }
		        err |= direct(endExp, crsETangent, dist3E, &pt3, EPS);
		        err |= inverse(pt3, endExp, &crs3E, &crsE3, &distE3, EPS);

		        //Compute course from end of Arc to pt3 at pt3
		        crs3E += M_PI;
		        if (crs3E >= M_2PI)
		        {
		            crs3E -= M_2PI;
		        }

		        err |= crsIntx(pt1,crs1S, &crsX1,&dist1X,
		                                 pt3,crs3E, &crsX3,&dist3X,
		                                 &intx,TOL,EPS);

		        if ( getMaskedError(err, getMaskAll()) )
		        {
		            printf("%s - Error occurred in pre-arcTanToTwoGeos err=0x%lx\n",testname, err);
		            errorCount++;
		            setupFailureCount++;
		            testCaseCount++;
		            continue;
		        }

		        //Set up each case
		        // 1. pt1, crs1S, pt3, crs3E
		        // 2. pt1, crs1S, intx, crsX3
		        // 3. intx, crsX1+pi, pt3, crs3E
		        // 4. intx, crsX1+pi, intx, crsX3
		        pt1Arr[0] = pt1;
		        crs1SArr[0] = crs1S;
		        pt3Arr[0] = pt3;
		        crs3EArr[0] = crs3E;

		        pt1Arr[1] = pt1;
		        crs1SArr[1] = crs1S;
		        pt3Arr[1] = intx;
		        crs3EArr[1] = crsX3;

		        pt1Arr[2] = intx;
		        crs1SArr[2] = crsX1+M_PI;
		        pt3Arr[2] = pt3;
		        crs3EArr[2] = crs3E;

		        pt1Arr[3] = intx;
		        crs1SArr[3] = crsX1+M_PI;
		        pt3Arr[3] = intx;
		        crs3EArr[3] = crsX3;

		        //for (ic=0; ic<icMax; ic++) {
			    for (ic=0; ic<4; ic++) {

		            err = 0;
		        	if (crs1SArr[ic] > M_2PI)
		        	{
		        		crs1SArr[ic] -= M_2PI;
		        	}
			        sprintf(testname,"TEST%-d-%-d-%-d-%-d",testNum,ir+1,iz+1,ic+1);

			        center.latitude = 0.0;
			        center.longitude = 0.0;
			        start.latitude = 0.0;
			        start.longitude = 0.0;
			        end.latitude = 0.0;
			        end.longitude = 0.0;
			        //dir = 0;
			        err |= arcTanToTwoGeos(pt1Arr[ic],crs1SArr[ic],pt3Arr[ic],crs3EArr[ic],r1,
			                                          &center,&start,&end,&dir, TOL, EPS);

			        if ( err )
			        {
			            printf("%s - Error occurred in arcTanToTwoGeos err=0x%lx\n",testname,err);
			            errorCount++;
			            failedCount++;
			            testCaseCount++;
			            continue;
			        }

			        err |= invDist(center,centerExp,&centerDistErr,EPS);
			        err |= invDist(start,startExp,&startDistErr,EPS);
			        err |= invDist(end,endExp,&endDistErr,EPS);

			        if ( getMaskedError(err, getMaskAll()) )
			        {
			            printf("%s - Error occurred in post-arcTanToTwoGeos err=0x%lx\n",testname,err);
			            errorCount++;
			            failedCount++;
			            testCaseCount++;
			            continue;
			        }

                                if (ptsAreSame(center, centerExp, TESTTOL) && ptsAreSame(start, startExp, TESTTOL) && ptsAreSame(end, endExp, TESTTOL) && (dirExp == dir))
			        {
			        	if (PRINT_PASSED_CASES) printf("%s passed  dirExp= %d dir= %d  centerDistErr=%16.8e  startDistErr=%16.8e  endDistErr=%16.8e err=%ld\n",
				               testname,dirExp,dir,centerDistErr,startDistErr,endDistErr,err);
			        	if (PRINT_PASSED_CASES) printf("%s Center: %14.8f %14.8f  az1=%14.8f dist1S=%14.8f dist3E=%14.8f\n",
			                   testname,latxExp,lonxExp,az1,dist1S,dist3E);
			        	if (PRINT_PASSED_CASES) printf("%s r1=%14.8f  crsCS=%14.8f  crsCE=%14.8f\n",
				               testname,r1,crsCS/DEG2RAD,crsCE/DEG2RAD);
			        	if (PRINT_PASSED_CASES) printLLPoint("   startExp ",startExp);
			        	if (PRINT_PASSED_CASES) printLLPoint("   centerExp",centerExp);
			        	if (PRINT_PASSED_CASES) printLLPoint("   endExp   ",endExp);
			        	if (PRINT_PASSED_CASES) printLLPoint("   pt1      ",pt1);
			        	if (PRINT_PASSED_CASES) printLLPoint("   pt3      ",pt3);
			        	if (PRINT_PASSED_CASES) printLLPoint("   intx     ",intx);
			        	if (PRINT_PASSED_CASES) printf("%s crs1S=%14.8f  crs3E=%14.8f  crsX1=%14.8f  crsX3=%14.8f\n",
				               testname,crs1S/DEG2RAD,crs3E/DEG2RAD,crsX1/DEG2RAD,crsX3/DEG2RAD);
			        	if (PRINT_PASSED_CASES) printf("%s crs1SArr=%14.8f  crs3EArr=%14.8f\n",
				               testname,crs1SArr[ic]/DEG2RAD,crs3EArr[ic]/DEG2RAD);
			        	if (PRINT_PASSED_CASES) printLLPoint("   pt1Arr   ",pt1Arr[ic]);
			        	if (PRINT_PASSED_CASES) printLLPoint("   pt3Arr   ",pt3Arr[ic]);
			        	if (PRINT_PASSED_CASES) printLLPoint("   start    ",start);
			        	if (PRINT_PASSED_CASES) printLLPoint("   center   ",center);
			        	if (PRINT_PASSED_CASES) printLLPoint("   end      ",end);


			            passedCount++;
			        }
			        else
			        {
				        printf("%s failed  dirExp= %d dir= %d  centerDistErr=%16.8e  startDistErr=%16.8e  endDistErr=%16.8e err=%ld\n",
				               testname,dirExp,dir,centerDistErr,startDistErr,endDistErr,err);
			            printf("%s Center: %.20lf %.20lf  az1=%.20lf dist1S=%.20lf dist3E=%.20lf\n",
			                   testname,latxExp,lonxExp,az1,dist1S,dist3E);
				        printf("%s r1=%.20lf  crsCS=%.20lf  crsCE=%.20lf\n",
				               testname,r1,crsCS/DEG2RAD,crsCE/DEG2RAD);
				        printLLPoint("   startExp ",startExp);
				        printLLPoint("   centerExp",centerExp);
				        printLLPoint("   endExp   ",endExp);
				        printLLPoint("   pt1      ",pt1);
				        printLLPoint("   pt3      ",pt3);
				        printLLPoint("   intx     ",intx);
				        printf("%s crs1S=%.20lf  crs3E=%.20lf  crsX1=%.20lf  crsX3=%.20lf\n",
				               testname,crs1S/DEG2RAD,crs3E/DEG2RAD,crsX1/DEG2RAD,crsX3/DEG2RAD);
				        printf("%s crs1SArr=%.20lf  crs3EArr=%.20lf\n",
				               testname,crs1SArr[ic]/DEG2RAD,crs3EArr[ic]/DEG2RAD);
				        printLLPoint("   pt1Arr   ",pt1Arr[ic]);
				        printLLPoint("   pt3Arr   ",pt3Arr[ic]);
				        printLLPoint("   start    ",start);
				        printLLPoint("   center   ",center);
				        printLLPoint("   end      ",end);


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

    printf("\nFinish testArcTanToTwoGeos_Set2\n\n\n");

    return set;

}

/*
 * NAME: testArcTanToTwoGeos_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the arcTanToTwoGeos function.
 *
 * 		This function runs all the test cases from all the test sets for testing
 *      the above function and returns the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testArcTanToTwoGeos_AllSets(TestSuite) - A test suite with the folling metrics:
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
TestSuite testArcTanToTwoGeos_AllSets()
{
	TestSuite suite;
	TestSet set1;
	TestSet set2;

    printf("\nStart testArcTanToTwoGeos_AllSets\n");

    suite = newTestSuite("testArcTanToTwoGeos_AllSets");

    set1 = testArcTanToTwoGeos_Set1();
    addTestSet(set1,&suite);

    set2 = testArcTanToTwoGeos_Set2();
    addTestSet(set2,&suite);

    displayTestSuite(suite);

    printf("\nFinish testArcTanToTwoGeos_AllSets\n\n\n");

    return suite;
}

/*
 * NAME: testPtIsOnGeo_Set1
 *
 * DESCRIPTION:
 * 		This function is used to test the ptIsOnGeo function.
 *
 * 		This function runs the test data created by Balakrishna Babu and uses
 * 		his original test plan to execute the actual tests.
 *
 *     	The approach for testing ptIsOnGeo is as follows:
 *     	1. Select a start point, and a course and a distance to the end point.
 *     	2. Compute the end point using direct.
 *     	3. Compute test points along the line using distances of
 *     	   [-1000,-10,-NMTOL,0,0.1*dist,0.5*dist,0.8*dist,dist,dist+NMTOL,dist+10,dist+1000] from start.
 *     	4. Test each of these points with LengthCode equal to 0,1,2.
 *     	5. For LengthCode 0, [-10,dist+10,dist+1000] should result in 0 (not on line) and the rest
 *     	   should result in 1 (on line)
 *     	6. For LengthCode 1, [-10] should result in 0 and the rest 1.
 *     	7. For LengthCode 2, all should result in 1.
 *     	8. For each of the test points, compute the course to the end point. Add 90 deg to find the
 *     	   perpendicular course and find points -1000,-10,-NMTOL,NMTOL,10,1000 nm along this course.
 *     	9. Repeat the test for each of these points. [-1000,-10,10,1000] should result in 0.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testPtIsOnGeo_Set1(TestSet) - A test set with the folling metrics:
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

TestSet testPtIsOnGeo_Set1()
{
	double DEG2RAD = M_PI / 180.0;
	double NMTOL = 0.03 / 100.0 / 1852.0;  //0.03 cm or ~1.62e-7 nm


    int testNum=0;
    char testname[80];
    double latStart, lonStart, latEnd, lonEnd, latTest, lonTest;
    double az1, d1;
    double d0[11] = {-1000.0,-10.0,-NMTOL,0.0,0.1,0.5,0.8,1.0,NMTOL,10.0,1000.0};
    double d[11];
    double d90[7] = {0.0,-1000.0,-10.0,-NMTOL,NMTOL,10.0,1000.0};
    int onLineExp[3][11][7] = {
                                  {		  {0,0,0,0,0,0,0},   {0,0,0,0,0,0,0},   {0,0,0,0,0,0,0},
                                		  {1,0,0,0,0,0,0},   {1,0,0,0,0,0,0},   {1,0,0,0,0,0,0},   {1,0,0,0,0,0,0},
                                		  {1,0,0,0,0,0,0},   {0,0,0,0,0,0,0},   {0,0,0,0,0,0,0},   {0,0,0,0,0,0,0} },
                                  {		  {0,0,0,0,0,0,0},   {0,0,0,0,0,0,0},   {0,0,0,0,0,0,0},
                                		  {1,0,0,0,0,0,0},   {1,0,0,0,0,0,0},   {1,0,0,0,0,0,0},   {1,0,0,0,0,0,0},
                                		  {1,0,0,0,0,0,0},   {1,0,0,0,0,0,0},   {1,0,0,0,0,0,0},   {1,0,0,0,0,0,0} },
                                  {		  {1,0,0,0,0,0,0},   {1,0,0,0,0,0,0},   {1,0,0,0,0,0,0},
                                		  {1,0,0,0,0,0,0},   {1,0,0,0,0,0,0},   {1,0,0,0,0,0,0},   {1,0,0,0,0,0,0},
                                		  {1,0,0,0,0,0,0},   {1,0,0,0,0,0,0},   {1,0,0,0,0,0,0},   {1,0,0,0,0,0,0} } };
    int onLine = 0;
    int i, j, lengthCode;

    LLPoint start, end, testBase, test[10];
    double crsSE, crsES, crstemp, disttemp;
    double crsTE, crsET, distTE, crsTE90;
    int passedCount=0, failedCount=0;
    int errorCount = 0;
    int unverifiedCount = 0;
    int setupFailureCount = 0;
    int testCaseCount = 0;
    ErrorSet err = 0;

    TestSet set;

    long newSeed=20080110;

    printf("\nStart testPtIsOnGeo_Set1\n");

    set = newTestSet("testPtIsOnGeo_Set1");

    srand(newSeed);  //Initialize the random number generator

    while (testNum<100)
    {
        //Select start point, azimuth and distance to end point
        latStart = randLat();
        lonStart = randLon();
        az1 = randAzimuth();
        d1 = 0.5 * randDist();  //about 2700 nm
        testNum++;

        sprintf(testname,"TEST%-d",testNum);

        start.latitude = latStart * DEG2RAD;
        start.longitude = lonStart * DEG2RAD;
        crsSE = az1 * DEG2RAD;

        //Compute end point
        err |= direct(start, crsSE, d1, &end, EPS);

        latEnd = end.latitude / DEG2RAD;
        lonEnd = end.longitude / DEG2RAD;

        //Set distances along geodesic
        d[0] = d0[0];
        d[1] = d0[1];
        d[2] = d0[2];

        d[3] = d1 * d0[3];
        d[4] = d1 * d0[4];
        d[5] = d1 * d0[5];
        d[6] = d1 * d0[6];
        d[7] = d1 * d0[7];

        d[8] = d1 + d0[8];
        d[9] = d1 + d0[9];
        d[10] = d1 + d0[10];

        //For each point along geodesic
        for (i=0; i<11; i++)
        {

            err |= direct(start,crsSE,d[i],&testBase,EPS);

            if (i == 3) //testBase is start
              crsTE = crsSE;
            else if (i == 7) //testBase is end
            {
              err |= inverse(testBase,start,&crsES,&crstemp,&disttemp,EPS);
              crsTE = crsES + M_PI;
            }
            else
              err |= inverse(testBase,end,&crsTE,&crsET,&distTE,EPS);
            crsTE90 = crsTE + M_PI_2;
            if (crsTE90 >= M_2PI)
            {
                crsTE90 -= M_2PI;
            }

            //For each point along perpendicular course
            for (j=0; j<7; j++)
            {

                err |= direct(testBase,crsTE90,d90[j],&test[j],EPS);
                latTest = test[j].latitude / DEG2RAD;
                lonTest = test[j].longitude / DEG2RAD;

                for (lengthCode=0; lengthCode<3; lengthCode++)
                {

                    if( getMaskedError(err, getMaskAll()) )
                    {
                        printf("Error occurred in pre-ptIsOnGeo err=0x%lx i=%d j=%d lc=%d\n",
                               err,i,j,lengthCode);
                        setupFailureCount++;
                        errorCount++;
                        testCaseCount++;
                        continue;
                    }

                    onLine = ptIsOnGeo(start,end,test[j],static_cast<LineType>(lengthCode),&err,TOL,EPS);

                    if( getMaskedError(err, getMaskAll()) )
                    {
                        printf("Error occurred in ptIsOnGeo err=0x%lx i=%d j=%d lc=%d\n",
                               err,i,j,lengthCode);
                        failedCount++;
                        errorCount++;
                        testCaseCount++;
                        continue;
                    }

                    if ( onLine == onLineExp[lengthCode][i][j] )
                    {
                    	if (PRINT_PASSED_CASES) printf("%s Start: %14.8f %14.8f  az1=%14.8f d1=%14.8f\n",
                               testname, latStart, lonStart, az1, d1);
                    	if (PRINT_PASSED_CASES) printf("%s End: %14.8f %14.8f\n",testname, latEnd, lonEnd);

                    	if (PRINT_PASSED_CASES) printf("%s passed lengthCode=%d  test: %14.8f %14.8f  dist=%14.8f d90=%14.8f  onLine=%d\n",
                               testname,lengthCode,latTest,lonTest,d[i],d90[j],onLine);
                        passedCount++;
                    }
                    else
                    {
                        printf("%s Start: %14.8f %14.8f  az1=%14.8f d1=%14.8f\n",
                               testname, latStart, lonStart, az1, d1);
                        printf("%s End: %14.8f %14.8f\n",testname, latEnd, lonEnd);

                        printf("%s failed lengthCode=%d  test: %14.8f %14.8f  dist=%14.8f d90=%14.8f  onLine=%d\n",
                               testname,lengthCode,latTest,lonTest,d[i],d90[j],onLine);
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

    printf("\nFinish testPtIsOnGeo_Set1\n\n\n");

    return set;

}

/*
 * NAME: testPtIsOnGeo_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the ptIsOnGeo function.
 *
 * 		This function runs all the test cases from all the test sets for testing
 *      the above function and returns the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testPtIsOnGeo_AllSets(TestSuite) - A test suite with the folling metrics:
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
TestSuite testPtIsOnGeo_AllSets()
{
	TestSuite suite;
	TestSet set1;

    printf("\nStart testPtIsOnGeo_AllSets\n");

    suite = newTestSuite("testPtIsOnGeo_AllSets");

    set1 = testPtIsOnGeo_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("\nFinish testPtIsOnGeo_AllSets\n\n\n");

    return suite;
}

/*
 * NAME: testGeoCrs_Set1
 *
 * DESCRIPTION:
 * 		This function is used to test the geoCrs function.
 *
 * 		This function runs the test data created by Balakrishna Babu and uses
 * 		his original test plan to execute the actual tests.
 *
 *     	The approach for testing geoCrs is as follows:
 *     	1. Select a start point, and a course and a distance to the end point.
 *     	2. Compute the end point using direct and reverse course using inverse.
 *     	3. Compute test points along the line using distances of
 *     	   [-1000,-10,-NMTOL,0,0.1*dist,0.5*dist,0.8*dist,dist,dist+NMTOL,dist+10,dist+1000] from start.
 *     	4. For each test point, compute the course and distance to both start and end points
 *     	   using inverse.
 *     	5. Forward course at the test point is the course to the end point if the test point is before
 *     	   the end point, and 180+course to start point if the test point is after the end point.
 *     	6. Test each of these points with LengthCode equal to 0,1,2.
 *     	7. For LengthCode 0, [-1000,-10,-NMTOL,dist+NMTOL,dist+10,dist+1000] should result in invalid course and the rest
 *     	   should result in values that match values calculated in step 5.
 *     	8. For LengthCode 1, [-1000,-10,-NMTOL] should result in invalid course and the rest should provide valid results.
 *     	9. For LengthCode 2, all test points should provide valid results and match valus from step 5.
 *     	10. For each of the test points, compute the course to the end point. Add 90 deg to find the
 *     	   perpendicular course and find points -10,-NMTOL,0,NMTOL,10 nm along this course.
 *     	11. Repeat the test for each of these points. [-10,-NMTOL,NMTOL,10] should result in invalid course.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testGeoCrs_Set1(TestSet) - A test set with the folling metrics:
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
TestSet testGeoCrs_Set1()
{
	double DEG2RAD = M_PI / 180.0;
	double NMTOL = 0.03 / 100.0 / 1852.0;  //0.03 cm or ~1.62e-7 nm
	double AZDEGTOL = 5.0e-7; //5e-7 deg or 8.72665e-9 rad or 0.162 cm at 100 nm

    int testNum=0;
    char testname[80];
    double latStart, lonStart, latEnd, lonEnd, latTest, lonTest;
    double az1, d1, d12;
    double d0[11] = {-1000.0,-10.0,-NMTOL,0.0,0.1,0.5,0.8,1.0,NMTOL,10.0,1000.0};
    double d[11];
    double d90[5] = {0.0,-10.0,-NMTOL,NMTOL,10.0};
    int onLineExp[3][11][5] = {
                                  {		  {0,0,0,0,0},   {0,0,0,0,0},   {0,0,0,0,0},
                                		  {1,0,0,0,0},   {1,0,0,0,0},   {1,0,0,0,0},   {1,0,0,0,0},
                                		  {1,0,0,0,0},   {0,0,0,0,0},   {0,0,0,0,0},   {0,0,0,0,0} },
                                  {		  {0,0,0,0,0},   {0,0,0,0,0},   {0,0,0,0,0},
                                		  {1,0,0,0,0},   {1,0,0,0,0},   {1,0,0,0,0},   {1,0,0,0,0},
                                		  {1,0,0,0,0},   {1,0,0,0,0},   {1,0,0,0,0},   {1,0,0,0,0} },
                                  {	      {1,0,0,0,0},   {1,0,0,0,0},   {1,0,0,0,0},
                                		  {1,0,0,0,0},   {1,0,0,0,0},   {1,0,0,0,0},   {1,0,0,0,0},
                                		  {1,0,0,0,0},   {1,0,0,0,0},   {1,0,0,0,0},   {1,0,0,0,0} } };
    int onLine = 0;
    int i, j;
    int lengthCode;

    LLPoint start, end, testBase, test[10];
    double crsSE, crsES, crsSE2;
    double crsTS, crsST, distST;
    double crsTE, crsET, distTE, crsTE90;
    double startCrs, revCrs, distToPt, crsAtPt;
    int passedCount=0, failedCount=0;
    int errorCount = 0;
    int unverifiedCount = 0;
    int setupFailureCount = 0;
    int testCaseCount = 0;
    ErrorSet err = 0;
    Geodesic tempGeo;

    TestSet set;

    long newSeed=20080117;

    printf("\nStart testGeoCrs_Set1\n");

    set = newTestSet("testGeoCrs_Set1");

    srand(newSeed);  //Initialize the random number generator

    while (testNum<100)
    {
        //Select start point, azimuth and distance to end point
        latStart = randLat();
        lonStart = randLon();
        az1 = randAzimuth();
        d1 = 0.5 * randDist();  //about 2700 nm
        testNum++;

        sprintf(testname,"TEST%-d",testNum);

        start.latitude = latStart * DEG2RAD;
        start.longitude = lonStart * DEG2RAD;
        crsSE = az1 * DEG2RAD;

        //Compute end point
        err |= direct(start, crsSE, d1, &end, EPS);
        err |= inverse(end,start,&crsES,&crsSE2,&d12,EPS);

        latEnd = end.latitude / DEG2RAD;
        lonEnd = end.longitude / DEG2RAD;

        //Set distances along geodesic
        d[0] = d0[0];
        d[1] = d0[1];
        d[2] = d0[2];

        d[3] = d1 * d0[3];
        d[4] = d1 * d0[4];
        d[5] = d1 * d0[5];
        d[6] = d1 * d0[6];
        d[7] = d1 * d0[7];

        d[8] = d1 + d0[8];
        d[9] = d1 + d0[9];
        d[10] = d1 + d0[10];

        //For each point along geodesic
        for (i=0; i<11; i++)
        {

            err |= direct(start,crsSE,d[i],&testBase,EPS);

            if (i == 3) //testBase is start
              crsTE = crsSE;
            else if (i == 7) //testBase is end
              crsTE = crsES + M_PI;
            else
              err |= inverse(testBase,end,&crsTE,&crsET,&distTE,EPS);
            crsTE90 = crsTE + M_PI_2;
            if (crsTE90 >= M_2PI)
            {
                crsTE90 -= M_2PI;
            }

            //For each point along perpendicular course
            for (j=0; j<5; j++)
            {

                err |= direct(testBase,crsTE90,d90[j],&test[j],EPS);
                latTest = test[j].latitude / DEG2RAD;
                lonTest = test[j].longitude / DEG2RAD;

                //Calculate expected values for the test point
                err |= inverse(test[j],start,&crsTS,&crsST,&distST,EPS);
                err |= inverse(test[j],end,&crsTE,&crsET,&distTE,EPS);
                if (distST > distTE)
                {  //test point is after end
                    crsTE = crsTS + M_PI;
                    if (crsTE > M_2PI)
                    {
                        crsTE -= M_2PI;
                    }
                }

                for (lengthCode=0; lengthCode<3; lengthCode++)
                {
                	err |= createGeo(&tempGeo, start, end, static_cast<LineType>(lengthCode), EPS);

                    if ( getMaskedError(err, getMaskAll()) )
                    {
                        printf("Error occurred in pre-geoCrs err=0x%lx i=%d j=%d lc=%d\n",
                               err,i,j,lengthCode);
                        setupFailureCount++;
                        errorCount++;
                        testCaseCount++;
                        continue;
                    }

                	crsAtPt = geoCrs(tempGeo, test[j], &startCrs, &revCrs, &distToPt, &err, TOL, EPS);

                    if ( getMaskedError(err, getMaskAll()) )
                    {
                        printf("Error occurred in geoCrs err=0x%lx i=%d j=%d lc=%d\n",
                               err,i,j,lengthCode);
                        failedCount++;
                        errorCount++;
                        testCaseCount++;
                        continue;
                    }

                    onLine = onLineExp[lengthCode][i][j];
                    if ( onLineExp[lengthCode][i][j] == 1 )
                    {
                        if ( //(fabs(startCrs-crsSE) < AZDEGTOL*DEG2RAD) &&
                            //(fabs(revCrs-crsES) < AZDEGTOL*DEG2RAD) &&
                            //(fabs(distToPt-distST) < NMTOL) &&
                            (fabs(crsAtPt-crsTE) < AZDEGTOL*DEG2RAD) )
                        {
                        	if (PRINT_PASSED_CASES) printf("%s Start: %14.8f %14.8f  az1=%14.8f d1=%14.8f\n",
                                   testname, latStart, lonStart, az1, d1);
                        	if (PRINT_PASSED_CASES) printf("%s End: %14.8f %14.8f  crsSE2=%14.8f  crsES=%14.8f  d12=%14.8f \n",
                                   testname, latEnd, lonEnd, crsSE2/DEG2RAD, crsES/DEG2RAD, d12);

                        	if (PRINT_PASSED_CASES) printf("%s passed lc=%d  test: %14.8f %14.8f  dist=%14.8f d90=%14.8f OL=%d\n",
                                   testname,lengthCode,latTest,lonTest,d[i],d90[j],onLine);
                        	if (PRINT_PASSED_CASES) printf("%s passed crsAtPt=%14.8f crsTE=%14.8f  distToPt=%14.8f  distST=%14.8f\n",
                                   testname,crsAtPt/DEG2RAD,crsTE/DEG2RAD,distToPt,distST);
                            passedCount++;
                        }
                        else
                        {
                            printf("%s Start: %14.8f %14.8f  az1=%14.8f d1=%14.8f\n",
                                   testname, latStart, lonStart, az1, d1);
                            printf("%s End: %14.8f %14.8f  crsSE2=%14.8f  crsES=%14.8f  d12=%14.8f \n",
                                   testname, latEnd, lonEnd, crsSE2/DEG2RAD, crsES/DEG2RAD, d12);

                            printf("%s failed lc=%d  test: %14.8f %14.8f  dist=%14.8f d90=%14.8f OL=%d\n",
                                   testname,lengthCode,latTest,lonTest,d[i],d90[j],onLine);
                            printf("%s failed crsAtPt=%14.8f crsTE=%14.8f  distToPt=%14.8f  distST=%14.8f\n",
                                   testname,crsAtPt/DEG2RAD,crsTE/DEG2RAD,distToPt,distST);
                            failedCount++;
                        }
                    }
                    else if ( onLineExp[lengthCode][i][j] == 0 )
                    {
                        if (crsAtPt == -1.0)
                        {
                        	if (PRINT_PASSED_CASES) printf("%s Start: %14.8f %14.8f  az1=%14.8f d1=%14.8f\n",
                                   testname, latStart, lonStart, az1, d1);
                        	if (PRINT_PASSED_CASES) printf("%s End: %14.8f %14.8f  crsSE2=%14.8f  crsES=%14.8f  d12=%14.8f \n",
                                   testname, latEnd, lonEnd, crsSE2/DEG2RAD, crsES/DEG2RAD, d12);

                        	if (PRINT_PASSED_CASES) printf("%s passed lc=%d  test: %14.8f %14.8f  dist=%14.8f d90=%14.8f OL=%d\n",
                                   testname,lengthCode,latTest,lonTest,d[i],d90[j],onLine);
                            passedCount++;
                        }
                        else
                        {
                            printf("%s Start: %14.8f %14.8f  az1=%14.8f d1=%14.8f\n",
                                   testname, latStart, lonStart, az1, d1);
                            printf("%s End: %14.8f %14.8f  crsSE2=%14.8f  crsES=%14.8f  d12=%14.8f \n",
                                   testname, latEnd, lonEnd, crsSE2/DEG2RAD, crsES/DEG2RAD, d12);

                            printf("%s failed lc=%d  test: %14.8f %14.8f  dist=%14.8f d90=%14.8f OL=%d\n",
                                   testname,lengthCode,latTest,lonTest,d[i],d90[j],onLine);
                            failedCount++;
                        }
                    }
                    else
                    {
                        printf("%s Start: %14.8f %14.8f  az1=%14.8f d1=%14.8f\n",
                               testname, latStart, lonStart, az1, d1);
                        printf("%s End: %14.8f %14.8f  crsSE2=%14.8f  crsES=%14.8f  d12=%14.8f \n",
                               testname, latEnd, lonEnd, crsSE2/DEG2RAD, crsES/DEG2RAD, d12);

                        printf("%s failed PROBLEM WITH TEST lc=%d  test: %14.8f %14.8f  dist=%14.8f d90=%14.8f\n",
                               testname,lengthCode,latTest,lonTest,d[i],d90[j]);
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

    printf("\nFinish testGeoCrs_Set1\n\n\n");

    return set;
}

/*
 * NAME: testGeoCrs_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the geoCrs function.
 *
 * 		This function runs all the test cases from all the test sets for testing
 *      the above function and returns the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testGeoCrs_AllSets(TestSuite) - A test suite with the folling metrics:
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
TestSuite testGeoCrs_AllSets()
{
	TestSuite suite;
	TestSet set1;

    printf("\nStart testGeoCrs_AllSets\n");

    suite = newTestSuite("testGeoCrs_AllSets");

    set1 = testGeoCrs_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("\nFinish testGeoCrs_AllSets\n\n\n");

    return suite;
}

/*
 * NAME: testProjectArcTanPtsToGeo_Set1
 *
 * DESCRIPTION:
 *
 * 		This function runs the test data created by Balakrishna Babu and uses
 * 		his test plan to execute the actual tests.
 *
 *      The approach for testing projectArcTanPtsToGeo is as follows:
 *
 *      1. Select a center point (C) and radius for Arc.
 *      2. Select an azimuth and distance from C to compute point (P) and select
 *         a course for the geodesic at P.
 *      3. Compute the perpendicular intercept on geodesic from center using
 *         projectToGeo.  Compute two points on the Arc (T1,T2) which are
 *         90 degrees from the perpendicular intercept line.
 *      4. Compute perpendicular intercept from T1 to the geodesic using
 *         projectToGeo (L1 on the geodesic).
 *      5. Compute the intersection points (U,V) between geodesic T1L1 and the Arc
 *         using geoArcIntx.
 *      6. Calculate the mid point of U and V (say W).
 *      7. Calculate a point (new T1) on the Arc along a geodesic CW.
 *      8. Repeat steps 4 through 7 until point T1 converges.  This will be one tangent point.
 *         Corresponding L1 will be the point on the geodesic.
 *      9. Repeat the process for T2.
 *      10.	The results from the function should match the points (T1, T2, L1 and L2)
 *         from steps 8 & 9.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testProjectArcTanPtsToGeo_Set1(TestSet) - A test set with the folling metrics:
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

TestSet testProjectArcTanPtsToGeo_Set1()
{

    double DEG2RAD = M_PI / 180.0;
    double NMTOL = 0.03 / 100.0 / 1852.0;  //0.03 cm or ~1.62e-7 nm

    int testCaseCount = 0, setupFailureCount = 0, errorCount = 0, unverifiedCount = 0;
    char testname[80];
    double latC, lonC, latP, lonP;
    double r1;
    double azCP, distCP, azLineP;

    double crsCP, crsLineP; //course from center to P, course of geodesic at P
    LLPoint p1, c1; //Point on geodesic and center of Arc
    LLPoint tanP1; //Tangent points on Arc from lines perpendicular to geodesic
    LLPoint lineP1; //Projection of tangent points on geodesic
    LLPoint lineP3; //Projection of C on geodesic through P
    LLPointPair lineP, tanP;//Solution from function under test
    LLPointPair intx;
    LLPoint midP, tanP1new;

    int nX;
    double crsC3, distC3;
    double latTP1, lonTP1, latTP2, lonTP2;
    double latLP1, lonLP1, latLP2, lonLP2;
    double crsCTP1, crsCTP2; //course from center to initial guess tangent points
    double crsTL, distTL; //course and distance from tangent to line
    double crs12, crs21, dist12;
    double crsCT, crsTC, distCT;

    double distConverge;
    int converged = 0, nIter;
    int passedCount = 0, failedCount = 0;
    ErrorSet err = 0;

    TestSet set;

    long newSeed = 20080124;

    printf("Start testProjectArcTanPtsToGeo_Set1\n");

    set = newTestSet("testProjectArcTanPtsToGeo_Set1");

    srand(newSeed); //Initialize the random number generator

    while (testCaseCount < 10000)
    {
        err = 0;

        //Select center, radius for Arc
        latC = randLat();
        lonC = randLon();
        r1 = 0.002 * randDist(); //about 10.80 nm

        //Select azimuth and diatance from center to locate P
        azCP = randAzimuth();
        distCP = 0.02 * randDist(); //about 108.0 nm

        //Select azimuth for geodesic at P
        azLineP = randAzimuth();

        sprintf(testname, "TEST%-d", testCaseCount);

        c1.latitude = latC * DEG2RAD;
        c1.longitude = lonC * DEG2RAD;
        crsCP = azCP * DEG2RAD;
        crsLineP = azLineP * DEG2RAD;

        //Calculate point P
        err |= direct(c1, crsCP, distCP, &p1, EPS);
        latP = p1.latitude / DEG2RAD;
        lonP = p1.longitude / DEG2RAD;

        err |= projectToGeo(p1, crsLineP, c1, &lineP3, &crsC3, &distC3,
                TOL, EPS);
        crsCTP1 = crsC3 + M_PI_2;
        crsCTP2 = crsC3 - M_PI_2;

        err |= direct(c1, crsCTP1, r1, &tanP1, EPS);

        converged = 0;
        nIter = 0;
        while (!converged)
        { //0=not converged, 1=converged

            //Calculate the intersection points between the Arc and the projection
            //geodesic from tangent point to the given geodesic through P.
            //There should be two intersection points T and U.
            err |= projectToGeo(p1, crsLineP, tanP1, &lineP1, &crsTL,
                    &distTL, TOL, EPS);
            err |= geoArcIntx(tanP1, crsTL, c1, r1, intx, &nX, TOL,
                    EPS);
            nIter++;

            if (nX == 1)
            {
                tanP1 = intx[0];
                converged = 1;
            }

            else
            {
                //Calculate the mid point of T and U (say V).
                err |= inverse(intx[0], intx[1], &crs12, &crs21, &dist12,
                        EPS);
                err |= direct(intx[0], crs12, 0.5 * dist12, &midP, EPS);

                //Compute the course from center of Arc to V and calculate the
                // new point (T) on the Arc.
                err |= inverse(c1, midP, &crsCT, &crsTC, &distCT, EPS);
                err |= direct(c1, crsCT, r1, &tanP1new, EPS);
                err |= invDist(tanP1, tanP1new, &distConverge, EPS);
                tanP1 = tanP1new;
                if (distConverge < NMTOL)
                {
                    converged = 1;
                }

                if (getMaskedError(err, getMaskAll()))
                {
                    break;
                }
            }
        }

        if (getMaskedError(err, getMaskAll()))
        {
            printf("Error occurred in pre-projectArcTanPtsToGeo err=0x%lx\n",
                    err);
            testCaseCount++;
            setupFailureCount++;
            errorCount++;
            continue;
        }

        //Initialize to zeros
        tanP[0].latitude = 0.0;
        tanP[0].longitude = 0.0;
        tanP[1].latitude = 0.0;
        tanP[1].longitude = 0.0;
        lineP[0].latitude = 0.0;
        lineP[0].longitude = 0.0;
        lineP[1].latitude = 0.0;
        lineP[1].longitude = 0.0;

        err |= projectArcTanPtsToGeo(p1, crsLineP, c1, r1, lineP, tanP, TOL,
                EPS);

        if (getMaskedError(err, getMaskAll()))
        {
            printf("Error occurred in projectArcTanPtsToGeo err=0x%lx\n", err);
            failedCount++;
            errorCount++;
            testCaseCount++;
            continue;
        }

        latTP1 = tanP[0].latitude / DEG2RAD;
        lonTP1 = tanP[0].longitude / DEG2RAD;
        latTP2 = tanP[1].latitude / DEG2RAD;
        lonTP2 = tanP[1].longitude / DEG2RAD;
        latLP1 = lineP[0].latitude / DEG2RAD;
        lonLP1 = lineP[0].longitude / DEG2RAD;
        latLP2 = lineP[1].latitude / DEG2RAD;
        lonLP2 = lineP[1].longitude / DEG2RAD;

        if (ptsAreSame(tanP1, tanP[0], TESTTOL) || ptsAreSame(tanP1, tanP[1], TESTTOL))
        {
        	if (PRINT_PASSED_CASES) printf("%s C: %14.8f %14.8f  r1: %14.8f  az: %14.8f  dist: %14.8f\n",
                    testname, latC, lonC, r1, azCP, distCP);
        	if (PRINT_PASSED_CASES) printf("%s P: %14.8f %14.8f  azLineP: %14.8f\n", testname, latP, lonP,
                    azLineP);

        	if (PRINT_PASSED_CASES) printf("%s passed tanP1: %14.8f %14.8f  tanP2: %14.8f %14.8f\n",
                    testname, latTP1, lonTP1, latTP2, lonTP2);
        	if (PRINT_PASSED_CASES) printf("%s passed lineP1: %14.8f %14.8f  lineP2: %14.8f %14.8f\n",
                    testname, latLP1, lonLP1, latLP2, lonLP2);
            passedCount++;
        }
        else
        {
            printf("%s C: %14.8f %14.8f  r1: %14.8f  az: %14.8f  dist: %14.8f\n",
                    testname, latC, lonC, r1, azCP, distCP);
            printf("%s P: %14.8f %14.8f  azLineP: %14.8f\n", testname, latP, lonP,
                    azLineP);

            printf("%s failed tanP1: %14.8f %14.8f  tanP2: %14.8f %14.8f\n",
                    testname, latTP1, lonTP1, latTP2, lonTP2);
            printf("%s failed lineP1: %14.8f %14.8f  lineP2: %14.8f %14.8f\n",
                    testname, latLP1, lonLP1, latLP2, lonLP2);
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

    printf("\nFinish testProjectArcTanPtsToGeo_Set1\n");

    return set;
}

/*
 * NAME: testProjectArcTanPtsToGeo_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the projectArcTanPtsToGeo function.
 *
 * 		This function runs all the test cases from all the test sets for testing
 *      the above function and returns the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testProjectArcTanPtsToGeo_AllSets(TestSuite) - A test suite with the folling metrics:
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
TestSuite testProjectArcTanPtsToGeo_AllSets()
{
	TestSuite suite;
	TestSet set1;

    printf("\nStart testProjectArcTanPtsToGeo_AllSets\n");

    suite = newTestSuite("testProjectArcTanPtsToGeo_AllSets");

    set1 = testProjectArcTanPtsToGeo_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("\nFinish testProjectArcTanPtsToGeo_AllSets\n\n\n");

    return suite;
}

/*
 * NAME: testLocusGeoIntx_Set1
 *
 * DESCRIPTION:
 * 		This function is used to test the locusGeoIntx function.
 *
 * 		This function runs the test data created by Balakrishna Babu and uses
 * 		his test plan to execute the actual tests.
 *
 *      1. Select an intersection point using random.
 *      2. Construct a locus that passes through the selected intersection point.
 *         2.1. Azimuth and a distance to find the point on the defining geodesic
 *         2.2. Compute azimuth to the point on locus from point on geodesic
 *         2.3. Compute two azimuths 90 degress on either side of azimuth to point
 *             and select two distances to locate the start and end points of the
 *             defining geodesic
 *         2.4. Select startDist and compute endDist
 *         2.5. Create locus that passed through the selected intersection point
 *      3. Construct a geodesic that passes through the same intersection point.
 *      4. Find intersection using all lineTypes for both locus and geodesic.
 *      5. Verify if NO_INTERSECTION_ERR is expected.
 *      6. Verify if the intescetion point lies on both locus and geodesic.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testLocusGeoIntx_Set1(TestSet) - A test set with the folling metrics:
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

TestSet testLocusGeoIntx_Set1()
{

    int testNum=0;
    char testname[256];
    double DEG2RAD = M_PI / 180.0;

    double latX, lonX;
    double az1, d1, d2, d3, slope1;
    double startDist1, endDist1;

    LLPoint proj, proj1;
    LLPoint geoStart1, geoEnd1;
    LLPoint geoStart2, geoEnd2;
    LineType lineType1, lineType2;  //0=SEGMENT,1=SEMIINFINITE,2=INFINITE
    Locus locus1;
    Locus* locusPtr1;

    double fcrs,bcrs;
    LLPoint intxExp, intx, intxINF;
    double crsC1, crsC2, crsC3;
    int expectValidIntx1, expectValidIntx2;
    int passedCount=0, failedCount=0;
    int testCaseCount = 0, setupFailureCount = 0, errorCount = 0, unverifiedCount = 0;
    ErrorSet err=0;

    TestSet set;

    long newSeed=20080501;

    printf("Start testLocusGeoIntx_Set1\n");

    set = newTestSet("testLocusGeoIntx_Set1");

	srand(newSeed);  //Initialize the random number generator

	locusPtr1 = &locus1;

	while (testNum < 1000)
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
    	d2 = 0.01*randDist();
    	d3 = 0.01*randDist();
        err |= direct(proj1, crsC2, d2, &geoStart1, EPS);
        err |= direct(proj1, crsC3, d3, &geoEnd1, EPS);

        //4. Select startDist and compute endDist
    	startDist1 = 0.001*randDist();
    	slope1 = (d1-startDist1)/d2;
    	endDist1 = d1 + d3 * (d1-startDist1)/d2;
    	lineType1 = INFINITE;

    	//5. Create locus that passed through the selected intersection point
        err |= createLocus(locusPtr1,geoStart1,geoEnd1,startDist1,endDist1,lineType1,TOL,EPS);
    	if (ptIsOnLocus(locus1,intxExp,&proj,&err,TOL,EPS)) {
    		expectValidIntx1 = 1;  //true
    	}
    	else {
    		expectValidIntx1 = 0;  //false
    	}

    	//Construct a geodesic that passes through the selected intersection point
        //Compute two azimuths and select two distances to locate the start and
        //end points of the geodesic
    	az1 = randAzimuth();
        crsC1 = az1 * DEG2RAD;

        crsC2 = crsC1 + M_PI_2;
        if (crsC2 > M_2PI) {
        	crsC2 -= M_2PI;
        }
        crsC3 = crsC1 - M_PI_2;
        if (crsC3 < 0.0) {
        	crsC3 += M_2PI;
        }
    	d2 = 0.01*randDist();
    	d3 = 0.01*randDist();
        err |= direct(intxExp, crsC2, d2, &geoStart2, EPS);
        err |= direct(intxExp, crsC3, d3, &geoEnd2, EPS);
    	lineType2 = INFINITE;

        if ( getMaskedError(err, getMaskAll()) )
        {
        	printf("Error occurred in pre-locusGeoIntx INFINITE err=0x%lx\n",err);
        	testCaseCount+=3;//This test function actually tests each line type which is a total of 3 (SEGMENT, SEMIINFINITE, and INFINITE)
            errorCount++;
            setupFailureCount++;
        	continue;
        }

        intxINF.latitude = 0.0;
        intxINF.longitude = 0.0;
        err |= locusGeoIntx(geoStart2,geoEnd2,locus1,&intxINF,TOL,EPS);

        if ( getMaskedError(err, getMaskAll()) )
        {
        	printf("Error occurred in locusGeoIntx INFINITE testNum %d err=0x%lx\n",testNum,err);
        	testCaseCount+=3;//This test function actually tests each line type which is a total of 3 (SEGMENT, SEMIINFINITE, and INFINITE)
        	failedCount++;
            errorCount++;
        	continue;
        }

        lineType2 = INFINITE;
    	if (ptIsOnGeo(geoStart2,geoEnd2,intxINF,lineType2,&err,TOL,EPS)) {
    		expectValidIntx2 = 1;  //true
    	}
    	else {
    		expectValidIntx2 = 0;  //false
    	}

        for (int locus1lineType=SEGMENT; locus1lineType<=INFINITE; locus1lineType++) {
         locus1.lineType = static_cast<LineType>(locus1lineType);
        	if (ptIsOnLocus(locus1,intxINF,&proj,&err,TOL,EPS)) {
        		expectValidIntx1 = 1;  //true
        	}
        	else {
        		expectValidIntx1 = 0;  //false
        	}

            if ( getMaskedError(err, getMaskAll()) )
            {
            	printf("Error occurred in pre-locusIntx testNum %d err= 0x%lx\n",testNum,err);
            	failedCount++;
                errorCount++;
                testCaseCount++;//only count once since we're inside the for loop that iterates through line types
            	continue;
            }

        	intx.latitude = 0.0;
            intx.longitude = 0.0;
            err |= locusGeoIntx(geoStart2,geoEnd2,locus1,&intx,TOL,EPS);

            if ( getMaskedError(err, getMaskAll()) )
            {
            	printf("Error occurred in locusIntx err= 0x%lx\n",err);
            	failedCount++;
                errorCount++;
                testCaseCount++;//only count once since we're inside the for loop that iterates through line types
            	continue;
            }

            if ( (expectValidIntx1 == 0) || (expectValidIntx2 == 0) ) {
            	if (err == NO_INTERSECTION_ERR) {
            		if (PRINT_PASSED_CASES) printf("%s latX=%14.8f lonX=%14.8f\n",
                      		testname,latX,lonX);
            		if (PRINT_PASSED_CASES) printLLPoint("Expected Intersection point",intxExp);
            		if (PRINT_PASSED_CASES) printLLPoint("point on defining geodesic",proj1);
            		if (PRINT_PASSED_CASES) printLLPoint("geoStart1",geoStart1);
            		if (PRINT_PASSED_CASES) printLLPoint("geoEnd1",geoEnd1);
            		if (PRINT_PASSED_CASES) printf("%s X: %14.8f %14.8f  az1=%14.8f d1=%14.8f d2=%14.8f d3=%14.8f startDist1=%14.8f endDist1=%14.8f slope1=%14.8f\n",
                      		testname,latX,lonX,az1,d1,d2,d3,startDist1,endDist1,slope1);
            		if (PRINT_PASSED_CASES) printf("%s crsC1=%14.8f crsC2=%14.8f crsC3=%14.8f\n",
                			testname,crsC1/DEG2RAD,crsC2/DEG2RAD,crsC3/DEG2RAD);
            		if (PRINT_PASSED_CASES) printLLPoint("locusStart1",locus1.locusStart);
            		if (PRINT_PASSED_CASES) printLLPoint("locusEnd1",locus1.locusEnd);
            		if (PRINT_PASSED_CASES) printf("%s lineType1=%d expIntx1=%d\n",
                			testname,locus1.lineType,expectValidIntx1);
            		if (PRINT_PASSED_CASES) printf("%s X: %14.8f %14.8f  crsC2=%14.8f crsC3=%14.8f d2=%14.8f d3=%14.8f\n",
                      		testname,latX,lonX,crsC2/DEG2RAD,crsC3/DEG2RAD,d2,d3);
            		if (PRINT_PASSED_CASES) printLLPoint("geoStart2",geoStart2);
            		if (PRINT_PASSED_CASES)  printLLPoint("geoEnd2",geoEnd2);
            		if (PRINT_PASSED_CASES)  printLLPoint("intxINF",intxINF);


            		if (PRINT_PASSED_CASES) printf("%s passed lineType1=%d lineType2=%d expIntx1=%d expIntx2=%d\n",
                			testname,locus1.lineType,lineType2,expectValidIntx1,expectValidIntx2);
                	passedCount++;
            	}
            	else {
                    printf("%s latX=%14.8f lonX=%14.8f\n",
                      		testname,latX,lonX);
                    printLLPoint("Expected Intersection point",intxExp);
                    printLLPoint("point on defining geodesic",proj1);
                    printLLPoint("geoStart1",geoStart1);
                    printLLPoint("geoEnd1",geoEnd1);
                    printf("%s X: %14.8f %14.8f  az1=%14.8f d1=%14.8f d2=%14.8f d3=%14.8f startDist1=%14.8f endDist1=%14.8f slope1=%14.8f\n",
                      		testname,latX,lonX,az1,d1,d2,d3,startDist1,endDist1,slope1);
                	printf("%s crsC1=%14.8f crsC2=%14.8f crsC3=%14.8f\n",
                			testname,crsC1/DEG2RAD,crsC2/DEG2RAD,crsC3/DEG2RAD);
                    printLLPoint("locusStart1",locus1.locusStart);
                    printLLPoint("locusEnd1",locus1.locusEnd);
                	printf("%s lineType1=%d expIntx1=%d\n",
                			testname,locus1.lineType,expectValidIntx1);
                	printf("%s X: %14.8f %14.8f  crsC2=%14.8f crsC3=%14.8f d2=%14.8f d3=%14.8f\n",
                      		testname,latX,lonX,crsC2/DEG2RAD,crsC3/DEG2RAD,d2,d3);
                    printLLPoint("geoStart2",geoStart2);
                    printLLPoint("geoEnd2",geoEnd2);
                    printLLPoint("intxINF",intxINF);

                	printf("%s failed lineType1=%d lineType2=%d expIntx1=%d expIntx2=%d\n",
                			testname,locus1.lineType,lineType2,expectValidIntx1,expectValidIntx2);
                	failedCount++;
            	}
            }
            else {
                if ( ptIsOnLocus(locus1,intx,&proj,&err,TESTTOL,EPS) &&
                     ptIsOnGeo(geoStart2,geoEnd2,intx,lineType2,&err,TESTTOL,EPS) ) {

                	if ( getMaskedError(err, getMaskAll()) )
                	{
                		printf("Error occurred in post-locusIntx err= 0x%lx\n",err);
                		failedCount++;
                        errorCount++;
                        testCaseCount++;//only count once since we're inside the for loop that iterates through line types
                		continue;
                	}
            		if (PRINT_PASSED_CASES) printf("%s latX=%14.8f lonX=%14.8f\n",
                      		testname,latX,lonX);
            		if (PRINT_PASSED_CASES) printLLPoint("Expected Intersection point",intxExp);
            		if (PRINT_PASSED_CASES) printLLPoint("point on defining geodesic",proj1);
            		if (PRINT_PASSED_CASES) printLLPoint("geoStart1",geoStart1);
            		if (PRINT_PASSED_CASES) printLLPoint("geoEnd1",geoEnd1);
            		if (PRINT_PASSED_CASES) printf("%s X: %14.8f %14.8f  az1=%14.8f d1=%14.8f d2=%14.8f d3=%14.8f startDist1=%14.8f endDist1=%14.8f slope1=%14.8f\n",
                      		testname,latX,lonX,az1,d1,d2,d3,startDist1,endDist1,slope1);
            		if (PRINT_PASSED_CASES) printf("%s crsC1=%14.8f crsC2=%14.8f crsC3=%14.8f\n",
                			testname,crsC1/DEG2RAD,crsC2/DEG2RAD,crsC3/DEG2RAD);
            		if (PRINT_PASSED_CASES) printLLPoint("locusStart1",locus1.locusStart);
            		if (PRINT_PASSED_CASES) printLLPoint("locusEnd1",locus1.locusEnd);
            		if (PRINT_PASSED_CASES) printf("%s lineType1=%d expIntx1=%d\n",
                			testname,locus1.lineType,expectValidIntx1);
            		if (PRINT_PASSED_CASES) printf("%s X: %14.8f %14.8f  crsC2=%14.8f crsC3=%14.8f d2=%14.8f d3=%14.8f\n",
                      		testname,latX,lonX,crsC2/DEG2RAD,crsC3/DEG2RAD,d2,d3);
            		if (PRINT_PASSED_CASES) printLLPoint("geoStart2",geoStart2);
            		if (PRINT_PASSED_CASES) printLLPoint("geoEnd2",geoEnd2);
            		if (PRINT_PASSED_CASES) printLLPoint("intxINF",intxINF);

                	if (PRINT_PASSED_CASES) printf("%s passed %14.8f %14.8f\n",
                			testname,intx.latitude/DEG2RAD,intx.longitude/DEG2RAD);
                	passedCount++;
                }
                else {
                    printf("%s latX=%14.8f lonX=%14.8f\n",
                      		testname,latX,lonX);
                    printLLPoint("Expected Intersection point",intxExp);
                    printLLPoint("point on defining geodesic",proj1);
                    printLLPoint("geoStart1",geoStart1);
                    printLLPoint("geoEnd1",geoEnd1);
                    printf("%s X: %14.8f %14.8f  az1=%14.8f d1=%14.8f d2=%14.8f d3=%14.8f startDist1=%14.8f endDist1=%14.8f slope1=%14.8f\n",
                      		testname,latX,lonX,az1,d1,d2,d3,startDist1,endDist1,slope1);
                	printf("%s crsC1=%14.8f crsC2=%14.8f crsC3=%14.8f\n",
                			testname,crsC1/DEG2RAD,crsC2/DEG2RAD,crsC3/DEG2RAD);
                    printLLPoint("locusStart1",locus1.locusStart);
                    printLLPoint("locusEnd1",locus1.locusEnd);
                	printf("%s lineType1=%d expIntx1=%d\n",
                			testname,locus1.lineType,expectValidIntx1);
                	printf("%s X: %14.8f %14.8f  crsC2=%14.8f crsC3=%14.8f d2=%14.8f d3=%14.8f\n",
                      		testname,latX,lonX,crsC2/DEG2RAD,crsC3/DEG2RAD,d2,d3);
                    printLLPoint("geoStart2",geoStart2);
                    printLLPoint("geoEnd2",geoEnd2);
                    printLLPoint("intxINF",intxINF);

                	printf("%s failed %14.8f %14.8f\n",
                			testname,intx.latitude/DEG2RAD,intx.longitude/DEG2RAD);
                	failedCount++;
                }
            }
            testCaseCount++;
        }
    }
    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("Finish testLocusGeoIntx_Set1\n");
    return set;
}

/*
 * NAME: testLocusGeoIntx_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the locusGeoIntx function.
 *
 * 		This function runs all the test cases from all the test sets for testing
 *      the above function and returns the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testLocusGeoIntx_AllSets(TestSuite) - A test suite with the folling metrics:
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
TestSuite testLocusGeoIntx_AllSets()
{
	TestSuite suite;
	TestSet set1;

    printf("\nStart testLocusGeoIntx_AllSets\n");

    suite = newTestSuite("testLocusGeoIntx_AllSets");

    set1 = testLocusGeoIntx_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("\nFinish testLocusGeoIntx_AllSets\n\n\n");

    return suite;
}


/*
 * NAME: testProjectToGeoAtAngle_Set1
 *
 * DESCRIPTION:
 * 		This function is used to test the projToGeoAtAngle function.
 *
 * 		This function runs the test data created by Balakrishna Babu and uses
 * 		his original test plan to execute the actual tests.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testProjectToGeoAtAngle_Set1(TestSet) - A test set with the folling metrics:
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
TestSet testProjectToGeoAtAngle_Set1()
{
	double DEG2RAD = M_PI / 180.0;

    int testNum=0;
    char testname[80];
    double latExp, lonExp;
    double az1, d1, az2, d2;
    LLPoint p3Exp;

    double latx, lonx;
    LLPoint p1, p2, p3;
    double crs13, crs31, dist13;
    double crs32, dist23;
    int passedCount=0, failedCount=0;
    int errorCount = 0;
    int unverifiedCount = 0;
    int setupFailureCount = 0;
    int testCaseCount = 0;
    double intAngle, angleMin = 0.01;
    ErrorSet err=0;

    TestSet set;

    long newSeed=12262007;

    printf("\nStart testProjectToGeoAtAngle_Set1\n");

    set = newTestSet("testProjectToGeoAtAngle_Set1");

    srand(newSeed);  //Initialize the random number generator

    while (testNum < 10000)
    {
        err = 0;

        latExp = randLat();
        lonExp = randLon();
        az1 = randAzimuth();
        d1 = 0.1*randDist();
        az2 = fmod(az1 + angleMin + (180.0 - 2*angleMin)*(double)rand()/RAND_MAX, 360.0);
        d2 = 0.1*randDist();
        testNum++;

        sprintf(testname,"TEST%-d",testNum);

        p3Exp.latitude = latExp * DEG2RAD;
        p3Exp.longitude = lonExp * DEG2RAD;
        crs31 = az1 * DEG2RAD;
        crs32 = az2 * DEG2RAD;

        err |= minSubtendedAngle(crs31,crs32,&intAngle);

        err |= direct(p3Exp, crs31, d1, &p1, EPS);
        err |= direct(p3Exp, crs32, d2, &p2, EPS);

        err |= inverse(p1, p3Exp, &crs13, &crs31, &dist13, EPS);

        if ( getMaskedError(err, getMaskAll()) )
        {
            printf("Error occurred in pre-projectToGeoAtAngle err=0x%lx\n",err);
            testCaseCount++;
            errorCount++;
            setupFailureCount++;
            continue;
        }

        err |= projectToGeoAtAngle(p1, crs13, p2, intAngle, &p3, &crs32, &dist23, TOL, EPS);

        if ( getMaskedError(err, getMaskAll()) )
        {
            printf("Error occurred in projectToGeoAtAngle err=0x%lx\n",err);
            printf("%s p3: %3.14lf %3.14lf  az1=%3.14lf d1=%3.14lf  az2=%3.14lf d2=%3.14lf\n",
                   testname, latExp, lonExp, az1, d1, az2, d2);
            printf("%s p1: %3.14lf =%3.14lf  p2: %3.14lf %3.14lf\n",
                   testname,p1.latitude/DEG2RAD,p1.longitude/DEG2RAD,p2.latitude/DEG2RAD,p2.longitude/DEG2RAD);

            printf("%s crs13=%3.14lf crs31=%3.14lf d13=%3.14lf crs32=%3.14lf d23=%3.14lf\n",
                   testname,crs13/DEG2RAD,crs31/DEG2RAD,dist13,crs32/DEG2RAD,dist23);
            printf("%s failed %3.14lf %3.14lf     %3.14lf %3.14lf     %3.14lf\n",
                   testname,latExp,latx,lonExp,lonx, intAngle/DEG2RAD);
            errorCount++;
            testCaseCount++;
            failedCount++;
            continue;
        }

        latx = p3.latitude / DEG2RAD;
        lonx = p3.longitude / DEG2RAD;

            if (ptsAreSame(p3, p3Exp, TESTTOL))
            {
            	if (PRINT_PASSED_CASES) printf("%s p3: %3.14lf %3.14lf  az1=%3.14lf d1=%3.14lf  az2=%3.14lf d2=%3.14lf\n",
                       testname, latExp, lonExp, az1, d1, az2, d2);
            	if (PRINT_PASSED_CASES) printf("%s p1: %3.14lf =%3.14lf  p2: %3.14lf %3.14lf\n",
                       testname,p1.latitude/DEG2RAD,p1.longitude/DEG2RAD,p2.latitude/DEG2RAD,p2.longitude/DEG2RAD);

            	if (PRINT_PASSED_CASES) printf("%s crs13=%3.14lf crs31=%3.14lf d13=%3.14lf crs32=%3.14lf d23=%3.14lf\n",
                       testname,crs13/DEG2RAD,crs31/DEG2RAD,dist13, crs32/DEG2RAD,dist23);
            	if (PRINT_PASSED_CASES) printf("%s passed %3.14lf %14.8e     %3.14lf %14.8e\n",
                       testname,latx,fabs(latExp-latx),lonx,fabs(lonExp-lonx));
                passedCount++;
            }
            else
            {
                printf("%s p3: %3.14lf %3.14lf  az1=%3.14lf d1=%3.14lf  az2=%3.14lf d2=%3.14lf\n",
                       testname, latExp, lonExp, az1, d1, az2, d2);
                printf("%s p1: %3.14lf =%3.14lf  p2: %3.14lf %3.14lf\n",
                       testname,p1.latitude/DEG2RAD,p1.longitude/DEG2RAD,p2.latitude/DEG2RAD,p2.longitude/DEG2RAD);
                printf("%s crs13=%3.14lf crs31=%3.14lf d13=%3.14lf crs32=%3.14lf d23=%3.14lf\n",
                       testname,crs13/DEG2RAD,crs31/DEG2RAD,dist13, crs32/DEG2RAD,dist23);
                printf("%s failed %3.14lf %3.14lf     %3.14lf %3.14lf     %3.14lf\n",
                       testname,latExp,latx,lonExp,lonx, intAngle/DEG2RAD);
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

    printf("\nFinish testProjectToGeoAtAngle_Set1\n\n\n");

    return set;

}

/*
 * NAME: testProjectToGeoAtAngle_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the projectToGeoAtAngle function.
 *
 * 		This function runs all the test cases from all the test sets for testing
 *      the above function and returns the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testProjectToGeoAtAngle_AllSets(TestSuite) - A test suite with the folling metrics:
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
TestSuite testProjectToGeoAtAngle_AllSets()
{
	TestSuite suite;
	TestSet set1;

    printf("\nStart testProjectToGeoAtAngle_AllSets\n");

    suite = newTestSuite("testProjectToGeoAtAngle_AllSets");

    set1 = testProjectToGeoAtAngle_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("\nFinish testProjectToGeoAtAngle_AllSets\n\n\n");

    return suite;
}

/*
 * NAME: testPtIsOnCrs_Set1
 *
 * DESCRIPTION:
 * 		This function is used to test the ptIsOnCrs function.
 *
 * 		This function runs the test data created by Rich Snow.
 *
 *     	The approach for testing ptIsOnCrs is as follows:
 *     	1. Select a start point, and a course and a distance to the end point.
 *     	2. Compute the end point using direct.
 *     	3. Compute test points along the line using distances of
 *     	   [-1000,-10,-NMTOL,0,0.1*dist,0.5*dist,0.8*dist,dist,dist+NMTOL,dist+10,dist+1000] from start.
 *     	4. For each of the test points, compute the course to the end point. Add 90 deg to find the
 *     	   perpendicular course and find points at dist 0,-1,-NMTOL,NMTOL,1 nm along this course.
 *     	5. Test each of these points. For distance 0 in 4., all the distances in 3. should return 1.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testPtIsOnCrs_Set1(TestSet) - A test set with the folling metrics:
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

TestSet testPtIsOnCrs_Set1()
{
	double DEG2RAD = M_PI / 180.0;
	double NMTOL = 0.03 / 100.0 / 1852.0;  //0.03 cm or ~1.62e-7 nm


    int testNum=0;
    char testname[80];
    double latStart, lonStart, latEnd, lonEnd, latTest, lonTest;
    double az1, d1;
    double d0[11] = {-1000.0,-10.0,-NMTOL,0.0,0.1,0.5,0.8,1.0,NMTOL,10.0,1000.0};
    double d[11];
    double d90[5] = {0.0,-1.0,-NMTOL,NMTOL,1.0};
    int onLineExp[11][5] = {
                                  {1,0,0,0,0},   {1,0,0,0,0},   {1,0,0,0,0},
                                  {1,0,0,0,0},   {1,0,0,0,0},   {1,0,0,0,0},   {1,0,0,0,0},
                                  {1,0,0,0,0},   {1,0,0,0,0},   {1,0,0,0,0},   {1,0,0,0,0} };
    int onLine = 0;
    int i, j;

    LLPoint start, end, testBase, test[10];
    double crsSE, crsTest, tempdbl, crsTest1, dist1Test;
    double crsTE, crsET, distTE, crsTE90;
    double crsES, crstemp, disttemp;
    int passedCount=0, failedCount=0;
    int errorCount = 0;
    int unverifiedCount = 0;
    int setupFailureCount = 0;
    int testCaseCount = 0;
    ErrorSet err = 0;

    TestSet set;

    long newSeed=20080110;

    printf("\nStart testPtIsOncrs_Set1\n");

    set = newTestSet("testPtIsOncrs_Set1");

    srand(newSeed);  //Initialize the random number generator

    while (testNum<100)
    {
        //Select start point, azimuth and distance to end point
        latStart = randLat();
        lonStart = randLon();
        az1 = randAzimuth();
        d1 = 0.5 * randDist();  //about 2700 nm
        testNum++;

        sprintf(testname,"TEST%-d",testNum);

        start.latitude = latStart * DEG2RAD;
        start.longitude = lonStart * DEG2RAD;
        crsSE = az1 * DEG2RAD;

        //Compute end point
        err |= direct(start, crsSE, d1, &end, EPS);

        latEnd = end.latitude / DEG2RAD;
        lonEnd = end.longitude / DEG2RAD;

        //Set distances along geodesic
        d[0] = d0[0];
        d[1] = d0[1];
        d[2] = d0[2];

        d[3] = d1 * d0[3];
        d[4] = d1 * d0[4];
        d[5] = d1 * d0[5];
        d[6] = d1 * d0[6];
        d[7] = d1 * d0[7];

        d[8] = d1 + d0[8];
        d[9] = d1 + d0[9];
        d[10] = d1 + d0[10];

        //For each point along geodesic
        for (i=0; i<11; i++)
        {

            err |= direct(start,crsSE,d[i],&testBase,EPS);

            if (i == 3) //testBase is start
              crsTE = crsSE;  
            else if (i == 7) //testBase is end
            {
              err |= inverse(testBase,start,&crsES,&crstemp,&disttemp,EPS);
              crsTE = crsES + M_PI;
            }
            else
              err |= inverse(testBase,end,&crsTE,&crsET,&distTE,EPS);
            crsTE90 = crsTE + M_PI_2;
            if (crsTE90 >= M_2PI)
            {
                crsTE90 -= M_2PI;
            }

            //For each point along perpendicular course
            for (j=0; j<5; j++)
            {

                err |= direct(testBase,crsTE90,d90[j],&test[j],EPS);
                latTest = test[j].latitude / DEG2RAD;
                lonTest = test[j].longitude / DEG2RAD;

                if( getMaskedError(err, getMaskAll()) )
                {
                        printf("Error occurred in pre-ptIsOnCrs err=0x%lx i=%d j=%d\n", err,i,j);
                        setupFailureCount++;
                        errorCount++;
                        testCaseCount++;
                        continue;
                }

                err |= invCrs(start,end,&crsTest,&tempdbl,EPS);
                if (d[i] < 0.0)
                    crsTest = modcrs(crsTest + M_PI);
                onLine = ptIsOnCrs(start,crsTest,test[j],&crsTest1,&dist1Test,&err,TOL,EPS);

                    if( getMaskedError(err, getMaskAll()) )
                    {
                        printf("Error occurred in ptIsOnCrs err=0x%lx i=%d j=%d\n", err,i,j);
                        failedCount++;
                        errorCount++;
                        testCaseCount++;
                        continue;
                    }

                    if ( onLine == onLineExp[i][j] )
                    {
                    	if (PRINT_PASSED_CASES) printf("%s Start: %14.8f %14.8f  az1=%14.8f d1=%14.8f\n",
                               testname, latStart, lonStart, az1, d1);
                    	if (PRINT_PASSED_CASES) printf("%s End: %14.8f %14.8f\n",testname, latEnd, lonEnd);

                    	if (PRINT_PASSED_CASES) printf("%s passed  test: %14.8f %14.8f  dist=%14.8f d90=%14.8f  onLine=%d\n",
                               testname,latTest,lonTest,d[i],d90[j],onLine);
                        passedCount++;
                    }
                    else
                    {
                        printf("%s Start: %14.8f %14.8f  az1=%14.8f d1=%14.8f\n",
                               testname, latStart, lonStart, az1, d1);
                        printf("%s End: %14.8f %14.8f\n",testname, latEnd, lonEnd);

                        printf("%s failed  test: %14.8f %14.8f  dist=%14.8f d90=%14.8f  onLine=%d\n",
                               testname,latTest,lonTest,d[i],d90[j],onLine);
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

    printf("\nFinish testPtIsOnCrs_Set1\n\n\n");

    return set;
}

/*
 * NAME: testPtIsOnCrs_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the ptIsOnCrs function.
 *
 * 		This function runs all the test cases from all the test sets for testing
 *      the above function and returns the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testPtIsOnCrs_AllSets(TestSuite) - A test suite with the folling metrics:
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
TestSuite testPtIsOnCrs_AllSets()
{
    TestSuite suite;
    TestSet set1;

    printf("\nStart testPtIsOnCrs_AllSets\n");

    suite = newTestSuite("testPtIsOnCrs_AllSets");

    set1 = testPtIsOnCrs_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("\nFinish testPtIsOnCrs_AllSets\n\n\n");

    return suite;

}
} //namespace
