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
 * NAME: testPtsAreSame_Set1
 *
 * DESCRIPTION:
 * 		This function is used to test the ptsAreSame function.
 *
 *		Test data is created by generating a random point then using the direct function
 *		to project a second point a small distance away.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testPtsAreSame_Set1(TestSet) - A test set with the following metrics:
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
TestSet testPtsAreSame_Set1()
{
	double DEG2RAD = M_PI / 180.0;
	double RAD2DEG = 180.0 / M_PI;

    LLPoint point1, point2;

    int expected, areSame; //0=false,1=true

    int passedCount=0, failedCount=0;
    int errorCount = 0;
    int unverifiedCount = 0;
    int setupFailureCount = 0;
    int testCaseCount = 0;
    ErrorSet err = 0;
    TestSet set;

    double latRanges[5][2] = {{ 0, 30 }, {30, 45}, {45, 45}, {45, 60}, {60, 90}};
    double lonRanges[5][2] = {{0, 60}, {60, 90}, {90, 90}, {90, 120}, {120, 180}};
    double azimuths[8] = {0.0, 45.0, 90.0, 135.0, 180.0, 225.0, 270.0, 315.0};
    double latsign[] = {-1, 1};
    double lonsign[] = {-1, 1};
    double distances[3] = {0.5*TOL, 0.99*TOL, 2.0*TOL}; //Due to round-off error we cannot use 1.0*TOL. Failures first appear at .999*TOL
    double randLat, randLon;
    int i,j,k,l,m,n;
    long newSeed = 10080111;

    srand(newSeed);

    printf("\n\nStart testPtsAreSame_Set1\n");

    set = newTestSet("testPtsAreSame_Set1");

    for(i = 0; i < 5; i++){
    	for(j = 0; j < 5; j++){
    		randLat = (latRanges[i][0] + (latRanges[i][1] - latRanges[i][0]) * (double) (rand())  / (double) RAND_MAX);
    		randLon = (lonRanges[j][0] + (lonRanges[j][1] - lonRanges[j][0]) * (double) (rand())  / (double) RAND_MAX);
    		for(k = 0; k < 2; k++){
    			for(l = 0; l < 2; l++){
    	    		point1.latitude = latsign[k]*randLat*DEG2RAD;
    	    		point1.longitude = lonsign[l]*randLon*DEG2RAD;
    	    		for(m = 0; m < 8; m++){
    	    			for(n = 0; n < 3; n++){
    	    				testCaseCount++;

    	    				err |= direct(point1,azimuths[m]*DEG2RAD,distances[n],&point2,EPS);
//	    					printf("point1: %3.20lf, %3.20lf \n", point1.latitude, point1.longitude);
//	    					printf("point2: %3.20lf, %3.20lf \n", point2.latitude, point2.longitude);
    	    				if(err){
    	    					errorCount++;
    	    					err = 0;
    	    				}

    	    				expected = (distances[n] < TOL) ? 1 : 0;
    	    				areSame = ptsAreSame(point1,point2,TOL);
    	    				if(areSame == expected){
    	    					passedCount++;
    	    				} else {
    	    					failedCount++;
    	    					printf("Failed: %d", testCaseCount);
    	    					printf("point1: %3.20lf, %3.20lf \n", point1.latitude*RAD2DEG, point1.longitude*RAD2DEG);
    	    					printf("point2: %3.20lf, %3.20lf \n", point2.latitude*RAD2DEG, point2.longitude*RAD2DEG);
    	    					printf("expected: %d, areSame: %d", expected, areSame);
    	    				}

    	    			}
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

    printf("\nFinish testPtsAreSame_Set1\n\n\n");

    return set;
}

/*
 * NAME: testPtsAreSame_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the ptsAreSame function.
 *
 * 		This function runs all the test cases from all the sets that test ptsAreSame
 * 		and returns the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testPtsAreSame_AllSets(TestSuite) - A test suite with the following metrics:
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
TestSuite testPtsAreSame_AllSets()
{

    TestSet set1;
    TestSuite suite;

    printf("\nStart testPtsAreSame_AllSets\n");

    suite = newTestSuite("testPtsAreSame_AllSets");

    set1 = testPtsAreSame_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);


    printf("Finish testPtsAreSame_AllSets\n\n\n");

    return suite;
}

/*
 * NAME: testsphereInverse
 *
 * DESCRIPTION:
 * 		This function is used to test the sphereInverse function.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testsphereInverse(TestSet) - A test set with the following metrics:
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
TestSet testSphereInverse_Set1()
{
    double DEG2RAD = M_PI / 180.0;
//    double RAD2DEG = 180.0 / M_PI;
    double NMTOL = 0.03 / 100.0 / 1852.0; //0.03 cm or ~1.62e-7 nm
//    double LLDEGTOL = 5.0e-8; //5e-8 degrees or 0.00018 sec
    double AZDEGTOL = 5.0e-7; //5e-7 deg or 8.72665e-9 rad or 0.162 cm at 100 nm
    double NMSHORTDIST = SPHERICAL_AZIMUTH_CUTOFF_DISTANCE; //5e-4 nm or 0.926 m
    double NMLARGEDIST = 10000.0; //10000 nm
    //LLPoint org, LLPoint dest,double* crs, double* dist, double eps);
    FILE *inputfile;
    int mm;
    char inputFileName[] = FILEROOT "/test/testLLPoint/data/testSphereInverse/sphereTestDataWithExpectedResults.txt";
    char testname[80];
    double lat1, lon1, lat2, lon2;
    double crs12Exp, crs21Exp, dist12Exp;
    LLPoint p1, p2;
    double crs12, crs21, dist12;
//    double inputROC;
    int passedCount = 0, failedCount = 0;
    int testCaseCount = 0, setupFailureCount = 0, errorCount = 0, unverifiedCount = 0;

    TestSet set;

    printf("Start testSphereInverse\n");

    set = newTestSet("testSphereInverse_Set1");

    inputfile = fopen(inputFileName, "r");

    if (inputfile == NULL)
    {
        printf("No input file %s found", inputFileName);
        set.testingError = 1;
        return set;
    }

    while (1)
    {
        //Input data: string, deg, deg, deg, deg, deg, deg, nm
        mm = fscanf(inputfile, "%s %lf %lf %lf %lf %lf %lf %lf", testname,
                &lat1, &lon1, &lat2, &lon2, &dist12Exp, &crs12Exp, &crs21Exp);
        if (mm == EOF)
        {
            fclose(inputfile);
            break;
        }

        p1.latitude = (lat1 * DEG2RAD);
        p1.longitude = lon1 * DEG2RAD;
        p2.latitude = (lat2 * DEG2RAD);
        p2.longitude = lon2 * DEG2RAD;

        sphereInverse(p1, p2, &crs12, &dist12, NULL, EPS);
        sphereInverse(p2, p1, &crs21, &dist12, NULL, EPS);

        if (((fabs(crs12Exp - crs12 / DEG2RAD) < AZDEGTOL) || (fabs(360.0
                - fabs(crs12Exp - crs12 / DEG2RAD)) < AZDEGTOL)) && ((fabs(
                crs21Exp - crs21 / DEG2RAD) < AZDEGTOL) || (fabs(360.0 - fabs(
                crs21Exp - crs21 / DEG2RAD)) < AZDEGTOL)) && (fabs(dist12Exp
                - dist12) < NMTOL))
        {
        	if (PRINT_PASSED_CASES) printf("%s passed\n", testname);
            passedCount++;
            testCaseCount++;
        }

        else if ((dist12Exp < NMSHORTDIST) || (dist12Exp > NMLARGEDIST))
        {
            if (fabs(dist12Exp - dist12) < NMTOL)
            {
            	if (PRINT_PASSED_CASES) printf("%s passed\n", testname);
            	if (PRINT_PASSED_CASES) printf("%s az1diff= %14.8e az2diff= %14.8e distdiff= %14.8e\n",
                        testname, fabs(crs12Exp - crs12 / DEG2RAD), fabs(
                                crs21Exp - crs21 / DEG2RAD), fabs(dist12Exp
                                - dist12));
                passedCount++;
                testCaseCount++;
            }
            else
            {
                printf(
                        "%s p1: %14.8f %14.8f  p2: %14.8f %14.8f  crs: %12.6f %12.6f  dist: %14.8f\n",
                        testname, lat1, lon1, lat2, lon2, crs12Exp, crs21Exp,
                        dist12Exp);
                printf("%s failed %14.8f %14.8f %14.8f %14.8f %14.8f %14.8f\n",
                        testname, crs12Exp, crs12 / DEG2RAD, crs21Exp, crs21
                                / DEG2RAD, dist12Exp, dist12);
                printf("%s az1diff= %14.8e az2diff= %14.8e distdiff= %14.8e\n",
                        testname, fabs(crs12Exp - crs12 / DEG2RAD), fabs(
                                crs21Exp - crs21 / DEG2RAD), fabs(dist12Exp
                                - dist12));
                failedCount++;
                testCaseCount++;
            }
        }

        else
        {
            printf(
                    "%s p1: %14.8f %14.8f  p2: %14.8f %14.8f  crs: %12.6f %12.6f  dist: %14.8f\n",
                    testname, lat1, lon1, lat2, lon2, crs12Exp, crs21Exp,
                    dist12Exp);
            printf("%s failed %14.8f %14.8f %14.8f %14.8f %14.8f %14.8f\n",
                    testname, crs12Exp, crs12 / DEG2RAD, crs21Exp, crs21
                            / DEG2RAD, dist12Exp, dist12);
            printf("%s az1diff= %14.8e az2diff= %14.8e distdiff= %14.8e\n",
                    testname, fabs(crs12Exp - crs12 / DEG2RAD), fabs(crs21Exp
                            - crs21 / DEG2RAD), fabs(dist12Exp - dist12));
            failedCount++;
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

    printf("Finish testSphereInverse_Set1\n");

    return set;

}

/*
 * NAME: testSphereInverse_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the testsphereInverse function.
 *
 * 		This function runs all the test cases from all the sets that test testsphereInverse
 * 		and returns the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testSphereInverse_AllSets(TestSuite) - A test set with the following metrics:
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
TestSuite testSphereInverse_AllSets()
{
    TestSuite suite;
    TestSet set1;

    printf("\nStart testSphereInverse_AllSets\n");

    suite = newTestSuite("testSphereInverse_AllSets");

    set1 = testSphereInverse_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("Finish testSphereInverse_AllSets\n\n\n");

    return suite;
}

/*
 * NAME: testsphereInvDist_Set1
 *
 * DESCRIPTION:
 * 		This function is used to test the sphereInvDist function.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testsphereInvDist(TestSet) - A test set with the following metrics:
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
TestSet testsphereInvDist_Set1()
{
    double DEG2RAD = M_PI / 180.0;
//    double RAD2DEG = 180.0 / M_PI;
    double NMTOL = 0.03 / 100.0 / 1852.0; //0.03 cm or ~1.62e-7 nm
//    double LLDEGTOL = 5.0e-8; //5e-8 degrees or 0.00018 sec
//    double AZDEGTOL = 5.0e-7; //5e-7 deg or 8.72665e-9 rad or 0.162 cm at 100 nm
//    double NMSHORTDIST = SPHERICAL_AZIMUTH_CUTOFF_DISTANCE; //5e-4 nm or 0.926 m
//    double NMLARGEDIST = 10000.0; //10000 nm
    //LLPoint org, LLPoint dest);
    FILE *inputfile;
    int mm;
    char inputFileName[] = FILEROOT "/test/testLLPoint/data/testSphereInvDist/sphereTestDataWithExpectedResults.txt";
    char testname[80];
    double lat1, lon1, lat2, lon2;
    double crs12Exp, crs21Exp, dist12Exp;
    LLPoint p1, p2;
    double dist12;
//    double inputROC;
    int passedCount = 0, failedCount = 0;
    int testCaseCount = 0, setupFailureCount = 0, errorCount = 0, unverifiedCount = 0;

    TestSet set;

    printf("Start testsphereInvDist\n");

    set = newTestSet("testsphereInvDist_Set1");

    inputfile = fopen(inputFileName, "r");

    if (inputfile == NULL)
    {
        printf("No input file %s found", inputFileName);
        set.testingError = 1;
        return set;
    }

    while (1)
    {
        //Input data: string, deg, deg, deg, deg, deg, deg, nm
        mm = fscanf(inputfile, "%s %lf %lf %lf %lf %lf %lf %lf", testname,
                &lat1, &lon1, &lat2, &lon2, &dist12Exp, &crs12Exp, &crs21Exp);
        if (mm == EOF)
        {
            fclose(inputfile);
            break;
        }

        p1.latitude = (lat1 * DEG2RAD);
        p1.longitude = lon1 * DEG2RAD;
        p2.latitude = (lat2 * DEG2RAD);
        p2.longitude = lon2 * DEG2RAD;

        dist12 = sphereInvDist(p1, p2, NULL);

        if (fabs(dist12Exp - dist12) < NMTOL)
        {
        	if (PRINT_PASSED_CASES) printf("%s passed\n", testname);
            passedCount++;
            testCaseCount++;
        }
        else
        {
            printf(
                    "%s p1: %14.8f %14.8f  p2: %14.8f %14.8f  crs: %12.6f %12.6f  dist: %14.8f\n",
                    testname, lat1, lon1, lat2, lon2, crs12Exp, crs21Exp,
                    dist12Exp);
            printf("%s failed %14.8f %14.8f distdiff= %14.8e\n", testname,
                    dist12Exp, dist12, fabs(dist12Exp - dist12));
            failedCount++;
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

    printf("Finish testsphereInvDist_Set1\n");

    return set;

}

/*
 * NAME: testsphereInvDist_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the testsphereInvDist function.
 *
 * 		This function runs all the test cases from all the sets that test testsphereInvDist
 * 		and returns the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testsphereInvDist_AllSets(TestSuite) - A test set with the following metrics:
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
TestSuite testsphereInvDist_AllSets()
{
    TestSuite suite;
    TestSet set1;

    printf("\nStart testsphereInvDist_AllSets\n");

    suite = newTestSuite("testsphereInvDist_AllSets");

    set1 = testsphereInvDist_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("Finish testsphereInvDist_AllSets\n\n\n");

    return suite;
}

/*
 * NAME: testsphereInvCrs_Set1
 *
 * DESCRIPTION:
 * 		This function is used to test the sphereInvCrs function.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testsphereInvCrs_Set1(TestSet) - A test set with the following metrics:
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
TestSet testsphereInvCrs_Set1()
{
    double DEG2RAD = M_PI / 180.0;
//    double RAD2DEG = 180.0 / M_PI;
//    double NMTOL = 0.03 / 100.0 / 1852.0; //0.03 cm or ~1.62e-7 nm
//    double LLDEGTOL = 5.0e-8; //5e-8 degrees or 0.00018 sec
    double AZDEGTOL = 5.0e-7; //5e-7 deg or 8.72665e-9 rad or 0.162 cm at 100 nm
//    double NMSHORTDIST = SPHERICAL_AZIMUTH_CUTOFF_DISTANCE; //5e-4 nm or 0.926 m
//    double NMLARGEDIST = 10000.0; //10000 nm
    //LLPoint org, LLPoint dest, double eps);
    FILE *inputfile;
    int mm;
    char inputFileName[] = FILEROOT "/test//testLLPoint/data/testSphereInvCrs/sphereTestDataWithExpectedResults.txt";
    char testname[80];
    double lat1, lon1, lat2, lon2;
    double crs12Exp, crs21Exp, dist12Exp;
    LLPoint p1, p2;
    double crs12, crs21;
    int passedCount = 0, failedCount = 0;
    int testCaseCount = 0, setupFailureCount = 0, errorCount = 0, unverifiedCount = 0;

    TestSet set;

    printf("Start testsphereInvCrs\n");

    set = newTestSet("testsphereInvCrs_Set1");

    inputfile = fopen(inputFileName, "r");

    if (inputfile == NULL)
    {
        printf("No input file %s found", inputFileName);
        set.testingError = 1;
        return set;
    }

    while (1)
    {
        //Input data: string, deg, deg, deg, deg, deg, deg, nm
        mm = fscanf(inputfile, "%s %lf %lf %lf %lf %lf %lf %lf", testname,
                &lat1, &lon1, &lat2, &lon2, &dist12Exp, &crs12Exp, &crs21Exp);
        if (mm == EOF)
        {
            fclose(inputfile);
            break;
        }

        p1.latitude = (lat1 * DEG2RAD);
        p1.longitude = lon1 * DEG2RAD;
        p2.latitude = (lat2 * DEG2RAD);
        p2.longitude = lon2 * DEG2RAD;

        crs12 = sphereInvCrs(p1, p2, EPS);
        crs21 = sphereInvCrs(p2, p1, EPS);

        if (((fabs(crs12Exp - crs12 / DEG2RAD) < AZDEGTOL) || (fabs(360.0
                - fabs(crs12Exp - crs12 / DEG2RAD)) < AZDEGTOL)) && ((fabs(
                crs21Exp - crs21 / DEG2RAD) < AZDEGTOL) || (fabs(360.0 - fabs(
                crs21Exp - crs21 / DEG2RAD)) < AZDEGTOL)))
        {
        	if (PRINT_PASSED_CASES) printf("%s passed\n", testname);
            passedCount++;
            testCaseCount++;
        }
        else
        {
            printf(
                    "%s p1: %14.8f %14.8f  p2: %14.8f %14.8f  crs: %12.6f %12.6f  dist: %14.8f\n",
                    testname, lat1, lon1, lat2, lon2, crs12Exp, crs21Exp,
                    dist12Exp);
            printf("%s failed %14.8f %14.8f %14.8f %14.8f\n", testname,
                    crs12Exp, crs12 / DEG2RAD, crs21Exp, crs21 / DEG2RAD);
            printf("%s az1diff= %14.8e az2diff= %14.8e\n", testname, fabs(
                    crs12Exp - crs12 / DEG2RAD), fabs(crs21Exp - crs21
                    / DEG2RAD));
            failedCount++;
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

    return set;
}

/*
 * NAME: testsphereInvCrs_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the testsphereInvCrs function.
 *
 * 		This function runs all the test cases from all the sets that test testsphereInvCrs
 * 		and returns the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testsphereInvCrs_AllSets(TestSuite) - A test set with the following metrics:
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
TestSuite testsphereInvCrs_AllSets()
{
    TestSuite suite;
    TestSet set1;

    printf("\nStart testsphereInvCrs_AllSets\n");

    suite = newTestSuite("testsphereInvCrs_AllSets");

    set1 = testsphereInvCrs_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    return suite;
}

/*
 * NAME: testInverse_Set1
 *
 * DESCRIPTION:
 * 		This function is used to test the inverse function
 *
 * The test data comes from the GeographicLib project by Charles Karney. Which uses an extension of Bessel's and Helmert's method of solving
 * the direct problem. Data and code can be obtained from http://geographiclib.sourceforge.net/. A description of the data follows.
 *
 * "Each line of the test set gives 9 space delimited numbers
 *
 * latitude for point 1, lat1 (degrees, exact)
 * longitude for point 1, lon1 (degrees, always 0)
 * azimuth for point 1, azi1 (clockwise from north in degrees, exact)
 * latitude for point 2, lat2 (degrees, accurate to 10^-18 deg)
 * longitude for point 2, lon2 (degrees, accurate to 10^-18 deg)
 * azimuth for point 2, azi2 (degrees, accurate to 10^-18 deg)
 * geodesic distance from point 1 to point 2, s12 (meters, exact)
 * arc distance on the auxiliary sphere, a12 (degrees, accurate to 10^-18 deg)
 * reduced length of the geodesic, m12 (meters, accurate to 0.1 pm)
 * the area under the geodesic, S12 (meters^2, accurate to 1 mm^2)
 *
 * These are computed using as direct geodesic calculations with the given lat1, lon1, azi1, and s12. The distance
 * s12 always corresponds to an arc length a12 <= 180 deg, so the given geodesics give the shortest paths from
 * point 1 to point 2. For simplicity and without loss of generality, lat1 is chosen in [0 deg, 90 deg], lon1 is
 * taken to be zero, azi1 is chosen in [0 deg, 180 deg]. Furthermore, lat1 and azi1 are taken to be multiples
 * of 10^-12 deg and s12 is a multiple of 0.1 um in [0 m, 20003931.4586254 m]. This results lon2 in [0 deg, 180 deg] and
 * azi2 in [0 deg, 180 deg].
 *
 * The direct calculation uses an expansion of the geodesic equations accurate to f^20 (approximately 1 part in 10^50) and is computed with with Maxima's bfloats and fpprec set to 100 (so the errors in the data are probably 1/2 of the values quoted above).
 *
 * The contents of [GeodTest.dat] are as follows [GeodTest-short.dat consists of 1/50 of these tests]:
 *
 * 100000 entries randomly distributed
 * 50000 entries which are nearly antipodal
 * 50000 entries with short distances
 * 50000 entries with one end near a pole
 * 50000 entries with both ends near opposite poles
 * 50000 entries which are nearly meridional
 * 50000 entries which are nearly equatorial
 * 50000 entries running between vertices (azi1 = azi2 = 90 deg)
 * 50000 entries ending close to vertices
 *
 * (a total of 500000 entries). The values for s12 for the geodesics running between vertices are truncated
 * to a multiple of 0.1 pm and this is used to determine point 2"
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testInverse_Set1(TestSet) - A test set with the following metrics:
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
TestSet testInverse_Set1()
{
    double DEG2RAD = M_PI / 180.0;
    double RAD2DEG = 180.0 / M_PI;
    double NMTOL = 0.03 / 100.0 / 1852.0; //0.03 cm or ~1.62e-7 nm
    double AZDEGTOL = 5.0e-7; //5e-7 deg or 8.72665e-9 rad or 0.162 cm at 100 nm
    FILE *inputfile;
    int mm;
    char inputFileName[] = FILEROOT "/geolib/src/test/resources/testLLPoint/data/GeodTest.dat";
    /* For a longer test set (500,000 cases) use GeodTest.dat. For a shorter test set (10,000 cases) use GeodTest-short.dat. */
    char testname[80];
    double lat1, lon1, lat2, lon2;
    double crs12Exp, crs21Exp, dist12Exp;
    double arcDistance12, reducedLength12, area12;
    LLPoint p1, p2;
    double crs12, crs21, dist12;
    int passedCount = 0, failedCount = 0, unverifiedCount = 0;
    int setupFailureCount = 0, errorCount = 0, testCaseCount = 0;
    ErrorSet err = 0;
    int testNum = 0;
    int failed = 0;
    double crs12Diff, crs21Diff, distDiff;

    TestSet set;

    printf("Start testInverse_Set1\n");

    set = newTestSet("testInverse_Set1");

    inputfile = fopen(inputFileName, "r");

    if (inputfile == NULL)
    {
        printf("No input file %s found", inputFileName);
        set.testingError = 1;
        return set;
    }

    int skipCount = 0;

    while (1)
    {
        //Input data: deg, deg, deg, deg, deg, deg, m, deg, m, m^2
        mm = fscanf(inputfile, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &lat1,
                &lon1, &crs12Exp, &lat2, &lon2, &crs21Exp, &dist12Exp, &arcDistance12, &reducedLength12, &area12);
        if (mm == EOF)
        {
            fclose(inputfile);
            break;
        }
        if(dist12Exp < 1 || dist12Exp > 10740.0*1852.0) {
           skipCount++;
           continue; //Vincenty's algorithm will fail in the nearly anti-podal case
        }
        testNum++;

        p1.latitude = lat1 * DEG2RAD;
        p1.longitude = lon1 * DEG2RAD;
        p2.latitude = lat2 * DEG2RAD;
        p2.longitude = lon2 * DEG2RAD;
        sprintf(testname,"Test%i",testNum);

        err = inverse(p1, p2, &crs12, &crs21, &dist12, EPS);

        if (getMaskedError(err, getMaskAll()))
        {
            printf("Error occurred in inverse err=0x%lx\n", err);
            failedCount++;
            errorCount++;
            testCaseCount++;
            continue;
        }

        failed = 0;

        crs12Diff = fabs(crs12 - crs12Exp*DEG2RAD);
        crs21Diff = fabs(crs21 - fmod(crs21Exp + 180.0,360.0)*DEG2RAD);
        distDiff = fabs(dist12 - dist12Exp/1852.0);

        if ((crs12Diff*RAD2DEG) > AZDEGTOL)
        {
        	failed = 1;
        	printf("\n%s failed: crs12=%.20lf crs12Exp=%.20lf crs12diff=%le", testname, crs12*RAD2DEG, crs12Exp, crs12Diff*RAD2DEG);
        }
        if ((crs21Diff*RAD2DEG) > AZDEGTOL)
        {
        	printf("\n%s failed: crs21=%.20lf crs21Exp=%.20lf crs21diff=%le", testname, crs21*RAD2DEG, fmod(crs21Exp + 180.0,360.0), crs21Diff*RAD2DEG);
        	failed = 1;
        }
        if (distDiff > NMTOL)
        {
        	printf("\n%s failed: dist12=%.20lf dist12Exp=%.20lf distDiff=%le", testname, dist12, dist12Exp, distDiff);
        	failed = 1;
        }

        if (failed)
        {
			printf("\n%s failed p1: %.20lf %.20lf  p2: %.20lf %.20lf\n",
					testname, lat1, lon1, lat2, lon2);
			printf("DistExpected: %.20lf Dist:%.20lf\n", dist12Exp/1852.0, dist12);
			printf("crs12: %.20lf crs12Exp%.20lf\n",crs12/DEG2RAD,crs12Exp);
			printf("crs21: %.20lf crs21Exp%.20lf\n",crs21/DEG2RAD,fmod(crs21Exp + 180.0,360.0));
			failedCount++;
        } else {
        	passedCount++;
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

    printf("Finish testInverse_Set1, skipped %d cases.\n", skipCount);

    return set;

}


/*
 * NAME: testInverse_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the testInverse function.
 *
 * 		This function runs all the test cases from all the sets that test testInverse
 * 		and returns the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testInverse_AllSets(TestSuite) - A test set with the following metrics:
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
TestSuite testInverse_AllSets()
{
    TestSuite suite;
    TestSet set1;

    printf("\nStart testInverse_AllSets\n");

    suite = newTestSuite("testInverse_AllSets");

    set1 = testInverse_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("Finish testInverse_AllSets\n\n\n");

    return suite;
}

/*
 * NAME: testInvCrs_Set1
 *
 * DESCRIPTION:
 * 		This function is used to test the invCrs function
 *
 * The test data comes from the GeographicLib project by Charles Karney. Which uses an extension of Bessel's and Helmert's method of solving
 * the direct problem. Data and code can be obtained from http://geographiclib.sourceforge.net/. A description of the data follows.
 *
 * "Each line of the test set gives 9 space delimited numbers
 *
 * latitude for point 1, lat1 (degrees, exact)
 * longitude for point 1, lon1 (degrees, always 0)
 * azimuth for point 1, azi1 (clockwise from north in degrees, exact)
 * latitude for point 2, lat2 (degrees, accurate to 10^-18 deg)
 * longitude for point 2, lon2 (degrees, accurate to 10^-18 deg)
 * azimuth for point 2, azi2 (degrees, accurate to 10^-18 deg)
 * geodesic distance from point 1 to point 2, s12 (meters, exact)
 * arc distance on the auxiliary sphere, a12 (degrees, accurate to 10^-18 deg)
 * reduced length of the geodesic, m12 (meters, accurate to 0.1 pm)
 * the area under the geodesic, S12 (meters^2, accurate to 1 mm^2)
 *
 * These are computed using as direct geodesic calculations with the given lat1, lon1, azi1, and s12. The distance
 * s12 always corresponds to an arc length a12 <= 180 deg, so the given geodesics give the shortest paths from
 * point 1 to point 2. For simplicity and without loss of generality, lat1 is chosen in [0 deg, 90 deg], lon1 is
 * taken to be zero, azi1 is chosen in [0 deg, 180 deg]. Furthermore, lat1 and azi1 are taken to be multiples
 * of 10^-12 deg and s12 is a multiple of 0.1 um in [0 m, 20003931.4586254 m]. This results lon2 in [0 deg, 180 deg] and
 * azi2 in [0 deg, 180 deg].
 *
 * The direct calculation uses an expansion of the geodesic equations accurate to f^20 (approximately 1 part in 10^50) and is computed with with Maxima's bfloats and fpprec set to 100 (so the errors in the data are probably 1/2 of the values quoted above).
 *
 * The contents of [GeodTest.dat] are as follows [GeodTest-short.dat consists of 1/50 of these tests]:
 *
 * 100000 entries randomly distributed
 * 50000 entries which are nearly antipodal
 * 50000 entries with short distances
 * 50000 entries with one end near a pole
 * 50000 entries with both ends near opposite poles
 * 50000 entries which are nearly meridional
 * 50000 entries which are nearly equatorial
 * 50000 entries running between vertices (azi1 = azi2 = 90 deg)
 * 50000 entries ending close to vertices
 *
 * (a total of 500000 entries). The values for s12 for the geodesics running between vertices are truncated
 * to a multiple of 0.1 pm and this is used to determine point 2"
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testInvCrs_Set1(TestSet) - A test set with the following metrics:
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
TestSet testInvCrs_Set1()
{
    double DEG2RAD = M_PI / 180.0;
    double RAD2DEG = 180.0 / M_PI;
    double AZDEGTOL = 5.0e-7; //5e-7 deg or 8.72665e-9 rad or 0.162 cm at 100 nm
    FILE *inputfile;
    int mm;
    char inputFileName[] = FILEROOT "/geolib/src/test/resources/testLLPoint/data/GeodTest-short.dat";
    /* For a longer test set (500,000 cases) use GeodTest.dat. For a shorter test set (10,000 cases) use GeodTest-short.dat. */
    char testname[80];
    double lat1, lon1, lat2, lon2;
    double crs12Exp, crs21Exp, dist12Exp;
    double arcDistance12, reducedLength12, area12;
    LLPoint p1, p2;
    double crs12, crs21;
    int passedCount = 0, failedCount = 0, unverifiedCount = 0;
    int setupFailureCount = 0, errorCount = 0, testCaseCount = 0;
    ErrorSet err = 0;
    int testNum = 0;
    int failed = 0;
    double crs12Diff, crs21Diff;

    TestSet set;

    printf("Start testInvCrs_Set1\n");

    set = newTestSet("testInvCrs_Set1");

    inputfile = fopen(inputFileName, "r");

    if (inputfile == NULL)
    {
        printf("No input file %s found", inputFileName);
        set.testingError = 1;
        return set;
    }

    while (1)
    {
        //Input data: deg, deg, deg, deg, deg, deg, m, deg, m, m^2
        mm = fscanf(inputfile, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &lat1,
                &lon1, &crs12Exp, &lat2, &lon2, &crs21Exp, &dist12Exp, &arcDistance12, &reducedLength12, &area12);
        if (mm == EOF)
        {
            fclose(inputfile);
            break;
        }
        if(dist12Exp > 10740.0*1852.0) continue; //Vincenty's algorithm will fail in the nearly anti-podal case
        testNum++;

        p1.latitude = lat1 * DEG2RAD;
        p1.longitude = lon1 * DEG2RAD;
        p2.latitude = lat2 * DEG2RAD;
        p2.longitude = lon2 * DEG2RAD;
        sprintf(testname,"Test%i",testNum);

        err = invCrs(p1, p2, &crs12, &crs21,EPS);

        if (getMaskedError(err, getMaskAll()))
        {
            printf("Error occurred in invCrs err=0x%lx\n", err);
            failedCount++;
            errorCount++;
            testCaseCount++;
            continue;
        }

        failed = 0;

        crs12Diff = fabs(crs12 - crs12Exp*DEG2RAD);
        crs21Diff = fabs(crs21 - fmod(crs21Exp + 180.0,360.0)*DEG2RAD);

        if ((crs12Diff*RAD2DEG) > AZDEGTOL)
        {
        	failed = 1;
        	printf("\n%s failed: crs12=%.20lf crs12Exp=%.20lf crs12diff=%le", testname, crs12*RAD2DEG, crs12Exp, crs12Diff*RAD2DEG);
        }
        if ((crs21Diff*RAD2DEG) > AZDEGTOL)
        {
        	printf("\n%s failed: crs21=%.20lf crs21Exp=%.20lf crs21diff=%le", testname, crs21*RAD2DEG, crs21Exp, crs21Diff*RAD2DEG);
        	failed = 1;
        }

        if (failed)
        {
			printf("\n%s failed p1: %.20lf %.20lf  p2: %.20lf %.20lf\n",
					testname, lat1, lon1, lat2, lon2);
			printf("crs12: %.20lf crs12Exp%.20lf\n",crs12/DEG2RAD,crs12Exp);
			printf("crs12: %.20lf crs12Exp%.20lf\n",crs21/DEG2RAD,crs21Exp);
			failedCount++;
        } else {
        	passedCount++;
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

    printf("Finish testInvCrs_Set1\n");

    return set;

}
/*
 * NAME: testInvCrs_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the testInvCrs function.
 *
 * 		This function runs all the test cases from all the sets that test testInvCrs
 * 		and returns the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testInvCrs_AllSets(TestSuite) - A test set with the following metrics:
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
TestSuite testInvCrs_AllSets()
{
    TestSuite suite;
    TestSet set1;

    printf("\nStart testInvCrs_AllSets\n");

    suite = newTestSuite("testInvCrs_AllSets");

    set1 = testInvCrs_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("Finish testInvCrs_AllSets\n\n\n");

    return suite;
}

/*
 * NAME: testInvDist_Set1
 *
 * DESCRIPTION:
 * 		This function is used to test the inverse function
 *
 * The test data comes from the GeographicLib project by Charles Karney. Which uses an extension of Bessel's and Helmert's method of solving
 * the direct problem. Data and code can be obtained from http://geographiclib.sourceforge.net/. A description of the data follows.
 *
 * "Each line of the test set gives 9 space delimited numbers
 *
 * latitude for point 1, lat1 (degrees, exact)
 * longitude for point 1, lon1 (degrees, always 0)
 * azimuth for point 1, azi1 (clockwise from north in degrees, exact)
 * latitude for point 2, lat2 (degrees, accurate to 10^-18 deg)
 * longitude for point 2, lon2 (degrees, accurate to 10^-18 deg)
 * azimuth for point 2, azi2 (degrees, accurate to 10^-18 deg)
 * geodesic distance from point 1 to point 2, s12 (meters, exact)
 * arc distance on the auxiliary sphere, a12 (degrees, accurate to 10^-18 deg)
 * reduced length of the geodesic, m12 (meters, accurate to 0.1 pm)
 * the area under the geodesic, S12 (meters^2, accurate to 1 mm^2)
 *
 * These are computed using as direct geodesic calculations with the given lat1, lon1, azi1, and s12. The distance
 * s12 always corresponds to an arc length a12 <= 180 deg, so the given geodesics give the shortest paths from
 * point 1 to point 2. For simplicity and without loss of generality, lat1 is chosen in [0 deg, 90 deg], lon1 is
 * taken to be zero, azi1 is chosen in [0 deg, 180 deg]. Furthermore, lat1 and azi1 are taken to be multiples
 * of 10^-12 deg and s12 is a multiple of 0.1 um in [0 m, 20003931.4586254 m]. This results lon2 in [0 deg, 180 deg] and
 * azi2 in [0 deg, 180 deg].
 *
 * The direct calculation uses an expansion of the geodesic equations accurate to f^20 (approximately 1 part in 10^50) and is computed with with Maxima's bfloats and fpprec set to 100 (so the errors in the data are probably 1/2 of the values quoted above).
 *
 * The contents of [GeodTest.dat] are as follows [GeodTest-short.dat consists of 1/50 of these tests]:
 *
 * 100000 entries randomly distributed
 * 50000 entries which are nearly antipodal
 * 50000 entries with short distances
 * 50000 entries with one end near a pole
 * 50000 entries with both ends near opposite poles
 * 50000 entries which are nearly meridional
 * 50000 entries which are nearly equatorial
 * 50000 entries running between vertices (azi1 = azi2 = 90 deg)
 * 50000 entries ending close to vertices
 *
 * (a total of 500000 entries). The values for s12 for the geodesics running between vertices are truncated
 * to a multiple of 0.1 pm and this is used to determine point 2"
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testInvDist_Set1(TestSet) - A test set with the following metrics:
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
TestSet testInvDist_Set1()
{
    double DEG2RAD = M_PI / 180.0;
    double NMTOL = 0.03 / 100.0 / 1852.0; //0.03 cm or ~1.62e-7 nm
    FILE *inputfile;
    int mm;
    char inputFileName[] = FILEROOT "/geolib/src/test/resources/testLLPoint/data/GeodTest-short.dat";
    /* For a longer test set (500,000 cases) use GeodTest.dat. For a shorter test set (10,000 cases) use GeodTest-short.dat. */
    char testname[80];
    double lat1, lon1, lat2, lon2;
    double crs12Exp, crs21Exp, dist12Exp;
    double arcDistance12, reducedLength12, area12;
    LLPoint p1, p2;
    double dist12;
    int passedCount = 0, failedCount = 0, unverifiedCount = 0;
    int setupFailureCount = 0, errorCount = 0, testCaseCount = 0;
    ErrorSet err = 0;
    int testNum = 0;
    int failed = 0;
    double distDiff;

    TestSet set;

    printf("Start testInvDist_Set1\n");

    set = newTestSet("testInvDist_Set1");

    inputfile = fopen(inputFileName, "r");

    if (inputfile == NULL)
    {
        printf("No input file %s found", inputFileName);
        set.testingError = 1;
        return set;
    }

    while (1)
    {
        //Input data: deg, deg, deg, deg, deg, deg, m, deg, m, m^2
        mm = fscanf(inputfile, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &lat1,
                &lon1, &crs12Exp, &lat2, &lon2, &crs21Exp, &dist12Exp, &arcDistance12, &reducedLength12, &area12);
        if (mm == EOF)
        {
            fclose(inputfile);
            break;
        }
        if(dist12Exp > 10740.0*1852.0) continue; //Vincenty's algorithm will fail in the nearly anti-podal case
        testNum++;

        p1.latitude = lat1 * DEG2RAD;
        p1.longitude = lon1 * DEG2RAD;
        p2.latitude = lat2 * DEG2RAD;
        p2.longitude = lon2 * DEG2RAD;
        sprintf(testname,"Test%i",testNum);

        err = invDist(p1, p2, &dist12, EPS);

        if (getMaskedError(err, getMaskAll()))
        {
            printf("Error occurred in invDist err=0x%lx\n", err);
            failedCount++;
            errorCount++;
            testCaseCount++;
            continue;
        }

        failed = 0;

        distDiff = fabs(dist12 - dist12Exp/1852.0);

        if (distDiff > NMTOL)
        {
        	printf("\n%s failed: dist12=%.20lf dist12Exp=%.20lf distDiff=%le", testname, dist12, dist12Exp, distDiff);
        	failed = 1;
        }

        if (failed)
        {
			printf("\n%s failed p1: %.20lf %.20lf  p2: %.20lf %.20lf\n",
					testname, lat1, lon1, lat2, lon2);
			printf("DistExpected: %.20lf Dist:%.20lf\n", dist12Exp/1852.0, dist12);
			failedCount++;
        } else {
        	passedCount++;
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

    printf("Finish testInvDist_Set1\n");

    return set;

}

/*
 * NAME: testInvDist_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the testInvDist function.
 *
 * 		This function runs all the test cases from all the sets that test testInvDist
 * 		and returns the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testInvDist_AllSets(TestSuite) - A test set with the following metrics:
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
TestSuite testInvDist_AllSets()
{
    TestSuite suite;
    TestSet set1;

    printf("\nStart testInvDist_AllSets\n");

    suite = newTestSuite("testInvDist_AllSets");

    set1 = testInvDist_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("Finish testInvDist_AllSets\n\n\n");

    return suite;
}

/*
 * NAME: testDirectLat_Set1
 *
 * DESCRIPTION:
 * This function is used to test the directLat function
 *
 * The test data comes from the GeographicLib project by Charles Karney. Which uses an extension of Bessel's and Helmert's method of solving
 * the direct problem. Data and code can be obtained from http://geographiclib.sourceforge.net/. A description of the data follows.
 *
 * "Each line of the test set gives 9 space delimited numbers
 *
 * latitude for point 1, lat1 (degrees, exact)
 * longitude for point 1, lon1 (degrees, always 0)
 * azimuth for point 1, azi1 (clockwise from north in degrees, exact)
 * latitude for point 2, lat2 (degrees, accurate to 10^-18 deg)
 * longitude for point 2, lon2 (degrees, accurate to 10^-18 deg)
 * azimuth for point 2, azi2 (degrees, accurate to 10^-18 deg)
 * geodesic distance from point 1 to point 2, s12 (meters, exact)
 * arc distance on the auxiliary sphere, a12 (degrees, accurate to 10^-18 deg)
 * reduced length of the geodesic, m12 (meters, accurate to 0.1 pm)
 * the area under the geodesic, S12 (meters^2, accurate to 1 mm^2)
 *
 * These are computed using as direct geodesic calculations with the given lat1, lon1, azi1, and s12. The distance
 * s12 always corresponds to an arc length a12 <= 180 deg, so the given geodesics give the shortest paths from
 * point 1 to point 2. For simplicity and without loss of generality, lat1 is chosen in [0 deg, 90 deg], lon1 is
 * taken to be zero, azi1 is chosen in [0 deg, 180 deg]. Furthermore, lat1 and azi1 are taken to be multiples
 * of 10^-12 deg and s12 is a multiple of 0.1 um in [0 m, 20003931.4586254 m]. This results lon2 in [0 deg, 180 deg] and
 * azi2 in [0 deg, 180 deg].
 *
 * The direct calculation uses an expansion of the geodesic equations accurate to f^20 (approximately 1 part in 10^50) and is computed with with Maxima's bfloats and fpprec set to 100 (so the errors in the data are probably 1/2 of the values quoted above).
 *
 * The contents of [GeodTest.dat] are as follows [GeodTest-short.dat consists of 1/50 of these tests]:
 *
 * 100000 entries randomly distributed
 * 50000 entries which are nearly antipodal
 * 50000 entries with short distances
 * 50000 entries with one end near a pole
 * 50000 entries with both ends near opposite poles
 * 50000 entries which are nearly meridional
 * 50000 entries which are nearly equatorial
 * 50000 entries running between vertices (azi1 = azi2 = 90 deg)
 * 50000 entries ending close to vertices
 *
 * (a total of 500000 entries). The values for s12 for the geodesics running between vertices are truncated
 * to a multiple of 0.1 pm and this is used to determine point 2"
 *
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testDirectLat_Set1(TestSet) - A test set with the following metrics:
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
TestSet testDirectLat_Set1()
{
    double DEG2RAD = M_PI / 180.0;
//    double RAD2DEG = 180.0 / M_PI;
//    double NMTOL = 0.03 / 100.0 / 1852.0; //0.03 cm or ~1.62e-7 nm
    double LLDEGTOL = 5.0e-8; //5e-8 degrees or 0.00018 sec
//    double AZDEGTOL = 5.0e-7; //5e-7 deg or 8.72665e-9 rad or 0.162 cm at 100 nm
//    double NMSHORTDIST = SPHERICAL_AZIMUTH_CUTOFF_DISTANCE; //5e-4 nm or 0.926 m
//    double NMLARGEDIST = 10000.0; //10000 nm
//    LLPoint origin, double course, double distance,double eps);

    FILE *inputfile;
    char inputFileName[] = FILEROOT "/geolib/src/test/resources/testLLPoint/data/GeodTest-short.dat";
    /* For a longer test set (500,000 cases) use GeodTest.dat. For a shorter test set (10,000 cases) use GeodTest-short.dat. */
    int mm;
    char testname[80];

    double lat1, lon1;
    double lat2;
    double lat2Exp, lon2Exp;
    LLPoint p1;
    double crs12;
    double dist12;
    double crs21, arcDistance12, reducedLength12, area12;

    int passedCount = 0, failedCount = 0, unverifiedCount = 0;
    int testCaseCount = 0, errorCount = 0, setupFailureCount = 0;
    ErrorSet err = 0;
    long testNum = 0;
    TestSet set;

    printf("Start testDirectLat\n");

    set = newTestSet("testDirectLat_Set1");

    inputfile = fopen(inputFileName, "r");

    if (inputfile == NULL)
    {
        printf("No input file %s found\n", inputFileName);
        set.testingError = 1;
        return set;
    }

    while (1)
    {
        //Input data: deg, deg, deg, deg, deg, deg, m, deg, m, m^2
        mm = fscanf(inputfile, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                &lat1, &lon1, &crs12, &lat2Exp, &lon2Exp, &crs21, &dist12, &arcDistance12, &reducedLength12, &area12);
        if (mm == EOF)
        {
            fclose(inputfile);
            break;
        }

        testNum++;
        sprintf(testname,"test%ld",testNum);

        p1.latitude = lat1 * DEG2RAD;
        p1.longitude = lon1 * DEG2RAD;
        crs12 *= DEG2RAD;
        dist12 /= 1852.0;

        err = directLat(p1, crs12, dist12, &lat2, EPS);

        if (getMaskedError(err, getMaskAll()))
        {
            printf("Error occurred in DirectLat err=0x%lx\n", err);
            errorCount++;
            failedCount++;
            testCaseCount++;
            continue;
        }

        lat2 /= DEG2RAD;

        if ((fabs(lat2Exp - lat2) < LLDEGTOL) || (fabs(360.0 - fabs(lat2Exp
                - lat2)) < LLDEGTOL))
        {
        	if (PRINT_PASSED_CASES) printf("%s passed\n", testname);
            passedCount++;
            testCaseCount++;
        }
        else
        {
            printf( "%s p1: %14.15lf %14.15lf  p2: %14.15lf %14.15lf  crs: %14.15lf %14.15lf  dist: %14.15lf\n",
                    testname, lat1, lon1, lat2Exp, lon2Exp, crs12, crs21, dist12);

            printf("%s failed %14.8f %14.8f\n", testname, lat2Exp, lat2);

            failedCount++;
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

    printf("Finish testDirectLat_Set1\n");

    return set;

}

/*
 * NAME: testDirectLat_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the testDirectLat function.
 *
 * 		This function runs all the test cases from all the sets that test testDirectLat
 * 		and returns the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testDirectLat_AllSets(TestSuite) - A test set with the following metrics:
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
TestSuite testDirectLat_AllSets()
{
    TestSuite suite;
    TestSet set1;

    printf("\nStart testDirectLat_AllSets\n");

    suite = newTestSuite("testDirectLat_AllSets");

    set1 = testDirectLat_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("Finish testDirectLat_AllSets\n\n\n");

    return suite;
}

/*
 * NAME: testDirectLon_Set1
 *
 * DESCRIPTION:
 * This function is used to test the directLon function
 *
 * The test data comes from the GeographicLib project by Charles Karney. Which uses an extension of Bessel's and Helmert's method of solving
 * the direct problem. Data and code can be obtained from http://geographiclib.sourceforge.net/. A description of the data follows.
 *
 * "Each line of the test set gives 9 space delimited numbers
 *
 * latitude for point 1, lat1 (degrees, exact)
 * longitude for point 1, lon1 (degrees, always 0)
 * azimuth for point 1, azi1 (clockwise from north in degrees, exact)
 * latitude for point 2, lat2 (degrees, accurate to 10^-18 deg)
 * longitude for point 2, lon2 (degrees, accurate to 10^-18 deg)
 * azimuth for point 2, azi2 (degrees, accurate to 10^-18 deg)
 * geodesic distance from point 1 to point 2, s12 (meters, exact)
 * arc distance on the auxiliary sphere, a12 (degrees, accurate to 10^-18 deg)
 * reduced length of the geodesic, m12 (meters, accurate to 0.1 pm)
 * the area under the geodesic, S12 (meters^2, accurate to 1 mm^2)
 *
 * These are computed using as direct geodesic calculations with the given lat1, lon1, azi1, and s12. The distance
 * s12 always corresponds to an arc length a12 <= 180 deg, so the given geodesics give the shortest paths from
 * point 1 to point 2. For simplicity and without loss of generality, lat1 is chosen in [0 deg, 90 deg], lon1 is
 * taken to be zero, azi1 is chosen in [0 deg, 180 deg]. Furthermore, lat1 and azi1 are taken to be multiples
 * of 10^-12 deg and s12 is a multiple of 0.1 um in [0 m, 20003931.4586254 m]. This results lon2 in [0 deg, 180 deg] and
 * azi2 in [0 deg, 180 deg].
 *
 * The direct calculation uses an expansion of the geodesic equations accurate to f^20 (approximately 1 part in 10^50) and is computed with with Maxima's bfloats and fpprec set to 100 (so the errors in the data are probably 1/2 of the values quoted above).
 *
 * The contents of [GeodTest.dat] are as follows [GeodTest-short.dat consists of 1/50 of these tests]:
 *
 * 100000 entries randomly distributed
 * 50000 entries which are nearly antipodal
 * 50000 entries with short distances
 * 50000 entries with one end near a pole
 * 50000 entries with both ends near opposite poles
 * 50000 entries which are nearly meridional
 * 50000 entries which are nearly equatorial
 * 50000 entries running between vertices (azi1 = azi2 = 90 deg)
 * 50000 entries ending close to vertices
 *
 * (a total of 500000 entries). The values for s12 for the geodesics running between vertices are truncated
 * to a multiple of 0.1 pm and this is used to determine point 2"
 *
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testDirect_Set1(TestSet) - A test set with the following metrics:
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
TestSet testDirectLon_Set1()
{
    double DEG2RAD = M_PI / 180.0;
//    double RAD2DEG = 180.0 / M_PI;
//    double NMTOL = 0.03 / 100.0 / 1852.0; //0.03 cm or ~1.62e-7 nm
    double LLDEGTOL = 5.0e-8; //5e-8 degrees or 0.00018 sec
//    double AZDEGTOL = 5.0e-7; //5e-7 deg or 8.72665e-9 rad or 0.162 cm at 100 nm
//    double NMSHORTDIST = SPHERICAL_AZIMUTH_CUTOFF_DISTANCE; //5e-4 nm or 0.926 m
//    double NMLARGEDIST = 10000.0; //10000 nm
//    LLPoint origin, double course, double distance,double eps);

    FILE *inputfile;
    char inputFileName[] = FILEROOT "/geolib/src/test/resources/testLLPoint/data/GeodTest-short.dat";
    /* For a longer test set (500,000 cases) use GeodTest.dat. For a shorter test set (10,000 cases) use GeodTest-short.dat. */
    int mm;
    char testname[80];

    double lat1, lon1;
    double lon2;
    double lat2Exp, lon2Exp;
    LLPoint p1;
    double crs12;
    double dist12;
    double crs21, arcDistance12, reducedLength12, area12;

    int passedCount = 0, failedCount = 0, unverifiedCount = 0;
    int testCaseCount = 0, errorCount = 0, setupFailureCount = 0;
    ErrorSet err = 0;
    long testNum = 0;
    TestSet set;

    printf("Start testDirectLon\n");

    set = newTestSet("testDirectLon_Set1");

    inputfile = fopen(inputFileName, "r");

    if (inputfile == NULL)
    {
        printf("No input file %s found\n", inputFileName);
        set.testingError = 1;
        return set;
    }

    while (1)
    {
        //Input data: deg, deg, deg, deg, deg, deg, m, deg, m, m^2
        mm = fscanf(inputfile, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                &lat1, &lon1, &crs12, &lat2Exp, &lon2Exp, &crs21, &dist12, &arcDistance12, &reducedLength12, &area12);
        if (mm == EOF)
        {
            fclose(inputfile);
            break;
        }
        if((fabs(lat1) > 89.9) || (fabs(lat2Exp) > 89.9)) continue;
        //The distance between longitudes becomes smaller near the poles so Vincenty's algorithm is only accurate to 10^-5 deg there

        testNum++;
        sprintf(testname,"test%ld",testNum);

        p1.latitude = lat1 * DEG2RAD;
        p1.longitude = lon1 * DEG2RAD;
        crs12 *= DEG2RAD;
        dist12 /= 1852.0;

        err = directLon(p1, crs12, dist12, &lon2, EPS);

        if (getMaskedError(err, getMaskAll()))
        {
            printf("Error occurred in directLon err=0x%lx\n", err);
            errorCount++;
            failedCount++;
            testCaseCount++;
            continue;
        }

        lon2 /= DEG2RAD;

        if ((fabs(lon2Exp - lon2) < LLDEGTOL) || (fabs(360.0 - fabs(lon2Exp
                - lon2)) < LLDEGTOL))
        {
        	if (PRINT_PASSED_CASES) printf("%s passed\n", testname);
            passedCount++;
            testCaseCount++;
        }
        else
        {
            printf( "%s p1: %14.15lf %14.15lf  p2: %14.15lf %14.15lf  crs: %14.15lf %14.15lf  dist: %14.15lf\n",
                    testname, lat1, lon1, lat2Exp, lon2Exp, crs12, crs21, dist12);

            printf("%s failed %14.8f %14.8f\n", testname, lon2Exp, lon2);

            failedCount++;
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

    printf("Finish testDirectLon_Set1\n");

    return set;

}

/*
 * NAME: testDirectLon_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the testDirectLon function.
 *
 * 		This function runs all the test cases from all the sets that test testDirectLon
 * 		and returns the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testDirectLon_AllSets(TestSuite) - A test set with the following metrics:
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
TestSuite testDirectLon_AllSets()
{
    TestSuite suite;
    TestSet set1;

    printf("\nStart testDirectLon_AllSets\n");

    suite = newTestSuite("testDirectLon_AllSets");

    set1 = testDirectLon_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("Finish testDirectLon_AllSets\n\n\n");

    return suite;
}

/*
 * NAME: testDirect_Set1
 *
 * DESCRIPTION:
 * 		This function is used to test the direct function
 *
 * The test data comes from the GeographicLib project by Charles Karney. Which uses an extension of Bessel's and Helmert's method of solving
 * the direct problem. Data and code can be obtained from http://geographiclib.sourceforge.net/. A description of the data follows.
 *
 * "Each line of the test set gives 9 space delimited numbers
 *
 * latitude for point 1, lat1 (degrees, exact)
 * longitude for point 1, lon1 (degrees, always 0)
 * azimuth for point 1, azi1 (clockwise from north in degrees, exact)
 * latitude for point 2, lat2 (degrees, accurate to 10^-18 deg)
 * longitude for point 2, lon2 (degrees, accurate to 10^-18 deg)
 * azimuth for point 2, azi2 (degrees, accurate to 10^-18 deg)
 * geodesic distance from point 1 to point 2, s12 (meters, exact)
 * arc distance on the auxiliary sphere, a12 (degrees, accurate to 10^-18 deg)
 * reduced length of the geodesic, m12 (meters, accurate to 0.1 pm)
 * the area under the geodesic, S12 (meters^2, accurate to 1 mm^2)
 *
 * These are computed using as direct geodesic calculations with the given lat1, lon1, azi1, and s12. The distance
 * s12 always corresponds to an arc length a12 <= 180 deg, so the given geodesics give the shortest paths from
 * point 1 to point 2. For simplicity and without loss of generality, lat1 is chosen in [0 deg, 90 deg], lon1 is
 * taken to be zero, azi1 is chosen in [0 deg, 180 deg]. Furthermore, lat1 and azi1 are taken to be multiples
 * of 10^-12 deg and s12 is a multiple of 0.1 um in [0 m, 20003931.4586254 m]. This results lon2 in [0 deg, 180 deg] and
 * azi2 in [0 deg, 180 deg].
 *
 * The direct calculation uses an expansion of the geodesic equations accurate to f^20 (approximately 1 part in 10^50) and is computed with with Maxima's bfloats and fpprec set to 100 (so the errors in the data are probably 1/2 of the values quoted above).
 *
 * The contents of [GeodTest.dat] are as follows [GeodTest-short.dat consists of 1/50 of these tests]:
 *
 * 100000 entries randomly distributed
 * 50000 entries which are nearly antipodal
 * 50000 entries with short distances
 * 50000 entries with one end near a pole
 * 50000 entries with both ends near opposite poles
 * 50000 entries which are nearly meridional
 * 50000 entries which are nearly equatorial
 * 50000 entries running between vertices (azi1 = azi2 = 90 deg)
 * 50000 entries ending close to vertices
 *
 * (a total of 500000 entries). The values for s12 for the geodesics running between vertices are truncated
 * to a multiple of 0.1 pm and this is used to determine point 2"
 *
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testDirect_Set1(TestSet) - A test set with the following metrics:
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
TestSet testDirect_Set1()
{
    double DEG2RAD = M_PI / 180.0;
//    double RAD2DEG = 180.0 / M_PI;
//    double NMTOL = 0.03 / 100.0 / 1852.0; //0.03 cm or ~1.62e-7 nm
//    double LLDEGTOL = 5.0e-8; //5e-8 degrees or 0.00018 sec
//    double AZDEGTOL = 5.0e-7; //5e-7 deg or 8.72665e-9 rad or 0.162 cm at 100 nm
//    double NMSHORTDIST = SPHERICAL_AZIMUTH_CUTOFF_DISTANCE; //5e-4 nm or 0.926 m
//    double NMLARGEDIST = 10000.0; //10000 nm
//    LLPoint origin, double course, double distance,double eps);

    FILE *inputfile;
    char inputFileName[] = FILEROOT "/geolib/src/test/resources/testLLPoint/data/GeodTest.dat";
    /* For a longer test set (500,000 cases) use GeodTest.dat. For a shorter test set (10,000 cases) use GeodTest-short.dat. */
    int mm;
    char testname[80];

    double lat1, lon1;
    double lat2, lon2;
    double lat2Exp, lon2Exp;
    LLPoint p1, p2, p2Exp;
    double crs12;
    double dist12;
    double crs21, arcDistance12, reducedLength12, area12;

    double distErr;

    int passedCount = 0, failedCount = 0, unverifiedCount = 0;
    int testCaseCount = 0, errorCount = 0, setupFailureCount = 0;
    ErrorSet err = 0;
    long testNum = 0;
    TestSet set;

    printf("Start testDirect\n");

    set = newTestSet("testDirect_Set1");

    inputfile = fopen(inputFileName, "r");

    if (inputfile == NULL)
    {
        printf("No input file %s found\n", inputFileName);
        set.testingError = 1;
        return set;
    }

    while (1)
    {
        //Input data: deg, deg, deg, deg, deg, deg, m, deg, m, m^2
        mm = fscanf(inputfile, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                &lat1, &lon1, &crs12, &lat2Exp, &lon2Exp, &crs21, &dist12, &arcDistance12, &reducedLength12, &area12);
        if (mm == EOF)
        {
            fclose(inputfile);
            break;
        }

        testNum++;
        sprintf(testname,"test%ld",testNum);

        p1.latitude = lat1 * DEG2RAD;
        p1.longitude = lon1 * DEG2RAD;
        p2Exp.latitude = lat2Exp * DEG2RAD;
        p2Exp.longitude = lon2Exp * DEG2RAD;
        crs12 *= DEG2RAD;
        dist12 /= 1852.0;

        err = direct(p1, crs12, dist12, &p2, EPS);

        if (getMaskedError(err, getMaskAll()))
        {
            printf("Error occurred in direct err=0x%lx\n", err);
            errorCount++;
            testCaseCount++;
            failedCount++;
            continue;
        }

        lat2 = p2.latitude / DEG2RAD;
        lon2 = p2.longitude / DEG2RAD;
        //        printf("%s p1: %13.8f %13.8f crs: %11.6f dist: %10.5f p2: %13.8f %13.8f\n",
        //        		testname, lat1, lon1, crs, dist, lat2, lon2);


        if (isEqualLLPoint(p2, p2Exp))
        {
        	if (PRINT_PASSED_CASES) printf("%s passed\n", testname);
            passedCount++;
            testCaseCount++;
        }
        else
        {
            printf( "%s p1: %14.15lf %14.15lf  p2: %14.15lf %14.15lf  crs: %14.15lf %14.15lf  dist: %14.15lf\n",
                    testname, lat1, lon1, lat2Exp, lon2Exp, crs12, crs21, dist12);

            err |= invDist(p2Exp, p2, &distErr, EPS);

            if (getMaskedError(err, getMaskAll()))
            {
                printf("Error occurred in post-direct err=0x%lx\n", err);
            }

            printf("%s failed p2Exp: %14.15f %14.15f  p2: %14.15f %14.15f  distErr: %18.15e\n",
                    testname, lat2Exp, lon2Exp, lat2, lon2, distErr);

            failedCount++;
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

    printf("Finish testDirect_Set1\n");

    return set;

}

/*
 * NAME: testDirect_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the direct function.
 *
 * 		This function runs all the test cases from all the sets that test direct
 * 		and returns the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testDirect_AllSets(TestSuite) - A test set with the following metrics:
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
TestSuite testDirect_AllSets()
{
    TestSuite suite;
    TestSet set1;

    printf("\nStart testDirect_AllSets\n");

    suite = newTestSuite("testDirect_AllSets");

    set1 = testDirect_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("Finish testDirect_AllSets\n\n\n");

    return suite;
}

static void smallDistInverse(LLPoint p1, LLPoint p2, double* az, double* dist)
{
    double dlon;
    minSubtendedAngle(p2.longitude,p1.longitude, &dlon);
    double dlat = p2.latitude - p1.latitude;
    double N = findN(p1);
    double M = findM(p1);
    double coslat = cos(p1.latitude);
    //    double tanAlpha = (N / M) * (cos(p1.latitude) * dlon / dlat);

    //    NdivM = coslat*coslat/BdivA/BdivA + sinlat*sinlat;
    //    tanAlpha = NdivM*(cosl((long double) p1.latitude)*dlambda/dphi);
    //    printf("tan(az') = %.16e\n",(double)tanAlpha);
    if (NULL != az)
    {

        *az = atan2(N * coslat * dlon, M * dlat);

        if (*az < 0)
        {
            *az += M_2PI;
        }
    }
    if (NULL != dist)
    {
        *dist = sqrt(pow(M * dlat, 2.0) + pow(N * coslat * dlon, 2.0));
    }

}

TestSet testDirectInverseConsistency()
{
    double DEG2RAD = M_PI / 180.0;
    FILE *outputfile;
    /* For a longer test set (500,000 cases) use GeodTest.dat. For a shorter test set (10,000 cases) use GeodTest-short.dat. */
	char outputFileName[] = FILEROOT "/geolib/src/test/resources/testLLPoint/data/consistencyTestNGSVincenty.csv";


    double lat1Exp, lon1Exp;
    double crs12Exp, dist12Exp;
    LLPoint p1Exp, p2, p1;
    double crs12, crs21, dist12;
    int passedCount = 0, failedCount = 0, unverifiedCount = 0;
    int setupFailureCount = 0, errorCount = 0, testCaseCount = 0;
    ErrorSet err = 0;
    double crs12Diff, distDiff, distErr;

    TestSet set;
    set = newTestSet("testDirectInverseConsistency");


    outputfile = fopen(outputFileName, "w");

    int i, j;
    long newSeed = 13487209;
    srand(newSeed);  //Initialize the random number generator

    int NGSErr = 0;

    for(i = 0; i < 1e5; i++)
    {
    	err = 0;
    	NGSErr = 0;
        //generate input data for Arc
        lat1Exp = randLat()*DEG2RAD;
        lon1Exp = randLon()*DEG2RAD;
        dist12Exp = (double) ((rand() % 10000));
        crs12Exp = randAzimuth()*DEG2RAD;

        p1Exp.latitude = lat1Exp;
        p1Exp.longitude = lon1Exp;
        crs12 = crs12Exp;
        dist12 = dist12Exp;
        p1 = p1Exp;

        for(j=0;j < 10000; j++){
        	err |= direct(p1, crs12, dist12, &p2, EPS);
        	err |= inverse(p2, p1, &crs21, &crs12, &dist12, EPS);
        	if(dist12 == 0 && crs12 == 0 && crs21 == 0)
        		NGSErr = 1;
        	err |= direct(p2, crs21, dist12, &p1, EPS);
        	err |= inverse(p1, p2, &crs12, &crs21, &dist12, EPS);
        	if(dist12 == 0 && crs12 == 0 && crs21 == 0)
        		NGSErr = 1;

        }
        p1.longitude = modlon(p1.longitude);
    	smallDistInverse(p1, p1Exp, NULL, &distErr);
    	minSubtendedAngle(crs12, crs12Exp, &crs12Diff);
        distDiff = fabs(dist12 - dist12Exp);


//        printf("%.18lf %.18lf, %.18lf, %.18lf\n", lat1Exp/DEG2RAD, lon1Exp/DEG2RAD, crs12Exp/DEG2RAD, dist12Exp);
//        printf("%.18lf %.18lf, %.18lf, %.18lf\n", p1.latitude/DEG2RAD, p1.longitude/DEG2RAD, crs12/DEG2RAD, dist12);
//        printf("%.18e\n", distErr);

        fprintf(outputfile,"%.16e, %.16e, %.16e, %.16e, %.16e, %.16e, %.16e, %.16e, %.16e, %.16e, %.16e, %d\n",
        		lat1Exp/DEG2RAD, lon1Exp/DEG2RAD, crs12Exp/DEG2RAD, dist12Exp,
        		p1.latitude/DEG2RAD, p1.longitude/DEG2RAD, crs12/DEG2RAD, dist12,
        		distErr, crs12Diff/DEG2RAD, distDiff, NGSErr);

        testCaseCount++;

    }


    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("Finish testDirectInverseConsistency\n");

    return set;

}

TestSet testDirectInverseMathematicaData()
{
    double DEG2RAD = M_PI / 180.0;
    FILE *inputfile;
    FILE *outputfileInv;
    FILE *outputfileDir;
    int mm;
    char inputFileName[] = FILEROOT "/geolib/src/test/resources/testLLPoint/data/mathematicaData.dat";
    /* For a longer test set (500,000 cases) use GeodTest.dat. For a shorter test set (10,000 cases) use GeodTest-short.dat. */
	char outputFileNameInv[] = FILEROOT "/geolib/src/test/resources/testLLPoint/data/MathematicaDataOutputInvNGSVincenty.csv";
	char outputFileNameDir[] = FILEROOT "/geolib/src/test/resources/testLLPoint/data/MathematicaDataOutputDirNGSVincenty.csv";

    double lat1, lon1, lat2, lon2;
    double lat2Exp, lon2Exp;
    double crs12Exp, crs21Exp, dist12Exp;
    LLPoint p1, p2, p2Exp;
    double crs12, crs21, dist12;
    int passedCount = 0, failedCount = 0, unverifiedCount = 0;
    int setupFailureCount = 0, errorCount = 0, testCaseCount = 0;
    ErrorSet err = 0;
    double crs12Diff, crs21Diff, distDiff, distErr;

    TestSet set;
    set = newTestSet("testDirectInverseMathematicaData");


    outputfileInv = fopen(outputFileNameInv, "w");
    outputfileDir = fopen(outputFileNameDir, "w");



    inputfile = fopen(inputFileName, "r");
    if (inputfile == NULL)
    {
        printf("No input file %s found", inputFileName);
        set.testingError = 1;
        return set;
    }
    while (1)
    {
    	err = 0;
        //Mathematica Input data: deg, deg, deg, deg, deg, deg, m
        mm = fscanf(inputfile, "%lf %lf %lf %lf %lf %lf %lf",
        		&lat1, &lon1, &crs12Exp, &lat2, &lon2, &crs21Exp, &dist12Exp);
//        //Karney Input data: deg, deg, deg, deg, deg, deg, m, deg, m, m^2
//        mm = fscanf(inputfile, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
//        		&lat1, &lon1, &crs12Exp, &lat2, &lon2, &crs21Exp, &dist12Exp, &arcDistance12, &reducedLength12, &area12);
        if (mm == EOF)
        {
            fclose(inputfile);
            break;
        }

        p1.latitude = lat1 * DEG2RAD;
        p1.longitude = lon1 * DEG2RAD;
        p2.latitude = lat2 * DEG2RAD;
        p2.longitude = lon2 * DEG2RAD;

        err = inverse(p1, p2, &crs12, &crs21, &dist12, EPS);

        crs12Diff = fabs(crs12 - crs12Exp*DEG2RAD);
        crs21Diff = fabs(crs21 - fmod(crs21Exp + 180.0,360.0)*DEG2RAD);
        distDiff = fabs(dist12 - dist12Exp/1852.0);

        fprintf(outputfileInv,"%.16e, %.16e, %.16e, %.16e, %.16e, %.16e, %.16e, %.16e, %.16e, %.16e, %.16e, %.16e, %ld\n",
        		lat1, lon1, crs12/DEG2RAD, lat2, lon2, dist12, dist12Exp/1852, crs12Exp, crs21Exp, crs12Diff/DEG2RAD, crs21Diff/DEG2RAD, distDiff, (long)err);

        testCaseCount++;

    }

    inputfile = fopen(inputFileName, "r");
    if (inputfile == NULL)
    {
        printf("No input file %s found", inputFileName);
        set.testingError = 1;
        return set;
    }
    while (1)
    {
    	err = 0;
        //Input data: deg, deg, deg, deg, deg, deg, m, deg, m, m^2
        mm = fscanf(inputfile, "%lf %lf %lf %lf %lf %lf %lf",
        		&lat1, &lon1, &crs12Exp, &lat2Exp, &lon2Exp, &crs21Exp, &dist12);
        if (mm == EOF)
        {
            fclose(inputfile);
            break;
        }

        p1.latitude = lat1 * DEG2RAD;
        p1.longitude = lon1 * DEG2RAD;

        err = direct(p1, crs12Exp*DEG2RAD, dist12/1852, &p2, EPS);

        lat2 = p2.latitude/DEG2RAD;
        lon2 = p2.longitude/DEG2RAD;

        p2Exp.latitude = lat2Exp*DEG2RAD;
        p2Exp.longitude = lon2Exp*DEG2RAD;

        smallDistInverse(p2, p2Exp, NULL, &distErr);
//        invDist(p2, p2Exp, &distErr, EPS);

        fprintf(outputfileDir,"%.16e, %.16e, %.16e, %.16e, %.16e, %.16e, %.16e, %.16e, %.16e, %ld\n",
        		lat1, lon1, crs12Exp, lat2, lon2, dist12/1852, lat2Exp, lon2Exp, distErr, (long)err);

        testCaseCount++;

    }

    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("Finish testDirectInverseMathematicaData\n");

    return set;

}

} // namespace
