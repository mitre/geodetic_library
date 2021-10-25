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
#include "testGeolib.h"

namespace geolib_idealab {



/*
 * NAME: testFindSetMaxAndMin_Set1
 *
 * DESCRIPTION:
 * 		This function is used to test the findSetMaxAndMin function.
 *
 * 		This function runs the test data created by Juan Amezcua and uses
 * 		his original test plan to execute the actual tests.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testFindSetMaxAndMin_Set1(TestSet) - A test set with the following metrics:
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
TestSet testFindSetMaxAndMin_Set1()
{
//	double myTol = 1.37e-9;//nautical miles
//	double myEps = 1.0e-20;


    int passedCount=0, failedCount=0;
    int errorCount = 0;
    int unverifiedCount = 0;
    int setupFailureCount = 0;
    int testCaseCount = 0;
    ErrorSet  err=0;

    int size = 10;
    double values[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    double max;
    double min;
    int i = 0;

    TestSet set;

    set = newTestSet("testFindSetMaxAndMin_Set1");

    printf("\nStart %s\n", set.name.c_str());


    err |= findSetMaxAndMin(values, size, &max, &min);

	printf("Set = {  ");
    for (i = 0; i < size; i++)
    {
    	printf("%lf ", values[i]);
    }
    printf("}\n");
    printf("max=%lf\n", max);
    printf("min=%lf\n", min);

    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("\nFinish %s\n\n\n", set.name.c_str());

    return set;

}

/*
 * NAME: testFindSetMaxAndMin_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the complex boundary's findSetMaxAndMin function.
 *
 * 		This function runs all the test cases and returns
 * 		the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testFindSetMaxAndMin_AllSets(TestSuite) - A test suite with the folling metrics:
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
TestSuite testFindSetMaxAndMin_AllSets()
{
	TestSuite suite;
	TestSet set1;

    suite = newTestSuite("testFindSetMaxAndMin_AllSets");

    printf("\nStart %s\n", suite.name.c_str());

    set1 = testFindSetMaxAndMin_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("\nFinish %s\n\n\n", suite.name.c_str());

    return suite;
}

} //namespace
