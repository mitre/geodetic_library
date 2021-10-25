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
#include <stdarg.h>
#include <float.h>
#include <limits.h>
#include "Geolib.h"
#include "testGeolib.h"

namespace geolib_idealab {





/*
 * NAME: testAddPointToLLPointSet_Set1
 *
 * DESCRIPTION:
 * 		This function is used to test the addPtToPtSet function.
 *
 * 		This function runs the test data created by Juan Amezcua and uses
 * 		his original test plan to execute the actual tests.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testAddPointToLLPointSet_Set1(TestSet) - A test set with the folling metrics:
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
TestSet testAddPointToLLPointSet_Set1()
{


    int passedCount=0, failedCount=0;
    int errorCount = 0;
    int unverifiedCount = 0;
    int setupFailureCount = 0;
    int testCaseCount = 0;
    ErrorSet err=0;
    TestSet set;
    char testName[100];

    LLPointSet llpSet;
    LLPoint pt1,pt2;


    set = newTestSet("testAddPointToLLPointSet_Set1");
    printf("\n\n%s\n", set.name.c_str());

    //test 1
    err = 0;
    testCaseCount++;
    sprintf(testName, "TEST %i", testCaseCount);

    llpSet = createPtSet();
    err |= createPt(&pt1, 1.0, 2.0);
    err |= createPt(&pt2, 1.0, 2.0);

    if (err) {
    	printf("\n%s - Error occurred during setup: 0x%lx", testName, err);
    	setupFailureCount++;
    	errorCount++;
    } else {
    	err |= addPtToPtSet(&llpSet, &pt1);

    	if (err){
    		printf("\n%s - Error occurred with addPointToLLPointSet: 0x%lx", testName, err);
    		failedCount++;
    		errorCount++;
    	} else {
    		//exact comparsion, don't use tolerances (i.e. neighborhoods) here
    		if (llpSet.elements[0]->latitude != pt2.latitude || llpSet.elements[0]->longitude != pt2.longitude){
    			//failure - point in LLPointSet doesn't exactly match
    			printf("\n\n%s FAILED", testName);
    			printf("\nexpected point = (%.20lf,%.20lf)", pt2.latitude, pt2.longitude);
    			printf("\nactual point = (%.20lf,%.20lf)", llpSet.elements[0]->latitude, llpSet.elements[0]->longitude);
    			failedCount++;
    		} else {
    			passedCount++;
    		}
    	}
    }

    //test 2
    err = 0;
    testCaseCount++;
    sprintf(testName, "TEST %i", testCaseCount);

    llpSet = createPtSet();
    err |= createPt(&pt1, 1.0, 2.0);
    err |= addPtToPtSet(&llpSet, &pt1);

    err |= createPt(&pt2, 3.0, 4.0);

    if (err) {
    	printf("\n%s - Error occurred during setup: 0x%lx", testName, err);
    	setupFailureCount++;
    	errorCount++;
    } else {
    	err |= addPtToPtSet(&llpSet, &pt2);

    	if (err){
    		printf("\n%s - Error occurred with addPointToLLPointSet: 0x%lx", testName, err);
    		failedCount++;
    		errorCount++;
    	} else {
    		//exact comparsion, don't use tolerances (i.e. neighborhoods) here
    		if (llpSet.elements[1]->latitude != pt2.latitude || llpSet.elements[1]->longitude != pt2.longitude){
    			//failure - point in LLPointSet doesn't exactly match
    			printf("\n\n%s FAILED", testName);
    			printf("\nexpected point = (%.20lf,%.20lf)", pt2.latitude, pt2.longitude);
    			printf("\nactual point = (%.20lf,%.20lf)", llpSet.elements[0]->latitude, llpSet.elements[0]->longitude);
    			failedCount++;
    		} else {
    			passedCount++;
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

    printf("Finish %s\n\n\n", set.name.c_str());

    return set;

}


} //namespace
