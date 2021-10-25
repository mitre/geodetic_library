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

#ifndef TESTARC_C_
#define TESTARC_C_
#endif

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

#define MAX_LINE_SIZE 1000
#define MAX_LINE_SIZE2 100
#define nInputs 8
#define nInputs2 21
#define nInputs3 15

#ifndef FILEROOT
#define FILEROOT "."
#endif



/*
 * NAME: testInitArcIntx_Set1
 *
 * DESCRIPTION:
 * 		This function is used to test the initArcIntx function.
 *
 * 		This function runs the test data created by Balakrishna Babu and uses
 * 		his original test plan to execute the actual tests.
 *
 *   	The approach for testing initArcIntx is as follows:
 *
 *     	1. Select center (c1) and radius (r1) for Arc 1.
 *    	2. Select center (c2) for Arc 2.
 *    	3. Calculate the sphere distance (d) between the two centers using sphereInvDist.
 *    	4. Choose radius (r2) for Arc 2 based on r1 and d as follows:
 *     	5. Compute abs(d-r1) say (dm) and d+r1 say (dp)
 *     	5. Choose r2 = {0.5*dm,dm-eps,dm,dm+eps,d,d+r1/2,dp-eps,dp,dp+eps,1.1*dp} where eps
 *     	   is a small distance (use NMTOL)
 *     	6. Two center points and their radii form the input.
 *     	7. Execute the function and using the results perform validation as follows:
 *     	8. For r2 = 0.5*dm and 1.1*dp result should have zero intersection points.
 *     	9. For r2 = dm-eps,dm,dm+eps,dp-eps,dp,dp+eps there should be one intersection
 *     	   point.
 *     	10. For r2 = d and d < r1, if r2+d < r1-eps then result should be zero intersection
 *     	    points, if r1-eps < r2+d < r1+eps then result should be one intersection point,
 *     	    if r2+d > r1+eps then result should be two intersection points.
 *     	11. For r2 = d and d > r1, then result should be two intersection points.
 *     	12. For r2 = d+r1/2, there should be two intersection points.
 *     	13. Calculate the distances from each intersection point to both centers and they
 *     	    should match the radii r1 and r2.
 *
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testInitArcIntx_Set1(TestSet) - A test set with the folling metrics:
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
TestSet testInitArcIntx_Set1()
{
	double DEG2RAD = M_PI / 180.0;
	double NMTOL = 0.03 / 100.0 / 1852.0;  //0.03 cm or ~1.62e-7 nm

    int testCase = 0;
    char testname[80];
    double latC1, lonC1, latC2, lonC2;
    double r1, r2, d12;  //radii of Arcs 1 and 2, distance between centers
    double r2A[10],dm,dp;
    LLPoint c1, c2;
    LLPointPair intx;
    int nX, i, j;
    double latx1, lonx1, latx2, lonx2;
    double distC1X1, distC1X2, distC2X1, distC2X2, distX1X2;
    double usedROC = 0;
    int passedCount=0, failedCount=0;
    int errorCount = 0;
    int unverifiedCount = 0;
    int setupFailureCount = 0;
    int testCaseCount = 0;
    ErrorSet err = 0;
    long seed = 20080128;
	double crs12;
	TestSet set;

	srand(seed);

    printf("Start testInitArcIntx_Set1\n");

    set = newTestSet("testInitArcIntx_Set1");

    while( testCaseCount <= 10000 )
    {
	err = 0;
        testCase = testCaseCount;
        //Select center and radius for Arc 1
        latC1 = randLat();
        lonC1 = randLon();
        r1 = 0.2 * randDist();  //about 1080 nm

        //Select center for Arc 2 in the same hemisphere
        latC2 = latC1 + randLat();
        lonC2 = lonC1 + 0.5*randLon();

        //Normalize lonC2 to within [-180,180]
        if (lonC2 > 180.0)
        {
            lonC2 -= 360.0;
        }
        else if (lonC2 <= -180.0)
        {
            lonC2 += 360.0;
        }

        //Normalize latC2 to within [-90,90] and adjust lonC2
        if (latC2 > 90.0)
        {
            latC2 = 180.0 - latC2;
            lonC2 += 180.0;
        }
        else if (latC2 < -90.0)
        {
            latC2 = -180.0 - latC2;
            lonC2 += 180.0;
        }

        //Once again Normalize lonC2 to within [-180,180]
        if (lonC2 > 180.0)
        {
            lonC2 -= 360.0;
        }
        else if (lonC2 <= -180.0)
        {
            lonC2 += 360.0;
        }

        //Set c1 and c2 lat/lon values in radians
        c1.latitude = latC1 * DEG2RAD;
        c1.longitude = lonC1 * DEG2RAD;
        c2.latitude = latC2 * DEG2RAD;
        c2.longitude = lonC2 * DEG2RAD;

        //Calculate spheroidal distance between center points
        d12 = sphereInvDist(c1,c2,NULL);

        //Set up r2 array
        //{0.5*dm,dm-eps,dm,dm+eps,d,d+r1/2,dp-eps,dp,dp+eps,1.1*dp}
        dm = fabs(d12-r1);
        dp = d12+r1;
        r2A[0] = 0.5*dm;
        r2A[1] = dm-0.5*NMTOL;
        r2A[2] = dm;
        r2A[3] = dm+0.5*NMTOL;
        r2A[4] = d12;
        r2A[5] = d12+0.5*r1;
        r2A[6] = dp-0.5*NMTOL;
        r2A[7] = dp;
        r2A[8] = dp+0.5*NMTOL;
        r2A[9] = 1.1*dp;

        for (j=0; j<10; j++)
        {
            err = 0;
            r2 = r2A[j];


	    //Convert test case to Babu's original numbering format
	    sprintf(testname, "TEST%i", (testCaseCount / 10) + 1);
	    i = testCaseCount % 10;

            //Run the tests starting here

            d12 = sphereInvDist(c1,c2,NULL);
            crs12 = sphereInvCrs(c1,c2,1.0e-12);

            intx[0].latitude = 0.0;
            intx[0].longitude = 0.0;
            intx[1].latitude = 0.0;
            intx[1].longitude = 0.0;
            err |= initArcIntx(c1, r1, c2, r2, intx, &nX, NULL, NMTOL);

            if ( getMaskedError(err, getMaskAll()) )
            {
            	printf("Error occurred in initArcIntx err=0x%lx\n",err);

            	printf("%s-%-d failed nX=%d  X1: %14.8f %14.8f  X2: %14.8f %14.8f\n",
                       testname,i,nX,latx1,lonx1,latx2,lonx2);
                printf("       arc1: (%20.15f, %20.15f, %20.15f);  arc2: (%20.15f, %20.15f, %20.15f)\n",latC1*180./M_PI,
                       lonC1*180./M_PI, r1, latC2*180./M_PI, lonC2*180./M_PI, r2);
                {
                    double sphDist, num, den, distp;
                    double r1p = r1/SPHERE_RADIUS_NMI, r2p = r2/SPHERE_RADIUS_NMI;
                    sphDist = sphereInvDist(c1,c2,NULL);
                    distp = sphDist/SPHERE_RADIUS_NMI;
                    num = cos(r2p)-cos(distp)*cos(r1p);
                    den = sin(distp)*sin(r1p);

                    printf("       dist = %.8e, num = %.8e, den = %.8e, frac = %.8e\n",sphDist,num,den,num/den);
                    printf("       outDist = %.15e, inDist = %.15e\n", sphDist-(r1+r2), r1>r2?r1-(sphDist+r2):r2-(sphDist+r1));
                }
                errorCount++;
                failedCount++;
                testCaseCount++;
                continue;
            }

            latx1 = intx[0].latitude / DEG2RAD;
            lonx1 = intx[0].longitude / DEG2RAD;
            latx2 = intx[1].latitude / DEG2RAD;
            lonx2 = intx[1].longitude / DEG2RAD;

            if (nX == 2)
            {
                distC1X1 = sphereInvDist(c1,intx[0],NULL);
                distC2X1 = sphereInvDist(c2,intx[0],NULL);
                distC1X2 = sphereInvDist(c1,intx[1],NULL);
                distC2X2 = sphereInvDist(c2,intx[1],NULL);
                distX1X2 = sphereInvDist(intx[0],intx[1],NULL);
                if ( (fabs(distC1X1-r1)<=TESTTOL) &&
                        (fabs(distC2X1-r2)<=TESTTOL) &&
                        (fabs(distC1X2-r1)<=TESTTOL) &&
                        (fabs(distC2X2-r2)<=TESTTOL) &&
                        (distX1X2>TESTTOL) )  //two distinct points
                {
                	if (PRINT_PASSED_CASES) printf("%s-%-d nX=%d  r1=%14.8f  r2=%14.8f  d12=%14.8f  usedROC=%14.8f\n",
                           testname,i ,nX,r1,r2,d12,usedROC);
                    if (PRINT_PASSED_CASES) printf("%s-%-d passed nX=%d  X1: %14.8f %14.8f  X2: %14.8f %14.8f\n",
                           testname,i,nX,latx1,lonx1,latx2,lonx2);
                    if (PRINT_PASSED_CASES) printf("%s-%-d passed C1X1=%14.8f C2X1=%14.8f C1X2=%14.8f C2X2=%14.8f X1X2=%14.8f\n",
                           testname,i,distC1X1,distC2X1,distC1X2,distC2X2,distX1X2);
                    if (PRINT_PASSED_CASES) printf("%s-%-d passed C1X1Err=%16.8e C2X1Err=%16.8e C1X2Err=%16.8e C2X2Err=%16.8e\n",
                           testname,i,fabs(distC1X1-r1),fabs(distC2X1-r2),fabs(distC1X2-r1),fabs(distC2X2-r2));
                    passedCount++;
                }
                else
                {
                   printf("%s-%-d C1(%14.8f,%14.8f), C2(%14.8f,%14.8f)\n",
                          testname,i ,nX,latC1,lonC1,latC2,lonC2);
                   printf("%s-%-d nX=%d  r1=%14.8f  r2=%14.8f  d12=%14.8f  usedROC=%14.8f\n",
                          testname,i ,nX,r1,r2,d12,usedROC);
                    printf("%s-%-d failed nX=%d  X1: %14.8f %14.8f  X2: %14.8f %14.8f\n",
                           testname,i,nX,latx1,lonx1,latx2,lonx2);
                    printf("%s-%-d failed C1X1=%14.8f C2X1=%14.8f C1X2=%14.8f C2X2=%14.8f X1X2=%14.8f\n",
                           testname,i,distC1X1,distC2X1,distC1X2,distC2X2,distX1X2);
                    printf("%s-%-d failed C1X1Err=%16.8e C2X1Err=%16.8e C1X2Err=%16.8e C2X2Err=%16.8e\n",
                           testname,i,fabs(distC1X1-r1),fabs(distC2X1-r2),fabs(distC1X2-r1),fabs(distC2X2-r2));
                    failedCount++;
                }
            }
            if (nX == 1)
            {
                distC1X1 = sphereInvDist(c1,intx[0],NULL);
                distC2X1 = sphereInvDist(c2,intx[0],NULL);
                if ( (fabs(distC1X1-r1)<=TESTTOL) &&
                        (fabs(distC2X1-r2)<=TESTTOL) )
                {
                	if (PRINT_PASSED_CASES) printf("%s-%-d nX=%d  r1=%14.8f  r2=%14.8f  d12=%14.8f  usedROC=%14.8f\n",
                           testname,i ,nX,r1,r2,d12,usedROC);
                	if (PRINT_PASSED_CASES) printf("%s-%-d passed nX=%d  X1: %14.8f %14.8f\n",
                           testname,i,nX,latx1,lonx1);
                	if (PRINT_PASSED_CASES) printf("%s-%-d passed C1X1=%14.8f C1X1Err=%16.8e C2X1=%14.8f C2X1Err=%16.8e\n",
                           testname,i,distC1X1,fabs(distC1X1-r1),distC2X1,fabs(distC2X1-r2));
                    passedCount++;
                }
                else
                {
                    printf("%s-%-d nX=%d  r1=%14.8f  r2=%14.8f  d12=%14.8f  usedROC=%14.8f\n",
                           testname,i ,nX,r1,r2,d12,usedROC);
                    printf("%s-%-d failed nX=%d  X1: %14.8f %14.8f\n",
                           testname,i,nX,latx1,lonx1);
                    printf("%s-%-d failed C1X1=%14.8f C1X1Err=%16.8e C2X1=%14.8f C2X1Err=%16.8e",
                           testname,i,distC1X1,fabs(distC1X1-r1),distC2X1,fabs(distC2X1-r2));
                    displayMatlabPt(c1,"c1",0);
                    displayMatlabPt(c2,"c2",0);
                    printf("\ncrs12 = %20.15f\n",crs12*180.0/M_PI);
                    crs12 = sphereInvCrs(c1,intx[0],1.0e-12);
                    printf("crsC1X1 = %20.15f; distC1X1 = %20.15f\n", crs12*180.0/M_PI, distC1X1);
                    crs12 = sphereInvCrs(c2,intx[0],1.0e-12);
                    printf("crsC2X1 = %20.15f; distC2X1 = %20.15f\n", crs12*180.0/M_PI, distC2X1);
                    printf("--------------------------------------------------------\n");
                    failedCount++;
                }
            }
            if (nX == 0)
            {
                if ( (d12>(r1+r2+TESTTOL)) ||
                        (d12<(fabs(r1-r2)-TESTTOL)) )
                {
                	if (PRINT_PASSED_CASES) printf("%s-%-d nX=%d  r1=%14.8f  r2=%14.8f  d12=%14.8f  usedROC=%14.8f\n",
                           testname,i ,nX,r1,r2,d12,usedROC);
                	if (PRINT_PASSED_CASES) printf("%s-%-d passed nX=%d\n",
                           testname,i,nX);
                    passedCount++;
                }
                else
                {
                    printf("%s-%-d nX=%d  r1=%14.8f  r2=%14.8f  d12=%14.8f  usedROC=%14.8f\n",
                           testname,i ,nX,r1,r2,d12,usedROC);
                    printf("%s-%-d failed nX=%d\n",
                           testname,i,nX);
                    failedCount++;
                }
            }

            testCaseCount++;
        } //for j
    }

    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("Finish testInitArcIntx_Set1\n\n\n");

    return set;

}

/*
 * NAME: testInitArcIntx_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the initArcIntx function.
 *
 * 		This function runs all the test cases from set1 and set2 and returns
 * 		the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testInitArcIntx_AllSets(TestSuite) - A test suite with the folling metrics:
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
TestSuite testInitArcIntx_AllSets()
{
	TestSuite suite;
	TestSet set1;

    printf("\nStart testInitArcIntx_AllSets\n");

    suite = newTestSuite("testInitArcIntx_AllSets");

    set1 = testInitArcIntx_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("Finish testInitArcIntx_AllSets\n\n\n");

    return suite;
}

/* NAME: testPtIsOnArc_Set1
 *
 * DESCRIPTION:
 *
 * 		This function runs the test data created by Balakrishna Babu and uses
 * 		his original test plan to execute the actual tests.
 *
 *      The approach for testing ptIsOnArc is as follows:
 *
 *      1. Select a center point, a radius, a start course and an end course.
 *      2. Compute the start and end points using direct.
 *      3. Compute test points along the Arc using azimuths of
 *         [0,45,90,135,180,225,270,315,azStart-AZDEGTOL,azStart,azStart+AZDEGTOL,
 *         0.5*(azStart+azEnd),azEnd-AZDEGTOL,azEnd,azEnd+AZDEGTOL].
 *      4. Compute test points at the above azimuths with Arc radii of [0.5*r, r-NMTOL, r+NMTOL, 1.1*r].
 *      5. For test points along Arc with azimuths between azStart and azEnd the result should be 1,
 *         and and 0 (not on Arc) for the rest.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testPtIsOnArc_Set1(TestSet) - A test set with the folling metrics:
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
TestSet testPtIsOnArc_Set1()
{
    //LLPoint center, double radius, double startCrs,
    //double endCrs, int orientation, LLPoint testPt,
    //ErrorSet* err, double tol, double eps)

    double DEG2RAD = M_PI / 180.0;
    double NMTOL = 0.03 / 100.0 / 1852.0; //0.03 cm or ~1.62e-7 nm
//    double LLDEGTOL = 5.0e-8; //5e-8 degrees or 0.00018 sec
    double AZDEGTOL = 5.0e-7; //5e-7 deg or 8.72665e-9 rad or 0.162 cm at 100 nm
    int testCaseCount = 0, setupFailureCount = 0, errorCount = 0, unverifiedCount = 0;
    char testname[80];
    double latStart, lonStart, latEnd, lonEnd, latTest, lonTest;
    double latCenter, lonCenter;
    double az1, az2, r1;
    double
            az0[15] = { 0.0, 45.0, 90.0, 135.0, 180.0, 225.0, 270.0, 315.0,
                        -AZDEGTOL, 0.0, AZDEGTOL, 0.5, -AZDEGTOL, 0.0, AZDEGTOL };
    double crs[15];
//    double r0[5] = { 0.5, -NMTOL, 0.0, NMTOL, 1.1 };
    double r[5];
    ArcDirection orient; //-1=counter clockwise, +1=clockwise
    int onArcExp[2][15][5] = { { { 0 } } }; //Initialize to all 0's
    int onArc = 0;
    int i, j;
    int testNum = 0;

    LLPoint center, start, end, test[5];
    double crsCS, crsCE;
    int passedCount = 0, failedCount = 0;
    ErrorSet err = 0;

    TestSet set;

    long newSeed = 20080111;

    printf("Start testPtIsOnArc_Set1\n");

    set = newTestSet("testPtIsOnArc_Set1");

    srand(newSeed); //Initialize the random number generator

    while (testNum < 100)
    {
        //Select center point, radius, start and end azimuths
        latCenter = randLat();
        lonCenter = randLon();
        az1 = randAzimuth();
        az2 = randAzimuth();
        r1 = 0.5 * randDist(); //about 2700 nm
        testNum++;

        sprintf(testname, "TEST%-d", testCaseCount);

        center.latitude = latCenter * DEG2RAD;
        center.longitude = lonCenter * DEG2RAD;
        crsCS = az1 * DEG2RAD;
        crsCE = az2 * DEG2RAD;

        //Compute start and end points
        err |= direct(center, crsCS, r1, &start, EPS);
        err |= direct(center, crsCE, r1, &end, EPS);

        latStart = start.latitude / DEG2RAD;
        lonStart = start.longitude / DEG2RAD;
        latEnd = end.latitude / DEG2RAD;
        lonEnd = end.longitude / DEG2RAD;

        //Set azimuths along Arc
        for (i = 0; i < 8; i++)
        {
            crs[i] = az0[i] * DEG2RAD;
        }
        crs[8] = (az1 - AZDEGTOL) * DEG2RAD;
        crs[9] = az1 * DEG2RAD;
        crs[10] = (az1 + AZDEGTOL) * DEG2RAD;
        crs[11] = 0.5 * (az1 + az2) * DEG2RAD;
        crs[12] = (az2 - AZDEGTOL) * DEG2RAD;
        crs[13] = az2 * DEG2RAD;
        crs[14] = (az2 + AZDEGTOL) * DEG2RAD;
        for (i = 0; i < 15; i++)
        {
            crs[i] = fmod((crs[i] + M_2PI * 2.0), M_2PI);
        }

        //Set radii for each Arc
        r[0] = 0.5 * r1; //not on Arc outside TOL
        r[1] = r1 - NMTOL; //not on Arc outside TOL
        r[2] = r1; //on Arc within TOL
        r[3] = r1 + NMTOL; //not on Arc outside TOL
        r[4] = 1.1 * r1; //not on Arc outside TOL

        //For each point along Arc
        for (i = 0; i < 15; i++)
        {

            //For each Arc with different radii
            for (j = 0; j < 5; j++)
            {

                orient = COUNTERCLOCKWISE; //counter clockwise
                onArcExp[0][i][j] = 0;
                if (j == 2)
                {
                    if (crsCE < crsCS)
                    { //doesn't cross 0
                        if ((crs[i] >= crsCE - TOL) && (crs[i] <= crsCS + TOL))
                        {
                            onArcExp[0][i][j] = 1;
                        }
                    }
                    else
                    { //crosses 0
                        if (((crs[i] >= crsCE - TOL) && (crs[i] <= M_2PI))
                                || ((crs[i] >= 0.0) && (crs[i] <= crsCS + TOL)))
                        {
                            onArcExp[0][i][j] = 1;
                        }
                    }
                }

                err |= direct(center, crs[i], r[j], &test[j], EPS);
                latTest = test[j].latitude / DEG2RAD;
                lonTest = test[j].longitude / DEG2RAD;

                if (getMaskedError(err, getMaskAll()))
                {
                    printf(
                            "Error occurred in pre-ptIsOnArc err=0x%lx i=%i j=%i orient=%i\n",
                            err, i, j, orient);
                    testCaseCount++;
                    setupFailureCount++;
                    errorCount++;
                    continue;
                }

                onArc = ptIsOnArc(center, r1, crsCS, crsCE, orient,
                        test[j], &err, TOL, EPS);

                if (getMaskedError(err, getMaskAll()))
                {
                    printf(
                            "Error occurred in ptIsOnArc err=0x%lx i=%i j=%i orient=%d\n",
                            err, i, j, orient);
                    failedCount++;
                    errorCount++;
                    testCaseCount++;
                    continue;
                }

                if (onArc == onArcExp[0][i][j])
                {
                	if (PRINT_PASSED_CASES) printf(
                            "%s Center: %14.8f %14.8f  azStart=%14.8f  azEnd=%14.8f  r1=%14.8f\n",
                            testname, latCenter, lonCenter, az1, az2, r1);
                	if (PRINT_PASSED_CASES) printf("%s start: %14.8f %14.8f  End: %14.8f %14.8f\n", testname,
                            latStart, lonStart, latEnd, lonEnd);
                	if (PRINT_PASSED_CASES) printf(
                            "%s passed orient=%d  az1=%14.8f az2=%14.8f  i=%d j=%d crs=%14.8f r=%14.8f  onArc=%d onArcExp=%d\n",
                            testname, orient, az1, az2, i, j, crs[i] / DEG2RAD,
                            r[j], onArc, onArcExp[0][i][j]);
                    passedCount++;
                    testCaseCount++;
                }
                else
                {
                    printf(
                            "%s Center: %14.8f %14.8f  azStart=%14.8f  azEnd=%14.8f  r1=%14.8f\n",
                            testname, latCenter, lonCenter, az1, az2, r1);
                    printf("%s start: %14.8f %14.8f  End: %14.8f %14.8f\n", testname,
                            latStart, lonStart, latEnd, lonEnd);
                    printf(
                            "%s failed orient=%d  az1=%14.8f az2=%14.8f  i=%d j=%d crs=%14.8f r=%14.8f  onArc=%d onArcExp=%d\n",
                            testname, orient, az1, az2, i, j, crs[i] / DEG2RAD,
                            r[j], onArc, onArcExp[0][i][j]);
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

    printf("Finish testPtIsOnArc_Set1\n");

    return set;
}

/*
 * NAME: testPtIsOnArc_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the ptIsOnArc function.
 *
 * 		This function runs all the test cases from set1 and set2 and returns
 * 		the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testPtIsOnArc_AllSets(TestSuite) - A test suite with the folling metrics:
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
TestSuite testPtIsOnArc_AllSets()
{
	TestSuite suite;
	TestSet set1;

    printf("\nStart testPtIsOnArc_AllSets\n");

    suite = newTestSuite("testPtIsOnArc_AllSets");

    set1 = testPtIsOnArc_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("Finish testPtIsOnArc_AllSets\n\n\n");

    return suite;
}

/* NAME: testPtIsInsideArc_Set1
 *
 * DESCRIPTION:
 *      Runs tests to determine if a given point is inside the pie slice under an arc
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testPtIsInsideArc_Set1(TestSet) - A test set with the following metrics:
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
TestSet testPtIsInsideArc_Set1()
{

    double DEG2RAD = M_PI / 180.0;
    int passedCount = 0;
    int failedCount = 0;
    int errorCount = 0;
    int testCaseCount = 0;
    int setupFailureCount = 0;
    int unverifiedCount = 0;
    LLPoint pt;
    LLPoint c;
    LLPoint startPt, endPt;
    double startCrs, endCrs;
    double radius;
    double ptCrs, angle, ptRadius, distToPoint;
    ErrorSet err = 0;
    int i,j;
    int inside, insideexp;
    long newSeed = 20080523;
    ArcDirection dir = ArcDirection::CLOCKWISE;

    TestSet set;

    printf("\nStart testPtIsInsideArc_Set1\n");

    set = newTestSet("testPtIsInsideArc_Set1");

    srand(newSeed); //Initialize the random number generator

    for(i = 0; i < 110; i++){
    	err = 0;

    	// Create arc center
    	if((i >= 100) && (i <= 104)){
    		c.latitude = 90.0*DEG2RAD; //A few North Pole tests
    	} else if ((i >= 105) && (i <= 109)){
    		c.latitude = -90.0*DEG2RAD; //A few South Pole tests
    	} else {
    		c.latitude = (-90.0 + 180.0*rand()/RAND_MAX)*DEG2RAD;
    	}
    	c.longitude = (-180.0 + 360.0*rand()/RAND_MAX)*DEG2RAD;

    	// Choose random parameters for arc.
    	radius = 0.1 + 1000.0*rand()/RAND_MAX;
    	startCrs = (0 + 360.0*rand()/RAND_MAX)*DEG2RAD;
    	endCrs = (0 + 360.0*rand()/RAND_MAX)*DEG2RAD;
    	if(startCrs == endCrs) endCrs = startCrs + 1.0e-6;
    	angle = computeSubtendedAngle(startCrs,endCrs,COUNTERCLOCKWISE);

    	// Test random points both in and out of the arc.
    	for(j = 0; j < 22; j++){
    		err = 0;

    		// Choose course to pt (ptCrs)
    		if(j == 20){
    			ptCrs = startCrs; //Test on boundary
    		} else if(j == 21){
    			ptCrs = endCrs; //Test on boundary
    		} else {
    			// Random course between startCrs and endCrs not on boundary
    			ptCrs = fmod((startCrs + (1.0e-6 + (angle - 2.0e-6)*rand()/RAND_MAX)) + M_2PI,M_2PI);
    		}

    		// Choose random distance for pt (ptRadius) not zero less than radius
    	 	ptRadius = 1.1*TOL + (radius - 1.1*TOL)*rand()/RAND_MAX;

    	 	err |= direct(c,ptCrs,ptRadius,&pt,EPS);
    	 	if(err) printf("Error occurred in direct %#lx\n",err);

    	 	// Inside radius/Inside courses. Test should be inside
    	 	dir = COUNTERCLOCKWISE;
    	 	err = 0;
    	 	testCaseCount++;
    	 	inside = ptIsInsideArc(c,radius,startCrs,endCrs,dir,pt,&err,TOL,EPS);
        	if(err) printf("Error occurred in ptIsInsideArc: %#lx at Test:%d\n",err, testCaseCount);

    	 	if(inside){
    	 		passedCount++;
    	 	} else {
    	 		failedCount++;
    	 		printf("fail1 %d: {%.20lf %.20lf} %.20lf %.20lf %.20lf %.20lf %.20lf %d\n",
    	 				testCaseCount, c.latitude/DEG2RAD, c.longitude/DEG2RAD, ptRadius,radius,
    	 				ptCrs/DEG2RAD, startCrs/DEG2RAD, endCrs/DEG2RAD,dir);
    	 	}

    	 	// Inside radius/Outside courses. Test should be outside
    	 	dir = CLOCKWISE;
    	 	err = 0;
    	 	testCaseCount++;
    	 	inside = ptIsInsideArc(c,radius,startCrs,endCrs,dir,pt,&err,TOL,EPS);
        	if(err) printf("Error occurred in ptIsInsideArc: %#lx at Test:%d\n",err, testCaseCount);

    	 	if(!inside){
    	 		if(j == 20 || j == 21){
        	 		// On boundary pt should be inside regardless of dir.
    	 			failedCount++;
        	 		printf("fail2 %d: {%.20lf %.20lf} %.20lf %.20lf %.20lf %.20lf %.20lf %d\n",
        	 				testCaseCount, c.latitude/DEG2RAD, c.longitude/DEG2RAD, ptRadius,radius,
        	 				ptCrs/DEG2RAD, startCrs/DEG2RAD, endCrs/DEG2RAD,dir);
    	 		} else {
    	 			passedCount++;
    	 		}
    	 	} else {
    	 		if(j == 20 || j == 21){
        	 		// On boundary pt should be inside regardless of dir.
    	 			passedCount++;
    	 		} else {
    	 			failedCount++;
    	 			printf("fail2 %d: {%.20lf %.20lf} %.20lf %.20lf %.20lf %.20lf %.20lf %d\n",
    	 				testCaseCount, c.latitude/DEG2RAD, c.longitude/DEG2RAD, ptRadius,radius,
    	 				ptCrs/DEG2RAD, startCrs/DEG2RAD, endCrs/DEG2RAD,dir);
    	 		}
    	 	}

    	 	// Choose random distance for pt (ptRadius) larger than radius
    	 	ptRadius = radius + 1.1*TOL + 1000*rand()/RAND_MAX;

    	 	direct(c,ptCrs,ptRadius,&pt,EPS);
    	 	if(err) printf("Error occurred in direct %#lx\n",err);

    	 	// Outside radius/Inside courses. Test should be outside.
    	 	dir = COUNTERCLOCKWISE;
    	 	err = 0;
    	 	testCaseCount++;
    	 	inside = ptIsInsideArc(c,radius,startCrs,endCrs,dir,pt,&err,TOL,EPS);
        	if(err) printf("Error occurred in ptIsInsideArc: %#lx at Test:%d\n",err, testCaseCount);

    	 	if(!inside){
    	 		passedCount++;
    	 	} else {
    	 		failedCount++;
    	 		printf("fail3 %d: {%.20lf %.20lf} %.20lf %.20lf %.20lf %.20lf %.20lf %d\n",
    	 				testCaseCount, c.latitude/DEG2RAD, c.longitude/DEG2RAD, ptRadius,radius,
    	 				ptCrs/DEG2RAD, startCrs/DEG2RAD, endCrs/DEG2RAD,dir);
    	 	}

    	 	// Outside radius/Outside courses. Test should be outside
    	 	dir = CLOCKWISE;
    	 	err = 0;
    	 	testCaseCount++;
    	 	inside = ptIsInsideArc(c,radius,startCrs,endCrs,dir,pt,&err,TOL,EPS);
        	if(err) printf("Error occurred in ptIsInsideArc: %#lx at Test:%d\n",err, testCaseCount);

    	 	if(!inside){
    	 		passedCount++;
    	 	} else {
    	 		failedCount++;
    	 		printf("fail4 %d: {%.20lf %.20lf} %.20lf %.20lf %.20lf %.20lf %.20lf %d\n",
    	 				testCaseCount, c.latitude/DEG2RAD, c.longitude/DEG2RAD, ptRadius,radius,
    	 				ptCrs/DEG2RAD, startCrs/DEG2RAD, endCrs/DEG2RAD,dir);
    	 	}

    	}

    	// For each circle do one pt at the center. Test should be inside.
    	ptCrs = (0.0 + 360.0*rand()/RAND_MAX)*DEG2RAD;
	 	err |= direct(c,ptCrs,0.0,&pt,EPS);
	 	if(err) printf("Error occurred in direct %#lx\n",err);

	 	err = 0;
	 	testCaseCount++;
	 	inside = ptIsInsideArc(c,radius,startCrs,endCrs,dir,pt,&err,TOL,EPS);
    	if(err) printf("Error occurred in ptIsInsideArc: %#lx at Test:%d\n",err, testCaseCount);

	 	if(inside){
	 		passedCount++;
	 	} else {
	 		failedCount++;
	 		printf("fail5 %d: {%.20lf %.20lf} %.20lf %.20lf %.20lf %.20lf %.20lf %d\n",
	 				testCaseCount, c.latitude/DEG2RAD, c.longitude/DEG2RAD, 0.0 ,radius,
	 				ptCrs/DEG2RAD, startCrs/DEG2RAD, endCrs/DEG2RAD,dir);
	 	}

    }

    //Test for test pt at pole
    for(i = 0;  i < 5; i++){
    	err = 0;
        insideexp = 0;

    	pt.latitude = 90.0*DEG2RAD;
    	pt.longitude = 0.0;
    	ptCrs = (0.0 + 360.0*rand()/RAND_MAX)*DEG2RAD;
	 	ptRadius = 1.1*TOL + 1000*rand()/RAND_MAX;

	 	err |= direct(pt,ptCrs,ptRadius,&c,EPS);
                err |= direct(c,10.0*DEG2RAD,ptRadius/2.0,&startPt,EPS);
                err |= direct(c,350.0*DEG2RAD,ptRadius/2.0,&endPt,EPS);
                err |= invDist(pt,startPt,&distToPoint,EPS);
                if (distToPoint <= TOL)
                  insideexp = 1;
                err |= invDist(pt,endPt,&distToPoint,EPS);
                if (distToPoint <= TOL)
                  insideexp = 1;
	 	if(err) printf("Error occurred in direct %#lx\n",err);

	 	err = 0;
	 	testCaseCount++;
    	inside = ptIsInsideArc(c,ptRadius/2.0,10.0*DEG2RAD,350*DEG2RAD,COUNTERCLOCKWISE,pt,&err,TOL,EPS);
    	if(err) printf("Error occurred in ptIsInsideArc: %#lx at Test:%d\n",err, testCaseCount);

	 	if(inside == insideexp){
	 		passedCount++;
	 	} else {
	 		failedCount++;
	 		printf("fail6 %d: {%.20lf %.20lf} %.20lf %.20lf %.20lf %.20lf %.20lf %d\n",
	 				testCaseCount, c.latitude/DEG2RAD, c.longitude/DEG2RAD, ptRadius,radius,
	 				ptCrs/DEG2RAD, startCrs/DEG2RAD, endCrs/DEG2RAD,dir);
	 	}

	 	err = 0;
	 	testCaseCount++;
    	inside = ptIsInsideArc(c,2.0*ptRadius,10.0*DEG2RAD,350*DEG2RAD,COUNTERCLOCKWISE,pt,&err,TOL,EPS);
	 	if(err) printf("Error occurred in ptIsInsideArc: %#lx at Test:%d\n",err, testCaseCount);

	 	if(inside){
	 		passedCount++;
	 	} else {
	 		failedCount++;
	 		printf("fail6 %d: {%.20lf %.20lf} %.20lf %.20lf %.20lf %.20lf %.20lf %d\n",
	 				testCaseCount, c.latitude/DEG2RAD, c.longitude/DEG2RAD, ptRadius,radius,
	 				ptCrs/DEG2RAD, startCrs/DEG2RAD, endCrs/DEG2RAD,dir);
	 	}

    }


    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("\nFinish testPtIsInsideArc_Set1\n");

    return set;
}

/*
 * NAME: testPtIsInsideArc_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the ptIsInsideArc function.
 *
 * 		This function runs all the test cases from set1 and set2 and returns
 * 		the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testPtIsInsideArc_AllSets(TestSuite) - A test suite with the folling metrics:
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
TestSuite testPtIsInsideArc_AllSets()
{
	TestSuite suite;
	TestSet set1;

    printf("\nStart testPtIsInsideArc_AllSets\n");

    suite = newTestSuite("testPtIsInsideArc_AllSets");

    set1 = testPtIsInsideArc_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("Finish testPtIsInsideArc_AllSets\n\n\n");

    return suite;
}

/*
 * NAME: testArcIntx_Set1
 *
 * DESCRIPTION:
 *      Tests arcIntx by creating two arcs that should have two intersections.
 *
 * 		This function runs the test data created by Balakrishna Babu and uses
 * 		his original test plan to execute the actual tests.
 *
 *      The approach for testing arcIntx is as follows:
 *
 *   1. Start with one expected intersection point.
 *   2. Select an azimuth and radius to compute center point of first Arc using direct.
 *   3. Select an azimuth and radius to compute center point of second Arc using direct.
 *   4. Two center points and their radii form the input.
 *   5. One of the resulting two intersection points should match the expected intersection point.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testArcIntx_Set1(TestSet) - A test set with the folling metrics:
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
 */
TestSet testArcIntx_Set1()
{
    double DEG2RAD = M_PI / 180.0;
    char testname[80];
    double latxExp, lonxExp;
    double az1, r1, az2, r2;
    LLPoint intxExp;

    LLPoint c1, c2;
    LLPointPair intx;
    int nX;
    double latx1, lonx1, latx2, lonx2;
    double crs12, crs21, dist12;
    double crs13, crs31, dist13;
    double crs23, crs32, dist23;
    int passedCount = 0, failedCount = 0;
    int testCaseCount = 0, setupFailureCount = 0, errorCount = 0, unverifiedCount = 0;
    ErrorSet err = 0;

    TestSet set;

    long newSeed = 20080101;

    printf("Start testArcIntx_Set1\n");

    set = newTestSet("testArcIntx_Set1");

    srand(newSeed); //Initialize the random number generator

    while (testCaseCount < 1000)
    {
        err = 0;
        latxExp = randLat(); //intersection latitude, deg
        lonxExp = randLon(); //intersection longitude, deg
        az1 = randAzimuth(); //0-360, course to centerPoint1
        r1 = 0.2 * randDist(); //about 1080 nm
        az2 = randAzimuth(); //0-360, course to centerPoint2
        r2 = 0.2 * randDist();

        sprintf(testname, "TEST%-d", testCaseCount);

        intxExp.latitude = latxExp * DEG2RAD;
        intxExp.longitude = lonxExp * DEG2RAD;

        crs31 = az1 * DEG2RAD; //course intersection to centerPoint1
        crs32 = az2 * DEG2RAD; //course intersection to centerPoint2

        err |= direct(intxExp, crs31, r1, &c1, EPS); //populates centerPoint1
        err |= direct(intxExp, crs32, r2, &c2, EPS); //populates centerPoint2

        //populates courses between intersection and centerPoint1, and distance
        err |= inverse(c1, intxExp, &crs13, &crs31, &dist13, EPS);
        //populates courses between intersection and centerPoint2, and distance
        err |= inverse(c2, intxExp, &crs23, &crs32, &dist23, EPS);
        //populates courses between two centerPoints, and the distance
        err |= inverse(c1, c2, &crs12, &crs21, &dist12, EPS);

        if (getMaskedError(err, getMaskAll()))
        {
            printf("Error occurred in pre-arcIntx err=0x%lx\n", err);
            testCaseCount++;
            setupFailureCount++;
            errorCount++;
            continue;
        }

        //updates intx (array of LLPoint structs), nX with # intersections
        err |= arcIntx(c1, r1, c2, r2, intx, &nX, TOL, EPS);

        if (getMaskedError(err, getMaskAll()))
        {
            printf("Error occurred in arcIntx err=0x%lx\n", err);
            failedCount++;
            errorCount++;
            testCaseCount++;
            continue;
        }

        latx1 = intx[0].latitude / DEG2RAD;
        lonx1 = intx[0].longitude / DEG2RAD;
        latx2 = intx[1].latitude / DEG2RAD;
        lonx2 = intx[1].longitude / DEG2RAD;

        if ( (nX > 0) && (
                     ptsAreSame(intx[0], intxExp, TESTTOL) || ptsAreSame(intx[1], intxExp, TESTTOL)) )
        {

        	if (PRINT_PASSED_CASES) printf("%s,%14.8f,%14.8f,", testname, latxExp, lonxExp);
            //print the center points
        	if (PRINT_PASSED_CASES) printf("%14.8f,%14.8f,%14.8f,%14.8f,", c1.latitude / DEG2RAD,
                    c1.longitude / DEG2RAD, c2.latitude / DEG2RAD, c2.longitude
                            / DEG2RAD);
            //print the radii
        	if (PRINT_PASSED_CASES) printf("%14.8f,%14.8f,", r1, r2);
            //print the courses
        	if (PRINT_PASSED_CASES) printf("%14.8f,%14.8f,%14.8f,%14.8f,%14.8f,%14.8f,", crs13, crs31,
                    crs23, crs32, crs12, crs21);

            //print true for passed, print the two intersection results
        	if (PRINT_PASSED_CASES) printf("true,%14.8f,%14.8f,%14.8f,%14.8f\n", latx1, lonx1, latx2,
                    lonx2);
            passedCount++;
        }
        else
        {

        	printf("%s,%14.8f,%14.8f,", testname, latxExp, lonxExp);
            //print the center points
            printf("%14.8f,%14.8f,%14.8f,%14.8f,", c1.latitude / DEG2RAD,
                    c1.longitude / DEG2RAD, c2.latitude / DEG2RAD, c2.longitude
                            / DEG2RAD);
            //print the radii
            printf("%14.8f,%14.8f,", r1, r2);
            //print the courses
            printf("%14.8f,%14.8f,%14.8f,%14.8f,%14.8f,%14.8f,", crs13, crs31,
                    crs23, crs32, crs12, crs21);

            //print false for passed, print N/A for all intersection values
            //printf("false,N/A,N/A,N/A,N/A\n");
            printf("false,%14.8f,%14.8f,%14.8f,%14.8f\n", latx1, lonx1, latx2, lonx2);
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

    printf("\nFinish testArcIntx_Set1\n");

    return set;

}

/*
 * NAME:
 * 		testArcIntx_Set2
 *
 * DESCRIPTION:
 *      Tests arcIntx by creating two arcs that should have no intersections.
 *
 * 		This function runs the test data created by Balakrishna Babu and uses
 * 		his original test plan to execute the actual tests.
 *
 * 		This method was originally named testArcIntxZero.
 *
 *      The approach for testing arcIntx is as follows:
 *
 *      1. Start with one expected intersection point.
 *      2. Select an azimuth and radius to compute center point of first Arc using direct.
 *      3. Select an azimuth and radius to compute center point of second Arc using direct.
 *      4. Compute the distance between the two center points using inverse.
 *      5. Set the radius of the second Arc as the difference between the radius of the first Arc and
 *          the distance between the two centers minus the distance tolerance.  So, zero intersection points
 *          are expected.
 *      6. Two center points and their radii form the input.
 *      7. Result should be zero intersection points.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testArcIntx_Set1(TestSet) - A test set with the folling metrics:
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
TestSet testArcIntx_Set2()
{
    //LLPoint center1, double r1, LLPoint center2, double r2,
    //LLPointPair intx, int* n, double tol, double eps);
    double DEG2RAD = M_PI / 180.0;
//    double RAD2DEG = 180.0 / M_PI;
    double NMTOL = 0.03 / 100.0 / 1852.0;  //0.03 cm or ~1.62e-7 nm
    char testname[80];
    double latxExp, lonxExp;
    double az1, r1, az2, r2;
    LLPoint intxExp;

    LLPoint c1, c2;
    LLPointPair intx;
    int nX;
    double latx1, lonx1, latx2, lonx2;
    double crs12, crs21, dist12;
    double crs31;
    double crs32;
    int passedCount = 0, failedCount = 0;
    int testCaseCount = 0, setupFailureCount = 0, errorCount = 0, unverifiedCount = 0;
    ErrorSet err = 0;

    TestSet set;

    long newSeed = 20080101;

    printf("Start testArcIntx_Set2\n");

    set = newTestSet("testArcIntx_Set2");

    srand(newSeed); //Initialize the random number generator

    while (testCaseCount < 1000)
    {
    	err = 0;
        latxExp = randLat();
        lonxExp = randLon();
        az1 = randAzimuth();
        r1 = 0.2 * randDist(); //about 1080 nm
        az2 = randAzimuth();
        r2 = 0.2 * randDist();

        sprintf(testname, "TEST%-d", testCaseCount);

        intxExp.latitude = latxExp * DEG2RAD;
        intxExp.longitude = lonxExp * DEG2RAD;
        crs31 = az1 * DEG2RAD;
        crs32 = az2 * DEG2RAD;

        err |= direct(intxExp, crs31, r1, &c1, EPS);
        err |= direct(intxExp, crs32, r2, &c2, EPS);

        err |= inverse(c1, c2, &crs12, &crs21, &dist12, EPS);
        err |= direct(c1, crs12, r1, &intxExp, EPS);

        r2 = fabs(r1 - dist12) - NMTOL; //0.03 cm
        latxExp = intxExp.latitude / DEG2RAD;
        lonxExp = intxExp.longitude / DEG2RAD;

        if (getMaskedError(err, getMaskAll()))
        {
            printf("Error occurred in pre-arcIntx Zero err=0x%lx\n",
                    err);
            testCaseCount++;
            setupFailureCount++;
            errorCount++;
            continue;
        }

        err |= arcIntx(c1, r1, c2, r2, intx, &nX, TOL, EPS);

        if (getMaskedError(err, getMaskAll()))
        {
            printf("Error occurred in arcIntx Zero err=0x%lx\n", err);
            failedCount++;
            errorCount++;
            testCaseCount++;
            continue;
        }

        latx1 = intx[0].latitude / DEG2RAD;
        lonx1 = intx[0].longitude / DEG2RAD;
        latx2 = intx[1].latitude / DEG2RAD;
        lonx2 = intx[1].longitude / DEG2RAD;

        if (nX == 0)
        {
        	if (PRINT_PASSED_CASES) printf(
                    "%s Intx: %14.8f %14.8f  az1=%14.8f r1=%14.8f  az2=%14.8f r2=%14.8f\n",
                    testname, latxExp, lonxExp, az1, r1, az2, r2);
        	if (PRINT_PASSED_CASES) printf("%s c1:   %14.8f %14.8f   r1: %14.8f\n", testname, c1.latitude
                    / DEG2RAD, c1.longitude / DEG2RAD, r1);

        	if (PRINT_PASSED_CASES) printf("%s c2:   %14.8f %14.8f   r2: %14.8f\n", testname, c2.latitude
                    / DEG2RAD, c2.longitude / DEG2RAD, r2);

        	if (PRINT_PASSED_CASES) printf("%s intx: %14.8f %14.8f   dist12: %14.8f\n", testname,
                    intxExp.latitude / DEG2RAD, intxExp.longitude / DEG2RAD, dist12);
        	if (PRINT_PASSED_CASES) printf("%s passed nX: %d\n", testname, nX);
            passedCount++;
        }
        else
        {
            printf(
                    "%s Intx: %14.8f %14.8f  az1=%14.8f r1=%14.8f  az2=%14.8f r2=%14.8f\n",
                    testname, latxExp, lonxExp, az1, r1, az2, r2);
            printf("%s c1:   %14.8f %14.8f   r1: %14.8f\n", testname, c1.latitude
                    / DEG2RAD, c1.longitude / DEG2RAD, r1);

            printf("%s c2:   %14.8f %14.8f   r2: %14.8f\n", testname, c2.latitude
                    / DEG2RAD, c2.longitude / DEG2RAD, r2);

            printf("%s intx: %14.8f %14.8f   dist12: %14.8f\n", testname,
                    intxExp.latitude / DEG2RAD, intxExp.longitude / DEG2RAD, dist12);
            printf(
                    "%s failed nX: %d  Exp: %14.8f %14.8f   intx1: %14.8f %14.8f  intx2: %14.8f %14.8f\n",
                    testname, nX, latxExp, lonxExp, latx1, lonx1, latx2, lonx2);
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

    printf("Finish testArcIntx_Set2\n");

    return set;

}
/*
 * NAME: testArcIntx_Set3
 *
 * DESCRIPTION:
 *      Tests arcIntx by creating two arcs that should have one intersection.
 *
 * 		This function runs the test data created by Balakrishna Babu and uses
 * 		his original test plan to execute the actual tests.
 *
 * 		This method was originally named testArcIntxOne.
 *
 *   1. Start with one expected intersection point.
 *   2. Select an azimuth and radius to compute center point of first Arc using direct.
 *   3. Select an azimuth and radius to compute center point of second Arc using direct.
 *   4. Compute the distance between the two center points using inverse.
 *   5. Set the radius of the second Arc as the difference between the radius of the first Arc and
 *      the distance between the two centers.  So, the two Arcs should touch and not intersect.
 *   6. Compute the course from center of first Arc to center of second Arc.
 *   7. Using this course, center and radius of first Arc compute the touching point using inverse.
 *   8. Two center points and their radii form the input.
 *   9. Result should be one intersection point, which is the point where two Arcs touch.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testArcIntx_Set3(TestSet) - A test set with the folling metrics:
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
TestSet testArcIntx_Set3()
{
    double DEG2RAD = M_PI / 180.0;
    int testCaseCount = 0, setupFailureCount = 0, errorCount = 0, unverifiedCount = 0;
    char testname[80];
    double latxExp, lonxExp;
    double az1, r1, az2, r2;
    LLPoint intxExp;

    LLPoint c1, c2;
    LLPointPair intx;
    int nX;
    double latx1, lonx1, latx2, lonx2;
    double crs12, crs21, dist12;
    double crs31;
    double crs32;
    int passedCount = 0, failedCount = 0;
    ErrorSet err = 0;

    TestSet set;

    long newSeed = 20080101;

    printf("Start testArcIntx_Set3\n");

    set = newTestSet("testArcIntx_Set3");

    srand(newSeed); //Initialize the random number generator

    while (testCaseCount < 1000)
    {
    	err = 0;
        latxExp = randLat();
        lonxExp = randLon();
        az1 = randAzimuth();
        r1 = 0.2 * randDist(); //about 1080 nm
        az2 = randAzimuth();
        r2 = 0.2 * randDist();

        sprintf(testname, "TEST%-d", testCaseCount);

        intxExp.latitude = latxExp * DEG2RAD;
        intxExp.longitude = lonxExp * DEG2RAD;
        crs31 = az1 * DEG2RAD;
        crs32 = az2 * DEG2RAD;

        err |= direct(intxExp, crs31, r1, &c1, EPS);
        err |= direct(intxExp, crs32, r2, &c2, EPS);

        err |= inverse(c1, c2, &crs12, &crs21, &dist12, EPS);
        err |= direct(c1, crs12, r1, &intxExp, EPS);

        r2 = fabs(r1 - dist12);
        latxExp = intxExp.latitude / DEG2RAD;
        lonxExp = intxExp.longitude / DEG2RAD;

        if (getMaskedError(err, getMaskAll()))
        {
            printf("Error occurred in pre-arcIntx One err=%s",
            		formatErrorMessage(err));
            testCaseCount++;
            setupFailureCount++;
            errorCount++;
            continue;
        }

        err |= arcIntx(c1, r1, c2, r2, intx, &nX, TOL, EPS);

        if (getMaskedError(err, getMaskAll()))
        {
            printf("Error occurred in arcIntx One err=0x%lx\n", err);
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
        	if (PRINT_PASSED_CASES) printf(
                    "%s Intx: %14.8f %14.8f  az1=%14.8f r1=%14.8f  az2=%14.8f r2=%14.8f   ",
                    testname, latxExp, lonxExp, az1, r1, az2, r2);
        	if (PRINT_PASSED_CASES) printf("%s c1:   %14.8f %14.8f   r1: %14.8f    ", testname, c1.latitude
                    / DEG2RAD, c1.longitude / DEG2RAD, r1);

        	if (PRINT_PASSED_CASES) printf("%s c2:   %14.8f %14.8f   r2: %14.8f    ", testname, c2.latitude
                    / DEG2RAD, c2.longitude / DEG2RAD, r2);

        	if (PRINT_PASSED_CASES) printf("%s intx: %14.8f %14.8f   dist12: %14.8f    ", testname,
                    intxExp.latitude / DEG2RAD, intxExp.longitude / DEG2RAD, dist12);
        	if (PRINT_PASSED_CASES) printf(
                    "%s passed nX: %d  latx1: %14.8f %14.8e    lonx1: %14.8f %14.8e\n",
                    testname, nX, latx1, fabs(latxExp - latx1), lonx1, fabs(
                            lonxExp - lonx1));
            passedCount++;
        }
        else
        {
            printf(
                    "%s Intx: %14.8f %14.8f  az1=%14.8f r1=%14.8f  az2=%14.8f r2=%14.8f   ",
                    testname, latxExp, lonxExp, az1, r1, az2, r2);
            printf("%s c1:   %14.8f %14.8f   r1: %14.8f   ", testname, c1.latitude
                    / DEG2RAD, c1.longitude / DEG2RAD, r1);

            printf("%s c2:   %14.8f %14.8f   r2: %14.8f   ", testname, c2.latitude
                    / DEG2RAD, c2.longitude / DEG2RAD, r2);

            printf("%s intx: %14.8f %14.8f   dist12: %14.8f   ", testname,
                    intxExp.latitude / DEG2RAD, intxExp.longitude / DEG2RAD, dist12);
            printf(
                    "%s failed nX: %d  Exp: %14.8f %14.8f   intx1: %14.8f %14.8f  intx2: %14.8f %14.8f\n",
                    testname, nX, latxExp, lonxExp, latx1, lonx1, latx2, lonx2);
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


    printf("Finish testArcIntx_Set3\n");

    return set;

}

/*
 * NAME: testArcIntx_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the ArcIntersect function.
 *
 * 		This function runs all the test cases and returns
 * 		the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testArcIntx_AllSets(TestSuite) - A test suite with the folling metrics:
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
TestSuite testArcIntx_AllSets()
{
	TestSuite suite;
	TestSet set1;
	TestSet set2;
	TestSet set3;


    printf("\nStart testArcIntx_AllSets\n");

    suite = newTestSuite("testArcIntx_AllSets");

    set1 = testArcIntx_Set1();
    addTestSet(set1,&suite);

    set2 = testArcIntx_Set2();
    addTestSet(set2,&suite);

    set3 = testArcIntx_Set3();
    addTestSet(set3,&suite);

    displayTestSuite(suite);

    printf("Finish testArcIntx_AllSets\n\n\n");

    return suite;
}

/*
 * NAME: testArcLength_Set1
 *
 * 		This function is used to test the arcLength function.
 *
 * DESCRIPTION:
 *      The approach for testing arcLength is as follows:
 *
 *      1. Select a center, a radius, start course and end course.
 *      2. Using the selected data as input, compute the discretized Arc
 *         length with orient = -1, and +1.
 *      3. Divide the azimuth between start and courses into 100 and compute
 *         points along the Arc at those azimuths using direct.
 *      4. Compute the distances between the consecutive 101 points along the
 *         Arc using inverse.
 *      5. Sum of the distances should approximately equal the Arc length.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testArcLength_Set1(TestSet) - A test set with the folling metrics:
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
TestSet testArcLength_Set1()
{
    double DEG2RAD = M_PI / 180.0;
    double RAD2DEG = 180.0 / M_PI;
    char testname[80];
    double latCenter, lonCenter;
    double az1, az2, r1;
    ArcDirection orient; //-1=counter clockwise, 1=clockwise
    int steps, i;

    LLPoint center, start, end, pt, ptNext;
    double crsCS, crsCE, deltacrs, crs, crsIJ, crsJI;
    double arcLenExp, arcLen;
    double distI, arcLenApprx, meanR;
    int passedCount = 0, failedCount = 0;
    int testCaseCount = 0, setupFailureCount = 0, errorCount = 0, unverifiedCount = 0;
    ErrorSet err = 0;
    long stepsForExp = 100000;
    long seed = 20080115;
    int testCase;

    TestSet set;

    srand(seed);
    printf("Start testArcLength_Set1\n");

    set = newTestSet("testArcLength_Set1");

    while(testCaseCount < 100)
    {
	err = 0;
	steps = 0;
        testCase = testCaseCount;
        //Select center point, radius, start and end azimuths
        latCenter = randLat();
        lonCenter = randLon();
        az1 = randAzimuth();
        az2 = randAzimuth();
        r1 = 0.01 * randDist(); //about 54 nm

        sprintf(testname, "TEST%-d", testCase);

        center.latitude = latCenter * DEG2RAD;
        center.longitude = lonCenter * DEG2RAD;
        crsCS = az1 * DEG2RAD;
        crsCE = az2 * DEG2RAD;

        //Compute start and end points
        err |= direct(center, crsCS, r1, &start, EPS);
        err |= direct(center, crsCE, r1, &end, EPS);

        orient = CLOCKWISE; //clockwise
        arcLenApprx = approxArcLength(center, r1, crsCS, crsCE, orient, &meanR);
        stepsForExp = (int) ceil((1.0 + arcLenApprx) * 100.0);
        if (r1 < 3.0) //more steps needed for small arcs
          stepsForExp = 2000;

        //Compute delta course along clockwise Arc
        if (crsCE < crsCS)
        {
            deltacrs = (M_2PI + crsCE - crsCS) / stepsForExp;
        }
        else
        {
            if (fabs(crsCE - crsCS) > 0.0)
              deltacrs = (crsCE - crsCS) / stepsForExp;
            else
              deltacrs = M_2PI / stepsForExp;
        }

        //Compute points along clockwise Arc
        arcLenExp = 0.0;
        pt = start;
        crs = crsCS;
        for (i = 0; i < stepsForExp; i++)
        {
            crs += deltacrs;
            if (crs >= M_2PI)
            {
                crs -= M_2PI;
            }
            err |= direct(center, crs, r1, &ptNext, EPS);
            err |= inverse(pt, ptNext, &crsIJ, &crsJI, &distI, EPS);
            arcLenExp += distI;
            pt = ptNext;
        }// for i

	if (getMaskedError(err, getMaskAll()))
	{
	     printf("Error occurred in pre-arcLength err=0x%lx\n",err);
	     testCaseCount++;
	     setupFailureCount++;
	     errorCount++;
	     continue;
	}

	arcLen = arcLength(center, r1, crsCS, crsCE, orient, &steps, &err, TOL, EPS);

	if (getMaskedError(err, getMaskAll()))
	{
	    printf("Error occurred in arcLength err=0x%lx\n", err);
	    failedCount++;
	    errorCount++;
	    testCaseCount++;
	    continue;
	}

	if (fabs(arcLen - arcLenExp) < TESTTOL * 100.0)
	{
	        	//if (PRINT_PASSED_CASES) printf("\nPass,%s",line);
            passedCount++;
	}
	else
	{
            printf("\nFail,%s latC %6.3f lonC %7.3f rad %7.3f crsCS %6.3f crsCE %6.3f arcLen %14.8f arcLenExp %14.8f",testname,latCenter,lonCenter,r1,crsCS*RAD2DEG,crsCE*RAD2DEG,arcLen,arcLenExp);
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

    printf("Finish testArcLength_Set1\n");

    return set;
}

/*
 * NAME: testArcLength_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the arcLength function.
 *
 * 		This function runs all the test cases and returns the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testArcLength_AllSets(TestSuite) - A test suite with the folling metrics:
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
TestSuite testArcLength_AllSets()
{
	TestSuite suite;
	TestSet set1;

    printf("\nStart testArcLength_AllSets\n");

    suite = newTestSuite("testArcLength_AllSets");

    set1 = testArcLength_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("Finish testArcLength_AllSets\n\n\n");

    return suite;
}

/*
 * NAME: testArcFromLength_Set1
 *
 * DESCRIPTION:
 *      The approach for testing arcFromLength is as follows:
 *
 *      1. Select a center, a radius, start course and end course.
 *      2. Using direct compute start and end points of the Arc.
 *      3. Using the selected data as input, compute the discretized Arc
 *         length with orient = -1, and +1.
 *      4. Now, using the descretized Arc length, center, start point, and orinet
 *         as input compute the end point.
 *      5. The computed end point should match the end point from step 2.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testArcFromLength_Set1(TestSet) - A test set with the folling metrics:
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
TestSet testArcFromLength_Set1()
{
    //LLPoint center, LLPoint givenPoint, ArcDirection dir,
    //double arcDist, LLPoint* endPoint, double *subtendedAngle,
    //double tol, double eps)

    double DEG2RAD = M_PI / 180.0;
    int testCaseCount = 0, setupFailureCount = 0, errorCount = 0, unverifiedCount = 0;
    char testname[80];
    double latStart, lonStart, latEnd, lonEnd;
    double latCenter, lonCenter;
    double az1, az2, r1;
    int steps;

    LLPoint center, start, endExp, endClock, endCC;
    double crsCS, crsCE;
    double subtendedClock, subtendedCC;
    double arcLenClock, arcLenCC;
    int passedCount = 0, failedCount = 0;
    ErrorSet err = 0;

    TestSet set;

    long newSeed = 20080208;

    printf("Start testArcFromLength_Set1\n");

    set = newTestSet("testArcFromLength_Set1");

    srand(newSeed); //Initialize the random number generator

    while (testCaseCount < 1000)
    {
        //Select center point, radius, start and end azimuths
        latCenter = randLat();
        lonCenter = randLon();
        az1 = randAzimuth();
        az2 = randAzimuth();
        r1 = 0.2 * randDist(); //about 1080 nm

        sprintf(testname, "TEST%-d", testCaseCount);

        center.latitude = latCenter * DEG2RAD;
        center.longitude = lonCenter * DEG2RAD;
        crsCS = az1 * DEG2RAD;
        crsCE = az2 * DEG2RAD;

        //Compute start and end points
        err |= direct(center, crsCS, r1, &start, EPS);
        err |= direct(center, crsCE, r1, &endExp, EPS);

        latStart = start.latitude / DEG2RAD;
        lonStart = start.longitude / DEG2RAD;
        latEnd = endExp.latitude / DEG2RAD;
        lonEnd = endExp.longitude / DEG2RAD;

        steps = 0;
        arcLenClock = arcLength(center, r1, crsCS, crsCE,
                CLOCKWISE, &steps, &err, TOL, EPS);

        steps = 0;
        arcLenCC = arcLength(center, r1, crsCS, crsCE,
                COUNTERCLOCKWISE, &steps, &err, TOL, EPS);

        if (getMaskedError(err, getMaskAll()))
        {
            printf("Error occurred in pre-arcFromLength err=0x%lx\n",
                    err);\

            testCaseCount++;
            setupFailureCount++;
            errorCount++;
            continue;
        }

        err |= arcFromLength(center, start, CLOCKWISE, arcLenClock,
                &endClock, &subtendedClock, TOL, EPS);
        err |= arcFromLength(center, start, COUNTERCLOCKWISE,
                arcLenCC, &endCC, &subtendedCC, TOL, EPS);

        if (getMaskedError(err, getMaskAll()))
        {
            printf("Error occurred in arcFromLength err=0x%lx\n", err);\

            failedCount++;
            errorCount++;
            testCaseCount++;
            continue;
        }

        if (ptsAreSame(endExp, endClock, TESTTOL) && ptsAreSame(endExp, endCC, TESTTOL))
        {
        	if (PRINT_PASSED_CASES) printf(
                    "%s Center: %14.8f %14.8f  azStart=%14.8f  azEnd=%14.8f  r1=%14.8f\n",
                    testname, latCenter, lonCenter, az1, az2, r1);
        	if (PRINT_PASSED_CASES) printf(
                    "%s start: %14.8f %14.8f  End: %14.8f %14.8f  ArcLenClock: %14.8f  ArcLenCC: %14.8f\n",
                    testname, latStart, lonStart, latEnd, lonEnd, arcLenClock,
                    arcLenCC);
        	if (PRINT_PASSED_CASES) printf(
                    "%s passed EndClock: %14.8f %14.8f  EndCC: %14.8f %14.8f  subtendedClock: %14.8f  subtendedCC: %14.8f\n",
                    testname, endClock.latitude / DEG2RAD, endClock.longitude
                            / DEG2RAD, endCC.latitude / DEG2RAD,
                    endCC.longitude / DEG2RAD, subtendedClock / DEG2RAD,
                    subtendedCC / DEG2RAD);
            passedCount++;
        }
        else
        {
            printf(
                    "%s Center: %14.8f %14.8f  azStart=%14.8f  azEnd=%14.8f  r1=%14.8f\n",
                    testname, latCenter, lonCenter, az1, az2, r1);
            printf(
                    "%s start: %14.8f %14.8f  End: %14.8f %14.8f  ArcLenClock: %14.8f  ArcLenCC: %14.8f\n",
                    testname, latStart, lonStart, latEnd, lonEnd, arcLenClock,
                    arcLenCC);
            printf(
                    "%s failed EndClock: %14.8f %14.8f  EndCC: %14.8f %14.8f  subtendedClock: %14.8f  subtendedCC: %14.8f\n",
                    testname, endClock.latitude / DEG2RAD, endClock.longitude
                            / DEG2RAD, endCC.latitude / DEG2RAD,
                    endCC.longitude / DEG2RAD, subtendedClock / DEG2RAD,
                    subtendedCC / DEG2RAD);
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

    printf("Finish testArcFromLength_Set1\n");

    return set;
}

/*
 * NAME: testArcFromLength_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the arcFromLength function.
 *
 * 		This function runs all the test cases and returns the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testArcFromLength_AllSets(TestSuite) - A test suite with the folling metrics:
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
TestSuite testArcFromLength_AllSets()
{
	TestSuite suite;
	TestSet set1;

    printf("\nStart testArcFromLength_AllSets\n");

    suite = newTestSuite("testArcFromLength_AllSets");

    set1 = testArcFromLength_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("Finish testArcFromLength_AllSets\n\n\n");

    return suite;
}

/*
 * NAME: testGetArcExtent_Set1
 *
 * DESCRIPTION:
 *
 * 		This function runs the test data created by Balakrishna Babu and uses
 * 		his original test plan to execute the actual tests.
 *
 *      The approach for testing getArcExtent is as follows:
 *      Note that this function under test returns a signed value.
 *      The sign of the returned value is the same as the sign of the arc's
 *      orientation.  In other words, if the orientation < 0 (clockwise),
 *      then arc extent < 0.  This is counterintuitive.
 *
 *      1. Select a start course and an end course.
 *      2. Calculate the difference between start and end courses moving clockwise
 *         from start and taking into account crossing 0 (or 360) azimuth.
 *      3. Repeat the previous step moving in the counter clockwise direction.
 *      4. Using the selected data in step 1 as input, compute the Arc Extent
 *         with orient = -1, and +1.
 *      5. The result from step 4 should match those calculated in steps 2 and 3.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testGetArcExtent_Set1(TestSet) - A test set with the folling metrics:
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
TestSet testGetArcExtent_Set1()
{
    double DEG2RAD = M_PI / 180.0;
    int testCaseCount = 0, setupFailureCount = 0, errorCount = 0, unverifiedCount = 0;
    ErrorSet err = 0;

    double crsCS, crsCE, arcExtentA, arcExtentB, arcExtentAExp, arcExtentBExp;
    int passedCount = 0, failedCount = 0;
    double radius;

    LLPoint center, startPt, endPt;

    TestSet set;

    long newSeed = 20080115;

    printf("Start testGetArcExtent_Set1\n");

    set = newTestSet("testGetArcExtent_Set1");

    srand(newSeed); //Initialize the random number generator

    while(testCaseCount < 10000){

    	//Create arc
    	center.latitude = randLat()*DEG2RAD;
    	center.longitude = randLon()*DEG2RAD;

    	crsCS = randAzimuth()*DEG2RAD;
    	crsCE = randAzimuth()*DEG2RAD;
    	radius = randDist();

    	//Set-up expected values
    	arcExtentAExp = computeSubtendedAngle(crsCS, crsCE, CLOCKWISE);
    	arcExtentBExp = computeSubtendedAngle(crsCS, crsCE, COUNTERCLOCKWISE);

    	err |= direct(center, crsCS, radius, &startPt, EPS);
    	err |= direct(center, crsCE, radius, &endPt, EPS);

    	if(err){
    		testCaseCount++;
    		errorCount++;
    		continue;
    	}

    	if(crsCS != crsCE){
			if(ptsAreSame(startPt,endPt,TOL)){
				arcExtentAExp = (arcExtentAExp > M_PI) ? M_2PI:0.0;
				arcExtentBExp = (arcExtentBExp > M_PI) ? M_2PI:0.0;
			}
    	}

    	//Test getArcExtent
    	err |= getArcExtent(center, radius, crsCS, crsCE, CLOCKWISE, &arcExtentA, TOL, EPS);
    	err |= getArcExtent(center, radius, crsCS, crsCE, COUNTERCLOCKWISE, &arcExtentB, TOL, EPS);

    	if(err){
    		testCaseCount++;
    		errorCount++;
    		continue;
    	}

    	testCaseCount++;
        if (arcExtentA == arcExtentAExp) {
    	/*if(fabs(arcExtentA - arcExtentAExp) < TOL){
*/
    		passedCount++;
    	} else {
    		failedCount++;
    		printf("Test-%d, center: (%3.15lf, %3.15lf), radius: %3.15lf, crsCS: %3.15lf, crsCE: %3.15lf\n",testCaseCount,center.latitude/DEG2RAD,center.longitude/DEG2RAD,radius,crsCS/DEG2RAD,crsCE/DEG2RAD);
    		printf("Test-%d, arcExtentAExp: %3.15lf, arcExtentA: %3.15lf (clockwise)\n", testCaseCount, arcExtentAExp/DEG2RAD, arcExtentA/DEG2RAD);
    	}

    	testCaseCount++;
    	if(fabs(arcExtentB - arcExtentBExp) < TOL){
    		passedCount++;
    	} else {
    		failedCount++;
    		printf("Test-%d, center: (%3.15lf, %3.15lf), radius: %3.15lf, crsCS: %3.15lf, crsCE: %3.15lf\n",testCaseCount,center.latitude/DEG2RAD,center.longitude/DEG2RAD,radius,crsCS/DEG2RAD,crsCE/DEG2RAD);
    		printf("Test-%d, arcExtentBExp: %3.15lf, arcExtentB: %3.15lf (counter-clockwise)\n", testCaseCount, arcExtentBExp/DEG2RAD, arcExtentB/DEG2RAD);
    	}

    }

    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("Finish testGetArcExtent_Set1\n");

    return set;
    }

/*
 * NAME: testGetArcExtent_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the getArcExtent function.
 *
 * 		This function runs all the test cases from set1 and set2 and returns
 * 		the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testGetArcExtent_AllSets(TestSuite) - A test suite with the folling metrics:
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
TestSuite testGetArcExtent_AllSets()
{
	TestSuite suite;
	TestSet set1;

    printf("\nStart testGetArcExtent_AllSets\n");

    suite = newTestSuite("testGetArcExtent_AllSets");

    set1 = testGetArcExtent_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("Finish testGetArcExtent_AllSets\n\n\n");

    return suite;
}
/*
 * NAME: testPtsOnArcOnTanThruPt_Set1
 *
 * DESCRIPTION:
 *
 * 		This function runs the test data created by Balakrishna Babu and uses
 * 		his original test plan to execute the actual tests.
 *
 *      The approach for testing ptsOnArcOnTanThruPt is as follows:
 *
 *      1. Select a center (C) and a radius.
 *      2. Select an azimuth and a distance from C, and compute point P using direct.
 *      3. If distance is less than radius, then there is no solution and expect
 *         zero tangent points.
 *      4. If distance equals radius, then the point itself is the one and only one
 *         tangent point.
 *      5. If distance is greater than radius, then there are two tangent points
 *         which are computed as follows.
 *      6. Compute course from center of Arc to point (crsCP) using inverse.
 *      7. Calculate a point (T) on the Arc along a geodesic from center with azimuth
 *         equal to crsCP + 90 using direct.
 *      8. Calculate the intersection points between the Arc and geodesic PT using
 *         geoArcIntx.  There should be two T and U.
 *      9. Calculate the mid point of T and U (say V).
 *      10. Compute the course from center of Arc to V and calculate the new point
 *         (T) on the Arc.
 *      11. Repeat steps 8 through 10 until point T converges.  This will be one tangent
 *         point.  Repeat with crsAP - 90 for the second tangent point.
 *      12. The results from the function should match the points from step 11.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testPtsOnArcOnTanThruPt_Set1(TestSet) - A test set with the folling metrics:
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
TestSet testPtsOnArcOnTanThruPt_Set1()
 {
    //LLPoint point, LLPoint center, double radius,
    //LLPointPair tanPt, int* n, double tol, double eps);
    double DEG2RAD = M_PI / 180.0;
//    double RAD2DEG = 180.0 / M_PI;
    double NMTOL = 0.03 / 100.0 / 1852.0;  //0.03 cm or ~1.62e-7 nm
    int testCaseCount = 0, setupFailureCount = 0, errorCount =0, unverifiedCount = 0;
    char testname[80];
    double latC, lonC, latP, lonP;
    double r1;
    double azCP;

    LLPoint p1, c1, tanP1, tanP2, midP, tanP1new, tanP2new;
    LLPointPair intx, tanP, tanPExp;
    int nExp, nX, nP;
    double latTP1, lonTP1, latTP2, lonTP2;
    double crsCP;
//    double crsPC;
    double distCP;
    double crsPT, crsTP, distTP;
    double crsCTP1, crsCTP2;
    double crs12, crs21, dist12;
    double crsCT, crsTC, distCT;
    double distConverge;
    int converged = 0, nIter;
    int passedCountZero = 0, failedCountZero = 0;
    int passedCountOne = 0, failedCountOne = 0;
    int passedCountTwo = 0, failedCountTwo = 0;
    int passedCount = 0, failedCount = 0;
    ErrorSet err = 0;

    TestSet set;

    long newSeed = 20080122;

    printf("\nStart testPtsOnArcOnTanThruPt_Set1\n");

    set = newTestSet("testPtsOnArcOnTanThruPt_Set1");

    srand(newSeed); //Initialize the random number generator

    while (testCaseCount <= 1000)
    {
        err = 0;

        //Select center and radius
        latC = randLat();
        lonC = randLon();
        r1 = 0.001 * randDist(); //about 5.40 nm

        //Select azimuth and distance from center to locate P
        azCP = randAzimuth();
        distCP = 0.002 * randDist(); //about 10.80 nm

        sprintf(testname, "TEST%-d", testCaseCount);

        c1.latitude = latC * DEG2RAD;
        c1.longitude = lonC * DEG2RAD;
        crsCP = azCP * DEG2RAD;

        //Calculate point P
        err |= direct(c1, crsCP, distCP, &p1, EPS);
        latP = p1.latitude / DEG2RAD;
        lonP = p1.longitude / DEG2RAD;

        //Initialize the expected tangent points
        tanPExp[0].latitude = 0.0;
        tanPExp[0].longitude = 0.0;
        tanPExp[1].latitude = 0.0;
        tanPExp[1].longitude = 0.0;

        //Initialize the resulting tangent points
        tanP[0].latitude = 0.0;
        tanP[0].longitude = 0.0;
        tanP[1].latitude = 0.0;
        tanP[1].longitude = 0.0;

        if (distCP < r1)
        { //no solution
            nExp = 0;
        }

        else if ((distCP - r1) < NMTOL)
        { //P is the tangent point
            nExp = 1;
            tanPExp[0] = p1;
        }

        else
        { //two points
            nExp = 2;
            crsCTP1 = crsCP + M_PI_2;
            crsCTP2 = crsCP - M_PI_2;

            //Initial guess for first tangent point
            err |= direct(c1, crsCTP1, r1, &tanP1, EPS);

            converged = 0;
            nIter = 0;
            while (!converged)
            { //0=not converged, 1=converged
                //Calculate the intersection points between the Arc and geodesic PT
                // using geoArcIntx.  There should be two T and U.
                err |= inverse(tanP1, p1, &crsTP, &crsPT, &distTP, EPS);
                err |= geoArcIntx(tanP1, crsTP, c1, r1, intx, &nX,
                        TOL, EPS);
                if (nX == 1)
                { //Found 1 point, implies tangent
                    tanP1 = intx[0];
                    converged = 1;
                }

                else
                {
                    //Calculate the mid point of T and U (say V).
                    err |= inverse(intx[0], intx[1], &crs12, &crs21,
                            &dist12, EPS);
                    err |= direct(intx[0], crs12, 0.5 * dist12, &midP, EPS);

                    //Compute the course from center of Arc to V and calculate the
                    // new point (T) on the Arc.
                    err |= inverse(c1, midP, &crsCT, &crsTC, &distCT, EPS);
                    err |= direct(c1, crsCT, r1, &tanP1new, EPS);
                    err |= invDist(tanP1, tanP1new, &distConverge, EPS);
                    tanP1 = tanP1new;
                    nIter++;
                    if (distConverge < NMTOL)
                    {
                        converged = 1;
                    }
                }
            }
            tanPExp[0] = tanP1;

            //Initial guess for second tangent point
            err |= direct(c1, crsCTP2, r1, &tanP2, EPS);

            converged = 0;
            nIter = 0;
            while (!converged)
            { //0=not converged, 1=converged
                //Calculate the intersection points between the Arc and geodesic PT
                // using geoArcIntx.  There should be two T and U.
                err |= inverse(tanP2, p1, &crsTP, &crsPT, &distTP, EPS);
                err |= geoArcIntx(tanP2, crsTP, c1, r1, intx, &nX,
                        TOL, EPS);
                if (nX == 1)
                { //Found 1 point, implies tangent
                    tanP2 = intx[0];
                    converged = 1;
                }
                else
                {
                    //Calculate the mid point of T and U (say V).
                    err |= inverse(intx[0], intx[1], &crs12, &crs21,
                            &dist12, EPS);
                    err |= direct(intx[0], crs12, 0.5 * dist12, &midP, EPS);

                    //Compute the course from center of Arc to V and calculate the
                    // new point (T) on the Arc.
                    err |= inverse(c1, midP, &crsCT, &crsTC, &distCT, EPS);
                    err |= direct(c1, crsCT, r1, &tanP2new, EPS);
                    err |= invDist(tanP2, tanP2new, &distConverge, EPS);
                    tanP2 = tanP2new;
                    nIter++;
                    if (distConverge < NMTOL)
                    {
                        converged = 1;
                    }
                }
            }
            tanPExp[1] = tanP2;

        }

        if (getMaskedError(err, getMaskAll()))
        {
            printf("Error occurred in pre-ptsOnArcOnTanThruPt err=0x%lx\n",
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

        err |= ptsOnArcOnTanThruPt(p1, c1, r1, tanP, &nP, TOL, EPS);

        if (getMaskedError(err, getMaskAll()))
        {
            printf("Error occurred in ptsOnArcOnTanThruPt err=0x%lx\n", err);
            failedCount++;
            errorCount++;
            testCaseCount++;
            continue;
        }

        latTP1 = tanP[0].latitude / DEG2RAD;
        lonTP1 = tanP[0].longitude / DEG2RAD;
        latTP2 = tanP[1].latitude / DEG2RAD;
        lonTP2 = tanP[1].longitude / DEG2RAD;

        switch (nP) {
        case 2:
            if ((nExp == nP) && ((ptsAreSame(tanP[0], tanPExp[0], TESTTOL)
                    && ptsAreSame(tanP[1], tanPExp[1], TESTTOL))
                    || (ptsAreSame(tanP[0], tanPExp[1], TESTTOL) && ptsAreSame( tanP[1], tanPExp[0], TESTTOL))))
            {
            	if (PRINT_PASSED_CASES) printf("%s C: %14.8f %14.8f  r1: %14.8f  az: %14.8f  dist: %14.8f\n",
                        testname, latC, lonC, r1, azCP, distCP);
            	if (PRINT_PASSED_CASES) printf("%s P: %14.8f %14.8f\n", testname, latP, lonP);
            	if (PRINT_PASSED_CASES) printf(
                        "%s passed nExp=%d  nP=%d  tanP1: %14.8f %14.8f  tanP2: %14.8f %14.8f\n",
                        testname, nExp, nP, latTP1, lonTP1, latTP2, lonTP2);
                passedCount++;
                passedCountTwo++;
            }
            else
            {
                printf("%s C: %14.8f %14.8f  r1: %14.8f  az: %14.8f  dist: %14.8f\n",
                        testname, latC, lonC, r1, azCP, distCP);
                printf("%s P: %14.8f %14.8f\n", testname, latP, lonP);
                printf(
                        "%s failed nExp=%d  nP=%d  tanP1: %14.8f %14.8f  tanP2: %14.8f %14.8f\n",
                        testname, nExp, nP, latTP1, lonTP1, latTP2, lonTP2);
                failedCount++;
                failedCountTwo++;
            }
            break;
        case 1:
            if ((nExp == nP) && ptsAreSame(tanP[0], tanPExp[0], TESTTOL))
            {
            	if (PRINT_PASSED_CASES) printf("%s C: %14.8f %14.8f  r1: %14.8f  az: %14.8f  dist: %14.8f\n",
                        testname, latC, lonC, r1, azCP, distCP);
            	if (PRINT_PASSED_CASES) printf("%s P: %14.8f %14.8f\n", testname, latP, lonP);
            	if (PRINT_PASSED_CASES) printf("%s passed nExp=%d  nP=%d  tanP1: %14.8f %14.8f\n",
                        testname, nExp, nP, latTP1, lonTP1);
                passedCount++;
                passedCountOne++;
            }
            else
            {
                printf("%s C: %14.8f %14.8f  r1: %14.8f  az: %14.8f  dist: %14.8f\n",
                        testname, latC, lonC, r1, azCP, distCP);
                printf("%s P: %14.8f %14.8f\n", testname, latP, lonP);
                printf("%s failed nExp=%d  nP=%d  tanP1: %14.8f %14.8f\n",
                        testname, nExp, nP, latTP1, lonTP1);
                failedCount++;
                failedCountOne++;
            }
            break;
        case 0:
            if (nExp == nP)
            {
            	if (PRINT_PASSED_CASES) printf("%s C: %14.8f %14.8f  r1: %14.8f  az: %14.8f  dist: %14.8f\n",
                        testname, latC, lonC, r1, azCP, distCP);
            	if (PRINT_PASSED_CASES) printf("%s P: %14.8f %14.8f\n", testname, latP, lonP);
            	if (PRINT_PASSED_CASES) printf("%s passed nExp=%d  nP=%d\n", testname, nExp, nP);
                passedCount++;
                passedCountZero++;
            }
            else
            {
                printf("%s C: %14.8f %14.8f  r1: %14.8f  az: %14.8f  dist: %14.8f\n",
                        testname, latC, lonC, r1, azCP, distCP);
                printf("%s P: %14.8f %14.8f\n", testname, latP, lonP);
                printf("%s failed nExp=%d  nP=%d\n", testname, nExp, nP);
                failedCount++;
                failedCountZero++;
            }
            break;
        default:
            printf("%s Unexpected Error nP=%d nExp=%d\n", testname, nP, nExp);
        }
        testCaseCount++;
    }

    /*printf("Random: %d passedZero, %d failedZero\n", passedCountZero,
            failedCountZero);
    printf("Random: %d passedOne, %d failedOne\n", passedCountOne,
            failedCountOne);
    printf("Random: %d passedTwo, %d failedTwo\n", passedCountTwo,
            failedCountTwo);
    printf("Random: %d passed, %d failed\n", passedCount, failedCount);
*/

    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("Finish testPtsOnArcOnTanThruPt_Set1\n");

    return set;
}

/*
 * NAME: testPtsOnArcOnTanThruPt_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the ptsOnArcOnTanThruPt function.
 *
 * 		This function runs all the test cases from set1 and set2 and returns
 * 		the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testGetArcExtent_AllSets(TestSuite) - A test suite with the folling metrics:
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
TestSuite testPtsOnArcOnTanThruPt_AllSets()
{
	TestSuite suite;
	TestSet set1;

    printf("\nStart testPtsOnArcOnTanThruPt_AllSets\n");

    suite = newTestSuite("testPtsOnArcOnTanThruPt_AllSets");

    set1 = testPtsOnArcOnTanThruPt_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("Finish testPtsOnArcOnTanThruPt_AllSets\n\n\n");

    return suite;
}

/*
 * NAME: testGeoTanToTwoCircles_Set1
 *
 * DESCRIPTION:
 *      Tests geoTanToTwoCircles function.
 *
 * 		This function runs the test cases created by Rich Snow
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testGeoTanToTwoCircles_Set1(TestSet) - A test set with the folling metrics:
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
TestSet testGeoTanToTwoCircles_Set1()
{
    double DEG2RAD = M_PI / 180.0;
    double ANGLE_TOL = 5.0e-8;//1.0e-5//degrees
    LLPoint center1, center2;
    Geodesic tanLines[2];
    ArcDirection dir1, dir2;

    double r1,r2;
    double az12,az21,minSubAng1,minSubAng2, raddist1, raddist2;
    double crsFromPoint;
    double startDist;
    double endDist;
    double distFromPoint1, distFromPoint2;
    double latcent1, loncent1, lineAz;
    double temp1, temp2, distcent;
    LLPoint projPoint1;
    LLPoint projPoint2;

    int j,jj;
    int passed = 1;

    ErrorSet err = 0;
    long newSeed = 20080523;

    int passedCount=0, failedCount=0;
    int errorCount = 0;
    int unverifiedCount = 0;
    int setupFailureCount = 0;
    int testCaseCount = 0;

    TestSet set;

    printf("\nStart testGeoTanToTwoCircles_Set1\n");

    set = newTestSet("testGeoTanToTwoCircles_Set1");

    srand(newSeed);

    for (jj = 0; jj < 200; jj++)
    {
      err = 0;
      latcent1 = DEG2RAD * randLat();
      loncent1 = DEG2RAD * randLon();
      center1.latitude = latcent1;
      center1.longitude = loncent1;
      if (jj < 50)
        r1 = 30.0 + 0.01 * randDist();
      else 
        r1 = 1.0 + 0.01 * randDist();
      lineAz = DEG2RAD * randAzimuth();
      if ((jj < 25) || (jj > 74))
        r2 = 30.0 + 0.01 * randDist();
      else 
        r2 = 1.0 + 0.01 * randDist();
      distcent =  0.01 * randDist();
      if ( (rand() % 2) == 0)
        dir1 = CLOCKWISE;
      else 
        dir1 = COUNTERCLOCKWISE;
      if ( (rand() % 2) == 0)
        dir2 = CLOCKWISE;
      else 
        dir2 = COUNTERCLOCKWISE;
      if (jj == 0)
        distcent = 1.0;
      else
        distcent = r1 + r2 + distcent;
      err |= direct(center1, lineAz, distcent, &center2, EPS);
      //printf("jj = %d r1 = %4.8f dir1 = %d r2 = %4.8f dir2 = %d distcent %4.8f\n",jj,r1,dir1,r2,dir2,distcent);

        if (getMaskedError(err, getMaskAll()))
        {
            printf("Error occurred in pre-geoTanToTwoCircles err=0x%lx\n", err);
            testCaseCount++;
            errorCount++;
            setupFailureCount++;
            continue;
        }

        err |= geoTanToTwoCircles(center1, r1, dir1, center2, r2, dir2, tanLines, TOL, EPS);

        if (err)
        {
            if ((jj == 0) && (err == CIRCLE_INSIDE_CIRCLE_ERR))
            {
              testCaseCount++;
              passedCount++;
              continue;
            }
            else
            {
              printf("Error: %ld\n%s",err,formatErrorMessage(err));
              testCaseCount++;
              errorCount++;
              failedCount++;
              continue;
            }
        }
        else
        {

            passed = 1;

            for (j=0;j<2;j++)
            {
                err |= invCrs(tanLines[j].startPoint, tanLines[j].endPoint, &az12,&az21,EPS);
                if (j == 0)
                err |= invCrs(tanLines[j].startPoint, center1, &temp1,&temp2,EPS);
                else
                err |= invCrs(tanLines[j].startPoint, center2, &temp1,&temp2,EPS);
                err |= minSubtendedAngle(az12, temp1, &minSubAng1);
                //printf("tan1 %4.12f\n",minSubAng1 / DEG2RAD);
                if (j == 0)
                err |= invCrs(tanLines[j].endPoint, center2, &temp1,&temp2,EPS);
                else
                err |= invCrs(tanLines[j].endPoint, center1, &temp1,&temp2,EPS);
                err |= minSubtendedAngle(az21, temp1, &minSubAng2);
              
                //printf("tan2 %4.12f\n",minSubAng2 / DEG2RAD);
                err |= projectToGeo(tanLines[j].startPoint,az12,center1,&projPoint1,&crsFromPoint,
                                          &distFromPoint1,TOL/10.0,EPS);
                err |= projectToGeo(tanLines[j].startPoint,az12,center2,&projPoint2,&crsFromPoint,
                                          &distFromPoint2,TOL/10.0,EPS);

                if (j==0)
                {
                    err |= invDist(projPoint1, tanLines[j].startPoint, &startDist, EPS);
                    err |= invDist(projPoint2, tanLines[j].endPoint, &endDist, EPS);
                    err |= invDist(tanLines[j].startPoint, center1, &raddist1, EPS);
                    err |= invDist(tanLines[j].endPoint, center2, &raddist2, EPS);
                }
                else
                {
                    err |= invDist(projPoint2, tanLines[j].startPoint, &startDist, EPS);
                    err |= invDist(projPoint1, tanLines[j].endPoint, &endDist, EPS);
                    err |= invDist(tanLines[j].startPoint, center2, &raddist2, EPS);
                    err |= invDist(tanLines[j].endPoint, center1, &raddist1, EPS);
                }

                if (fabs(minSubAng1 / DEG2RAD -90.0) >= ANGLE_TOL)
                  passed = 0;

                if (fabs(minSubAng2 / DEG2RAD -90.0) >= ANGLE_TOL)
                  passed = 0;

                if (fabs(r1 - raddist1) > TESTTOL)
                  passed = 0;
                 if (fabs(r2 - raddist2) > TESTTOL)
                  passed = 0;

                if  (startDist > TESTTOL)
                {
                    printf("StartDistance %d > TOL: %e\n",j,startDist);
                    passed = 0;
                }
                if ( endDist > TESTTOL )
                {
                    printf("EndDistance %d > TOL: %e\n",j,endDist);
                    passed = 0;
                }

            }

            if (getMaskedError(err, getMaskAll()))
            {
            	printf("Error occurred in post-geoTanToTwoCircles err=0x%lx\n", err);
                testCaseCount++;
                errorCount++;
                failedCount++;
                continue;
            }

            if (passed)
            {
                passedCount++;
            }
            else
            {
                failedCount++;
            }

            testCaseCount++;
        }

    } //for jj

    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("Finish testGeoTanToTwoCircles_Set1\n");

    return set;
}

/*
 * NAME: testGeoTanToTwoCircles_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the geoTanToTwoCircles function.
 *
 * 		This function runs all the test cases and returns the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testGeoTanToTwoCircles_AllSets(TestSuite) - A test suite with the folling metrics:
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
TestSuite testGeoTanToTwoCircles_AllSets()
{
	TestSuite suite;
	TestSet set1;

    printf("\nStart testGeoTanToTwoCircles_AllSets\n");

    suite = newTestSuite("testGeoTanToTwoCircles_AllSets");

    set1 = testGeoTanToTwoCircles_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("Finish testGeoTanToTwoCircles_AllSets\n\n\n");

    return suite;
}


/* NAME: testComputeSubtendedAngle_Set1
 *
 * DESCRIPTION:
 *      Runs tests to determine the subtended angle given a startCrs, endCrs, and orientation.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testComputeSubtendedAngle_Set1(TestSet) - A test set with the folling metrics:
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
TestSet testComputeSubtendedAngle_Set1()
{
//    double DEG2RAD = M_PI / 180.0;
    double RAD2DEG = 180.0 / M_PI;
//    double AZDEGTOL = 5.0e-7; //5e-7 deg or 8.72665e-9 rad or 0.162 cm at 100 nm

    double startCrs, endCrs;
    ArcDirection orientation; //-1=counter clockwise, 1=clockwise

    double expectedResult, expectedResultRad, subtendedAngleRad;

    int passedCount = 0, failedCount = 0, errorCount = 0;
    int testCaseCount = 0;
    int setupFailureCount = 0;
    int unverifiedCount = 0;


    TestSet set;

    long newSeed = 20080101;
    srand(newSeed); //Initialize the random number generator


    printf("Start testComputeSubtendedAngle_Set1\n");

    set = newTestSet("testComputeSubtendedAngle_Set1");

    for (testCaseCount = 0; testCaseCount < 10000; testCaseCount++) {
                    startCrs = randDouble(0.0, M_2PI);
                    endCrs = randDouble(0.0, M_2PI);
                    orientation = COUNTERCLOCKWISE;  // counter-clockwise
                    if (randDouble(-1.0, 1.0) < 0) {
                        orientation = CLOCKWISE;  // clockwise
                    }

                    if (startCrs == endCrs){
                    	expectedResult = orientation*M_2PI;
                    } else {
						if (orientation == COUNTERCLOCKWISE) {
							if (startCrs > endCrs) {
								expectedResult = startCrs - endCrs;
							}
							else {
								expectedResult = M_2PI - (endCrs - startCrs);
							}
							expectedResult *= -1;
						}
						else {
							if (startCrs > endCrs) {
								expectedResult = M_2PI - (startCrs - endCrs);
							}
							else {
								expectedResult = endCrs - startCrs;
							}
						}
                    }

        subtendedAngleRad = computeSubtendedAngle(startCrs,endCrs,orientation);
        expectedResultRad = expectedResult;
        expectedResult *= RAD2DEG;

        if (subtendedAngleRad == expectedResultRad)
//        if (fabs(subtendedAngleRad - expectedResultRad) < (AZDEGTOL * DEG2RAD))
        {
        	if (PRINT_PASSED_CASES)
        		printf("%d startCrs: %20.15f  endCrs: %20.15f  orientation: %d  expected: %20.15f  returned: %20.15f\n",
                    testCaseCount, startCrs, endCrs, orientation, expectedResult, subtendedAngleRad*RAD2DEG);

        	if (PRINT_PASSED_CASES)
        		printf("%d passed\n", testCaseCount);
            passedCount++;
        }
        else
        {
            printf("%d startCrs: %20.15f  endCrs: %20.15f  orientation: %d  expected: %20.15f  returned: %20.15f\n",
            		testCaseCount , startCrs, endCrs, orientation, expectedResult, subtendedAngleRad*RAD2DEG);
            printf("%d failed\n", testCaseCount);
            failedCount++;
        }
    }

    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("Finish testComputeSubtendedAngle_Set1\n");

    return set;
}


/*
 * NAME: testComputeSubtendedAngle_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the computeSubtendedAngle function.
 *
 * 		This function runs all the test cases and returns the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testComputeSubtendedAngle_AllSets(TestSuite) - A test suite with the folling metrics:
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
TestSuite testComputeSubtendedAngle_AllSets()
{
	TestSuite suite;
	TestSet set1;

    printf("\nStart testComputeSubtendedAngle_AllSets\n");

    suite = newTestSuite("testComputeSubtendedAngle_AllSets");

    set1 = testComputeSubtendedAngle_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("Finish testComputeSubtendedAngle_AllSets\n\n\n");

    return suite;
}

/* NAME: testGeoTanToArcAtAngleToGeo_Set1
 *
 * DESCRIPTION:
 *      Runs tests to determine the tangent line given a geodesic, arc, and angle.
 *      Tests developed by Joe Hopper
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testGeoTanToArcAtAngleToGeo_Set1(TestSet) - A test set with the following metrics:
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
TestSet testGeoTanToArcAtAngleToGeo_Set1()
{
    char testName[80];

    int passedCount = 0, failedCount = 0, errorCount = 0;
    int testCaseCount = 0;
    int setupFailureCount = 0;
    int unverifiedCount = 0;

	double DEG2RAD = M_PI/180.0;
	double tol = 1.37e-9, eps = 1e-20;

	LLPoint anglePtExp, tanPtExp;
	Arc testArc;
	Geodesic testGeo;
	double angle;
	double crsFromAnglePtToGeoStart, crsFromGeoStartToAnglePt;
	double crsFromAnglePtToTanPt, crsFromTanPtToAnglePt;
	double crsFromTanPtToArcCenter;
	double geoDist, tanDist;

	Geodesic tangentLine;
	int tangentLineLocation;

	double anglePtErr, tanPtErr;

	int k, i, j;
	long err = 0;
    long newSeed = 20080523;

    TestSet set;

    set = newTestSet("testGeoTanToArcAtAngleToGeo_Set1");

    srand(newSeed); //Initialize the random number generator

//	printf("anglePtErr,tanPtErr,anglePtExplat,anglePtExplon,anglePtlat,anglePtlon,tanPtExplat,tanPtExplon,tanPtlat,tanPtlon,angle,geoDist,tanDist,radius,err\n");
    for(k = 0; k < 2500; k++){
    	err = 0;

    	/* Set up random test data */
		anglePtExp.latitude = (-80.0 + 160.0*rand()/RAND_MAX)*DEG2RAD;
		anglePtExp.longitude = (-180.0 + 360.0*rand()/RAND_MAX)*DEG2RAD;
		crsFromAnglePtToGeoStart = (0.0 + (360.0 - 0.0)*rand()/RAND_MAX)*DEG2RAD;
		angle = (1.0 + (179.0 - 1.0)*rand()/RAND_MAX)*DEG2RAD;
		testArc.radius = 1.0e-8 + (200.0 - 1.0e-8)*rand()/RAND_MAX;
		geoDist = -2000.0 + (4000.0 - 0.0)*rand()/RAND_MAX;
		tanDist = 0.0 + (1500.0 - 0.0)*rand()/RAND_MAX;

		/* Set up Geo */
		err |= direct(anglePtExp, crsFromAnglePtToGeoStart, geoDist, &testGeo.startPoint, eps);
		if(geoDist == 0.0){
			testGeo.startAz = crsFromAnglePtToGeoStart + M_PI;
		} else {
			err |= inverse(anglePtExp, testGeo.startPoint, &crsFromAnglePtToGeoStart, &crsFromGeoStartToAnglePt, NULL, eps);
			if(geoDist < 0.0){
				crsFromGeoStartToAnglePt += M_PI;
				crsFromAnglePtToGeoStart += M_PI;
			}
			testGeo.startAz = crsFromGeoStartToAnglePt;
		}

		/* Loop through 4 casses */
		for(i = -1; i < 2; i += 2){
			for(j = -1; j < 2; j += 2){
				if(i == -1 )
			    sprintf(testName, "TEST%-d", testCaseCount);

				/* Set up Arc */
                                /*if (i * j > 0)
				crsFromAnglePtToTanPt = crsFromAnglePtToGeoStart + i*angle;
                                else*/
				crsFromAnglePtToTanPt = crsFromAnglePtToGeoStart - i*angle;
				err |= direct(anglePtExp, crsFromAnglePtToTanPt, tanDist, &tanPtExp, eps);
				if(tanDist == 0.0){
					crsFromTanPtToAnglePt = crsFromAnglePtToTanPt + M_PI;
				} else {
					err |= inverse(anglePtExp, tanPtExp, &crsFromAnglePtToTanPt, &crsFromTanPtToAnglePt, NULL, eps);
				}
				crsFromTanPtToArcCenter = crsFromTanPtToAnglePt + j*M_PI_2;
				testArc.dir = static_cast<ArcDirection>(j);
                                /*err |= crsIntx(testGeo.startPoint, testGeo.startAz, NULL, &dist, tanPtExp, crsFromTanPtToArcCenter, NULL, NULL, NULL, tol, eps);
                                if ((i * j < 0) && (testArc.radius < dist))
                                  testArc.radius = dist + 1.0;*/
				err |= direct(tanPtExp, crsFromTanPtToArcCenter, testArc.radius, &testArc.centerPoint, eps);

				/* Test function */
				err |= geoTanToArcAtAngleToGeo(testArc, testGeo, angle, &tangentLine, &tangentLineLocation, tol, eps);
				if(tangentLineLocation == 0){
//					printf("angle %lf\n", angle/DEG2RAD);
					continue;
				}
				testCaseCount++;

				if(err){
					errorCount++;
					printf(" %s Arcdir %d CentSide %d angle %f radius %f tanDist %f err=0x%lx\n",testName,testArc.dir,-i,angle/DEG2RAD,testArc.radius,tanDist,err);
				}

				err |= inverse(anglePtExp, tangentLine.endPoint, NULL, NULL, &anglePtErr, eps);
				err |= inverse(tanPtExp, tangentLine.startPoint, NULL, NULL, &tanPtErr, eps);

		//		if(err == 0){
		//			printf("%3.15e,%3.15e,%3.15lf,%3.15lf,%3.15lf,%3.15lf,%3.15lf,%3.15lf,%3.15lf,%3.15lf,%3.15lf,%3.15lf,%3.15lf,%3.15lf,%lx\n",
		//					anglePtErr, tanPtErr, anglePtExp.latitude/DEG2RAD, anglePtExp.longitude/DEG2RAD, tangentLine.endPoint.latitude/DEG2RAD, tangentLine.endPoint.longitude/DEG2RAD,
		//					tanPtExp.latitude/DEG2RAD, tanPtExp.longitude/DEG2RAD, tangentLine.startPoint.latitude/DEG2RAD, tangentLine.startPoint.longitude/DEG2RAD,
		//					angle/DEG2RAD, geoDist, tanDist, testArc.radius, err);
		//		}

				if(!(ptsAreSame(anglePtExp,tangentLine.endPoint,TESTTOL) && ptsAreSame(tanPtExp,tangentLine.startPoint,TESTTOL))){
					printf("%s\n", testName);
					printf("ArcDirection: %d, CenterSide: %d\n", testArc.dir, -i);

					printf("anglePtErr %3.15e\n", anglePtErr);
					printf("tanPtErr   %3.15e\n", tanPtErr);
					printf("\n");
					printf("anglePtExp %3.15lf, %3.15lf\n", anglePtExp.latitude/DEG2RAD, anglePtExp.longitude/DEG2RAD);
					printf("anglePt    %3.15lf, %3.15lf\n", tangentLine.endPoint.latitude/DEG2RAD, tangentLine.endPoint.longitude/DEG2RAD);

					printf("tanPtExp   %3.15lf, %3.15lf\n", tanPtExp.latitude/DEG2RAD, tanPtExp.longitude/DEG2RAD);
					printf("tanPt      %3.15lf, %3.15lf\n", tangentLine.startPoint.latitude/DEG2RAD, tangentLine.startPoint.longitude/DEG2RAD);
					printf("\n");
					printf("angle: %3.15lf, geoDist: %3.15lf, tanDist: %3.15lf, radius: %3.15lf\n", angle/DEG2RAD, geoDist, tanDist, testArc.radius);
					printf("\n");
					printf("testArc.centerPoint = {%3.18lf, %3.18lf}\n", testArc.centerPoint.latitude, testArc.centerPoint.longitude);
					printf("radius = %3.18lf\n", testArc.radius);
					printf("testArc.dir = %d\n", testArc.dir);
					printf("testGeo.startPoint = {%3.18lf, %3.18lf}\n", testGeo.startPoint.latitude, testGeo.startPoint.longitude);
					printf("testGeo.startAz = %3.18lf\n", testGeo.startAz);
					printf("angle = %3.18lf\n", angle);
					printf("\n");
					failedCount++;
				} else {
					passedCount++;
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

    printf("Finish %s\n", set.name.c_str());

    return set;

}

/* NAME: testGeoTanToArcAtAngleToGeo_Set2
 *
 * DESCRIPTION:
 *      Runs tests to determine the tangent line given a geodesic, arc, and angle.
 *      Tests developed by Rich Snow
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testGeoTanToArcAtAngleToGeo_Set3(TestSet) - A test set with the following metrics:
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
TestSet testGeoTanToArcAtAngleToGeo_Set2()
{

    double DEG2RAD = M_PI / 180.0;
    const double RAD2DEG = 180.0 / M_PI;
    double ANGLE_TOL = 5.0e-8;//1.0e-5;//degrees

    char testName[80];

    ArcDirection arcDir;
    LineType lineType;
    double angle;
    double tanstlatexp, tanstlonexp, crstscent, crscentts, temp;
    double geodist;


    LLPoint arcCenter, arcStart, arcEnd;
    Arc arc;
    LLPoint geoStart, geoEnd;
    LLPoint tanstexp, tanendexp, tempLLPoint;
    Geodesic geo,tangeoexp;


    Geodesic tangentGeo;
    int tangentGeoLocation, tangentGeoLocationexp;
    int tanLineOnGeo = 0;
    int tanLineOnArc = 0;
    LLPointPair tanLineAndArcIntersections;
    int tanLineAndArcIntersectionCnt = -1;
    double crsFromArcCenterToTanLineStart, crsFromTanLineStartToArcCenter;
    double tangencyAngle;
    double crsFromTanLineEndToGeoEnd, crsFromGeoEndToTanLineEnd;
    double interiorAngle;
    double latC, lonC, r1, az1, az2, crs1, crs2, tempcrs;
    double crs31, dist13, crs32, dist23;
    int outputMatlab = 0;//0 = don't output matlab code to the console

    int fail;
    int testNum = 0;

    int passedCount = 0, failedCount = 0, errorCount = 0;
    int testCaseCount = 0;
    int setupFailureCount = 0;
    int unverifiedCount = 0;
    long err = 0;
    long newSeed = 20080523;

    TestSet set;

    set = newTestSet("testGeoTanToArcAtAngleToGeo_Set2");

    printf("Start %s\n", set.name.c_str());

    srand(newSeed); //Initialize the random number generator

    while (testNum < 28)
    {
      err = 0;
      fail = 0;
      testNum++;
      sprintf(testName, "TEST%-d",testNum);

      if (testNum == 2 || testNum == 4 || testNum == 8 || testNum == 10 || testNum == 16 || testNum == 18 || testNum == 22 || testNum == 24)
      {
        //printf("%s not used\n",testName);
        continue; 
      }
      //generate input data for Arc
      latC = DEG2RAD*randLat();
      lonC = DEG2RAD*randLon();
      err |= createPt(&arcCenter, latC, lonC);
      r1 = (double) ((rand() % 100) + 1);
      if (testNum >= 1 && testNum <= 7)
      {
        crscentts = 315.0;
        //if (testNum == 7)
        //crscentts = 335.0;
        //angle = (double) ((rand() % 90) + 1); //degrees (0,90]
        angle = 45.0;
        arcDir = CLOCKWISE;
        tangentGeoLocationexp = -1;
      }
      else if (testNum >= 8 && testNum <=14)
      {
        crscentts = 45.0;
        //angle = (double) ((rand() % 90) + 1); //degrees (0,90]
        angle = 45.0; 
        arcDir = COUNTERCLOCKWISE;
        tangentGeoLocationexp = 1;
      }
      else if (testNum >= 15 && testNum <= 21)
      {
        crscentts = 225.0;
        //angle = (double) ((rand() % 90) + 1); //degrees (0,90]
        angle = 45.0; 
        arcDir = COUNTERCLOCKWISE;
        tangentGeoLocationexp = 1;
      }
      else if (testNum >= 22 && testNum <= 28)
      {
        crscentts = 135.0;
        //angle = (double) ((rand() % 90) + 1); //degrees (0,90]
        angle = 45.0; 
        arcDir = CLOCKWISE;
        tangentGeoLocationexp = -1;
      }
      crscentts = crscentts * DEG2RAD;
      err |= direct(arcCenter, crscentts, r1, &tanstexp, EPS);
      tanstlatexp = tanstexp.latitude;
      tanstlonexp = tanstexp.longitude;
      err |= invCrs(tanstexp, arcCenter, &crstscent, &temp, EPS);
      lineType = SEGMENT;
      if (testNum == 3 || testNum == 11 || testNum == 17 || testNum == 25)
      {
        if (testNum == 3 || testNum == 17)
          err |= direct(arcCenter, M_PI_2, r1, &tempLLPoint, EPS);
        else if (testNum == 11 || testNum == 25)
          err |= direct(arcCenter, 3.0 * M_PI_2, r1, &tempLLPoint, EPS);
        err |= invCrs(tempLLPoint, arcCenter, &tempcrs, &temp, EPS);
        if (testNum == 3 || testNum == 25)
          err |= crsIntx(tanstexp, crstscent - M_PI_2, &crs31, &dist13, tempLLPoint, tempcrs + M_PI_2, &crs32, &dist23, &tanendexp, TOL, EPS);
        else if (testNum == 11 || testNum == 17)
          err |= crsIntx(tanstexp, crstscent + M_PI_2, &crs31, &dist13, tempLLPoint, tempcrs - M_PI_2, &crs32, &dist23, &tanendexp, TOL, EPS);

        geoStart.latitude = tempLLPoint.latitude;
        geoStart.longitude = tempLLPoint.longitude;
        geoEnd.latitude = tanendexp.latitude;
        geoEnd.longitude = tanendexp.longitude;  
        if (testNum == 11)
        {
          err |= direct(geoStart, tempcrs + M_PI_2, 1.0, &tempLLPoint, EPS);
          geoStart.latitude = tempLLPoint.latitude;
          geoStart.longitude = tempLLPoint.longitude;
        }

        //Extend the geo 1 nmi, so tanendexp and geoEnd are not same point.
          err |= direct(geoEnd, crs32 + M_PI, 1.0, &tempLLPoint, EPS);
          geoEnd.latitude = tempLLPoint.latitude;
          geoEnd.longitude = tempLLPoint.longitude;
          err |= invCrs(tanendexp, geoEnd, &crs32, &temp, EPS);
          err |= minSubtendedAngle(crs32, crs31 + M_PI, &angle); 
          angle = angle * RAD2DEG;
        err |= createGeo(&tangeoexp, tanstexp, tanendexp, lineType, EPS);
      }

      if (testNum == 7 || testNum == 14 || testNum == 21 || testNum == 28)
      {
        if (testNum == 7)  
          err |= crsIntx(tanstexp, crstscent - M_PI_2, &crs31, &dist13, arcCenter, 5.0 * DEG2RAD, &crs32, &dist23, &tanendexp, TOL, EPS);
        else if (testNum == 14)
          err |= crsIntx(tanstexp, crstscent + M_PI_2, &crs31, &dist13, arcCenter, 5.0 * DEG2RAD, &crs32, &dist23, &tanendexp, TOL, EPS);
        else if (testNum == 21)
            err |= crsIntx(tanstexp, crstscent + M_PI_2, &crs31, &dist13, arcCenter, 175.0 * DEG2RAD, &crs32, &dist23, &tanendexp, TOL, EPS);
        else if (testNum == 28)
            err |= crsIntx(tanstexp, crstscent - M_PI_2, &crs31, &dist13, arcCenter, 175.0 * DEG2RAD, &crs32, &dist23, &tanendexp, TOL, EPS);
        geoStart.latitude = arcCenter.latitude;
        geoStart.longitude = arcCenter.longitude;
        geoEnd.latitude = tanendexp.latitude;
        geoEnd.longitude = tanendexp.longitude;  

        //Extend the geo 1 nmi, so tanendexp and geoEnd are not same point.
          err |= direct(geoEnd, crs32 + M_PI, 1.0, &tempLLPoint, EPS);
          geoEnd.latitude = tempLLPoint.latitude;
          geoEnd.longitude = tempLLPoint.longitude;
          err |= invCrs(tanendexp, geoEnd, &crs32, &temp, EPS);
          err |= minSubtendedAngle(crs32, crs31 + M_PI, &angle); 
          angle = angle * RAD2DEG;
        err |= createGeo(&tangeoexp, tanstexp, tanendexp, lineType, EPS);
      }

      if (testNum == 1 || testNum == 5 || testNum == 6 || testNum == 9 || testNum == 12 || testNum == 13 || testNum == 15 || testNum == 19 || testNum == 20 || testNum == 23 || testNum == 26 || testNum == 27)
      {
        if (testNum == 1 || testNum == 23)
         err |= direct(tanstexp, crstscent - M_PI_2, 3.0 * r1, &tanendexp, EPS);
        else if (testNum == 5 || testNum == 27)
          err |= direct(tanstexp, crstscent - M_PI_2, 1.2 * r1, &tanendexp, EPS);
        else if (testNum == 6 || testNum == 26)
          err |= direct(tanstexp, crstscent - M_PI_2, 0.5 * r1, &tanendexp, EPS);
        else if (testNum == 9 || testNum == 15)
         err |= direct(tanstexp, crstscent + M_PI_2, 3.0 * r1, &tanendexp, EPS);
        else if (testNum == 13 || testNum == 19)
          err |= direct(tanstexp, crstscent + M_PI_2, 1.2 * r1, &tanendexp, EPS);
        else if (testNum == 12 || testNum == 20)
          err |= direct(tanstexp, crstscent + M_PI_2, 0.5 * r1, &tanendexp, EPS);
        err |= createGeo(&tangeoexp, tanstexp, tanendexp, lineType, EPS);
        geodist = (double) ((rand() % 50) + 50);
        geoStart.latitude = tanendexp.latitude;
        geoStart.longitude = tanendexp.longitude;

        if (testNum == 1 || testNum == 5 || testNum == 6 || testNum == 23 || testNum == 26 || testNum == 27)
          err |= direct(tanendexp, tangeoexp.endAz - angle*DEG2RAD, geodist, &geoEnd, EPS);
        else if (testNum == 9 || testNum == 12 || testNum == 13 || testNum == 15 || testNum == 19 || testNum == 20)
          err |= direct(tanendexp, tangeoexp.endAz + angle*DEG2RAD, geodist, &geoEnd, EPS);
      }

      az1 = (double) ((rand() % 90) + 1); //degrees (0,90]
      az2 = (double) ((rand() % 90) + 1); //degrees (0,90]
      crs1 = az1 * DEG2RAD;
      crs2 = az2 * DEG2RAD;


      crs1 = modpos(crscentts - crs1, M_2PI);
      crs2 = modpos(crscentts + crs2, M_2PI);

      if (arcDir == CLOCKWISE)
      {
        err |= direct(arcCenter, crs1, r1, &arcStart, EPS);
        err |= direct(arcCenter, crs2, r1, &arcEnd, EPS);
      }
      else
      {
        err |= direct(arcCenter, crscentts + crs1, r1, &arcStart, EPS);
        err |= direct(arcCenter, crscentts - crs2, r1, &arcEnd, EPS);
      }


      err |= createGeo(&geo, geoStart, geoEnd, lineType, EPS);
      err |= createArc(&arc, arcCenter, arcStart, arcEnd, arcDir, TOL, EPS);
    if (outputMatlab) displayMatlabGeo(geo, "geo", 0);
    if (outputMatlab) displayMatlabArc(arc, "arc", 0);


	if (getMaskedError(err, getMaskAll()))
	{
		printf("Error occurred in pre-WGS84FindTangentLineAtAngle err=0x%lx\n", err);
		errorCount++;
		setupFailureCount++;
		testCaseCount++;
	} else {

        //RUN TEST
        err = geoTanToArcAtAngleToGeo(arc, geo, angle*DEG2RAD, &tangentGeo, &tangentGeoLocation, TOL, EPS);
        //The following line is used when checking that the test setup is correct
	//err |= createGeo(&tangentGeo, tanstexp, tanendexp, lineType, EPS);
        if (outputMatlab) displayMatlabGeo(tangentGeo, "tangeo", 0);

        if (getMaskedError(err, getMaskAll()))
        {
            printf("Error occurred in geoTanToArcAtAngleToGeo err=0x%lx\n", err);
            errorCount++;
            failedCount++;
            testCaseCount++;
        } else {

        	err = 0;

        	//treat the input geo as infinite
        	tanLineOnGeo = ptIsOnGeo(geo.startPoint, geo.endPoint, tangentGeo.endPoint, INFINITE, &err, TOL, EPS);

        	//treat the input arc as a circle
        	tanLineOnArc = ptIsOnArc(arc.centerPoint, arc.radius, arc.startAz, arc.startAz, arc.dir, tangentGeo.startPoint, &err, TOL, EPS);

        	//find the intersections of the tangent line and arc
        	err |= geoArcIntx(tangentGeo.endPoint, tangentGeo.endAz + M_PI, arc.centerPoint, arc.radius, tanLineAndArcIntersections, &tanLineAndArcIntersectionCnt, TOL, EPS);

        	//find the angle at the tangency of the arc with the tangent line
        	err |= invCrs(arc.centerPoint, tangentGeo.startPoint, &crsFromArcCenterToTanLineStart, &crsFromTanLineStartToArcCenter, EPS);
        	err |= minSubtendedAngle(crsFromTanLineStartToArcCenter, tangentGeo.startAz, &tangencyAngle);

        	//find the interior angle formed
        	err |= invCrs(tangentGeo.endPoint, geo.endPoint, &crsFromTanLineEndToGeoEnd, &crsFromGeoEndToTanLineEnd, EPS);
        	err |= minSubtendedAngle(crsFromTanLineEndToGeoEnd, tangentGeo.endAz + M_PI, &interiorAngle);
                if (interiorAngle > M_PI_2)
                  interiorAngle = M_PI - interiorAngle;

        	//check if errors occurred while testing result
        	if (getMaskedError(err, getMaskAll()))
        	{
        		printf("Error occurred in pre-geoTanToArcAtAngleToGeo err=0x%lx\n", err);
        		errorCount++;
        		setupFailureCount++;
        		testCaseCount++;
        	} else {

        		//test that tangent line lies on input geo
        		if (!tanLineOnGeo)
        			fail = 1;

        		//test that tangent line lies on arc
        		if (!tanLineOnArc)
        			fail = 1;

        		//test tangency of line with respect to intersection points
        		if (tanLineAndArcIntersectionCnt == 0)
        		{
        			fail = 1;
        		}
        		else if ( tanLineAndArcIntersectionCnt == 2 && (!ptsAreSame(tanLineAndArcIntersections[0], tanLineAndArcIntersections[1], TESTTOL)) )
        		{
					fail = 1;
        		}

        		//test tangency of line with respect to angle formed with the arc
        		if (fabs(fabs(tangencyAngle*RAD2DEG) - M_PI_2*RAD2DEG) > ANGLE_TOL)
        		{
        			fail = 1;
        		}

        		//test interior angle formed with geo
        		if (fabs(fabs(interiorAngle*RAD2DEG) - angle) > ANGLE_TOL)
        		{
        			fail = 1;
        		}
                    
                        if (!ptsAreSame(tanstexp, tangentGeo.startPoint, TESTTOL) || !ptsAreSame(tanendexp, tangentGeo.endPoint, TESTTOL))
                        {
                          fail = 1;
                        }

                        if (tangentGeoLocation != tangentGeoLocationexp)
                        {
                          fail = 1;
                        }

				if (fail)
				{
                                  printf("%s tanLineOnGeo %d tanLineOnArc %d IntCnt %d arcDir %d\n",testName,tanLineOnGeo,tanLineOnArc,tanLineAndArcIntersectionCnt,arcDir);              
      printf("geo.startAz %16.12f\n",geo.startAz*RAD2DEG);
      printf("crscentts %16.12f crs1 %16.12f crs2 %16.12f\n",crscentts*RAD2DEG,crs1*RAD2DEG,crs2*RAD2DEG);
      printf("crstscent %16.12f tangeoexp.startAz %16.12f tangeoexp.endAz %16.12f\n",crstscent*RAD2DEG,tangeoexp.startAz*RAD2DEG,tangeoexp.endAz*RAD2DEG); 
      printf("crstscentcalc %16.12f tangeocalc.startAz %16.12f tangeocalc.endAz %16.12f crsFromTanLineEndToGeoEndcalc %16.12f\n",crsFromTanLineStartToArcCenter*RAD2DEG,tangentGeo.startAz*RAD2DEG,tangentGeo.endAz*RAD2DEG,crsFromTanLineEndToGeoEnd*RAD2DEG); 
      printf("tanstlatexp %16.12f tanstlonexp %16.12f\n",tanstlatexp*RAD2DEG,tanstlonexp*RAD2DEG);
      printf("tanendlatexp %16.12f tanendlonexp %16.12f\n",tanendexp.latitude*RAD2DEG,tanendexp.longitude*RAD2DEG);
                                  printf("tanAngle %16.12f intAngle %16.12f angle %16.12f\n",tangencyAngle*RAD2DEG,interiorAngle*RAD2DEG,angle);
      printf("tangentGeoLocation %d tangentGeoLocationexp %d\n",tangentGeoLocation,tangentGeoLocationexp);
					printf("\n%s Failed\n", testName);
					failedCount++;
				} else
                                {
                                  /*printf("%s tanLineOnGeo %d tanLineOnArc %d IntCnt %d arcDir %d\n",testName,tanLineOnGeo,tanLineOnArc,tanLineAndArcIntersectionCnt,arcDir);              
      printf("tanstlatexp %16.12f tanstlonexp %16.12f\n",tanstlatexp*RAD2DEG,tanstlonexp*RAD2DEG);
      printf("tanendlatexp %16.12f tanendlonexp %16.12f\n",tanendexp.latitude*RAD2DEG,tanendexp.longitude*RAD2DEG);
                                  printf("tanAngle %16.12f intAngle %16.12f angle %16.12f\n",tangencyAngle*RAD2DEG,interiorAngle*RAD2DEG,angle);
					printf("\n%s Passed\n", testName)*/;
					passedCount++;
                                 }

				testCaseCount++;
        	}
        }
    }
    } //while

    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("Finish %s\n", set.name.c_str());

    return set;

}

/*
 * NAME: testGeoTanToArcAtAngleToGeo_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the geoTanToArcAtAngleToGeo function.
 *
 * 		This function runs all the test cases and returns the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testGeoTanToArcAtAngleToGeo_AllSets(TestSuite) - A test suite with the following metrics:
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
TestSuite testGeoTanToArcAtAngleToGeo_AllSets()
{
	TestSuite suite;
	TestSet set1;
	TestSet set2;

    printf("\nStart testGeoTanToArcAtAngleToGeo_AllSets\n");

    suite = newTestSuite("testGeoTanToArcAtAngleToGeo_AllSets");

    set1 = testGeoTanToArcAtAngleToGeo_Set1();
    addTestSet(set1,&suite);

    set2 = testGeoTanToArcAtAngleToGeo_Set2();
    addTestSet(set2,&suite);

    displayTestSuite(suite);

    printf("Finish testGeoTanToArcAtAngleToGeo_AllSets\n\n\n");

    return suite;
}

/* NAME: testArcTanToArcAndGeo_Set1
 *
 * DESCRIPTION:
 *      Runs tests to determine the tangent arc given a geodesic, arc, and radius.
 *      for the arcTanToArcAndGeo function.
 *
 *      Original set of test data created by Rich Snow, testing implementation
 *      done by Rich Snow.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testArcTanToArcAndGeo_Set1(TestSet) - A test set with the following metrics:
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
TestSet testArcTanToArcAndGeo_Set1()
{

    double DEG2RAD = M_PI / 180.0;
    const double RAD2DEG = 180.0 / M_PI;
    double ANGLE_TOL = 5.0e-8;//1.0e-5;//degrees

    char testName[80];

    ArcDirection arcDir;
    LineType lineType;


    LLPoint arcCenter, arcStart, arcEnd;
    LLPoint arcCenterexp, arcStartexp, arcEndexp;
    Arc arc;
    LLPoint geoStart, geoEnd;
    LLPoint tempLLPoint;
    Geodesic geo;
    double radius, startaz, r2;
    double temp, temp1, tempcrs; 
    double latC, lonC, r1, az1, az2, crs1, crs2;

    Arc tangentArc;
    int tangentArcLocation;
    int tanArcOnGeo = 0;
    int tanArcOnArc = 0;
    LLPointPair tanArcAndArcIntersections;
    int tanArcAndArcIntersectionCnt = -1;
    LLPointPair tanArcAndGeoIntersections;
    int tanArcAndGeoIntersectionCnt = -1;
    double tangencyAngleAtGeo;
    double crsFromTanArcEndToGeoEnd, crsFromGeoEndToTanArcEnd;
    int outputMatlab = 0;//0 = don't output matlab code to the console
    int tempint;


    int fail;
    int testNum = 0;

    int passedCount = 0, failedCount = 0, errorCount = 0;
    int testCaseCount = 0;
    int setupFailureCount = 0;
    int unverifiedCount = 0;
    ErrorSet err = 0;
    long newSeed = 20080523;

    TestSet set;

    set = newTestSet("testArcTanToArcAndGeo_Set1");

    printf("Start %s\n", set.name.c_str());

    srand(newSeed);  //Initialize the random number generator

    while (testNum < 1000)
    {
      err = 0;
      fail = 0;

      //generate input data for Arc
      latC = randLat();
      lonC = randLon();
      r1 = (double) ((rand() % 100) + 1);
      r2 = r1 / 5.0;
      az1 = (double) ((rand() % 45) + 1); //degrees (0,90]
      az2 = (double) ((rand() % 45) + 1); //degrees (0,90]
      testNum++;
      sprintf(testName, "TEST%-d",testNum);

      tempint = (rand() % 4);

      if ((latC > 89.0) && (tempint == 3 || tempint == 1))
        latC = latC - 1.0;
      if ((latC < -89.0) && (tempint == 0 || tempint == 2))
        latC = latC + 1.0;
      arcCenter.latitude = latC * DEG2RAD;
      arcCenter.longitude = lonC * DEG2RAD;
      crs1 = az1 * DEG2RAD;
      crs2 = az2 * DEG2RAD;


      if ((latC > 87.9) && (tempint == 3 || tempint == 1))
      {
        r1 = r1 / 2.0;
        r2 = r1 / 5.0;
      }
      if ((latC < -88.0) && (tempint == 0 || tempint == 2))
      {
        r1 = r1 / 2.0;
        r2 = r1 / 5.0;
      }

      if (tempint < 2)
      {
        arcDir = CLOCKWISE;
        if (tempint == 0)
          startaz = 105.0 * DEG2RAD;
        else
          startaz = 315.0 * DEG2RAD;
        err |= direct(arcCenter, startaz - crs1, r1, &arcStart, EPS);
        err |= direct(arcCenter, startaz + crs2, r1, &arcEnd, EPS);
      }
      else  
      {
        arcDir = COUNTERCLOCKWISE;
        if (tempint == 2)
          startaz = 255.0 * DEG2RAD;
        else
          startaz = 45.0 * DEG2RAD;
        err |= direct(arcCenter, startaz + crs1, r1, &arcStart, EPS);
        err |= direct(arcCenter, startaz - crs2, r1, &arcEnd, EPS);
      } 
      err |= createArc(&arc, arcCenter, arcStart, arcEnd, arcDir, TOL, EPS);

      temp = (double) ((rand() % 21) - 10);
      if (latC > 88.0)
      {
        if (tempint == 1)
          temp = -20.0;
        else if (tempint == 3)
          temp = 20.0;
      }
      if (latC < -88.0)
      {
        if (tempint == 0)
          temp = -20.0;
        else if (tempint == 2)
          temp = 20.0;
      }
      if (tempint == 0)
        err |= direct(arcCenter, (225.0+temp)*DEG2RAD, r1, &arcStartexp, EPS);
      else if (tempint == 1)
        err |= direct(arcCenter, (45.0+temp)*DEG2RAD, r1, &arcStartexp, EPS);
      else if (tempint == 2)
        err |= direct(arcCenter, (135.0+temp)*DEG2RAD, r1, &arcStartexp, EPS);
      else
        err |= direct(arcCenter, (315.0+temp)*DEG2RAD, r1, &arcStartexp, EPS);
      err |= invCrs(arcStartexp, arcCenter, &tempcrs, &temp, EPS);
      err |= direct(arcStartexp, tempcrs + M_PI, r2, &arcCenterexp, EPS);
        temp = (double) ((rand() % 11) - 5);
      if (tempint == 0 || tempint == 2) 
        err |= direct(arcCenterexp, temp*DEG2RAD, r2, &arcEndexp, EPS);
      else
        err |= direct(arcCenterexp, temp*DEG2RAD + M_PI, r2, &arcEndexp, EPS);
      err |= invCrs(arcEndexp, arcCenterexp, &tempcrs, &temp, EPS);
      if (tempint < 2)
      {
        err |= direct(arcEndexp, tempcrs - M_PI_2, 2.0 * r1, &geoStart, EPS);
        err |= direct(arcEndexp, tempcrs + M_PI_2, 1.0, &tempLLPoint, EPS);
      }
      else 
      {
        err |= direct(arcEndexp, tempcrs + M_PI_2, 2.0 * r1, &geoStart, EPS);
        err |= direct(arcEndexp, tempcrs - M_PI_2, 1.0, &tempLLPoint, EPS);
      }
  
      //geoEnd has been moved 1 nmi so geoEnd and arcEndexp are not the same point
      geoEnd.latitude = tempLLPoint.latitude;
      geoEnd.longitude = tempLLPoint.longitude;


      lineType = SEGMENT;

	err |= createGeo(&geo, geoStart, geoEnd, lineType, EPS);


        radius = r2;

	if (outputMatlab) displayMatlabGeo(geo, "geo", 0);
	if (outputMatlab) displayMatlabArc(arc, "arc", 0);

	if (getMaskedError(err, getMaskAll()))
	{
		printf("\nError occurred in pre-arcTanToArcAndGeo err=0x%lx\n", err);
		errorCount++;
		setupFailureCount++;
		testCaseCount++;
	} else {

        //RUN TEST
        err = arcTanToArcAndGeo(arc, geo, radius, &tangentArc, &tangentArcLocation, TOL, EPS);
        //The following block is used when checking that the test setup is correct. Also, the line above is commented out when doing that.
        /*if (arcDir == CLOCKWISE)
          err |= createArc(&tangentArc, arcCenterexp, arcStartexp, arcEndexp, COUNTERCLOCKWISE, TOL, EPS);
        else
          err |= createArc(&tangentArc, arcCenterexp, arcStartexp, arcEndexp, CLOCKWISE, TOL, EPS);*/

        if (outputMatlab) displayMatlabArc(tangentArc, "tangentArc", 0);

        if (getMaskedError(err, getMaskAll()))
        {
            printf("\nError occurred in arcTanToArcAndGeo err=0x%lx\n", err);
            errorCount++;
            failedCount++;
            testCaseCount++;
        } else {

        	err = 0;

        	//treat the input geo as infinite
        	tanArcOnGeo = ptIsOnGeo(geo.startPoint, geo.endPoint, tangentArc.endPoint, INFINITE, &err, TESTTOL, EPS);

        	//treat the input arc as a circle
        	tanArcOnArc = ptIsOnArc(arc.centerPoint, arc.radius, arc.startAz, arc.startAz, arc.dir, tangentArc.startPoint, &err, TESTTOL, EPS);

        	//find the intersections of the tangent arc and arc
        	err |= arcIntx(arc.centerPoint, arc.radius, tangentArc.centerPoint, tangentArc.radius, tanArcAndArcIntersections, &tanArcAndArcIntersectionCnt, TOL, EPS);

        	//find the intersections of the tangent arc and geo
        	err |= geoArcIntx(geo.startPoint, geo.startAz, tangentArc.centerPoint, tangentArc.radius, tanArcAndGeoIntersections, &tanArcAndGeoIntersectionCnt, TOL, EPS);

        	//find the angle at the tangency of the arc with the geodesic
        	err |= invCrs(tangentArc.endPoint, geo.endPoint, &crsFromTanArcEndToGeoEnd, &crsFromGeoEndToTanArcEnd, EPS);
        	err |= invCrs(tangentArc.endPoint, tangentArc.centerPoint, &temp, &temp1, EPS);
        	//err |= minSubtendedAngle(crsFromTanArcEndToGeoEnd, tangentArc.endAz, &tangencyAngleAtGeo);
        	err |= minSubtendedAngle(crsFromTanArcEndToGeoEnd, temp, &tangencyAngleAtGeo);

        	//check if errors occurred while testing result
        	if (getMaskedError(err, getMaskAll()))
        	{
        		printf("\nError occurred in pre-arcTanToArcAndGeo err=0x%lx\n", err);
        		errorCount++;
        		setupFailureCount++;
        		testCaseCount++;
        	} else {

        		//test that tangent arc lies on input geo
        		if (!tanArcOnGeo)
        			fail = 1;

        		//test that tangent arc lies on arc
        		if (!tanArcOnArc)
        			fail = 1;

        		//test tangency of tangent arc and arc with respect to intersection points
        		if (tanArcAndArcIntersectionCnt == 0)
        		{
        			fail = 1;
        		}
        		else if ( tanArcAndArcIntersectionCnt == 2 && (!ptsAreSame(tanArcAndArcIntersections[0], tanArcAndArcIntersections[1], TESTTOL)) )
        		{
					fail = 1;
        		}

        		//test tangency of tangent arc and geo with respect to intersection points
        		if (tanArcAndGeoIntersectionCnt == 0)
        		{
        			fail = 1;
        		}
        		else if ( tanArcAndGeoIntersectionCnt == 2 && (!ptsAreSame(tanArcAndGeoIntersections[0], tanArcAndGeoIntersections[1], TESTTOL)) )
        		{
					fail = 1;
        		}

        		//test tangency of arc with respect to angle formed with the geo
        		if (fabs(fabs(tangencyAngleAtGeo*RAD2DEG) - M_PI_2*RAD2DEG) > ANGLE_TOL)
        		{
        			fail = 1;
        		}

                        if (!ptsAreSame(tangentArc.startPoint, arcStartexp, TESTTOL) || !ptsAreSame(tangentArc.endPoint, arcEndexp, TESTTOL) || !ptsAreSame(tangentArc.centerPoint, arcCenterexp, TESTTOL))
                        {
                          fail = 1;
                        }

                        if (fabs(tangentArc.radius - r2) >= TESTTOL)
                        {
                          fail = 1;
                        }

				if (fail)
				{
                                  printf("%s tanArcOnGeo %d tanArcOnArc %d ArcArcIntCnt %d ArcGeoIntCnt %d arcDir %d\n",testName,tanArcOnGeo,tanArcOnArc,tanArcAndArcIntersectionCnt,tanArcAndGeoIntersectionCnt,arcDir);              
                                  printf("startaz %16.12f crs1 %16.12f crs2 %16.12f\n",startaz*RAD2DEG,az1,az2);
                                  printf("geo.startAz %10.6f geometry %d\n",geo.startAz*RAD2DEG, tempint);
                                  printf("tanAngle %16.12f\n",tangencyAngleAtGeo*RAD2DEG);
					printf("\n%s Failed\n", testName);
					failedCount++;
				} else {
                                  //printf("%s tanArcOnGeo %d tanArcOnArc %d ArcArcIntCnt %d ArcGeoIntCnt %d arcDir %d\n",testName,tanArcOnGeo,tanArcOnArc,tanArcAndArcIntersectionCnt,tanArcAndGeoIntersectionCnt,arcDir);              
                                  //printf("startaz %16.12f crs1 %16.12f crs2 %16.12f\n",startaz*RAD2DEG,az1,az2);
                                  //printf("geo.startAz %10.6f geometry %d\n",geo.startAz*RAD2DEG, tempint);
                                  //printf("tanAngle %16.12f\n",tangencyAngleAtGeo*RAD2DEG);
					passedCount++;
                                }

				testCaseCount++;
        	}
        }
    }
    } //while

    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("Finish %s\n", set.name.c_str());

    return set;

}

/*
 * NAME: testArcTanToArcAndGeo_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the arcTanToArcAndGeo function.
 *
 * 		This function runs all the test cases and returns the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testArcTanToArcAndGeo_AllSets(TestSuite) - A test suite with the following metrics:
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
TestSuite testArcTanToArcAndGeo_AllSets()
{
	TestSuite suite;
	TestSet set1;

    suite = newTestSuite("testArcTanToArcAndGeo_AllSets");

    printf("\nStart %s\n", suite.name.c_str());

    set1 = testArcTanToArcAndGeo_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("Finish %s\n\n\n", suite.name.c_str());

    return suite;
}

/* NAME: testArcFromStartAndEnd_Set1
 *
 * DESCRIPTION:
 *      Runs tests to create an arc given a start point, start course, and end
 *      point for the arcFromStartAndEnd function.
 *
 *      Original set of test data created by Rich Snow, testing implementation
 *      done by Rich Snow.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testArcFromStartAndEnd_Set1(TestSet) - A test set with the following metrics:
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
TestSet testArcFromStartAndEnd_Set1()
{

    double DEG2RAD = M_PI / 180.0;
    const double RAD2DEG = 180.0 / M_PI;
    double ANGLE_TOL = 5.0e-8;//1.0e-5;//degrees

    char testName[80];

    ArcDirection arcDirexp;
    LLPoint arcCenterexp, arcStartexp, arcEndexp;
    Arc arc, arcexp;
    double startaz, startcrs, startcrsgen;
    double temp, tempcrs; 
    double latC, lonC, r1, az1, az2, crs1, crs2;

    int startPtOnArc = 0;
    int endPtOnArc = 0;
    int outputMatlab = 0;//0 = don't output matlab code to the console


    int fail;
    int testNum = 0;

    int passedCount = 0, failedCount = 0, errorCount = 0;
    int testCaseCount = 0;
    int setupFailureCount = 0;
    int unverifiedCount = 0;
    ErrorSet err = 0;
    long newSeed = 20080523;

    TestSet set;

    set = newTestSet("testArcFromStartAndEnd_Set1");

    printf("Start %s\n", set.name.c_str());

    srand(newSeed);  //Initialize the random number generator

    while (testNum < 100)
    {
      err = 0;
      fail = 0;

      //generate the arc (arcexp) to be compared with arc constructed by this function 
      latC = randLat();
      lonC = randLon();
      r1 = (double) ((rand() % 100) + 1);
      az1 = (double) ((rand() % 180) + 1); //degrees (0,180]
      az2 = (double) ((rand() % 180) + 1); //degrees (0,180]
      startaz = DEG2RAD*randAzimuth();
      testNum++;
      sprintf(testName, "TEST%-d",testNum);

      arcCenterexp.latitude = latC * DEG2RAD;
      arcCenterexp.longitude = lonC * DEG2RAD;
      crs1 = az1 * DEG2RAD;
      crs2 = az2 * DEG2RAD;

      if ((rand() % 2) == 0)
        arcDirexp = CLOCKWISE;
      else
        arcDirexp = COUNTERCLOCKWISE;

      if (arcDirexp == CLOCKWISE)
      {
        err |= direct(arcCenterexp, startaz - crs1, r1, &arcStartexp, EPS);
        err |= direct(arcCenterexp, startaz + crs2, r1, &arcEndexp, EPS);
        err |= invCrs(arcStartexp, arcCenterexp, &tempcrs, &temp, EPS);
        startcrs = modcrs(tempcrs - M_PI_2);
      }
      else
      {
        err |= direct(arcCenterexp, startaz + crs1, r1, &arcStartexp, EPS);
        err |= direct(arcCenterexp, startaz - crs2, r1, &arcEndexp, EPS);
        err |= invCrs(arcStartexp, arcCenterexp, &tempcrs, &temp, EPS);
        startcrs = modcrs(tempcrs + M_PI_2);
      }

      err |= createArc(&arcexp, arcCenterexp, arcStartexp, arcEndexp, arcDirexp, TOL, EPS);
	if (outputMatlab) displayMatlabArc(arcexp, "arcexp", 0);

	if (getMaskedError(err, getMaskAll()))
	{
		printf("\nError occurred in pre-arcFromStartAndEnd err=0x%lx\n", err);
		errorCount++;
		setupFailureCount++;
		testCaseCount++;
	} else {

        //RUN TEST
        err = arcFromStartAndEnd(arcStartexp, startcrs, arcEndexp, &arc, TOL, EPS);

        if (outputMatlab) displayMatlabArc(arc, "arc", 0);

        if (getMaskedError(err, getMaskAll()))
        {
            printf("\nError occurred in arcFromStartAndEnd err=0x%lx\n", err);
            errorCount++;
            failedCount++;
            testCaseCount++;
        } else {

        	err = 0;

        	startPtOnArc = ptIsOnArc(arc.centerPoint, arc.radius, arc.startAz, arc.startAz, arc.dir, arcStartexp, &err, TESTTOL, EPS);

        	endPtOnArc = ptIsOnArc(arc.centerPoint, arc.radius, arc.startAz, arc.startAz, arc.dir, arcEndexp, &err, TESTTOL, EPS);

                err |= invCrs(arc.startPoint, arc.centerPoint, &startcrsgen, &temp, EPS);
                if (arc.dir == CLOCKWISE)
                  startcrsgen = modcrs(startcrsgen - M_PI_2);
                else
                  startcrsgen = modcrs(startcrsgen + M_PI_2);

        	//check if errors occurred while testing result
        	if (getMaskedError(err, getMaskAll()))
        	{
        		printf("\nError occurred in post-arcFromStartAndEnd err=0x%lx\n", err);

        		errorCount++;
        		setupFailureCount++;
        		testCaseCount++;
        	} else {

        		//test that start and end points lie on arc
        		if (!startPtOnArc || !endPtOnArc)
        			fail = 1;

        		//test that the arc has the right radius
        		if (fabs(arc.radius - r1) >= TESTTOL)
        			fail = 1;


        		if ( !ptsAreSame(arc.centerPoint, arcCenterexp, TESTTOL)) 
        		{
					fail = 1;
        		}

        		if ( !ptsAreSame(arc.startPoint, arcStartexp, TESTTOL)) 
        		{
					fail = 1;
        		}

        		if ( !ptsAreSame(arc.endPoint, arcEndexp, TESTTOL)) 
        		{
					fail = 1;
        		}

                        if (arcDirexp != arc.dir)
        		{
					fail = 1;
        		}

                        if (fabs(startcrs - startcrsgen)*RAD2DEG > ANGLE_TOL)
                        {
                          fail = 1;
                        }


				if (fail)
				{
                                  printf("%s startPtOnArc %d endPtOnArc %d arcDirexp %d arcDircalc %d\n",testName,startPtOnArc,endPtOnArc,arcDirexp,arc.dir);              
                                  printf("startaz %16.12f crs1 %16.12f crs2 %16.12f\n",startaz*RAD2DEG,az1,az2);
                                  printf("startcrs %16.12f tempcrs %16.12f\n",startcrs*RAD2DEG,tempcrs*RAD2DEG);
					printf("\n%s Failed \n", testName);
					failedCount++;
				} else
					passedCount++;

				testCaseCount++;
        	}
        }
    }
    } //while

    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("Finish %s\n", set.name.c_str());

    return set;

}

/*
 * NAME: testArcFromStartAndEnd_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the arcFromStartAndEnd function.
 *
 * 		This function runs all the test cases and returns the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testArcFromStartAndEnd_AllSets(TestSuite) - A test suite with the following metrics:
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
TestSuite testArcFromStartAndEnd_AllSets()
{
	TestSuite suite;
	TestSet set1;

    suite = newTestSuite("testArcFromStartAndEnd_AllSets");

    printf("\nStart %s\n", suite.name.c_str());

    set1 = testArcFromStartAndEnd_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("Finish %s\n\n\n", suite.name.c_str());

    return suite;
}

/* NAME: testArcTanToCrs_Set1
 *
 * DESCRIPTION:
 *      Runs tests to create an arc given a start point, tangent start course, 
 *      end point, and end course for the 
 *      arcTanToCrs function.
 *
 *      Original set of test data created by Rich Snow, testing implementation
 *      done by Rich Snow.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testArcTanToCrs_Set1(TestSet) - 
 *              A test set with the following metrics:
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
TestSet testArcTanToCrs_Set1()
{

    double DEG2RAD = M_PI / 180.0;
    const double RAD2DEG = 180.0 / M_PI;
    double ANGLE_TOL = 5.0e-8;//1.0e-5;//degrees

    char testName[80];

    ArcDirection arcDirexp;
    LLPoint arcCenterexp, arcStartexp, arcEndexp;
    Arc arc, arcexp;
    double startaz, startcrs, outcrs, startcrsgen, outcrsgen;
    double temp, tempcrs; 
    double latC, lonC, r1, az1, az2, crs1, crs2;

    int startPtOnArc = 0;
    int endPtOnArc = 0;
    int outputMatlab = 0;//0 = don't output matlab code to the console


    int fail;
    int testNum = 0;

    int passedCount = 0, failedCount = 0, errorCount = 0;
    int testCaseCount = 0;
    int setupFailureCount = 0;
    int unverifiedCount = 0;
    ErrorSet err = 0;
    long newSeed = 20080523;

    TestSet set;

    set = newTestSet("testArcTanToCrs_Set1");

    printf("Start %s\n", set.name.c_str());

    srand(newSeed);  //Initialize the random number generator

    while (testNum < 100)
    {
      err = 0;
      fail = 0;

      //generate the arc (arcexp) to be compared with arc constructed by this function 
      latC = randLat();
      lonC = randLon();
      r1 = (double) ((rand() % 100) + 1);
      az1 = (double) ((rand() % 180) + 1); //degrees (0,180]
      az2 = (double) ((rand() % 180) + 1); //degrees (0,180]
      startaz = DEG2RAD*randAzimuth();
      testNum++;
      sprintf(testName, "TEST%-d",testNum);

      arcCenterexp.latitude = latC * DEG2RAD;
      arcCenterexp.longitude = lonC * DEG2RAD;
      crs1 = az1 * DEG2RAD;
      crs2 = az2 * DEG2RAD;

      if ((rand() % 2) == 0)
        arcDirexp = CLOCKWISE;
      else
        arcDirexp = COUNTERCLOCKWISE;

      if (arcDirexp == CLOCKWISE)
      {
        err |= direct(arcCenterexp, startaz - crs1, r1, &arcStartexp, EPS);
        err |= direct(arcCenterexp, startaz + crs2, r1, &arcEndexp, EPS);
        err |= invCrs(arcStartexp, arcCenterexp, &tempcrs, &temp, EPS);
        startcrs = modcrs(tempcrs - M_PI_2);
        err |= invCrs(arcEndexp, arcCenterexp, &tempcrs, &temp, EPS);
        outcrs = modcrs(tempcrs - M_PI_2);
      }
      else
      {
        err |= direct(arcCenterexp, startaz + crs1, r1, &arcStartexp, EPS);
        err |= direct(arcCenterexp, startaz - crs2, r1, &arcEndexp, EPS);
        err |= invCrs(arcStartexp, arcCenterexp, &tempcrs, &temp, EPS);
        startcrs = modcrs(tempcrs + M_PI_2);
        err |= invCrs(arcEndexp, arcCenterexp, &tempcrs, &temp, EPS);
        outcrs = modcrs(tempcrs + M_PI_2);
      }

      err |= createArc(&arcexp, arcCenterexp, arcStartexp, arcEndexp, arcDirexp, TOL, EPS);
	if (outputMatlab) displayMatlabArc(arcexp, "arcexp", 0);

	if (getMaskedError(err, getMaskAll()))
	{
		printf("\nError occurred in pre-arcTanToCrs err=0x%lx\n", err);
		errorCount++;
		setupFailureCount++;
		testCaseCount++;
	} else {

        //RUN TEST
        err = arcTanToCrs(arcStartexp, startcrs, arcEndexp, outcrs, &arc, TOL, EPS);

        if (outputMatlab) displayMatlabArc(arc, "arc", 0);

        if (getMaskedError(err, getMaskAll()))
        {
            printf("\nError occurred in arcTanToCrs err=0x%lx\n", err);
            errorCount++;
            failedCount++;
            testCaseCount++;
        } else {

        	err = 0;

        	startPtOnArc = ptIsOnArc(arc.centerPoint, arc.radius, arc.startAz, arc.startAz, arc.dir, arcStartexp, &err, TESTTOL, EPS);

        	endPtOnArc = ptIsOnArc(arc.centerPoint, arc.radius, arc.startAz, arc.startAz, arc.dir, arcEndexp, &err, TESTTOL, EPS);

                err |= invCrs(arc.startPoint, arc.centerPoint, &startcrsgen, &temp, EPS);
                err |= invCrs(arc.endPoint, arc.centerPoint, &outcrsgen, &temp, EPS);
                if (arc.dir == CLOCKWISE)
                {
                  startcrsgen = modcrs(startcrsgen - M_PI_2);
                  outcrsgen = modcrs(outcrsgen - M_PI_2);
                }
                else
                {
                  startcrsgen = modcrs(startcrsgen + M_PI_2);
                  outcrsgen = modcrs(outcrsgen + M_PI_2);
                }

        	//check if errors occurred while testing result
        	if (getMaskedError(err, getMaskAll()))
        	{
        		printf("\nError occurred in post-arcTanToCrs err=0x%lx\n", err);
        		errorCount++;
        		setupFailureCount++;
        		testCaseCount++;
        	} else {

        		//test that start and end points lie on arc
        		if (!startPtOnArc || !endPtOnArc)
        			fail = 1;

        		//test that the arc has the right radius
        		if (fabs(arc.radius - r1) >= TESTTOL)
        			fail = 1;


        		if ( !ptsAreSame(arc.centerPoint, arcCenterexp, TESTTOL)) 
        		{
					fail = 1;
        		}

        		if ( !ptsAreSame(arc.startPoint, arcStartexp, TESTTOL)) 
        		{
					fail = 1;
        		}

        		if ( !ptsAreSame(arc.endPoint, arcEndexp, TESTTOL)) 
        		{
					fail = 1;
        		}

                        if (arcDirexp != arc.dir)
        		{
					fail = 1;
        		}

                        if (fabs(startcrs - startcrsgen)*RAD2DEG > ANGLE_TOL)
                        {
                          fail = 1;
                        }

                        if (fabs(outcrs - outcrsgen)*RAD2DEG > ANGLE_TOL)
                        {
                          fail = 1;
                        }


				if (fail)
				{
                                  printf("%s startPtOnArc %d endPtOnArc %d arcDirexp %d arcDircalc %d\n",testName,startPtOnArc,endPtOnArc,arcDirexp,arc.dir);              
                                  printf("startaz %16.12f crs1 %16.12f crs2 %16.12f\n",startaz*RAD2DEG,az1,az2);
                                  printf("startcrs %16.12f tempcrs %16.12f\n",startcrs*RAD2DEG,tempcrs*RAD2DEG);
					printf("\n%s Failed \n", testName);
					failedCount++;
				} else
					passedCount++;

				testCaseCount++;
        	}
        }
    }
    } //while

    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("Finish %s\n", set.name.c_str());

    return set;

}

/*
 * NAME: testArcTanToCrs_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the arcTanToCrs function.
 *
 * 		This function runs all the test cases and returns the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		arcTanToCrs(TestSuite) - A test suite with the following metrics:
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
TestSuite testArcTanToCrs_AllSets()
{
	TestSuite suite;
	TestSet set1;

    suite = newTestSuite("testArcTanToCrs_AllSets");

    printf("\nStart %s\n", suite.name.c_str());

    set1 = testArcTanToCrs_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("Finish %s\n\n\n", suite.name.c_str());

    return suite;
}

/* NAME: testArcEndFromStartAndRadius_Set1
 *
 * DESCRIPTION:
 *      Runs tests to create an arc given a start point, tangent start course, 
 *      and radius for the 
 *      arcEndFromStartAndRadius function.
 *
 *      Original set of test data created by Rich Snow, testing implementation
 *      done by Rich Snow.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testArcEndFromStartAndRadius_Set1(TestSet) - 
 *              A test set with the following metrics:
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
TestSet testArcEndFromStartAndRadius_Set1()
{

    double DEG2RAD = M_PI / 180.0;
    const double RAD2DEG = 180.0 / M_PI;
    double ANGLE_TOL = 5.0e-8;//1.0e-5;//degrees

    char testName[80];

    ArcDirection arcDirexp;
    LLPoint arcCenterexp, arcStartexp, arcEndexp, nextEndPoint;
    Arc arc, arcexp;
    double startaz, startcrs, outcrs, startcrsgen, outcrsgen;
    double temp, tempcrs; 
    double latC, lonC, r1, az1, az2, crs1, crs2;

    int startPtOnArc = 0;
    int endPtOnArc = 0;
    int outputMatlab = 0;//0 = don't output matlab code to the console


    int fail;
    int testNum = 0;

    int passedCount = 0, failedCount = 0, errorCount = 0;
    int testCaseCount = 0;
    int setupFailureCount = 0;
    int unverifiedCount = 0;
    ErrorSet err = 0;
    long newSeed = 20080523;

    TestSet set;

    set = newTestSet("testArcEndFromStartAndRadius_Set1");

    printf("Start %s\n", set.name.c_str());

    srand(newSeed);  //Initialize the random number generator

    while (testNum < 500)
    {
      err = 0;
      fail = 0;

      //generate the arc (arcexp) to be compared with arc constructed by this function 
      latC = randLat();
      lonC = randLon();
      r1 = (double) ((rand() % 100) + 1);
      az1 = (double) ((rand() % 180) + 1); //degrees (0,180]
      az2 = (double) ((rand() % 180) + 1); //degrees (0,180]
      startaz = DEG2RAD*randAzimuth();
      testNum++;
      sprintf(testName, "TEST%-d",testNum);

      arcCenterexp.latitude = latC * DEG2RAD;
      arcCenterexp.longitude = lonC * DEG2RAD;
      crs1 = az1 * DEG2RAD;
      crs2 = az2 * DEG2RAD;

      if ((rand() % 2) == 0)
        arcDirexp = CLOCKWISE;
      else
        arcDirexp = COUNTERCLOCKWISE;

      if (arcDirexp == CLOCKWISE)
      {
        err |= direct(arcCenterexp, startaz - crs1, r1, &arcStartexp, EPS);
        err |= direct(arcCenterexp, startaz + crs2, r1, &arcEndexp, EPS);
        err |= invCrs(arcStartexp, arcCenterexp, &tempcrs, &temp, EPS);
        startcrs = modcrs(tempcrs - M_PI_2);
        err |= invCrs(arcEndexp, arcCenterexp, &tempcrs, &temp, EPS);
        outcrs = modcrs(tempcrs - M_PI_2);
        err |= direct(arcEndexp, outcrs, 2.0 * r1, &nextEndPoint, EPS);
      }
      else
      {
        err |= direct(arcCenterexp, startaz + crs1, r1, &arcStartexp, EPS);
        err |= direct(arcCenterexp, startaz - crs2, r1, &arcEndexp, EPS);
        err |= invCrs(arcStartexp, arcCenterexp, &tempcrs, &temp, EPS);
        startcrs = modcrs(tempcrs + M_PI_2);
        err |= invCrs(arcEndexp, arcCenterexp, &tempcrs, &temp, EPS);
        outcrs = modcrs(tempcrs + M_PI_2);
        err |= direct(arcEndexp, outcrs, 2.0 * r1, &nextEndPoint, EPS);
      }

      err |= createArc(&arcexp, arcCenterexp, arcStartexp, arcEndexp, arcDirexp, TOL, EPS);
	if (outputMatlab) displayMatlabArc(arcexp, "arcexp", 0);

	if (getMaskedError(err, getMaskAll()))
	{
		printf("\nError occurred in pre-arcEndFromStartAndRadius err=0x%lx\n", err);
		errorCount++;
		setupFailureCount++;
		testCaseCount++;
	} else {

        //RUN TEST
        err = arcEndFromStartAndRadius(arcStartexp, startcrs, r1, arcDirexp, nextEndPoint, &arc, TOL, EPS);

        if (outputMatlab) displayMatlabArc(arc, "arc", 0);

        if (getMaskedError(err, getMaskAll()))
        {
            printf("\nError occurred in arcEndFromStartAndRadius err=0x%lx\n", err);
            errorCount++;
            failedCount++;
            testCaseCount++;
        } else {

        	err = 0;

        	startPtOnArc = ptIsOnArc(arc.centerPoint, arc.radius, arc.startAz, arc.startAz, arc.dir, arcStartexp, &err, TESTTOL, EPS);

        	endPtOnArc = ptIsOnArc(arc.centerPoint, arc.radius, arc.startAz, arc.startAz, arc.dir, arcEndexp, &err, TESTTOL, EPS);

                err |= invCrs(arc.startPoint, arc.centerPoint, &startcrsgen, &temp, EPS);
                err |= invCrs(arc.endPoint, arc.centerPoint, &outcrsgen, &temp, EPS);
                if (arc.dir == CLOCKWISE)
                {
                  startcrsgen = modcrs(startcrsgen - M_PI_2);
                  outcrsgen = modcrs(outcrsgen - M_PI_2);
                }
                else
                {
                  startcrsgen = modcrs(startcrsgen + M_PI_2);
                  outcrsgen = modcrs(outcrsgen + M_PI_2);
                }

        	//check if errors occurred while testing result
        	if (getMaskedError(err, getMaskAll()))
        	{
        		printf("\nError occurred in post-arcEndFromStartAndRadius err=0x%lx\n", err);
        		errorCount++;
        		setupFailureCount++;
        		testCaseCount++;
        	} else {

        		//test that start and end points lie on arc
        		if (!startPtOnArc || !endPtOnArc)
        			fail = 1;

        		//test that the arc has the right radius
        		if (fabs(arc.radius - r1) >= TESTTOL)
        			fail = 1;


        		if ( !ptsAreSame(arc.centerPoint, arcCenterexp, TESTTOL)) 
        		{
					fail = 1;
        		}

        		if ( !ptsAreSame(arc.startPoint, arcStartexp, TESTTOL)) 
        		{
					fail = 1;
        		}

        		if ( !ptsAreSame(arc.endPoint, arcEndexp, TESTTOL)) 
        		{
					fail = 1;
        		}

                        if (arcDirexp != arc.dir)
        		{
					fail = 1;
        		}

                        if (fabs(startcrs - startcrsgen)*RAD2DEG > ANGLE_TOL)
                        {
                          fail = 1;
                        }

                        if (fabs(outcrs - outcrsgen)*RAD2DEG > ANGLE_TOL)
                        {
                          fail = 1;
                        }


				if (fail)
				{
                                  printf("%s startPtOnArc %d endPtOnArc %d arcDirexp %d arcDircalc %d\n",testName,startPtOnArc,endPtOnArc,arcDirexp,arc.dir);              
                                  printf("startaz %16.12f crs1 %16.12f crs2 %16.12f\n",startaz*RAD2DEG,az1,az2);
                                  printf("startcrs %16.12f tempcrs %16.12f\n",startcrs*RAD2DEG,tempcrs*RAD2DEG);
					printf("\n%s Failed \n", testName);
					failedCount++;
				} else {
                                  //printf("%s startPtOnArc %d endPtOnArc %d arcDirexp %d arcDircalc %d\n",testName,startPtOnArc,endPtOnArc,arcDirexp,arc.dir);              
                                 
					passedCount++;
                                }

				testCaseCount++;
        	}
        }
    }
    } //while

    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("Finish %s\n", set.name.c_str());

    return set;

}

/*
 * NAME: testArcEndFromStartAndRadius_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the arcEndFromStartAndRadius function.
 *
 * 		This function runs all the test cases and returns the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		arcEndFromStartAndRadius(TestSuite) - A test suite with the following metrics:
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
TestSuite testArcEndFromStartAndRadius_AllSets()
{
	TestSuite suite;
	TestSet set1;

    suite = newTestSuite("testArcEndFromStartAndRadius_AllSets");

    printf("\nStart %s\n", suite.name.c_str());

    set1 = testArcEndFromStartAndRadius_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("Finish %s\n\n\n", suite.name.c_str());

    return suite;
}

/* NAME: testArcEndFromStartAndCenter_Set1
 *
 * DESCRIPTION:
 *      Runs tests to create an arc given a start point, tangent start course, 
 *      and center for the 
 *      arcEndFromStartAndCenter function.
 *
 *      Original set of test data created by Rich Snow, testing implementation
 *      done by Rich Snow.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testArcEndFromStartAndCenter_Set1(TestSet) - 
 *              A test set with the following metrics:
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
TestSet testArcEndFromStartAndCenter_Set1()
{

    double DEG2RAD = M_PI / 180.0;
    const double RAD2DEG = 180.0 / M_PI;
    double ANGLE_TOL = 5.0e-8;//1.0e-5;//degrees

    char testName[80];

    ArcDirection arcDirexp;
    LLPoint arcCenterexp, arcStartexp, arcEndexp, nextEndPoint;
    Arc arc, arcexp;
    double startaz, startcrs, outcrs, startcrsgen, outcrsgen;
    double temp, tempcrs; 
    double latC, lonC, r1, az1, az2, crs1, crs2;

    int startPtOnArc = 0;
    int endPtOnArc = 0;
    int outputMatlab = 0;//0 = don't output matlab code to the console


    int fail;
    int testNum = 0;

    int passedCount = 0, failedCount = 0, errorCount = 0;
    int testCaseCount = 0;
    int setupFailureCount = 0;
    int unverifiedCount = 0;
    ErrorSet err = 0;
    long newSeed = 20080523;

    TestSet set;

    set = newTestSet("testArcEndFromStartAndCenter_Set1");

    printf("Start %s\n", set.name.c_str());

    srand(newSeed);  //Initialize the random number generator

    while (testNum < 500)
    {
      err = 0;
      fail = 0;

      //generate the arc (arcexp) to be compared with arc constructed by this function 
      latC = randLat();
      lonC = randLon();
      r1 = (double) ((rand() % 100) + 1);
      az1 = (double) ((rand() % 180) + 1); //degrees (0,180]
      az2 = (double) ((rand() % 180) + 1); //degrees (0,180]
      startaz = DEG2RAD*randAzimuth();
      testNum++;
      sprintf(testName, "TEST%-d",testNum);

      arcCenterexp.latitude = latC * DEG2RAD;
      arcCenterexp.longitude = lonC * DEG2RAD;
      crs1 = az1 * DEG2RAD;
      crs2 = az2 * DEG2RAD;

      if ((rand() % 2) == 0)
        arcDirexp = CLOCKWISE;
      else
        arcDirexp = COUNTERCLOCKWISE;

      if (arcDirexp == CLOCKWISE)
      {
        err |= direct(arcCenterexp, startaz - crs1, r1, &arcStartexp, EPS);
        err |= direct(arcCenterexp, startaz + crs2, r1, &arcEndexp, EPS);
        err |= invCrs(arcStartexp, arcCenterexp, &tempcrs, &temp, EPS);
        startcrs = modcrs(tempcrs - M_PI_2);
        err |= invCrs(arcEndexp, arcCenterexp, &tempcrs, &temp, EPS);
        outcrs = modcrs(tempcrs - M_PI_2);
        err |= direct(arcEndexp, outcrs, 2.0 * r1, &nextEndPoint, EPS);
      }
      else
      {
        err |= direct(arcCenterexp, startaz + crs1, r1, &arcStartexp, EPS);
        err |= direct(arcCenterexp, startaz - crs2, r1, &arcEndexp, EPS);
        err |= invCrs(arcStartexp, arcCenterexp, &tempcrs, &temp, EPS);
        startcrs = modcrs(tempcrs + M_PI_2);
        err |= invCrs(arcEndexp, arcCenterexp, &tempcrs, &temp, EPS);
        outcrs = modcrs(tempcrs + M_PI_2);
        err |= direct(arcEndexp, outcrs, 2.0 * r1, &nextEndPoint, EPS);
      }

      err |= createArc(&arcexp, arcCenterexp, arcStartexp, arcEndexp, arcDirexp, TOL, EPS);
	if (outputMatlab) displayMatlabArc(arcexp, "arcexp", 0);

	if (getMaskedError(err, getMaskAll()))
	{
		printf("\nError occurred in pre-arcEndFromStartAndCenter err=0x%lx\n", err);
		errorCount++;
		setupFailureCount++;
		testCaseCount++;
	} else {

        //RUN TEST
        err = arcEndFromStartAndCenter(arcStartexp, startcrs, arcCenterexp, arcDirexp, nextEndPoint, &arc, TOL, EPS);

        if (outputMatlab) displayMatlabArc(arc, "arc", 0);

        if (getMaskedError(err, getMaskAll()))
        {
            printf("\nError occurred in arcEndFromStartAndCenter err=0x%lx\n", err);
            errorCount++;
            failedCount++;
            testCaseCount++;
        } else {

        	err = 0;

        	startPtOnArc = ptIsOnArc(arc.centerPoint, arc.radius, arc.startAz, arc.startAz, arc.dir, arcStartexp, &err, TESTTOL, EPS);

        	endPtOnArc = ptIsOnArc(arc.centerPoint, arc.radius, arc.startAz, arc.startAz, arc.dir, arcEndexp, &err, TESTTOL, EPS);

                err |= invCrs(arc.startPoint, arc.centerPoint, &startcrsgen, &temp, EPS);
                err |= invCrs(arc.endPoint, arc.centerPoint, &outcrsgen, &temp, EPS);
                if (arc.dir == CLOCKWISE)
                {
                  startcrsgen = modcrs(startcrsgen - M_PI_2);
                  outcrsgen = modcrs(outcrsgen - M_PI_2);
                }
                else
                {
                  startcrsgen = modcrs(startcrsgen + M_PI_2);
                  outcrsgen = modcrs(outcrsgen + M_PI_2);
                }

        	//check if errors occurred while testing result
        	if (getMaskedError(err, getMaskAll()))
        	{
        		printf("\nError occurred in post-arcEndFromStartAndCenter err=0x%lx\n", err);
        		errorCount++;
        		setupFailureCount++;
        		testCaseCount++;
        	} else {

        		//test that start and end points lie on arc
        		if (!startPtOnArc || !endPtOnArc)
        			fail = 1;

        		//test that the arc has the right radius
        		if (fabs(arc.radius - r1) >= TESTTOL)
        			fail = 1;


        		if ( !ptsAreSame(arc.centerPoint, arcCenterexp, TESTTOL)) 
        		{
					fail = 1;
        		}

        		if ( !ptsAreSame(arc.startPoint, arcStartexp, TESTTOL)) 
        		{
					fail = 1;
        		}

        		if ( !ptsAreSame(arc.endPoint, arcEndexp, TESTTOL)) 
        		{
					fail = 1;
        		}

                        if (arcDirexp != arc.dir)
        		{
					fail = 1;
        		}

                        if (fabs(startcrs - startcrsgen)*RAD2DEG > ANGLE_TOL)
                        {
                          fail = 1;
                        }

                        if (fabs(outcrs - outcrsgen)*RAD2DEG > ANGLE_TOL)
                        {
                          fail = 1;
                        }


				if (fail)
				{
                                  printf("%s startPtOnArc %d endPtOnArc %d arcDirexp %d arcDircalc %d\n",testName,startPtOnArc,endPtOnArc,arcDirexp,arc.dir);              
                                  printf("startaz %16.12f crs1 %16.12f crs2 %16.12f\n",startaz*RAD2DEG,az1,az2);
                                  printf("startcrs %16.12f tempcrs %16.12f\n",startcrs*RAD2DEG,tempcrs*RAD2DEG);
					printf("\n%s Failed \n", testName);
					failedCount++;
				} else {
                                  //printf("%s startPtOnArc %d endPtOnArc %d arcDirexp %d arcDircalc %d\n",testName,startPtOnArc,endPtOnArc,arcDirexp,arc.dir);              
                                 
					passedCount++;
                                }

				testCaseCount++;
        	}
        }
    }
    } //while

    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("Finish %s\n", set.name.c_str());

    return set;

}

/*
 * NAME: testArcEndFromStartAndCenter_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the arcEndFromStartAndCenter function.
 *
 * 		This function runs all the test cases and returns the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		arcEndFromStartAndCenter(TestSuite) - A test suite with the following metrics:
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
TestSuite testArcEndFromStartAndCenter_AllSets()
{
	TestSuite suite;
	TestSet set1;

    suite = newTestSuite("testArcEndFromStartAndCenter_AllSets");

    printf("\nStart %s\n", suite.name.c_str());

    set1 = testArcEndFromStartAndCenter_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("Finish %s\n\n\n", suite.name.c_str());

    return suite;
}
} //namespace
