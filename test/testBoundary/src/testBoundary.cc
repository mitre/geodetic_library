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
 * NAME: testBndryCircleIntx_Set1
 *
 * DESCRIPTION:
 * 		This function is used to test the bndryCircleIntx function.
 *
 * 		This function runs the test data created by Rich Snow.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testBndryCircleIntx_Set1(TestSet) - A test set with the folling metrics:
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
TestSet testBndryCircleIntx_Set1()
{
  double DEG2RAD = M_PI / 180.0;
  int passedCount=0, failedCount=0;
  int errorCount = 0;
  int unverifiedCount = 0;
  int setupFailureCount = 0;
  int testCaseCount = 0;
  int outputMatlab = 0; //0 = don't output matlab code to the console
  ErrorSet err=0;
  double eps = 1.0e-20;
  double tol = 1.37e-9;
  TestSet set;

  LLPoint line1Start, line1End, line2Start, line2End;
  LLPoint locus1Start, locus1End, locus2Start, locus2End, locus3Start, locus3End;
  LLPoint arc1Center, arc1Start, arc1End, arc2Center, arc2Start, arc2End;
  double locus1Dist, locus2Dist, locus3Dist;
  Geodesic line1, line2;
  Locus locus1, locus2, locus3;
  Arc arc1, arc2;
  LLPoint obCenter, obPoint, geoStartExp, locusIntxExp;
  LLPointPair intx;
  
  Arc ob;
  Boundary b, cCommonBoundary;
  int i, j, n, temppass = 0;
  Shape* shape = NULL;
  ShapeType shapeClass;
  Arc* shapeArc = NULL;
  Geodesic* shapeGeo = NULL;
  Locus* shapeLocus = NULL;
  int shapeCount = -1;

  set = newTestSet("testBndryCircleIntx_Set1");
  printf("\n\nStart testBndryCircleIntx_Set1\n");

  //Bug 29391
  b = createBndry();
  err |= createPt(&line1Start, 37.364296535439436*DEG2RAD, -118.35575321139568*DEG2RAD);
  err |= createPt(&line1End, 37.365836493163414*DEG2RAD, -118.35378991460784*DEG2RAD);
  err |= createPt(&line2Start, 37.35507130404578*DEG2RAD, -118.34712186236166*DEG2RAD);
  err |= createPt(&line2End, 37.35878598896732*DEG2RAD, -118.34238551200936*DEG2RAD);
  err |= createPt(&locus1Start, 37.365066518425145*DEG2RAD, -118.35477157303242*DEG2RAD);
  err |= createPt(&locus1End, 37.35692867050427*DEG2RAD, -118.34475374553963*DEG2RAD);
  err |= createPt(&locus2Start, 37.365066518425145*DEG2RAD, -118.35477157303242*DEG2RAD);
  err |= createPt(&locus2End, 37.35692867050427*DEG2RAD, -118.34475374553963*DEG2RAD);
  err |= createPt(&obCenter, 37.36585*DEG2RAD, -118.35376111111111*DEG2RAD);
  err |= createPt(&obPoint, 37.36590492658283*DEG2RAD, -118.35376111111107*DEG2RAD);
  locus1Dist = 0.06583153347732182;
  locus2Dist = 0.15881570690152078;
  if (err) printf("\\nError(s) occurred during setup: 0x%lx", err);
  err = 0;
  err = createGeo(&line1, line1Start, line1End, SEGMENT, eps);
  err |= createGeo(&line2, line2Start, line2End, SEGMENT, eps);
  err |= createLocus(&locus1, locus1Start, locus1End, -locus1Dist, -locus2Dist, SEGMENT, tol, eps);
  err |= createLocus(&locus2, locus2Start, locus2End, locus1Dist, locus2Dist, SEGMENT, tol, eps);
  err |= addGeoToBndry(&b, &line1);
  err |= addGeoToBndry(&b, &line2);
  err |= addLocusToBndry(&b, &locus1);
  err |= addLocusToBndry(&b, &locus2);
  err |= createArc(&ob, obCenter, obPoint, obPoint, ArcDirection::CLOCKWISE, tol, eps);
  if (err) printf("\\nError(s) occurred during construction: 0x%lx", err);
        if (outputMatlab) displayMatlabGeo(line1, "line1", 0);
        if (outputMatlab) displayMatlabGeo(line2, "line2", 0);
        if (outputMatlab) displayMatlabLocus(locus1, "locus1", 0);
        if (outputMatlab) displayMatlabLocus(locus2, "locus2", 0);
        if (outputMatlab) displayMatlabArc(ob, "ob", 0);
  err |= geoArcIntx(line1Start, line1.startAz, obCenter, ob.radius, intx, &n, tol, eps);
  for (j = 0; j < n; j++)
  {
     if (ptIsOnGeo(line1Start, line1End, intx[j], SEGMENT, &err, tol, eps))
        geoStartExp = intx[j];
  }
  err |= locusArcIntx(locus1, obCenter, ob.radius, intx, &n, tol, eps);
  locusIntxExp = intx[0];
  cCommonBoundary = createBndry();
  err |= getMaskedError(bndryCircleIntx(b, ob, &cCommonBoundary, tol, eps), getMask(0, 1, 0, 0, 0, 0, 1));

  if (err)
  {
    errorCount++;
    failedCount++;
  }
  else
  {
     shapeCount = cCommonBoundary.length;
     for (i=0; i < shapeCount; i++) {
       //Get the shape
       shape = &cCommonBoundary.elements[i];
       //Get the shape type
       shapeClass = cCommonBoundary.elements[i].type;
       switch (shapeClass) {
       case ARC:
         shapeArc = (Arc*)shape;
         if (ptsAreSame(ob.centerPoint, shapeArc->centerPoint, TESTTOL) && ptsAreSame(shapeArc->startPoint, locusIntxExp, TESTTOL) && ptsAreSame(shapeArc->endPoint, geoStartExp, TESTTOL))
           temppass = 1;
         else
           temppass = 0; 
         break;
       case LOCUS:
         shapeLocus = (Locus*)shape;
         if (ptsAreSame(shapeLocus->locusStart, line1End, TESTTOL) && ptsAreSame(shapeLocus->locusEnd, locusIntxExp, TESTTOL))
           temppass = 1;
         else
           temppass = 0;
         break;
       case GEODESIC:
         shapeGeo = (Geodesic*)shape; 
         if (ptsAreSame(shapeGeo->startPoint, geoStartExp, TESTTOL) && ptsAreSame(shapeGeo->endPoint, line1End, TESTTOL))
           temppass = 1;
         else
           temppass = 0; 
         break;
         default:
           temppass = 0;
       }
       if (!temppass)
         break;
     }//for i
     if (temppass)
       passedCount++;
     else
       failedCount++;
  }
  testCaseCount++;

  //Bug 25861
   b = createBndry();
   err |= createPt(&line1Start, 33.993536492523305*DEG2RAD, -118.8194712971946*DEG2RAD);
   err |= createPt(&line1End, 33.92901255874419*DEG2RAD, -118.7826601300351*DEG2RAD);
   err |= createPt(&line2Start, 33.971446487765*DEG2RAD, -118.58932391049225*DEG2RAD);
   err |= createPt(&line2End, 33.90602519670383*DEG2RAD, -118.60544028274295*DEG2RAD);
   err |= createPt(&locus1Start, 33.96337929221117*DEG2RAD, -118.802258742832*DEG2RAD);
   err |= createPt(&locus1End, 33.97054652*DEG2RAD, -118.78416011*DEG2RAD);
   err |= createPt(&locus2Start, 33.97054652*DEG2RAD, -118.78416011*DEG2RAD);
   err |= createPt(&locus2End, 33.938736192963965*DEG2RAD, -118.5973851782385*DEG2RAD);
   err |= createPt(&locus3Start, 33.96081269363565*DEG2RAD, -118.72681588046969*DEG2RAD);
   err |= createPt(&locus3End, 33.938736192963965*DEG2RAD, -118.5973851782385*DEG2RAD);
   err |= createPt(&arc1Center, 33.81766053528569*DEG2RAD, -118.76180431286852*DEG2RAD);
   err |= createPt(&arc1Start, 33.92901255874419*DEG2RAD, -118.7826601300351*DEG2RAD);
   err |= createPt(&arc1End, 33.92809337072112*DEG2RAD, -118.734823550297*DEG2RAD);
   err |= createPt(&arc2Center, 33.97054651999997*DEG2RAD, -118.78416010999997*DEG2RAD);
   err |= createPt(&arc2Start, 34.000706215528545*DEG2RAD, -118.80136771950963*DEG2RAD);
   err |= createPt(&arc2End, 34.00326882417624*DEG2RAD, -118.7761673606076*DEG2RAD);
   err |= createPt(&obCenter, 34.00273028992501*DEG2RAD, -118.77629895394571*DEG2RAD);
   err |= createPt(&obPoint, 34.00327986312095*DEG2RAD, -118.77629895394571*DEG2RAD);
   locus1Dist = -2.0;
   locus2Dist = -2.0;
   locus3Dist = 2.0;
   if (err) printf("\\nError(s) occurred during setup: 0x%lx", err);
   err = 0;
   err |= createGeo(&line1, line1Start, line1End, SEGMENT, eps);
   err |= createGeo(&line2, line2Start, line2End, SEGMENT, eps);
   err |= createLocus(&locus1, locus1Start, locus1End, locus1Dist, locus1Dist, SEGMENT, tol, eps);
   err |= createLocus(&locus2, locus2Start, locus2End, locus2Dist, locus2Dist, SEGMENT, tol, eps);
   err |= createLocus(&locus3, locus3Start, locus3End, locus3Dist, locus3Dist, SEGMENT, tol, eps);
   err |= createArc(&arc1, arc1Center, arc1Start, arc1End, CLOCKWISE, tol, eps);
   err |= createArc(&arc2, arc2Center, arc2Start, arc2End, CLOCKWISE, tol, eps);
   err |= addGeoToBndry(&b, &line1);
   err |= addGeoToBndry(&b, &line2);
   err |= addLocusToBndry(&b, &locus1);
   err |= addLocusToBndry(&b, &locus2);
   err |= addLocusToBndry(&b, &locus3);
   err |= addArcToBndry(&b, &arc1);
   err |= addArcToBndry(&b, &arc2);
   err |= createArc(&ob, obCenter, obPoint, obPoint, COUNTERCLOCKWISE, tol, eps);
   if (err) printf("\\nError(s) occurred during construction: 0x%lx", err);
   cCommonBoundary = createBndry();
   err |= getMaskedError(bndryCircleIntx(b, ob, &cCommonBoundary, tol, eps), getMask(0, 1, 0, 0, 0, 0, 1));

   if (err)
   {
     errorCount++;
     failedCount++;
   }
   else
   {
     shapeCount = cCommonBoundary.length;
     if (shapeCount == 0)
       failedCount++;
     for (i=0; i < shapeCount; i++) {
       //Get the shape
       shape = &cCommonBoundary.elements[i];
       //Get the shape type
       shapeClass = cCommonBoundary.elements[i].type;
       if ((i > 1) || (shapeClass != ARC))
       {
         failedCount++;
         break;
       }
       else
       {
         shapeArc = (Arc*)shape;
         if (ptsAreSame(obCenter, shapeArc->centerPoint, TESTTOL) && ptsAreSame(obPoint, shapeArc->startPoint, TESTTOL) && ptsAreSame(obPoint, shapeArc->endPoint, TESTTOL))
           passedCount++;
         else
           failedCount++;
       }
     }//for i
   }
   testCaseCount++;

   set.testCases = testCaseCount;
   set.pass = passedCount;
   set.fail = failedCount;
   set.unverified = unverifiedCount;
   set.setupFailures = setupFailureCount;
   set.errors = errorCount;

   displayTestSet(set);

   printf("Finish testBndryCircleIntx_Set1\n\n\n");

   return set;
}

/*
 * NAME: testBndryCircleIntx_Set2
 *
 * DESCRIPTION:
 * 		This function is used to test the bndryCircleIntx function.
 *
 * 		This function runs the test data created by Rich Snow.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testBndryCircleIntx_Set2(TestSet) - A test set with the folling metrics:
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
TestSet testBndryCircleIntx_Set2()
{
    double DEG2RAD = M_PI / 180.0;
    //const double RAD2DEG = 180.0 / M_PI;
    int passedCount=0, failedCount=0;
    int errorCount = 0;
    int unverifiedCount = 0;
    int setupFailureCount = 0;
    int testCaseCount = 0;
    int outputMatlab = 0; //0 = don't output matlab code to the console
    int i, j, ii, jj, temppass = 0, checkGeoOnly;
    ErrorSet err=0;
    long newSeed = 20080523;
    double latS, lonS, crs1, crs2, crs3, crs4, radius, disti;
    double delta, startcrs;
    LineType lineType;
    ShapeType thisType;
    Geodesic* geo1;
    Geodesic* geo2;
    Geodesic* geo3;
    Geodesic* geo4; 
    Arc* thisArc = NULL;
    Geodesic* thisGeo = NULL;
    Shape* thisShape = NULL;
    LLPoint nbArcStartexp, nbArcEndexp, nbArcCenterexp;
    LLPoint nbGeo1Startexp, nbGeo1Endexp; 
    LLPoint nbGeo2Startexp, nbGeo2Endexp; 
    LLPoint center, geo1Start, geo2Start, geo3Start, geo4Start;
    LLPoint arcStart, arcCenter, arcEnd;
    Arc testArc;
    Boundary b1, newBndry;

    TestSet set;
    set = newTestSet("testBndryCircleIntx_Set2");

    printf("\n\nStart testBndryCircleIntx_Set2\n");

    srand(newSeed);  //Initialize the random number generator
 
    lineType = SEGMENT;

    for (jj = 0; jj < 100; jj++)
    {
    err = 0;
    b1 = createBndry();

    latS = randLat();
    if (latS >= 89.0)
      latS = latS - 1.0;
    latS = DEG2RAD * latS;
    lonS = DEG2RAD * randLon();
    crs1 = DEG2RAD * ((double)((rand() % 89) + 1));
    crs2 = DEG2RAD * ((double)((rand() % 89) + 91));
    crs3 = DEG2RAD * ((double)((rand() % 89) + 181));
    crs4 = DEG2RAD * ((double)((rand() % 89) + 271));
    radius = (double)((rand() % 100) + 1);

    //printf("\nlatS %4.8f lonS %4.8f crs1 %4.8f crs2 %4.8f crs3 %4.8f crs4 %4.8f rad %4.8f\n",latS*RAD2DEG,lonS*RAD2DEG,crs1*RAD2DEG,crs2*RAD2DEG,crs3*RAD2DEG,crs4*RAD2DEG,radius);
    //boundary 1
    center.latitude = latS;
    center.longitude = lonS;
    err |= constructGeoBoundary(center, radius, crs1, crs2, crs3, crs4, &b1);
    geo1 = (Geodesic*) &b1.elements[0];
    geo2 = (Geodesic*) &b1.elements[1];
    geo3 = (Geodesic*) &b1.elements[2];
    geo4 = (Geodesic*) &b1.elements[3];
    geo1Start = geo1->startPoint;
    geo2Start = geo2->startPoint;
    geo3Start = geo3->startPoint;
    geo4Start = geo4->startPoint;
    if (err) printf("\nError(s) occurred during test %d setup: 0x%lx", jj,err);

        if (outputMatlab) displayMatlabGeo(*geo1, "geo1", 0);
        if (outputMatlab) displayMatlabGeo(*geo2, "geo2", 0);
        if (outputMatlab) displayMatlabGeo(*geo3, "geo3", 0);
        if (outputMatlab) displayMatlabGeo(*geo4, "geo4", 0);

    //testArc (which is a circle) has center at vertex of b1
    for (ii = 0; ii < 4; ii++) //each geodesic
    {
      newBndry = createBndry();
      temppass = 0;
      err = 0;
      if (ii == 0)
      {
        if (geo1->length > geo4->length)
          disti = geo4->length/4.0;
        else
          disti = geo1->length/4.0;
        err |= direct(geo1Start, geo1->startAz, disti, &arcStart, EPS);
        err |= direct(geo1Start, geo4->endAz + M_PI, disti, &arcEnd, EPS);
        arcCenter = geo1Start;
      }
      else if (ii == 1)
      {
        if (geo1->length > geo2->length)
          disti = geo2->length/4.0;
        else
          disti = geo1->length/4.0;
        err |= direct(geo2Start, geo2->startAz, disti, &arcStart, EPS);
        err |= direct(geo2Start, geo1->endAz + M_PI, disti, &arcEnd, EPS);
        arcCenter = geo2Start;
      }
      else if (ii == 2)
      {
        if (geo3->length > geo2->length)
          disti = geo2->length/4.0;
        else
          disti = geo3->length/4.0;
        err |= direct(geo3Start, geo3->startAz, disti, &arcStart, EPS);
        err |= direct(geo3Start, geo2->endAz + M_PI, disti, &arcEnd, EPS);
        arcCenter = geo3Start;
      }
      else if (ii == 3)
      {
        if (geo3->length > geo4->length)
          disti = geo4->length/4.0;
        else
          disti = geo3->length/4.0;
        err |= direct(geo4Start, geo4->startAz, disti, &arcStart, EPS);
        err |= direct(geo4Start, geo3->endAz + M_PI, disti, &arcEnd, EPS);
        arcCenter = geo4Start;
      }
        nbArcStartexp = arcStart;
        nbArcEndexp = arcEnd;
        nbArcCenterexp = arcCenter;

        if (ii == 0)
        {
          nbGeo1Startexp = arcCenter;
          nbGeo1Endexp = arcStart;
          nbGeo2Startexp = arcEnd;
          nbGeo2Endexp = arcCenter;
        }
        else
        {
          nbGeo1Startexp = arcEnd;
          nbGeo1Endexp = arcCenter;
          nbGeo2Startexp = arcCenter;
          nbGeo2Endexp = arcStart;
        }
       
        err |= createArc(&testArc, arcCenter, arcStart, arcStart, CLOCKWISE, TOL, EPS);
        err |= bndryCircleIntx(b1, testArc, &newBndry, TOL, EPS);
        if (err) {
          errorCount++;
          failedCount++;
        }
        else {
          for (j = 0; j < newBndry.length; j++)
          {
            thisShape = &newBndry.elements[j];
            if (thisShape == NULL)
            {
              printf("Shape not defined test %d\n",jj);
              temppass = 0; 
            }
            thisType = newBndry.elements[j].type;
            switch (thisType) {
            case ARC:
              thisArc = (Arc*) thisShape;
              if (ptsAreSame(nbArcStartexp, thisArc->startPoint, TESTTOL) && ptsAreSame(nbArcEndexp, thisArc->endPoint, TESTTOL) && ptsAreSame(nbArcCenterexp, thisArc->centerPoint, TESTTOL)) 
                temppass = 1;
               else
               {
                 printf("Failed arc jj = %d ii = %d\n",jj,ii);
                 temppass = 0;
               }
               break;
            case GEODESIC:
              thisGeo = (Geodesic*)thisShape;
              if (j == 0)
              {
                if (ptsAreSame(nbGeo1Startexp, thisGeo->startPoint, TESTTOL) && ptsAreSame(nbGeo1Endexp, thisGeo->endPoint, TESTTOL)) 
                  temppass = 1;
                else
                {
                  printf("Failed vertex geo1 jj = %d ii = %d\n",jj,ii); 
                  temppass = 0;
                } 
              }
              else if (j == 1)
              {
                if (ptsAreSame(nbGeo2Startexp, thisGeo->startPoint, TESTTOL) && ptsAreSame(nbGeo2Endexp, thisGeo->endPoint, TESTTOL)) 
                  temppass = 1;
                else
                {
                  printf("Failed vertex geo2 jj = %d ii = %d\n",jj,ii); 
                  temppass = 0;
                } 
              }
              break;
            default:
              printf("invalid type test %d\n",jj);
              temppass  = 0; 
            }
            if (!temppass)
              break;
          } //for j
        }
        testCaseCount++;
        if (temppass)
          passedCount++;
        else
          failedCount++;
    } //for ii

    //testArc (which is a circle) intersects geodesic between end points
    for (ii = 0; ii < 4; ii++) //each geodesic
    {
      if (ii == 0)
      {
        if (.1 * geo1->length < radius/4.0)
          delta = .1 * geo1->length;
        else
          delta = radius/4.0;
      }
      else if (ii == 1)
      {
        if (.1 * geo2->length < radius/4.0)
          delta = .1 * geo2->length;
        else
          delta = radius/4.0;
      }
      else if (ii == 2)
      {
        if (.1 * geo3->length < radius/4.0)
          delta = .1 * geo3->length;
        else
          delta = radius/4.0;
      }
      else if (ii == 3)
      {
        if (.1 * geo4->length < radius/4.0)
          delta = .1 * geo4->length;
        else
          delta = radius/4.0;
      }
      for (i = 1; i < 5; i++) 
      {
        newBndry = createBndry();
        if (ii == 0)
        {
          disti = geo1->length * .2 * ((double)i) - delta;
          err |= direct(geo1Start, geo1->startAz, disti, &arcStart, EPS);
          startcrs = modcrs(geo1->startAz + M_PI_2);
          disti = geo1->length * .2 * ((double)i) + delta;
          err |= direct(geo1Start, geo1->startAz, disti, &arcEnd, EPS);
        }
        else if (ii == 1)
        {
          disti = geo2->length * .2 * ((double)i) - delta;
          err |= direct(geo2Start, geo2->startAz, disti, &arcStart, EPS);
          startcrs = modcrs(geo2->startAz + M_PI_2);
          disti = geo2->length * .2 * ((double)i) + delta;
          err |= direct(geo2Start, geo2->startAz, disti, &arcEnd, EPS);
        }
        else if (ii == 2)
        {
          disti = geo3->length * .2 * ((double)i) - delta;
          err |= direct(geo3Start, geo3->startAz, disti, &arcStart, EPS);
          startcrs = modcrs(geo3->startAz + M_PI_2);
          disti = geo3->length * .2 * ((double)i) + delta;
          err |= direct(geo3Start, geo3->startAz, disti, &arcEnd, EPS);
        }
        else if (ii == 3)
        {
          disti = geo4->length * .2 * ((double)i) - delta;
          err |= direct(geo4Start, geo4->startAz, disti, &arcStart, EPS);
          startcrs = modcrs(geo4->startAz + M_PI_2);
          disti = geo4->length * .2 * ((double)i) + delta;
          err |= direct(geo4Start, geo4->startAz, disti, &arcEnd, EPS);
        }
        err |= arcFromStartAndEnd(arcStart, startcrs, arcEnd, &testArc, TOL, EPS);
        arcCenter = testArc.centerPoint;
        nbArcStartexp = arcEnd;
        nbArcEndexp = arcStart;
        nbArcCenterexp = arcCenter;
        nbGeo1Startexp = arcStart;
        nbGeo1Endexp = arcEnd;
        err |= createArc(&testArc, arcCenter, arcStart, arcStart, CLOCKWISE, TOL, EPS);
        err |= bndryCircleIntx(b1, testArc, &newBndry, TOL, EPS);
        if (err) {
          errorCount++;
          failedCount++;
        }
        else {
          checkGeoOnly = 0;  
          if (newBndry.length > 2)
          {
            //In this case we just check the side with the current value of ii
            checkGeoOnly = 1;
            printf("In test %d circle intersects more than one geo\n",jj);
          }
          for (j = 0; j < newBndry.length; j++)
          {
            if (checkGeoOnly && (j != ii))
              continue;
            thisShape = &newBndry.elements[j];
            if (thisShape == NULL)
            {
              printf("Shape not defined test %d\n",jj);
              temppass = 0; 
            }
            thisType = newBndry.elements[j].type;
            switch (thisType) {
            case ARC:
              thisArc = (Arc*) thisShape;
              if (ptsAreSame(nbArcStartexp, thisArc->startPoint, TESTTOL) && ptsAreSame(nbArcEndexp, thisArc->endPoint, TESTTOL) && ptsAreSame(nbArcCenterexp, thisArc->centerPoint, TESTTOL)) 
                temppass = 1;
               else
               {
                 printf("Failed arc jj = %d ii = %d i = %d\n",jj,ii,i);
                 temppass = 0;
               }
               break;
            case GEODESIC:
              thisGeo = (Geodesic*)thisShape;
              if (j == 0)
              {
                if (ptsAreSame(nbGeo1Startexp, thisGeo->startPoint, TESTTOL) && ptsAreSame(nbGeo1Endexp, thisGeo->endPoint, TESTTOL)) 
                  temppass = 1;
                else
                {
                  printf("Failed geo1 jj = %d ii = %d i = %d\n",jj,ii,i); 
                  temppass = 0;
                } 
              }
              break;
            default:
              printf("invalid type test %d\n",jj);
              temppass = 0; 
            }
            if (!temppass)
              break;
          } //for j
        }
        testCaseCount++;
        if (temppass)
          passedCount++;
        else
          failedCount++;
    } //for i
    } //for ii
    } //for jj

    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("Finish testBndryCircleIntx_Set2\n\n\n");

    return set;
}

/*
 * NAME: testBndryCircleIntx_Set3
 *
 * DESCRIPTION:
 * 		This function is used to test the bndryCircleIntx function.
 *
 * 		This function runs the test data created by Rich Snow.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testBndryCircleIntx_Set3(TestSet) - A test set with the folling metrics:
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
TestSet testBndryCircleIntx_Set3()
{
    double DEG2RAD = M_PI / 180.0;
    const double RAD2DEG = 180.0 / M_PI;
    int passedCount=0, failedCount=0;
    int errorCount = 0;
    int unverifiedCount = 0;
    int setupFailureCount = 0;
    int testCaseCount = 0;
    int outputMatlab = 0; //0 = don't output matlab code to the console
    int i, j, jj, n1, n2, temppass = 0;
    ErrorSet err=0;
    long newSeed = 20080523;
    double latS, lonS, crs, geolen, disti, perpCrs, locDist;
    double startcrs, delta;
    double fcrs, bcrs;
    LineType lineType;
    ShapeType thisType;
    LLPoint geoStart, geoEnd, geoPt, arcStart, arcEnd, arcCenter, tempLLPoint;
    LLPoint nbArcStartexp, nbArcEndexp, nbArcCenterexp;
    LLPoint nbArc1Startexp, nbArc1Endexp, nbArc1Centerexp;
    LLPoint nbLoc3Startexp, nbLoc3Endexp;
    Arc* thisArc = NULL;
    Locus* thisLocus = NULL;
    Shape* thisShape = NULL;
    Locus* loc3; 
    Locus* loc4;
    Arc* arc1;
    Arc* arc2;
    Arc testArc;
    Boundary b2, newBndry;
    LLPointPair tmpIntx1, tmpIntx2;

    TestSet set;
    set = newTestSet("testBndryCircleIntx_Set3");

    printf("\n\nStart testBndryCircleIntx_Set3\n");

    srand(newSeed);  //Initialize the random number generator

    lineType = SEGMENT;
    for (jj = 0; jj < 100; jj++)
    {
    err = 0;
    b2 = createBndry();

    latS = randLat();
    if (latS >= 89.3)
      latS = 88.0;
    latS = DEG2RAD * latS;
    lonS = DEG2RAD * randLon();
    crs = DEG2RAD * randAzimuth();
    geolen = (double)((rand() % 100) + 1);
    locDist = (double)((rand() % 20) + 1);

    //printf("\nlatS %4.8f lonS %4.8f crs %4.8f geolen %4.8f locDist %4.8f\n",latS*RAD2DEG,lonS*RAD2DEG,crs*RAD2DEG,geolen,locDist);
    geoStart.latitude = latS;
    geoStart.longitude = lonS;

    err |= direct(geoStart, crs, geolen, &geoEnd, EPS);
    //boundary 2
    err |= constructLocusArcBoundary(geoStart, crs, geolen, locDist, &b2);
    loc3 = (Locus*) &b2.elements[0];
    arc1 = (Arc*) &b2.elements[1];
    loc4  = (Locus*) &b2.elements[2];
    arc2 = (Arc*) &b2.elements[3];

    if (err) printf("\nError(s) occurred during test %d setup: 0x%lx", jj,err);

        if (outputMatlab) displayMatlabLocus(*loc3, "loc3", 0);
        if (outputMatlab) displayMatlabArc(*arc1, "arc1", 0);
        if (outputMatlab) displayMatlabLocus(*loc4, "loc4", 0);
        if (outputMatlab) displayMatlabArc(*arc2, "arc2", 0);

    //testArc (which is a circle) intersects loc3 between endpoints of loc3
    for (i = 1; i < 5; i++) //increment disti along geodesic of loci
    {
      err = 0;
      newBndry = createBndry();
      if (.1 * geolen < locDist/3.0)
        delta = .1 * geolen;
      else
        delta = locDist/3.0;
      disti = geolen * .2 * ((double)i) - delta;
      err |= direct(geoStart, crs, disti, &geoPt, EPS);
      err |= ptOnLocusFromGeoPt(*loc3, geoPt, &arcStart, &perpCrs, TOL, EPS);
      startcrs = modcrs(perpCrs + M_PI );
      disti = geolen * .2 * ((double)i) + delta;
      err |= direct(geoStart, crs, disti, &geoPt, EPS);
      err |= ptOnLocusFromGeoPt(*loc3, geoPt, &arcEnd, &perpCrs, TOL, EPS);
      err |= arcFromStartAndEnd(arcStart, startcrs, arcEnd, &testArc, TOL, EPS);
      arcCenter = testArc.centerPoint;
      nbArcStartexp = arcStart;
      nbArcEndexp = arcEnd; 
      nbArcCenterexp = arcCenter;
      nbLoc3Startexp = arcStart;
      nbLoc3Endexp = arcEnd;
        err |= createArc(&testArc, arcCenter, arcStart, arcStart, CLOCKWISE, TOL, EPS);
        if (outputMatlab) displayMatlabArc(testArc, "testArc", 0);
        err |= bndryCircleIntx(b2, testArc, &newBndry, TOL, EPS);

        if (err) {
          printf("\nError(s) occurred during test %d error: 0x%lx", jj,err);
          errorCount++;
          failedCount++;
        }
        else {
          for (j = 0; j < newBndry.length; j++)
          {
            thisShape = &newBndry.elements[j];
            if (thisShape == NULL)
            {
              printf("Shape not defined test %d\n",jj);
              temppass = 0; 
              break;
            }
            thisType = newBndry.elements[j].type;
            switch (thisType) {
            case ARC:
              thisArc = (Arc*) thisShape;
              if (ptsAreSame(nbArcStartexp, thisArc->startPoint, TESTTOL) && ptsAreSame(nbArcEndexp, thisArc->endPoint, TESTTOL) && ptsAreSame(nbArcCenterexp, thisArc->centerPoint, TESTTOL)) 
                temppass = 1;
               else
               {
                 printf("Failed arc jj = %d i = %d\n",jj,i);
                 temppass = 0;
                 break;
               }
               break;
            case LOCUS:
              thisLocus = (Locus*)thisShape;
              if (ptsAreSame(nbLoc3Startexp, thisLocus->locusStart, TESTTOL) && ptsAreSame(nbLoc3Endexp, thisLocus->locusEnd, TESTTOL) && (fabs(locDist - thisLocus->startDist) < TESTTOL)) 
                temppass = 1;
              else
              {
                printf("Failed loc3 jj = %d i = %d\n",jj,i); 
                 if (!ptsAreSame(nbLoc3Startexp, thisLocus->locusStart, TESTTOL))
                 {
                   err |= invCrs(nbLoc3Startexp, thisLocus->locusStart, &fcrs, &bcrs, EPS);
                   printf("loc start pt failed crs exp to calc %8.3f\n",fcrs*RAD2DEG);
                 }
                 if (!ptsAreSame(nbLoc3Endexp, thisLocus->locusEnd, TESTTOL))
                 {
                   err |= invCrs(nbLoc3Endexp, thisLocus->locusEnd, &fcrs, &bcrs, EPS);
                   printf("loc end pt failed crs exp to calc %8.3f\n",fcrs*RAD2DEG);
                 } 
                temppass = 0;
                break;
              } 
              break;
            default:
              printf("invalid type test %d\n",jj);
              temppass = 0; 
            }
            if (!temppass)
              break;
          } //for j
        }
        testCaseCount++;
        if (temppass)
          passedCount++;
        else
          failedCount++;
    } //for i

    //testArc (which is a circle) centered at start point of arc1
    newBndry = createBndry();
    err = 0;
    if (geolen < locDist)
      delta = geolen/4.0;
    else
      delta = locDist/4.0;
     
    if (delta > arc1->radius/4.0)
      delta = arc1->radius/4.0;

    err |= direct(arc1->startPoint, 0.0, delta, &tempLLPoint, EPS);
    err |= createArc(&testArc, arc1->startPoint, tempLLPoint, tempLLPoint, CLOCKWISE, TOL, EPS);
    err |= locusArcIntx(*loc3, testArc.centerPoint, delta, tmpIntx1, &n1, TOL, EPS);
    nbLoc3Endexp = tmpIntx1[0];
    nbLoc3Startexp = arc1->startPoint;
    nbArcCenterexp = testArc.centerPoint;
    nbArcEndexp = tmpIntx1[0];
    err |= arcIntx(arc1->centerPoint, arc1->radius, testArc.centerPoint, testArc.radius, tmpIntx2, &n2, TOL, EPS);
    if (n2 == 2)
    {
      if (ptIsOnArc(arc1->centerPoint, arc1->radius, arc1->startAz, arc1->endAz, arc1->dir, tmpIntx2[0], &err, TOL, EPS))
        nbArcStartexp = tmpIntx2[0];
      else
        nbArcStartexp = tmpIntx2[1];
    }
    else if (n2 == 1)
        nbArcStartexp = tmpIntx2[0];

    nbArc1Startexp = arc1->startPoint;
    nbArc1Endexp = nbArcStartexp;
    nbArc1Centerexp = arc1->centerPoint;

    err |= bndryCircleIntx(b2, testArc, &newBndry, TOL, EPS);

    if (err) {
      printf("\nError(s) occurred during vertex test %d error: 0x%lx", jj,err);
      errorCount++;
      failedCount++;
    }
    else {
          for (j = 0; j < newBndry.length; j++)
          {
            thisShape = &newBndry.elements[j];
            if (thisShape == NULL)
            {
              printf("Shape not defined test %d\n",jj);
              temppass = 0; 
              break;
            }
            thisType = newBndry.elements[j].type;
            switch (thisType) {
            case ARC:
              thisArc = (Arc*) thisShape;
              if (j == 1)
              {
                if (ptsAreSame(nbArc1Startexp, thisArc->startPoint, TESTTOL) && ptsAreSame(nbArc1Endexp, thisArc->endPoint, TESTTOL) && ptsAreSame(nbArc1Centerexp, thisArc->centerPoint, TESTTOL)) 
                  temppass = 1;
                else
                {
                  printf("Failed vertex arc1 jj = %d i = %d\n",jj,i);
                  temppass = 0;
                  break;
                }
              }
              else if (j == 2)
              {
                if (ptsAreSame(nbArcStartexp, thisArc->startPoint, TESTTOL) && ptsAreSame(nbArcEndexp, thisArc->endPoint, TESTTOL) && ptsAreSame(nbArcCenterexp, thisArc->centerPoint, TESTTOL)) 
                  temppass = 1;
                else
                {
                  printf("Failed vertex arc jj = %d i = %d\n",jj,i);
                  temppass = 0;
                  break;
                }
              }
              break;
            case LOCUS:
              thisLocus = (Locus*)thisShape;
              if (ptsAreSame(nbLoc3Startexp, thisLocus->locusStart, TESTTOL) && ptsAreSame(nbLoc3Endexp, thisLocus->locusEnd, TESTTOL) && (fabs(locDist - thisLocus->startDist) < TESTTOL)) 
                temppass = 1;
              else
              {
                printf("Failed vertex loc3 jj = %d i = %d\n",jj,i); 
                temppass = 0;
                break;
              } 
              break;
            default:
              printf("invalid type test %d\n",jj);
              temppass = 0; 
            }
            if (!temppass)
              break;
          } //for j
    }
    testCaseCount++;
        if (temppass)
          passedCount++;
        else
          failedCount++;

    } //for jj

    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("Finish testBndryCircleIntx_Set3\n\n\n");

    return set;
}

/*
 * NAME: testBndryCircleIntx_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the bndryCircleIntx function.
 *
 * 		This function runs all the test cases and returns
 * 		the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testBndryCircleIntx_AllSets(TestSuite) - A test suite with the folling metrics:
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
TestSuite testBndryCircleIntx_AllSets()
{
	TestSuite suite;
	TestSet set1;
	TestSet set2;
	TestSet set3;

    printf("\nStart testBndryCircleIntx_AllSets\n");

    suite = newTestSuite("testBndryCircleIntx_AllSets");

    set1 = testBndryCircleIntx_Set1();
    addTestSet(set1,&suite);

    set2 = testBndryCircleIntx_Set2();
    addTestSet(set2,&suite);

    set3 = testBndryCircleIntx_Set3();
    addTestSet(set3,&suite);

    displayTestSuite(suite);

    printf("\nFinish testBndryCircleIntx_AllSets\n\n\n");

    return suite;
}

/*
 * NAME: testBndryCircleIntxExists_Set1
 *
 * DESCRIPTION:
 * 		This function is used to test the bndryCircleIntxExists function.
 *
 * 		This function runs the test data created by Rich Snow.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testBndryCircleIntxExists_Set1(TestSet) - A test set with the folling metrics:
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
TestSet testBndryCircleIntxExists_Set1()
{
    double DEG2RAD = M_PI / 180.0;
    //const double RAD2DEG = 180.0 / M_PI;
    int passedCount=0, failedCount=0;
    int errorCount = 0;
    int unverifiedCount = 0;
    int setupFailureCount = 0;
    int testCaseCount = 0;
    int outputMatlab = 0; //0 = don't output matlab code to the console
    int i, j, ii, jj;
    int checkSurface= 0, intersect, intersectexp;
    ErrorSet err=0;
    long newSeed = 20080523;
    double latS, lonS, crs1, crs2, crs3, crs4, radius, fcrs, bcrs;
    double circrad, temp, crsStart, distcent;
    double distArray[6] = {.4, .8, .999999, 1.000001, 1.2, 1.5}; 
    double dist1, dist2, dist3, dist4, distmin;
    LineType lineType;
    Geodesic* geo1;
    Geodesic* geo2;
    Geodesic* geo3;
    Geodesic* geo4; 
    LLPoint center, geo1Start, geo2Start, geo3Start, geo4Start;
    LLPoint arcStart, circCenter, tempLLPoint, pt2;
    Arc testArc;
    Boundary b1;

    TestSet set;
    set = newTestSet("testBndryCircleIntxExists_Set1");

    printf("\n\nStart testBndryCircleIntxExists_Set1\n");

    srand(newSeed);  //Initialize the random number generator
 
    lineType = SEGMENT;

    for (jj = 0; jj < 100; jj++)
    {
    err = 0;
    b1 = createBndry();

    latS = DEG2RAD * randLat();
    lonS = DEG2RAD * randLon();
    crs1 = DEG2RAD * ((double)((rand() % 89) + 1));
    crs2 = DEG2RAD * ((double)((rand() % 89) + 91));
    crs3 = DEG2RAD * ((double)((rand() % 89) + 181));
    crs4 = DEG2RAD * ((double)((rand() % 89) + 271));
    radius = (double)((rand() % 100) + 1);
    circrad = 30.0 + 0.01 * randDist();

    //printf("\nlatS %4.8f lonS %4.8f crs1 %4.8f crs2 %4.8f crs3 %4.8f crs4 %4.8f rad %4.8f circrad %4.8f\n",latS*RAD2DEG,lonS*RAD2DEG,crs1*RAD2DEG,crs2*RAD2DEG,crs3*RAD2DEG,crs4*RAD2DEG,radius,circrad);
    //boundary 1
    center.latitude = latS;
    center.longitude = lonS;
    err |= constructGeoBoundary(center, radius, crs1, crs2, crs3, crs4, &b1);
    geo1 = (Geodesic*) &b1.elements[0];
    geo2 = (Geodesic*) &b1.elements[1];
    geo3 = (Geodesic*) &b1.elements[2];
    geo4 = (Geodesic*) &b1.elements[3];
    geo1Start = geo1->startPoint;
    geo2Start = geo2->startPoint;
    geo3Start = geo3->startPoint;
    geo4Start = geo4->startPoint;
    if (err) printf("\nError(s) occurred during test %d setup: 0x%lx", jj,err);

        if (outputMatlab) displayMatlabGeo(*geo1, "geo1", 0);
        if (outputMatlab) displayMatlabGeo(*geo2, "geo2", 0);
        if (outputMatlab) displayMatlabGeo(*geo3, "geo3", 0);
        if (outputMatlab) displayMatlabGeo(*geo4, "geo4", 0);

    //testArc (which is a circle) has center relative to vertex of b1
    for (ii = 0; ii < 4; ii++) //each vertex
    {
      err = 0;
      if (ii == 0)
      {
        tempLLPoint = geo1Start;
        err |= invCrs(tempLLPoint, center, &fcrs, &bcrs, EPS);
        err |= direct(tempLLPoint, fcrs + M_PI, circrad, &circCenter, EPS);
      }
      else if (ii == 1)
      {
        tempLLPoint = geo2Start;
        err |= invCrs(tempLLPoint, center, &fcrs, &bcrs, EPS);
        err |= direct(tempLLPoint, fcrs + M_PI, circrad, &circCenter, EPS);
      }
      else if (ii == 2)
      {
        tempLLPoint = geo3Start;
        err |= invCrs(tempLLPoint, center, &fcrs, &bcrs, EPS);
        err |= direct(tempLLPoint, fcrs + M_PI, circrad, &circCenter, EPS);
      }
      else if (ii == 3)
      {
        tempLLPoint = geo4Start;
        err |= invCrs(tempLLPoint, center, &fcrs, &bcrs, EPS);
        err |= direct(tempLLPoint, fcrs + M_PI, circrad, &circCenter, EPS);
      }

      err |= inverse(circCenter, tempLLPoint, &crsStart, &temp, &distcent, EPS);

      for (j = 0; j < 6; j++)
      { 
        if (j <= 3) 
         err |= direct(circCenter, crsStart, distArray[j] * distcent, &arcStart, EPS);
        else
         err |= direct(circCenter, crsStart, distcent + distArray[j] * 0.5 * radius, &arcStart, EPS);
        err |= createArc(&testArc, circCenter, arcStart, arcStart, CLOCKWISE, TOL, EPS);
        err |= bndryCircleIntxExists(b1, testArc, checkSurface, &intersect, TOL, EPS);

        if (j <= 2)
          intersectexp = 0;
        else
          intersectexp = 1;

        if (err) {
          errorCount++;
          failedCount++;
        }
        else {
          if (intersect == intersectexp)
            passedCount++;
          else
          {
            //printf("Failed vertex jj = %d ii = %d j = %d\n",jj,ii,j);
            failedCount++;
          }
        }
        testCaseCount++;
      } //for j
    } //for ii

    //testArc (which is a circle) intersects interior of boundary 
    err |= projectToGeo(geo1Start, geo1->startAz, center, &pt2, &temp, &dist1, TOL, EPS);
    err |= projectToGeo(geo2Start, geo2->startAz, center, &pt2, &temp, &dist2, TOL, EPS);
    err |= projectToGeo(geo3Start, geo3->startAz, center, &pt2, &temp, &dist3, TOL, EPS);
    err |= projectToGeo(geo4Start, geo4->startAz, center, &pt2, &temp, &dist4, TOL, EPS);
    distmin = dist1;
    if (distmin > dist2)
      distmin = dist2;
    if (distmin > dist3)
      distmin = dist3;
    if (distmin > dist4)
      distmin = dist4;

    err |= direct(center, crs1, distmin / 2.0, &arcStart, EPS);
    err |= createArc(&testArc, center, arcStart, arcStart, CLOCKWISE, TOL, EPS);
    for (i = 0; i < 2; i++)
    {
      if (i == 0)
      {
        checkSurface = 0;
        intersectexp = 0;
      }
      else
      {
        checkSurface = 1;
        intersectexp = 1;
      }
      err |= bndryCircleIntxExists(b1, testArc, checkSurface, &intersect, TOL, EPS);

        if (err) {
          errorCount++;
          failedCount++;
        }
        else {
          if (intersect == intersectexp)
            passedCount++;
          else
          {
            printf("Failed surface jj = %d ii = %d j = %d\n",jj,ii,j);
            failedCount++;
          }
        }
        testCaseCount++;
    } //for i
    } //for jj

    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("Finish testBndryCircleIntxExists_Set1\n\n\n");

    return set;
}

/*
 * NAME: testBndryCircleIntxExists_Set2
 *
 * DESCRIPTION:
 * 		This function is used to test the bndryCircleIntxExists function.
 *
 * 		This function runs the test data created by Rich Snow.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testBndryCircleIntxExists_Set2(TestSet) - A test set with the folling metrics:
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
TestSet testBndryCircleIntxExists_Set2()
{
    double DEG2RAD = M_PI / 180.0;
    //const double RAD2DEG = 180.0 / M_PI;
    int passedCount=0, failedCount=0;
    int errorCount = 0;
    int unverifiedCount = 0;
    int setupFailureCount = 0;
    int testCaseCount = 0;
    int outputMatlab = 0; //0 = don't output matlab code to the console
    int i, j, jj, ii;
    int checkSurface = 0, intersect, intersectexp;
    ErrorSet err=0;
    long newSeed = 20080523;
    double latS, lonS, crs, geolen, disti, perpCrs, locDist;
    double circrad, temp, crsStart, distcent;
    double distArray[6] = {.4, .8, .999999, 1.000001, 1.2, 1.5}; 
    double delta, fcrs, bcrs, tempcrs;
    LineType lineType;
    LLPoint geoStart, geoEnd, geoPt, arcStart, tempLLPoint, circCenter;
    Locus* loc3; 
    Locus* loc4;
    Arc* arc1;
    Arc* arc2;
    Arc testArc;
    Boundary b2;

    TestSet set;
    set = newTestSet("testBndryCircleIntxExists_Set2");

    printf("\n\nStart testBndryCircleIntxExists_Set2\n");

    srand(newSeed);  //Initialize the random number generator

    lineType = SEGMENT;
    for (jj = 0; jj < 100; jj++)
    {
    err = 0;
    b2 = createBndry();

    latS = randLat();
    if (latS >= 89.3)
      latS = 88.0;
    if (latS <= -89.3)
      latS = -88.0;
    latS = DEG2RAD * latS;
    lonS = DEG2RAD * randLon();
    crs = DEG2RAD * randAzimuth();
    geolen = (double)((rand() % 100) + 1);
    locDist = (double)((rand() % 20) + 1);
    circrad = 30.0 + 0.01 * randDist();

    //printf("\nlatS %4.8f lonS %4.8f crs %4.8f geolen %4.8f locDist %4.8f\n",latS*RAD2DEG,lonS*RAD2DEG,crs*RAD2DEG,geolen,locDist);
    geoStart.latitude = latS;
    geoStart.longitude = lonS;

    err |= direct(geoStart, crs, geolen, &geoEnd, EPS);
    //boundary 2
    err |= constructLocusArcBoundary(geoStart, crs, geolen, locDist, &b2);
    loc3 = (Locus*) &b2.elements[0];
    arc1 = (Arc*) &b2.elements[1];
    loc4  = (Locus*) &b2.elements[2];
    arc2 = (Arc*) &b2.elements[3];

    if (err) printf("\nError(s) occurred during test %d setup: 0x%lx", jj,err);

        if (outputMatlab) displayMatlabLocus(*loc3, "loc3", 0);
        if (outputMatlab) displayMatlabArc(*arc1, "arc1", 0);
        if (outputMatlab) displayMatlabLocus(*loc4, "loc4", 0);
        if (outputMatlab) displayMatlabArc(*arc2, "arc2", 0);

    //testArc (which is a circle) has center relative to loc3 
    if (.1 * geolen < locDist/3.0)
      delta = .1 * geolen;
    else
      delta = locDist/3.0;
    circrad = delta;
    for (i = 1; i < 5; i++) //increment disti along geodesic of loci
    {
      err = 0;
      disti = geolen * .2 * ((double)i);
      err |= direct(geoStart, crs, disti, &geoPt, EPS);
      err |= ptOnLocusFromGeoPt(*loc3, geoPt, &tempLLPoint, &perpCrs, TOL, EPS);
      err |= invCrs(tempLLPoint, geoPt, &tempcrs, &temp, EPS);
      tempcrs += M_PI;
      err |= direct(tempLLPoint, tempcrs, circrad, &circCenter, EPS);
      err |= inverse(circCenter, tempLLPoint, &crsStart, &temp, &distcent, EPS);
      for (j = 0; j < 6; j++)
      { 
        if (j <= 3) 
         err |= direct(circCenter, crsStart, distArray[j] * distcent, &arcStart, EPS);
        else
         err |= direct(circCenter, crsStart, distcent + distArray[j] * 0.5 * delta, &arcStart, EPS);
        err |= createArc(&testArc, circCenter, arcStart, arcStart, CLOCKWISE, TOL, EPS);
        err |= bndryCircleIntxExists(b2, testArc, checkSurface, &intersect, TOL, EPS);

        if (j <= 2)
          intersectexp = 0;
        else
          intersectexp = 2;

        if (err) {
          errorCount++;
          failedCount++;
        }
        else {
          if (intersect == intersectexp)
            passedCount++;
          else
          {
            //printf("Failed locus jj = %d i = %d j = %d intsect = %d intsectexp = %d delta = %4.8f\n",jj,i,j,intersect,intersectexp,circrad);
            failedCount++;
          }
        }
        testCaseCount++;
      } //for j

    } //for i

    //testArc (which is a circle) has center relative to start point of arc1
    err = 0;
    if (geolen < locDist)
      delta = geolen/4.0;
    else
      delta = locDist/4.0;
     
    if (delta > arc1->radius/4.0)
      delta = arc1->radius/4.0;

    err |= invCrs(arc1->centerPoint, arc1->startPoint, &fcrs, &bcrs, EPS);
    err |= direct(arc1->startPoint, bcrs + M_PI, delta, &circCenter, EPS);
      err |= inverse(circCenter, arc1->startPoint, &crsStart, &temp, &distcent, EPS);

      for (j = 0; j < 6; j++)
      { 
        if (j <= 3) 
         err |= direct(circCenter, crsStart, distArray[j] * distcent, &arcStart, EPS);
        else
         err |= direct(circCenter, crsStart, distcent + distArray[j] * 0.5 * delta, &arcStart, EPS);
        err |= createArc(&testArc, circCenter, arcStart, arcStart, CLOCKWISE, TOL, EPS);
        err |= bndryCircleIntxExists(b2, testArc, checkSurface, &intersect, TOL, EPS);

        if (j <= 2)
          intersectexp = 0;
        else
          intersectexp = 1;

        if (err) {
          errorCount++;
          failedCount++;
        }
        else {
          if (intersect == intersectexp)
            passedCount++;
          else
          {
            //printf("Failed vertex jj = %d j = %d intsect = %d intsectexp = %d delta = %4.8f\n",jj,j,intersect,intersectexp,delta);
            failedCount++;
          }
        }
        testCaseCount++;
      } //for j

    //testArc (which is a circle) intersects interior of boundary 
    err |= direct(geoStart, crs, 0.5 * geolen, &geoPt, EPS);

    err |= direct(geoStart, crs, 0.5 * geolen + circrad, &arcStart, EPS);
    err |= createArc(&testArc, geoPt, arcStart, arcStart, CLOCKWISE, TOL, EPS);
    for (ii = 0; ii < 2; ii++)
    {
      if (ii == 0)
      {
        checkSurface = 0;
        intersectexp = 0;
      }
      else
      {
        checkSurface = 1;
        intersectexp = 1;
      }
      err |= bndryCircleIntxExists(b2, testArc, checkSurface, &intersect, TOL, EPS);

        if (err) {
          errorCount++;
          failedCount++;
        }
        else {
          if (intersect == intersectexp)
            passedCount++;
          else
          {
            printf("Failed surface jj = %d ii = %d\n",jj,ii);
            failedCount++;
          }
        }
        testCaseCount++;
    } //for ii

    } //for jj

    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("Finish testBndryCircleIntxExists_Set2\n\n\n");

    return set;
}

/*
 * NAME: testBndryCircleIntxExists_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the bndryCircleIntxExists function.
 *
 * 		This function runs all the test cases and returns
 * 		the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testBndryCircleIntxExists_AllSets(TestSuite) - A test suite with the folling metrics:
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
TestSuite testBndryCircleIntxExists_AllSets()
{
	TestSuite suite;
	TestSet set1;
	TestSet set2;

    printf("\nStart testBndryCircleIntxExists_AllSets\n");

    suite = newTestSuite("testBndryCircleIntxExists_AllSets");

    set1 = testBndryCircleIntxExists_Set1();
    addTestSet(set1,&suite);

    set2 = testBndryCircleIntxExists_Set2();
    addTestSet(set2,&suite);


    displayTestSuite(suite);

    printf("\nFinish testBndryCircleIntxExists_AllSets\n\n\n");

    return suite;
}

TestSet testBndryIntxExists_Set1()
{

    double DEG2RAD = M_PI / 180.0;
	TestSet set;
    int passedCount=0, failedCount=0;
    int errorCount = 0;
    int unverifiedCount = 0;
    int setupFailureCount = 0;
    int testCaseCount = 0;

    LLPoint center1, center2, center3;
    double radius;
    double crs;

//    Boundary geoTestBndry = createBndry();
//    Boundary arcTestBndry = createBndry();
//    Boundary locusTestBndry = createBndry();
//    Boundary spiralTestBndry = createBndry();

    Boundary testBndryArray1[4], testBndryArray2[4], testBndryArray3[4];

    int i, j, k;
    int intxExist;

    ErrorSet err=0;

    long newSeed = 20080523;

    set = newTestSet("testBndryIntxExists_Set1");

    printf("\n\nStart testBndryIntxExists_Set1\n");

    srand(newSeed);  //Initialize the random number generator

    //Initialize boundary sets
    for(i = 0; i < 4; i++){
    	testBndryArray1[i] = createBndry();
    	testBndryArray2[i] = createBndry();
    	testBndryArray3[i] = createBndry();
    }


    for(i = 0; i < 1; i++){

    	// Create randomly located set of test boundaries
		center1.latitude = randLat()*DEG2RAD;
		center1.longitude = randLon()*DEG2RAD;

		radius = 50;

		err |= createGeoTestBndry(center1, radius, &testBndryArray1[0]);
		err |= createLocusTestBndry(center1, radius, &testBndryArray1[1]);
		err |= createArcTestBndry(center1, radius, &testBndryArray1[2]);
		err |= createSpiralTestBndry(center1, radius, radius, &testBndryArray1[3]);

		if(err){
			setupFailureCount++;
			err = 0;
			continue;
		}

		// Create second set of test boundaries that will intersect the first set
//		crs = randAzimuth()*DEG2RAD;
		crs = 90*DEG2RAD;
		err |= direct(center1, crs, .5*radius, &center2, EPS);

		err |= createGeoTestBndry(center2, radius, &testBndryArray2[0]);
		err |= createLocusTestBndry(center2, radius, &testBndryArray2[1]);
		err |= createArcTestBndry(center2, radius, &testBndryArray2[2]);
		err |= createSpiralTestBndry(center2, radius, radius, &testBndryArray2[3]);

		if(err){
			setupFailureCount++;
			err = 0;
			continue;
		}

		// Create third set of test boundaries that will not intersect first set
		err |= direct(center1, crs, 5.0*radius, &center3, EPS);

		err |= createGeoTestBndry(center3, radius, &testBndryArray3[0]);
		err |= createLocusTestBndry(center3, radius, &testBndryArray3[1]);
		err |= createArcTestBndry(center3, radius, &testBndryArray3[2]);
		err |= createSpiralTestBndry(center3, radius, radius, &testBndryArray3[3]);

		if(err){
			setupFailureCount++;
			err = 0;
			continue;
		}

		// Test all combinations of intersecting boundary sets
		for(j = 0; j < 4; j++){
			for(k = 0; k < 4; k++){
				testCaseCount++;
				intxExist = -1;

				err |= bndryIntxExists(testBndryArray1[j],testBndryArray2[k], &intxExist, TOL, EPS);

				if(err) {
					errorCount++;
					printf("\nError in running test %d: %#lx\n",testCaseCount, err);
					printf("center1: %lf, %lf\n", center1.latitude/DEG2RAD, center1.longitude/DEG2RAD);
					printf("center3: %lf, %lf\n", center2.latitude/DEG2RAD, center2.longitude/DEG2RAD);
//					displayMatlabBndry(testBndryArray1[j],"testBndryArray1",0);
//					displayMatlabBndry(testBndryArray2[k],"testBndryArray2",0);
					err = 0;
					continue;
				}

				if(intxExist == 1){
					passedCount++;
				} else {
					failedCount++;
					printf("\nFailure running test %d\n",testCaseCount);
					printf("center1: %lf, %lf\n", center1.latitude/DEG2RAD, center1.longitude/DEG2RAD);
					printf("center3: %lf, %lf\n", center2.latitude/DEG2RAD, center2.longitude/DEG2RAD);
					printf("IntersectionExists:%d\n",intxExist);
//					displayMatlabBndry(testBndryArray1[j],"testBndryArray1",0);
//					displayMatlabBndry(testBndryArray2[k],"testBndryArray2",0);
				}

			}
		}

		// Test all combinations of non-intersecting boundary sets
		for(j = 0; j < 4; j++){
			for(k = 0; k < 4; k++){
				testCaseCount++;
				intxExist = -1;

				err |= bndryIntxExists(testBndryArray1[j],testBndryArray3[k], &intxExist, TOL, EPS);

				if(err) {
					errorCount++;
					printf("\nError in running test %d: %#lx\n",testCaseCount, err);
					printf("center1: %lf, %lf\n", center1.latitude, center1.longitude);
					printf("center3: %lf, %lf\n", center3.latitude, center3.longitude);
//					displayMatlabBndry(testBndryArray1[j],"testBndryArray1",0);
//					displayMatlabBndry(testBndryArray3[k],"testBndryArray3",0);
					err = 0;
					continue;
				}

				if(intxExist == 0){
					passedCount++;
				} else {
					failedCount++;
					printf("\nFailure running test %d\n",testCaseCount);
					printf("center1: %lf, %lf\n", center1.latitude, center1.longitude);
					printf("center3: %lf, %lf\n", center3.latitude, center3.longitude);
					printf("IntersectionExists:%d\n",intxExist);
//					displayMatlabBndry(testBndryArray1[j],"testBndryArray1",0);
//					displayMatlabBndry(testBndryArray3[k],"testBndryArray3",0);
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

    printf("Finish testBndryCircleIntxExists_Set2\n\n\n");

    return set;

}

/*
 * NAME: testBndryIntxExists_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the bndryCircleIntxExists function.
 *
 * 		This function runs all the test cases and returns
 * 		the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testBndryIntxExists_AllSets(TestSuite) - A test suite with the folling metrics:
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
TestSuite testBndryIntxExists_AllSets()
{
	TestSuite suite;
	TestSet set1;
//	TestSet set2;

    printf("\nStart testBndryIntxExists_AllSets\n");

    suite = newTestSuite("testBndryIntxExists_AllSets");

    set1 = testBndryIntxExists_Set1();
    addTestSet(set1,&suite);

//    set2 = testBndryIntxExists_Set2();
//    addTestSet(set2,&suite);


    displayTestSuite(suite);

    printf("\nFinish testBndryIntxExists_AllSets\n\n\n");

    return suite;
}

/*
 * NAME: testArcsCoincide_Set1
 *
 * DESCRIPTION:
 * 		This function is used to test the arc coincide function.
 *
 * 		This function runs the test data created by Rich Snow.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testArcsCoincide_Set1(TestSet) - A test set with the folling metrics:
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
TestSet testArcsCoincide_Set1()
{
    double DEG2RAD = M_PI / 180.0;
    int passedCount=0, failedCount=0;
    int errorCount = 0;
    int unverifiedCount = 0;
    int setupFailureCount = 0;
    int testCaseCount = 0;
    ErrorSet err=0;
    LLPoint* comPt1;
    LLPoint* comPt2;
    Arc* myArc;
    LLPoint arc1_centerPoint, arc2_centerPoint, arc1_startPoint, arc2_startPoint, arc1_endPoint, arc2_endPoint;
    LLPoint arc1_centerPoint0, arc1_startPoint0, arc1_endPoint0;
    LLPoint exArc1_center, exArc1_start, exArc1_end, expectedPoint, expectedPoint2;
    LLPoint exArc2_center, exArc2_start, exArc2_end;
    ArcDirection arc1_direction, arc2_direction, expectedDirection, expectedDirection2;
	Arc arc1;
	Arc arc2;
	Arc nullArc;
	Arc expectedArc;
	Arc expectedArc2;
	int coincide;
    int shapeCount = 0;  //The number of shapes that coincide between a circle and arc
    Shape commonShapes[2];  //Common shapes if arc and circle coincide
    Shape nullCommonShapes[0];
	double tol = TOL;
    double radius;

    TestSet set;

    printf("\n\nStart testArcsCoincide_Set1\n");

    set = newTestSet("testArcsCoincide_Set1");
    srand(06042011);

    for (int jj = 0; jj < 10; jj++)
    { 
    radius = 10.0 + 0.01*randDist();
    arc1_centerPoint0.latitude = DEG2RAD*randLat();
    arc1_centerPoint0.longitude = DEG2RAD*randLon();
    err |= direct(arc1_centerPoint0, 7.0*M_PI_4, radius, &arc1_startPoint0, EPS);
    err |= direct(arc1_centerPoint0, M_PI_4, radius, &arc1_endPoint0, EPS);

    // Test Case 1 -- Arcs have different center points
    arc1_centerPoint = arc1_centerPoint0;
    arc1_startPoint = arc1_startPoint0;
    arc1_endPoint = arc1_endPoint0;
    arc1_direction = CLOCKWISE;

    err |= direct(arc1_centerPoint0, 0.0, 10.0, &arc2_centerPoint, EPS);
    err |= direct(arc2_centerPoint, 7.0*M_PI_4, radius, &arc2_startPoint, EPS);
    err |= direct(arc2_centerPoint, M_PI_4, radius, &arc2_endPoint, EPS);
    arc2_direction = CLOCKWISE;

	// create the Arcs
	err = createArc(&arc1, arc1_centerPoint, arc1_startPoint, arc1_endPoint, arc1_direction, tol, EPS);
	err = createArc(&arc2, arc2_centerPoint, arc2_startPoint, arc2_endPoint, arc1_direction, tol, EPS);

	if (err) {
		errorCount++;
		failedCount++;
	}
	else {
		coincide = arcsCoincide(arc1, arc2, commonShapes, &shapeCount, &err, tol, EPS);
		if (err) {
			printf(formatErrorMessage(err));
			errorCount++;
			failedCount++;
		}
		else {
			if (coincide == 0) {
				passedCount++;
				//printf("Test Case 1: Passed\n");
			}
			else {
				printf("Test Case 1: Failed Boolean Test\n");
				failedCount++;
			}
		}
	}
	testCaseCount++;

    // Test Case 2 -- Arcs have different center points that are within tolerance

    err |= direct(arc1_centerPoint, 0.0, 0.5*tol, &arc2_centerPoint, EPS);
    arc2_startPoint = arc1_startPoint; 
    arc2_endPoint = arc1_endPoint; 
    arc2_direction = CLOCKWISE;

    exArc1_center = arc1_centerPoint;
    exArc1_start = arc1_startPoint;
    exArc1_end = arc1_endPoint;
    expectedDirection = CLOCKWISE;

	// create the Arc
    err= 0;
	err = createArc(&arc1, arc1_centerPoint, arc1_startPoint, arc1_endPoint, arc1_direction, tol, EPS);
	err = createArc(&arc2, arc2_centerPoint, arc2_startPoint, arc2_endPoint, arc2_direction, tol, EPS);

	if (err) {
		errorCount++;
		failedCount++;
		printf("Error making arc in Test2\n");
	}
	else {
		coincide = arcsCoincide(arc1, arc2, commonShapes, &shapeCount, &err, tol, EPS);
		if (err) {
			printf("Error in running test 2: %ld\n", err);
			printf("%s",formatErrorMessage(err));
			errorCount++;
			failedCount++;
		}
		else {
			if (coincide == 1) {
				if (shapeCount == 1) {
					err = createArc(&expectedArc, exArc1_center, exArc1_start, exArc1_end, expectedDirection, tol, EPS);
					myArc = (Arc*) &commonShapes[0];
					if (compareArcs(expectedArc, myArc, TESTTOL) == 1) {
				           //printf("Test Case 2: Passed\n");
                                
                                           passedCount++;
					}
					else {
						printf("Test Case 2: Failed Shape Compare Test\n");
						failedCount++;
					}
				}
				else {
					printf("Test Case 2: Failed Shape Count Test.  Expected 2, returned %d\n", shapeCount);
					failedCount++;
				}
			}
			else {
				printf("Test Case 2: Failed Boolean Test\n");
				failedCount++;
			}
		}
	}
	testCaseCount++;

    // Test Case 3 -- Arcs have different radii

    arc2_centerPoint = arc1_centerPoint0;
    err |= direct(arc2_centerPoint, 7.0*M_PI_4, radius + 10.0, &arc2_startPoint, EPS);
    err |= direct(arc2_centerPoint, M_PI_4, radius + 10.0, &arc2_endPoint, EPS);
    arc2_direction = CLOCKWISE;

	// create the Arc
    err= 0;
	err = createArc(&arc1, arc1_centerPoint, arc1_startPoint, arc1_endPoint, arc1_direction, tol, EPS);
	if (err) {
		errorCount++;
		failedCount++;
	}
	else {
		err = 0;
		err = createArc(&arc2, arc2_centerPoint, arc2_startPoint, arc2_endPoint, arc1_direction, tol, EPS);

		if (err) {
			errorCount++;
			failedCount++;
		}
		else {
			err = 0;
			coincide = arcsCoincide(arc1, arc2, commonShapes, &shapeCount, &err, tol, EPS);
			if (err) {
				printf(formatErrorMessage(err));
				errorCount++;
				failedCount++;
			}
			else {
				if (coincide == 0) {
					passedCount++;
					//printf("Test Case 3: Passed\n");
				}
				else {
					printf("Test Case 3: Failed Boolean Test\n");
					failedCount++;
				}
			}
		}
	}
	testCaseCount++;

    // Test Case 4 -- Arcs have different radii that are within tolerance

    arc2_centerPoint = arc1_centerPoint0;
    err |= direct(arc2_centerPoint, 7.0*M_PI_4, radius + 0.5*tol, &arc2_startPoint, EPS);
    err |= direct(arc2_centerPoint, M_PI_4, radius + 0.5*tol, &arc2_endPoint, EPS);
    arc2_direction = CLOCKWISE;

    exArc1_center = arc1_centerPoint0;
    exArc1_start = arc1_startPoint0;
    exArc1_end = arc1_endPoint0;
    expectedDirection = CLOCKWISE;

	// create the Arc
    err= 0;
	err = createArc(&arc1, arc1_centerPoint, arc1_startPoint, arc1_endPoint, arc1_direction, tol, EPS);
	err = createArc(&arc2, arc2_centerPoint, arc2_startPoint, arc2_endPoint, arc1_direction, tol, EPS);

	if (err) {
		errorCount++;
		failedCount++;
	}
	else {
		err = 0;
		coincide = arcsCoincide(arc1, arc2, commonShapes, &shapeCount, &err, tol, EPS);
		if (err) {
    		printf("Error in test 4: %s",formatErrorMessage(err));
			errorCount++;
			failedCount++;
		}
		else {
			if (coincide == 1) {
				if (shapeCount == 1) {
					err = createArc(&expectedArc, exArc1_center, exArc1_start, exArc1_end, expectedDirection, tol, EPS);
					myArc = (Arc*) &commonShapes[0];
					if (compareArcs(expectedArc, myArc, TESTTOL) == 1) {
						passedCount++;
						//printf("Test Case 4: Passed\n");
					}
					else {
						printf("Test Case 4: Failed Shape Compare Test\n");
						failedCount++;
					}
				}
				else {
					printf("Test Case 4: Failed Shape Count Test.  Expected 2, returned %d\n", shapeCount);
					failedCount++;
				}
			}
			else {
				printf("Test Case 4: Failed Boolean Test\n");
				failedCount++;
			}
		}
	}
	testCaseCount++;

    // Test Case 5 -- Arc radius is zero
    arc1_centerPoint = arc1_centerPoint0;
    arc1_startPoint = arc1_centerPoint0;
    arc1_endPoint = arc1_centerPoint0;
    arc1_direction = CLOCKWISE;

    arc2_centerPoint = arc1_centerPoint;
    arc2_startPoint = arc1_startPoint;
    arc2_endPoint = arc1_endPoint;
    arc2_direction = CLOCKWISE;

	// create the Arc
    err= 0;
	err = createArc(&arc1, arc1_centerPoint, arc1_startPoint, arc1_endPoint, arc1_direction, tol, EPS);
	if (err) {
                if (err == RADIUS_OUT_OF_RANGE_ERR) 
                {
		  passedCount++;
		  //printf("Test Case 5: Passed. Expected Error: %s",formatErrorMessage(err));
                }
                else
                {
		  printf("Test Case 5: Failed Test\n");
		  failedCount++;
                }
	}
	else {
		printf("Test Case 5: Failed Test\n");
		failedCount++;
	}
	testCaseCount++;

    // Test Case 6 -- Arc radius is smaller than tolerance
    err |= direct(arc1_centerPoint, M_PI_4, 0.5*tol, &arc1_startPoint, EPS);
    err |= direct(arc1_centerPoint, M_PI_2,  0.5*tol, &arc1_endPoint, EPS);
    arc1_direction = CLOCKWISE;

    arc2_centerPoint = arc1_centerPoint;
    arc2_startPoint = arc1_startPoint;
    arc2_endPoint = arc1_endPoint;
    arc2_direction = CLOCKWISE;

	// create the Arc
    err= 0;
	err = createArc(&arc1, arc1_centerPoint, arc1_startPoint, arc1_endPoint, arc1_direction, tol, EPS);
	err = createArc(&arc2, arc2_centerPoint, arc2_startPoint, arc2_endPoint, arc2_direction, tol, EPS);

	if (err) {
		errorCount++;
		failedCount++;
	}
	else {
		coincide = arcsCoincide(arc1, arc2, commonShapes, &shapeCount, &err, tol, EPS);
		if (err) {
                       if (err == SHAPE_NOT_DEFINED_ERR)
                       {
			 //printf("Test Case 6: Passed. Expected Error: %s", formatErrorMessage(err));
			 passedCount++;
                       }
                       else
                       {
			 printf("Test Case 6: Failed\n");
                         failedCount++;
                       }
		}
		else {
				printf("Test Case 6: Failed to detect error condition\n");
				failedCount++;
		}
	}
	testCaseCount++;

    // Test Case 7 -- Arc radius is larger than 5000 nm
    err |= direct(arc1_centerPoint, M_PI_4, 5010.0, &arc1_startPoint, EPS);
    err |= direct(arc1_centerPoint, M_PI_2, 5010.0, &arc1_endPoint, EPS);
    arc1_direction = CLOCKWISE;

    arc2_centerPoint = arc1_centerPoint;
    arc2_startPoint = arc1_startPoint;
    arc2_endPoint = arc1_endPoint;
    arc2_direction = CLOCKWISE;

    exArc1_center = arc1_centerPoint;
    exArc1_start = arc1_startPoint;
    exArc1_end = arc1_endPoint;
    expectedDirection = CLOCKWISE;

	// create the Arc
    err= 0;
	err = createArc(&arc1, arc1_centerPoint, arc1_startPoint, arc1_endPoint, arc1_direction, tol, EPS);
	err = createArc(&arc2, arc2_centerPoint, arc2_startPoint, arc2_endPoint, arc2_direction, tol, EPS);

	if (err) {
		errorCount++;
		failedCount++;
	}
	else {
		coincide = arcsCoincide(arc1, arc2, commonShapes, &shapeCount, &err, tol, EPS);
		if (err) {
			printf("Error in running test 7: %ld\n", err);
			printf("%s",formatErrorMessage(err));
			errorCount++;
			failedCount++;
		}
		else {
			if (coincide == 1) {
				if (shapeCount == 1) {
					err = createArc(&expectedArc, exArc1_center, exArc1_start, exArc1_end, expectedDirection, tol, EPS);
					myArc = (Arc*) &commonShapes[0];
					if (compareArcs(expectedArc, myArc, TESTTOL) == 1) {
						//printf("Test Case 7: Passed\n");
						passedCount++;
					}
					else {
						printf("Test Case 7: Failed Shape Compare Test\n");
						failedCount++;
					}
				}
				else {
					printf("Test Case 7: Failed Shape Count Test.  Expected 2, returned %d\n", shapeCount);
					failedCount++;
				}
			}
			else {
				printf("Test Case 7: Failed Boolean Test\n");
				failedCount++;
			}
		}
	}
	testCaseCount++;

    // Test Case 8 -- Arc1 is circle
    arc1_startPoint = arc1_startPoint0;
    arc1_endPoint = arc1_startPoint;
    arc1_direction = CLOCKWISE;

    arc2_centerPoint = arc1_centerPoint;
    arc2_startPoint = arc1_startPoint;
    arc2_endPoint = arc1_endPoint0;
    arc2_direction = CLOCKWISE;

	// create the Arc
    err= 0;
	err = createArc(&arc1, arc1_centerPoint, arc1_startPoint, arc1_endPoint, arc1_direction, tol, EPS);
	err = createArc(&arc2, arc2_centerPoint, arc2_startPoint, arc2_endPoint, arc2_direction, tol, EPS);

	if (err) {
		errorCount++;
		failedCount++;
	}
	else {
		coincide = arcsCoincide(arc1, arc2, commonShapes, &shapeCount, &err, tol, EPS);
		if (err) {
			printf("Error in running test 8: %ld\n", err);
			printf("%s",formatErrorMessage(err));
			errorCount++;
			failedCount++;
		}
		else {
			if (coincide == 1) {
				if (shapeCount == 1) {

					myArc = (Arc*) &commonShapes[0];
					if (compareArcs(arc2, myArc, TESTTOL) == 1) {
						passedCount++;
						//printf("Test Case 8: Passed\n");
					}
					else {
						printf("Test Case 8: Failed Shape Compare Test\n");
						failedCount++;
					}
				}
				else {
					printf("Test Case 8: Failed Shape Count Test\n");
					failedCount++;
				}
			}
			else {
				printf("Test Case 8: Failed Boolean Test\n");
				failedCount++;
			}
		}
	}
	testCaseCount++;

    // Test Case 9 -- Arc2 is circle
    arc1_startPoint = arc1_startPoint0;
    arc1_endPoint = arc1_endPoint0;
    arc1_direction = CLOCKWISE;

    arc2_centerPoint = arc1_centerPoint;
    arc2_startPoint = arc1_startPoint;
    arc2_endPoint = arc2_startPoint;
    arc2_direction = CLOCKWISE;

	// create the Arc
    err= 0;
	err = createArc(&arc1, arc1_centerPoint, arc1_startPoint, arc1_endPoint, arc1_direction, tol, EPS);
	err = createArc(&arc2, arc2_centerPoint, arc2_startPoint, arc2_endPoint, arc2_direction, tol, EPS);

	if (err) {
		errorCount++;
		failedCount++;
	}
	else {
		err = 0;
		coincide = arcsCoincide(arc1, arc2, commonShapes, &shapeCount, &err, tol, EPS);
		if (err) {
			printf("Error in running test 9: %ld\n", err);
			printf("%s",formatErrorMessage(err));
			errorCount++;
			failedCount++;
		}
		else {
			if (coincide == 1) {
				if (shapeCount == 1) {

					myArc = (Arc*) &commonShapes[0];
					if (compareArcs(arc1, myArc, TESTTOL) == 1) {
						passedCount++;
						//printf("Test Case 9: Passed\n");
					}
					else {
						printf("Test Case 9: Failed Shape Compare Test\n");
						failedCount++;
					}
				}
				else {
					printf("Test Case 9: Failed Shape Count Test\n");
					failedCount++;
				}
			}
			else {
				printf("Test Case 9: Failed Boolean Test\n");
				failedCount++;
			}
		}
	}
	testCaseCount++;

    // Test Case 10 -- Both arcs are circles
    arc1_startPoint = arc1_startPoint0;
    arc1_endPoint = arc1_startPoint;
    arc1_direction = CLOCKWISE;

    arc2_centerPoint = arc1_centerPoint;
    arc2_startPoint = arc1_startPoint;
    arc2_endPoint = arc1_endPoint;
    arc2_direction = CLOCKWISE;

	// create the Arc
    err= 0;
	err = createArc(&arc1, arc1_centerPoint, arc1_startPoint, arc1_endPoint, arc1_direction, tol, EPS);
	err = createArc(&arc2, arc2_centerPoint, arc2_startPoint, arc2_endPoint, arc2_direction, tol, EPS);

	if (err) {
	  errorCount++;
	  failedCount++;
	}
	else {
		coincide = arcsCoincide(arc1, arc2, commonShapes, &shapeCount, &err, tol, EPS);
		if (err) {
			printf("Error in running test 10: %ld\n", err);
			printf("%s",formatErrorMessage(err));
			errorCount++;
			failedCount++;
		}
		else {
			if (coincide == 1) {
			    if (shapeCount == 1) {
				  myArc = (Arc*) &commonShapes[0];
				  if (compareArcs(arc1, myArc, TESTTOL) == 1) {
					passedCount++;
					//printf("Test Case 10: Passed\n");
			          }
			          else {
				     printf("Test Case 10: Failed Shape Compare Test\n");
				     failedCount++;
			          }
			    }
			    else {
				printf("Test Case 10: Failed Shape Count Test.  Expected 1, returned %d\n", shapeCount);
				failedCount++;
			    }
			}
			else {
				printf("Test Case 10: Failed Boolean Test\n");
				failedCount++;
			}
		}
	}
	testCaseCount++;

    // Test Case 11 -- Arc1s start and end point are within tolerance
    err |= direct(arc1_startPoint, 3.0*M_PI_4, 0.5*tol, &arc1_endPoint, EPS);
    arc1_direction = CLOCKWISE;

    arc2_centerPoint = arc1_centerPoint;
    arc2_startPoint = arc1_startPoint;
    arc2_endPoint = arc1_startPoint;
    arc2_direction = CLOCKWISE;

	// create the Arc
    err= 0;
	err = createArc(&arc1, arc1_centerPoint, arc1_startPoint, arc1_endPoint, arc1_direction, tol, EPS);
	err = createArc(&arc2, arc2_centerPoint, arc2_startPoint, arc2_endPoint, arc2_direction, tol, EPS);

	if (err) {
		errorCount++;
		failedCount++;
	}
	else {
		coincide = arcsCoincide(arc1, arc2, commonShapes, &shapeCount, &err, tol, EPS);
		if (err) {
			printf("Error in running test 11: %ld\n", err);
			printf("%s",formatErrorMessage(err));
			errorCount++;
			failedCount++;
		}
		else {
			if (coincide == 1) {
				if (shapeCount == 1) {
					myArc = (Arc*) &commonShapes[0];
					if (compareArcs(arc1, myArc, TESTTOL) == 1) {
						passedCount++;
						//printf("Test Case 11: Passed\n");
					}
					else {
						printf("Test Case 11: Failed Shape Compare Test\n");
						failedCount++;
					}
				}
				else {
					printf("Test Case 11: Failed Shape Count Test\n");
					failedCount++;
				}
			}
			else {
				printf("Test Case 11: Failed Boolean Test\n");
				failedCount++;
			}
		}
	}
	testCaseCount++;

    // Test Case 12 -- Arc1s start and end point are just outside of tolerance
    err |= direct(arc1_startPoint, 3.0*M_PI_4, 1.1*tol, &arc1_endPoint, EPS);
    arc1_direction = CLOCKWISE;

    arc2_centerPoint = arc1_centerPoint;
    arc2_startPoint = arc1_startPoint;
    arc2_endPoint = arc1_startPoint;
    arc2_direction = CLOCKWISE;

	// create the Arc
    err= 0;
	err = createArc(&arc1, arc1_centerPoint, arc1_startPoint, arc1_endPoint, arc1_direction, tol, EPS);
	err = createArc(&arc2, arc2_centerPoint, arc2_startPoint, arc2_endPoint, arc2_direction, tol, EPS);

	if (err) {
		errorCount++;
		failedCount++;
	}
	else {
		coincide = arcsCoincide(arc1, arc2, commonShapes, &shapeCount, &err, tol, EPS);
		if (err) {
			printf("Error in running test 12: %ld\n", err);
			printf("%s",formatErrorMessage(err));
			errorCount++;
			failedCount++;
		}
		else {
			if (coincide == 1) {
				if (shapeCount == 1) {
					myArc = (Arc*) &commonShapes[0];
					if (compareArcs(arc1, myArc, TESTTOL) == 1) {
						passedCount++;
						//printf("Test Case 12: Passed\n");
					}
					else {
						printf("Test Case 12: Failed Shape Compare Test\n");
						failedCount++;
					}
				}
				else {
					printf("Test Case 12: Failed Shape Count Test\n");
					failedCount++;
				}
			}
			else {
				printf("Test Case 12: Failed Boolean Test\n");
				failedCount++;
			}
		}
	}
	testCaseCount++;

    // Test Case 13 -- Arcs have same start and end points but different orientation
    arc1_startPoint = arc1_startPoint0;
    arc1_endPoint = arc1_endPoint0;
    arc1_direction = CLOCKWISE;

    arc2_centerPoint = arc1_centerPoint;
    arc2_startPoint = arc1_startPoint;
    arc2_endPoint = arc1_endPoint;
    arc2_direction = COUNTERCLOCKWISE;

    expectedPoint = arc1_startPoint;
    expectedPoint2 = arc1_endPoint;

	// create the Arc
    err= 0;
	err = createArc(&arc1, arc1_centerPoint, arc1_startPoint, arc1_endPoint, arc1_direction, tol, EPS);
	err = createArc(&arc2, arc2_centerPoint, arc2_startPoint, arc2_endPoint, arc2_direction, tol, EPS);

	if (err) {
		errorCount++;
		failedCount++;
	}
	else {
		coincide = arcsCoincide(arc1, arc2, commonShapes, &shapeCount, &err, tol, EPS);
		if (err) {
			printf("Error in running test 13: %ld\n", err);
			printf("%s",formatErrorMessage(err));
			errorCount++;
			failedCount++;
		}
		else {
			if (coincide == 1) {
				if (shapeCount == 2) {

					comPt1 = (LLPoint*) &commonShapes[0];
					comPt2 = (LLPoint*) &commonShapes[1];
	        	    if (comPt1 != NULL && comPt2 != NULL &&
	        	    		((ptsAreSame(expectedPoint, *comPt1, TESTTOL) && ptsAreSame(expectedPoint2, *comPt2, TESTTOL)) || (ptsAreSame(expectedPoint, *comPt2, TESTTOL) && ptsAreSame(expectedPoint2, *comPt1, TESTTOL)))) {
						passedCount++;
						//printf("Test Case 13: Passed\n");
					}
					else {
						printf("Test Case 13: Failed Shape Compare Test\n");
						failedCount++;
					}
				}
				else {
					printf("Test Case 13: Failed Shape Count Test\n");
					failedCount++;
				}
			}
			else {
				printf("Test Case 13: Failed Boolean Test\n");
				failedCount++;
			}
		}
	}
	testCaseCount++;

    // Test Case 14 -- Arc2 overlaps part of Arc1
    err |= direct(arc1_centerPoint, 0.0, radius, &arc1_startPoint, EPS);
    err |= direct(arc1_centerPoint, M_PI, radius, &arc1_endPoint, EPS);
    arc1_direction = CLOCKWISE;

    arc2_centerPoint = arc1_centerPoint;
    err |= direct(arc2_centerPoint, M_PI_2, radius, &arc2_startPoint, EPS);
    err |= direct(arc2_centerPoint, 3.0*M_PI_2, radius, &arc2_endPoint, EPS);
    arc2_direction = CLOCKWISE;

    exArc1_center = arc1_centerPoint;
    exArc1_start = arc2_startPoint;
    exArc1_end = arc1_endPoint;
    expectedDirection = CLOCKWISE;

	// create the Arc
    err= 0;
	err = createArc(&arc1, arc1_centerPoint, arc1_startPoint, arc1_endPoint, arc1_direction, tol, EPS);
	err = createArc(&arc2, arc2_centerPoint, arc2_startPoint, arc2_endPoint, arc2_direction, tol, EPS);

	if (err) {
		errorCount++;
		failedCount++;
	}
	else {
		coincide = arcsCoincide(arc1, arc2, commonShapes, &shapeCount, &err, tol, EPS);
		if (err) {
			printf("Error in running test 14: %ld\n", err);
			printf("%s",formatErrorMessage(err));
			errorCount++;
			failedCount++;
		}
		else {
			if (coincide == 1) {
				if (shapeCount == 1) {
					err = createArc(&expectedArc, exArc1_center, exArc1_start, exArc1_end, expectedDirection, tol, EPS);
					myArc = (Arc*) &commonShapes[0];
					if (compareArcs(expectedArc, myArc, TESTTOL) == 1) {
						passedCount++;
						//printf("Test Case 14: Passed\n");
					}
					else {
						printf("Test Case 14: Failed Shape Compare Test\n");
						failedCount++;
					}
				}
				else {
					printf("Test Case 14: Failed Shape Count Test.  Expected 1, returned %d\n", shapeCount);
					failedCount++;
				}
			}
			else {
				printf("Test Case 14: Failed Boolean Test\n");
				failedCount++;
			}
		}
	}
	testCaseCount++;

    // Test Case 15 -- Arc2s start and end points are opposite of Arc1s start and end points
    arc1_startPoint = arc1_startPoint0;
    arc1_endPoint = arc1_endPoint0;
    arc1_direction = CLOCKWISE;

    arc2_centerPoint = arc1_centerPoint;
    arc2_startPoint = arc1_endPoint;
    arc2_endPoint = arc1_startPoint;
    arc2_direction = CLOCKWISE;

    expectedPoint = arc1_startPoint;
    expectedPoint2 = arc1_endPoint;

	// create the Arc
    err= 0;
	err = createArc(&arc1, arc1_centerPoint, arc1_startPoint, arc1_endPoint, arc1_direction, tol, EPS);
	err = createArc(&arc2, arc2_centerPoint, arc2_startPoint, arc2_endPoint, arc2_direction, tol, EPS);
	if (err) {
		printf("Error making Arc1 in Test case 15\n");
		errorCount++;
		failedCount++;
	}
	else {
		coincide = arcsCoincide(arc1, arc2, commonShapes, &shapeCount, &err, tol, EPS);
		if (err) {
			printf("Error in running test 15: %ld\n", err);
			printf("%s",formatErrorMessage(err));
			errorCount++;
			failedCount++;
		}
		else {
			if (coincide == 1) {
				if (shapeCount == 2) {
					comPt1 = (LLPoint*) &commonShapes[0];
					comPt2 = (LLPoint*) &commonShapes[1];
	        	    if (comPt1 != NULL && comPt2 != NULL && ((ptsAreSame(expectedPoint, *comPt1, TESTTOL) &&  ptsAreSame(expectedPoint2, *comPt2, TESTTOL)) || (ptsAreSame(expectedPoint, *comPt2, TESTTOL) &&  ptsAreSame(expectedPoint2, *comPt1, TESTTOL)))) {
						passedCount++;
						//printf("Test Case 15: Passed\n");
					}
					else {
						printf("Test Case 15: Failed Shape Compare Test\n");
						failedCount++;
					}
				}
				else {
					printf("Test Case 15: Failed Shape Count Test\n");
					failedCount++;
				}
			}
			else {
				printf("Test Case 15: Failed Boolean Test\n");
				failedCount++;
			}
		}
	}
	testCaseCount++;

    // Test Case 16 -- Arc1 and Arc2 only have one point in common
    err |= direct(arc1_centerPoint, M_PI_2, radius, &arc1_startPoint, EPS);
    err |= direct(arc1_centerPoint, 0.0, radius, &arc1_endPoint, EPS);
    arc1_direction = COUNTERCLOCKWISE;

    arc2_centerPoint = arc1_centerPoint;
    arc2_startPoint = arc1_startPoint;
    err |= direct(arc2_centerPoint, M_PI, radius, &arc2_endPoint, EPS);
    arc2_direction = CLOCKWISE;

    expectedPoint = arc1_startPoint;

	// create the Arc
    err= 0;
	err = createArc(&arc1, arc1_centerPoint, arc1_startPoint, arc1_endPoint, arc1_direction, tol, EPS);
	err = createArc(&arc2, arc2_centerPoint, arc2_startPoint, arc2_endPoint, arc2_direction, tol, EPS);

	if (err) {
		errorCount++;
		failedCount++;
	}
	else {
		coincide = arcsCoincide(arc1, arc2, commonShapes, &shapeCount, &err, tol, EPS);
		if (err) {
			printf("Error in running test 16: %ld\n", err);
			printf("%s",formatErrorMessage(err));
			errorCount++;
			failedCount++;
		}
		else {
			if (coincide == 1) {
				if (shapeCount == 1) {

					comPt1 = (LLPoint*) &commonShapes[0];

		        	        if (comPt1 != NULL && ptsAreSame(expectedPoint, *comPt1, TESTTOL)) {
						passedCount++;
						//printf("Test Case 16: Passed\n");
					}
					else {
						printf("Test Case 16: Failed Shape Compare Test\n");
						failedCount++;
					}
				}
				else {
					printf("Test Case 16: Failed Shape Count Test.  Expected 1, returned %d\n", shapeCount);
					failedCount++;
				}
			}
			else {
				printf("Test Case 16: Failed Boolean Test\n");
				failedCount++;
			}
		}
	}
	testCaseCount++;

    // Test Case 17 -- Arc1 and Arc2 have an arc and point in common
    err |= direct(arc1_centerPoint, 0.0, radius, &arc1_startPoint, EPS);
    err |= direct(arc1_centerPoint, M_PI, radius, &arc1_endPoint, EPS);
    arc1_direction = CLOCKWISE;

    arc2_centerPoint = arc1_centerPoint;
    arc2_startPoint = arc1_endPoint;
    err |= direct(arc2_centerPoint, M_PI_2, radius, &arc2_endPoint, EPS);
    arc2_direction = CLOCKWISE;

    exArc1_center = arc1_centerPoint;
    exArc1_start = arc1_startPoint;
    exArc1_end = arc2_endPoint;
    expectedDirection = CLOCKWISE;

    expectedPoint = arc1_endPoint;

	// create the Arc
    err= 0;
	err = createArc(&arc1, arc1_centerPoint, arc1_startPoint, arc1_endPoint, arc1_direction, tol, EPS);
	err = createArc(&arc2, arc2_centerPoint, arc2_startPoint, arc2_endPoint, arc2_direction, tol, EPS);

	if (err) {
		errorCount++;
		failedCount++;
	}
	else {
		coincide = arcsCoincide(arc1, arc2, commonShapes, &shapeCount, &err, tol, EPS);
		if (err) {
			printf("Error in running test 17: %ld\n", err);
			printf("%s",formatErrorMessage(err));
			errorCount++;
			failedCount++;
		}
		else {
			if (coincide == 1) {
				if (shapeCount == 2) {
                                        if (commonShapes[0].type == LLPOINT)
                                        {
					  comPt1 = (LLPoint*) &commonShapes[0];
                                          myArc = (Arc*) &commonShapes[1];
                                        }
                                        else
                                        {
					  comPt1 = (LLPoint*) &commonShapes[1];
                                          myArc = (Arc*) &commonShapes[0];
                                        }
					err = createArc(&expectedArc, exArc1_center, exArc1_start, exArc1_end, expectedDirection, tol, EPS);
					if ((compareArcs(expectedArc, myArc, TESTTOL) == 1) && (comPt1 != NULL) && ptsAreSame(expectedPoint, *comPt1, TESTTOL)) {
						passedCount++;
						//printf("Test Case 17: Passed\n");
					}
					else {
						printf("Test Case 17: Failed Shape Compare Test\n");
						failedCount++;
					}
				}
				else {
					printf("Test Case 17: Failed Shape Count Test.  Expected 2, returned %d\n", shapeCount);
					failedCount++;

				}
			}
			else {
				printf("Test Case 17: Failed Boolean Test\n");
				failedCount++;
			}
		}
	}
	testCaseCount++;

    // Test Case 18 -- Arc1 and Arc2 have two arcs in common
    err |= direct(arc1_centerPoint, M_PI_4, radius, &arc1_startPoint, EPS);
    err |= direct(arc1_centerPoint, 3.0*M_PI_2, radius, &arc1_endPoint, EPS);
    arc1_direction = CLOCKWISE;

    arc2_centerPoint = arc1_centerPoint;
    err |= direct(arc2_centerPoint, M_PI, radius, &arc2_startPoint, EPS);
    err |= direct(arc2_centerPoint, M_PI_2, radius, &arc2_endPoint, EPS);
    arc2_direction = CLOCKWISE;

    exArc1_center = arc1_centerPoint;
    exArc1_start = arc1_startPoint;
    exArc1_end =  arc2_endPoint;
    expectedDirection = CLOCKWISE;

    exArc2_center = arc1_centerPoint;
    exArc2_start = arc2_startPoint;
    exArc2_end = arc1_endPoint;
    expectedDirection2 = CLOCKWISE;

	// create the Arc
    err= 0;
	err = createArc(&arc1, arc1_centerPoint, arc1_startPoint, arc1_endPoint, arc1_direction, tol, EPS);
	err = createArc(&arc2, arc2_centerPoint, arc2_startPoint, arc2_endPoint, arc2_direction, tol, EPS);

	if (err) {
		errorCount++;
		failedCount++;
	}
	else {
		coincide = arcsCoincide(arc1, arc2, commonShapes, &shapeCount, &err, tol, EPS);
		if (err) {
			printf("Error in running test 18: %ld\n", err);
			printf("%s",formatErrorMessage(err));
			errorCount++;
			failedCount++;
		}
		else {
			if (coincide == 1) {
				if (shapeCount == 2) {
					err = createArc(&expectedArc, exArc1_center, exArc1_start, exArc1_end, expectedDirection, tol, EPS);
					err = createArc(&expectedArc2, exArc2_center, exArc2_start, exArc2_end, expectedDirection2, tol, EPS);
		        		Arc* myArc1;
		        		Arc* myArc2;
					myArc1 = (Arc*) &commonShapes[0];
					myArc2  = (Arc*) &commonShapes[1];
					if (((compareArcs(expectedArc, myArc1, TESTTOL) == 1) && (compareArcs(expectedArc2, myArc2, TESTTOL) == 1) )|| ((compareArcs(expectedArc2, myArc1, TESTTOL) == 1) && (compareArcs(expectedArc, myArc2, TESTTOL) == 1))) {
						passedCount++;
						//printf("Test Case 18: Passed\n");
					}
					else {
						printf("Test Case 18: Failed Shape Compare Test\n");
						failedCount++;
					}
				}
				else {
					printf("Test Case 18: Failed Shape Count Test.  Expected 2, returned %d\n", shapeCount);
					failedCount++;
				}
			}
			else {
				printf("Test Case 18: Failed Boolean Test\n");
				failedCount++;
			}
		}
	}
	testCaseCount++;

    // Test Case 19 -- Arcs are on the same circle, but do not have any common shapes
    err |= direct(arc1_centerPoint, 0.0, radius, &arc1_startPoint, EPS);
    err |= direct(arc1_centerPoint, M_PI_2, radius, &arc1_endPoint, EPS);
    arc1_direction = CLOCKWISE;

    arc2_centerPoint = arc1_centerPoint;
    err |= direct(arc2_centerPoint, M_PI, radius, &arc2_startPoint, EPS);
    err |= direct(arc2_centerPoint, 3.0*M_PI_2, radius, &arc2_endPoint, EPS);
    arc2_direction = CLOCKWISE;

	// create the Arc
    err= 0;
	err = createArc(&arc1, arc1_centerPoint, arc1_startPoint, arc1_endPoint, arc1_direction, tol, EPS);
	err = createArc(&arc2, arc2_centerPoint, arc2_startPoint, arc2_endPoint, arc2_direction, tol, EPS);

	if (err) {
		errorCount++;
		failedCount++;
	}
	else {
		err = 0;
		coincide = arcsCoincide(arc1, arc2, commonShapes, &shapeCount, &err, tol, EPS);
		if (err) {
			printf("Error in running test 19: %ld\n", err);
			printf("%s",formatErrorMessage(err));
			errorCount++;
			failedCount++;
		}
		else {
			if (coincide == 1 && shapeCount == 0) {
				passedCount++;
				//printf("Test Case 19: Passed\n");
			}
			else {
				printf("Test Case 19: Failed Boolean Test\n");
				failedCount++;
			}
		}
	}
	testCaseCount++;

    // Test Case 20 -- Arc1 is null
	err = 0;
	coincide = arcsCoincide(nullArc, arc2, commonShapes, &shapeCount, &err, tol, EPS);
	if (err) {
		printf("Error in running test 20: %ld\n", err);
		printf("%s",formatErrorMessage(err));
		errorCount++;
		failedCount++;
	}
	else {
		if (coincide == 0) {
			passedCount++;
			//printf("Test Case 20: Passed\n");
		}
		else {
			printf("Test Case 20: Failed Boolean Test\n");
			failedCount++;
		}
	}
	testCaseCount++;

    // Test Case 21 -- Arc2 is null
	err = 0;
	coincide = arcsCoincide(arc1, nullArc, commonShapes, &shapeCount, &err, tol, EPS);
	if (err) {
		printf("Error in running test 21: %ldn", err);
		printf("%s",formatErrorMessage(err));
		errorCount++;
		failedCount++;
	}
	else {
		if (coincide == 0) {
			passedCount++;
			//printf("Test Case 21: Passed\n");
		}
		else {
			printf("Test Case 21: Failed Boolean Test\n");
			failedCount++;
		}
	}
	testCaseCount++;

    // Test Case 22 -- nullcommonShapes[0] is used
    err |= direct(arc1_centerPoint, M_PI_2, radius, &arc1_startPoint, EPS);
    err |= direct(arc1_centerPoint, 0.0, radius, &arc1_endPoint, EPS);
    arc1_direction = COUNTERCLOCKWISE;

    arc2_centerPoint = arc1_centerPoint;
    arc2_startPoint = arc1_startPoint;
    err |= direct(arc2_centerPoint, M_PI, radius, &arc2_endPoint, EPS);
    arc2_direction = CLOCKWISE;

    expectedPoint = arc1_startPoint;

	// create the Arc
    err= 0;
	err = createArc(&arc1, arc1_centerPoint, arc1_startPoint, arc1_endPoint, arc1_direction, tol, EPS);
	err = createArc(&arc2, arc2_centerPoint, arc2_startPoint, arc2_endPoint, arc2_direction, tol, EPS);

	if (err) {
		errorCount++;
		failedCount++;
	}
	else {
		coincide = arcsCoincide(arc1, arc2, nullCommonShapes, &shapeCount, &err, tol, EPS);
		if (err) {
			printf("Error in running test 22: %ld\n", err);
			printf("%s",formatErrorMessage(err));
			errorCount++;
			failedCount++;
		}
		else {
			if (coincide == 1) {
				if (shapeCount == 1) {

					comPt1 = (LLPoint*) &nullCommonShapes[0];
		        	    if (comPt1 != NULL && ptsAreSame(expectedPoint, *comPt1, TESTTOL)) {
						passedCount++;
						//printf("Test Case 22: Passed\n");
					}
					else {
						printf("Test Case 22: Failed Shape Compare Test\n");
						failedCount++;
					}
				}
				else {
					printf("Test Case 22: Failed Shape Count Test.  Expected 1, returned %d\n", shapeCount);
					failedCount++;
				}
			}
			else {
				printf("Test Case 22: Failed Boolean Test\n");
				failedCount++;
			}
		}
	}
	testCaseCount++;

    // Test Case 23 -- shapeCount is null
    err |= direct(arc1_centerPoint, M_PI_2, radius, &arc1_startPoint, EPS);
    err |= direct(arc1_centerPoint, 0.0, radius, &arc1_endPoint, EPS);
    arc1_direction = COUNTERCLOCKWISE;

    arc2_centerPoint = arc1_centerPoint;
    arc2_startPoint = arc1_startPoint;
    err |= direct(arc2_centerPoint, M_PI, radius, &arc2_endPoint, EPS);
    arc2_direction = CLOCKWISE;

    expectedPoint = arc1_startPoint;

	// create the Arc
    err= 0;
	err = createArc(&arc1, arc1_centerPoint, arc1_startPoint, arc1_endPoint, arc1_direction, tol, EPS);
	err = createArc(&arc2, arc2_centerPoint, arc2_startPoint, arc2_endPoint, arc2_direction, tol, EPS);

	if (err) {
		errorCount++;
		failedCount++;
	}
	else {
		coincide = arcsCoincide(arc1, arc2, commonShapes, NULL, &err, tol, EPS);
		if (err) {
                        if (err == NO_MEMORY_ALLOCATED_ERR)
                        {
			  passedCount++;
		          //printf("Test Case 23: Passed. Expected Error: %s",formatErrorMessage(err));
                        }
                        else
                        {
                          printf("Test Case 23: Failed Test\n");
                          failedCount++;
                        }
		}
		else {
				printf("Test Case 23: Failed Test\n");
				failedCount++;
		}
	}
	testCaseCount++;

    // Test Case 24 -- Same arc, Clockwise, subang < pi
    err = direct(arc1_centerPoint, 100.0 * DEG2RAD, radius, &arc1_startPoint, EPS);
    err = direct(arc1_centerPoint, 200.0 * DEG2RAD, radius, &arc1_endPoint, EPS);
    arc1_direction = CLOCKWISE;

    arc2_centerPoint = arc1_centerPoint;
    arc2_startPoint = arc1_startPoint;
    arc2_endPoint = arc1_endPoint;
    arc2_direction = CLOCKWISE;

    exArc1_center = arc1_centerPoint;
    exArc1_start = arc2_startPoint;
    exArc1_end = arc1_endPoint;
    expectedDirection = CLOCKWISE;

	// create the Arc
    err= 0;
	err = createArc(&arc1, arc1_centerPoint, arc1_startPoint, arc1_endPoint, arc1_direction, tol, EPS);
	err = createArc(&arc2, arc2_centerPoint, arc2_startPoint, arc2_endPoint, arc2_direction, tol, EPS);

	if (err) {
		errorCount++;
		failedCount++;
	}
	else {
		coincide = arcsCoincide(arc1, arc2, commonShapes, &shapeCount, &err, tol, EPS);
		if (err) {
			printf("Error in running test 24: %ld\n", err);
			printf("%s",formatErrorMessage(err));
			errorCount++;
			failedCount++;
		}
		else {
			if (coincide == 1) {
				if (shapeCount == 1) {
					err = createArc(&expectedArc, exArc1_center, exArc1_start, exArc1_end, expectedDirection, tol, EPS);
					myArc = (Arc*) &commonShapes[0];
					if (compareArcs(expectedArc, myArc, TESTTOL) == 1) {
						passedCount++;
						//printf("Test Case 24: Passed\n");
					}
					else {
						printf("Test Case 24: Failed Shape Compare Test\n");
						failedCount++;
					}
				}
				else {
					printf("Test Case 24: Failed Shape Count Test.  Expected 1, returned %d\n", shapeCount);
					failedCount++;
				}
			}
			else {
				printf("Test Case 24: Failed Boolean Test\n");
				failedCount++;
			}
		}
	}
	testCaseCount++;

    // Test Case 25 -- Same arc, Clockwise, subang > pi
    err = direct(arc1_centerPoint, 100.0 * DEG2RAD, radius, &arc1_startPoint, EPS);
    err = direct(arc1_centerPoint, 300.0 * DEG2RAD, radius, &arc1_endPoint, EPS);
    arc1_direction = CLOCKWISE;

    arc2_centerPoint = arc1_centerPoint;
    arc2_startPoint = arc1_startPoint;
    arc2_endPoint = arc1_endPoint;
    arc2_direction = CLOCKWISE;

    exArc1_center = arc1_centerPoint;
    exArc1_start = arc2_startPoint;
    exArc1_end = arc1_endPoint;
    expectedDirection = CLOCKWISE;

	// create the Arc
    err= 0;
	err = createArc(&arc1, arc1_centerPoint, arc1_startPoint, arc1_endPoint, arc1_direction, tol, EPS);
	err = createArc(&arc2, arc2_centerPoint, arc2_startPoint, arc2_endPoint, arc2_direction, tol, EPS);

	if (err) {
		errorCount++;
		failedCount++;
	}
	else {
		coincide = arcsCoincide(arc1, arc2, commonShapes, &shapeCount, &err, tol, EPS);
		if (err) {
			printf("Error in running test 25: %ld\n", err);
			printf("%s",formatErrorMessage(err));
			errorCount++;
			failedCount++;
		}
		else {
			if (coincide == 1) {
				if (shapeCount == 1) {
					err = createArc(&expectedArc, exArc1_center, exArc1_start, exArc1_end, expectedDirection, tol, EPS);
					myArc = (Arc*) &commonShapes[0];
					if (compareArcs(expectedArc, myArc, TESTTOL) == 1) {
						passedCount++;
						//printf("Test Case 25: Passed\n");
					}
					else {
						printf("Test Case 25: Failed Shape Compare Test\n");
						failedCount++;
					}
				}
				else {
					printf("Test Case 25: Failed Shape Count Test.  Expected 1, returned %d\n", shapeCount);
					failedCount++;
				}
			}
			else {
				printf("Test Case 25: Failed Boolean Test\n");
				failedCount++;
			}
		}
	}
	testCaseCount++;

    // Test Case 26 -- Same arc, counterClock, subang < pi
    err = direct(arc1_centerPoint, 200.0 * DEG2RAD, radius, &arc1_startPoint, EPS);
    err = direct(arc1_centerPoint, 100.0 * DEG2RAD, radius, &arc1_endPoint, EPS);
    arc1_direction = COUNTERCLOCKWISE;

    arc2_centerPoint = arc1_centerPoint;
    arc2_startPoint = arc1_startPoint;
    arc2_endPoint = arc1_endPoint;
    arc2_direction = COUNTERCLOCKWISE;

    exArc1_center = arc1_centerPoint;
    exArc1_start = arc2_startPoint;
    exArc1_end = arc1_endPoint;
    expectedDirection = COUNTERCLOCKWISE;

	// create the Arc
    err= 0;
	err = createArc(&arc1, arc1_centerPoint, arc1_startPoint, arc1_endPoint, arc1_direction, tol, EPS);
	err = createArc(&arc2, arc2_centerPoint, arc2_startPoint, arc2_endPoint, arc2_direction, tol, EPS);

	if (err) {
		errorCount++;
		failedCount++;
	}
	else {
		coincide = arcsCoincide(arc1, arc2, commonShapes, &shapeCount, &err, tol, EPS);
		if (err) {
			printf("Error in running test 26: %ld\n", err);
			printf("%s",formatErrorMessage(err));
			errorCount++;
			failedCount++;
		}
		else {
			if (coincide == 1) {
				if (shapeCount == 1) {
					err = createArc(&expectedArc, exArc1_center, exArc1_start, exArc1_end, expectedDirection, tol, EPS);
					myArc = (Arc*) &commonShapes[0];
					if (compareArcs(expectedArc, myArc, TESTTOL) == 1) {
						passedCount++;
						//printf("Test Case 26: Passed\n");
					}
					else {
						printf("Test Case 26: Failed Shape Compare Test\n");
						failedCount++;
					}
				}
				else {
					printf("Test Case 26: Failed Shape Count Test.  Expected 1, returned %d\n", shapeCount);
					failedCount++;
				}
			}
			else {
				printf("Test Case 26: Failed Boolean Test\n");
				failedCount++;
			}
		}
	}
	testCaseCount++;

    // Test Case 27 -- Same arc, counterClock, subang > pi
    err = direct(arc1_centerPoint, 300.0 * DEG2RAD, radius, &arc1_startPoint, EPS);
    err = direct(arc1_centerPoint, 100.0 * DEG2RAD, radius, &arc1_endPoint, EPS);
    arc1_direction = COUNTERCLOCKWISE;

    arc2_centerPoint = arc1_centerPoint;
    arc2_startPoint = arc1_startPoint;
    arc2_endPoint = arc1_endPoint;
    arc2_direction = COUNTERCLOCKWISE;

    exArc1_center = arc1_centerPoint;
    exArc1_start = arc2_startPoint;
    exArc1_end = arc1_endPoint;
    expectedDirection = COUNTERCLOCKWISE;

	// create the Arc
    err= 0;
	err = createArc(&arc1, arc1_centerPoint, arc1_startPoint, arc1_endPoint, arc1_direction, tol, EPS);
	err = createArc(&arc2, arc2_centerPoint, arc2_startPoint, arc2_endPoint, arc2_direction, tol, EPS);

	if (err) {
		errorCount++;
		failedCount++;
	}
	else {
		coincide = arcsCoincide(arc1, arc2, commonShapes, &shapeCount, &err, tol, EPS);
		if (err) {
			printf("Error in running test 27: %ld\n", err);
			printf("%s",formatErrorMessage(err));
			errorCount++;
			failedCount++;
		}
		else {
			if (coincide == 1) {
				if (shapeCount == 1) {
					err = createArc(&expectedArc, exArc1_center, exArc1_start, exArc1_end, expectedDirection, tol, EPS);
					myArc = (Arc*) &commonShapes[0];
					if (compareArcs(expectedArc, myArc, TESTTOL) == 1) {
						passedCount++;
						//printf("Test Case 27: Passed\n");
					}
					else {
						printf("Test Case 27: Failed Shape Compare Test\n");
						failedCount++;
					}
				}
				else {
					printf("Test Case 27: Failed Shape Count Test.  Expected 1, returned %d\n", shapeCount);
					failedCount++;
				}
			}
			else {
				printf("Test Case 27: Failed Boolean Test\n");
				failedCount++;
			}
		}
	}
	testCaseCount++;

    // Test Case 28 -- Same arc, Clockwise, subang > pi, cross 2pi
    err = direct(arc1_centerPoint, 300.0 * DEG2RAD, radius, &arc1_startPoint, EPS);
    err = direct(arc1_centerPoint, 200.0 * DEG2RAD, radius, &arc1_endPoint, EPS);
    arc1_direction = CLOCKWISE;

    arc2_centerPoint = arc1_centerPoint;
    arc2_startPoint = arc1_startPoint;
    arc2_endPoint = arc1_endPoint;
    arc2_direction = CLOCKWISE;

    exArc1_center = arc1_centerPoint;
    exArc1_start = arc2_startPoint;
    exArc1_end = arc1_endPoint;
    expectedDirection = CLOCKWISE;

	// create the Arc
    err= 0;
	err = createArc(&arc1, arc1_centerPoint, arc1_startPoint, arc1_endPoint, arc1_direction, tol, EPS);
	err = createArc(&arc2, arc2_centerPoint, arc2_startPoint, arc2_endPoint, arc2_direction, tol, EPS);

	if (err) {
		errorCount++;
		failedCount++;
	}
	else {
		coincide = arcsCoincide(arc1, arc2, commonShapes, &shapeCount, &err, tol, EPS);
		if (err) {
			printf("Error in running test 28: %ld\n", err);
			printf("%s",formatErrorMessage(err));
			errorCount++;
			failedCount++;
		}
		else {
			if (coincide == 1) {
				if (shapeCount == 1) {
					err = createArc(&expectedArc, exArc1_center, exArc1_start, exArc1_end, expectedDirection, tol, EPS);
					myArc = (Arc*) &commonShapes[0];
					if (compareArcs(expectedArc, myArc, TESTTOL) == 1) {
						passedCount++;
						//printf("Test Case 28: Passed\n");
					}
					else {
						printf("Test Case 28: Failed Shape Compare Test\n");
						failedCount++;
					}
				}
				else {
					printf("Test Case 28: Failed Shape Count Test.  Expected 1, returned %d\n", shapeCount);
					failedCount++;
				}
			}
			else {
				printf("Test Case 28: Failed Boolean Test\n");
				failedCount++;
			}
		}
	}
	testCaseCount++;

    // Test Case 29 -- Same arc, counterClock, subang < pi, cross 2pi
    err = direct(arc1_centerPoint, 100.0 * DEG2RAD, radius, &arc1_startPoint, EPS);
    err = direct(arc1_centerPoint, 300.0 * DEG2RAD, radius, &arc1_endPoint, EPS);
    arc1_direction = COUNTERCLOCKWISE;

    arc2_centerPoint = arc1_centerPoint;
    arc2_startPoint = arc1_startPoint;
    arc2_endPoint = arc1_endPoint;
    arc2_direction = COUNTERCLOCKWISE;

    exArc1_center = arc1_centerPoint;
    exArc1_start = arc2_startPoint;
    exArc1_end = arc1_endPoint;
    expectedDirection = COUNTERCLOCKWISE;

	// create the Arc
    err= 0;
	err = createArc(&arc1, arc1_centerPoint, arc1_startPoint, arc1_endPoint, arc1_direction, tol, EPS);
	err = createArc(&arc2, arc2_centerPoint, arc2_startPoint, arc2_endPoint, arc2_direction, tol, EPS);

	if (err) {
		errorCount++;
		failedCount++;
	}
	else {
		coincide = arcsCoincide(arc1, arc2, commonShapes, &shapeCount, &err, tol, EPS);
		if (err) {
			printf("Error in running test 29: %ld\n", err);
			printf("%s",formatErrorMessage(err));
			errorCount++;
			failedCount++;
		}
		else {
			if (coincide == 1) {
				if (shapeCount == 1) {
					err = createArc(&expectedArc, exArc1_center, exArc1_start, exArc1_end, expectedDirection, tol, EPS);
					myArc = (Arc*) &commonShapes[0];
					if (compareArcs(expectedArc, myArc, TESTTOL) == 1) {
						passedCount++;
						//printf("Test Case 29: Passed\n");
					}
					else {
						printf("Test Case 29: Failed Shape Compare Test\n");
						failedCount++;
					}
				}
				else {
					printf("Test Case 29: Failed Shape Count Test.  Expected 1, returned %d\n", shapeCount);
					failedCount++;
				}
			}
			else {
				printf("Test Case 29: Failed Boolean Test\n");
				failedCount++;
			}
		}
	}
	testCaseCount++;

    // Test Case 30 -- Same arc, counterClock, subang > pi, cross 2pi
    err = direct(arc1_centerPoint, 100.0 * DEG2RAD, radius, &arc1_startPoint, EPS);
    err = direct(arc1_centerPoint, 200.0 * DEG2RAD, radius, &arc1_endPoint, EPS);
    arc1_direction = COUNTERCLOCKWISE;

    arc2_centerPoint = arc1_centerPoint;
    arc2_startPoint = arc1_startPoint;
    arc2_endPoint = arc1_endPoint;
    arc2_direction = COUNTERCLOCKWISE;

    exArc1_center = arc1_centerPoint;
    exArc1_start = arc2_startPoint;
    exArc1_end = arc1_endPoint;
    expectedDirection = COUNTERCLOCKWISE;

	// create the Arc
    err= 0;
	err = createArc(&arc1, arc1_centerPoint, arc1_startPoint, arc1_endPoint, arc1_direction, tol, EPS);
	err = createArc(&arc2, arc2_centerPoint, arc2_startPoint, arc2_endPoint, arc2_direction, tol, EPS);

	if (err) {
		errorCount++;
		failedCount++;
	}
	else {
		coincide = arcsCoincide(arc1, arc2, commonShapes, &shapeCount, &err, tol, EPS);
		if (err) {
			printf("Error in running test 30: %ld\n", err);
			printf("%s",formatErrorMessage(err));
			errorCount++;
			failedCount++;
		}
		else {
			if (coincide == 1) {
				if (shapeCount == 1) {
					err = createArc(&expectedArc, exArc1_center, exArc1_start, exArc1_end, expectedDirection, tol, EPS);
					myArc = (Arc*) &commonShapes[0];
					if (compareArcs(expectedArc, myArc, TESTTOL) == 1) {
						passedCount++;
						//printf("Test Case 30: Passed\n");
					}
					else {
						printf("Test Case 30: Failed Shape Compare Test\n");
						failedCount++;
					}
				}
				else {
					printf("Test Case 30: Failed Shape Count Test.  Expected 1, returned %d\n", shapeCount);
					failedCount++;
				}
			}
			else {
				printf("Test Case 30: Failed Boolean Test\n");
				failedCount++;
			}
		}
	}
	testCaseCount++;
    }//for jj

    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("Finish testArcsCoincide_Set1\n\n\n");

    return set;
}

/*
 * NAME: testArcsCoincide_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the boundary's ArcsCoincide function.
 *
 * 		This function runs all the test cases and returns
 * 		the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testArcsCoincide_AllSets(TestSuite) - A test suite with the folling metrics:
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
TestSuite testArcsCoincide_AllSets()
{
	TestSuite suite;
	TestSet set1;

    printf("\nStart testArcsCoincide_AllSets\n");

    suite = newTestSuite("testArcsCoincide_AllSets");

    set1 = testArcsCoincide_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("\nFinish testArcsCoincide_AllSets\n\n\n");

    return suite;
}

int compareArcs(Arc arc1, Arc* arc2, double tol) {

	if (arc2 == NULL)
	{
		return 0;
	}
    if (!ptsAreSame(arc1.centerPoint, arc2->centerPoint, tol))
    {
    	printf("different center points\n");
    	return 0;
    }
    if (!ptsAreSame(arc1.startPoint, arc2->startPoint, tol))
    {
    	printf("different start points\n");
    	return 0;
    }
    if (!ptsAreSame(arc1.endPoint, arc2->endPoint, tol))
    {
    	printf("different end points\n");
    	return 0;
    }
    if (arc1.dir * arc2->dir < 0)
    {
    	printf("different directions\n");
    	return 0;
    }

	return 1;
}

/*
 * NAME: testPtIsInsideBndry_Set1
 *
 * DESCRIPTION:
 * 		This function is used to test the ptIsInsideBndry function.
 *
 * 		This function runs the test data created by Rich Snow
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testPtIsInsideBndry_Set1(TestSet) - A test set with the folling metrics:
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
TestSet testPtIsInsideBndry_Set1()
{
    double DEG2RAD = M_PI / 180.0;
    //const double RAD2DEG = 180.0 / M_PI;
    int passedCount=0, failedCount=0;
    int errorCount = 0;
    int unverifiedCount = 0;
    int setupFailureCount = 0;
    int testCaseCount = 0;
    int outputMatlab = 0; //0 = don't output matlab code to the console
    int inside, insideexp, i, j, ii, jj;
    ErrorSet err=0;
    long newSeed = 20080523;
    double latS, lonS, crs, geolen, distp, disti, distj, perpCrs;
    double tempd1, tempd2, tempd3, crsAtPt, locDist;
    double crsst, crsi, arcradius;
    LineType lineType;
    Geodesic line1, line2;
    LLPoint geoStart, geoEnd, testPoint, geoPt, ptOnLoc;
    Locus loc1, loc2;
    Locus* loc3; 
    Locus* loc4;
    Arc* arc1;
    Arc* arc2;
    Boundary b1, b2;

    TestSet set;

    printf("\n\nStart testPtIsInsideBndry_Set1\n");

    set = newTestSet("testPtIsInsideBndry_Set1");
 
    srand(newSeed);  //Initialize the random number generator

    lineType = SEGMENT;
    for (jj = 0; jj < 100; jj++)
    {
    b1 = createBndry();
    b2 = createBndry();

    latS = randLat();
    if (latS == 90.0)
      latS = 89.0; 
    latS = DEG2RAD * latS;
    lonS = DEG2RAD * randLon();
    crs = DEG2RAD * randAzimuth();
    geolen = (double)((rand() % 100) + 1);
    locDist = (double)((rand() % 20) + 1);

    //printf("\ntest %d latS %4.8f lonS %4.8f crs %4.8f geolen %4.8f locDist %4.8f\n",jj,latS*RAD2DEG,lonS*RAD2DEG,crs*RAD2DEG,geolen,locDist);
    //boundary 1
    geoStart.latitude = latS;
    geoStart.longitude = lonS;
    err |= direct(geoStart, crs, geolen, &geoEnd, EPS);
    err |= createLocus(&loc1, geoStart, geoEnd, 12.0, 2.0, lineType, TOL, EPS);
    err |= createLocus(&loc2, geoStart, geoEnd, -2.0, -12.0, lineType, TOL, EPS);
    err |= createGeo(&line1, loc1.locusStart, loc2.locusStart, lineType, EPS);
    err |= createGeo(&line2, loc1.locusEnd, loc2.locusEnd, lineType, EPS);
    err = addLocusToBndry(&b1, &loc1);
    err = addGeoToBndry(&b1, &line1);
    err = addLocusToBndry(&b1, &loc2);
    err = addGeoToBndry(&b1, &line2);

    //boundary 2
    err |= constructLocusArcBoundary(geoStart, crs, geolen, locDist, &b2);
    loc3 = (Locus*) &b2.elements[0];
    arc1 = (Arc*) &b2.elements[1];
    loc4  = (Locus*) &b2.elements[2];
    arc2 = (Arc*) &b2.elements[3];

    if (err) printf("\nError(s) occurred during setup: 0x%lx", err);

        if (outputMatlab) displayMatlabLocus(*loc3, "loc3", 0);
        if (outputMatlab) displayMatlabArc(*arc1, "arc1", 0);
        if (outputMatlab) displayMatlabLocus(*loc4, "loc4", 0);
        if (outputMatlab) displayMatlabArc(*arc2, "arc2", 0);

    if (jj < 0) //This test is not being used
    //Test points relative to loc1 and loc2 of boundary 1
    for (ii = 0; ii < 2; ii++) //each locus
    {
    for (i = 0; i < 6; i++) //increment disti along geodesic
    {
      insideexp = 1;
      disti = geolen * .2 * ((double)i);
      if (i == 5)
        disti = .99 * geolen;
      err |= direct(geoStart, crs, disti, &geoPt, EPS);
      if (ii == 0)
        err |= ptOnLocusFromGeoPt(loc1, geoPt, &ptOnLoc, &perpCrs, TOL,EPS);
      else
        err |= ptOnLocusFromGeoPt(loc2, geoPt, &ptOnLoc, &perpCrs, TOL,EPS);
      err |= invDist(geoPt, ptOnLoc, &distp, EPS);
      for (j = 0;  j < 9; j++) //increment distj from geoPt
      {
        err = 0;
        distj = distp * .2 * ((double)j);
        if (j == 6)
          distj = distp - 2.0 * TOL;
        else if (j == 8)
          distj = distp + 2.0 * TOL;
        err |= direct(geoPt, perpCrs, distj, &testPoint, EPS);
        if (j >= 7)
          insideexp = 0;
        inside = ptIsInsideBndry(b1, testPoint, &err, TOL, EPS);

        if (err) {
          errorCount++;
          failedCount++;
        }
        else {
          if (inside == insideexp) {
            passedCount++;
          }
          else {
            printf("Failed bndy1 loc jj = %d ii = %d i = %d j = %d\n",jj,ii,i,j);
            failedCount++;
          }
        }
        testCaseCount++;
      } //for j
    } //for i
    } //for ii

    //Test points relative to line1 and line2 of boundary 1
    for (ii = 0; ii < 2; ii++) //each line
    {
    for (i = 0; i < 5; i++) //increment disti along line
    {
      insideexp = 0;
      if (ii == 0)
      {
        disti = line1.length * .2 * ((double)i);
        err |= direct(line1.startPoint, line1.startAz, disti, &geoPt, EPS);
        crsAtPt = geoCrs(line1, geoPt, &tempd1, &tempd2, &tempd3, &err, TOL, EPS);
        perpCrs = modcrs(crsAtPt - M_PI_2);
      }
      else
      {
        disti = line2.length * .2 * ((double)i);
        err |= direct(line2.startPoint, line2.startAz, disti, &geoPt, EPS);
        crsAtPt = geoCrs(line2, geoPt, &tempd1, &tempd2, &tempd3, &err, TOL, EPS);
        perpCrs = modcrs(crsAtPt + M_PI_2);
      }
      for (j = 0; j < 6; j++) //increment distj from geoPt
      {
        err = 0;
        distj = 1.0 + ((double)j) * 3.0;
        if (j == 4)
          distj = 1.1 * TOL;
        else if (j == 5)
        {
          distj = -1.1 * TOL;
          if (i == 0 && ii == 0)
            distj = 0.0;
          insideexp = 1;
        }
        err |= direct(geoPt, perpCrs, distj, &testPoint, EPS);
        inside = ptIsInsideBndry(b1, testPoint, &err, TOL, EPS);
        if (err) {
          printf("Error(s) occurred line test %d ii = %d i = %d j = %d: 0x%lx", jj,ii,i,j,err);
          errorCount++;
          failedCount++;
        }
        else {
          if (inside == insideexp) {
            passedCount++;
          }
          else {
            printf("Failed line jj = %d ii = %d i = %d j = %d\n",jj,ii,i,j);
            failedCount++;
          }
        }
        testCaseCount++;
      } //for j
    } //for i
    } //for ii

    //Test points relative to loc3 and loc4 of boundary 2
    for (ii = 0; ii < 2; ii++)  // each locus
    {
    for (i = 0; i < 6; i++) //increment disti along geodesic
    {
      insideexp = 1;
      disti = geolen * .2 * ((double)i);
      if (i == 5)
        disti = .99 * geolen;
      err |= direct(geoStart, crs, disti, &geoPt, EPS);
      if (ii == 0)
        err |= ptOnLocusFromGeoPt(*loc3, geoPt, &ptOnLoc, &perpCrs, TOL,EPS);
      else
        err |= ptOnLocusFromGeoPt(*loc4, geoPt, &ptOnLoc, &perpCrs, TOL,EPS);
      err |= invDist(geoPt, ptOnLoc, &distp, EPS);
      for (j = 0;  j < 9; j++) //increment distj from geoPt
      {
        err = 0;
        distj = distp * .2 * ((double)j);
        if (j == 6)
          distj = distp - 1.1 * TOL;
        else if (j == 8)
          distj = distp + 1.1 * TOL;
        err |= direct(geoPt, perpCrs, distj, &testPoint, EPS);
        if (j >= 7)
          insideexp = 0;
        inside = ptIsInsideBndry(b2, testPoint, &err, TOL, EPS);

        if (err) {
          printf("Error(s) occurred bndy2 loc test %d ii = %d i = %d j = %d: 0x%lx", jj,ii,i,j,err);
          errorCount++;
          failedCount++;
        }
        else {
          if (inside == insideexp) {
            passedCount++;
          }
          else {
            printf("Failed bndy2 loc jj = %d ii = %d i = %d j = %d\n",jj,ii,i,j);
            failedCount++;
          }
        }
        testCaseCount++;
      } //for j
    } //for i
    } //for ii

    //Test points relative to the arcs of boundary 2
    for (ii = 0; ii < 2; ii++) //each arc
    {
      if (ii == 0)
      {
        err |= invCrs(geoStart, loc3->locusStart, &crsst, &tempd1, EPS);
        arcradius = arc1->radius;
      }
      else
      {
        err |= invCrs(geoEnd, loc4->locusEnd, &crsst, &tempd1, EPS);
        arcradius = arc2->radius;
      }
    for (i = 0; i < 5; i++) //increment crsi from arc center
    {
      insideexp = 1;
      crsi = crsst  + M_PI * .2 * ((double)i);
      for (j = 0;  j < 9; j++) //increment distj along crsi
      {
        err = 0;
        distj = arcradius * .2 * ((double)j);
        if (j == 6)
          distj = arcradius - 1.1 * TOL;
        else if (j == 8)
          distj = arcradius + 1.1 * TOL;
        if (ii == 0)
          err |= direct(geoStart, crsi, distj, &testPoint, EPS);
        else
          err |= direct(geoEnd, crsi, distj, &testPoint, EPS);
        if (j >= 7)
          insideexp = 0;
        inside = ptIsInsideBndry(b2, testPoint, &err, TOL, EPS);

        if (err) {
          printf("Error(s) occurred arc test %d ii = %d i = %d j = %d: 0x%lx", jj,ii,i,j,err);
          errorCount++;
          failedCount++;
        }
        else {
          if (inside == insideexp) {
            passedCount++;
          }
          else {
            failedCount++;
            printf("Failed arc jj = %d ii = %d i = %d j = %d\n",jj,ii,i,j);
          }
        }
        testCaseCount++;
      }//for j
    } //for i
    } //for ii
    } //for jj
    
    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("Finish testPtIsInsideBndry_Set1\n\n\n");

    return set;
}

/*
 * NAME: testPtIsInsideBndry_Set2
 *
 * DESCRIPTION:
 * 		This function is used to test the ptIsInsideBndry function.
 *
 * 		This function runs the test data created by Rich Snow
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testPtIsInsideBndry_Set2(TestSet) - A test set with the folling metrics:
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
TestSet testPtIsInsideBndry_Set2()
{
    double DEG2RAD = M_PI / 180.0;
    //const double RAD2DEG = 180.0 / M_PI;
    int passedCount=0, failedCount=0;
    int errorCount = 0;
    int unverifiedCount = 0;
    int setupFailureCount = 0;
    int testCaseCount = 0;
    int outputMatlab = 0; //0 = don't output matlab code to the console
    int inside, insideexp, i, j, ii, jj;
    ErrorSet err=0;
    long newSeed = 20080523;
    double latS, lonS, crs1, crs2, crs3, crs4, crsi, radius, disti, distj;
    double tempd1, tempd2, tempd3, crsAtPt, perpCrs;
    LineType lineType;
    Geodesic* geo1;
    Geodesic* geo2;
    Geodesic* geo3;
    Geodesic* geo4; 
    LLPoint center, geo1Start, geo2Start, geo3Start, geo4Start, testPoint, intx, geoPt;
    Boundary b1;

    TestSet set;
    set = newTestSet("testPtIsInsideBndry_Set2");

    printf("\n\nStart testPtIsInsideBndry_Set2\n");

    srand(newSeed);  //Initialize the random number generator
 
    lineType = SEGMENT;

    for (jj = 0; jj < 100; jj++)
    {
    err = 0;
    b1 = createBndry();

    latS = DEG2RAD * randLat();
    lonS = DEG2RAD * randLon();
    crs1 = DEG2RAD * ((double)((rand() % 89) + 1));
    crs2 = DEG2RAD * ((double)((rand() % 89) + 91));
    crs3 = DEG2RAD * ((double)((rand() % 89) + 181));
    crs4 = DEG2RAD * ((double)((rand() % 89) + 271));
    radius = (double)((rand() % 100) + 1);

    //printf("\nlatS %4.8f lonS %4.8f crs1 %4.8f crs2 %4.8f crs3 %4.8f crs4 %4.8f rad %4.8f\n",latS*RAD2DEG,lonS*RAD2DEG,crs1*RAD2DEG,crs2*RAD2DEG,crs3*RAD2DEG,crs4*RAD2DEG,radius);
    //boundary 1
    center.latitude = latS;
    center.longitude = lonS;
    err |= constructGeoBoundary(center, radius, crs1, crs2, crs3, crs4, &b1);
    geo1 = (Geodesic*) &b1.elements[0];
    geo2 = (Geodesic*) &b1.elements[1];
    geo3 = (Geodesic*) &b1.elements[2];
    geo4 = (Geodesic*) &b1.elements[3];
    geo1Start = geo1->startPoint;
    geo2Start = geo2->startPoint;
    geo3Start = geo3->startPoint;
    geo4Start = geo4->startPoint;
    if (err) printf("Error(s) occurred during setup test %d: 0x%lx", jj, err);

        if (outputMatlab) displayMatlabGeo(*geo1, "geo1", 0);
        if (outputMatlab) displayMatlabGeo(*geo2, "geo2", 0);
        if (outputMatlab) displayMatlabGeo(*geo3, "geo3", 0);
        if (outputMatlab) displayMatlabGeo(*geo4, "geo4", 0);
    //Test points relative to (on and inside) each geodesic of boundary 1
    insideexp = 1;
    for (ii = 0; ii < 4; ii++) //each geodesic
    {
    for (i = 0; i < 5; i++) //crsi increments
    {
      if (ii == 0)
      {
        crsi = crs1 + (crs2 - crs1) *.2 * ((double)i);
        if (i == 0)
          err |= invDist(center, geo1Start, &disti, EPS);
        else
          err |= crsIntx(center, crsi, &tempd1, &disti, geo1Start, geo1->startAz, &tempd2, &tempd3, &intx, TOL, EPS);
      }
      else if (ii == 1)
      {
        crsi = crs2 + (crs3 - crs2) *.2 * ((double)i);
        if (i == 0)
          err |= invDist(center, geo2Start, &disti, EPS);
        else
          err |= crsIntx(center, crsi, &tempd1, &disti, geo2Start, geo2->startAz, &tempd2, &tempd3, &intx, TOL, EPS);
      }
      else if (ii == 2)
      {
        crsi = crs3 + (crs4 - crs3) *.2 * ((double)i);
        if (i == 0)
          err |= invDist(center, geo3Start, &disti, EPS);
        else
          err |= crsIntx(center, crsi, &tempd1, &disti, geo3Start, geo3->startAz, &tempd2, &tempd3, &intx, TOL, EPS);
      }
      else if (ii == 3)
      {
        crsi = crs4 + modcrs(crs1 - crs4) *.2 * ((double)i);
        if (i == 0)
          err |= invDist(center, geo4Start, &disti, EPS);
        else
          err |= crsIntx(center, crsi, &tempd1, &disti, geo4Start, geo4->startAz, &tempd2, &tempd3, &intx, TOL, EPS);
      }

      for (j = 1; j < 7; j++) //distj increments along each crsi
      {
        err = 0;
        distj = disti * .2 * ((double)j);
        if (j == 6)
          distj = disti - 1.2 * TOL;
        err |= direct(center, crsi, distj, &testPoint, EPS);
        inside = ptIsInsideBndry(b1, testPoint, &err, TOL, EPS);
        if (err) {
          printf("Error(s) occurred during inside %d ii = %d i = %d j = %d: 0x%lx", jj,ii,i,j,err);
          errorCount++;
          failedCount++;
        }
        else {
          if (inside == insideexp) {
            passedCount++;
          }
          else {
            failedCount++;
            printf("Failed inside jj = %d ii = %d i = %d j = %d\n",jj,ii,i,j);
          }
        }
        testCaseCount++;
      } //for j
    } //for i
    } //for ii

    //Test points relative to (outside) each geodesic of boundary 1
    insideexp = 0;
    for (ii = 0; ii < 4; ii++) //each geodesic
    {
    for (i = 0; i < 5; i++) //disti increments along geodesic
    {
      if (ii == 0)
      {
        disti = geo1->length * .2 * ((double)i);
        err |= direct(geo1Start, geo1->startAz, disti, &geoPt, EPS);
        crsAtPt = geoCrs(*geo1, geoPt, &tempd1, &tempd2, &tempd3, &err, TOL, EPS);
        perpCrs = modcrs(crsAtPt - M_PI_2);
      }
      else if (ii == 1)
      {
        disti = geo2->length * .2 * ((double)i);
        err |= direct(geo2Start, geo2->startAz, disti, &geoPt, EPS);
        crsAtPt = geoCrs(*geo2, geoPt, &tempd1, &tempd2, &tempd3, &err, TOL, EPS);
        perpCrs = modcrs(crsAtPt - M_PI_2);
      }
      else if (ii == 2)
      {
        disti = geo3->length * .2 * ((double)i);
        err |= direct(geo3Start, geo3->startAz, disti, &geoPt, EPS);
        crsAtPt = geoCrs(*geo3, geoPt, &tempd1, &tempd2, &tempd3, &err, TOL, EPS);
        perpCrs = modcrs(crsAtPt - M_PI_2);
      }
      else if (ii == 3)
      {
        disti = geo4->length * .2 * ((double)i);
        err |= direct(geo4Start, geo4->startAz, disti, &geoPt, EPS);
        crsAtPt = geoCrs(*geo4, geoPt, &tempd1, &tempd2, &tempd3, &err, TOL, EPS);
        perpCrs = modcrs(crsAtPt - M_PI_2);
      }
      for (j = 0; j < 6; j++) //distj increments from geoPt
      {
        err = 0;
        distj = 1.0 + ((double)j) * 3.0;
        if (j == 5)
          distj = 1.1 * TOL;
        err |= direct(geoPt, perpCrs, distj, &testPoint, EPS);
        inside = ptIsInsideBndry(b1, testPoint, &err, TOL, EPS);
        if (err) {
          printf("Error(s) occurred during outside %d ii = %d i = %d j = %d: 0x%lx", jj,ii,i,j,err);
          errorCount++;
          failedCount++;
        }
        else {
          if (inside == insideexp) {
            passedCount++;
          }
          else {
            printf("Failed outside jj = %d ii = %d i = %d j = %d\n",jj,ii,i,j);
            failedCount++;
          }
        }
        testCaseCount++;
      } //for j
    } //for i
    } //for ii
    } //for jj

    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("Finish testPtIsInsideBndry_Set2\n\n\n");

    return set;
}

/*
 * NAME: testPtIsInsideBndry_Set3
 *
 * DESCRIPTION:
 * 		This function is used to test the ptIsInsideBndry function.
 *
 * 		This function runs the test data created by Rich Snow
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testPtIsInsideBndry_Set3(TestSet) - A test set with the folling metrics:
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
TestSet testPtIsInsideBndry_Set3()
{
    double DEG2RAD = M_PI / 180.0;
    //const double RAD2DEG = 180.0 / M_PI;
    int passedCount=0, failedCount=0;
    int errorCount = 0;
    int unverifiedCount = 0;
    int setupFailureCount = 0;
    int testCaseCount = 0;
    int inside, insideexp, i, j, ii, jj;
    ErrorSet err=0;
    long newSeed = 20080523;
    double latS, lonS, innerRadius, outerRadius, startSubAngle, endSubAngle;
    double length, disti, distj, distp, crsi;
    double tempd1, tempd2, tempd3, crsAtPt, perpCrs, subtendedAngle, raddiff;
    Geodesic* beginGeo;
    Geodesic* endGeo;
    Arc* innerArc;
    Arc* outerArc;
    Locus* startLocLeft;
    Locus* startLocRight;
    Locus* endLocLeft;
    Locus* endLocRight;
    LLPoint center, testPoint,  geoPt, ptOnLoc;
    Boundary b;

    TestSet set;
    set = newTestSet("testPtIsInsideBndry_Set3");

    printf("\n\nStart testPtIsInsideBndry_Set3\n");

    srand(newSeed);  //Initialize the random number generator

    for (jj = 0; jj < 100; jj++)
    {
      b = createBndry();
      latS = DEG2RAD * randLat();
      lonS = DEG2RAD * randLon();
      startSubAngle = DEG2RAD * randAzimuth();
      endSubAngle = startSubAngle + DEG2RAD * ((double)((rand() % 60) + 30)); 
      innerRadius = (double)((rand() % 50) + 10);
      outerRadius = innerRadius + (double)((rand() % 10) + 10);
      length = (double)((rand() % 50) + 10);
      //printf("\nlatS %4.8f lonS %4.8f innerRad %4.8f outerRad %4.8f startAngle %4.8f endAngle %4.8f length %4.8f\n",latS*RAD2DEG,lonS*RAD2DEG,innerRadius,outerRadius,startSubAngle*RAD2DEG,endSubAngle*RAD2DEG,length);

      center.latitude = latS;
      center.longitude = lonS;
      err |= construct2Arc2LocusBoundary(center, innerRadius, outerRadius, startSubAngle, endSubAngle, length, &b);
      beginGeo = (Geodesic*) &b.elements[0];
      startLocLeft = (Locus*) &b.elements[1];
      innerArc = (Arc*) &b.elements[2];
      endLocRight = (Locus*) &b.elements[3];
      endGeo = (Geodesic*) &b.elements[4];
      endLocLeft = (Locus*) &b.elements[5];
      outerArc = (Arc*) &b.elements[6];
      startLocRight = (Locus*) &b.elements[7];

    //Test points relative to beginGeo and endGeo of boundary b
    for (ii = 0; ii < 2; ii++) //each geodesic
    {
    for (i = 0; i < 6; i++) //increment disti along geodesic
    {
      insideexp = 0;
      if (ii == 0)
      {
        disti = beginGeo->length * .2 * ((double)i);
        err |= direct(beginGeo->startPoint, beginGeo->startAz, disti, &geoPt, EPS);
        crsAtPt = geoCrs(*beginGeo, geoPt, &tempd1, &tempd2, &tempd3, &err, TOL, EPS);
        perpCrs = modcrs(crsAtPt - M_PI_2);
      }
      else
      {
        disti = endGeo->length * .2 * ((double)i);
        err |= direct(endGeo->startPoint, endGeo->startAz, disti, &geoPt, EPS);
        crsAtPt = geoCrs(*endGeo, geoPt, &tempd1, &tempd2, &tempd3, &err, TOL, EPS);
        perpCrs = modcrs(crsAtPt - M_PI_2);
      }
      for (j = 0; j < 6; j++) //increment distj from geoPt
      {
        err = 0;
        distj = 1.0 + ((double)j) * 3.0;
        if (j == 4)
          distj = 1.1 * TOL;
        else if (j == 5)
        {
          distj = -1.1 * TOL;
          insideexp = 1;
        }
        err |= direct(geoPt, perpCrs, distj, &testPoint, EPS);
        inside = ptIsInsideBndry(b, testPoint, &err, TOL, EPS);
        if (err) {
          printf("\nError(s) occurred during line %d ii = %d i = %d j = %d: 0x%lx", jj,ii,i,j,err);
          errorCount++;
          failedCount++;
        }
        else {
          if (inside == insideexp) {
            passedCount++;
          }
          else {
            printf("Failed line jj = %d ii = %d i = %d j = %d\n",jj,ii,i,j);
            failedCount++;
          }
        }
        testCaseCount++;
      } //for j
    } //for i
    } //for ii

    //Test points relative to startLocLeft and startLocRight
    for (ii = 0; ii < 2; ii++)  // each locus
    {
    for (i = 0; i < 6; i++) //increment disti along geodesic
    {
      insideexp = 1;
      disti = startLocLeft->geoLength * .2 * ((double)i);
      if (i == 5)
        disti = .99 * startLocLeft->geoLength;
      err |= direct(startLocLeft->geoStart, startLocLeft->geoAz, disti, &geoPt, EPS);
      if (ii == 0)
        err |= ptOnLocusFromGeoPt(*startLocLeft, geoPt, &ptOnLoc, &perpCrs, TOL,EPS);
      else
        err |= ptOnLocusFromGeoPt(*startLocRight, geoPt, &ptOnLoc, &perpCrs, TOL,EPS);
      err |= invDist(geoPt, ptOnLoc, &distp, EPS);
      for (j = 0;  j < 9; j++) //increment distj from geoPt
      {
        err = 0;
        distj = distp * .2 * ((double)j);
        if (j == 6)
          distj = distp - 1.1 * TOL;
        else if (j == 8)
          distj = distp + 1.1 * TOL;
        err |= direct(geoPt, perpCrs, distj, &testPoint, EPS);
        if (j >= 7)
          insideexp = 0;
        inside = ptIsInsideBndry(b, testPoint, &err, TOL, EPS);

        if (err) {
          printf("\nError(s) occurred during startloc %d ii = %d i = %d j = %d: 0x%lx", jj,ii,i,j,err);
          errorCount++;
          failedCount++;
        }
        else {
          if (inside == insideexp) {
            passedCount++;
          }
          else {
            printf("Failed startLoc jj = %d ii = %d i = %d j = %d\n",jj,ii,i,j);
            failedCount++;
          }
        }
        testCaseCount++;
      } //for j
    } //for i
    } //for ii

    //Test points relative to endLocLeft and endLocRight
    for (ii = 0; ii < 2; ii++)  // each locus
    {
    for (i = 0; i < 6; i++) //increment disti along geodesic
    {
      insideexp = 1;
      disti = endLocLeft->geoLength * .2 * ((double)i);
      if (i == 5)
        disti = .99 * endLocLeft->geoLength;
      err |= direct(endLocLeft->geoStart, endLocLeft->geoAz, disti, &geoPt, EPS);
      if (ii == 0)
        err |= ptOnLocusFromGeoPt(*endLocLeft, geoPt, &ptOnLoc, &perpCrs, TOL,EPS);
      else
        err |= ptOnLocusFromGeoPt(*endLocRight, geoPt, &ptOnLoc, &perpCrs, TOL,EPS);
      err |= invDist(geoPt, ptOnLoc, &distp, EPS);
      for (j = 0;  j < 9; j++) //increment distj from geoPt
      {
        err = 0;
        distj = distp * .2 * ((double)j);
        if (j == 6)
          distj = distp -1.1 * TOL;
        else if (j == 8)
          distj = distp + 1.1 * TOL;
        err |= direct(geoPt, perpCrs, distj, &testPoint, EPS);
        if (j >= 7)
          insideexp = 0;
        inside = ptIsInsideBndry(b, testPoint, &err, TOL, EPS);

        if (err) {
          printf("\nError(s) occurred during endloc %d ii = %d i = %d j = %d: 0x%lx", jj,ii,i,j,err);
          errorCount++;
          failedCount++;
        }
        else {
          if (inside == insideexp) {
            passedCount++;
          }
          else {
            printf("Failed endLoc jj = %d ii = %d i = %d j = %d\n",jj,ii,i,j);
            failedCount++;
          }
        }
        testCaseCount++;
      } //for j
    } //for i
    } //for ii

    //Test points relative to the arcs 
      subtendedAngle = computeSubtendedAngle(startSubAngle, endSubAngle, CLOCKWISE);
      raddiff = outerRadius - innerRadius;
    for (i = 0; i < 5; i++) //increment crsi from arc center
    {
      insideexp = 1;
      crsi = startSubAngle  + subtendedAngle * .2 * ((double)i);
      for (j = 0;  j < 10; j++) //increment distj along crsi
      {
        err = 0;
        distj = innerRadius + raddiff * .2 * ((double)j);
        if (j == 6)
          distj = outerRadius - 1.1 * TOL;
        else if (j == 7)
          distj = innerRadius + 1.1 * TOL;
        else if (j == 8)
          distj = outerRadius + 1.1 * TOL;
        else if (j == 9)
          distj = innerRadius - 1.1 * TOL;
        err |= direct(center, crsi, distj, &testPoint, EPS);
        if (j >= 8)
          insideexp = 0;
        inside = ptIsInsideBndry(b, testPoint, &err, TOL, EPS);

        if (err) {
          errorCount++;
          failedCount++;
        }
        else {
          if (inside == insideexp) {
            passedCount++;
          }
          else {
            failedCount++;
            printf("Failed arc jj = %d i = %d j = %d\n",jj,i,j);
          }
        }
        testCaseCount++;
      }//for j
    } //for i

    } //for jj

    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("Finish testPtIsInsideBndry_Set3\n\n\n");

    return set;
}

/*
 * NAME: testPtIsInsideBndry_Set4
 *
 * DESCRIPTION:
 * 		This function is used to test the ptIsInsideBndry function.
 *
 * 		This function runs the test data created by Joe Heidrich
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testPtIsInsideBndry_Set4(TestSet) - A test set with the folling metrics:
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
TestSet testPtIsInsideBndry_Set4()
{
    int passedCount=0, failedCount=0;
    int errorCount = 0;
    int unverifiedCount = 0;
    int setupFailureCount = 0;
    int testCaseCount = 0;
    double latxExp, lonxExp, len, lineAz, startRad, endRad, az, rad, subAngle;
    LLPoint sp, tempPt1, testPt;
    Geodesic line1, line2;
    double az12, dist;
    int i, j, inside, insideexp;
    ErrorSet err = 0;

    Boundary b;

    double DEG2RAD = M_PI / 180;

    Locus* outLoc1;
    Locus* outLoc2;
    Spiral* spiral1;
    Spiral* spiral2;

    double eps = 1.0e-20;
    double tol = 1.37e-9;

    //Percentage of distances to step through
    double distArray[7] = {0, .5, .99, .999999, 1.000001, 1.01, 1.5};
    TestSet set;
    set = newTestSet("testPtIsInsideBndry_Set4");

    printf("\n\nStart testPtIsInsideBndry_Set4\n");

    srand(04012011);

	for (i=0;i<100;i++) {

                err = 0;
		b = createBndry();

		latxExp = randLat();
		lonxExp = randLon();

		err |= createPt(&sp, latxExp*DEG2RAD, lonxExp*DEG2RAD);
		startRad = 30 + 0.01 * randDist() / 2; //0-54 nm
		endRad = 30 + 0.01 * randDist() / 2; //0-54 nm
		lineAz = randAzimuth();
		len = 120 + .01 * randDist();

                err |= createLocusSpiralBndry(sp, lineAz, len, startRad, endRad, &b);
                
                outLoc1 = (Locus*)&b.elements[0];
                outLoc2 = (Locus*)&b.elements[1];
                spiral1 = (Spiral*)&b.elements[2];
                spiral2 = (Spiral*)&b.elements[3];

		//Create Points that are inside/outside
		err |= createGeo(&line1, outLoc1->geoStart, outLoc1->geoEnd, SEGMENT, eps);
		err |= createGeo(&line2, outLoc2->geoStart, outLoc2->geoEnd, SEGMENT, eps);

                if (err) {
                  printf("\nError(s) occurred during setup i = %d: 0x%lx", i,err);
                  continue;
                }

		for (j=0;j<7;j++) {
                        err = 0;
			dist = rand() % (int) floor(line1.length);
			err |= direct(line1.startPoint, line1.startAz, dist, &tempPt1, eps);
			err |= ptOnLocusFromGeoPt(*outLoc1, tempPt1, &testPt, &az12, tol, eps);
			err |= direct(tempPt1, az12, distArray[j] * outLoc1->startDist, &testPt, eps);
			inside = ptIsInsideBndry(b, testPt, &err, tol, eps);
                        if (j <= 3)
                          insideexp = 1;
                        else
                          insideexp = 0;

                        if (err) {
                          errorCount++;
                          failedCount++;
                         }
                         else {
                           if (inside == insideexp) {
                             passedCount++;
                           }
                           else {
                             failedCount++;
                             printf("Failed loc1 i = %d j = %d\n",i,j);
                           }
                         }
                         testCaseCount++;
		}// for j
		
		for (j=0;j<7;j++) {
                        err = 0;
			dist = rand() % (int) floor(line2.length);
			err |= direct(line2.startPoint, line2.startAz, dist, &tempPt1, eps);
			err |= ptOnLocusFromGeoPt(*outLoc2, tempPt1, &testPt, &az12, tol, eps);
			err |= direct(tempPt1, az12, distArray[j] * fabs(outLoc2->startDist), &testPt, eps);
			inside = ptIsInsideBndry(b, testPt, &err, tol, eps);
                        if (j <= 3)
                          insideexp = 1;
                        else
                          insideexp = 0;

                        if (err) {
                          errorCount++;
                          failedCount++;
                         }
                         else {
                           if (inside == insideexp) {
                             passedCount++;
                           }
                           else {
                             failedCount++;
                             printf("Failed loc2 i = %d j = %d\n",i,j);
                           }
                         }
                         testCaseCount++;
		}// for j
		
		for (j=0;j<7;j++) {
                        err = 0;
			subAngle = spiral1->subtendedAngle / DEG2RAD;
			az = spiral1->startAz + (spiral1->dir * (rand() % (int) floor(subAngle)) * DEG2RAD);
			err |= spiralRadius(*spiral1, az, &rad);
			err |= direct(spiral1->centerPoint, az, distArray[j] * rad, &testPt, eps);
			inside = ptIsInsideBndry(b, testPt, &err, tol, eps);
                        if (j <= 3)
                          insideexp = 1;
                        else
                          insideexp = 0;

                        if (err) {
                          errorCount++;
                          failedCount++;
                         }
                         else {
                           if (inside == insideexp) {
                             passedCount++;
                           }
                           else {
                             failedCount++;
                             printf("Failed spiral1 i = %d j = %d\n",i,j);
                           }
                         }
                         testCaseCount++;
		}// for j
		
		for (j=0;j<7;j++) {
                        err = 0;
			subAngle = spiral2->subtendedAngle / DEG2RAD;
			az = spiral2->startAz + (spiral2->dir * (rand() % (int) floor(subAngle))) * DEG2RAD;
			err |= spiralRadius(*spiral2, az, &rad);
			err |= direct(spiral2->centerPoint, az, distArray[j] * rad, &testPt, eps);
			inside = ptIsInsideBndry(b, testPt, &err, tol, eps);
                        if (j <= 3)
                          insideexp = 1;
                        else
                          insideexp = 0;

                        if (err) {
                          errorCount++;
                          failedCount++;
                         }
                         else {
                           if (inside == insideexp) {
                             passedCount++;
                           }
                           else {
                             failedCount++;
                             printf("Failed spiral2 i = %d j = %d\n",i,j);
                           }
                         }
                         testCaseCount++;
		} //for j
	} //for i

    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("Finish testPtIsInsideBndry_Set4\n\n\n");

    return set;
}

/*
 * NAME: testPtIsInsideBndry_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the boundary's ptIsInsideBndry function.
 *
 * 		This function runs all the test cases and returns
 * 		the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testPtIsInsideBndry_AllSets(TestSuite) - A test suite with the folling metrics:
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

TestSuite testPtIsInsideBndry_AllSets()
{
	TestSuite suite;
	TestSet set1;
	TestSet set2;
	TestSet set3;
	TestSet set4;

    printf("\nStart testPtIsInsideBndry_AllSets\n");

    suite = newTestSuite("testPtIsInsideBndry_AllSets");

    set1 = testPtIsInsideBndry_Set1();
    addTestSet(set1,&suite);

    set2 = testPtIsInsideBndry_Set2();
    addTestSet(set2,&suite);

    set3 = testPtIsInsideBndry_Set3();
    addTestSet(set3,&suite);

    set4 = testPtIsInsideBndry_Set4();
    addTestSet(set4,&suite);

    displayTestSuite(suite);

    printf("\nFinish testPtIsInsideBndry_AllSets\n\n\n");

    return suite;
}

/*
 * NAME: testPtsAreInsideBndry_Set1
 *
 * DESCRIPTION:
 * 		This function is used to test the ptsAreInsideBndry function.
 *
 * 		This function runs the test data created by Rich Snow
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testPtsAreInsideBndry_Set1(TestSet) - A test set with the folling metrics:
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
TestSet testPtsAreInsideBndry_Set1()
{
    double DEG2RAD = M_PI / 180.0;
    //const double RAD2DEG = 180.0 / M_PI;
    int passedCount=0, failedCount=0;
    int errorCount = 0;
    int unverifiedCount = 0;
    int setupFailureCount = 0;
    int testCaseCount = 0;
    int outputMatlab = 0; //0 = don't output matlab code to the console
    int inside, insideexp, i, j, ii, jj;
    ErrorSet err=0;
    long newSeed = 20080523;
    double latS, lonS, crs, geolen, distp, disti, distj, perpCrs;
    double tempd1, tempd2, tempd3, crsAtPt, locDist;
    double crsst, crsi, arcradius;
    LineType lineType;
    Geodesic line1, line2;
    LLPoint geoStart, geoEnd, testPoint, geoPt, ptOnLoc;
    Locus loc1, loc2;
    Locus* loc3;
    Locus* loc4;
    Arc* arc1;
    Arc* arc2;
    Boundary b1, b2;
    double latList[1];
    double lonList[1];
    int insideArr[1];
    insideArr[0] = 0;

    TestSet set;

    printf("\n\nStart testPtsAreInsideBndry_Set1\n");

    set = newTestSet("testPtsAreInsideBndry_Set1");

    srand(newSeed);  //Initialize the random number generator

    lineType = SEGMENT;
    for (jj = 0; jj < 100; jj++)
    {
		b1 = createBndry();
		b2 = createBndry();

		latS = randLat();
		if (latS == 90.0)
		  latS = 89.0;

		latS = DEG2RAD * latS;
		lonS = DEG2RAD * randLon();
		crs = DEG2RAD * randAzimuth();
		geolen = (double)((rand() % 100) + 1);
		locDist = (double)((rand() % 20) + 1);

		//printf("\ntest %d latS %4.8f lonS %4.8f crs %4.8f geolen %4.8f locDist %4.8f\n",jj,latS*RAD2DEG,lonS*RAD2DEG,crs*RAD2DEG,geolen,locDist);
		//boundary 1
		geoStart.latitude = latS;
		geoStart.longitude = lonS;
		err |= direct(geoStart, crs, geolen, &geoEnd, EPS);
		err |= createLocus(&loc1, geoStart, geoEnd, 12.0, 2.0, lineType, TOL, EPS);
		err |= createLocus(&loc2, geoStart, geoEnd, -2.0, -12.0, lineType, TOL, EPS);
		err |= createGeo(&line1, loc1.locusStart, loc2.locusStart, lineType, EPS);
		err |= createGeo(&line2, loc1.locusEnd, loc2.locusEnd, lineType, EPS);
		err = addLocusToBndry(&b1, &loc1);
		err = addGeoToBndry(&b1, &line1);
		err = addLocusToBndry(&b1, &loc2);
		err = addGeoToBndry(&b1, &line2);

		//boundary 2
		err |= constructLocusArcBoundary(geoStart, crs, geolen, locDist, &b2);
		loc3 = (Locus*) &b2.elements[0];
		arc1 = (Arc*) &b2.elements[1];
		loc4  = (Locus*) &b2.elements[2];
		arc2 = (Arc*) &b2.elements[3];

		if (err) printf("\nError(s) occurred during setup: 0x%lx", err);

			if (outputMatlab) displayMatlabLocus(*loc3, "loc3", 0);
			if (outputMatlab) displayMatlabArc(*arc1, "arc1", 0);
			if (outputMatlab) displayMatlabLocus(*loc4, "loc4", 0);
			if (outputMatlab) displayMatlabArc(*arc2, "arc2", 0);

		if (jj < 0) //This test is not being used
		//Test points relative to loc1 and loc2 of boundary 1
		for (ii = 0; ii < 2; ii++) //each locus
		{
			for (i = 0; i < 6; i++) //increment disti along geodesic
			{
				insideexp = 1;
				disti = geolen * .2 * ((double)i);

				if (i == 5)
					disti = .99 * geolen;

				err |= direct(geoStart, crs, disti, &geoPt, EPS);

				if (ii == 0)
					err |= ptOnLocusFromGeoPt(loc1, geoPt, &ptOnLoc, &perpCrs, TOL,EPS);
				else
					err |= ptOnLocusFromGeoPt(loc2, geoPt, &ptOnLoc, &perpCrs, TOL,EPS);

				err |= invDist(geoPt, ptOnLoc, &distp, EPS);
				for (j = 0;  j < 9; j++) //increment distj from geoPt
				{
					err = 0;
					distj = distp * .2 * ((double)j);

					if (j == 6)
						distj = distp - 2.0 * TOL;
					else if (j == 8)
						distj = distp + 2.0 * TOL;

					err |= direct(geoPt, perpCrs, distj, &testPoint, EPS);

					if (j >= 7)
						insideexp = 0;

//					inside = ptsAreInsideBndry(b1, testPoint, &err, TOL, EPS);
//					pointsArr[0]  = testPoint;
//					err |= ptsAreInsideBndry(b1, pointsArr, insideArr, 1, TOL, EPS);
					latList[0] = testPoint.latitude;
					lonList[0] = testPoint.longitude;
					err |= ptsAreInsideBndry(b1, latList, lonList, insideArr, 1, TOL, EPS);
					inside = insideArr[0];

					if (err) {
						errorCount++;
						failedCount++;
					} else {
						if (inside == insideexp) {
							passedCount++;
						} else {
							printf("Failed bndy1 loc jj = %d ii = %d i = %d j = %d\n",jj,ii,i,j);
							failedCount++;
						}
					}
					testCaseCount++;
				} //for j
			} //for i
		} //for ii

			//Test points relative to line1 and line2 of boundary 1
		for (ii = 0; ii < 2; ii++) //each line
		{
			for (i = 0; i < 5; i++) //increment disti along line
			{
				insideexp = 0;
				if (ii == 0)
				{
					disti = line1.length * .2 * ((double)i);
					err |= direct(line1.startPoint, line1.startAz, disti, &geoPt, EPS);
					crsAtPt = geoCrs(line1, geoPt, &tempd1, &tempd2, &tempd3, &err, TOL, EPS);
					perpCrs = modcrs(crsAtPt - M_PI_2);
				}
				else
				{
					disti = line2.length * .2 * ((double)i);
					err |= direct(line2.startPoint, line2.startAz, disti, &geoPt, EPS);
					crsAtPt = geoCrs(line2, geoPt, &tempd1, &tempd2, &tempd3, &err, TOL, EPS);
					perpCrs = modcrs(crsAtPt + M_PI_2);
				}

				for (j = 0; j < 6; j++) //increment distj from geoPt
				{
					err = 0;
					distj = 1.0 + ((double)j) * 3.0;
					if (j == 4)
						distj = 1.1 * TOL;
					else if (j == 5)
					{
						distj = -1.1 * TOL;
						if (i == 0 && ii == 0)
							distj = 0.0;

						insideexp = 1;
					}
					err |= direct(geoPt, perpCrs, distj, &testPoint, EPS);

//					inside = ptsAreInsideBndry(b1, testPoint, &err, TOL, EPS);
//					pointsArr[0]  = testPoint;
//					err |= ptsAreInsideBndry(b1, pointsArr, insideArr, 1, TOL, EPS);
					latList[0] = testPoint.latitude;
					lonList[0] = testPoint.longitude;
					err |= ptsAreInsideBndry(b1, latList, lonList, insideArr, 1, TOL, EPS);
					inside = insideArr[0];

					if (err) {
						printf("Error(s) occurred line test %d ii = %d i = %d j = %d: 0x%lx", jj,ii,i,j,err);
						errorCount++;
						failedCount++;
					}
					else {
						if (inside == insideexp) {
							passedCount++;
						} else {
							printf("Failed line jj = %d ii = %d i = %d j = %d\n",jj,ii,i,j);

							failedCount++;
						}
					}
					testCaseCount++;
				} //for j
			} //for i
		} //for ii

			//Test points relative to loc3 and loc4 of boundary 2
		for (ii = 0; ii < 2; ii++)  // each locus
		{
			for (i = 0; i < 6; i++) //increment disti along geodesic
			{
				insideexp = 1;
				disti = geolen * .2 * ((double)i);
				if (i == 5)
					disti = .99 * geolen;

				err |= direct(geoStart, crs, disti, &geoPt, EPS);

				if (ii == 0)
					err |= ptOnLocusFromGeoPt(*loc3, geoPt, &ptOnLoc, &perpCrs, TOL,EPS);
				else
					err |= ptOnLocusFromGeoPt(*loc4, geoPt, &ptOnLoc, &perpCrs, TOL,EPS);

				err |= invDist(geoPt, ptOnLoc, &distp, EPS);
				for (j = 0;  j < 9; j++) //increment distj from geoPt
				{
					err = 0;
					distj = distp * .2 * ((double)j);
					if (j == 6)
					distj = distp - 1.1 * TOL;
					else if (j == 8)
					distj = distp + 1.1 * TOL;

					err |= direct(geoPt, perpCrs, distj, &testPoint, EPS);

					if (j >= 7)
					insideexp = 0;

//					inside = ptsAreInsideBndry(b2, testPoint, &err, TOL, EPS);
//					pointsArr[0]  = testPoint;
//					err |= ptsAreInsideBndry(b2, pointsArr, insideArr, 1, TOL, EPS);
					latList[0] = testPoint.latitude;
					lonList[0] = testPoint.longitude;
					err |= ptsAreInsideBndry(b2, latList, lonList, insideArr, 1, TOL, EPS);
					inside = insideArr[0];

					if (err) {
						printf("Error(s) occurred bndy2 loc test %d ii = %d i = %d j = %d: 0x%lx", jj,ii,i,j,err);
						errorCount++;
						failedCount++;
					}
					else {
						if (inside == insideexp) {
							passedCount++;
						} else {
							printf("Failed bndy2 loc jj = %d ii = %d i = %d j = %d\n",jj,ii,i,j);
							failedCount++;
						}
					}
					testCaseCount++;
				} //for j
			} //for i
		} //for ii

		//Test points relative to the arcs of boundary 2
		for (ii = 0; ii < 2; ii++) //each arc
		{
			if (ii == 0)
			{
				err |= invCrs(geoStart, loc3->locusStart, &crsst, &tempd1, EPS);
				arcradius = arc1->radius;
			}
			else
			{
				err |= invCrs(geoEnd, loc4->locusEnd, &crsst, &tempd1, EPS);
				arcradius = arc2->radius;
			}

			for (i = 0; i < 5; i++) //increment crsi from arc center
			{
				insideexp = 1;
				crsi = crsst  + M_PI * .2 * ((double)i);
				for (j = 0;  j < 9; j++) //increment distj along crsi
				{
					err = 0;
					distj = arcradius * .2 * ((double)j);
					if (j == 6)
						distj = arcradius - 1.1 * TOL;
					else if (j == 8)
						distj = arcradius + 1.1 * TOL;

					if (ii == 0)
						err |= direct(geoStart, crsi, distj, &testPoint, EPS);
					else
						err |= direct(geoEnd, crsi, distj, &testPoint, EPS);

					if (j >= 7)
						insideexp = 0;

//					inside = ptsAreInsideBndry(b2, testPoint, &err, TOL, EPS);
//					pointsArr[0]  = testPoint;
//					err |= ptsAreInsideBndry(b2, pointsArr, insideArr, 1, TOL, EPS);
					latList[0] = testPoint.latitude;
					lonList[0] = testPoint.longitude;
					err |= ptsAreInsideBndry(b2, latList, lonList, insideArr, 1, TOL, EPS);
					inside = insideArr[0];

					if (err) {
						printf("Error(s) occurred arc test %d ii = %d i = %d j = %d: 0x%lx", jj,ii,i,j,err);
						errorCount++;
						failedCount++;
					} else {
						if (inside == insideexp) {
							passedCount++;
						} else {
							failedCount++;

							printf("Failed arc jj = %d ii = %d i = %d j = %d\n",jj,ii,i,j);
//							if(jj == 54 && ii == 0 && i == 0 && j == 5)
//							{
//								displayDataBndry(b2,"b2",1);
//								displayDataPt(testPoint,"testPt",1);
//								displayMatlabBndry(b2,"b2",0);
//								displayMatlabPt(testPoint,"testPt",0);
//								displayBndry(b2,"b2",0);
//								displayPt(testPoint,"testPt",0);
//								printf("\n%d\n",inside);
//							}
						}
					}

					testCaseCount++;
				}//for j
			} //for i
		} //for ii
    } //for jj

    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("Finish testPtsAreInsideBndry_Set1\n\n\n");

    return set;
}

/*
 * NAME: testPtsAreInsideBndry_Set2
 *
 * DESCRIPTION:
 * 		This function is used to test the ptsAreInsideBndry function.
 *
 * 		This function runs the test data created by Rich Snow
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testPtsAreInsideBndry_Set2(TestSet) - A test set with the folling metrics:
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
TestSet testPtsAreInsideBndry_Set2()
{
    double DEG2RAD = M_PI / 180.0;
    //const double RAD2DEG = 180.0 / M_PI;
    int passedCount=0, failedCount=0;
    int errorCount = 0;
    int unverifiedCount = 0;
    int setupFailureCount = 0;
    int testCaseCount = 0;
    int outputMatlab = 0; //0 = don't output matlab code to the console
    int inside, insideexp, i, j, ii, jj;
    ErrorSet err=0;
    long newSeed = 20080523;
    double latS, lonS, crs1, crs2, crs3, crs4, crsi, radius, disti, distj;
    double tempd1, tempd2, tempd3, crsAtPt, perpCrs;
    LineType lineType;
    Geodesic* geo1;
    Geodesic* geo2;
    Geodesic* geo3;
    Geodesic* geo4;
    LLPoint center, geo1Start, geo2Start, geo3Start, geo4Start, testPoint, intx, geoPt;
    Boundary b1;
    double latList[1];
    double lonList[1];
    int insideArr[1];

    TestSet set;
    set = newTestSet("testPtsAreInsideBndry_Set2");

    printf("\n\nStart testPtsAreInsideBndry_Set2\n");

    srand(newSeed);  //Initialize the random number generator

    lineType = SEGMENT;

    for (jj = 0; jj < 100; jj++)
    {
    err = 0;
    b1 = createBndry();

    latS = DEG2RAD * randLat();
    lonS = DEG2RAD * randLon();
    crs1 = DEG2RAD * ((double)((rand() % 89) + 1));
    crs2 = DEG2RAD * ((double)((rand() % 89) + 91));
    crs3 = DEG2RAD * ((double)((rand() % 89) + 181));
    crs4 = DEG2RAD * ((double)((rand() % 89) + 271));
    radius = (double)((rand() % 100) + 1);

    //printf("\nlatS %4.8f lonS %4.8f crs1 %4.8f crs2 %4.8f crs3 %4.8f crs4 %4.8f rad %4.8f\n",latS*RAD2DEG,lonS*RAD2DEG,crs1*RAD2DEG,crs2*RAD2DEG,crs3*RAD2DEG,crs4*RAD2DEG,radius);
    //boundary 1
    center.latitude = latS;
    center.longitude = lonS;
    err |= constructGeoBoundary(center, radius, crs1, crs2, crs3, crs4, &b1);
    geo1 = (Geodesic*) &b1.elements[0];
    geo2 = (Geodesic*) &b1.elements[1];
    geo3 = (Geodesic*) &b1.elements[2];
    geo4 = (Geodesic*) &b1.elements[3];
    geo1Start = geo1->startPoint;
    geo2Start = geo2->startPoint;
    geo3Start = geo3->startPoint;
    geo4Start = geo4->startPoint;
    if (err) printf("Error(s) occurred during setup test %d: 0x%lx", jj, err);

        if (outputMatlab) displayMatlabGeo(*geo1, "geo1", 0);
        if (outputMatlab) displayMatlabGeo(*geo2, "geo2", 0);
        if (outputMatlab) displayMatlabGeo(*geo3, "geo3", 0);
        if (outputMatlab) displayMatlabGeo(*geo4, "geo4", 0);
    //Test points relative to (on and inside) each geodesic of boundary 1

        if(jj == 23) continue; //Polar Case

    insideexp = 1;
    for (ii = 0; ii < 4; ii++) //each geodesic
    {
    for (i = 0; i < 5; i++) //crsi increments
    {
      if (ii == 0)
      {
        crsi = crs1 + (crs2 - crs1) *.2 * ((double)i);
        if (i == 0)
          err |= invDist(center, geo1Start, &disti, EPS);
        else
          err |= crsIntx(center, crsi, &tempd1, &disti, geo1Start, geo1->startAz, &tempd2, &tempd3, &intx, TOL, EPS);
      }
      else if (ii == 1)
      {
        crsi = crs2 + (crs3 - crs2) *.2 * ((double)i);
        if (i == 0)
          err |= invDist(center, geo2Start, &disti, EPS);
        else
          err |= crsIntx(center, crsi, &tempd1, &disti, geo2Start, geo2->startAz, &tempd2, &tempd3, &intx, TOL, EPS);
      }
      else if (ii == 2)
      {
        crsi = crs3 + (crs4 - crs3) *.2 * ((double)i);
        if (i == 0)
          err |= invDist(center, geo3Start, &disti, EPS);
        else
          err |= crsIntx(center, crsi, &tempd1, &disti, geo3Start, geo3->startAz, &tempd2, &tempd3, &intx, TOL, EPS);
      }
      else if (ii == 3)
      {
        crsi = crs4 + modcrs(crs1 - crs4) *.2 * ((double)i);
        if (i == 0)
          err |= invDist(center, geo4Start, &disti, EPS);
        else
          err |= crsIntx(center, crsi, &tempd1, &disti, geo4Start, geo4->startAz, &tempd2, &tempd3, &intx, TOL, EPS);
      }

      for (j = 1; j < 7; j++) //distj increments along each crsi
      {
        err = 0;
        distj = disti * .2 * ((double)j);
        if (j == 6)
          distj = disti - 1.2 * TOL;
        err |= direct(center, crsi, distj, &testPoint, EPS);

//        inside = ptsAreInsideBndry(b1, testPoint, &err, TOL, EPS);
//        pointsArr[0]  = testPoint;
//        err |= ptsAreInsideBndry(b1, pointsArr, insideArr, 1, TOL, EPS);
        latList[0] = testPoint.latitude;
        lonList[0] = testPoint.longitude;
		err |= ptsAreInsideBndry(b1, latList, lonList, insideArr, 1, TOL, EPS);
        inside = insideArr[0];

        if (err) {
          printf("Error(s) occurred during inside %d ii = %d i = %d j = %d: 0x%lx", jj,ii,i,j,err);
          errorCount++;
          failedCount++;
        }
        else {
          if (inside == insideexp) {
            passedCount++;
          }
          else {
            failedCount++;
            printf("Failed inside jj = %d ii = %d i = %d j = %d\n",jj,ii,i,j);
          }
        }
        testCaseCount++;
      } //for j
    } //for i
    } //for ii

    //Test points relative to (outside) each geodesic of boundary 1
    insideexp = 0;
    for (ii = 0; ii < 4; ii++) //each geodesic
    {
    for (i = 0; i < 5; i++) //disti increments along geodesic
    {
      if (ii == 0)
      {
        disti = geo1->length * .2 * ((double)i);
        err |= direct(geo1Start, geo1->startAz, disti, &geoPt, EPS);
        crsAtPt = geoCrs(*geo1, geoPt, &tempd1, &tempd2, &tempd3, &err, TOL, EPS);
        perpCrs = modcrs(crsAtPt - M_PI_2);
      }
      else if (ii == 1)
      {
        disti = geo2->length * .2 * ((double)i);
        err |= direct(geo2Start, geo2->startAz, disti, &geoPt, EPS);
        crsAtPt = geoCrs(*geo2, geoPt, &tempd1, &tempd2, &tempd3, &err, TOL, EPS);
        perpCrs = modcrs(crsAtPt - M_PI_2);
      }
      else if (ii == 2)
      {
        disti = geo3->length * .2 * ((double)i);
        err |= direct(geo3Start, geo3->startAz, disti, &geoPt, EPS);
        crsAtPt = geoCrs(*geo3, geoPt, &tempd1, &tempd2, &tempd3, &err, TOL, EPS);
        perpCrs = modcrs(crsAtPt - M_PI_2);
      }
      else if (ii == 3)
      {
        disti = geo4->length * .2 * ((double)i);
        err |= direct(geo4Start, geo4->startAz, disti, &geoPt, EPS);
        crsAtPt = geoCrs(*geo4, geoPt, &tempd1, &tempd2, &tempd3, &err, TOL, EPS);
        perpCrs = modcrs(crsAtPt - M_PI_2);
      }
      for (j = 0; j < 6; j++) //distj increments from geoPt
      {
        err = 0;
        distj = 1.0 + ((double)j) * 3.0;
        if (j == 5)
          distj = 1.1 * TOL;
        err |= direct(geoPt, perpCrs, distj, &testPoint, EPS);

//        inside = ptsAreInsideBndry(b1, testPoint, &err, TOL, EPS);
//        pointsArr[0]  = testPoint;
//        err |= ptsAreInsideBndry(b1, pointsArr, insideArr, 1, TOL, EPS);
        latList[0] = testPoint.latitude;
        lonList[0] = testPoint.longitude;
		err |= ptsAreInsideBndry(b1, latList, lonList, insideArr, 1, TOL, EPS);
        inside = insideArr[0];

        if (err) {
          printf("Error(s) occurred during outside %d ii = %d i = %d j = %d: 0x%lx", jj,ii,i,j,err);
          errorCount++;
          failedCount++;
        }
        else {
          if (inside == insideexp) {
            passedCount++;
          }
          else {
            printf("Failed outside jj = %d ii = %d i = %d j = %d\n",jj,ii,i,j);
//			if(jj == 0 && ii == 1 && i == 0 && j == 5)
//			{
//				displayDataBndry(b1,"b1",1);
//				displayDataPt(testPoint,"testPt",1);
//				displayMatlabBndry(b1,"b1",0);
//				displayMatlabPt(testPoint,"testPt",0);
//	//				displayBndry(b,"b2",0);
//	//				displayPt(testPoint,"testPt",0);
//				printf("\n%d\n",inside);
//			}
            failedCount++;
          }
        }
        testCaseCount++;
      } //for j
    } //for i
    } //for ii
    } //for jj

    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("Finish testPtsAreInsideBndry_Set2\n\n\n");

    return set;
}

/*
 * NAME: testPtsAreInsideBndry_Set3
 *
 * DESCRIPTION:
 * 		This function is used to test the ptsAreInsideBndry function.
 *
 * 		This function runs the test data created by Rich Snow
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testPtsAreInsideBndry_Set3(TestSet) - A test set with the folling metrics:
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
TestSet testPtsAreInsideBndry_Set3()
{
    double DEG2RAD = M_PI / 180.0;
    //const double RAD2DEG = 180.0 / M_PI;
    int passedCount=0, failedCount=0;
    int errorCount = 0;
    int unverifiedCount = 0;
    int setupFailureCount = 0;
    int testCaseCount = 0;
    int inside, insideexp, i, j, ii, jj;
    ErrorSet err=0;
    long newSeed = 20080523;
    double latS, lonS, innerRadius, outerRadius, startSubAngle, endSubAngle;
    double length, disti, distj, distp, crsi;
    double tempd1, tempd2, tempd3, crsAtPt, perpCrs, subtendedAngle, raddiff;
    Geodesic* beginGeo;
    Geodesic* endGeo;
    Arc* innerArc;
    Arc* outerArc;
    Locus* startLocLeft;
    Locus* startLocRight;
    Locus* endLocLeft;
    Locus* endLocRight;
    LLPoint center, testPoint,  geoPt, ptOnLoc;
    Boundary b;
    double latList[1];
    double lonList[1];
    int insideArr[1];

    TestSet set;
    set = newTestSet("testPtsAreInsideBndry_Set3");

    printf("\n\nStart testPtsAreInsideBndry_Set3\n");

    srand(newSeed);  //Initialize the random number generator

    for (jj = 0; jj < 100; jj++)
    {
      b = createBndry();
      latS = DEG2RAD * randLat();
      lonS = DEG2RAD * randLon();
      startSubAngle = DEG2RAD * randAzimuth();
      endSubAngle = startSubAngle + DEG2RAD * ((double)((rand() % 60) + 30));
      innerRadius = (double)((rand() % 50) + 10);
      outerRadius = innerRadius + (double)((rand() % 10) + 10);
      length = (double)((rand() % 50) + 10);
      //printf("\nlatS %4.8f lonS %4.8f innerRad %4.8f outerRad %4.8f startAngle %4.8f endAngle %4.8f length %4.8f\n",latS*RAD2DEG,lonS*RAD2DEG,innerRadius,outerRadius,startSubAngle*RAD2DEG,endSubAngle*RAD2DEG,length);

      center.latitude = latS;
      center.longitude = lonS;
      err |= construct2Arc2LocusBoundary(center, innerRadius, outerRadius, startSubAngle, endSubAngle, length, &b);
      beginGeo = (Geodesic*) &b.elements[0];
      startLocLeft = (Locus*) &b.elements[1];
      innerArc = (Arc*) &b.elements[2];
      endLocRight = (Locus*) &b.elements[3];
      endGeo = (Geodesic*) &b.elements[4];
      endLocLeft = (Locus*) &b.elements[5];
      outerArc = (Arc*) &b.elements[6];
      startLocRight = (Locus*) &b.elements[7];

    //Test points relative to beginGeo and endGeo of boundary b
    for (ii = 0; ii < 2; ii++) //each geodesic
    {
    for (i = 0; i < 6; i++) //increment disti along geodesic
    {
      insideexp = 0;
      if (ii == 0)
      {
        disti = beginGeo->length * .2 * ((double)i);
        err |= direct(beginGeo->startPoint, beginGeo->startAz, disti, &geoPt, EPS);
        crsAtPt = geoCrs(*beginGeo, geoPt, &tempd1, &tempd2, &tempd3, &err, TOL, EPS);
        perpCrs = modcrs(crsAtPt - M_PI_2);
      }
      else
      {
        disti = endGeo->length * .2 * ((double)i);
        err |= direct(endGeo->startPoint, endGeo->startAz, disti, &geoPt, EPS);
        crsAtPt = geoCrs(*endGeo, geoPt, &tempd1, &tempd2, &tempd3, &err, TOL, EPS);
        perpCrs = modcrs(crsAtPt - M_PI_2);
      }
      for (j = 0; j < 6; j++) //increment distj from geoPt
      {
        err = 0;
        distj = 1.0 + ((double)j) * 3.0;
        if (j == 4)
          distj = 1.1 * TOL;
        else if (j == 5)
        {
          distj = -1.1 * TOL;
          insideexp = 1;
        }
        err |= direct(geoPt, perpCrs, distj, &testPoint, EPS);

//        inside = ptsAreInsideBndry(b, testPoint, &err, TOL, EPS);
//        pointsArr[0]  = testPoint;
//        err |= ptsAreInsideBndry(b, pointsArr, insideArr, 1, TOL, EPS);
		latList[0] = testPoint.latitude;
		lonList[0] = testPoint.longitude;
		err |= ptsAreInsideBndry(b, latList, lonList, insideArr, 1, TOL, EPS);
        inside = insideArr[0];

        if (err) {
          printf("\nError(s) occurred during line %d ii = %d i = %d j = %d: 0x%lx", jj,ii,i,j,err);
          errorCount++;
          failedCount++;
        }
        else {
          if (inside == insideexp) {
            passedCount++;
          }
          else {
            printf("Failed line jj = %d ii = %d i = %d j = %d\n",jj,ii,i,j);
            failedCount++;
          }
        }
        testCaseCount++;
      } //for j
    } //for i
    } //for ii

    //Test points relative to startLocLeft and startLocRight
    for (ii = 0; ii < 2; ii++)  // each locus
    {
    for (i = 0; i < 6; i++) //increment disti along geodesic
    {
      insideexp = 1;
      disti = startLocLeft->geoLength * .2 * ((double)i);
      if (i == 5)
        disti = .99 * startLocLeft->geoLength;
      err |= direct(startLocLeft->geoStart, startLocLeft->geoAz, disti, &geoPt, EPS);
      if (ii == 0)
        err |= ptOnLocusFromGeoPt(*startLocLeft, geoPt, &ptOnLoc, &perpCrs, TOL,EPS);
      else
        err |= ptOnLocusFromGeoPt(*startLocRight, geoPt, &ptOnLoc, &perpCrs, TOL,EPS);
      err |= invDist(geoPt, ptOnLoc, &distp, EPS);
      for (j = 0;  j < 9; j++) //increment distj from geoPt
      {
        err = 0;
        distj = distp * .2 * ((double)j);
        if (j == 6)
          distj = distp - 1.1 * TOL;
        else if (j == 8)
          distj = distp + 1.1 * TOL;
        err |= direct(geoPt, perpCrs, distj, &testPoint, EPS);
        if (j >= 7)
          insideexp = 0;

//        inside = ptsAreInsideBndry(b, testPoint, &err, TOL, EPS);
//        pointsArr[0]  = testPoint;
//        err |= ptsAreInsideBndry(b, pointsArr, insideArr, 1, TOL, EPS);
        latList[0] = testPoint.latitude;
        lonList[0] = testPoint.longitude;
		err |= ptsAreInsideBndry(b, latList, lonList, insideArr, 1, TOL, EPS);
        inside = insideArr[0];

        if (err) {
          printf("\nError(s) occurred during startloc jj = %d ii = %d i = %d j = %d: 0x%lx", jj,ii,i,j,err);
//          			if(jj == 97 && ii == 1 && i == 0 && j == 6)
//          			{
//          				displayDataBndry(b,"b1",1);
//          				displayDataPt(testPoint,"testPt",1);
//          				displayMatlabBndry(b,"b1",0);
//          				displayMatlabPt(testPoint,"testPt",0);
//          	//				displayBndry(b,"b2",0);
//          	//				displayPt(testPoint,"testPt",0);
//          				printf("\n%d\n",inside);
//          			}
          errorCount++;
          failedCount++;
        }
        else {
          if (inside == insideexp) {
            passedCount++;
          }
          else {
            printf("Failed startLoc jj = %d ii = %d i = %d j = %d\n",jj,ii,i,j);
            failedCount++;
          }
        }
        testCaseCount++;
      } //for j
    } //for i
    } //for ii

    //Test points relative to endLocLeft and endLocRight
    for (ii = 0; ii < 2; ii++)  // each locus
    {
    for (i = 0; i < 6; i++) //increment disti along geodesic
    {
      insideexp = 1;
      disti = endLocLeft->geoLength * .2 * ((double)i);
      if (i == 5)
        disti = .99 * endLocLeft->geoLength;
      err |= direct(endLocLeft->geoStart, endLocLeft->geoAz, disti, &geoPt, EPS);
      if (ii == 0)
        err |= ptOnLocusFromGeoPt(*endLocLeft, geoPt, &ptOnLoc, &perpCrs, TOL,EPS);
      else
        err |= ptOnLocusFromGeoPt(*endLocRight, geoPt, &ptOnLoc, &perpCrs, TOL,EPS);
      err |= invDist(geoPt, ptOnLoc, &distp, EPS);
      for (j = 0;  j < 9; j++) //increment distj from geoPt
      {
        err = 0;
        distj = distp * .2 * ((double)j);
        if (j == 6)
          distj = distp -1.1 * TOL;
        else if (j == 8)
          distj = distp + 1.1 * TOL;
        err |= direct(geoPt, perpCrs, distj, &testPoint, EPS);
        if (j >= 7)
          insideexp = 0;

//        inside = ptsAreInsideBndry(b, testPoint, &err, TOL, EPS);
//        pointsArr[0]  = testPoint;
//        err |= ptsAreInsideBndry(b, pointsArr, insideArr, 1, TOL, EPS);
        latList[0] = testPoint.latitude;
        lonList[0] = testPoint.longitude;
		err |= ptsAreInsideBndry(b, latList, lonList, insideArr, 1, TOL, EPS);
        inside = insideArr[0];

        if (err) {
          printf("\nError(s) occurred during endloc %d ii = %d i = %d j = %d: 0x%lx", jj,ii,i,j,err);
          errorCount++;
          failedCount++;
        }
        else {
          if (inside == insideexp) {
            passedCount++;
          }
          else {
            printf("Failed endLoc jj = %d ii = %d i = %d j = %d\n",jj,ii,i,j);
            failedCount++;
          }
        }
        testCaseCount++;
      } //for j
    } //for i
    } //for ii

    //Test points relative to the arcs
      subtendedAngle = computeSubtendedAngle(startSubAngle, endSubAngle, CLOCKWISE);
      raddiff = outerRadius - innerRadius;
    for (i = 0; i < 5; i++) //increment crsi from arc center
    {
      insideexp = 1;
      crsi = startSubAngle  + subtendedAngle * .2 * ((double)i);
      for (j = 0;  j < 10; j++) //increment distj along crsi
      {
        err = 0;
        distj = innerRadius + raddiff * .2 * ((double)j);
        if (j == 6)
          distj = outerRadius - 1.1 * TOL;
        else if (j == 7)
          distj = innerRadius + 1.1 * TOL;
        else if (j == 8)
          distj = outerRadius + 1.1 * TOL;
        else if (j == 9)
          distj = innerRadius - 1.1 * TOL;
        err |= direct(center, crsi, distj, &testPoint, EPS);
        if (j >= 8)
          insideexp = 0;

//        inside = ptsAreInsideBndry(b, testPoint, &err, TOL, EPS);
//        pointsArr[0]  = testPoint;
//        err |= ptsAreInsideBndry(b, pointsArr, insideArr, 1, TOL, EPS);
		latList[0] = testPoint.latitude;
		lonList[0] = testPoint.longitude;
		err |= ptsAreInsideBndry(b, latList, lonList, insideArr, 1, TOL, EPS);
        inside = insideArr[0];

        if (err) {
          errorCount++;
          failedCount++;
        }
        else {
          if (inside == insideexp) {
            passedCount++;
          }
          else {
            failedCount++;
            printf("Failed arc jj = %d i = %d j = %d\n",jj,i,j);
          }
        }
        testCaseCount++;
      }//for j
    } //for i

    } //for jj

    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("Finish testPtsAreInsideBndry_Set3\n\n\n");

    return set;
}

/*
 * NAME: testPtsAreInsideBndry_Set4
 *
 * DESCRIPTION:
 * 		This function is used to test the ptsAreInsideBndry function.
 *
 * 		This function runs the test data created by Joe Heidrich
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testPtsAreInsideBndry_Set4(TestSet) - A test set with the folling metrics:
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
TestSet testPtsAreInsideBndry_Set4()
{
    int passedCount=0, failedCount=0;
    int errorCount = 0;
    int unverifiedCount = 0;
    int setupFailureCount = 0;
    int testCaseCount = 0;
    double latxExp, lonxExp, len, lineAz, startRad, endRad, az, rad, subAngle;
    LLPoint sp, tempPt1, testPt;
    Geodesic line1, line2;
    double az12, dist;
    int i, j, inside, insideexp;
    ErrorSet err = 0;

    Boundary b;
    double latList[1];
    double lonList[1];
    int insideArr[1];

    double DEG2RAD = M_PI / 180;

    Locus* outLoc1;
    Locus* outLoc2;
    Spiral* spiral1;
    Spiral* spiral2;

    double eps = 1.0e-20;
    double tol = 1.37e-9;

    //Percentage of distances to step through
    double distArray[7] = {0, .5, .99, .999999, 1.000001, 1.01, 1.5};
    TestSet set;
    set = newTestSet("testPtsAreInsideBndry_Set4");

    printf("\n\nStart testPtsAreInsideBndry_Set4\n");

    srand(04012011);

	for (i=0;i<100;i++) {

                err = 0;
		b = createBndry();

		latxExp = randLat();
		lonxExp = randLon();

		err |= createPt(&sp, latxExp*DEG2RAD, lonxExp*DEG2RAD);
		startRad = 30 + 0.01 * randDist() / 2; //0-54 nm
		endRad = 30 + 0.01 * randDist() / 2; //0-54 nm
		lineAz = randAzimuth();
		len = 120 + .01 * randDist();

                err |= createLocusSpiralBndry(sp, lineAz, len, startRad, endRad, &b);

//                displayMatlabBndry(b,"b",0);

                outLoc1 = (Locus*)&b.elements[0];
                outLoc2 = (Locus*)&b.elements[1];
                spiral1 = (Spiral*)&b.elements[2];
                spiral2 = (Spiral*)&b.elements[3];

		//Create Points that are inside/outside
		err |= createGeo(&line1, outLoc1->geoStart, outLoc1->geoEnd, SEGMENT, eps);
		err |= createGeo(&line2, outLoc2->geoStart, outLoc2->geoEnd, SEGMENT, eps);

                if (err) {
                  printf("\nError(s) occurred during setup i = %d: 0x%lx", i,err);
                  continue;
                }

                if(i==99) continue; //exclude polar case

		for (j=0;j<7;j++) {
                        err = 0;
			dist = rand() % (int) floor(line1.length);
			err |= direct(line1.startPoint, line1.startAz, dist, &tempPt1, eps);
			err |= ptOnLocusFromGeoPt(*outLoc1, tempPt1, &testPt, &az12, tol, eps);
			err |= direct(tempPt1, az12, distArray[j] * outLoc1->startDist, &testPt, eps);

//			inside = ptsAreInsideBndry(b, testPt, &err, tol, eps);
//	        pointsArr[0]  = testPt;
//	        err |= ptsAreInsideBndry(b, pointsArr, insideArr, 1, TOL, EPS);
			latList[0] = testPt.latitude;
			lonList[0] = testPt.longitude;
			err |= ptsAreInsideBndry(b, latList, lonList, insideArr, 1, TOL, EPS);
	        inside = insideArr[0];

                        if (j <= 3)
                          insideexp = 1;
                        else
                          insideexp = 0;

                        if (err) {
                        	printf("Failed loc1 i = %d j = %d  ",i,j);
                        	printf("Error: %s", formatErrorMessage(err));
                          errorCount++;
                          failedCount++;
                         }
                         else {
                           if (inside == insideexp) {
                             passedCount++;
                           }
                           else {
                             failedCount++;
                             printf("Failed loc1 i = %d j = %d\n",i,j);
                           }
                         }
                         testCaseCount++;
		}// for j

		for (j=0;j<7;j++) {
                        err = 0;
			dist = rand() % (int) floor(line2.length);
			err |= direct(line2.startPoint, line2.startAz, dist, &tempPt1, eps);
			err |= ptOnLocusFromGeoPt(*outLoc2, tempPt1, &testPt, &az12, tol, eps);
			err |= direct(tempPt1, az12, distArray[j] * fabs(outLoc2->startDist), &testPt, eps);

//			inside = ptsAreInsideBndry(b, testPt, &err, tol, eps);
//			pointsArr[0]  = testPt;
//			err |= ptsAreInsideBndry(b, pointsArr, insideArr, 1, TOL, EPS);
			latList[0] = testPt.latitude;
			lonList[0] = testPt.longitude;
			err |= ptsAreInsideBndry(b, latList, lonList, insideArr, 1, TOL, EPS);
			inside = insideArr[0];

                        if (j <= 3)
                          insideexp = 1;
                        else
                          insideexp = 0;

                        if (err) {
                        	printf("Failed loc2 i = %d j = %d  ",i,j);
                        	printf("Error: %s", formatErrorMessage(err));
                          errorCount++;
                          failedCount++;
                         }
                         else {
                           if (inside == insideexp) {
                             passedCount++;
                           }
                           else {
                             failedCount++;
                             printf("Failed loc2 i = %d j = %d\n",i,j);
                           }
                         }
                         testCaseCount++;
		}// for j

		for (j=0;j<7;j++) {
                        err = 0;
			subAngle = spiral1->subtendedAngle / DEG2RAD;
			az = spiral1->startAz + (spiral1->dir * (rand() % (int) floor(subAngle)) * DEG2RAD);
			err |= spiralRadius(*spiral1, az, &rad);
			err |= direct(spiral1->centerPoint, az, distArray[j] * rad, &testPt, eps);

//			inside = ptsAreInsideBndry(b, testPt, &err, tol, eps);
//			pointsArr[0]  = testPt;
//			err |= ptsAreInsideBndry(b, pointsArr, insideArr, 1, TOL, EPS);
			latList[0] = testPt.latitude;
			lonList[0] = testPt.longitude;
			err |= ptsAreInsideBndry(b, latList, lonList, insideArr, 1, TOL, EPS);
			inside = insideArr[0];


//			if(i==99 && j==2){
//				displayDataBndry(b,"b",1);
//				displayDataPt(testPt,"testPt",1);
////				displayMatlabBndry(b,"b",0);
////				displayMatlabPt(testPt,"testPt",0);
////				displayBndry(b,"b",0);
////				displayPt(testPt,"testPt",0);
//				printf("\n%d\n",inside);
//			}

                        if (j <= 3)
                          insideexp = 1;
                        else
                          insideexp = 0;

                        if (err) {
                        	printf("Failed spiral1 i = %d j = %d  ",i,j);
                        	printf("Error: %s", formatErrorMessage(err));
                          errorCount++;
                          failedCount++;
                         }
                         else {
                           if (inside == insideexp) {
                             passedCount++;
                           }
                           else {
                             failedCount++;
                             printf("Failed spiral1 i = %d j = %d\n",i,j);
                           }
                         }
                         testCaseCount++;
		}// for j

		for (j=0;j<7;j++) {
                        err = 0;
			subAngle = spiral2->subtendedAngle / DEG2RAD;
			az = spiral2->startAz + (spiral2->dir * (rand() % (int) floor(subAngle))) * DEG2RAD;
			err |= spiralRadius(*spiral2, az, &rad);
			err |= direct(spiral2->centerPoint, az, distArray[j] * rad, &testPt, eps);

//			inside = ptsAreInsideBndry(b, testPt, &err, tol, eps);
//			pointsArr[0]  = testPt;
//			err |= ptsAreInsideBndry(b, pointsArr, insideArr, 1, TOL, EPS);
			latList[0] = testPt.latitude;
			lonList[0] = testPt.longitude;
			err |= ptsAreInsideBndry(b, latList, lonList, insideArr, 1, TOL, EPS);
			inside = insideArr[0];

                        if (j <= 3)
                          insideexp = 1;
                        else
                          insideexp = 0;

                        if (err) {
                        	printf("Failed spiral2 i = %d j = %d  ",i,j);
                        	printf("Error: %s", formatErrorMessage(err));
                          errorCount++;
                          failedCount++;
                         }
                         else {
                           if (inside == insideexp) {
                             passedCount++;
                           }
                           else {
                             failedCount++;
                             printf("Failed spiral2 i = %d j = %d\n",i,j);
                           }
                         }
                         testCaseCount++;
		} //for j
	} //for i

    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("Finish testPtsAreInsideBndry_Set4\n\n\n");

    return set;
}

TestSuite testPtsAreInsideBndry_AllSets()
{
	TestSuite suite;
	TestSet set1;
	TestSet set2;
	TestSet set3;
	TestSet set4;

    printf("\nStart testPtsAreInsideBndry_AllSets\n");

    suite = newTestSuite("testPtsAreInsideBndry_AllSets");

    set1 = testPtsAreInsideBndry_Set1();
    addTestSet(set1,&suite);

    set2 = testPtsAreInsideBndry_Set2();
    addTestSet(set2,&suite);

    set3 = testPtsAreInsideBndry_Set3();
    addTestSet(set3,&suite);

    set4 = testPtsAreInsideBndry_Set4();
    addTestSet(set4,&suite);

    displayTestSuite(suite);

    printf("\nFinish testPtsAreInsideBndry_AllSets\n\n\n");

    return suite;
}


/*
 * NAME: testBoundaryGeoIntersection_Set1
 *
 * DESCRIPTION:
 * 		This function is used to test the boundaryGeoIntersection function.
 *
 * 		This function runs the test data created by Blythe Debenport and uses
 * 		her original test plan to execute the actual tests.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testBoundaryGeoIntersection_Set1(TestSet) - A test set with the folling metrics:
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
TestSet testBndryGeoIntx_Set1()
{
    int passedCount=0, failedCount=0;
    int errorCount = 0;
    int unverifiedCount = 0;
    int setupFailureCount = 0;
    int testCaseCount = 0;
    long err=0;

    TestSet set;

    double locus1_startDist, locus1_endDist, locus2_startDist, locus2_endDist;
	LLPoint line1_startPoint, line1_endPoint, line2_startPoint, line2_endPoint;
	LLPoint locus1_startPoint, locus1_endPoint,  locus2_startPoint, locus2_endPoint;
	Geodesic line1, line2;
	Locus locus1, locus2;
	LLPoint testPoint;
	LLPoint intx1;
	double az;

    LLPoint* intxList = NULL;
    int intxCount = 0;
    int intxCode = 0;
    int outputMatlab = 0; //0 = don't output matlab code to the console

	Boundary b1 =  createBndry();

    // Boundary 1
    line1_startPoint.latitude = 39.325230354683285*M_PI/180;
    line1_startPoint.longitude = -77.71817718753329*M_PI/180;
    line1_endPoint.latitude = 39.28975340430264*M_PI/180;
    line1_endPoint.longitude = -77.79092092779719*M_PI/180;
    line2_startPoint.latitude = 39.11317792490197*M_PI/180;
    line2_startPoint.longitude = -77.54755357485374*M_PI/180;
    line2_endPoint.latitude = 39.07780615758011*M_PI/180;
    line2_endPoint.longitude = -77.62016459538935*M_PI/180;
    locus1_startPoint.latitude = 39.30749758814272*M_PI/180;
    locus1_startPoint.longitude = -77.7545582397233*M_PI/180;
    locus1_endPoint.latitude = 39.09549772052336*M_PI/180;
    locus1_endPoint.longitude = -77.58386815425973*M_PI/180;
    locus2_startPoint.latitude = 39.30749758814272*M_PI/180;
    locus2_startPoint.longitude = -77.7545582397233*M_PI/180;
    locus2_endPoint.latitude = 39.09549772052336*M_PI/180;
    locus2_endPoint.longitude = -77.58386815425973*M_PI/180;

    locus1_startDist = -2.0;
    locus1_endDist = -2.0;
    locus2_startDist = 2.0;
    locus2_endDist = 2.0;

    printf("\n\nStart testBoundaryGeoIntersection_Set1\n");

    set = newTestSet("testBndryGeoIntx_Set1");

	// create the Geodesic
	err = createGeo(&line1,  line1_startPoint, line1_endPoint, SEGMENT, EPS);
	err = createGeo(&line2,  line2_startPoint, line2_endPoint, SEGMENT, EPS);

	//  create the loci
	err = createLocus(&locus1,  locus1_startPoint, locus1_endPoint, locus1_startDist, locus1_endDist, SEGMENT, TOL, EPS);
	err = createLocus(&locus2,  locus2_startPoint, locus2_endPoint, locus2_startDist, locus2_endDist, SEGMENT, TOL, EPS);
        if (outputMatlab) displayMatlabGeo(line1, "line1", 0);
        if (outputMatlab) displayMatlabGeo(line2, "line2", 0);
        if (outputMatlab) displayMatlabLocus(locus1, "locus1", 0);
        if (outputMatlab) displayMatlabLocus(locus2, "locus2", 0);

	// create  boundry1
	err = addGeoToBndry(&b1, &line1);
	err = addGeoToBndry(&b1, &line2);
	err = addLocusToBndry(&b1, &locus1);
	err = addLocusToBndry(&b1, &locus2);

        //if (outputMatlab) displayMatlabBndry(b1, "bndy1", 0);

	// Test Case 9 -- Point is inside boundary and test geo goes
        // through a vertex of the boundary.
	testPoint.latitude = 39.28875771541921*M_PI/180;
	testPoint.longitude = -77.71038254006952*M_PI/180;
	az = 270.9369627534785*M_PI/180;

	intx1.latitude = 39.28975340430235*M_PI/180;
	intx1.longitude = -77.79092092779698*M_PI/180;

	err = 0;
	err = bndryGeoIntx(b1, testPoint, az, &intxList, &intxCount, &intxCode, TOL, EPS);

	if (err) {
		errorCount++;
		failedCount++;
	}
	else {
                //intxCode 2 means the test geo passes through a vertex
		if (intxCount == 2 && intxCode == 2) {
			if (ptsAreSame(intx1, intxList[0], TOL) && ptsAreSame(intx1, intxList[1], TESTTOL)) {
				//printf("Test 9: PASSED\n");
				passedCount++;
			}
			else {
				printf("Test 9: FAILED Intersection\n");
				failedCount++;
			}
		}
		else {
			printf("Test 9: FAILED Count\n");
			failedCount++;
		}
	}
	testCaseCount++;

    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("Finish testBndryGeoIntx_Set1\n\n\n");

    return set;
}

/*
 * NAME: testBndryGeoIntx_Set2
 *
 * DESCRIPTION:
 * 		This function is used to test the bndryGeoIntx function.
 *
 * 		This function runs the test data created by Richard Snow.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testBndryGeoIntx_Set2(TestSet) - A test set with the folling metrics:
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
TestSet testBndryGeoIntx_Set2()
{
    double DEG2RAD = M_PI / 180.0;
    //const double RAD2DEG = 180.0 / M_PI;
    int passedCount=0, failedCount=0;
    int errorCount = 0;
    int unverifiedCount = 0;
    int setupFailureCount = 0;
    int testCaseCount = 0;
    int outputMatlab = 0; //0 = don't output matlab code to the console
    int i, j, ii, jj;
    ErrorSet err=0;
    long newSeed = 20080523;
    double latS, lonS, crs1, crs2, crs3, crs4, crsi, radius, disti;
    double tempd1, tempd2, tempd3, crsj, tempcrs;
    LineType lineType;
    Geodesic* geo1;
    Geodesic* geo2;
    Geodesic* geo3;
    Geodesic* geo4; 
    LLPoint center, geo1Start, geo2Start, geo3Start, geo4Start, intx, intx1;
    Boundary b1;
    LLPoint* intxList = NULL;
    int intxCount = 0, intxCountExp = 0;
    int intxCode = 0;

    TestSet set;
    set = newTestSet("testBndryGeoIntx_Set2");

    printf("\n\nStart testBndryGeoIntx_Set2\n");

    srand(newSeed);  //Initialize the random number generator
 
    lineType = SEGMENT;

    for (jj = 0; jj < 100; jj++)
    {
    b1 = createBndry();

    latS = DEG2RAD * randLat();
    lonS = DEG2RAD * randLon();
    crs1 = DEG2RAD * ((double)((rand() % 89) + 1));
    crs2 = DEG2RAD * ((double)((rand() % 89) + 91));
    crs3 = DEG2RAD * ((double)((rand() % 89) + 181));
    crs4 = DEG2RAD * ((double)((rand() % 89) + 271));
    radius = (double)((rand() % 100) + 1);

    //printf("\nlatS %4.8f lonS %4.8f crs1 %4.8f crs2 %4.8f crs3 %4.8f crs4 %4.8f rad %4.8f\n",latS*RAD2DEG,lonS*RAD2DEG,crs1*RAD2DEG,crs2*RAD2DEG,crs3*RAD2DEG,crs4*RAD2DEG,radius);
    //boundary 1
    center.latitude = latS;
    center.longitude = lonS;
    err |= constructGeoBoundary(center, radius, crs1, crs2, crs3, crs4, &b1);
    geo1 = (Geodesic*) &b1.elements[0];
    geo2 = (Geodesic*) &b1.elements[1];
    geo3 = (Geodesic*) &b1.elements[2];
    geo4 = (Geodesic*) &b1.elements[3];
    geo1Start = geo1->startPoint;
    geo2Start = geo2->startPoint;
    geo3Start = geo3->startPoint;
    geo4Start = geo4->startPoint;
    if (err) printf("\nError(s) occurred during setup: 0x%lx", err);

        if (outputMatlab) displayMatlabGeo(*geo1, "geo1", 0);
        if (outputMatlab) displayMatlabGeo(*geo2, "geo2", 0);
        if (outputMatlab) displayMatlabGeo(*geo3, "geo3", 0);
        if (outputMatlab) displayMatlabGeo(*geo4, "geo4", 0);

    //Test line with 1 intersection on side and two at vertex
    for (ii = 0; ii < 4; ii++) //each geodesic
    {
    for (i = 0; i < 5; i++) //crsi increments
    {
      err = 0;
      intxCountExp = 1;
      if (ii == 0)
      {
        crsi = crs1 + (crs2 - crs1) *.2 * ((double)i);
        if (i == 0)
        {
          intx.latitude = geo1Start.latitude;
          intx.longitude = geo1Start.longitude;
          intxCountExp = 2;
        }
        else
          err |= crsIntx(center, crsi, &tempd1, &disti, geo1Start, geo1->startAz, &tempd2, &tempd3, &intx, TOL, EPS);
      }
      else if (ii == 1)
      {
        crsi = crs2 + (crs3 - crs2) *.2 * ((double)i);
        if (i == 0)
        {
          intx.latitude = geo2Start.latitude;
          intx.longitude = geo2Start.longitude;
          intxCountExp = 2;
        }
        else
          err |= crsIntx(center, crsi, &tempd1, &disti, geo2Start, geo2->startAz, &tempd2, &tempd3, &intx, TOL, EPS);
      }
      else if (ii == 2)
      {
        crsi = crs3 + (crs4 - crs3) *.2 * ((double)i);
        if (i == 0)
        {
          intx.latitude = geo3Start.latitude;
          intx.longitude = geo3Start.longitude;
          intxCountExp = 2;
        }
        else
          err |= crsIntx(center, crsi, &tempd1, &disti, geo3Start, geo3->startAz, &tempd2, &tempd3, &intx, TOL, EPS);
      }
      else if (ii == 3)
      {
        crsi = crs4 + modcrs(crs1 - crs4) *.2 * ((double)i);
          err |= invDist(center, geo4Start, &disti, EPS);
        if (i == 0)
        {
          intx.latitude = geo4Start.latitude;
          intx.longitude = geo4Start.longitude;
          intxCountExp = 2;
        }
        else
          err |= crsIntx(center, crsi, &tempd1, &disti, geo4Start, geo4->startAz, &tempd2, &tempd3, &intx, TOL, EPS);
      }

        err |= bndryGeoIntx(b1, center, crsi, &intxList, &intxCount, &intxCode, TOL, EPS);
        if (err) {
          errorCount++;
          failedCount++;
        }
        else {
          if (intxCount == intxCountExp) {
            if ((intxCount == 2) && ptsAreSame(intx, intxList[0], TESTTOL) && ptsAreSame(intx, intxList[1], TESTTOL)) {
              passedCount++;
            }
            else if ((intxCount == 1) && ptsAreSame(intx, intxList[0], TESTTOL)) {
              passedCount++;
            }
            else {
              failedCount++;
              printf("Failed intx point jj = %d ii = %d i = %d\n",jj,ii,i);
            }
          }
          else {
            failedCount++;
            printf("Failed intx count jj = %d ii = %d i = %d\n",jj,ii,i);
          }
        }
        testCaseCount++;
    } //for i
    } //for ii

    //Test line intersecting two sides
    for (i = 0; i < 5; i++) //crsi increments
    {
      err = 0;
      intxCountExp = 2;
      crsi = crs1 + (crs2 - crs1) *.2 * ((double)i);
      if (i == 0)
      {
        intx.latitude = geo1Start.latitude;
        intx.longitude = geo1Start.longitude;
        intxCountExp = 3;
      }
      else
        err |= crsIntx(center, crsi, &tempd1, &disti, geo1Start, geo1->startAz, &tempd2, &tempd3, &intx, TOL, EPS);

      for (j = 1; j < 6; j++)
      {
        crsj = crs2 + (crs3 - crs2) *.2 * ((double)j);
        if (j == 5)
        {
          intx1.latitude = geo3Start.latitude;
          intx1.longitude = geo3Start.longitude;
          intxCountExp = 3;
          if (i == 0)
            intxCountExp = 4;
        }
        else
          err |= crsIntx(center, crsj, &tempd1, &disti, geo2Start, geo2->startAz, &tempd2, &tempd3, &intx1, TOL, EPS);
          err |= invCrs(intx, intx1, &tempcrs, &tempd1, EPS);

        err |= bndryGeoIntx(b1, intx, tempcrs, &intxList, &intxCount, &intxCode, TOL, EPS);
        if (err) {
          errorCount++;
          failedCount++;
        }
        else {
          if (intxCount == intxCountExp) {
            if ((intxCount == 4) && ptsAreSame(intx, intxList[0], TESTTOL) && ptsAreSame(intx1, intxList[1], TESTTOL) && ptsAreSame(intx1, intxList[2], TESTTOL) && ptsAreSame(intx, intxList[3],TESTTOL)) {
              passedCount++;
            }
            else if (intxCount == 3) {
              if (i == 0  && j < 5 && ptsAreSame(intx, intxList[0], TESTTOL) && ptsAreSame(intx1, intxList[1], TESTTOL) && ptsAreSame(intx, intxList[2], TESTTOL)) 
                passedCount++;
              else if (i > 0  && j == 5 && ptsAreSame(intx, intxList[0], TESTTOL) && ptsAreSame(intx1, intxList[1], TESTTOL) && ptsAreSame(intx1, intxList[2], TESTTOL)) 
                passedCount++;
            }
            else if ((intxCount == 2) && ptsAreSame(intx, intxList[0], TESTTOL) && ptsAreSame(intx1, intxList[1], TESTTOL)) {
              passedCount++;
            }
            else {
              failedCount++;
              printf("Failed intx point jj = %d ii = %d i = %d\n",jj,ii,i);
            }
          }
          else {
            failedCount++;
            printf("Failed intx count jj = %d ii = %d i = %d\n",jj,ii,i);
          }
        }
        testCaseCount++;
     }//for j
    } //for i

    } //for jj

    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("Finish testBndryGeoIntx_Set2\n\n\n");

    return set;
}

/*
 * NAME: testBndryGeoIntx_Set3
 *
 * DESCRIPTION:
 * 		This function is used to test the bndryGeoIntx function.
 *
 * 		This function runs the test data created by Richard Snow.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testBndryGeoIntx_Set3(TestSet) - A test set with the folling metrics:
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
TestSet testBndryGeoIntx_Set3()
{
    double DEG2RAD = M_PI / 180.0;
    //const double RAD2DEG = 180.0 / M_PI;
    int passedCount=0, failedCount=0;
    int errorCount = 0;
    int unverifiedCount = 0;
    int setupFailureCount = 0;
    int testCaseCount = 0;
    int outputMatlab = 0; //0 = don't output matlab code to the console
    int i, j, jj;
    ErrorSet err=0;
    long newSeed = 20080523;
    double latS, lonS, crs, geolen, disti, distj, perpCrs;
    double tempd1, locDist;
    double crsst, crsi, arcradius, crs34, crsarc, crsj;
    LineType lineType;
    LLPoint geoStart, geoEnd, geoPt, intx, intx1;
    LLPoint ptOnLoc3, ptOnLoc4;
    Locus* loc3; 
    Locus* loc4;
    Arc* arc1;
    Arc* arc2;
    Boundary b2;
    LLPoint* intxList = NULL;
    int intxCount = 0, intxCountExp = 0, intxCode = 0;

    TestSet set;
    set = newTestSet("testBndryGeoIntx_Set3");

    printf("\n\nStart testBndryGeoIntx_Set3\n");

    srand(newSeed);  //Initialize the random number generator

    lineType = SEGMENT;
    for (jj = 0; jj < 100; jj++)
    {
    b2 = createBndry();

    latS = DEG2RAD * randLat();
    lonS = DEG2RAD * randLon();
    crs = DEG2RAD * randAzimuth();
    geolen = (double)((rand() % 100) + 1);
    locDist = (double)((rand() % 20) + 1);

    //printf("\nlatS %4.8f lonS %4.8f crs %4.8f geolen %4.8f locDist %4.8f\n",latS*RAD2DEG,lonS*RAD2DEG,crs*RAD2DEG,geolen,locDist);
    geoStart.latitude = latS;
    geoStart.longitude = lonS;

    err |= direct(geoStart, crs, geolen, &geoEnd, EPS);
    //boundary 2
    err |= constructLocusArcBoundary(geoStart, crs, geolen, locDist, &b2);
    loc3 = (Locus*) &b2.elements[0];
    arc1 = (Arc*) &b2.elements[1];
    loc4  = (Locus*) &b2.elements[2];
    arc2 = (Arc*) &b2.elements[3];

    if (err) printf("\nError(s) occurred during setup: 0x%lx", err);

        if (outputMatlab) displayMatlabLocus(*loc3, "loc3", 0);
        if (outputMatlab) displayMatlabArc(*arc1, "arc1", 0);
        if (outputMatlab) displayMatlabLocus(*loc4, "loc4", 0);
        if (outputMatlab) displayMatlabArc(*arc2, "arc2", 0);

    //Test lines relative to loc3 and loc4 of boundary 2
    for (i = 0; i < 6; i++) //increment disti along geodesic of loci
    {
      err = 0;
      intxCountExp = 2;
      if (i == 0)
        intxCountExp = 3;
      disti = geolen * .2 * ((double)i);
      if (i == 5)
        disti = .99 * geolen;
      err |= direct(geoStart, crs, disti, &geoPt, EPS);
      err |= ptOnLocusFromGeoPt(*loc3, geoPt, &ptOnLoc3, &perpCrs, TOL,EPS);
      for (j = 1; j < 6; j++) //increment distj along geodesic of loci
      {
        distj = geolen * .2 * ((double)j);
        if (j == 5)
          distj = .99 * geolen;
        err |= direct(geoStart, crs, distj, &geoPt, EPS);
        err |= ptOnLocusFromGeoPt(*loc4, geoPt, &ptOnLoc4, &perpCrs, TOL,EPS);
        err |= invCrs(ptOnLoc3, ptOnLoc4, &crs34, &tempd1, EPS);
        intx.latitude = ptOnLoc3.latitude;
        intx.longitude = ptOnLoc3.longitude;
        intx1.latitude = ptOnLoc4.latitude;
        intx1.longitude = ptOnLoc4.longitude;
        err |= bndryGeoIntx(b2, ptOnLoc3, crs34, &intxList, &intxCount, &intxCode, TOL, EPS);

        if (err) {
          errorCount++;
          failedCount++;
        }
        else {
          if (intxCount == intxCountExp) {
            if (intxCount == 3 && ptsAreSame(intx, intxList[0], TESTTOL) && ptsAreSame(intx, intxList[1], TESTTOL) && ptsAreSame(intx1, intxList[2], TESTTOL)) 
              passedCount++;
            else if (intxCount == 2 && ptsAreSame(intx, intxList[0], TESTTOL) && ptsAreSame(intx1, intxList[1], TESTTOL)) 
              passedCount++;
            else if (intxCount == 1 && ptsAreSame(intx, intxList[0], TESTTOL)) 
              passedCount++;
            else {
              printf("Failed loc intx point jj = %d i = %d j = %d\n",jj,i,j);
              failedCount++;
            }  
          }
          else {
            printf("Failed loc intx count jj = %d i = %d j = %d intxCount = %d intxCountExp %d\n",jj,i,j,intxCount,intxCountExp);
            failedCount++;
          }
        }
        testCaseCount++;
     } //for j
    } //for i

    //Test lines relative to the arc1 of boundary 2
    err |= invCrs(geoStart, loc3->locusStart, &crsst, &tempd1, EPS);
    arcradius = arc1->radius;
    for (i = 0; i < 6; i++) //increment crsi from arc center
    {
      err = 0;
      intxCountExp = 2;
      if (i == 0)
        intxCountExp = 3;
      crsi = crsst  + M_PI_2 * .2 * ((double)i);
      err |= direct(geoStart, crsi, arcradius, &intx, EPS);
      for (j = 1; j < 6; j++) //increment crsj from arc center
      {
        crsj = crsst + M_PI  - M_PI_2 * .2 * ((double)j);
        err |= direct(geoStart, crsj, arcradius, &intx1, EPS);
        err |= invCrs(intx, intx1, &crsarc, &tempd1, EPS);
        if (i == 5 && j == 5)
        {
          err |= invCrs(intx, geoStart, &crsarc, &tempd1, EPS);
          crsarc = crsarc + M_PI_2;
          intxCountExp = 1;
        }
        err |= bndryGeoIntx(b2, intx, crsarc, &intxList, &intxCount, &intxCode, TOL, EPS);

        if (err) {
          errorCount++;
          failedCount++;
        }
        else {
          if (intxCount == intxCountExp) {
            if (intxCount == 3 && ptsAreSame(intx, intxList[0], TESTTOL) && ptsAreSame(intx, intxList[1], TESTTOL) && ptsAreSame(intx1, intxList[2], TESTTOL)) 
              passedCount++;
            else if (intxCount == 2 && ptsAreSame(intx, intxList[0], TESTTOL) && ptsAreSame(intx1, intxList[1], TESTTOL)) 
              passedCount++;
            else if (intxCount == 1 && ptsAreSame(intx, intxList[0], TESTTOL)) 
              passedCount++;
            else {
              failedCount++;
              printf("Failed arc intx point jj = %d i = %d j = %d\n",jj,i,j);
            }
          }
          else {
            failedCount++;
            printf("Failed arc intx count jj = %d i = %d j = %d intxCount = %d intxCountExp %d\n",jj,i,j,intxCount,intxCountExp);
          }
        }
        testCaseCount++;
      } //for j
    } //for i
    } //for jj

    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("Finish testBndryGeoIntx_Set3\n\n\n");

    return set;
}

/*
 * NAME: testBndryGeoIntx_Set4
 *
 * DESCRIPTION:
 * 		This function is used to test the bndryGeoIntx function.
 *
 * 		This function runs the test data created by Richard Snow.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testBndryGeoIntx_Set4(TestSet) - A test set with the folling metrics:
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
TestSet testBndryGeoIntx_Set4()
{
    double DEG2RAD = M_PI / 180.0;
    //const double RAD2DEG = 180.0 / M_PI;
    int passedCount=0, failedCount=0;
    int errorCount = 0;
    int unverifiedCount = 0;
    int setupFailureCount = 0;
    int testCaseCount = 0;
    int outputMatlab = 0; //0 = don't output matlab code to the console
    int i, j, jj;
    ErrorSet err=0;
    long newSeed = 20080523;
    double latS, lonS, crs, len, disti, distj, perpCrs;
    double tempd1, geo3Az, geo4Az, geo3len, geo4len, startRad, endRad;
    double crsst, crsi, crs34, crsarc, crsj;
    LineType lineType;
    LLPoint sp, geo3Start, geo4Start, geoPt, intx, intx1;
    LLPoint ptOnLoc3, ptOnLoc4;
    Locus* loc3; 
    Locus* loc4;
    Spiral* spiral1;
    Spiral* spiral2;
    Boundary b2;
    LLPoint* intxList = NULL;
    int intxCount = 0, intxCountExp = 0, intxCode = 0;

    TestSet set;
    set = newTestSet("testBndryGeoIntx_Set4");

    printf("\n\nStart testBndryGeoIntx_Set4\n");

    srand(newSeed);  //Initialize the random number generator

    lineType = SEGMENT;
    for (jj = 0; jj < 100; jj++)
    {
    b2 = createBndry();

    latS = DEG2RAD * randLat();
    lonS = DEG2RAD * randLon();
    crs = DEG2RAD * randAzimuth();
    len = 120.0 + 0.01 * randDist();
    startRad = 30.0 + 0.01 * randDist() / 2.0; 
    endRad = 30.0 + 0.01 * randDist() / 2.0; 

    //printf("\nlatS %4.8f lonS %4.8f crs %4.8f len %4.8f startRad %4.8f endRad %4.8f\n",latS*RAD2DEG,lonS*RAD2DEG,crs*RAD2DEG,len,startRad,endRad);
    sp.latitude = latS;
    sp.longitude = lonS;

    //boundary 2
    err |= createLocusSpiralBndry(sp, crs, len, startRad, endRad, &b2);
    loc3 = (Locus*) &b2.elements[0];
    loc4  = (Locus*) &b2.elements[1];
    spiral1 = (Spiral*) &b2.elements[2];
    spiral2 = (Spiral*) &b2.elements[3];

    if (err) printf("\nError(s) occurred during setup: 0x%lx", err);

    geo3Start = loc3->geoStart;
    geo3Az = loc3->geoAz;
    geo3len = loc3->geoLength;
    geo4Start = loc4->geoStart;
    geo4Az = loc4->geoAz;
    geo4len = loc4->geoLength;

        if (outputMatlab) displayMatlabLocus(*loc3, "loc3", 0);
        //if (outputMatlab) displayMatlabArc(*spiral1, "spiral1", 0);
        if (outputMatlab) displayMatlabLocus(*loc4, "loc4", 0);
        //if (outputMatlab) displayMatlabArc(*spiral2, "spiral2", 0);

    //Test lines relative to loc3 and loc4 of boundary 2
    for (i = 0; i < 6; i++) //increment disti along geodesic of loc3
    {
      err = 0;
      intxCountExp = 2;
      if (i == 0)
        intxCountExp = 3;
      disti = geo3len * .2 * ((double)i);
      if (i == 5)
        disti = .99 * geo3len;
      err |= direct(geo3Start, geo3Az, disti, &geoPt, EPS);
      err |= ptOnLocusFromGeoPt(*loc3, geoPt, &ptOnLoc3, &perpCrs, TOL,EPS);
      for (j = 1; j < 6; j++) //increment distj along geodesic of loc4
      {
        distj = geo4len * .2 * ((double)j);
        if (j == 5)
          distj = .99 * geo4len;
        err |= direct(geo4Start, geo4Az, distj, &geoPt, EPS);
        err |= ptOnLocusFromGeoPt(*loc4, geoPt, &ptOnLoc4, &perpCrs, TOL,EPS);
        err |= invCrs(ptOnLoc3, ptOnLoc4, &crs34, &tempd1, EPS);
        intx = ptOnLoc3;
        intx1 = ptOnLoc4;
        err |= bndryGeoIntx(b2, ptOnLoc3, crs34, &intxList, &intxCount, &intxCode, TOL, EPS);

        if (err) {
          errorCount++;
          failedCount++;
        }
        else {
          if (intxCount == intxCountExp) {
            if (intxCount == 3 && ptsAreSame(intx, intxList[0], TESTTOL) && ptsAreSame(intx1, intxList[1], TESTTOL) && ptsAreSame(intx, intxList[2], TESTTOL)) 
              passedCount++;
            else if (intxCount == 2 && ptsAreSame(intx, intxList[0], TESTTOL) && ptsAreSame(intx1, intxList[1], TESTTOL)) 
              passedCount++;
            else if (intxCount == 1 && ptsAreSame(intx, intxList[0], TESTTOL)) 
              passedCount++;
            else {
              printf("Failed loc intx point jj = %d i = %d j = %d\n",jj,i,j);
              failedCount++;
            }  
          }
          else {
            printf("Failed loc intx count jj = %d i = %d j = %d intxCount = %d intxCountExp %d\n",jj,i,j,intxCount,intxCountExp);
            failedCount++;
          }
        }
        testCaseCount++;
     } //for j
    } //for i

    //Test lines relative to the spiral1 of boundary 2
    crsst = spiral1->startAz;
    for (i = 0; i < 5; i++) //increment crsi from spiral center
    {
      err = 0;
      intxCountExp = 2;
      if (i == 0)
        intxCountExp = 3;
      crsi = crsst + spiral1->subtendedAngle * .1 * ((double)i);
      err |= ptOnSpiral(*spiral1, crsi, &intx, EPS);
      for (j = 1; j < 5; j++) //increment crsj from spiral center
      {
        crsj = crsst + spiral1->subtendedAngle * (1.0  - .1 * ((double)j));
        err |= ptOnSpiral(*spiral1, crsj, &intx1, EPS);
        err |= invCrs(intx, intx1, &crsarc, &tempd1, EPS);
        err |= bndryGeoIntx(b2, intx, crsarc, &intxList, &intxCount, &intxCode, TOL, EPS);

        if (err) {
          errorCount++;
          failedCount++;
        }
        else {
          if (intxCount == intxCountExp) {
            if (intxCount == 3 && ptsAreSame(intx, intxList[0], TESTTOL) && ptsAreSame(intx, intxList[1], TESTTOL) && ptsAreSame(intx1, intxList[2], TESTTOL)) 
              passedCount++;
            else if (intxCount == 2 && ptsAreSame(intx, intxList[0], TESTTOL) && ptsAreSame(intx1, intxList[1], TESTTOL)) 
              passedCount++;
            else if (intxCount == 1 && ptsAreSame(intx, intxList[0], TESTTOL)) 
              passedCount++;
            else {
              failedCount++;
              printf("Failed spiral intx point jj = %d i = %d j = %d\n",jj,i,j);
            }
          }
          else {
            failedCount++;
            printf("Failed spiral intx count jj = %d i = %d j = %d intxCount = %d intxCountExp %d\n",jj,i,j,intxCount,intxCountExp);
          }
        }
        testCaseCount++;
      } //for j
    } //for i
    } //for jj

    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("Finish testBndryGeoIntx_Set4\n\n\n");

    return set;
}

/*
 * NAME: testBndryGeoIntx_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the boundary's bndryGeoIntx function.
 *
 * 		This function runs all the test cases and returns
 * 		the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testBndryGeoIntx_AllSets(TestSuite) - A test suite with the folling metrics:
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
TestSuite testBndryGeoIntx_AllSets()
{
	TestSuite suite;
	TestSet set1;
	TestSet set2;
	TestSet set3;
	TestSet set4;

    printf("\nStart testBndryGeoIntx_AllSets\n");

    suite = newTestSuite("testBndryGeoIntx_AllSets");

    set1 = testBndryGeoIntx_Set1();
    addTestSet(set1,&suite);

    set2 = testBndryGeoIntx_Set2();
    addTestSet(set2,&suite);

    set3 = testBndryGeoIntx_Set3();
    addTestSet(set3,&suite);

    set4 = testBndryGeoIntx_Set4();
    addTestSet(set4,&suite);

    displayTestSuite(suite);

    printf("\nFinish testBndryGeoIntx_AllSets\n\n\n");

    return suite;
}

/*
 * NAME: testBndryLocusIntx_Set1
 *
 * DESCRIPTION:
 * 		This function is used to test the bndryLocusIntx function.
 *
 * 		This function runs the test data created by Richard Snow.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testBndryLocusIntx_Set1(TestSet) - A test set with the folling metrics:
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
TestSet testBndryLocusIntx_Set1()
{
    double DEG2RAD = M_PI / 180.0;
    //const double RAD2DEG = 180.0 / M_PI;
    int passedCount=0, failedCount=0;
    int errorCount = 0;
    int unverifiedCount = 0;
    int setupFailureCount = 0;
    int testCaseCount = 0;
    int outputMatlab = 0; //0 = don't output matlab code to the console
    int i, ii, jj;
    ErrorSet err=0;
    long newSeed = 20080523;
    double latS, lonS, crs1, crs2, crs3, crs4, radius;
    double distgeoStart, distgeoEnd, distlocStart,distlocEnd;
    double crsX1XP1, crsX2XP2, crslocgeo1X1, crslocgeo2X2;
    double h1,h2,fcrs1,fcrs2,bcrs1,bcrs2,d1,d2,slope1,hs,he;
    LineType lineType;
    Geodesic locgeo;
    Locus loc;
    LLPoint locgeoStart,locgeoEnd,X1,X2,XP1,XP2;
    Geodesic* geo1;
    Geodesic* geo2;
    Geodesic* geo3;
    Geodesic* geo4; 
    LLPoint center, geo1Start, geo2Start, geo3Start, geo4Start;
    Boundary b1;
    LLPoint* intxList = NULL;
    int intxCount = 0, intxCountExp = 0;

    TestSet set;
    set = newTestSet("testBndryLocusIntx_Set1");

    printf("\n\nStart testBndryLocusIntx_Set1\n");

    srand(newSeed);  //Initialize the random number generator
 
    lineType = SEGMENT;

    for (jj = 0; jj < 100; jj++)
    {
    b1 = createBndry();

    latS = DEG2RAD * randLat();
    lonS = DEG2RAD * randLon();
    crs1 = DEG2RAD * ((double)((rand() % 89) + 1));
    crs2 = DEG2RAD * ((double)((rand() % 89) + 91));
    crs3 = DEG2RAD * ((double)((rand() % 89) + 181));
    crs4 = DEG2RAD * ((double)((rand() % 89) + 271));
    radius = (double)((rand() % 100) + 1);

    //printf("\nlatS %4.8f lonS %4.8f crs1 %4.8f crs2 %4.8f crs3 %4.8f crs4 %4.8f rad %4.8f\n",latS*RAD2DEG,lonS*RAD2DEG,crs1*RAD2DEG,crs2*RAD2DEG,crs3*RAD2DEG,crs4*RAD2DEG,radius);
    //boundary 1
    center.latitude = latS;
    center.longitude = lonS;
    err |= constructGeoBoundary(center, radius, crs1, crs2, crs3, crs4, &b1);
    geo1 = (Geodesic*) &b1.elements[0];
    geo2 = (Geodesic*) &b1.elements[1];
    geo3 = (Geodesic*) &b1.elements[2];
    geo4 = (Geodesic*) &b1.elements[3];
    geo1Start = geo1->startPoint;
    geo2Start = geo2->startPoint;
    geo3Start = geo3->startPoint;
    geo4Start = geo4->startPoint;
    if (err) printf("\nError(s) occurred during setup: 0x%lx", err);

        if (outputMatlab) displayMatlabGeo(*geo1, "geo1", 0);
        if (outputMatlab) displayMatlabGeo(*geo2, "geo2", 0);
        if (outputMatlab) displayMatlabGeo(*geo3, "geo3", 0);
        if (outputMatlab) displayMatlabGeo(*geo4, "geo4", 0);

    //Test locus intersecting two sides
    for (ii = 0; ii < 4; ii++) //each successive pair of geos
    { 
    for (i = 1; i < 5; i++) //4 increments on each geo
    {
      err = 0;
      intxCountExp = 2;
      if (ii == 0)
      {
        distgeoStart = geo1->length *.2 * ((double)(i-1));
        distgeoEnd = geo2->length * .2 * ((double)(i+1));
        distlocStart = geo1->length *.2 * ((double)i);
        distlocEnd = geo2->length * .2 * ((double)i);
        err |= direct(geo1Start, geo1->startAz, distgeoStart, &locgeoStart, EPS);
        err |= direct(geo2Start, geo2->startAz, distgeoEnd, &locgeoEnd, EPS);
        err |= createGeo(&locgeo, locgeoStart, locgeoEnd, INFINITE, EPS);
        err |= direct(geo1Start, geo1->startAz, distlocStart, &X1, EPS);
        err |= direct(geo2Start, geo2->startAz, distlocEnd, &X2, EPS);
      }
      else if (ii == 1)
      {
        distgeoStart = geo2->length *.2 * ((double)(i-1));
        distgeoEnd = geo3->length * .2 * ((double)(i+1));
        distlocStart = geo2->length *.2 * ((double)i);
        distlocEnd = geo3->length * .2 * ((double)i);
        err |= direct(geo2Start, geo2->startAz, distgeoStart, &locgeoStart, EPS);
        err |= direct(geo3Start, geo3->startAz, distgeoEnd, &locgeoEnd, EPS);
        err |= createGeo(&locgeo, locgeoStart, locgeoEnd, INFINITE, EPS);
        err |= direct(geo2Start, geo2->startAz, distlocStart, &X1, EPS);
        err |= direct(geo3Start, geo3->startAz, distlocEnd, &X2, EPS);
      }
      else if (ii == 2)
      {
        distgeoStart = geo3->length *.2 * ((double)(i-1));
        distgeoEnd = geo4->length * .2 * ((double)(i+1));
        distlocStart = geo3->length *.2 * ((double)i);
        distlocEnd = geo4->length * .2 * ((double)i);
        err |= direct(geo3Start, geo3->startAz, distgeoStart, &locgeoStart, EPS);
        err |= direct(geo4Start, geo4->startAz, distgeoEnd, &locgeoEnd, EPS);
        err |= createGeo(&locgeo, locgeoStart, locgeoEnd, INFINITE, EPS);
        err |= direct(geo3Start, geo3->startAz, distlocStart, &X1, EPS);
        err |= direct(geo4Start, geo4->startAz, distlocEnd, &X2, EPS);
      }
      else if (ii == 3)
      {
        distgeoStart = geo4->length *.2 * ((double)(i-1));
        distgeoEnd = geo1->length * .2 * ((double)(i+1));
        distlocStart = geo4->length *.2 * ((double)i);
        distlocEnd = geo1->length * .2 * ((double)i);
        err |= direct(geo4Start, geo4->startAz, distgeoStart, &locgeoStart, EPS);
        err |= direct(geo1Start, geo1->startAz, distgeoEnd, &locgeoEnd, EPS);
        err |= createGeo(&locgeo, locgeoStart, locgeoEnd, INFINITE, EPS);
        err |= direct(geo4Start, geo4->startAz, distlocStart, &X1, EPS);
        err |= direct(geo1Start, geo1->startAz, distlocEnd, &X2, EPS);
      }

      //Project X1 and X2 on to the locus geodesic
      err |= projectToGeo(locgeoStart, locgeo.startAz, X1, &XP1, &crsX1XP1, &h1, TOL, EPS);
      err |= projectToGeo(locgeoStart, locgeo.startAz, X2, &XP2, &crsX2XP2, &h2, TOL, EPS);
      //Compute the distances (d1, d2) along the locus geodesic of points XP1 and XP2 from geo start
      err |= inverse(locgeoStart, XP1, &fcrs1, &bcrs1, &d1, EPS);
      err |= inverse(locgeoStart, XP2, &fcrs2, &bcrs2, &d2, EPS);
      //Correct the sign of h1, h2, d1 and d2 based on the courses locgeo.startAz, fcrs1, fcrs2...
      err |= inverse(locgeoStart, X1, &crslocgeo1X1, &bcrs1, NULL, EPS);
      if (fabs(locgeo.startAz - crslocgeo1X1) > M_PI)
      {
        if (crslocgeo1X1 > locgeo.startAz)
          h1 = -h1;
      }
      else if (locgeo.startAz > crslocgeo1X1)
        h1 = -h1;

      err |= inverse(locgeoEnd, X2, &crslocgeo2X2, &bcrs1, NULL, EPS);
      if (fabs(crslocgeo2X2 - modcrs(locgeo.endAz + M_PI)) > M_PI)
      {
        if (modcrs(locgeo.endAz + M_PI) > crslocgeo2X2)
          h2 = -h2;
      }
      else if (crslocgeo2X2 > modcrs(locgeo.endAz + M_PI))
        h2 = -h2;

      if (fabs(locgeo.startAz - fcrs1) > 0.1)
        d1 = -d1;
      if (fabs(locgeo.startAz - fcrs2) > 0.1)
        d2 = -d2;

      //Compute slope for locus (h2-h1)/(d2-d1)
      slope1 = (h2 -h1)/(d2-d1);

      //Compute start and end distances for the locus
      hs = h1 -d1*slope1;
      he = h1 + (locgeo.length - d1)*slope1;

      //Create the locus
      err |= createLocus(&loc,locgeoStart,locgeoEnd,hs,he,INFINITE, TOL, EPS);

        err |= bndryLocusIntx(b1, loc, &intxList, &intxCount, TOL, EPS);
        if (err) {
          errorCount++;
          failedCount++;
        }
        else {
          if (intxCount == intxCountExp) {
            if (intxCount == 2) {
              if ((ii < 3) && ptsAreSame(X1, intxList[0], TESTTOL) && ptsAreSame(X2, intxList[1], TESTTOL)) 
                passedCount++;
              else if ((ii == 3) && ptsAreSame(X1, intxList[1], TESTTOL) && ptsAreSame(X2, intxList[0], TESTTOL)) 
                passedCount++;
            }
            else {
              failedCount++;
              printf("Failed intx point jj = %d ii = %d i = %d\n",jj,ii,i);
            }
          }
          else {
            failedCount++;
            printf("Failed intx count jj = %d ii = %d i = %d count %d expcount %d\n",jj,ii,i,intxCount,intxCountExp);
          }
        }
        testCaseCount++;
    } //for i
    } //for ii
    } //for jj

    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("Finish testBndryLocusIntx_Set1\n\n\n");

    return set;
}

/*
 * NAME: testBndryLocusIntx_Set2
 *
 * DESCRIPTION:
 * 		This function is used to test the bndryLocusIntx function.
 *
 * 		This function runs the test data created by Richard Snow.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testBndryLocusIntx_Set2(TestSet) - A test set with the folling metrics:
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
TestSet testBndryLocusIntx_Set2()
{
    double DEG2RAD = M_PI / 180.0;
    //const double RAD2DEG = 180.0 / M_PI;
    int passedCount=0, failedCount=0;
    int errorCount = 0;
    int unverifiedCount = 0;
    int setupFailureCount = 0;
    int testCaseCount = 0;
    int outputMatlab = 0; //0 = don't output matlab code to the console
    int i, jj;
    ErrorSet err=0;
    long newSeed = 20080523;
    double latS, lonS, crs, geolen, perpCrs;
    double tempd1, locDist;
    double crsst, crsi, arcradius;
    double distgeoStart, distgeoEnd, distlocStart,distlocEnd;
    double crsX1XP1, crsX2XP2, crslocgeo1X1, crslocgeo2X2;
    double h1,h2,fcrs1,fcrs2,bcrs1,bcrs2,d1,d2,slope1,hs,he;
    LineType lineType;
    LLPoint geoStart, geoEnd, geoPt, intx, intx1;
    LLPoint locgeoStart,locgeoEnd,X1,X2,XP1,XP2;
    Geodesic locgeo;
    Locus loc;
    Locus* loc3; 
    Locus* loc4;
    Arc* arc1;
    Arc* arc2;
    Boundary b2;
    LLPoint* intxList = NULL;
    int intxCount = 0, intxCountExp = 0;

    TestSet set;
    set = newTestSet("testBndryLocusIntx_Set2");

    printf("\n\nStart testBndryLocusIntx_Set2\n");

    srand(newSeed);  //Initialize the random number generator

    lineType = SEGMENT;
    for (jj = 0; jj < 100; jj++)
    {
    err = 0;
    b2 = createBndry();

    latS = DEG2RAD * randLat();
    lonS = DEG2RAD * randLon();
    crs = DEG2RAD * randAzimuth();
    geolen = (double)((rand() % 100) + 1);
    locDist = (double)((rand() % 20) + 1);

    //printf("\nlatS %4.8f lonS %4.8f crs %4.8f geolen %4.8f locDist %4.8f\n",latS*RAD2DEG,lonS*RAD2DEG,crs*RAD2DEG,geolen,locDist);
    geoStart.latitude = latS;
    geoStart.longitude = lonS;

    err |= direct(geoStart, crs, geolen, &geoEnd, EPS);
    //boundary 2
    err |= constructLocusArcBoundary(geoStart, crs, geolen, locDist, &b2);
    loc3 = (Locus*) &b2.elements[0];
    arc1 = (Arc*) &b2.elements[1];
    loc4  = (Locus*) &b2.elements[2];
    arc2 = (Arc*) &b2.elements[3];

    if (err) printf("\nError(s) occurred during setup: 0x%lx", err);

        if (outputMatlab) displayMatlabLocus(*loc3, "loc3", 0);
        if (outputMatlab) displayMatlabArc(*arc1, "arc1", 0);
        if (outputMatlab) displayMatlabLocus(*loc4, "loc4", 0);
        if (outputMatlab) displayMatlabArc(*arc2, "arc2", 0);

    //Test locus relative to loc3 and loc4 of boundary 2
    for (i = 1; i < 5; i++) //increment along geodesic of loci
    {
      err = 0;
      intxCountExp = 2;
      distgeoStart = geolen * .18 * ((double)(i-1));
      distgeoEnd = geolen * .18 * ((double)i);
      distlocStart = geolen * .18 * ((double)i);
      distlocEnd = geolen * .18 * ((double)(i+1));
      err |= direct(geoStart, crs, distgeoStart, &geoPt, EPS);
      err |= ptOnLocusFromGeoPt(*loc3, geoPt, &locgeoStart, &perpCrs, TOL,EPS);
      err |= direct(geoStart, crs, distgeoEnd, &geoPt, EPS);
      err |= ptOnLocusFromGeoPt(*loc4, geoPt, &locgeoEnd, &perpCrs, TOL,EPS);
      err |= createGeo(&locgeo, locgeoStart, locgeoEnd, INFINITE, EPS);
      err |= direct(geoStart, crs, distlocStart, &geoPt, EPS);
      err |= ptOnLocusFromGeoPt(*loc3, geoPt, &X1, &perpCrs, TOL,EPS);
      err |= direct(geoStart, crs, distlocEnd, &geoPt, EPS);
      err |= ptOnLocusFromGeoPt(*loc4, geoPt, &X2, &perpCrs, TOL,EPS);

      //Project X1 and X2 on to the locus geodesic
      err |= projectToGeo(locgeoStart, locgeo.startAz, X1, &XP1, &crsX1XP1, &h1, TOL, EPS);
      err |= projectToGeo(locgeoStart, locgeo.startAz, X2, &XP2, &crsX2XP2, &h2, TOL, EPS);
      //Compute the distances (d1, d2) along the locus geodesic of points XP1 and XP2 from geo start
      err |= inverse(locgeoStart, XP1, &fcrs1, &bcrs1, &d1, EPS);
      err |= inverse(locgeoStart, XP2, &fcrs2, &bcrs2, &d2, EPS);
      //Correct the sign of h1, h2, d1 and d2 based on the courses locgeo.startAz, fcrs1, fcrs2...
      err |= inverse(locgeoStart, X1, &crslocgeo1X1, &bcrs1, NULL, EPS);
      if (fabs(locgeo.startAz - crslocgeo1X1) > M_PI)
      {
        if (crslocgeo1X1 > locgeo.startAz)
          h1 = -h1;
      }
      else if (locgeo.startAz > crslocgeo1X1)
        h1 = -h1;

      err |= inverse(locgeoEnd, X2, &crslocgeo2X2, &bcrs1, NULL, EPS);
      if (fabs(crslocgeo2X2 - modcrs(locgeo.endAz + M_PI)) > M_PI)
      {
        if (modcrs(locgeo.endAz + M_PI) > crslocgeo2X2)
          h2 = -h2;
      }
      else if (crslocgeo2X2 > modcrs(locgeo.endAz + M_PI))
        h2 = -h2;

      if (fabs(locgeo.startAz - fcrs1) > 0.1)
        d1 = -d1;
      if (fabs(locgeo.startAz - fcrs2) > 0.1)
        d2 = -d2;

      //Compute slope for locus (h2-h1)/(d2-d1)
      slope1 = (h2 -h1)/(d2-d1);

      //Compute start and end distances for the locus
      hs = h1 -d1*slope1;
      he = h1 + (locgeo.length - d1)*slope1;

      //Create the locus
      err |= createLocus(&loc,locgeoStart,locgeoEnd,hs,he,INFINITE, TOL, EPS);

        err |= bndryLocusIntx(b2, loc, &intxList, &intxCount, TOL, EPS);
        if (err) {
          errorCount++;
          failedCount++;
        }
        else {
          if (intxCount == intxCountExp) {
            if (intxCount == 2 && ptsAreSame(X1, intxList[0], TESTTOL) && ptsAreSame(X2, intxList[1], TESTTOL)) 
              passedCount++;
            else {
              printf("Failed loc intx point jj = %d i = %d\n",jj,i);
              failedCount++;
            }  
          }
          else {
            printf("Failed loc intx count jj = %d i = %d intxCount = %d intxCountExp %d\n",jj,i,intxCount,intxCountExp);
            failedCount++;
          }
        }
        testCaseCount++;
    } //for i

    //Test locus relative to the arc1 of boundary 2
    err |= invCrs(geoStart, loc3->locusStart, &crsst, &tempd1, EPS);
    arcradius = arc1->radius;
    for (i = 1; i < 4; i++) //increment crsi from arc center
    {
      err = 0;
      intxCountExp = 2;
      crsi = crsst  + M_PI_2 * .2 * ((double)(i-1));
      err |= direct(geoStart, crsi, arcradius, &locgeoStart, EPS);
      crsi = crsst  + M_PI_2 * .2 * ((double)i);
      err |= direct(geoStart, crsi, arcradius, &X1, EPS);
      crsi = crsst + M_PI  - M_PI_2 * .2 * ((double)(i+1));
      err |= direct(geoStart, crsi, arcradius, &X2, EPS);
      crsi = crsst + M_PI  - M_PI_2 * .2 * ((double)i);
      err |= direct(geoStart, crsi, arcradius, &locgeoEnd, EPS);
      err |= createGeo(&locgeo, locgeoStart, locgeoEnd, INFINITE, EPS);

      //Project X1 and X2 on to the locus geodesic
      err |= projectToGeo(locgeoStart, locgeo.startAz, X1, &XP1, &crsX1XP1, &h1, TOL, EPS);
      err |= projectToGeo(locgeoStart, locgeo.startAz, X2, &XP2, &crsX2XP2, &h2, TOL, EPS);
      //Compute the distances (d1, d2) along the locus geodesic of points XP1 and XP2 from geo start
      err |= inverse(locgeoStart, XP1, &fcrs1, &bcrs1, &d1, EPS);
      err |= inverse(locgeoStart, XP2, &fcrs2, &bcrs2, &d2, EPS);
      //Correct the sign of h1, h2, d1 and d2 based on the courses locgeo.startAz, fcrs1, fcrs2...
      err |= inverse(locgeoStart, X1, &crslocgeo1X1, &bcrs1, NULL, EPS);
      if (fabs(locgeo.startAz - crslocgeo1X1) > M_PI)
      {
        if (crslocgeo1X1 > locgeo.startAz)
          h1 = -h1;
      }
      else if (locgeo.startAz > crslocgeo1X1)
        h1 = -h1;

      err |= inverse(locgeoEnd, X2, &crslocgeo2X2, &bcrs1, NULL, EPS);
      if (fabs(crslocgeo2X2 - modcrs(locgeo.endAz + M_PI)) > M_PI)
      {
        if (modcrs(locgeo.endAz + M_PI) > crslocgeo2X2)
          h2 = -h2;
      }
      else if (crslocgeo2X2 > modcrs(locgeo.endAz + M_PI))
        h2 = -h2;

      if (fabs(locgeo.startAz - fcrs1) > 0.1)
        d1 = -d1;
      if (fabs(locgeo.startAz - fcrs2) > 0.1)
        d2 = -d2;

      //Compute slope for locus (h2-h1)/(d2-d1)
      slope1 = (h2 -h1)/(d2-d1);

      //Compute start and end distances for the locus
      hs = h1 -d1*slope1;
      he = h1 + (locgeo.length - d1)*slope1;

      //Create the locus
      err |= createLocus(&loc,locgeoStart,locgeoEnd,hs,he,INFINITE, TOL, EPS);

        err |= bndryLocusIntx(b2, loc, &intxList, &intxCount, TOL, EPS);

        if (err) {
          errorCount++;
          failedCount++;
        }
        else {
          if (intxCount == intxCountExp) {
            if (intxCount == 3 && ptsAreSame(intx, intxList[0], TESTTOL) && ptsAreSame(intx, intxList[1], TESTTOL) && ptsAreSame(intx1, intxList[2], TESTTOL)) 
              passedCount++;
            else if (intxCount == 2 && ((ptsAreSame(X1, intxList[0], TESTTOL) && ptsAreSame(X2, intxList[1], TESTTOL)) || (ptsAreSame(X2, intxList[0], TESTTOL) && ptsAreSame(X1, intxList[1], TESTTOL))) ) 
              passedCount++;
            else if (intxCount == 1 && ptsAreSame(intx, intxList[0], TESTTOL)) 
              passedCount++;
            else {
              failedCount++;
              printf("Failed arc intx point jj = %d i = %d\n",jj,i);
            }
          }
          else {
            failedCount++;
            printf("Failed arc intx count jj = %d i = %d intxCount = %d intxCountExp %d\n",jj,i,intxCount,intxCountExp);
          }
        }
        testCaseCount++;
    } //for i
    } //for jj

    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("Finish testBndryLocusIntx_Set2\n\n\n");

    return set;
}

/*
 * NAME: testBndryLocusIntx_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the boundary's bndryLocusIntx function.
 *
 * 		This function runs all the test cases and returns
 * 		the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testBndryLocusIntx_AllSets(TestSuite) - A test suite with the folling metrics:
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
TestSuite testBndryLocusIntx_AllSets()
{
	TestSuite suite;
	TestSet set1;
	TestSet set2;

    printf("\nStart testBndryLocusIntx_AllSets\n");

    suite = newTestSuite("testBndryLocusIntx_AllSets");

    set1 = testBndryLocusIntx_Set1();
    addTestSet(set1,&suite);

    set2 = testBndryLocusIntx_Set2();
    addTestSet(set2,&suite);

    displayTestSuite(suite);

    printf("\nFinish testBndryLocusIntx_AllSets\n\n\n");

    return suite;
}

/*
 * NAME: testBndryArcIntx_Set2
 *
 * DESCRIPTION:
 * 		This function is used to test the bndryArcIntx function.
 *
 * 		This function runs the test data created by Rich Snow.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testBndryArcIntx_Set2(TestSet) - A test set with the folling metrics:
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
TestSet testBndryArcIntx_Set2()
{
    double DEG2RAD = M_PI / 180.0;
    //const double RAD2DEG = 180.0 / M_PI;
    int passedCount=0, failedCount=0;
    int errorCount = 0;
    int unverifiedCount = 0;
    int setupFailureCount = 0;
    int testCaseCount = 0;
    int outputMatlab = 0; //0 = don't output matlab code to the console
    int i, ii, jj;
    ErrorSet err=0;
    long newSeed = 20080523;
    double latS, lonS, crs1, crs2, crs3, crs4, radius, disti;
    double delta, startcrs, geolength;
    LineType lineType;
    Geodesic* geo1;
    Geodesic* geo2;
    Geodesic* geo3;
    Geodesic* geo4; 
    LLPoint center, geo1Start, geo2Start, geo3Start, geo4Start, intx, intx1;
    LLPoint arcStart, arcCenter, arcEnd;
    Arc testArc;
    Boundary b1;
    LLPoint* intxList = NULL;
    int intxCount = 0, intxCountExp = 0;

    TestSet set;
    set = newTestSet("testBndryArcIntx_Set2");

    printf("\n\nStart testBndryArcIntx_Set2\n");

    srand(newSeed);  //Initialize the random number generator
 
    lineType = SEGMENT;

    for (jj = 0; jj < 100; jj++)
    {
    err = 0;
    b1 = createBndry();

    latS = randLat();
    if (latS > 89.0)
      latS = latS - 1.0;
    latS = DEG2RAD * latS;
    lonS = DEG2RAD * randLon();
    crs1 = DEG2RAD * ((double)((rand() % 89) + 1));
    crs2 = DEG2RAD * ((double)((rand() % 89) + 91));
    crs3 = DEG2RAD * ((double)((rand() % 89) + 181));
    crs4 = DEG2RAD * ((double)((rand() % 89) + 271));
    radius = (double)((rand() % 100) + 1);

    //printf("\nlatS %4.8f lonS %4.8f crs1 %4.8f crs2 %4.8f crs3 %4.8f crs4 %4.8f rad %4.8f\n",latS*RAD2DEG,lonS*RAD2DEG,crs1*RAD2DEG,crs2*RAD2DEG,crs3*RAD2DEG,crs4*RAD2DEG,radius);
    //boundary 1
    center.latitude = latS;
    center.longitude = lonS;
    err |= constructGeoBoundary(center, radius, crs1, crs2, crs3, crs4, &b1);
    geo1 = (Geodesic*) &b1.elements[0];
    geo2 = (Geodesic*) &b1.elements[1];
    geo3 = (Geodesic*) &b1.elements[2];
    geo4 = (Geodesic*) &b1.elements[3];
    geo1Start = geo1->startPoint;
    geo2Start = geo2->startPoint;
    geo3Start = geo3->startPoint;
    geo4Start = geo4->startPoint;
    if (err) printf("\nError(s) occurred during setup: 0x%lx", err);

        if (outputMatlab) displayMatlabGeo(*geo1, "geo1", 0);
        if (outputMatlab) displayMatlabGeo(*geo2, "geo2", 0);
        if (outputMatlab) displayMatlabGeo(*geo3, "geo3", 0);
        if (outputMatlab) displayMatlabGeo(*geo4, "geo4", 0);

    //Arc has center at vertex
    for (ii = 0; ii < 4; ii++) //each geodesic
    {
      err = 0;
      intxCountExp = 2;
      if (ii == 0)
      {
        if (geo1->length > geo4->length)
          disti = geo4->length/4.0;
        else
          disti = geo1->length/4.0;
        err |= direct(geo1Start, geo1->startAz, disti, &arcStart, EPS);
        err |= direct(geo1Start, geo4->endAz + M_PI, disti, &arcEnd, EPS);
        arcCenter.latitude = geo1Start.latitude;
        arcCenter.longitude = geo1Start.longitude;
      }
      else if (ii == 1)
      {
        if (geo1->length > geo2->length)
          disti = geo2->length/4.0;
        else
          disti = geo1->length/4.0;
        err |= direct(geo2Start, geo2->startAz, disti, &arcStart, EPS);
        err |= direct(geo2Start, geo1->endAz + M_PI, disti, &arcEnd, EPS);
        arcCenter.latitude = geo2Start.latitude;
        arcCenter.longitude = geo2Start.longitude;
      }
      else if (ii == 2)
      {
        if (geo3->length > geo2->length)
          disti = geo2->length/4.0;
        else
          disti = geo3->length/4.0;
        err |= direct(geo3Start, geo3->startAz, disti, &arcStart, EPS);
        err |= direct(geo3Start, geo2->endAz + M_PI, disti, &arcEnd, EPS);
        arcCenter.latitude = geo3Start.latitude;
        arcCenter.longitude = geo3Start.longitude;
      }
      else if (ii == 3)
      {
        if (geo3->length > geo4->length)
          disti = geo4->length/4.0;
        else
          disti = geo3->length/4.0;
        err |= direct(geo4Start, geo4->startAz, disti, &arcStart, EPS);
        err |= direct(geo4Start, geo3->endAz + M_PI, disti, &arcEnd, EPS);
        arcCenter.latitude = geo4Start.latitude;
        arcCenter.longitude = geo4Start.longitude;
      }
       
        intx.latitude = arcStart.latitude;
        intx.longitude = arcStart.longitude;
        intx1.latitude = arcEnd.latitude;
        intx1.longitude = arcEnd.longitude;
        err |= createArc(&testArc, arcCenter, arcStart, arcEnd, CLOCKWISE, TOL, EPS);
        err |= bndryArcIntx(b1, testArc, &intxList, &intxCount,  TOL, EPS);
        if (err) {
          errorCount++;
          failedCount++;
        }
        else {
          if (intxCount == intxCountExp) {
            if (intxCount == 2 && ((ptsAreSame(intx, intxList[0], TESTTOL) && ptsAreSame(intx1, intxList[1], TESTTOL)) || (ptsAreSame(intx1, intxList[0], TESTTOL) && ptsAreSame(intx, intxList[1], TESTTOL))) ) 
              passedCount++;
            else if ((intxCount == 1) && ptsAreSame(intx, intxList[0], TESTTOL)) {
              passedCount++;
            }
            else {
              failedCount++;
              printf("Failed vertex intx point jj = %d ii = %d\n",jj,ii);
            }
          }
          else {
            failedCount++;
            printf("Failed vertex intx count jj = %d ii = %d intxCount = %d intxCountExp %d\n",jj,ii,intxCount,intxCountExp);
          }
        }
        testCaseCount++;
    } //for ii

    //Arc intersects geodesic between end points
    intxCountExp = 2;
    for (ii = 0; ii < 4; ii++)
    {
      if (ii == 0)
        geolength = geo1->length;
      else if (ii == 1)
        geolength = geo2->length;
      else if (ii == 2)
        geolength = geo3->length;
      else if (ii == 3)
        geolength = geo4->length;

      if (.1 * geolength < radius/4.0)
        delta = .1 * geolength;
      else
        delta = radius/4.0;

      for (i = 1; i < 5; i++)
      {
        if (ii == 0)
        {
          disti = geo1->length * .2 * ((double)i) - delta;
          err |= direct(geo1Start, geo1->startAz, disti, &arcStart, EPS);
          startcrs = modcrs(geo1->startAz + M_PI_4);
          disti = geo1->length * .2 * ((double)i) + delta;
          err |= direct(geo1Start, geo1->startAz, disti, &arcEnd, EPS);
        }
        else if (ii == 1) 
        {
          disti = geo2->length * .2 * ((double)i) - delta;
          err |= direct(geo2Start, geo2->startAz, disti, &arcStart, EPS);
          startcrs = modcrs(geo2->startAz + M_PI_4);
          disti = geo2->length * .2 * ((double)i) + delta;
          err |= direct(geo2Start, geo2->startAz, disti, &arcEnd, EPS);
        }
        else if (ii == 2)
        {
          disti = geo3->length * .2 * ((double)i) - delta;
          err |= direct(geo3Start, geo3->startAz, disti, &arcStart, EPS);
          startcrs = modcrs(geo3->startAz + M_PI_4);
          disti = geo3->length * .2 * ((double)i) + delta;
          err |= direct(geo3Start, geo3->startAz, disti, &arcEnd, EPS);
        }
        else if (ii == 3)
        {
          disti = geo4->length * .2 * ((double)i) - delta;
          err |= direct(geo4Start, geo4->startAz, disti, &arcStart, EPS);
          startcrs = modcrs(geo4->startAz + M_PI_4);
          disti = geo4->length * .2 * ((double)i) + delta;
          err |= direct(geo4Start, geo4->startAz, disti, &arcEnd, EPS);
        }

      err |= arcFromStartAndEnd(arcStart, startcrs, arcEnd, &testArc, TOL, EPS);
        intx.latitude = arcStart.latitude;
        intx.longitude = arcStart.longitude;
        intx1.latitude = arcEnd.latitude;
        intx1.longitude = arcEnd.longitude;

        err |= bndryArcIntx(b1, testArc, &intxList, &intxCount, TOL, EPS);
        if (err) {
          errorCount++;
          failedCount++;
        }
        else {
          if (intxCount == intxCountExp) {
            if (intxCount == 2 && ((ptsAreSame(intx, intxList[0], TESTTOL) && ptsAreSame(intx1, intxList[1], TESTTOL)) || (ptsAreSame(intx1, intxList[0], TESTTOL) && ptsAreSame(intx, intxList[1], TESTTOL))) ) {
              passedCount++;
            }
            else {
              failedCount++;
              printf("Failed geo intx point jj = %d ii = %d i = %d\n",jj,ii,i);
            }
          }
          else {
            failedCount++;
            printf("Failed geo intx count jj = %d ii = %d i = %d intxCount = %d intxCountExp %d\n",jj,ii,i,intxCount,intxCountExp);
          }
        }
        testCaseCount++;
    } //for i
    } //for ii

    } //for jj

    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("Finish testBndryArcIntx_Set2\n\n\n");

    return set;
}

/*
 * NAME: testBndryArcIntx_Set3
 *
 * DESCRIPTION:
 * 		This function is used to test the bndryArcIntx function.
 *
 * 		This function runs the test data created by Rich Snow.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testBndryArcIntx_Set3(TestSet) - A test set with the folling metrics:
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
TestSet testBndryArcIntx_Set3()
{
    double DEG2RAD = M_PI / 180.0;
    //const double RAD2DEG = 180.0 / M_PI;
    int passedCount=0, failedCount=0;
    int errorCount = 0;
    int unverifiedCount = 0;
    int setupFailureCount = 0;
    int testCaseCount = 0;
    int outputMatlab = 0; //0 = don't output matlab code to the console
    int i, jj, n1, n2;
    ErrorSet err=0;
    long newSeed = 20080523;
    double latS, lonS, crs, geolen, disti, perpCrs, locDist;
    double startcrs, delta;
    LineType lineType;
    LLPoint geoStart, geoEnd, geoPt, intx, intx1, arcStart, arcEnd, tempLLPoint;
    Locus* loc3; 
    Locus* loc4;
    Arc* arc1;
    Arc* arc2;
    Arc testArc;
    Boundary b2;
    LLPoint* intxList = NULL;
    LLPointPair tmpIntx1, tmpIntx2;
    int intxCount = 0, intxCountExp = 0;

    TestSet set;
    set = newTestSet("testBndryArcIntx_Set3");

    printf("\n\nStart testBndryArcIntx_Set3\n");

    srand(newSeed);  //Initialize the random number generator

    lineType = SEGMENT;
    for (jj = 0; jj < 100; jj++)
    {
    b2 = createBndry();

    latS = randLat();
    if (latS >= 89.3)
      latS = 88.0;
    latS = DEG2RAD * latS;
    lonS = DEG2RAD * randLon();
    crs = DEG2RAD * randAzimuth();
    geolen = (double)((rand() % 100) + 1);
    locDist = (double)((rand() % 20) + 1);

    //printf("\nlatS %4.8f lonS %4.8f crs %4.8f geolen %4.8f locDist %4.8f\n",latS*RAD2DEG,lonS*RAD2DEG,crs*RAD2DEG,geolen,locDist);
    geoStart.latitude = latS;
    geoStart.longitude = lonS;

    err |= direct(geoStart, crs, geolen, &geoEnd, EPS);
    //boundary 2
    err |= constructLocusArcBoundary(geoStart, crs, geolen, locDist, &b2);
    loc3 = (Locus*) &b2.elements[0];
    arc1 = (Arc*) &b2.elements[1];
    loc4  = (Locus*) &b2.elements[2];
    arc2 = (Arc*) &b2.elements[3];

    if (err) printf("\nError(s) occurred during setup: 0x%lx", err);

        if (outputMatlab) displayMatlabLocus(*loc3, "loc3", 0);
        if (outputMatlab) displayMatlabArc(*arc1, "arc1", 0);
        if (outputMatlab) displayMatlabLocus(*loc4, "loc4", 0);
        if (outputMatlab) displayMatlabArc(*arc2, "arc2", 0);

    //Test arcs relative to loc3 of boundary 2
    for (i = 1; i < 5; i++) //increment disti along geodesic of loci
    {
      err = 0;
      intxCountExp = 2;
      if (.1 * geolen < locDist/4.0)
        delta = .1 * geolen;
      else
        delta = locDist/4.0;
      disti = geolen * .2 * ((double)i) - delta;
      err |= direct(geoStart, crs, disti, &geoPt, EPS);
      err |= ptOnLocusFromGeoPt(*loc3, geoPt, &arcStart, &perpCrs, TOL, EPS);
      //startcrs = modcrs(perpCrs + M_PI - M_PI_4);
      startcrs = modcrs(perpCrs + M_PI);
      disti = geolen * .2 * ((double)i) + delta;
      err |= direct(geoStart, crs, disti, &geoPt, EPS);
      err |= ptOnLocusFromGeoPt(*loc3, geoPt, &arcEnd, &perpCrs, TOL, EPS);
      err |= arcFromStartAndEnd(arcStart, startcrs, arcEnd, &testArc, TOL, EPS);
        intx.latitude = arcStart.latitude;
        intx.longitude = arcStart.longitude;
        intx1.latitude = arcEnd.latitude;
        intx1.longitude = arcEnd.longitude;
        err |= bndryArcIntx(b2, testArc, &intxList, &intxCount, TOL, EPS);

        if (err) {
          errorCount++;
          failedCount++;
        }
        else {
          if (intxCount == intxCountExp) {
            if (intxCount == 3 && ptsAreSame(intx, intxList[0], TESTTOL) && ptsAreSame(intx, intxList[1], TESTTOL) && ptsAreSame(intx1, intxList[2], TESTTOL)) 
              passedCount++;
            else if (intxCount == 2 && ((ptsAreSame(intx, intxList[0], TESTTOL) && ptsAreSame(intx1, intxList[1], TESTTOL)) || (ptsAreSame(intx1, intxList[0], TESTTOL) && ptsAreSame(intx, intxList[1], TESTTOL))) ) 
              passedCount++;
            else if (intxCount == 1 && ptsAreSame(intx, intxList[0], TESTTOL)) 
              passedCount++;
            else {
              printf("Failed loc intx point jj = %d i = %d\n",jj,i);
              failedCount++;
            }  
          }
          else {
            printf("Failed loc intx count jj = %d i = %d intxCount = %d intxCountExp %d\n",jj,i,intxCount,intxCountExp);
            failedCount++;
          }
        }
        testCaseCount++;
    } //for i

    //Test arc centered at start point of arc1
    intxCountExp = 2;
    if (geolen < locDist)
      delta = geolen/4.0;
    else
      delta = locDist/4.0;
     
    if (delta > arc1->radius/4.0)
      delta = arc1->radius/4.0;

    err |= direct(arc1->startPoint, 0.0, delta, &tempLLPoint, EPS);
    err |= createArc(&testArc, arc1->startPoint, tempLLPoint, tempLLPoint, CLOCKWISE, TOL, EPS);
    err |= locusArcIntx(*loc3, testArc.centerPoint, delta, tmpIntx1, &n1, TOL, EPS);
      intx.latitude = tmpIntx1[0].latitude;
      intx.longitude = tmpIntx1[0].longitude;
    err |= arcIntx(arc1->centerPoint, arc1->radius, testArc.centerPoint, testArc.radius, tmpIntx2, &n2, TOL, EPS);
    if (n2 == 2)
    {
      if (ptIsOnArc(arc1->centerPoint, arc1->radius, arc1->startAz, arc1->endAz, arc1->dir, tmpIntx2[0], &err, TOL, EPS))
      {
        intx1.latitude = tmpIntx2[0].latitude;
        intx1.longitude = tmpIntx2[0].longitude;
      }
      else
      {
        intx1.latitude = tmpIntx2[1].latitude;
        intx1.longitude = tmpIntx2[1].longitude;
      }
    }
    else if (n2 == 1)
    {
      intx1.latitude = tmpIntx2[0].latitude;
      intx1.longitude = tmpIntx2[0].longitude;
    }

    err |= bndryArcIntx(b2, testArc, &intxList, &intxCount, TOL, EPS);

    if (err) {
      errorCount++;
      failedCount++;
    }
    else {
      if (intxCount == intxCountExp) {
        if (intxCount == 2 && ((ptsAreSame(intx, intxList[0], TESTTOL) && ptsAreSame(intx1, intxList[1], TESTTOL)) || (ptsAreSame(intx1, intxList[0], TESTTOL) && ptsAreSame(intx, intxList[1], TESTTOL))) ) 
          passedCount++;
        else if (intxCount == 1 && ptsAreSame(intx, intxList[0], TESTTOL)) 
          passedCount++;
        else {
          printf("Failed arc intx point jj = %d i = %d\n",jj,i);
          failedCount++;
        }  
      }
      else {
        printf("Failed arc intx count jj = %d i = %d intxCount = %d intxCountExp %d\n",jj,i,intxCount,intxCountExp);
        failedCount++;
      }
    }
    testCaseCount++;

    //The following code is commented out pending further review RES 9/30/11
    /*if (jj > 9)
      continue; 
    //Tests for CONCENTRIC_CIRCLE_ERR
    //These tests only cover the case when the intersection of the 
    //concentric arcs is one common arc
    intxCountExp = 2;
    for (i = 0; i < 4; i++)
    {
      if (i == 0)
      {
        //Clockwise arcs with arc1 containing testArc
        err |= direct(arc1->centerPoint, modcrs(arc1->startAz + M_PI/9.0), arc1->radius, &arcStart, EPS);
        err |= direct(arc1->centerPoint, modcrs(arc1->endAz - M_PI/9.0), arc1->radius, &arcEnd, EPS);
        err |= createArc(&testArc, arc1->centerPoint, arcStart, arcEnd, CLOCKWISE, TOL, EPS);
        intx = arcStart;
        intx1 = arcEnd;
      }
      else if (i == 1)
      {
        //Clockwise arcs with testArc containing arc1
        err |= direct(arc1->centerPoint, modcrs(arc1->startAz - M_PI/9.0), arc1->radius, &arcStart, EPS);
        err |= direct(arc1->centerPoint, modcrs(arc1->endAz + M_PI/9.0), arc1->radius, &arcEnd, EPS);
        err |= createArc(&testArc, arc1->centerPoint, arcStart, arcEnd, CLOCKWISE, TOL, EPS);
        intx = arc1->startPoint;
        intx1 = arc1->endPoint;
      }
      else if (i == 2)
      {
        //Counterclockwise arcs with arc2 containing testArc
        err |= direct(arc2->centerPoint, modcrs(arc2->startAz - M_PI/9.0), arc2->radius, &arcStart, EPS);
        err |= direct(arc2->centerPoint, modcrs(arc2->endAz + M_PI/9.0), arc2->radius, &arcEnd, EPS);
        err |= createArc(&testArc, arc2->centerPoint, arcStart, arcEnd, COUNTERCLOCKWISE, TOL, EPS);
        intx = arcStart;
        intx1 = arcEnd;
      }
      else if (i == 3)
      {
        //Counterclockwise arcs with testArc containing arc2
        err |= direct(arc2->centerPoint, modcrs(arc2->startAz + M_PI/9.0), arc2->radius, &arcStart, EPS);
        err |= direct(arc2->centerPoint, modcrs(arc2->endAz - M_PI/9.0), arc2->radius, &arcEnd, EPS);
        err |= createArc(&testArc, arc2->centerPoint, arcStart, arcEnd, COUNTERCLOCKWISE, TOL, EPS);
        intx = arc2->startPoint;
        intx1 = arc2->endPoint;
      }

      err |= bndryArcIntx(b2, testArc, &intxList, &intxCount, TOL, EPS);

      if (((i == 1) || (i == 3)) && (intxCount == 4))
      {
        if (ptsAreSame(intxList[0], intxList[1], TOL) && ptsAreSame(intxList[2], intxList[3], TOL))
        {
          intxList[1] = intxList[2];
          intxCount = 2;
        }
        else if (ptsAreSame(intxList[0], intxList[2], TOL) && ptsAreSame(intxList[1], intxList[3], TOL))
          intxCount = 2;
        else if (ptsAreSame(intxList[0], intxList[3], TOL) && ptsAreSame(intxList[1], intxList[2], TOL))
          intxCount = 2;
      }

      if ((err & CONCENTRIC_CIRCLE_ERR) && (intxCount == intxCountExp)) {
        if (intxCount == 2 && ((ptsAreSame(intx, intxList[0], TESTTOL) && ptsAreSame(intx1, intxList[1], TESTTOL)) || (ptsAreSame(intx1, intxList[0], TESTTOL) && ptsAreSame(intx, intxList[1], TESTTOL))) ) 
          passedCount++;
        else if (intxCount == 1 && ptsAreSame(intx, intxList[0], TESTTOL)) 
          passedCount++;
        else {
          printf("Failed arc intx point jj = %d i = %d\n",jj,i);
          failedCount++;
        }  
      }
      else {
        printf("Failed arc intx count jj = %d i = %d intxCount = %d intxCountExp %d\n",jj,i,intxCount,intxCountExp);
        failedCount++;
      }
    testCaseCount++;
    err = 0;
    } //for i
    */
    //end of commented out code

    } //for jj

    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("Finish testBndryArcIntx_Set3\n\n\n");

    return set;
}

/*
 * NAME: testBndryArcIntx_Set4
 *
 * DESCRIPTION:
 * 		This function is used to test the bndryArcIntx function.
 *
 * 		This function runs the test data created by Rich Snow.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testBndryArcIntx_Set4(TestSet) - A test set with the folling metrics:
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
TestSet testBndryArcIntx_Set4()
{
    double DEG2RAD = M_PI / 180.0;
    //const double RAD2DEG = 180.0 / M_PI;
    int passedCount=0, failedCount=0;
    int errorCount = 0;
    int unverifiedCount = 0;
    int setupFailureCount = 0;
    int testCaseCount = 0;
    int outputMatlab = 0; //0 = don't output matlab code to the console
    int i, j, jj;
    ErrorSet err=0;
    long newSeed = 20080523;
    double latS, lonS, crs, geolen, disti, perpCrs, locDist, len, startRad, endRad, geoAz;
    double startcrs, delta, crsi, crsj;
    LineType lineType;
    LLPoint geoStart, geoPt, intx, intx1, arcStart, arcEnd, sp;
    Locus* loc3; 
    Locus* loc4;
    Spiral* spiral1;
    Spiral* spiral2;
    Arc testArc;
    Boundary b2;
    LLPoint* intxList = NULL;
    int intxCount = 0, intxCountExp = 0;

    TestSet set;
    set = newTestSet("testBndryArcIntx_Set4");

    printf("\n\nStart testBndryArcIntx_Set4\n");

    srand(newSeed);  //Initialize the random number generator

    lineType = SEGMENT;
    for (jj = 0; jj < 100; jj++)
    {
    b2 = createBndry();

    latS = DEG2RAD * randLat();
    lonS = DEG2RAD * randLon();
    crs = DEG2RAD * randAzimuth();
    len = 120.0 + 0.01 * randDist();
    startRad = 30.0 + 0.01 * randDist() / 2.0;
    endRad = 30.0 + 0.01 * randDist() / 2.0;

    //printf("\nlatS %4.8f lonS %4.8f crs %4.8f len %4.8f startRad %4.8f endRad %4.8f\n",latS*RAD2DEG,lonS*RAD2DEG,crs*RAD2DEG,len,startRad,endRad );
    sp.latitude = latS;
    sp.longitude = lonS;

    //boundary 2
    err |= createLocusSpiralBndry(sp, crs, len, startRad, endRad, &b2);
    loc3 = (Locus*) &b2.elements[0];
    loc4  = (Locus*) &b2.elements[1];
    spiral1 = (Spiral*) &b2.elements[2];
    spiral2 = (Spiral*) &b2.elements[3];

    geoStart = loc3->geoStart;
    geolen = loc3->geoLength;
    geoAz = loc3->geoAz;
    if (loc3->startDist < loc3->endDist)
      locDist = loc3->startDist;
    else
      locDist = loc3->endDist;

    if (err) printf("\nError(s) occurred during setup: 0x%lx", err);

        if (outputMatlab) displayMatlabLocus(*loc3, "loc3", 1);
        if (outputMatlab) displayMatlabLocus(*loc4, "loc4", 1);
        if (outputMatlab) displayMatlabSpiral(*spiral1, "spiral1", 1);
        if (outputMatlab) displayMatlabSpiral(*spiral2, "spiral2", 1);
        /*if (outputMatlab) displayMatlabBndry(b2, "b2", 1);
         displayLocus(*loc3, "loc3", 1);
         displayLocus(*loc4, "loc4", 1);
         displaySpiral(*spiral1, "spiral1", 1);
         displaySpiral(*spiral2, "spiral2", 1);
         displayBndry(b2, "b2", 1);*/

    //Test arcs relative to loc3 of boundary 2
    for (i = 1; i < 5; i++) //increment disti along geodesic of loci
    {
      err = 0;
      intxCountExp = 2;
      if (.1 * geolen < locDist/4.0)
        delta = .1 * geolen;
      else
        delta = locDist/4.0;
      disti = geolen * .2 * ((double)i) - delta;
      err |= direct(geoStart, geoAz, disti, &geoPt, EPS);
      err |= ptOnLocusFromGeoPt(*loc3, geoPt, &arcStart, &perpCrs, TOL, EPS);
      startcrs = modcrs(perpCrs + M_PI);
      disti = geolen * .2 * ((double)i) + delta;
      err |= direct(geoStart, geoAz, disti, &geoPt, EPS);
      err |= ptOnLocusFromGeoPt(*loc3, geoPt, &arcEnd, &perpCrs, TOL, EPS);
      err |= arcFromStartAndEnd(arcStart, startcrs, arcEnd, &testArc, TOL, EPS);
        intx = arcStart;
        intx1 = arcEnd;
        err |= bndryArcIntx(b2, testArc, &intxList, &intxCount, TOL, EPS);

        if (err) {
          errorCount++;
          failedCount++;
        }
        else {
          if (intxCount == intxCountExp) {
            if (intxCount == 3 && ptsAreSame(intx, intxList[0], TESTTOL) && ptsAreSame(intx, intxList[1], TESTTOL) && ptsAreSame(intx1, intxList[2], TESTTOL)) 
              passedCount++;
            else if (intxCount == 2 && ((ptsAreSame(intx, intxList[0], TESTTOL) && ptsAreSame(intx1, intxList[1], TESTTOL)) || (ptsAreSame(intx1, intxList[0], TESTTOL) && ptsAreSame(intx, intxList[1], TESTTOL))) ) 
              passedCount++;
            else if (intxCount == 1 && ptsAreSame(intx, intxList[0], TESTTOL)) 
              passedCount++;
            else {
              printf("Failed loc intx point jj = %d i = %d\n",jj,i);
              failedCount++;
            }  
          }
          else {
            printf("Failed loc intx count jj = %d i = %d intxCount = %d intxCountExp %d\n",jj,i,intxCount,intxCountExp);
            failedCount++;
          }
        }
        testCaseCount++;
    } //for i

    //Test arc relative to the spiral1 of boundary 2
    for (i = 1; i < 4; i++) //increment crsi from spiral center
    {
      err = 0;
      intxCountExp = 2;
      crsi = spiral1->startAz + 40.0 *  ((double)i) * DEG2RAD;
      err |= ptOnSpiral(*spiral1, crsi, &arcStart, EPS);
      err |= spiralTanCrs(*spiral1, crsi, &startcrs);
      startcrs = modcrs(startcrs - M_PI_2);
      for (j = 1; j < 4; j++) //increment crsj from spiral center
      {
        crsj =  crsi + 20.0 * ((double)j) * DEG2RAD;
        if (i == 3 && j == 3) 
          crsj =  crsi + 15.0 * ((double)j) * DEG2RAD;
        err |= ptOnSpiral(*spiral1, crsj, &arcEnd, EPS);
        err |= arcFromStartAndEnd(arcStart, startcrs, arcEnd, &testArc, TOL, EPS);
        //printf("arcrad %f spiral startAz %f endAz %f crsi %f crsj %f\n",testArc.radius,spiral1->startAz,spiral1->endAz,crsi,crsj);
        intx = arcStart;
        intx1 = arcEnd;
        err |= bndryArcIntx(b2, testArc, &intxList, &intxCount, TOL, EPS);

        if (err) {
          errorCount++;
          failedCount++;
        }
        else {
          if (intxCount == intxCountExp) {
            if (intxCount == 3 && ptsAreSame(intx, intxList[0], TESTTOL) && ptsAreSame(intx, intxList[1], TESTTOL) && ptsAreSame(intx1, intxList[2], TESTTOL)) 
              passedCount++;
            else if (intxCount == 2 && ((ptsAreSame(intx, intxList[0], TESTTOL) && ptsAreSame(intx1, intxList[1], TESTTOL)) || (ptsAreSame(intx1, intxList[0], TESTTOL) && ptsAreSame(intx, intxList[1], TESTTOL))) ) 
              passedCount++;
            else if (intxCount == 1 && ptsAreSame(intx, intxList[0], TESTTOL)) 
              passedCount++;
            else {
              failedCount++;
              printf("Failed spiral intx point jj = %d i = %d j = %d\n",jj,i,j);
            }
          }
          else {
            failedCount++;
            printf("Failed spiral intx count jj = %d i = %d j = %d intxCount = %d intxCountExp %d\n",jj,i,j,intxCount,intxCountExp);
          }
        }
        testCaseCount++;
      } //for j
    } //for i
    } //for jj

    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("Finish testBndryArcIntx_Set4\n\n\n");

    return set;
}

/*
 * NAME: testBndryArcIntx_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the boundary's bndryArcIntx function.
 *
 * 		This function runs all the test cases and returns
 * 		the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testBndryArcIntx_AllSets(TestSuite) - A test suite with the folling metrics:
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
TestSuite testBndryArcIntx_AllSets()
{
	TestSuite suite;
	TestSet set2;
	TestSet set3;
	TestSet set4;

    printf("\nStart testBndryArcIntx_AllSets\n");

    suite = newTestSuite("testBndryArcIntx_AllSets");

    set2 = testBndryArcIntx_Set2();
    addTestSet(set2,&suite);

    set3 = testBndryArcIntx_Set3();
    addTestSet(set3,&suite);

    set4 = testBndryArcIntx_Set4();
    addTestSet(set4,&suite);

    displayTestSuite(suite);

    printf("\nFinish testBndryArcIntx_AllSets\n\n\n");

    return suite;
}

/*
 * NAME: testOrderBndry_Set1
 *
 * DESCRIPTION:
 * 		This function is used to test the orderBndry function.
 *
 * 		This function runs the test data created by Rich Snow.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testOrderBndry_Set1(TestSet) - A test set with the following metrics:
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
TestSet testOrderBndry_Set1()
{
	double DEG2RAD = M_PI / 180.0;
    int passedCount=0, failedCount=0;
    int errorCount = 0;
    int unverifiedCount = 0;
    int setupFailureCount = 0;
    int testCaseCount = 0;
    ErrorSet err=0;
    int jj, ii, failed = 0;

    Boundary bin, bout;
    Geodesic geo1, geo2, geo3, geo4;
    Locus locus1, locus2;
    Spiral spiral1, spiral2;
    Arc arc1, arc2;
    LLPoint geo1Start, geo2Start, geo3Start, geo4Start;
    LLPoint center, geoEnd, tempLLPoint;
    double latS, lonS, crs1, crs2, crs3, crs4, radius;
    double geocrs, geolen, locDist, locDist2; 
    double startRad, endRad, startAz, endAz, tempdbl;
    Geodesic* geo1out;
    Geodesic* geo2out;
    Geodesic* geo3out;
    Geodesic* geo4out;
    Arc* arc1out;
    Arc* arc2out;
    Locus* loc1out;
    Locus* loc2out;
    Spiral* spiral1out;
    Spiral* spiral2out;


    TestSet set;

    set = newTestSet("testOrderBndry_Set1");

    printf("\nStart testOrderBndry_Set1\n");
    srand(20110531);

    for (ii = 0; ii < 10; ii++)
    {
    //test case 1 -- 4 geos
    latS = DEG2RAD * randLat();
    lonS = DEG2RAD * randLon();
    center.latitude = latS;
    center.longitude = lonS;
    crs1 = DEG2RAD * ((double)((rand() % 89) + 1));
    crs2 = DEG2RAD * ((double)((rand() % 89) + 91));
    crs3 = DEG2RAD * ((double)((rand() % 89) + 181));
    crs4 = DEG2RAD * ((double)((rand() % 89) + 271));
    radius = (double)((rand() % 100) + 1);

    err |= direct(center, crs1, radius, &geo1Start, EPS);
    err |= direct(center, crs2, radius, &geo2Start, EPS);
    err |= direct(center, crs3, radius, &geo3Start, EPS);
    err |= direct(center, crs4, radius, &geo4Start, EPS);
    err |= createGeo(&geo1, geo1Start, geo2Start, SEGMENT, EPS);
    err |= createGeo(&geo2, geo3Start, geo2Start, SEGMENT, EPS);
    err |= createGeo(&geo3, geo3Start, geo4Start, SEGMENT, EPS);
    err |= createGeo(&geo4, geo4Start, geo1Start, SEGMENT, EPS);

    //add the shapes to the boundary out of connected order
    for (jj = 0; jj < 5; jj++)
    {
      err = 0;
      bin = createBndry();
      bout = createBndry();
      if (jj == 0)
      {
        err |= addGeoToBndry(&bin, &geo1);
        err |= addGeoToBndry(&bin, &geo3);
        err |= addGeoToBndry(&bin, &geo2);
        err |= addGeoToBndry(&bin, &geo4);
      }
      else if (jj == 1)
      {
        err |= addGeoToBndry(&bin, &geo1);
        err |= addGeoToBndry(&bin, &geo3);
        err |= addGeoToBndry(&bin, &geo4);
        err |= addGeoToBndry(&bin, &geo2);
      }
      else if (jj == 2)
      {
        err |= addGeoToBndry(&bin, &geo1);
        err |= addGeoToBndry(&bin, &geo4);
        err |= addGeoToBndry(&bin, &geo3);
        err |= addGeoToBndry(&bin, &geo2);
      }
      else if (jj == 3)
      {
        err |= addGeoToBndry(&bin, &geo1);
        err |= addGeoToBndry(&bin, &geo4);
        err |= addGeoToBndry(&bin, &geo2);
        err |= addGeoToBndry(&bin, &geo3);
      }
      else if (jj == 4)
      {
        err |= addGeoToBndry(&bin, &geo1);
        err |= addGeoToBndry(&bin, &geo2);
        err |= addGeoToBndry(&bin, &geo4);
        err |= addGeoToBndry(&bin, &geo3);
      }

    if (err) printf("\nError(s) occurred during setup: 0x%lx", err);

    //order the boundary
    err = orderBndry(bin, &bout, TOL, EPS);
    if (err) printf("\nError(s) occurred in orderBndry: 0x%lx", err);

    geo1out = (Geodesic*) &bout.elements[0];
    geo2out = (Geodesic*) &bout.elements[1];
    geo3out = (Geodesic*) &bout.elements[2];
    geo4out = (Geodesic*) &bout.elements[3];

    failed = 0;
    if (!ptsAreSame(geo1out->startPoint, geo1.startPoint, TESTTOL) || !ptsAreSame(geo1out->endPoint, geo1.endPoint, TESTTOL))
      failed = 1;
    //geo2 is reversed by orderBndry
    if (!ptsAreSame(geo2out->startPoint, geo2.endPoint, TESTTOL) || !ptsAreSame(geo2out->endPoint, geo2.startPoint, TESTTOL))
      failed = 1;
    if (!ptsAreSame(geo3out->startPoint, geo3.startPoint, TESTTOL) || !ptsAreSame(geo3out->endPoint, geo3.endPoint, TESTTOL))
      failed = 1;
    if (!ptsAreSame(geo4out->startPoint, geo4.startPoint, TESTTOL) || !ptsAreSame(geo4out->endPoint, geo4.endPoint, TESTTOL))
      failed = 1;

    if (err)
    {
      errorCount++;
      failedCount++;
    }
    else
    {
      if (failed)
        failedCount++;
      else
        passedCount++;
    }

    testCaseCount++;
    }// for jj

    //test case 2 -- locus - arc boundary
    geocrs = DEG2RAD * randAzimuth();
    geolen = (double)((rand() % 100) + 1);
    locDist = (double)((rand() % 20) + 1);

    err |= direct(center, geocrs, geolen, &geoEnd, EPS);
    err |= createLocus(&locus1, center, geoEnd, locDist, locDist, SEGMENT, TOL, EPS);
    err |= createLocus(&locus2, center, geoEnd, -locDist, -locDist, SEGMENT, TOL, EPS);
    err |= createArc(&arc1, center, locus1.locusStart, locus2.locusStart, CLOCKWISE, TOL, EPS);
    err |= createArc(&arc2, geoEnd, locus1.locusEnd, locus2.locusEnd, COUNTERCLOCKWISE, TOL, EPS);

    //add the shapes to the boundary out of connected order
    for (jj = 0; jj < 6; jj++)
    {
      err = 0;
      bin = createBndry();
      bout = createBndry();
      if (jj == 0)
      {
        err |= addLocusToBndry(&bin, &locus1);
        err |= addLocusToBndry(&bin, &locus2);
        err |= addArcToBndry(&bin, &arc1);
        err |= addArcToBndry(&bin, &arc2);
      }
      else if (jj == 1)
      {
        err |= addLocusToBndry(&bin, &locus1);
        err |= addLocusToBndry(&bin, &locus2);
        err |= addArcToBndry(&bin, &arc2);
        err |= addArcToBndry(&bin, &arc1);
      }
      else if (jj == 2)
      {
        err |= addLocusToBndry(&bin, &locus1);
        err |= addArcToBndry(&bin, &arc2);
        err |= addLocusToBndry(&bin, &locus2);
        err |= addArcToBndry(&bin, &arc1);
      }
      else if (jj == 3)
      {
        err |= addLocusToBndry(&bin, &locus1);
        err |= addArcToBndry(&bin, &arc2);
        err |= addArcToBndry(&bin, &arc1);
        err |= addLocusToBndry(&bin, &locus2);
      }
      else if (jj == 4)
      {
        err |= addLocusToBndry(&bin, &locus1);
        err |= addArcToBndry(&bin, &arc1);
        err |= addArcToBndry(&bin, &arc2);
        err |= addLocusToBndry(&bin, &locus2);
      }
      else if (jj == 5)
      {
        err |= addLocusToBndry(&bin, &locus1);
        err |= addArcToBndry(&bin, &arc1);
        err |= addLocusToBndry(&bin, &locus2);
        err |= addArcToBndry(&bin, &arc2);
      }

    if (err) printf("\nError(s) occurred during setup: 0x%lx", err);

    //order the boundary
    err = orderBndry(bin, &bout, TOL, EPS);
    if (err) printf("\nError(s) occurred in orderBndry: 0x%lx", err);

    loc1out = (Locus*) &bout.elements[0];
    arc2out = (Arc*) &bout.elements[1];
    loc2out = (Locus*) &bout.elements[2];
    arc1out = (Arc*) &bout.elements[3];

    //The ordered bndry is locus1, arc2, reversed locus2 (loc2out), and reversed arc1 (arc1out)
    failed = 0;
    if (!ptsAreSame(locus1.geoStart, loc1out->geoStart, TESTTOL) || !ptsAreSame(locus1.geoEnd, loc1out->geoEnd, TESTTOL) || (fabs(locus1.startDist - loc1out->startDist) > TESTTOL) || (fabs(locus1.endDist - loc1out->endDist) > TESTTOL))
      failed = 1;
    if (!ptsAreSame(arc2.centerPoint, arc2out->centerPoint, TESTTOL) || !ptsAreSame(arc2.startPoint, arc2out->startPoint, TESTTOL) || !ptsAreSame(arc2.endPoint, arc2out->endPoint, TESTTOL))
      failed = 1;
    //In the reversed locus2, locDist changes sign
    if (!ptsAreSame(locus2.geoEnd, loc2out->geoStart, TESTTOL) || !ptsAreSame(locus2.geoStart, loc2out->geoEnd, TESTTOL) || (fabs(locus2.startDist + loc2out->startDist) > TESTTOL) || (fabs(locus2.endDist + loc2out->endDist) > TESTTOL))
      failed = 1;
    if (!ptsAreSame(arc1.centerPoint, arc1out->centerPoint, TESTTOL) || !ptsAreSame(arc1.endPoint, arc1out->startPoint, TESTTOL) || !ptsAreSame(arc1.startPoint, arc1out->endPoint, TESTTOL))
      failed = 1;

    if (err)
    {
      errorCount++;
      failedCount++;
    }
    else
    {
      if (failed)
        failedCount++;
      else
        passedCount++;
    }

    testCaseCount++;
    }//for jj

    //test case 3 -- locus - spiral boundary
    //Use locus1 from above
    locDist2 = (double)((rand() % 20) + 1);
    err |= createLocus(&locus2, center, geoEnd, -locDist2, -locDist2, SEGMENT, TOL, EPS);
    err |= inverse(center, locus1.locusStart, &startAz, &tempdbl, &startRad, EPS);
    err |= inverse(center, locus2.locusStart, &endAz, &tempdbl, &endRad, EPS);
    err |= createSpiral(&spiral1, center, startRad, endRad, startAz, endAz, CLOCKWISE, EPS);
    err |= inverse(geoEnd, locus1.locusEnd, &startAz, &tempdbl, &startRad, EPS);
    err |= inverse(geoEnd, locus2.locusEnd, &endAz, &tempdbl, &endRad, EPS);
    err |= createSpiral(&spiral2, geoEnd, startRad, endRad, startAz, endAz, COUNTERCLOCKWISE, EPS);

    //add the shapes to the boundary out of connected order
    for (jj = 0; jj < 6; jj++)
    {
      err = 0;
      bin = createBndry();
      bout = createBndry();
      if (jj == 0)
      {
        err |= addLocusToBndry(&bin, &locus1);
        err |= addLocusToBndry(&bin, &locus2);
        err |= addSpiralToBndry(&bin, &spiral1);
        err |= addSpiralToBndry(&bin, &spiral2);
      }
      else if (jj == 1)
      {
        err |= addLocusToBndry(&bin, &locus1);
        err |= addSpiralToBndry(&bin, &spiral1);
        err |= addLocusToBndry(&bin, &locus2);
        err |= addSpiralToBndry(&bin, &spiral2);
      }
      else if (jj == 2)
      {
        err |= addLocusToBndry(&bin, &locus1);
        err |= addSpiralToBndry(&bin, &spiral1);
        err |= addSpiralToBndry(&bin, &spiral2);
        err |= addLocusToBndry(&bin, &locus2);
      }
      else if (jj == 3)
      {
        err |= addLocusToBndry(&bin, &locus1);
        err |= addSpiralToBndry(&bin, &spiral2);
        err |= addSpiralToBndry(&bin, &spiral1);
        err |= addLocusToBndry(&bin, &locus2);
      }
      else if (jj == 4)
      {
        err |= addLocusToBndry(&bin, &locus1);
        err |= addSpiralToBndry(&bin, &spiral2);
        err |= addLocusToBndry(&bin, &locus2);
        err |= addSpiralToBndry(&bin, &spiral1);
      }
      else if (jj == 5)
      {
        err |= addLocusToBndry(&bin, &locus1);
        err |= addLocusToBndry(&bin, &locus2);
        err |= addSpiralToBndry(&bin, &spiral2);
        err |= addSpiralToBndry(&bin, &spiral1);
      }

    if (err) printf("\nError(s) occurred during setup: 0x%lx", err);

    //order the boundary
    err = orderBndry(bin, &bout, TOL, EPS);
    if (err) printf("\nError(s) occurred in orderBndry: 0x%lx", err);

    loc1out = (Locus*) &bout.elements[0];
    spiral2out = (Spiral*) &bout.elements[1];
    loc2out = (Locus*) &bout.elements[2];
    spiral1out = (Spiral*) &bout.elements[3];

    //The ordered bndry is locus1, spiral2, reversed locus2 (loc2out), and reversed spiral1 (spiral1out)
    failed = 0;
    if (!ptsAreSame(locus1.geoStart, loc1out->geoStart, TESTTOL) || !ptsAreSame(locus1.geoEnd, loc1out->geoEnd, TESTTOL) || (fabs(locus1.startDist - loc1out->startDist) > TESTTOL) || (fabs(locus1.endDist - loc1out->endDist) > TESTTOL))
      failed = 1;
    if (!ptsAreSame(spiral2.centerPoint, spiral2out->centerPoint, TESTTOL) || !ptsAreSame(spiral2.startPoint, spiral2out->startPoint, TESTTOL) || !ptsAreSame(spiral2.endPoint, spiral2out->endPoint, TESTTOL))
      failed = 1;
    //In the reversed locus2, locDist changes sign
    if (!ptsAreSame(locus2.geoEnd, loc2out->geoStart, TESTTOL) || !ptsAreSame(locus2.geoStart, loc2out->geoEnd, TESTTOL) || (fabs(locus2.startDist + loc2out->startDist) > TESTTOL) || (fabs(locus2.endDist + loc2out->endDist) > TESTTOL))
      failed = 1;
    if (!ptsAreSame(spiral1.centerPoint, spiral1out->centerPoint, TESTTOL) || !ptsAreSame(spiral1.endPoint, spiral1out->startPoint, TESTTOL) || !ptsAreSame(spiral1.startPoint, spiral1out->endPoint, TESTTOL))
      failed = 1;

    if (err)
    {
      errorCount++;
      failedCount++;
    }
    else
    {
      if (failed)
        failedCount++;
      else
        passedCount++;
    }

    testCaseCount++;
    }//for jj

    //test case 4 -- spiral and geo boundary (3 shapes)
    //Use spiral2 from above
    err = 0;
    direct(spiral2.startPoint, spiral2.startAz + M_PI_4, (spiral2.startRadius + spiral2.endRadius)/4.0, &tempLLPoint, EPS);
    err |= createGeo(&geo1, spiral2.startPoint, tempLLPoint, SEGMENT, EPS);
    err |= createGeo(&geo2, spiral2.endPoint, tempLLPoint, SEGMENT, EPS);
    bin = createBndry();
    bout = createBndry();
    //add the shapes to the boundary 
    err |= addSpiralToBndry(&bin, &spiral2);
    err |= addGeoToBndry(&bin, &geo1);
    err |= addGeoToBndry(&bin, &geo2);
    if (err) printf("\nError(s) occurred during setup: 0x%lx", err);

    //order the boundary
    err = orderBndry(bin, &bout, TOL, EPS);
    if (err) printf("\nError(s) occurred in orderBndry: 0x%lx", err);
    
    spiral2out = (Spiral*) &bout.elements[0];
    geo2out = (Geodesic*) &bout.elements[1];
    geo1out = (Geodesic*) &bout.elements[2];

    failed = 0;
    if (!ptsAreSame(spiral2.centerPoint, spiral2out->centerPoint, TESTTOL) || !ptsAreSame(spiral2.startPoint, spiral2out->startPoint, TESTTOL) || !ptsAreSame(spiral2.endPoint, spiral2out->endPoint, TESTTOL))
      failed = 1;
    if (!ptsAreSame(geo2out->startPoint, geo2.startPoint, TESTTOL) || !ptsAreSame(geo2out->endPoint, geo2.endPoint, TESTTOL))
      failed = 1;
     //geo1 is reversed by orderBndry
    if (!ptsAreSame(geo1out->startPoint, geo1.endPoint, TESTTOL) || !ptsAreSame(geo1out->endPoint, geo1.startPoint, TESTTOL))
      failed = 1;

    if (err)
    {
      errorCount++;
      failedCount++;
    }
    else
    {
      if (failed)
        failedCount++;
      else
        passedCount++;
    }

    testCaseCount++;

    //test case 5 -- spiral and geo (only 2 shapes)
    //Use spiral2 from above
    err = 0;
    err |= createGeo(&geo1, spiral2.startPoint, spiral2.endPoint, SEGMENT, EPS);
    bin = createBndry();
    bout = createBndry();
    //add the shapes to the boundary 
    err |= addSpiralToBndry(&bin, &spiral2);
    err |= addGeoToBndry(&bin, &geo1);
    if (err) printf("\nError(s) occurred during setup: 0x%lx", err);

    //order the boundary
    err = orderBndry(bin, &bout, TOL, EPS);
    if (err) printf("\nError(s) occurred in orderBndry: 0x%lx", err);
    
    spiral2out = (Spiral*) &bout.elements[0];
    geo1out = (Geodesic*) &bout.elements[1];

    failed = 0;
    if (!ptsAreSame(spiral2.centerPoint, spiral2out->centerPoint, TESTTOL) || !ptsAreSame(spiral2.startPoint, spiral2out->startPoint, TESTTOL) || !ptsAreSame(spiral2.endPoint, spiral2out->endPoint, TESTTOL))
      failed = 1;
     //geo1 is reversed by orderBndry
    if (!ptsAreSame(geo1out->startPoint, geo1.endPoint, TESTTOL) || !ptsAreSame(geo1out->endPoint, geo1.startPoint, TESTTOL))
      failed = 1;

    if (err)
    {
      errorCount++;
      failedCount++;
    }
    else
    {
      if (failed)
        failedCount++;
      else
        passedCount++;
    }

    testCaseCount++;

    }//for ii

    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("\nFinish testOrderBndry_Set1\n\n\n");

    return set;
}

/*
 * NAME: testOrderBndry_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the boundary's orderBndry function.
 *
 * 		This function runs all the test cases and returns
 * 		the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testOrderBndry_AllSets(TestSuite) - A test suite with the folling metrics:
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
TestSuite testOrderBndry_AllSets()
{
	TestSuite suite;
	TestSet set1;

    suite = newTestSuite("testOrderBndry_AllSets");

    printf("\nStart testOrderBndry_AllSets\n");

    set1 = testOrderBndry_Set1();
    addTestSet(set1,&suite);

    displayTestSuite(suite);

    printf("\nFinish testOrderBndry_AllSets\n\n\n");

    return suite;
}

/*
 * NAME: testSeparateBndry_Set1
 *
 * DESCRIPTION:
 * 		This function is used to test the separateBndry function.
 *
 * 		This function runs the test data created by Rich Snow.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testSeparateBndry_Set1(TestSet) - A test set with the following metrics:
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
TestSet testSeparateBndry_Set1()
{
	double DEG2RAD = M_PI / 180.0;
    int passedCount=0, failedCount=0;
    int errorCount = 0;
    int unverifiedCount = 0;
    int setupFailureCount = 0;
    int testCaseCount = 0;
    ErrorSet err=0;
    int jj, ii, failed = 0;

    Boundary bin;
    Boundary bout[8];
    Geodesic geo1, geo2, geo3, geo4;
    Locus locus1, locus2;
    Spiral spiral1, spiral2;
    Arc arc1, arc2;
    LLPoint geo1Start, geo2Start, geo3Start, geo4Start;
    LLPoint center, geoEnd, tempLLPoint;
    double latS, lonS, crs1, crs2, crs3, crs4, radius;
    double geocrs, geolen, locDist, locDist2;
    double startRad, endRad, startAz, endAz, tempdbl;
    Geodesic* geo1out;
    Geodesic* geo2out;
    Geodesic* geo3out;
    Geodesic* geo4out;
    Arc* arc1out;
    Arc* arc2out;
    Locus* loc1out;
    Locus* loc2out;
    Spiral* spiral1out;
    Spiral* spiral2out;
    int numberOfBoundaries = 0;

    TestSet set;

    set = newTestSet("testSeparateBndry_Set1");

    printf("\nStart testSeparateBndry_Set1\n");
    srand(20110531);

    for (ii = 0; ii < 10; ii++)
    {
		//test case 1 -- 4 geos
		latS = DEG2RAD * randLat();
		lonS = DEG2RAD * randLon();
		center.latitude = latS;
		center.longitude = lonS;
		crs1 = DEG2RAD * ((double)((rand() % 89) + 1));
		crs2 = DEG2RAD * ((double)((rand() % 89) + 91));
		crs3 = DEG2RAD * ((double)((rand() % 89) + 181));
		crs4 = DEG2RAD * ((double)((rand() % 89) + 271));
		radius = (double)((rand() % 100) + 1);

		err |= direct(center, crs1, radius, &geo1Start, EPS);
		err |= direct(center, crs2, radius, &geo2Start, EPS);
		err |= direct(center, crs3, radius, &geo3Start, EPS);
		err |= direct(center, crs4, radius, &geo4Start, EPS);
		err |= createGeo(&geo1, geo1Start, geo2Start, SEGMENT, EPS);
		err |= createGeo(&geo2, geo3Start, geo2Start, SEGMENT, EPS);
		err |= createGeo(&geo3, geo3Start, geo4Start, SEGMENT, EPS);
		err |= createGeo(&geo4, geo4Start, geo1Start, SEGMENT, EPS);

		//add the shapes to the boundary out of connected order
		for (jj = 0; jj < 5; jj++)
		{
			err = 0;
			bin = createBndry();
			bout[0] = createBndry();
			if (jj == 0)
			{
				err |= addGeoToBndry(&bin, &geo1);
				err |= addGeoToBndry(&bin, &geo3);
				err |= addGeoToBndry(&bin, &geo2);
				err |= addGeoToBndry(&bin, &geo4);
			}
			else if (jj == 1)
			{
				err |= addGeoToBndry(&bin, &geo1);
				err |= addGeoToBndry(&bin, &geo3);
				err |= addGeoToBndry(&bin, &geo4);
				err |= addGeoToBndry(&bin, &geo2);
			}
			else if (jj == 2)
			{
				err |= addGeoToBndry(&bin, &geo1);
				err |= addGeoToBndry(&bin, &geo4);
				err |= addGeoToBndry(&bin, &geo3);
				err |= addGeoToBndry(&bin, &geo2);
			}
			else if (jj == 3)
			{
				err |= addGeoToBndry(&bin, &geo1);
				err |= addGeoToBndry(&bin, &geo4);
				err |= addGeoToBndry(&bin, &geo2);
				err |= addGeoToBndry(&bin, &geo3);
			}
			else if (jj == 4)
			{
				err |= addGeoToBndry(&bin, &geo1);
				err |= addGeoToBndry(&bin, &geo2);
				err |= addGeoToBndry(&bin, &geo4);
				err |= addGeoToBndry(&bin, &geo3);
			}

			if (err) printf("\nError(s) occurred during setup: 0x%lx", err);

			//order the boundary
			err = separateBndry(bin, bout, &numberOfBoundaries, TOL, EPS);
			if (err) printf("\nError(s) occurred in separateBndry: 0x%lx", err);

			geo1out = (Geodesic*) &bout[0].elements[0];
			geo2out = (Geodesic*) &bout[0].elements[1];
			geo3out = (Geodesic*) &bout[0].elements[2];
			geo4out = (Geodesic*) &bout[0].elements[3];

			failed = 0;
			if (!ptsAreSame(geo1out->startPoint, geo1.startPoint, TESTTOL) || !ptsAreSame(geo1out->endPoint, geo1.endPoint, TESTTOL))
				failed = 1;
			//geo2 is reversed by separateBndry
			if (!ptsAreSame(geo2out->startPoint, geo2.endPoint, TESTTOL) || !ptsAreSame(geo2out->endPoint, geo2.startPoint, TESTTOL))
				failed = 1;
			if (!ptsAreSame(geo3out->startPoint, geo3.startPoint, TESTTOL) || !ptsAreSame(geo3out->endPoint, geo3.endPoint, TESTTOL))
				failed = 1;
			if (!ptsAreSame(geo4out->startPoint, geo4.startPoint, TESTTOL) || !ptsAreSame(geo4out->endPoint, geo4.endPoint, TESTTOL))
				failed = 1;

			if (err)
			{
				errorCount++;
				failedCount++;
			}
			else
			{
				if (failed)
					failedCount++;
				else
					passedCount++;
			}

			testCaseCount++;
		}// for jj

		//test case 2 -- locus - arc boundary
		geocrs = DEG2RAD * randAzimuth();
		geolen = (double)((rand() % 100) + 1);
		locDist = (double)((rand() % 20) + 1);

		err |= direct(center, geocrs, geolen, &geoEnd, EPS);
		err |= createLocus(&locus1, center, geoEnd, locDist, locDist, SEGMENT, TOL, EPS);
		err |= createLocus(&locus2, center, geoEnd, -locDist, -locDist, SEGMENT, TOL, EPS);
		err |= createArc(&arc1, center, locus1.locusStart, locus2.locusStart, CLOCKWISE, TOL, EPS);
		err |= createArc(&arc2, geoEnd, locus1.locusEnd, locus2.locusEnd, COUNTERCLOCKWISE, TOL, EPS);

		//add the shapes to the boundary out of connected order
		for (jj = 0; jj < 6; jj++)
		{
			err = 0;
			bin = createBndry();
			bout[0] = createBndry();
			if (jj == 0)
			{
				err |= addLocusToBndry(&bin, &locus1);
				err |= addLocusToBndry(&bin, &locus2);
				err |= addArcToBndry(&bin, &arc1);
				err |= addArcToBndry(&bin, &arc2);
			}
			else if (jj == 1)
			{
				err |= addLocusToBndry(&bin, &locus1);
				err |= addLocusToBndry(&bin, &locus2);
				err |= addArcToBndry(&bin, &arc2);
				err |= addArcToBndry(&bin, &arc1);
			}
			else if (jj == 2)
			{
				err |= addLocusToBndry(&bin, &locus1);
				err |= addArcToBndry(&bin, &arc2);
				err |= addLocusToBndry(&bin, &locus2);
				err |= addArcToBndry(&bin, &arc1);
			}
			else if (jj == 3)
			{
				err |= addLocusToBndry(&bin, &locus1);
				err |= addArcToBndry(&bin, &arc2);
				err |= addArcToBndry(&bin, &arc1);
				err |= addLocusToBndry(&bin, &locus2);
			}
			else if (jj == 4)
			{
				err |= addLocusToBndry(&bin, &locus1);
				err |= addArcToBndry(&bin, &arc1);
				err |= addArcToBndry(&bin, &arc2);
				err |= addLocusToBndry(&bin, &locus2);
			}
			else if (jj == 5)
			{
				err |= addLocusToBndry(&bin, &locus1);
				err |= addArcToBndry(&bin, &arc1);
				err |= addLocusToBndry(&bin, &locus2);
				err |= addArcToBndry(&bin, &arc2);
			}

			if (err) printf("\nError(s) occurred during setup: 0x%lx", err);

			//order the boundary
			err = separateBndry(bin, bout, &numberOfBoundaries, TOL, EPS);
			if (err) printf("\nError(s) occurred in separateBndry: 0x%lx", err);

			loc1out = (Locus*) &bout[0].elements[0];
			arc2out = (Arc*) &bout[0].elements[1];
			loc2out = (Locus*) &bout[0].elements[2];
			arc1out = (Arc*) &bout[0].elements[3];

			//The ordered bndry is locus1, arc2, reversed locus2 (loc2out), and reversed arc1 (arc1out)
			failed = 0;
			if (!ptsAreSame(locus1.geoStart, loc1out->geoStart, TESTTOL) || !ptsAreSame(locus1.geoEnd, loc1out->geoEnd, TESTTOL) || (fabs(locus1.startDist - loc1out->startDist) > TESTTOL) || (fabs(locus1.endDist - loc1out->endDist) > TESTTOL))
				failed = 1;
			if (!ptsAreSame(arc2.centerPoint, arc2out->centerPoint, TESTTOL) || !ptsAreSame(arc2.startPoint, arc2out->startPoint, TESTTOL) || !ptsAreSame(arc2.endPoint, arc2out->endPoint, TESTTOL))
				failed = 1;
			//In the reversed locus2, locDist changes sign
			if (!ptsAreSame(locus2.geoEnd, loc2out->geoStart, TESTTOL) || !ptsAreSame(locus2.geoStart, loc2out->geoEnd, TESTTOL) || (fabs(locus2.startDist + loc2out->startDist) > TESTTOL) || (fabs(locus2.endDist + loc2out->endDist) > TESTTOL))
				failed = 1;
			if (!ptsAreSame(arc1.centerPoint, arc1out->centerPoint, TESTTOL) || !ptsAreSame(arc1.endPoint, arc1out->startPoint, TESTTOL) || !ptsAreSame(arc1.startPoint, arc1out->endPoint, TESTTOL))
				failed = 1;

			if (err)
			{
				errorCount++;
				failedCount++;
			}
			else
			{
				if (failed)
					failedCount++;
				else
					passedCount++;
			}

			testCaseCount++;
		}//for jj

		//test case 3 -- locus - spiral boundary
		//Use locus1 from above
		locDist2 = (double)((rand() % 20) + 1);
		err |= createLocus(&locus2, center, geoEnd, -locDist2, -locDist2, SEGMENT, TOL, EPS);
		err |= inverse(center, locus1.locusStart, &startAz, &tempdbl, &startRad, EPS);
		err |= inverse(center, locus2.locusStart, &endAz, &tempdbl, &endRad, EPS);
		err |= createSpiral(&spiral1, center, startRad, endRad, startAz, endAz, CLOCKWISE, EPS);
		err |= inverse(geoEnd, locus1.locusEnd, &startAz, &tempdbl, &startRad, EPS);
		err |= inverse(geoEnd, locus2.locusEnd, &endAz, &tempdbl, &endRad, EPS);
		err |= createSpiral(&spiral2, geoEnd, startRad, endRad, startAz, endAz, COUNTERCLOCKWISE, EPS);

		//add the shapes to the boundary out of connected order
		for (jj = 0; jj < 6; jj++)
		{
			err = 0;
			bin = createBndry();
			bout[0] = createBndry();
			if (jj == 0)
			{
				err |= addLocusToBndry(&bin, &locus1);
				err |= addLocusToBndry(&bin, &locus2);
				err |= addSpiralToBndry(&bin, &spiral1);
				err |= addSpiralToBndry(&bin, &spiral2);
			}
			else if (jj == 1)
			{
				err |= addLocusToBndry(&bin, &locus1);
				err |= addSpiralToBndry(&bin, &spiral1);
				err |= addLocusToBndry(&bin, &locus2);
				err |= addSpiralToBndry(&bin, &spiral2);
			}
			else if (jj == 2)
			{
				err |= addLocusToBndry(&bin, &locus1);
				err |= addSpiralToBndry(&bin, &spiral1);
				err |= addSpiralToBndry(&bin, &spiral2);
				err |= addLocusToBndry(&bin, &locus2);
			}
			else if (jj == 3)
			{
				err |= addLocusToBndry(&bin, &locus1);
				err |= addSpiralToBndry(&bin, &spiral2);
				err |= addSpiralToBndry(&bin, &spiral1);
				err |= addLocusToBndry(&bin, &locus2);
			}
			else if (jj == 4)
			{
				err |= addLocusToBndry(&bin, &locus1);
				err |= addSpiralToBndry(&bin, &spiral2);
				err |= addLocusToBndry(&bin, &locus2);
				err |= addSpiralToBndry(&bin, &spiral1);
			}
			else if (jj == 5)
			{
				err |= addLocusToBndry(&bin, &locus1);
				err |= addLocusToBndry(&bin, &locus2);
				err |= addSpiralToBndry(&bin, &spiral2);
				err |= addSpiralToBndry(&bin, &spiral1);
			}

			if (err) printf("\nError(s) occurred during setup: 0x%lx", err);

			//order the boundary
			err = separateBndry(bin, bout, &numberOfBoundaries, TOL, EPS);
			if (err) printf("\nError(s) occurred in separateBndry: 0x%lx", err);

			loc1out = (Locus*) &bout[0].elements[0];
			spiral2out = (Spiral*) &bout[0].elements[1];
			loc2out = (Locus*) &bout[0].elements[2];
			spiral1out = (Spiral*) &bout[0].elements[3];

			//The ordered bndry is locus1, spiral2, reversed locus2 (loc2out), and reversed spiral1 (spiral1out)
			failed = 0;
			if (!ptsAreSame(locus1.geoStart, loc1out->geoStart, TESTTOL) || !ptsAreSame(locus1.geoEnd, loc1out->geoEnd, TESTTOL) || (fabs(locus1.startDist - loc1out->startDist) > TESTTOL) || (fabs(locus1.endDist - loc1out->endDist) > TESTTOL))
				failed = 1;
			if (!ptsAreSame(spiral2.centerPoint, spiral2out->centerPoint, TESTTOL) || !ptsAreSame(spiral2.startPoint, spiral2out->startPoint, TESTTOL) || !ptsAreSame(spiral2.endPoint, spiral2out->endPoint, TESTTOL))
				failed = 1;
			//In the reversed locus2, locDist changes sign
			if (!ptsAreSame(locus2.geoEnd, loc2out->geoStart, TESTTOL) || !ptsAreSame(locus2.geoStart, loc2out->geoEnd, TESTTOL) || (fabs(locus2.startDist + loc2out->startDist) > TESTTOL) || (fabs(locus2.endDist + loc2out->endDist) > TESTTOL))
				failed = 1;
			if (!ptsAreSame(spiral1.centerPoint, spiral1out->centerPoint, TESTTOL) || !ptsAreSame(spiral1.endPoint, spiral1out->startPoint, TESTTOL) || !ptsAreSame(spiral1.startPoint, spiral1out->endPoint, TESTTOL))
				failed = 1;

			if (err)
			{
			  errorCount++;
			  failedCount++;
			}
			else
			{
				if (failed)
					failedCount++;
				else
					passedCount++;
			}

			testCaseCount++;
		}//for jj

		//test case 4 -- spiral and geo boundary (3 shapes)
		//Use spiral2 from above
		err = 0;
		direct(spiral2.startPoint, spiral2.startAz + M_PI_4, (spiral2.startRadius + spiral2.endRadius)/4.0, &tempLLPoint, EPS);
		err |= createGeo(&geo1, spiral2.startPoint, tempLLPoint, SEGMENT, EPS);
		err |= createGeo(&geo2, spiral2.endPoint, tempLLPoint, SEGMENT, EPS);
		bin = createBndry();
		bout[0] = createBndry();
		//add the shapes to the boundary
		err |= addSpiralToBndry(&bin, &spiral2);
		err |= addGeoToBndry(&bin, &geo1);
		err |= addGeoToBndry(&bin, &geo2);
		if (err) printf("\nError(s) occurred during setup: 0x%lx", err);

		//order the boundary
		err = separateBndry(bin, bout, &numberOfBoundaries, TOL, EPS);
		if (err) printf("\nError(s) occurred in separateBndry: 0x%lx", err);

		spiral2out = (Spiral*) &bout[0].elements[0];
		geo2out = (Geodesic*) &bout[0].elements[1];
		geo1out = (Geodesic*) &bout[0].elements[2];

		failed = 0;
		if (!ptsAreSame(spiral2.centerPoint, spiral2out->centerPoint, TESTTOL) || !ptsAreSame(spiral2.startPoint, spiral2out->startPoint, TESTTOL) || !ptsAreSame(spiral2.endPoint, spiral2out->endPoint, TESTTOL))
			failed = 1;
		if (!ptsAreSame(geo2out->startPoint, geo2.startPoint, TESTTOL) || !ptsAreSame(geo2out->endPoint, geo2.endPoint, TESTTOL))
			failed = 1;
		 //geo1 is reversed by separateBndry
		if (!ptsAreSame(geo1out->startPoint, geo1.endPoint, TESTTOL) || !ptsAreSame(geo1out->endPoint, geo1.startPoint, TESTTOL))
			failed = 1;

		if (err)
		{
			errorCount++;
			failedCount++;
		}
		else
		{
			if (failed)
				failedCount++;
			else
				passedCount++;
		}

		testCaseCount++;

		//test case 5 -- spiral and geo (only 2 shapes)
		//Use spiral2 from above
		err = 0;
		err |= createGeo(&geo1, spiral2.startPoint, spiral2.endPoint, SEGMENT, EPS);
		bin = createBndry();
		bout[0] = createBndry();
		//add the shapes to the boundary
		err |= addSpiralToBndry(&bin, &spiral2);
		err |= addGeoToBndry(&bin, &geo1);
		if (err) printf("\nError(s) occurred during setup: 0x%lx", err);

		//order the boundary
		err = separateBndry(bin, bout, &numberOfBoundaries, TOL, EPS);
		if (err) printf("\nError(s) occurred in separateBndry: 0x%lx", err);

		spiral2out = (Spiral*) &bout[0].elements[0];
		geo1out = (Geodesic*) &bout[0].elements[1];

		failed = 0;
		if (!ptsAreSame(spiral2.centerPoint, spiral2out->centerPoint, TESTTOL) || !ptsAreSame(spiral2.startPoint, spiral2out->startPoint, TESTTOL) || !ptsAreSame(spiral2.endPoint, spiral2out->endPoint, TESTTOL))
			failed = 1;
		 //geo1 is reversed by separateBndry
		if (!ptsAreSame(geo1out->startPoint, geo1.endPoint, TESTTOL) || !ptsAreSame(geo1out->endPoint, geo1.startPoint, TESTTOL))
			failed = 1;

		if (err)
		{
			errorCount++;
			failedCount++;
		}
		else
		{
			if (failed)
				failedCount++;
			else
				passedCount++;
		}

		testCaseCount++;

    }//for ii

    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("\nFinish testSeparateBndry_Set1\n\n\n");

    return set;
}

/*
 * NAME: testSeparateBndry_Set2
 *
 * DESCRIPTION:
 * 		This function is used to test the separateBndry function.
 *
 * 		This function runs the test data created by Rich Snow.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testSeparateBndry_Set1(TestSet) - A test set with the following metrics:
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
TestSet testSeparateBndry_Set2()
{
	double DEG2RAD = M_PI / 180.0;
    int passedCount = 0, failedCount = 0;
    int errorCount = 0;
    int unverifiedCount = 0;
    int setupFailureCount = 0;
    int testCaseCount = 0;
    ErrorSet err = 0;
    int ii, failed = 0;

    Boundary bin;
    Boundary bout[8];
    Geodesic geo11, geo12, geo13, geo14;
    Geodesic geo21, geo22, geo23, geo24;
    LLPoint center1;
    LLPoint center2;
    double latS1, lonS1;
    double latS2, lonS2;
    double crs11, crs12, crs13, crs14;
    double crs21, crs22, crs23, crs24;
    LLPoint geo11Start, geo12Start, geo13Start, geo14Start;
    LLPoint geo21Start, geo22Start, geo23Start, geo24Start;
    double radius1;
    double radius2;
    Geodesic* geo1out;
    Geodesic* geo2out;
    Geodesic* geo3out;
    Geodesic* geo4out;
    Geodesic* geo5out;
    Geodesic* geo6out;
    Geodesic* geo7out;
    Geodesic* geo8out;

    int numberOfBoundaries = 0;

    TestSet set;

    set = newTestSet("testSeparateBndry_Set2");

    printf("\nStart testSeparateBndry_Set2\n");
    srand(20110531);

    for (ii = 0; ii < 100; ii++)
    {

		err = 0;
		bin = createBndry();
		bout[0] = createBndry();
		bout[1] = createBndry();

		latS1 = DEG2RAD * randLat();
		lonS1 = DEG2RAD * randLon();
		center1.latitude = latS1;
		center1.longitude = lonS1;
		crs11 = DEG2RAD * ((double)((rand() % 89) + 1));
		crs12 = DEG2RAD * ((double)((rand() % 89) + 91));
		crs13 = DEG2RAD * ((double)((rand() % 89) + 181));
		crs14 = DEG2RAD * ((double)((rand() % 89) + 271));
		radius1 = (double)((rand() % 100) + 1);

		latS2 = DEG2RAD * randLat();
		lonS2 = DEG2RAD * randLon();
		center2.latitude = latS2;
		center2.longitude = lonS2;
		crs21 = DEG2RAD * ((double)((rand() % 89) + 1));
		crs22 = DEG2RAD * ((double)((rand() % 89) + 91));
		crs23 = DEG2RAD * ((double)((rand() % 89) + 181));
		crs24 = DEG2RAD * ((double)((rand() % 89) + 271));
		radius2 = (double)((rand() % 100) + 1);


		err |= direct(center1, crs11, radius1, &geo11Start, EPS);
		err |= direct(center1, crs12, radius1, &geo12Start, EPS);
		err |= direct(center1, crs13, radius1, &geo13Start, EPS);
		err |= direct(center1, crs14, radius1, &geo14Start, EPS);
		err |= createGeo(&geo11, geo11Start, geo12Start, SEGMENT, EPS);
		err |= createGeo(&geo12, geo13Start, geo12Start, SEGMENT, EPS); // Flipped
		err |= createGeo(&geo13, geo13Start, geo14Start, SEGMENT, EPS);
		err |= createGeo(&geo14, geo14Start, geo11Start, SEGMENT, EPS);

		err |= direct(center1, crs21, radius2, &geo21Start, EPS);
		err |= direct(center1, crs22, radius2, &geo22Start, EPS);
		err |= direct(center1, crs23, radius2, &geo23Start, EPS);
		err |= direct(center1, crs24, radius2, &geo24Start, EPS);
		err |= createGeo(&geo21, geo21Start, geo22Start, SEGMENT, EPS);
		err |= createGeo(&geo22, geo23Start, geo22Start, SEGMENT, EPS); // Flipped
		err |= createGeo(&geo23, geo23Start, geo24Start, SEGMENT, EPS);
		err |= createGeo(&geo24, geo24Start, geo21Start, SEGMENT, EPS);

		//add the shapes to the boundary out of connected order
		err |= addGeoToBndry(&bin, &geo11);
		err |= addGeoToBndry(&bin, &geo13);
		err |= addGeoToBndry(&bin, &geo12);
		err |= addGeoToBndry(&bin, &geo14);

		err |= addGeoToBndry(&bin, &geo21);
		err |= addGeoToBndry(&bin, &geo24);
		err |= addGeoToBndry(&bin, &geo22);
		err |= addGeoToBndry(&bin, &geo23);
		if (err) printf("\nError(s) occurred during setup: 0x%lx", err);

		//order the boundary
		err = separateBndry(bin, bout, &numberOfBoundaries, TOL, EPS);
		if (err) printf("\nError(s) occurred in separateBndry: 0x%lx", err);

		geo1out = (Geodesic*) &bout[0].elements[0];
		geo2out = (Geodesic*) &bout[0].elements[1];
		geo3out = (Geodesic*) &bout[0].elements[2];
		geo4out = (Geodesic*) &bout[0].elements[3];

		geo5out = (Geodesic*) &bout[1].elements[0];
		geo6out = (Geodesic*) &bout[1].elements[1];
		geo7out = (Geodesic*) &bout[1].elements[2];
		geo8out = (Geodesic*) &bout[1].elements[3];

		failed = 0;
		if (numberOfBoundaries != 2)
			failed = 1;

		if (!ptsAreSame(geo1out->startPoint, geo11.startPoint, TESTTOL) || !ptsAreSame(geo1out->endPoint, geo11.endPoint, TESTTOL))
			failed = 1;
		//geo2 is reversed by separateBndry
		if (!ptsAreSame(geo2out->startPoint, geo12.endPoint, TESTTOL) || !ptsAreSame(geo2out->endPoint, geo12.startPoint, TESTTOL))
			failed = 1;
		if (!ptsAreSame(geo3out->startPoint, geo13.startPoint, TESTTOL) || !ptsAreSame(geo3out->endPoint, geo13.endPoint, TESTTOL))
			failed = 1;
		if (!ptsAreSame(geo4out->startPoint, geo14.startPoint, TESTTOL) || !ptsAreSame(geo4out->endPoint, geo14.endPoint, TESTTOL))
			failed = 1;

		if (!ptsAreSame(geo5out->startPoint, geo21.startPoint, TESTTOL) || !ptsAreSame(geo5out->endPoint, geo21.endPoint, TESTTOL))
			failed = 1;
		//geo2 is reversed by separateBndry
		if (!ptsAreSame(geo6out->startPoint, geo22.endPoint, TESTTOL) || !ptsAreSame(geo6out->endPoint, geo22.startPoint, TESTTOL))
			failed = 1;
		if (!ptsAreSame(geo7out->startPoint, geo23.startPoint, TESTTOL) || !ptsAreSame(geo7out->endPoint, geo23.endPoint, TESTTOL))
			failed = 1;
		if (!ptsAreSame(geo8out->startPoint, geo24.startPoint, TESTTOL) || !ptsAreSame(geo8out->endPoint, geo24.endPoint, TESTTOL))
			failed = 1;

		if (err)
		{
			errorCount++;
			failedCount++;
		}
		else
		{
			if (failed)
				failedCount++;
			else
				passedCount++;
		}

		testCaseCount++;

//		//test case 2 -- locus - arc boundary
//		geocrs = DEG2RAD * randAzimuth();
//		geolen = (double)((rand() % 100) + 1);
//		locDist = (double)((rand() % 20) + 1);
//
//		err |= direct(center, geocrs, geolen, &geoEnd, EPS);
//		err |= createLocus(&locus1, center, geoEnd, locDist, locDist, SEGMENT, TOL, EPS);
//		err |= createLocus(&locus2, center, geoEnd, -locDist, -locDist, SEGMENT, TOL, EPS);
//		err |= createArc(&arc1, center, locus1.locusStart, locus2.locusStart, CLOCKWISE, TOL, EPS);
//		err |= createArc(&arc2, geoEnd, locus1.locusEnd, locus2.locusEnd, COUNTERCLOCKWISE, TOL, EPS);
//
//		//add the shapes to the boundary out of connected order
//		for (jj = 0; jj < 6; jj++)
//		{
//			err = 0;
//			bin = createBndry();
//			bout[0] = createBndry();
//			if (jj == 0)
//			{
//				err |= addLocusToBndry(&bin, &locus1);
//				err |= addLocusToBndry(&bin, &locus2);
//				err |= addArcToBndry(&bin, &arc1);
//				err |= addArcToBndry(&bin, &arc2);
//			}
//			else if (jj == 1)
//			{
//				err |= addLocusToBndry(&bin, &locus1);
//				err |= addLocusToBndry(&bin, &locus2);
//				err |= addArcToBndry(&bin, &arc2);
//				err |= addArcToBndry(&bin, &arc1);
//			}
//			else if (jj == 2)
//			{
//				err |= addLocusToBndry(&bin, &locus1);
//				err |= addArcToBndry(&bin, &arc2);
//				err |= addLocusToBndry(&bin, &locus2);
//				err |= addArcToBndry(&bin, &arc1);
//			}
//			else if (jj == 3)
//			{
//				err |= addLocusToBndry(&bin, &locus1);
//				err |= addArcToBndry(&bin, &arc2);
//				err |= addArcToBndry(&bin, &arc1);
//				err |= addLocusToBndry(&bin, &locus2);
//			}
//			else if (jj == 4)
//			{
//				err |= addLocusToBndry(&bin, &locus1);
//				err |= addArcToBndry(&bin, &arc1);
//				err |= addArcToBndry(&bin, &arc2);
//				err |= addLocusToBndry(&bin, &locus2);
//			}
//			else if (jj == 5)
//			{
//				err |= addLocusToBndry(&bin, &locus1);
//				err |= addArcToBndry(&bin, &arc1);
//				err |= addLocusToBndry(&bin, &locus2);
//				err |= addArcToBndry(&bin, &arc2);
//			}
//
//			if (err) printf("\nError(s) occurred during setup: 0x%lx", err);
//
//			//order the boundary
//			err = separateBndry(bin, bout, &numberOfBoundaries, TOL, EPS);
//			if (err) printf("\nError(s) occurred in separateBndry: 0x%lx", err);
//
//			loc1out = (Locus*) &bout[0].elements[0];
//			arc2out = (Arc*) &bout[0].elements[1];
//			loc2out = (Locus*) &bout[0].elements[2];
//			arc1out = (Arc*) &bout[0].elements[3];
//
//			//The ordered bndry is locus1, arc2, reversed locus2 (loc2out), and reversed arc1 (arc1out)
//			failed = 0;
//			if (!ptsAreSame(locus1.geoStart, loc1out->geoStart, TESTTOL) || !ptsAreSame(locus1.geoEnd, loc1out->geoEnd, TESTTOL) || (fabs(locus1.startDist - loc1out->startDist) > TESTTOL) || (fabs(locus1.endDist - loc1out->endDist) > TESTTOL))
//				failed = 1;
//			if (!ptsAreSame(arc2.centerPoint, arc2out->centerPoint, TESTTOL) || !ptsAreSame(arc2.startPoint, arc2out->startPoint, TESTTOL) || !ptsAreSame(arc2.endPoint, arc2out->endPoint, TESTTOL))
//				failed = 1;
//			//In the reversed locus2, locDist changes sign
//			if (!ptsAreSame(locus2.geoEnd, loc2out->geoStart, TESTTOL) || !ptsAreSame(locus2.geoStart, loc2out->geoEnd, TESTTOL) || (fabs(locus2.startDist + loc2out->startDist) > TESTTOL) || (fabs(locus2.endDist + loc2out->endDist) > TESTTOL))
//				failed = 1;
//			if (!ptsAreSame(arc1.centerPoint, arc1out->centerPoint, TESTTOL) || !ptsAreSame(arc1.endPoint, arc1out->startPoint, TESTTOL) || !ptsAreSame(arc1.startPoint, arc1out->endPoint, TESTTOL))
//				failed = 1;
//
//			if (err)
//			{
//				errorCount++;
//				failedCount++;
//			}
//			else
//			{
//				if (failed)
//					failedCount++;
//				else
//					passedCount++;
//			}
//
//			testCaseCount++;
//		}//for jj
//
//		//test case 3 -- locus - spiral boundary
//		//Use locus1 from above
//		locDist2 = (double)((rand() % 20) + 1);
//		err |= createLocus(&locus2, center, geoEnd, -locDist2, -locDist2, SEGMENT, TOL, EPS);
//		err |= inverse(center, locus1.locusStart, &startAz, &tempdbl, &startRad, EPS);
//		err |= inverse(center, locus2.locusStart, &endAz, &tempdbl, &endRad, EPS);
//		err |= createSpiral(&spiral1, center, startRad, endRad, startAz, endAz, CLOCKWISE, EPS);
//		err |= inverse(geoEnd, locus1.locusEnd, &startAz, &tempdbl, &startRad, EPS);
//		err |= inverse(geoEnd, locus2.locusEnd, &endAz, &tempdbl, &endRad, EPS);
//		err |= createSpiral(&spiral2, geoEnd, startRad, endRad, startAz, endAz, COUNTERCLOCKWISE, EPS);
//
//		//add the shapes to the boundary out of connected order
//		for (jj = 0; jj < 6; jj++)
//		{
//			err = 0;
//			bin = createBndry();
//			bout[0] = createBndry();
//			if (jj == 0)
//			{
//				err |= addLocusToBndry(&bin, &locus1);
//				err |= addLocusToBndry(&bin, &locus2);
//				err |= addSpiralToBndry(&bin, &spiral1);
//				err |= addSpiralToBndry(&bin, &spiral2);
//			}
//			else if (jj == 1)
//			{
//				err |= addLocusToBndry(&bin, &locus1);
//				err |= addSpiralToBndry(&bin, &spiral1);
//				err |= addLocusToBndry(&bin, &locus2);
//				err |= addSpiralToBndry(&bin, &spiral2);
//			}
//			else if (jj == 2)
//			{
//				err |= addLocusToBndry(&bin, &locus1);
//				err |= addSpiralToBndry(&bin, &spiral1);
//				err |= addSpiralToBndry(&bin, &spiral2);
//				err |= addLocusToBndry(&bin, &locus2);
//			}
//			else if (jj == 3)
//			{
//				err |= addLocusToBndry(&bin, &locus1);
//				err |= addSpiralToBndry(&bin, &spiral2);
//				err |= addSpiralToBndry(&bin, &spiral1);
//				err |= addLocusToBndry(&bin, &locus2);
//			}
//			else if (jj == 4)
//			{
//				err |= addLocusToBndry(&bin, &locus1);
//				err |= addSpiralToBndry(&bin, &spiral2);
//				err |= addLocusToBndry(&bin, &locus2);
//				err |= addSpiralToBndry(&bin, &spiral1);
//			}
//			else if (jj == 5)
//			{
//				err |= addLocusToBndry(&bin, &locus1);
//				err |= addLocusToBndry(&bin, &locus2);
//				err |= addSpiralToBndry(&bin, &spiral2);
//				err |= addSpiralToBndry(&bin, &spiral1);
//			}
//
//			if (err) printf("\nError(s) occurred during setup: 0x%lx", err);
//
//			//order the boundary
//			err = separateBndry(bin, bout, &numberOfBoundaries, TOL, EPS);
//			if (err) printf("\nError(s) occurred in separateBndry: 0x%lx", err);
//
//			loc1out = (Locus*) &bout[0].elements[0];
//			spiral2out = (Spiral*) &bout[0].elements[1];
//			loc2out = (Locus*) &bout[0].elements[2];
//			spiral1out = (Spiral*) &bout[0].elements[3];
//
//			//The ordered bndry is locus1, spiral2, reversed locus2 (loc2out), and reversed spiral1 (spiral1out)
//			failed = 0;
//			if (!ptsAreSame(locus1.geoStart, loc1out->geoStart, TESTTOL) || !ptsAreSame(locus1.geoEnd, loc1out->geoEnd, TESTTOL) || (fabs(locus1.startDist - loc1out->startDist) > TESTTOL) || (fabs(locus1.endDist - loc1out->endDist) > TESTTOL))
//				failed = 1;
//			if (!ptsAreSame(spiral2.centerPoint, spiral2out->centerPoint, TESTTOL) || !ptsAreSame(spiral2.startPoint, spiral2out->startPoint, TESTTOL) || !ptsAreSame(spiral2.endPoint, spiral2out->endPoint, TESTTOL))
//				failed = 1;
//			//In the reversed locus2, locDist changes sign
//			if (!ptsAreSame(locus2.geoEnd, loc2out->geoStart, TESTTOL) || !ptsAreSame(locus2.geoStart, loc2out->geoEnd, TESTTOL) || (fabs(locus2.startDist + loc2out->startDist) > TESTTOL) || (fabs(locus2.endDist + loc2out->endDist) > TESTTOL))
//				failed = 1;
//			if (!ptsAreSame(spiral1.centerPoint, spiral1out->centerPoint, TESTTOL) || !ptsAreSame(spiral1.endPoint, spiral1out->startPoint, TESTTOL) || !ptsAreSame(spiral1.startPoint, spiral1out->endPoint, TESTTOL))
//				failed = 1;
//
//			if (err)
//			{
//			  errorCount++;
//			  failedCount++;
//			}
//			else
//			{
//				if (failed)
//					failedCount++;
//				else
//					passedCount++;
//			}
//
//			testCaseCount++;
//		}//for jj
//
//		//test case 4 -- spiral and geo boundary (3 shapes)
//		//Use spiral2 from above
//		err = 0;
//		direct(spiral2.startPoint, spiral2.startAz + M_PI_4, (spiral2.startRadius + spiral2.endRadius)/4.0, &tempLLPoint, EPS);
//		err |= createGeo(&geo1, spiral2.startPoint, tempLLPoint, SEGMENT, EPS);
//		err |= createGeo(&geo2, spiral2.endPoint, tempLLPoint, SEGMENT, EPS);
//		bin = createBndry();
//		bout[0] = createBndry();
//		//add the shapes to the boundary
//		err |= addSpiralToBndry(&bin, &spiral2);
//		err |= addGeoToBndry(&bin, &geo1);
//		err |= addGeoToBndry(&bin, &geo2);
//		if (err) printf("\nError(s) occurred during setup: 0x%lx", err);
//
//		//order the boundary
//		err = separateBndry(bin, bout, &numberOfBoundaries, TOL, EPS);
//		if (err) printf("\nError(s) occurred in separateBndry: 0x%lx", err);
//
//		spiral2out = (Spiral*) &bout[0].elements[0];
//		geo2out = (Geodesic*) &bout[0].elements[1];
//		geo1out = (Geodesic*) &bout[0].elements[2];
//
//		failed = 0;
//		if (!ptsAreSame(spiral2.centerPoint, spiral2out->centerPoint, TESTTOL) || !ptsAreSame(spiral2.startPoint, spiral2out->startPoint, TESTTOL) || !ptsAreSame(spiral2.endPoint, spiral2out->endPoint, TESTTOL))
//			failed = 1;
//		if (!ptsAreSame(geo2out->startPoint, geo2.startPoint, TESTTOL) || !ptsAreSame(geo2out->endPoint, geo2.endPoint, TESTTOL))
//			failed = 1;
//		 //geo1 is reversed by separateBndry
//		if (!ptsAreSame(geo1out->startPoint, geo1.endPoint, TESTTOL) || !ptsAreSame(geo1out->endPoint, geo1.startPoint, TESTTOL))
//			failed = 1;
//
//		if (err)
//		{
//			errorCount++;
//			failedCount++;
//		}
//		else
//		{
//			if (failed)
//				failedCount++;
//			else
//				passedCount++;
//		}
//
//		testCaseCount++;
//
//		//test case 5 -- spiral and geo (only 2 shapes)
//		//Use spiral2 from above
//		err = 0;
//		err |= createGeo(&geo1, spiral2.startPoint, spiral2.endPoint, SEGMENT, EPS);
//		bin = createBndry();
//		bout[0] = createBndry();
//		//add the shapes to the boundary
//		err |= addSpiralToBndry(&bin, &spiral2);
//		err |= addGeoToBndry(&bin, &geo1);
//		if (err) printf("\nError(s) occurred during setup: 0x%lx", err);
//
//		//order the boundary
//		err = separateBndry(bin, bout, &numberOfBoundaries, TOL, EPS);
//		if (err) printf("\nError(s) occurred in separateBndry: 0x%lx", err);
//
//		spiral2out = (Spiral*) &bout[0].elements[0];
//		geo1out = (Geodesic*) &bout[0].elements[1];
//
//		failed = 0;
//		if (!ptsAreSame(spiral2.centerPoint, spiral2out->centerPoint, TESTTOL) || !ptsAreSame(spiral2.startPoint, spiral2out->startPoint, TESTTOL) || !ptsAreSame(spiral2.endPoint, spiral2out->endPoint, TESTTOL))
//			failed = 1;
//		 //geo1 is reversed by separateBndry
//		if (!ptsAreSame(geo1out->startPoint, geo1.endPoint, TESTTOL) || !ptsAreSame(geo1out->endPoint, geo1.startPoint, TESTTOL))
//			failed = 1;
//
//		if (err)
//		{
//			errorCount++;
//			failedCount++;
//		}
//		else
//		{
//			if (failed)
//				failedCount++;
//			else
//				passedCount++;
//		}
//
//		testCaseCount++;

    }//for ii

    set.testCases = testCaseCount;
    set.pass = passedCount;
    set.fail = failedCount;
    set.unverified = unverifiedCount;
    set.setupFailures = setupFailureCount;
    set.errors = errorCount;

    displayTestSet(set);

    printf("\nFinish testSeparateBndry_Set2\n\n\n");

    return set;
}

/*
 * NAME: testSeparateBndry_AllSets
 *
 * DESCRIPTION:
 * 		This function is used to test the boundary's separateBndry function.
 *
 * 		This function runs all the test cases and returns
 * 		the cumulative totals obtained.
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		testSeparateBndry_AllSets(TestSuite) - A test suite with the folling metrics:
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
TestSuite testSeparateBndry_AllSets()
{
	TestSuite suite;
	TestSet set1;
	TestSet set2;

    suite = newTestSuite("testSeparateBndry_AllSets");

    printf("\nStart testSeparateBndry_AllSets\n");

    set1 = testSeparateBndry_Set1();
    set2 = testSeparateBndry_Set2();
    addTestSet(set1,&suite);
    addTestSet(set2,&suite);

    displayTestSuite(suite);

    printf("\nFinish testSeparateBndry_AllSets\n\n\n");

    return suite;
}
/*
 * The following 3 routines are used to create Boundaries for testing.
*/
ErrorSet constructLocusArcBoundary(LLPoint geoStart, double geocrs, double geolen, double locusDist, Boundary *b)
{
  ErrorSet err = 0;
  LLPoint geoEnd;
  Locus loc1, loc2;
  Arc arc1, arc2;

  err |= direct(geoStart, geocrs, geolen, &geoEnd, EPS);
  err |= createLocus(&loc1, geoStart, geoEnd, locusDist, locusDist, SEGMENT, TOL, EPS);
  err |= createLocus(&loc2, geoStart, geoEnd, -locusDist, -locusDist, SEGMENT, TOL, EPS);
  err |= createArc(&arc1, geoStart, loc1.locusStart, loc2.locusStart, CLOCKWISE, TOL, EPS);
  err |= createArc(&arc2, geoEnd, loc1.locusEnd, loc2.locusEnd, COUNTERCLOCKWISE, TOL, EPS);
  err |= addLocusToBndry(b, &loc1);
  err |= addArcToBndry(b, &arc1);
  err |= addLocusToBndry(b, &loc2);
  err |= addArcToBndry(b, &arc2);
  return err;
}

ErrorSet constructGeoBoundary(LLPoint center, double radius, double crs1, double crs2, double crs3, double crs4, Boundary *b)
{
  ErrorSet err = 0;
  LLPoint geo1Start, geo2Start, geo3Start, geo4Start;
  Geodesic geo1, geo2, geo3, geo4;

  err |= direct(center, crs1, radius, &geo1Start, EPS);
  err |= direct(center, crs2, radius, &geo2Start, EPS);
  err |= direct(center, crs3, radius, &geo3Start, EPS);
  err |= direct(center, crs4, radius, &geo4Start, EPS);
  err |= createGeo(&geo1, geo1Start, geo2Start, SEGMENT, EPS);
  err |= createGeo(&geo2, geo2Start, geo3Start, SEGMENT, EPS);
  err |= createGeo(&geo3, geo3Start, geo4Start, SEGMENT, EPS);
  err |= createGeo(&geo4, geo4Start, geo1Start, SEGMENT, EPS);
  err |= addGeoToBndry(b, &geo1);
  err |= addGeoToBndry(b, &geo2);
  err |= addGeoToBndry(b, &geo3);
  err |= addGeoToBndry(b, &geo4);
  return err;
}

ErrorSet construct2Arc2LocusBoundary(LLPoint center, double innerRadius, double outerRadius, double startSubAngle, double endSubAngle, double length, Boundary *b)
{
  ErrorSet err = 0;
  double fcrs, bcrs, dist, locDist;
  Locus startLocLeft, startLocRight, endLocLeft, endLocRight;
  Arc innerArc, outerArc;
  Geodesic beginGeo, endGeo;
  LLPoint innerArcStart, innerArcEnd, outerArcStart, outerArcEnd;
  LLPoint tempLLPoint1, tempLLPoint2;

  err |= direct(center, startSubAngle, innerRadius, &innerArcStart, EPS);
  err |= direct(center, startSubAngle, outerRadius, &outerArcStart, EPS);
  err |= direct(center, endSubAngle, innerRadius, &innerArcEnd, EPS);
  err |= direct(center, endSubAngle, outerRadius, &outerArcEnd, EPS);
  err |= createArc(&innerArc, center, innerArcStart, innerArcEnd, CLOCKWISE, TOL, EPS);
  err |= createArc(&outerArc, center, outerArcStart, outerArcEnd, CLOCKWISE, TOL, EPS);
  err |= inverse(innerArcStart, outerArcStart, &fcrs, &bcrs, &dist, EPS);
  locDist = dist / 2.0;
  err |= direct(innerArcStart, fcrs, locDist, &tempLLPoint1, EPS);
  err |= invCrs(tempLLPoint1, innerArcStart, &fcrs, &bcrs, EPS);
  err |= direct(tempLLPoint1, fcrs + M_PI_2, length, &tempLLPoint2, EPS);
  err |= createLocus(&startLocLeft, tempLLPoint1, tempLLPoint2, -locDist, -locDist, SEGMENT, TOL, EPS);
  err |= createLocus(&startLocRight, tempLLPoint1, tempLLPoint2, locDist, locDist, SEGMENT, TOL, EPS);
  err |= createGeo(&beginGeo, startLocLeft.locusEnd, startLocRight.locusEnd, SEGMENT, EPS);
  err |= inverse(innerArcEnd, outerArcEnd, &fcrs, &bcrs, &dist, EPS);
  locDist = dist / 2.0;
  err |= direct(innerArcEnd, fcrs, locDist, &tempLLPoint1, EPS);
  err |= invCrs(tempLLPoint1, innerArcEnd, &fcrs, &bcrs, EPS);
  err |= direct(tempLLPoint1, fcrs - M_PI_2, length, &tempLLPoint2, EPS);
  err |= createLocus(&endLocLeft, tempLLPoint1, tempLLPoint2, -locDist, -locDist, SEGMENT, TOL, EPS);
  err |= createLocus(&endLocRight, tempLLPoint1, tempLLPoint2, locDist, locDist, SEGMENT, TOL, EPS);
  err |= createGeo(&endGeo, endLocLeft.locusEnd, endLocRight.locusEnd, SEGMENT, EPS);
  err |= addGeoToBndry(b, &beginGeo);
  err |= addLocusToBndry(b, &startLocLeft); 
  err |= addArcToBndry(b, &innerArc);
  err |= addLocusToBndry(b, &endLocRight); 
  err |= addGeoToBndry(b, &endGeo);
  err |= addLocusToBndry(b, &endLocLeft); 
  err |= addArcToBndry(b, &outerArc);
  err |= addLocusToBndry(b, &startLocRight); 
  return err;
}

//Creates a boundary composed of two loci capped by two spirals
ErrorSet createLocusSpiralBndry(LLPoint sp, double lineAz, double len, double startRad, double endRad, Boundary *b)
{
  ErrorSet err = 0;
  LLPoint ep, mc1, mc2, tempPt1, tempPt2, perpPt1, perpPt2;
  Geodesic line1, line2;
  Spiral spiral1, spiral2, mcSpiral, tempSp1, tempSp2;
  Locus outLoc1, outLoc2;
  double az12, az21, perpCrs1, perpDist1, perpCrs2, perpDist2,rad;  

  err |= direct(sp, lineAz, len, &ep, EPS);
  err |= invCrs(sp, ep, &az12, &az21, EPS);

  err |= createSpiral(&spiral1, sp, startRad, endRad, az12, az12, CLOCKWISE, EPS);
  err |= createSpiral(&spiral2, ep, startRad, endRad, az21, az21, COUNTERCLOCKWISE, EPS);

  err |= ptOnSpiral(spiral1, az12 + spiral1.dir * M_PI/4, &tempPt1, EPS);
  err |= ptOnSpiral(spiral2, az21 + spiral2.dir * M_PI/4, &tempPt2, EPS);
  err |= createGeo(&line1, tempPt1, tempPt2, INFINITE, EPS);
  err |= ptOnSpiral(spiral1, az12 - spiral1.dir * M_PI_4, &tempPt1, EPS);
  err |= ptOnSpiral(spiral2, az21 - spiral2.dir * M_PI_4, &tempPt2, EPS);
  err |= createGeo(&line2, tempPt1, tempPt2, INFINITE, EPS);

  //spiralMidChord calls determine where to truncate the spirals  // and also define the locus start and end points.
  err |= spiralRadius(spiral1, az12 + spiral1.dir * M_PI_2, &rad);
  err |= createSpiralSection(spiral1, az12 + spiral1.dir * M_PI_2, rad, &mcSpiral, EPS);
  err |= spiralMidChord(mcSpiral, line1, &mc1, TOL, EPS);

  err |= spiralRadius(spiral2, az21 + spiral2.dir * M_PI_2, &rad);
  err |= createSpiralSection(spiral2, az21 + spiral2.dir * M_PI_2, rad, &mcSpiral, EPS);
  err |= spiralMidChord(mcSpiral, line1, &mc2, TOL, EPS);

  err |= projectToGeo(sp, lineAz, mc1, &perpPt1, &perpCrs1, &perpDist1, TOL, EPS);
  err |= projectToGeo(sp, lineAz, mc2, &perpPt2, &perpCrs2, &perpDist2, TOL, EPS);

  err |= createLocus(&outLoc1, perpPt1, perpPt2, perpDist1, perpDist2, SEGMENT, TOL, EPS);

  err |= addLocusToBndry(b, &outLoc1);

  err |= moveSpiralStartToPt(spiral1, mc1, &tempSp1, TOL, EPS);
  err |= moveSpiralStartToPt(spiral2, mc2, &tempSp2, TOL, EPS);

  err |= spiralRadius(spiral1, az12 - spiral1.dir * M_PI_2, &rad);
  err |= createSpiralSection(spiral1, az12 - spiral1.dir * M_PI_2, rad, &mcSpiral, EPS);
  err |= spiralMidChord(mcSpiral, line2, &mc1, TOL, EPS);

  err |= spiralRadius(spiral2, az21 - spiral2.dir * M_PI_2, &rad);
  err |= createSpiralSection(spiral2, az21 - spiral2.dir * M_PI_2, rad, &mcSpiral, EPS);
  err |= spiralMidChord(mcSpiral, line2, &mc2, TOL, EPS);

  err |= projectToGeo(sp, lineAz, mc1, &perpPt1, &perpCrs1, &perpDist1, TOL, EPS);
  err |= projectToGeo(sp, lineAz, mc2, &perpPt2, &perpCrs2, &perpDist2, TOL, EPS);
  err |= createLocus(&outLoc2, perpPt1, perpPt2, -perpDist1, -perpDist2, SEGMENT, TOL, EPS);

  err |= addLocusToBndry(b, &outLoc2);

  err |= moveSpiralEndToPt(tempSp1, mc1, &spiral1, TOL, EPS);
  err |= moveSpiralEndToPt(tempSp2, mc2, &spiral2, TOL, EPS);

  err |= addSpiralToBndry(b, &spiral1);
  err |= addSpiralToBndry(b, &spiral2);
  return err;
}
} //namespace
