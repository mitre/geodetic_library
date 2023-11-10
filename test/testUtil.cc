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

#include <iostream>
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
#include "testUtil.h"

using namespace std;

namespace geolib_idealab {

//TODO this EPS declaration needs to go away in the long term, see isEqualLLPoint TODO
#ifndef EPS
#define EPS 0.5e-13
#endif

#ifdef TOL
#undef EPS
#endif
#define TOL 1.37e-9

#define MAX_LINE_SIZE 100 //array bounds should be set at compile time

const char* OUTPUT_DEC_PREC = "20";  //Output decimal precision
const long LATMIN = -90, LATMAX = 90; //degrees
const long LONMIN = -180, LONMAX = 180; //degrees
const long AZMIN = 0, AZMAX = 360; //degrees
const long DISTMIN = 0, DISTMAX = 5400; //nautical miles
const long MAXNUM = 3600; //non-zero modulus for generating a randon number
const long SLOPEMIN = 0, SLOPEMAX = 72;


/*
 * NAME: randLat
 *
 * DESCRIPTION:
 * 		Return a random latitude in the range [LATMIN,LATMAX]
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		randLat(double): a random latitude in degrees
 *
 */
double randLat()
{
   // To change to uniform distribution over area, uncomment these lines.  AAES-1098
   //double z = (rand() - RAND_MAX/2) / (double) RAND_MAX;
   //return asin(z) * 180 / M_PI;
    return (LATMIN + (LATMAX - LATMIN) * (double) (rand() % MAXNUM + 1)
            / (double) MAXNUM);
}

/*
 * NAME: randLon
 *
 * DESCRIPTION:
 * 		Return a random longitude in the range [LONMIN,LONMAX]
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		randLon(double): a random longitude in degrees
 *
 */
double randLon()
{
    return (LONMIN + (LONMAX - LONMIN) * (double) (rand() % MAXNUM + 1)
            / (double) MAXNUM);
}

/*
 * NAME: randAzimuth
 *
 * DESCRIPTION:
 * 		Return a random azimuth in the range [AZMIN,AZMAX]
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		randAzimuth(double): a random azimuth in degrees
 *
 */
double randAzimuth()
{
    return (AZMIN + (AZMAX - AZMIN) * (double) (rand() % MAXNUM + 1)
            / (double) MAXNUM);
}

/*
 * NAME: randDist
 *
 * DESCRIPTION:
 * 		Return a random distance in the range [DISTMIN,DISTMAX]
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		randDist(double): a random distance in nautical miles
 *
 */
double randDist()
{
    return (DISTMIN + (DISTMAX - DISTMIN) * (double) (rand() % MAXNUM + 1)
            / (double) MAXNUM);
}

/*
 * NAME: randDouble
 *
 * DESCRIPTION:
 * 		Return a random double value in the range [minDouble,maxDouble]
 *
 * INPUT(Type):
 * 		minDouble(double) = the lower bound of the desired range
 * 		maxDouble(double) = the upper bound of the desired range
 *
 * OUTPUT(Return Type):
 * 		randDouble(double): a random double
 *
 */
double randDouble(double minDouble, double maxDouble)
{
    return (minDouble + (maxDouble - minDouble)
            * (double) (rand() % MAXNUM + 1) / (double) MAXNUM);
}

/*
 * NAME: randDist
 *
 * DESCRIPTION:
 * 		Return a random signed distance
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		randSignedDist(double): a random signed distance in nautical miles
 *
 */
double randSignedDist()
{
    return (-DISTMAX + (DISTMAX + DISTMAX) * (double) (rand() % MAXNUM + 1)
            / (double) MAXNUM);
}

/*
 * NAME: randSlope
 *
 * DESCRIPTION:
 * 		Return a random slope
 *
 * INPUT(Type):
 * 		None
 *
 * OUTPUT(Return Type):
 * 		randSlope(double): a random slope
 *
 */
double randSlope()
{
    return (-SLOPEMAX + (SLOPEMAX + SLOPEMAX) * (double) (rand() % MAXNUM + 1)
            / (double) MAXNUM);
}


/*
 * NAME: testUtil.c
 *
 * DESCRIPTION:
 * 		This source file contains utility functions that are used for testing the
 * 		geolib project.  This source does not contain functions used to generate
 * 		test data, the functions to be tested, or the test functions used for
 * 		validating the geolib project.
 *
 */

/*
 * NAME: newTestSuite
 *
 * DESCRIPTION:
 * 		Constructor for a new test set.
 * 		Initialize the name of the test set
 * 		with the given name and initialize everything
 * 		else to zero.
 *
 * INPUT(Type):
 * 		name(string) = The name of the test set
 *
 * OUTPUT(Return Type):
 * 		newTestSet(TestSet) = A new instance of	TestSet.
 *
 */
TestSet newTestSet(std::string name)
{
    TestSet set = { name, 0, 0, 0, 0, 0, 0, 0 };
    return set;
}

/*
 * NAME: newTestSuite
 *
 * DESCRIPTION:
 * 		Constructor for a new testing suite.
 * 		Initialize the name of the testing suite
 * 		with the given name, initialize the array
 * 		of elements to NULL, initialize everything
 * 		else to zero.
 *
 * INPUT(Type):
 * 		name(string) = The name of the testing suite
 *
 * OUTPUT(Return Type):
 * 		newTestSuite(TestSuite) = A new instance of TestSuite.
 *
 */
TestSuite newTestSuite(std::string name)
{
    TestSuite suite = { name, 0, NULL, 0, 0, 0, 0, 0, 0, 0 };
    return suite;
}

/*
 * NAME: displayTestSet
 *
 * DESCRIPTION:
 * 		Display testing metrics for the given test set.
 *
 * INPUT(Type):
 * 		set(TestSet) = The test set to display
 *
 * OUTPUT(Return Type):
 * 		Nothing
 *
 */
void displayTestSet(TestSet set)
{
    cout << endl << "Test Function:" << set.name << endl;
    if (set.testingError)
        cout << "************ Unable to execute " << set.name << "******************" << endl;
    else
        cout << "Tests executed:" << set.testCases <<
                " Passed:" << set.pass << "Failed:" << set.fail <<
                " Setup Failures:" << set.setupFailures << " Unverified:" << set.unverified <<
                " Errors:" << set.errors << endl;
}

/*
 * NAME: displayTestSuite
 *
 * DESCRIPTION:
 * 		Display testing metrics for the given testing suite.  Note that
 *		the testing metrics for all the elements (test suite or test set)
 * 		in the suite will be displayed as well.  Summaries for the given
 * 		suite along with all children suite are displayed as well.
 *
 * INPUT(Type):
 * 		suite(TestSuite) = The testing suite to display
 *
 * OUTPUT(Return Type):
 * 		Nothing
 *
 */
void displayTestSuite(TestSuite suite)
{
    int n = suite.length;
    int i = 0;
    TestElement element;
    TestType testType;
    TestSet* mySet;
    TestSuite* mySuite;

    cout << "Test Suite Summary:" << suite.name << endl;

    if (suite.testingSuiteError)
       cout << "************ A problem was encountered with " << suite.name << " ******************" << endl;

    for (i = 0; i < n; i++)
    {
        element = suite.elements[i];
        testType = suite.elements[i].type;

        switch (testType) {
        case SET:
            mySet = (TestSet*) element.ptr;
            displayTestSet(*mySet);
            break;
        case SUITE:
            mySuite = (TestSuite*) element.ptr;
            printf(
                    "***************************************************************************************************************\n");
            displayTestSuite(*mySuite);
            printf(
                    "*****************************************************************************************************************\n");
            break;
        default:
            printf("\nUnexpected Test Element found.");
        }
    }

    if (!suite.testingSuiteError)
        printf(
                "\nTotal Tests executed:%i Total Passed:%i Total Failed:%i Total Setup Failures:%i Total Unverified:%i Total Errors:%i\n",
                suite.testCases, suite.pass, suite.fail, suite.setupFailures,
                suite.unverified, suite.errors);
}

/*
 * NAME: outputTestSuiteCSV
 *
 * DESCRIPTION:
 *
 * 		print the given test suite metrics to the given file in
 * 		a comma separated format (CSV).
 *
 * INPUT(Type):
 * 		suite(TestSuite) = The testing suite to output to file
 * 		file(*FILE) = The output file to write to
 *
 * OUTPUT(Return Type):
 * 		Nothing
 *
 */
void outputTestSuiteCSV(TestSuite suite, FILE* file)
{
    int n = suite.length;
    int i = 0;
    TestElement element;
    TestType testType;
    TestSet* mySet;
    TestSuite* mySuite;

    if (file == NULL)
    {
    	printf("\n********** No output file found ********************\n");
    	return;
    }

    fprintf(file,"%s,%i,%i,%i,%i,%i,%i,%i\n",
    		suite.name.c_str(),suite.testCases,suite.pass,suite.fail,
    		suite.setupFailures, suite.unverified, suite.errors,
    		suite.testingSuiteError);

    for (i = 0; i < n; i++)
    {
        element = suite.elements[i];
        testType = suite.elements[i].type;

        fprintf(file,"%s,",suite.name.c_str());

        switch (testType) {
        case SET:
            mySet = (TestSet*) element.ptr;
            fprintf(file,"%s,%i,%i,%i,%i,%i,%i,%i\n",
            		mySet->name.c_str(),mySet->testCases,mySet->pass,mySet->fail,
            		mySet->setupFailures, mySet->unverified, mySet->errors,
            		mySet->testingError);
            break;
        case SUITE:
            mySuite = (TestSuite*) element.ptr;
            outputTestSuiteCSV(*mySuite,file);
            break;
        default:
            fprintf(file,"%s\n","Unexpected Test Element found");
        }
    }
}


/*
 * NAME: outputTestSuiteXML
 *
 * DESCRIPTION:
 *
 *       print the given test suite metrics to the given file in
 *       a comma separated format (CSV).
 *
 * INPUT(Type):
 *       suite(TestSuite) = The testing suite to output to file
 *       level(int) = Level in heirarchy (0=testsuites, 1=testsuite, other=process children only)
 *       file(*FILE) = The output file to write to
 *
 * OUTPUT(Return Type):
 *       Nothing
 *
 */
void outputTestSuiteXML(const TestSuite &suite, const int level, FILE* file)
{
    const int n = suite.length;
    int i = 0;
    const int suiteErrors(suite.errors + suite.setupFailures + suite.testingSuiteError + suite.unverified);
    TestElement element;
    TestType testType;
    TestSet* mySet;
    TestSuite* mySuite;

    if (file == NULL)
    {
      printf("\n********** No output file found ********************\n");
      return;
    }

    string closingTag("");
    switch (level) {
    case 0:
       // top level goes to testsuites element
       fprintf(file, "<testsuites name=\"%s\" tests=\"%d\" failures=\"%d\" errors=\"%d\">\n",
             suite.name.c_str(),suite.testCases,suite.fail,suiteErrors);
       closingTag = "</testsuites>\n";
       break;
    case 1:
       // this level goes to a testsuite element
       fprintf(file, "  <testsuite name=\"%s\" tests=\"%d\" failures=\"%d\" errors=\"%d\">\n",
             suite.name.c_str(),suite.testCases,suite.fail,suiteErrors);
       closingTag = "  </testsuite>\n";
       break;
    default:
       // Flatten this level:  no tag and no closing tag
       // Children will be recursed to get all the TestSet members
       closingTag = "";
    }

    // handle each child
    for (i = 0; i < n; i++)
    {
        element = suite.elements[i];
        testType = suite.elements[i].type;

        switch (testType) {
        case SET:
            mySet = (TestSet*) element.ptr;
            if (level == 0) {
               cerr << "WARNING:  Found TestSet at first level -- " << mySet->name << endl;
            }
            else {
               fprintf(file,"    <testcase name=\"%s\" assertions=\"%d\">\n",
                  mySet->name.c_str(),mySet->testCases);
               // TODO Save failure and error messages in TestSet object.
               // For now, all we have are counts.
               // Report straight-up test failures in a failure element
               if (mySet->fail > 0) {
                  fprintf(file, "      <failure message=\"%d failures\"> </failure>\n", mySet->fail);
               }
               // use error elements for all 4 kinds of errors
               if (mySet->errors > 0) {
                  fprintf(file, "      <error message=\"%d errors\"> </error>\n", mySet->errors);
               }
               if (mySet->setupFailures > 0) {
                  fprintf(file, "      <error message=\"%d setup failures\"> </error>\n", mySet->setupFailures);
               }
               if (mySet->testingError > 0) {
                  fprintf(file, "      <error message=\"%d testing errors\"> </error>\n", mySet->testingError);
               }
               if (mySet->unverified > 0) {
                  fprintf(file, "      <error message=\"%d unverified\"> </error>\n", mySet->unverified);
               }
               fprintf(file,"    </testcase>\n");
            }
            break;
        case SUITE:
            mySuite = (TestSuite*) element.ptr;
            outputTestSuiteXML(*mySuite, level+1, file);
            break;
        default:
            fprintf(file,"%s\n","Unexpected Test Element found");
        }
    }
    fputs(closingTag.c_str(), file);
}


/*
 * NAME: outputTestSuiteMetrics
 *
 * DESCRIPTION:
 *
 * 		output the given test suite metrics to the given file in
 * 		a comma separated format (CSV).
 *
 * INPUT(Type):
 * 		suite(TestSuite) = The testing suite to display
 * 		outputFile(FILE*) = The output file path to write to.  If NULL,
 * 			then create an output file to in directory geolib/test.
 *
 * OUTPUT(Return Type):
 * 		Nothing
 *
 */
void outputTestSuiteMetrics(TestSuite suite, char* outputFile)
{
    const char* fileLocation = FILEROOT "/test/";
    char fileName[100];
    char timeStamp[100];
    char completeFilePath[500];
    time_t now;
    FILE* file;

    if (outputFile == NULL)
    {
/*		now = time(NULL);
		sprintf(fileName,"%s%s",suite.name.c_str(),"Metrics_");
		strftime(timeStamp,100,"%m%d%Y_%H%M%S%p",localtime(&now));
		sprintf(completeFilePath,"%s%s%s%s",fileLocation,fileName,timeStamp,".csv");*/
      sprintf(completeFilePath, "geolib-test-output.xml");
    }
    else
    {
    	sprintf(completeFilePath,"%s",outputFile);
    }

    //reportConstants(completeFilePath);
    //file = fopen(completeFilePath,"a");
    file = fopen(completeFilePath,"w");

	if (!file)
	{
		printf("\nCould not open output file:%s\n", completeFilePath);
		return;
	}

	printf("\nOutput file location:%s", completeFilePath);
    printf("\nOutputting Test Suite Summary:%s\n", suite.name.c_str());

    //fprintf(file,"%s","parent,child,testCaseCount,passCount,failCount,setupFailureCount,unverifiedCount,errorCount,testElementError\n");
    //fprintf(file,"%s,","root");
    fprintf(file, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");

    outputTestSuiteXML(suite, 0, file);

    fclose(file);
}

/*
 * NAME: growSuite
 *
 * DESCRIPTION:
 * 		Dynamically grow the array of TestElements for the given
 * 		testing suite by one element.
 *
 * INPUT(Type):
 * 		suite(TestSuite*) = The testing suite to grow
 *
 * OUTPUT(Return Type):
 * 		Nothing
 *
 */
void growSuite(TestSuite* suite)
{
    int n = suite->length;
    int i = 0;
    TestElement *newArray, *oldArray;

    n++;
    newArray = new TestElement[n];
    (newArray[n - 1]).ptr = NULL;
    (newArray[n - 1]).type = SET;

    for (i = 0; i < n - 1; i++)
    {
        /* Copy existing data pointers */
        newArray[i] = suite->elements[i];
    }

    oldArray = suite->elements;
    suite->elements = newArray;
    suite->length = n;
    delete[] oldArray;

}

/*
 * NAME: addTestSet
 *
 * DESCRIPTION:
 * 		Add the test set to the given testing suite
 *
 * INPUT(Type):
 * 		set(TestSet) = The test set being added to the test suite
 * 		suite(TestSuite*) = The testing suite
 *
 * OUTPUT(Return Type):
 * 		Nothing
 *
 */
void addTestSet(TestSet set, TestSuite* suite)
{
    int n;
    TestSet* newSet = NULL;

    growSuite(suite);
    n = suite->length;

    newSet = new TestSet();
    *newSet = set;

    suite->elements[n - 1].ptr = (TestSet*) newSet;
    suite->elements[n - 1].type = SET;

    suite->pass += set.pass;
    suite->fail += set.fail;
    suite->errors += set.errors;
    suite->unverified += set.unverified;
    suite->setupFailures += set.setupFailures;
    suite->testCases += set.testCases;
    suite->testingSuiteError |= set.testingError;

}

/*
 * NAME: addTestSuite
 *
 * DESCRIPTION:
 * 		Add the child testing suite to the parent testing suite
 *
 * INPUT(Type):
 * 		childSuite(TestSuite) = The suite to be added to the parent suite
 * 		parentSuite(TestSuite*) = The parent testing suite
 *
 * OUTPUT(Return Type):
 * 		Nothing
 *
 */
void addTestSuite(TestSuite childSuite, TestSuite* parentSuite)
{
    int n;
    TestSuite* newSuite = NULL;

    growSuite(parentSuite);
    n = parentSuite->length;

    newSuite = new TestSuite();
    *newSuite = childSuite;

    parentSuite->elements[n - 1].ptr = (TestSuite*) newSuite;
    parentSuite->elements[n - 1].type = SUITE;

    parentSuite->pass += childSuite.pass;
    parentSuite->fail += childSuite.fail;
    parentSuite->errors += childSuite.errors;
    parentSuite->unverified += childSuite.unverified;
    parentSuite->setupFailures += childSuite.setupFailures;
    parentSuite->testCases += childSuite.testCases;
    parentSuite->testingSuiteError |= childSuite.testingSuiteError;

}

/*
 * NAME: checkForNullAndDuplicatePointers
 *
 * DESCRIPTION:
 * 		Extend the checkForDuplicatePointers function to check a list of pointers
 * 		for duplicates and/or NULL pointers.
 *
 * INPUT(Type):
 * 		firstP(void*) = The first pointer to be verified
 * 		...(void*) = A list of additional pointers
 *
 * OUTPUT(Return Type):
 * 		checkForNullAndDuplicatePointers(int) =
 * 			Return 0 for no NULL pointers, and no duplicate pointers
 * 			Return 1 for NULL pointer found
 * 			Return 2 for duplicate pointers found
 * 			Return 3 for NULL and duplicate pointers found
 *
 */
int checkForNullAndDuplicatePointers(void* firstP, ...)
{
    const long NUMBERFORKNOWNPOINTER = 2008229;
    const void* KNOWNPOINTER = &NUMBERFORKNOWNPOINTER;

    va_list argp;
    int nullAndDuplicatePointersFound = 0;

    void* p;
    void* q[20];
    int i = 0, n = 0, j = 0, k = 0;
    int nullPointerFound = 0;
    int duplicatePointersFound = 0;
    int numberOfNotNULLPointers = 0;
    int numberOfPointers = 0;

    va_start(argp, firstP);
    p = firstP;
    n++;

    //Check for NULL
    while (p != KNOWNPOINTER)
    {
        if (p == NULL)
        {
            nullPointerFound = 1;
        }
        else
        {
            q[i] = p;
            i++;
        }
        p = va_arg(argp, void*);
        n++;
    }
    numberOfNotNULLPointers = i;
    numberOfPointers = n;

    //Check for non-NULL duplicates
    if (numberOfNotNULLPointers > 1)
    {
        for (j = 0; j < numberOfNotNULLPointers - 1; j++)
        {
            for (k = j + 1; k < numberOfNotNULLPointers; k++)
            {
                if (q[j] == q[k])
                {
                    duplicatePointersFound = 1;
                }
            }
        }
    }

    nullAndDuplicatePointersFound = nullPointerFound + 2
            * duplicatePointersFound;
    //	printf("nPtr=%d  nNotNullPtr=%d  nullPtrFound=%d  dupPtrFound=%d\n",
    //			numberOfPointers,numberOfNotNULLPointers,
    //			nullPointerFound,duplicatePointersFound);

    va_end(argp); //required to return normally
    return nullAndDuplicatePointersFound;
}

/*
 * NAME: isEqualLLPoint
 *
 * DESCRIPTION:
 * 		Compare the distance between two points.
 *      If the distance is less than 0.03 centimeters then the two points
 *		are considered to be equal.
 *
 * INPUT(Type):
 * 		pt1(LLPoint) = The first point
 * 		pt2(LLPoint) = The second point
 *
 * OUTPUT(Return Type):
 * 		isEqualLLPoint(int) =
 * 			Return 0 if the two points are not "equal"
 * 			Return 1 if the two points are "equal"
 *
 */
//TODO refactor and redesign this method so that EPS and NMTOL are passed in arguments
//TODO refactor to handle error codes properly; use masks and this should probably return the error code created by invDist
int isEqualLLPoint(LLPoint pt1, LLPoint pt2)
{
    ErrorSet err = 0;
    double dist12;
    double NMTOL = 0.03 / 100.0 / 1852.0; //0.03 cm or ~1.62e-7 nm

    err = invDist(pt1, pt2, &dist12, EPS);
    //printf("distance: %18.12f\n", fabs(dist12));
    if ((err == 0) && (fabs(dist12) < NMTOL))
    {
        return 1; //true
    }
    else
    {
        return 0; //false
    }
}

/*
 * NAME: convertStringToShapeType
 *
 * DESCRIPTION:
 * 		Convert a string to an object of type ShapeType.  First
 * 		checks the input string against the null delimter value
 * 		to determine if the two stings match.  If so then the
 * 		default value is returned.  If an invalid string is entered
 * 		then the default value is returned and a C error is raised.
 *
 * INPUT(Type):
 * 		str(char*) = The string value that is being converted.
 * 					 Valid string values: arc, geodesic, locus, spiral, llpoint
 * 		nullDelim(char*) = The string value that denotes a C NULL pointer.
 * 					 The input string value is first checked against this string.
 * 		dflt(int) =  The default type to be returned when an error occurs or
 * 					 if the input string matches the null delimiter string.
 * 					 Valid values:
 * 					 ARC = 0, GEODESIC = 1, LOCUS = 2, SPIRAL = 3, LLPOINT = 4
 *
 * OUTPUT(Return Type):
 * 		convertStringToShapeType(ShapeType) =
 * 			The ShapeType that corresponds to the input string.
 *
 */
//TODO Refactor so that the arguments names are more descriptive
//TODO Add valid input checks
ShapeType convertStringToShapeType(char* str, char* nullDelim, ShapeType dflt)
{
    ShapeType value = dflt;

    if (NULL != str)
    {
        if (!strcasecmp(str, nullDelim))
        {
            value = dflt;
        }
        else if (!strcasecmp(str, "ARC"))
        {
            value = ARC;
        }
        else if (!strcasecmp(str, "GEODESIC"))
        {
            value = GEODESIC;
        }
        else if (!strcasecmp(str, "LOCUS"))
        {
            value = LOCUS;
        }
        else if (!strcasecmp(str, "SPIRAL"))
        {
            value = SPIRAL;
        }
        else if (!strcasecmp(str, "LLPOINT"))
        {
            value = LLPOINT;
        }
        else
        {
            perror("convertStringToShapeType Error:  Unexpected shape type");
            value = dflt;
        }
    }
    return value;
}

/*
 * NAME: convertStringToShapeType
 *
 * DESCRIPTION:
 * 		Convert a string to an object of type LineType.  First
 * 		checks the input string against the null delimter value
 * 		to determine if the two stings match.  If so then the
 * 		default value is returned.  If an invalid string is entered
 * 		then the default value is returned and a C error is raised.
 *
 * INPUT(Type):
 * 		str(char*) = The string value that is being converted.
 * 					 Valid string values: segment, semiinfinite, infinite
 * 		nullDelim(char*) = The string value that denotes a C NULL pointer.
 * 					 The input string value is first checked against this string.
 * 		dflt(int) =  The default type to be returned when an error occurs or
 * 					 if the input string matches the null delimiter string.
 * 					 Valid values: 	SEGMENT = 0, SEMIINFINITE = 1, INFINITE = 2
 *
 * OUTPUT(Return Type):
 * 		convertStringToLineType(LineType) =
 * 			The LineType that corresponds to the input string.
 *
 */
//TODO Figure out why strcasecmp is causing build errors with non-gcc compilers, would help clean up this function if used
LineType convertStringToLineType(char* str, char* nullDelim, LineType dflt)
{
    LineType value;

    if (!strcmp(str, nullDelim))
    {
        value = dflt;
    }
    else if (!strcasecmp(str, "segment"))
    {
        value = SEGMENT;
    }
    else if (!strcasecmp(str, "semiinfinite"))
    {
        value = SEMIINFINITE;
    }
    else if (!strcasecmp(str, "infinite"))
    {
        value = INFINITE;
    }
    else
    {
        perror("convertStringToLineType Error:  Unexpected line type");
        value = dflt;
    }

    return value;
}

int quadrant(LLPoint pt, LLPoint ref)
{
    int quad = 0;
    if (pt.latitude > ref.latitude)
    {
        if (pt.longitude > ref.longitude)
            quad = 1;
        else
            quad = 2;
    }
    else
    {
        if (pt.longitude > ref.longitude)
            quad = 4;
        else
            quad = 3;
    }
    return quad;
}

//TODO This function is a temporary duplicate of the displayPt function.
//Once the refactoring of libWGS84.c is complete and displayPt is public
//all references to printLLPoint should be replaced with displayPt.
void printLLPoint(string ptName, LLPoint pt)
{
    double DEG2RAD = M_PI / 180.0;

    printf("%s:  LAT=%.20lf  LON=%.20lf\n", ptName.c_str(), pt.latitude / DEG2RAD,
            pt.longitude / DEG2RAD);
}

void newflush(FILE *in)
{
    int ch;

    do
        ch = fgetc(in); while (ch != EOF && ch != '\n');

    clearerr(in);
}

void newpause(void)
{
    printf("Press [Enter] to continue . . .");
    fflush(stdout);
    getchar();
}

//TODO This function is a temporary duplicate of the displayLocus function.
//Once the refactoring of libWGS84.c is complete and displayLocus is public
//all references to printLocus should be replaced with displayLocus.
void printLocus(const char* locusName, Locus locus)
{
    double DEG2RAD = M_PI / 180.0;
    double locslope = atan((locus.endDist - locus.startDist) / locus.geoLength);
    printf("%s: geoLength=%16.10f    geoAz=%16.10f     geoRevAz=%16.10f\n",
            locusName, locus.geoLength, locus.geoAz / DEG2RAD, locus.geoRevAz
                    / DEG2RAD);
    printf(
            "%s: startDist=%16.10f  endDist=%16.10f  locSlopeDeg=%16.10f  lineType=%d\n",
            locusName, locus.startDist, locus.endDist, locslope / DEG2RAD,
            locus.lineType);
    printLLPoint("   geoStart  ", locus.geoStart);
    printLLPoint("   geoEnd    ", locus.geoEnd);
    printLLPoint("   locusStart", locus.locusStart);
    printLLPoint("   locusEnd  ", locus.locusEnd);
}

double crsZeroToTwoPI(double crs)
{
    double crsWithinRange;
    crsWithinRange = crs;
    if (crs < 0.0)
        crsWithinRange = crs + M_2PI;
    if (crs >= M_2PI)
        crsWithinRange = crs - M_2PI;
    return crsWithinRange;
}

//TODO refactor and redesign this method so that EPS and TOL are passed in arguments
//TODO refactor to handle error codes properly; use masks and this should probably return the error code(s) generated
void createLocusThroughPoint(Locus *locusPtr, LLPoint pt, double crsLoc,
                             double slopeLocAngle, double d1, double d2,
                             double d3)
{
    /*
     This is a helper function to create a locus that passes through a given point
     at the specified azimuth (crsLoc).  Memory is allocated by calling function.

     a.	Compute the azimuth from pt to the corresponding point on the defining geodesic
     using course of locus at pt and the slope of locus
     b.	Compute the point on the defining geodesic using the azimuth and distance d1
     c.	Compute azimuth to pt from point on geodesic
     d.	Compute two azimuths 90 degress on either side of azimuth to point and
     using the two distances (d2, d3) locate the start and end points of the
     defining geodesic
     e.	Compute startDist and endDist of locus using slope and d1
     f.	Create locus using start and end points of defining geodesic, startDist and endDist
     */
    LLPoint geoPt, geoStart, geoEnd; //points on defining geodesic
    double startDist, endDist; //distances that are used to define the locus
    double crsPtGeo; //, crsGeoPt;       //courses from pt to geoPt and reverse
    double crsGeoStart, crsGeoEnd; //courses to startPt and endPt from geoPt
    double fcrs, bcrs;
    double slopeLoc = tan(slopeLocAngle);
    LineType lineType;
    ErrorSet err = 0;

    //Find the point on the defining geodesic
    crsPtGeo = crsZeroToTwoPI(crsLoc - slopeLocAngle - M_PI_2);
    err |= direct(pt, crsPtGeo, d1, &geoPt, EPS);
//    printLLPoint("   point on defining geodesic", geoPt);

    //Compute azimuth to the point on locus from point on geodesic
    err |= invCrs(geoPt, pt, &fcrs, &bcrs, EPS);

    //Compute two azimuths 90 degress on either side of azimuth to point and
    //using the two distances (d2, d3) locate the start and end points of the
    //defining geodesic
    crsGeoStart = crsZeroToTwoPI(fcrs + M_PI_2);
    crsGeoEnd = crsZeroToTwoPI(fcrs - M_PI_2);
    err |= direct(geoPt, crsGeoStart, d2, &geoStart, EPS);
    err |= direct(geoPt, crsGeoEnd, d3, &geoEnd, EPS);
    //    printLLPoint("   geoStart  ",geoStart);
    //    printLLPoint("   geoEnd    ",geoEnd);

    //Compute startDist and endDist using d1, d2, d3 and slope of locus
    startDist = d1 - d2 * slopeLoc;
    endDist = d1 + d3 * slopeLoc;
    lineType = INFINITE;
    //    printf("   startDist=%14.8f endDist=%14.8f\n",startDist,endDist);

    //Create locus using geoStart, geoEnd, startDist, endDist and lineType
    //Locus* loc, LLPoint geoStart, LLPoint geoEnd, double startDist,
    //double endDist, LineType lineType, double eps)
    err |= createLocus(locusPtr, geoStart, geoEnd, startDist, endDist,
            lineType, TOL, EPS);
    //    printLLPoint("   locusStart",locusPtr->locusStart);
    //    printLLPoint("   locusEnd  ",locusPtr->locusEnd);
    return;
}

//TODO refactor and redesign this method so that LLDEGTOL is a passed in argument
int isEqualLLPointWithinBox(LLPoint pt1, LLPoint pt2)
{
    double DEG2RAD = M_PI / 180.0;
    double LLDEGTOL = 5.0e-8; //5e-8 degrees or 0.00018 sec

    if ((fabs(pt1.latitude - pt2.latitude) < LLDEGTOL * DEG2RAD) && ((fabs(
            pt1.longitude - pt2.longitude) < LLDEGTOL * DEG2RAD) || (fabs(360.0
            * DEG2RAD - fabs(pt1.longitude - pt2.longitude)) < LLDEGTOL
            * DEG2RAD)))
    {
        return 1; //true
    }
    else
    {
        return 0; //false
    }
}

//TODO refactor and redesign this method so that LLDEGTOL is a passed in argument
int isEqualLatitude(double lat1, double lat2)
//latitudes are in radians

{
    double DEG2RAD = M_PI / 180.0;
    double LLDEGTOL = 5.0e-8; //5e-8 degrees or 0.00018 sec

    if ((fabs(lat1 - lat2) < LLDEGTOL * DEG2RAD))
    {
        return 1; //true
    }
    else
    {
        return 0; //false
    }
}

//TODO refactor and redesign this method so that LLDEGTOL is a passed in argument
int isEqualLongitude(double lon1, double lon2)
//longitudes are in radians

{
    double DEG2RAD = M_PI / 180.0;
    double LLDEGTOL = 5.0e-8; //5e-8 degrees or 0.00018 sec

    if ((fabs(lon1 - lon2) < LLDEGTOL * DEG2RAD) || (fabs(360.0 * DEG2RAD
            - fabs(lon1 - lon2)) < LLDEGTOL * DEG2RAD))
    {
        return 1; //true
    }
    else
    {
        return 0; //false
    }
}

//TODO This looks to be a modified version of the displayPt function, strictly
//used by testing at the moment.
//Once the refactoring of libWGS84.c is complete and displayPt is public,
//this function should be analyzed and refactored to see if it is redundant.
//If not, then it should be refactored to be more intuitive and moved into the LLPoint source file.
void display(LLPoint p)
{

    int d, m;
    double s;

    rad2dms(p.latitude, &d, &m, &s);
    printf("latitude: %20.15f = ", p.latitude * 180.0 / M_PI);
    if (d < 0.0 || m < 0.0 || s < 0.0)
    {
        printf("-%02d:%02d:%08.5f; ", abs(d), abs(m), fabs(s));
    }
    else
    {
        printf(" %02d:%02d:%08.5f; ", d, m, s);
    }
    rad2dms(p.longitude, &d, &m, &s);
    printf("longitude: %20.15f = ", p.longitude * 180.0 / M_PI);
    if (d < 0.0 || m < 0.0 || s < 0.0)
    {
        printf("-%03d:%02d:%08.5f\n", abs(d), abs(m), fabs(s));
    }
    else
    {
        printf(" %03d:%02d:%08.5f\n", d, m, s);
    }

}

void displayOld(LLPoint p)
{

    int d, m;
    double s;

    rad2dms(p.latitude, &d, &m, &s);
    printf("latitude: %20.15f = %03d:%02d:%8.5f;  ", p.latitude * 180.0 / M_PI,
            d, abs(m), fabs(s));
    rad2dms(p.longitude, &d, &m, &s);
    printf("longitude: %20.15f = %04d:%02d:%8.5f\n",
            p.longitude * 180.0 / M_PI, d, abs(m), fabs(s));

}

void displayOld2(LLPoint p)
{

    int d, m;
    double s;

    crsintrad2dms(p.latitude, &d, &m, &s);
    printf("latitude: %20.15f = %d d, %d m, %20.15f s; ", p.latitude * 180.0
            / M_PI, d, m, s);
    crsintrad2dms(p.longitude, &d, &m, &s);
    printf("longitude: %20.15f = %d d, %d m, %20.15f s\n", p.longitude * 180.0
            / M_PI, d, m, s);

}

void displayOld3(LLPoint p)
{

    int d, m;
    double s;

    crsintrad2dms(p.latitude, &d, &m, &s);
    printf("latitude: %20.15f = %03d:%02d:%8.5f;  ", p.latitude * 180.0 / M_PI,
            d, abs(m), fabs(s));
    crsintrad2dms(p.longitude, &d, &m, &s);
    printf("longitude: %20.15f = %04d:%02d:%8.5f\n",
            p.longitude * 180.0 / M_PI, d, abs(m), fabs(s));

}

void parseInputFindArcLength(char* str, LLPoint* center, double* radius,
                             double* startAz, double* endAz, int* dir)
{
    int d, m;
    double s;
    char hemi[2];

    char* tok;
    //    char* orgLatStr;

    /* Read degree minutes, seconds from file into variables */
    /* Origin latitude data */
    tok = strtok(str, " \t");
    //    printf("%s\n",tok);

    if (!center)
        center = (LLPoint*) malloc(sizeof(LLPoint));

    //    orgLatStr = strdup(tok);
    if (sscanf(tok, "%d:%d:%lf%s", &d, &m, &s, hemi))
    {
        //        printf("%s\n", tok);
        //        printf("%d -- %d -- %f -- %s\n", d,m,s,hemi);
        center->latitude = dms2radHemi(d, m, s, hemi);

    }
    else
    {
        printf("String %s could not be parsed\n", tok);
        exit(1);
    }

    tok = strtok(NULL, " \t");
    //    orgLonStr = strdup(tok);
    if (sscanf(tok, "%d:%d:%lf%s", &d, &m, &s, hemi))
    {
        //        printf("%s\n", tok);
        //       printf("%d -- %d -- %f -- %s\n", d,m,s,hemi);
        center->longitude = dms2radHemi(d, m, s, hemi);

    }
    else
    {
        printf("String %s could not be parsed\n", tok);
        exit(1);
    }

    tok = strtok(NULL, " \t");
    if (sscanf(tok, "%lf", radius))
    {
        //        printf("%f\n",*radius);
    }
    else
    {
        printf("String %s could not be parsed\n", tok);
        exit(1);
    }

    tok = strtok(NULL, " \t");
    if (sscanf(tok, "%lf", startAz))
    {
        //        printf("%f\n",*startAz);
        *startAz *= M_PI / 180.0;
    }
    else
    {
        //        printf("String %s could not be parsed\n",tok);
        exit(1);
    }

    tok = strtok(NULL, " \t");
    if (sscanf(tok, "%lf", endAz))
    {
        //       printf("%f\n",*endAz);
        *endAz *= M_PI / 180.0;
    }
    else
    {
        printf("String %s could not be parsed\n", tok);
        exit(1);
    }

    tok = strtok(NULL, " \t");
    if (sscanf(tok, "%d", dir))
    {
        //       printf("%d\n",*dir);
    }
    else
    {
        printf("String %s could not be parsed\n", tok);
        exit(1);
    }

}

void parseInputDiscretizedArcLength(char* str, LLPoint* center, double* radius,
                                    double* startAz, double* endAz, int* dir)
{
    int d, m;
    double s;
    char hemi[2];

    char* tok;
    //    char* orgLatStr;

    /* Read degree minutes, seconds from file into variables */
    /* Origin latitude data */
    tok = strtok(str, " \t");
    //    printf("%s\n",tok);

    if (!center)
        center = (LLPoint*) malloc(sizeof(LLPoint));

    //    orgLatStr = strdup(tok);
    if (sscanf(tok, "%d:%d:%lf%s", &d, &m, &s, hemi))
    {
        //        printf("%s\n", tok);
        //        printf("%d -- %d -- %f -- %s\n", d,m,s,hemi);
        center->latitude = dms2radHemi(d, m, s, hemi);

    }
    else
    {
        printf("String %s could not be parsed\n", tok);
        exit(1);
    }

    tok = strtok(NULL, " \t");
    //    orgLonStr = strdup(tok);
    if (sscanf(tok, "%d:%d:%lf%s", &d, &m, &s, hemi))
    {
        //        printf("%s\n", tok);
        //       printf("%d -- %d -- %f -- %s\n", d,m,s,hemi);
        center->longitude = dms2radHemi(d, m, s, hemi);

    }
    else
    {
        printf("String %s could not be parsed\n", tok);
        exit(1);
    }

    tok = strtok(NULL, " \t");
    if (sscanf(tok, "%lf", radius))
    {
        //        printf("%f\n",*radius);
    }
    else
    {
        printf("String %s could not be parsed\n", tok);
        exit(1);
    }

    tok = strtok(NULL, " \t");
    if (sscanf(tok, "%lf", startAz))
    {
        //        printf("%f\n",*startAz);
        *startAz *= M_PI / 180.0;
    }
    else
    {
        //        printf("String %s could not be parsed\n",tok);
        exit(1);
    }

    tok = strtok(NULL, " \t");
    if (sscanf(tok, "%lf", endAz))
    {
        //       printf("%f\n",*endAz);
        *endAz *= M_PI / 180.0;
    }
    else
    {
        printf("String %s could not be parsed\n", tok);
        exit(1);
    }

    tok = strtok(NULL, " \t");
    if (sscanf(tok, "%d", dir))
    {
        //       printf("%d\n",*dir);
    }
    else
    {
        printf("String %s could not be parsed\n", tok);
        exit(1);
    }

}

double parseDMSFormat(const char* text)
{
    int d, m;
    double s;
    double angle;
    char hemi[2];

    if (sscanf(text, "%d:%d:%lf%s", &d, &m, &s, hemi))
    {
        //        printf("%s\n", tok);
        //        printf("%d -- %d -- %f -- %s\n", d,m,s,hemi);
        angle = dms2radHemi(d, m, s, hemi);

    }
    else
    {
        printf("String %s could not be parsed\n", text);
        exit(1);
    }

    return angle;

}

void parseInputInterceptAtAngle(char* str, LLPoint* pt1, double* crs,
                                LLPoint* pt2, double* angle)
{

    char* tok;

    /* First Arc's Parameters */

    /* Read degree minutes, seconds from file into variables */
    /* Center1 Origin latitude data */
    tok = strtok(str, " \t");
    //    printf("%s\n",tok);
    pt1->latitude = parseDMSFormat(tok);

    tok = strtok(NULL, " \t");
    //    orgLonStr = strdup(tok);
    pt1->longitude = parseDMSFormat(tok);

    tok = strtok(NULL, " \t");
    if (!sscanf(tok, "%lf", crs))
    {
        printf("String %s could not be parsed\n", tok);
        exit(1);
    }
    *crs *= M_PI / 180.0;

    /* Second Arc's Parameters */
    tok = strtok(NULL, " \t");
    //    printf("%s\n",tok);
    pt2->latitude = parseDMSFormat(tok);

    tok = strtok(NULL, " \t");
    //    orgLonStr = strdup(tok);
    //   printf("%s\n",tok);
    pt2->longitude = parseDMSFormat(tok);

    tok = strtok(NULL, " \t");
    if (!sscanf(tok, "%lf", angle))
    {
        printf("String for dir2 '%s' could not be parsed\n", tok);
        exit(1);
    }
    *angle *= M_PI / 180.0;

}

/** @function displayBoundaryArc
 *
 * Print the given arc properties to the console in the desired units
 *
 * @param a : The Geolib Arc to be displayed
 * @param displayRadians: 0 = degrees, otherwise use radians
 *
 */
void displayBoundaryArc(Arc a, int displayRadians)
{
	char* output;

	output = createArcString(a, NULL, SYSTEM_OUT, displayRadians);
	printf("%s", output);

}

/** @function displayBoundaryLocus
 *
 * Print the given locus properties to the console in the desired units
 *
 * @param l : The Geolib Locus to be displayed
 * @param displayRadians: 0 = degrees, otherwise use radians
 *
 */
void displayBoundaryLocus(Locus l, int displayRadians)
{
	char* output;

	output = createLocusString(l, NULL, SYSTEM_OUT, displayRadians);
	printf("%s", output);

}

/** @function displayBoundaryGeodesic
 *
 * Print the given geodesic properties to the console in the desired units
 *
 * @param g : The Geolib Geodesic to be displayed
 * @param displayRadians: 0 = degrees, otherwise use radians
 *
 */
void displayBoundaryGeodesic(Geodesic g, int displayRadians)
{
	char* output;

	output = createGeoString(g, NULL, SYSTEM_OUT, displayRadians);
	printf("%s", output);
}

/** @function displayBoundaryLLPoint
 *
 * Print the given llpoint properties to the console in the desired units
 *
 * @param p : The Geolib LLPoint to be displayed
 * @param displayRadians: 0 = degrees, otherwise use radians
 *
 */
void displayBoundaryLLPoint(LLPoint p, int displayRadians)
{
	char* output;

	output = createPtString(p, NULL, SYSTEM_OUT, displayRadians);
	printf("%s", output);
}

/*
ErrorSet createGeoTestBndry(LLPoint center, double radius, Boundary* geoTestBndry){

	double DEG2RAD = M_PI/180;

	LLPoint v1, v2, v3, v4;
	double crs1 = 0*DEG2RAD;
	double crs2 = 90*DEG2RAD;
	double crs3 = 180*DEG2RAD;
	double crs4 = 270*DEG2RAD;
	Geodesic geo12, geo23, geo34, geo41;

	ErrorSet err = 0;

	err |= direct(center, crs1, radius, &v1, EPS);
	err |= direct(center, crs2, radius, &v2, EPS);
	err |= direct(center, crs3, radius, &v3, EPS);
	err |= direct(center, crs4, radius, &v4, EPS);

	err |= createGeo(&geo12, v1, v2, SEGMENT, EPS);
	err |= createGeo(&geo23, v2, v3, SEGMENT, EPS);
	err |= createGeo(&geo34, v3, v4, SEGMENT, EPS);
	err |= createGeo(&geo41, v4, v1, SEGMENT, EPS);

	err |= addGeoToBndry(geoTestBndry, &geo12);
	err |= addGeoToBndry(geoTestBndry, &geo23);
	err |= addGeoToBndry(geoTestBndry, &geo34);
	err |= addGeoToBndry(geoTestBndry, &geo41);

	return err;
}

ErrorSet createLocusTestBndry(LLPoint center, double radius, Boundary* locusTestBndry){

	double DEG2RAD = M_PI/180;

	LLPoint v1, v2, v3, v4;
	double dst12, dst23, dst34, dst41;
	double crs1 = 0*DEG2RAD;
	double crs2 = 90*DEG2RAD;
	double crs3 = 180*DEG2RAD;
	double crs4 = 270*DEG2RAD;
	Locus locus12, locus23, locus34, locus41;

	ErrorSet err = 0;

	err |= direct(center, crs1, radius, &v1, EPS);
	err |= direct(center, crs2, radius, &v2, EPS);
	err |= direct(center, crs3, radius, &v3, EPS);
	err |= direct(center, crs4, radius, &v4, EPS);

	err |= inverse(v1, v2, NULL, NULL, &dst12, EPS);
	err |= inverse(v2, v3, NULL, NULL, &dst23, EPS);
	err |= inverse(v3, v4, NULL, NULL, &dst34, EPS);
	err |= inverse(v4, v1, NULL, NULL, &dst41, EPS);

	err |= createLocus(&locus12, v4, v3, -dst41, -dst23, SEGMENT, TOL, EPS);
	err |= createLocus(&locus23, v1, v4, -dst12, -dst34, SEGMENT, TOL, EPS);
	err |= createLocus(&locus34, v2, v1, -dst23, -dst41, SEGMENT, TOL, EPS);
	err |= createLocus(&locus41, v3, v2, -dst34, -dst12, SEGMENT, TOL, EPS);

	err |= addLocusToBndry(locusTestBndry, &locus12);
	err |= addLocusToBndry(locusTestBndry, &locus23);
	err |= addLocusToBndry(locusTestBndry, &locus34);
	err |= addLocusToBndry(locusTestBndry, &locus41);

	return err;
}

ErrorSet createArcTestBndry(LLPoint center, double radius, Boundary* arcTestBndry){

	double DEG2RAD = M_PI/180;

	LLPoint v1, v2, v3, v4;
	double crs1 = 0*DEG2RAD;
	double crs2 = 90*DEG2RAD;
	double crs3 = 180*DEG2RAD;
	double crs4 = 270*DEG2RAD;
	Arc arc12, arc23, arc34, arc41;

	ErrorSet err = 0;

	err |= direct(center, crs1, radius, &v1, EPS);
	err |= direct(center, crs2, radius, &v2, EPS);
	err |= direct(center, crs3, radius, &v3, EPS);
	err |= direct(center, crs4, radius, &v4, EPS);

	err |= createArc(&arc12, center, v1, v2, CLOCKWISE, TOL, EPS);
	err |= createArc(&arc23, center, v2, v3, CLOCKWISE, TOL, EPS);
	err |= createArc(&arc34, center, v3, v4, CLOCKWISE, TOL, EPS);
	err |= createArc(&arc41, center, v4, v1, CLOCKWISE, TOL, EPS);

	err |= addArcToBndry(arcTestBndry, &arc12);
	err |= addArcToBndry(arcTestBndry, &arc23);
	err |= addArcToBndry(arcTestBndry, &arc34);
	err |= addArcToBndry(arcTestBndry, &arc41);

	return err;
}

ErrorSet createSpiralTestBndry(LLPoint center, double radius1, double radius2, Boundary* spiralTestBndry){

	double DEG2RAD = M_PI/180;

	double crs1 = 0*DEG2RAD;
	double crs2 = 90*DEG2RAD;
	double crs3 = 180*DEG2RAD;
	double crs4 = 270*DEG2RAD;
	Spiral spiral12, spiral23, spiral34, spiral41;

	ErrorSet err = 0;

	err |= createSpiral(&spiral12, center, radius1, radius2, crs1, crs2, CLOCKWISE, EPS);
	err |= createSpiral(&spiral23, center, radius2, radius1, crs2, crs3, CLOCKWISE, EPS);
	err |= createSpiral(&spiral34, center, radius1, radius2, crs3, crs4, CLOCKWISE, EPS);
	err |= createSpiral(&spiral41, center, radius2, radius1, crs4, crs1, CLOCKWISE, EPS);

	err |= addSpiralToBndry(spiralTestBndry, &spiral12);
	err |= addSpiralToBndry(spiralTestBndry, &spiral23);
	err |= addSpiralToBndry(spiralTestBndry, &spiral34);
	err |= addSpiralToBndry(spiralTestBndry, &spiral41);

	return err;
}
*/
} //namespace

