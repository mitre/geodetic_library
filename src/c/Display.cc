/*
 * Copyright 2007-2021 The MITRE Corporation.  All Rights reserved.
 *
 * This is the copyright work of The MITRE Corporation and was produced
 * for the U.S. Government under Contract Number DTFAWA-10-C-00080
 * and is subject to Federal Aviation Administration Acquisition
 * Management System Clause 3.5-13, Rights in Data-General, Alt. III
 * and Alt. IV (Oct. 1996).  No other use than that granted to the U.S.
 * Government, or to those acting on behalf of the U.S. Government, under
 * that Clause is authorized without the express written permission of
 * The MITRE Corporation.  
 * For further information, please contact The MITRE Corporation,
 * Contracts Office, 7515 Colshire Drive, McLean, VA 22102 (703) 983-6000.
 */

#include <sys/types.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#if REPLACE_WITH_AMDLIBM
#include "amdlibm.h"
#endif
#include "Util.h"

using namespace geolib_idealab;

/* Print out a list of constants that were set at compile time */
void geolib_idealab::reportConstants(const char* fileName)
{
    FILE *F = NULL;
    char remChar = '#';

    if (NULL != fileName)
    {
        F = fopen(fileName, "w");
    }
    /* If fopen fails or no name given, write to STDOUT */
    if (NULL == F)
    {
        F = stdout;
    }

#ifdef VERSION
    // VERSION macro should be defined at compile time by passing -D flag to compiler
    fprintf(F, "%c Repository Version: %s\n", remChar, strcmp(VERSION, "") ? VERSION : "Unknown");
#else
    fprintf(F,"%c No VERSION macro defined; repository version unknown\n",remChar);
#endif
/*
#ifdef TARGETS
    fprintf(F,"%c !! TARGETS flag set at compile time--LOW PRECISION ONLY !!\n",remChar);
#else
    fprintf(F,"%c Compiled for high precision\n",remChar);
#endif

    fprintf(F, "%c GEOLIB Constants\n", remChar);

#ifdef KEEP_USER_EPS
    fprintf(F,
            "%c  No DEFAULT_EPS imposed; user-supplied value will be used\n",
            remChar);
#else
    fprintf(F,"%c  DEFAULT_EPS used for Forward/Inverse computation: %.1e\n",remChar, DEFAULT_EPS);
#endif
*/
    fprintf(F, "%c  FLATTENING:             %0.17f\n", remChar, FLATTENING);
    fprintf(F, "%c  INVERSE_FLATTENING:     %-16f\n", remChar,
            INVERSE_FLATTENING);
    fprintf(F, "%c  SEMI_MAJOR_AXIS_FEET:   %-.7f\n", remChar,
            SEMI_MAJOR_AXIS_FEET);
    fprintf(F, "%c  SEMI_MAJOR_AXIS_METERS: %-16f\n", remChar,
            SEMI_MAJOR_AXIS_METERS);
    fprintf(F, "%c  SEMI_MAJOR_AXIS_NMI:    %-16f\n", remChar,
            SEMI_MAJOR_AXIS_NMI);
    fprintf(F, "%c  SEMI_MINOR_AXIS_NMI:    %-16f\n", remChar,
            SEMI_MINOR_AXIS_NMI);
    fprintf(F, "%c  FEET_PER_NMI:           %-16f\n", remChar, FEET_PER_NMI);
    fprintf(F, "%c  SPHERE_RADIUS_NMI:      %-16f\n", remChar, SPHERE_RADIUS);
    fprintf(F, "%c  FLAT_TOL:               %-.1e\n", remChar, FLAT_TOL);

    fprintf(F, "%c  MAX_ITERATIONS:         %d\n", remChar, MAX_ITERATIONS);
    fprintf(F, "%c  MAX_DISTANCE_ERROR:     %-16f\n", remChar,
            (double) MAX_DISTANCE_ERROR);
    fprintf(F, "%c  MAX_V_ITERATIONS:       %d\n", remChar, MAX_V_ITERATIONS);
    fprintf(F, "%c  ANTIPODAL_TOL:          %-16f\n", remChar, ANTIPODAL_TOL
            * 180.0 / M_PI);
    fprintf(F, "%c  M_PI:                   %0.24lf\n", remChar, M_PI);
    fprintf(F, "%c  M_PI_2:                 %0.24lf\n", remChar, M_PI_2);
    fprintf(F, "%c  M_2PI:                  %0.24lf\n", remChar, M_2PI);
    fprintf(F, "%c \n", remChar);
    fprintf(F, "%c  Double is %d bytes\n", remChar, sizeof(double));
    fprintf(F, "%c  Long Double is %d bytes\n", remChar, sizeof(long double));
    if (sizeof(double) == sizeof(SPECIAL_DOUBLE))
    {
        fprintf(F, "%c  Long doubles are not being used in iterateInverse\n",
                remChar);
    }

    if ((NULL != F) && (stdout != F))
    {
        fclose(F);
    }

}

/* Print out a description of Geolib's current version */
void geolib_idealab::aboutGeolib(char *result, const unsigned int size)
{
	char svnVer[128];
	char geoVer[128];
	char algVer[128];
	char contact[256];
	char cr[1024];
	char output[size];
	int i;

#ifdef VERSION
    // VERSION macro should be defined at compile time by passing -D flag to compiler
    snprintf(svnVer, sizeof svnVer, "Repository Version: %s\n", strcmp(VERSION, "") ? VERSION : "Unknown");
//#else
//    snprintf(svnVer, sizeof svnVer, "No VERSION macro defined; repository version unknown\n");
#endif
//    snprintf(geoVer, sizeof geoVer, "Current Version: Geolib%s", GEOLIB_VERSION_STR);
//    snprintf(algVer, sizeof algVer, "Current Algorithm: %s\n", ALGORITHM);
    snprintf(contact, sizeof contact, "Contacts: Mike Mills (mjmills@mitre.org)\n          Joe Hopper (jhopper@mitre.org)\n          Joe Heidrich (jheidrich@mitre.org)\n\n");
    snprintf(cr, sizeof cr, "Copyright 2007-2012 The MITRE Corporation.  All Rights reserved.\n\nThis is the copyright work of The MITRE Corporation and was produced\nfor the U.S. Government under Contract Number DTFAWA-10-C-00080\nand is subject to Federal Aviation Administration Acquisition Management System Clause\n3.5-13, Rights in Data-General, Alt. III and Alt. IV (Oct. 1996).  No other use than\nthat granted to the U.S. Government, or to those acting on behalf of\nthe U.S. Government, under that Clause is authorized without the\nexpress written permission of The MITRE Corporation.\n\nFor further information, please contact The MITRE Corporation, Contracts Office, 7515 Colshire Drive, McLean, VA 22102 (703) 983-6000.\n");
    snprintf(output, sizeof output, "%s\n%s\n%s\n%s\n%s", geoVer,svnVer,algVer,cr,contact);


    for (i=0; i < size; i++)
    {
        result[i]=output[i];
	}
}

char* geolib_idealab::createPtString(LLPoint p, char* pointName, OutputMode mode, int displayRadians)
{
	char* string = static_cast<char *>(calloc(sizeof(char), 500));

    double lat = p.latitude;
    double lon = p.longitude;
    char* angleUnits = "rad";
    char conversion[10] = "";
    const int maxNameLen = 100;
    char name[101];//+1 is for terminating character '\0'
    char strNameLen[4];
    char formatString[10];

    //check for conversion to degrees
    if (!displayRadians)
    {
        lat *= 180.0 / M_PI;
        lon *= 180.0 / M_PI;
        angleUnits = "deg";
        sprintf(conversion,"*(pi/180)");
    }

    //truncate variable name if necessary to max length
    if (pointName != NULL)
    {
		sprintf(strNameLen, "%i", maxNameLen);
		sprintf(formatString, "%%.%ss", strNameLen);
		sprintf(name, formatString, pointName);
    } else {
    	sprintf(name, "%s", "Point");
    }


    switch (mode)
    {
		case SYSTEM_OUT:
		    sprintf(string, "\n%s: (%.20lf, %.20lf) [%s]", name, lat, lon, angleUnits);
			break;
		case MATLAB_OUT:
			sprintf(string, "\n%s = CoordinatePoint([%.20lf %.20lf]%s); %%[%s]",
					name, lat, lon, conversion, angleUnits);
			break;
		case JAVA_OUT:
			break;
		case DATA_OUT:
			sprintf(string, "\nLLPoint %s = {%.20lf, %.20lf};", name, lat, lon);
			break;
    }


	return string;
}

char* geolib_idealab::createGeoString(Geodesic g, char* geoName, OutputMode mode, int displayRadians)
{
	char* string = static_cast<char *>(calloc(sizeof(char), 2000));
    char conversion[10] = "";
    const int maxNameLen = 40;
    char name[41];//+1 is for terminating character '\0'
    char strNameLen[3];
    char formatString[10];

    char startPtName[81];
    char endPtName[81];
    char* strStartPt;
    char* strEndPt;
    char strForAz[100];
    char strRevAz[100];
    char strEndAz[100];
    char strGeoLen[100];
    char strGeoType[100];
    char strGeo[100];
    double startAz = g.startAz;
    double endAz = g.endAz;
    double revAz = modpos(g.endAz + M_PI, M_2PI );
    char* angleUnits = "rad";

    if (!displayRadians)
    {
        startAz *= 180.0 / M_PI;
        endAz *= 180.0 / M_PI;
        revAz = modpos(((g.endAz) * 180.0 / M_PI) + 180.0, 360.0);
        angleUnits = "deg";
        sprintf(conversion,"*(pi/180)");
    }

    //truncate variable name if necessary to max length
    if (geoName != NULL)
    {
		sprintf(strNameLen, "%i", maxNameLen);
		sprintf(formatString, "%%.%ss", strNameLen);
		sprintf(name, formatString, geoName);
    } else {
    	sprintf(name, "%s", "Geodesic");
    }


    switch (mode)
    {
		case SYSTEM_OUT:
			sprintf(startPtName, "Start Point");
			strStartPt = geolib_idealab::createPtString(g.startPoint, startPtName, SYSTEM_OUT, displayRadians);

			sprintf(endPtName, "End Point");
			strEndPt = geolib_idealab::createPtString(g.endPoint, endPtName, SYSTEM_OUT, displayRadians);

		    sprintf(strForAz, "\nForward Azimuth: %.20lf [%s]", startAz, angleUnits);

		    sprintf(strRevAz, "\nReverse Azimuth: %.20lf [%s]", revAz, angleUnits);

		    sprintf(strEndAz, "\nEnd Azimuth: %.20lf [%s]", endAz, angleUnits);

		    sprintf(strGeoLen, "\nLength: %.20lf [nautical miles]", g.length);

		    sprintf(strGeoType, "\nLine Type [0=Segment, 1=Semi, 2=Inf]: %i", g.lineType);

		    sprintf(string, "\n%s:%s%s%s%s%s%s%s", name, strStartPt, strEndPt, strForAz, strRevAz, strEndAz, strGeoLen, strGeoType);
			break;
		case MATLAB_OUT:
			sprintf(startPtName, "%sStart", name);
			strStartPt = geolib_idealab::createPtString(g.startPoint, startPtName, MATLAB_OUT, displayRadians);

			sprintf(endPtName, "%sEnd", name);
			strEndPt = geolib_idealab::createPtString(g.endPoint, endPtName, MATLAB_OUT, displayRadians);

			sprintf(strGeo, "\n%s = LineSegment(%s, %s);", name, startPtName, endPtName);

			sprintf(string, "\n%s%s%s\n", strStartPt, strEndPt, strGeo);
			break;
		case JAVA_OUT:
			break;
		case DATA_OUT:
			sprintf(startPtName, "%sStart", name);
			strStartPt = geolib_idealab::createPtString(g.startPoint, startPtName, DATA_OUT, displayRadians);

			sprintf(endPtName, "%sEnd", name);
			strEndPt = geolib_idealab::createPtString(g.endPoint, endPtName, DATA_OUT, displayRadians);

			sprintf(strGeo, "\nGeodesic %s;\ncreateGeo(&%s, %s, %s, SEGMENT, eps);", name, name, startPtName, endPtName);

			sprintf(string, "\n%s%s%s\n", strStartPt, strEndPt, strGeo);
			break;
    }


	return string;
}

char* geolib_idealab::createArcString(Arc a, char* arcName, OutputMode mode, int displayRadians)
{
	char* string = static_cast<char *>(calloc(sizeof(char), 3000));
    char conversion[10] = "";
    const int maxNameLen = 40;
    char name[41];//+1 is for terminating character '\0'
    char strNameLen[3];
    char formatString[10];


    char centerPtName[81];
    char startPtName[81];
    char endPtName[81];
    char radiusName[81];
    char startAngleName[81];
    char strStartAngle[100];
    char subAngleName[81];
    char* strCenterPt;
    char* strStartPt;
    char* strEndPt;
    char strStartAz[100];
    char strEndAz[100];
    char strSubAngle[100];
    char strRadius[200];
    char strDir[100];
    char strArc[441];
    double startAz = a.startAz;
    double endAz = a.endAz;
    double subtendedAngle = a.subtendedAngle;
    double midPtCrs;
    char* angleUnits = "rad";

	//compute the course from the arc center to the arc mid point
	midPtCrs = modpos(a.startAz + (0.5*a.subtendedAngle), M_2PI);//subtended angle is signed already

    if (!displayRadians)
    {
        startAz *= 180.0 / M_PI;
        endAz *= 180.0 / M_PI;
        subtendedAngle *= 180.0 / M_PI;
        midPtCrs *= 180.0 / M_PI;
        angleUnits = "deg";
        sprintf(conversion,"*(pi/180)");
    }

    //truncate variable name if necessary to max length
    if (arcName != NULL)
    {
		sprintf(strNameLen, "%i", maxNameLen);
		sprintf(formatString, "%%.%ss", strNameLen);
		sprintf(name, formatString, arcName);
    } else {
    	sprintf(name, "%s", "Arc");
    }


    switch (mode)
    {
		case SYSTEM_OUT:
			sprintf(centerPtName, "Center Point");
			strCenterPt = geolib_idealab::createPtString(a.centerPoint, centerPtName, SYSTEM_OUT, displayRadians);

			sprintf(startPtName, "Start Point");
			strStartPt = geolib_idealab::createPtString(a.startPoint, startPtName, SYSTEM_OUT, displayRadians);

			sprintf(endPtName, "End Point");
			strEndPt = geolib_idealab::createPtString(a.endPoint, endPtName, SYSTEM_OUT, displayRadians);

			if(a.dir == COUNTERCLOCKWISE){
				sprintf(strDir, "\nDirection: %d [counterclockwise]", a.dir);
			} else {
				sprintf(strDir, "\nDirection: %d [clockwise]", a.dir);
			}

		    sprintf(strStartAz, "\nStart Azimuth: %.20lf [%s]", startAz, angleUnits);

		    sprintf(strEndAz, "\nEnd Azimuth: %.20lf [%s]", endAz, angleUnits);

		    sprintf(strRadius, "\nRadius: %.20lf [nm]", a.radius);

		    sprintf(strSubAngle, "\nSubtended Angle: %.20lf [%s]", subtendedAngle, angleUnits);

		    sprintf(string, "\n%s:%s%s%s%s%s%s%s%s", name, strCenterPt, strStartPt, strEndPt, strDir, strStartAz, strEndAz, strRadius, strSubAngle);
			break;
		case MATLAB_OUT:

			sprintf(centerPtName, "%sCenter", name);
			strCenterPt = geolib_idealab::createPtString(a.centerPoint, centerPtName, MATLAB_OUT, displayRadians);

			sprintf(radiusName, "%sRadius", name);
			sprintf(strRadius, "\n%s = %.20lf; %%[nm]", radiusName, a.radius);

			if (a.startPoint.latitude == a.endPoint.latitude && a.startPoint.longitude == a.endPoint.longitude)
			{

				sprintf(strArc, "\n%s = Arc(%s, %s);", name, centerPtName, radiusName);

				sprintf(string, "\n%s%s%s", strCenterPt, strRadius, strArc);
			}
			else
			{
			    sprintf(startAngleName, "%sStartAngle", name);
			    sprintf(strStartAngle, "\n%s = %.20lf%s; %%[%s]", startAngleName, startAz, conversion, angleUnits);

			    sprintf(subAngleName, "%sSubtendedAngle", name);
			    sprintf(strSubAngle, "\n%s = %.20lf%s; %%[%s]", subAngleName, subtendedAngle, conversion, angleUnits);

				sprintf(strArc, "\n%s = Arc(%s, %s, %s, %s);", name, centerPtName, radiusName, startAngleName, subAngleName);

				sprintf(string, "\n%s%s%s%s%s\n", strCenterPt, strRadius, strStartAngle, strSubAngle, strArc);
			}

			break;
		case JAVA_OUT:
			break;
		case DATA_OUT:
			sprintf(centerPtName, "%sCenter", name);
			strCenterPt = geolib_idealab::createPtString(a.centerPoint, centerPtName, DATA_OUT, displayRadians);
			sprintf(startPtName, "%sStartPoint", name);
			strStartPt = geolib_idealab::createPtString(a.startPoint, startPtName, DATA_OUT, displayRadians);
			sprintf(endPtName, "%sEndPoint", name);
			strEndPt = geolib_idealab::createPtString(a.endPoint, endPtName, DATA_OUT, displayRadians);

			sprintf(radiusName, "%sRadius", name);
			sprintf(strRadius, "\ndouble %s = %.20lf;", radiusName, a.radius);

			if(a.dir == CLOCKWISE)
				sprintf(strDir, "CLOCKWISE");
			else
				sprintf(strDir, "COUNTERCLOCKWISE");

				sprintf(strArc, "\nArc %s;\ncreateArc(&%s, %s, %s, %s, %s, tol, eps);",
						name, name, centerPtName, startPtName, endPtName, strDir);

				sprintf(string, "\n%s%s%s%s%s", strCenterPt, strStartPt, strEndPt, strRadius, strArc);

			break;
    }


	return string;
}

char* geolib_idealab::createSpiralString(Spiral sp, char* spiralName, OutputMode mode, int displayRadians)
{
	char* string = static_cast<char *>(calloc(sizeof(char), 3000));
    char conversion[10] = "";
    const int maxNameLen = 40;
    char name[41];//+1 is for terminating character '\0'
    char strNameLen[3];
    char formatString[10];

    char centerPtName[81];
    char startPtName[81];
    char endPtName[81];
    char orientationName[81];
    char startRadiusName[81];
    char endRadiusName[81];
    char startAzName[81];
    char endAzName[81];
    char* strCenterPt;
    char* strStartPt;
    char* strEndPt;
    char strStartAz[100];
    char strEndAz[100];
    char strSubAngle[100];
    char strStartRadius[100];
    char strEndRadius[100];
    char strDir[100];
    char strGrowthRate[100];
    char strOrientation[100];
    char strSpiral[441];
    double startAz = sp.startAz;
    double endAz = sp.endAz;
    double subtendedAngle = sp.subtendedAngle;
    double growthRate = sp.growthRate;
    char* angleUnits = "rad";

    if (!displayRadians)
    {
        startAz *= 180.0 / M_PI;
        endAz *= 180.0 / M_PI;
        subtendedAngle *= 180.0 / M_PI;
        growthRate *=  M_PI / 180.0;
        angleUnits = "deg";
        sprintf(conversion,"*(pi/180)");
    }

    //truncate variable name if necessary to max length
    if (spiralName != NULL)
    {
		sprintf(strNameLen, "%i", maxNameLen);
		sprintf(formatString, "%%.%ss", strNameLen);
		sprintf(name, formatString, spiralName);
    } else {
    	sprintf(name, "%s", "Spiral");
    }


    switch (mode)
    {
		case SYSTEM_OUT:
			sprintf(centerPtName, "Center Point");
			strCenterPt = geolib_idealab::createPtString(sp.centerPoint, centerPtName, SYSTEM_OUT, displayRadians);

			sprintf(startPtName, "Start Point");
			strStartPt = geolib_idealab::createPtString(sp.startPoint, startPtName, SYSTEM_OUT, displayRadians);

			sprintf(endPtName, "End Point");
			strEndPt = geolib_idealab::createPtString(sp.endPoint, endPtName, SYSTEM_OUT, displayRadians);

			if(sp.dir == COUNTERCLOCKWISE){
				sprintf(strDir, "\nDirection: %d [counterclockwise]", sp.dir);
			} else {
				sprintf(strDir, "\nDirection: %d [clockwise]", sp.dir);
			}

		    sprintf(strStartAz, "\nStart Azimuth: %.20lf [%s]", startAz, angleUnits);

		    sprintf(strEndAz, "\nEnd Azimuth: %.20lf [%s]", endAz, angleUnits);

		    sprintf(strStartRadius, "\nStart Radius: %.20lf [nm]", sp.startRadius);

		    sprintf(strEndRadius, "\nEnd Radius: %.20lf [nm]", sp.endRadius);

		    sprintf(strSubAngle, "\nSubtended Angle: %.20lf [%s]", subtendedAngle, angleUnits);

		    sprintf(strGrowthRate, "\nGrowth Rate: %.20lf [nm/%s]", growthRate, angleUnits);

		    sprintf(string, "\n%s:%s%s%s%s%s%s%s%s%s%s", name, strCenterPt, strStartPt, strEndPt, strDir, strStartAz, strEndAz, strStartRadius, strEndRadius, strSubAngle, strGrowthRate);
			break;
		case MATLAB_OUT:

			sprintf(centerPtName, "%sCenter", name);
			strCenterPt = geolib_idealab::createPtString(sp.centerPoint, centerPtName, MATLAB_OUT, displayRadians);

			sprintf(startPtName, "%sStartPt", name);
			strStartPt = geolib_idealab::createPtString(sp.startPoint, startPtName, MATLAB_OUT, displayRadians);

			sprintf(endPtName, "%sEndPt", name);
			strEndPt = geolib_idealab::createPtString(sp.endPoint, endPtName, MATLAB_OUT, displayRadians);

			sprintf(orientationName, "%sOrientation", name);
			if(sp.dir == COUNTERCLOCKWISE){
			        sprintf(strOrientation, "\n%s = %d;", orientationName, 1); // Reversed for MATLAB
			} else {
			        sprintf(strOrientation, "\n%s = %d;", orientationName, -1); // Reversed for MATLAB
			}

				sprintf(strSpiral, "\n%s = Spiral(%s, %s, %s, %s);", name, centerPtName, startPtName, endPtName, orientationName);

				sprintf(string, "\n%s%s%s%s%s\n", strCenterPt, strStartPt, strEndPt, strOrientation, strSpiral);

			break;
		case JAVA_OUT:
			break;
		case DATA_OUT:
			sprintf(centerPtName, "%sCenter", name);
			strCenterPt = geolib_idealab::createPtString(sp.centerPoint, centerPtName, DATA_OUT, displayRadians);

			sprintf(startAzName, "%sStartAzimuth", name);
		    sprintf(strStartAz, "\ndouble %s = %.20lf;", startAzName, startAz);
		    sprintf(endAzName, "%sEndAzimuth", name);
		    sprintf(strEndAz, "\ndouble %s = %.20lf;", endAzName, endAz);
		    sprintf(startRadiusName, "%sStartRadius", name);
		    sprintf(strStartRadius, "\ndouble %s = %.20lf;", startRadiusName, sp.startRadius);
		    sprintf(endRadiusName, "%sEndRadius", name);
		    sprintf(strEndRadius, "\ndouble %s = %.20lf;", endRadiusName, sp.endRadius);


			if(sp.dir == CLOCKWISE)
				sprintf(strDir, "CLOCKWISE");
			else
				sprintf(strDir, "COUNTERCLOCKWISE");

				sprintf(strSpiral, "\nSpiral %s;\ncreateSpiral(&%s, %s, %s, %s, %s, %s, %s, eps);",
						name, name, centerPtName, startRadiusName, endRadiusName, startAzName, endAzName, strDir);

				sprintf(string, "\n%s%s%s%s%s%s", strCenterPt, strStartRadius, strEndRadius, strStartAz, strEndAz, strSpiral);
			break;
    }

	return string;
}

char* geolib_idealab::createLocusString(Locus l, char* locusName, OutputMode mode, int displayRadians)
{
	char* string = static_cast<char *>(calloc(sizeof(char), 3000));
    char conversion[10] = "";
    const int maxNameLen = 40;
    char name[41];//+1 is for terminating character '\0'
    char strNameLen[3];
    char formatString[10];

    char geoStartPtName[81];
    char geoEndPtName[81];
    char locStartPtName[81];
    char locEndPtName[81];
    char locStartDistName[81];
    char locEndDistName[81];
    char* strGeoStartPt;
    char* strGeoEndPt = NULL;
    char* strLocStartPt;
    char* strLocEndPt;
    char strForAz[100];
    char strRevAz[100];
    char strEndAz[100];
    char strGeoLen[100];
    char strLocType[100];
    char strLocStartDist[100];
    char strLocEndDist[100];
    char strSlope[100];
    char strGeo[100];
    char strLocus[100];
    double startAz = l.geoAz;
    double endAz = modpos(l.geoRevAz + M_PI, M_2PI );
    double revAz = l.geoRevAz;
    char* angleUnits = "rad";

    if (!displayRadians)
    {
        startAz *= 180.0 / M_PI;
        revAz *= 180.0 / M_PI;
        endAz = modpos(((l.geoRevAz) * 180.0 / M_PI) + 180.0, 360.0);
        angleUnits = "deg";
        sprintf(conversion,"*(pi/180)");
    }

    //truncate variable name if necessary to max length
    if (locusName != NULL)
    {
		sprintf(strNameLen, "%i", maxNameLen);
		sprintf(formatString, "%%.%ss", strNameLen);
		sprintf(name, formatString, locusName);
    } else {
    	sprintf(name, "%s", "Locus");
    }


    switch (mode)
    {
		case SYSTEM_OUT:
			sprintf(geoStartPtName, "Geo Start Point");
			strGeoStartPt = geolib_idealab::createPtString(l.geoStart, geoStartPtName, SYSTEM_OUT, displayRadians);

			sprintf(geoEndPtName, "Geo End Point");
			strGeoEndPt = geolib_idealab::createPtString(l.geoEnd, geoEndPtName, SYSTEM_OUT, displayRadians);

		    sprintf(strForAz, "\nForward Azimuth: %.20lf [%s]", startAz, angleUnits);

		    sprintf(strRevAz, "\nReverse Azimuth: %.20lf [%s]", revAz, angleUnits);

		    sprintf(strEndAz, "\nEnd Azimuth: %.20lf [%s]", endAz, angleUnits);

		    sprintf(strGeoLen, "\nGeo Length: %.20lf [nautical miles]", l.geoLength);

		    sprintf(strLocType, "\nLine Type [0=Segment, 1=Semi, 2=Inf]: %i", l.lineType);

		    sprintf(strLocStartDist, "\nStart Distance: %.20lf [nm]", l.startDist);

		    sprintf(strLocEndDist, "\nEnd Distance: %.20lf [nm]", l.endDist);

		    sprintf(strSlope, "\nSlope: %.20lf", l.slope);

			sprintf(locStartPtName, "Locus Start Point");
			strLocStartPt = geolib_idealab::createPtString(l.locusStart, locStartPtName, SYSTEM_OUT, displayRadians);

			sprintf(locEndPtName, "Locus End Point");
			strLocEndPt = geolib_idealab::createPtString(l.locusEnd, locEndPtName, SYSTEM_OUT, displayRadians);

		    sprintf(string, "\n%s:%s%s%s%s%s%s%s%s%s%s%s%s", name, strGeoStartPt, strGeoEndPt, strForAz, strRevAz, strEndAz,
		    		strGeoLen, strLocType, strLocStartDist, strLocEndDist, strSlope, strLocStartPt, strLocEndPt);
			break;
		case MATLAB_OUT:
			sprintf(geoStartPtName, "%sGeoStart", name);
			strGeoStartPt = geolib_idealab::createPtString(l.geoStart, geoStartPtName, MATLAB_OUT, displayRadians);

			sprintf(geoEndPtName, "%sGeoEnd", name);
			strGeoEndPt = geolib_idealab::createPtString(l.geoEnd, geoEndPtName, MATLAB_OUT, displayRadians);

			sprintf(strGeo, "\n%sGeo = LineSegment(%s, %s);", name, geoStartPtName, geoEndPtName);

			sprintf(strLocus, "\n%s = Locus(%sGeo, %.20lf, %.20lf); %%[nm, nm]", name, name, l.startDist, l.endDist);

			sprintf(string, "\n%s%s%s%s\n", strGeoStartPt, strGeoEndPt, strGeo, strLocus);
			break;
		case JAVA_OUT:
			break;
		case DATA_OUT:
			sprintf(geoStartPtName, "%sGeoStart", name);
			strGeoStartPt = geolib_idealab::createPtString(l.geoStart, geoStartPtName, DATA_OUT, displayRadians);

			sprintf(geoEndPtName, "%sGeoEnd", name);
			strGeoEndPt = geolib_idealab::createPtString(l.geoEnd, geoEndPtName, DATA_OUT, displayRadians);

			sprintf(locStartDistName, "%sStartDistance", name);
		    sprintf(strLocStartDist, "\ndouble %s = %.20lf;", locStartDistName, l.startDist);
			sprintf(locEndDistName, "%sEndDistance", name);
		    sprintf(strLocEndDist, "\ndouble %s = %.20lf;", locEndDistName, l.endDist);

		    if(l.lineType == SEGMENT)
		    	sprintf(strLocType, "SEGMENT");
		    else if(l.lineType == SEMIINFINITE)
		    	sprintf(strLocType, "SEMIINFINITE");
		    else
		    	sprintf(strLocType, "INFINITE");

			sprintf(strLocus, "\nLocus %s;\ncreateLocus(&%s, %s, %s, %s, %s, %s, tol, eps);",
					name, name, geoStartPtName, geoEndPtName, locStartDistName, locEndDistName, strLocType);

			sprintf(string, "\n%s%s%s%s%s\n", strGeoStartPt, strGeoEndPt, strLocStartDist, strLocEndDist, strLocus);


			break;
    }


	return string;
}

char* geolib_idealab::createBndryString(Boundary b, char* boundaryName, OutputMode mode, int displayRadians)
{
    char* buff = NULL;
    Shape thisShape;
    char* shapeStringList[50];
    int shapeCount = 0;
    int charCount = 0;
    int i = 0;

#define maxNameLen 40
//    const int maxNameLen = 40;
    char name[41];//+1 is for terminating character '\0'
    char strNameLen[3];
    char formatString[10];
    char shapeName[2*maxNameLen];

    //truncate variable name if necessary to max length
    if (boundaryName != NULL)
    {
		sprintf(strNameLen, "%i", maxNameLen);
		sprintf(formatString, "\n%%.%ss\n", strNameLen);
		sprintf(name, formatString, boundaryName);
    } else {
    	sprintf(name, "%s", "\nBoundary\n");
    }

    for (i = 0; i < b.length; i++)
    {
        thisShape = b.elements[i];
        switch (thisShape.type) {
        case ARC:
        	sprintf(shapeName, "%s_Arc_%i", boundaryName, i);
            shapeStringList[shapeCount] = geolib_idealab::createArcString(*((Arc*) thisShape.this_shape), shapeName, mode, displayRadians);
            charCount = charCount + strlen(shapeStringList[shapeCount++]);
            break;
        case GEODESIC:
        	sprintf(shapeName, "%s_Geodesic_%i", boundaryName, i);
            shapeStringList[shapeCount] = geolib_idealab::createGeoString(*((Geodesic*) thisShape.this_shape), shapeName, mode, displayRadians);
            charCount = charCount + strlen(shapeStringList[shapeCount++]);
            break;
        case LOCUS:
        	sprintf(shapeName, "%s_Locus_%i", boundaryName, i);
            shapeStringList[shapeCount] = geolib_idealab::createLocusString(*((Locus*) thisShape.this_shape), shapeName, mode, displayRadians);
            charCount = charCount + strlen(shapeStringList[shapeCount++]);
            break;
        case SPIRAL:
        	sprintf(shapeName, "%s_Spiral_%i", boundaryName, i);
            shapeStringList[shapeCount] = geolib_idealab::createSpiralString(*((Spiral*) thisShape.this_shape), shapeName, mode, displayRadians);
            charCount = charCount + strlen(shapeStringList[shapeCount++]);
            break;
        case LLPOINT:
        	sprintf(shapeName, "%s_Point_%i", boundaryName, i);
            shapeStringList[shapeCount] = geolib_idealab::createPtString(*((LLPoint*) thisShape.this_shape), shapeName, mode, displayRadians);
            charCount = charCount + strlen(shapeStringList[shapeCount++]);
            break;
        }
    }

    buff = (char*) calloc(sizeof(char), charCount + 1 + maxNameLen);

    strcat(buff, name);
    for (i = 0; i < shapeCount; i++)
    {
        buff = strcat(buff, shapeStringList[i]);
    }

    return buff;

}

char* geolib_idealab::createComplexBndryString(ComplexBoundary c, char* ComplexBoundaryName, OutputMode mode, int displayRadians)
{
    char* buff = NULL;
    Boundary* boundary = NULL;
    char* boundaryStringList[50];
    int elementCount = 0;
    int charCount = 0;
    int i = 0;

//    const int maxNameLen = 40;
#define maxNameLen 40
    char name[41];//+1 is for terminating character '\0'
    char strNameLen[3];
    char formatString[10];
    char boundaryName[maxNameLen];

    //truncate variable name if necessary to max length
    if (ComplexBoundaryName != NULL)
    {
		sprintf(strNameLen, "%i", maxNameLen);
		sprintf(formatString, "\n%%.%ss\n", strNameLen);
		sprintf(name, formatString, ComplexBoundaryName);
    } else {
    	sprintf(name, "%s", "\nComplexBoundary\n");
    }

    for (i = 0; i < c.length; i++)
    {
        boundary = c.elements[i];
		sprintf(boundaryName, "Boundary_%i", i);
		boundaryStringList[elementCount] = geolib_idealab::createBndryString(*boundary, boundaryName, mode, displayRadians);
		charCount = charCount + strlen(boundaryStringList[elementCount++]);
    }

    buff = (char*) calloc(sizeof(char), charCount + 1 + maxNameLen);

    strcat(buff, name);
    for (i = 0; i < elementCount; i++)
    {
        buff = strcat(buff, boundaryStringList[i]);
    }

    return buff;

}

void geolib_idealab::displayPt(LLPoint p, char* pointName, int displayRadians)
{
	char* output;

	output = geolib_idealab::createPtString(p, pointName, SYSTEM_OUT, displayRadians);
	printf("%s", output);
}

void geolib_idealab::displayGeo(Geodesic g, char* geoName, int displayRadians)
{
	char* output;

	output = geolib_idealab::createGeoString(g, geoName, SYSTEM_OUT, displayRadians);
	printf("%s", output);
}

void geolib_idealab::displayArc(Arc a, char* arcName, int displayRadians)
{
	char* output;

	output = geolib_idealab::createArcString(a, arcName, SYSTEM_OUT, displayRadians);
	printf("%s", output);
}

void geolib_idealab::displayLocus(Locus l, char* locusName, int displayRadians)
{
	char* output;

	output = geolib_idealab::createLocusString(l, locusName, SYSTEM_OUT, displayRadians);
	printf("%s", output);
}

void geolib_idealab::displaySpiral(Spiral s, char* spiralName, int displayRadians)
{
	char* output;

	output = geolib_idealab::createSpiralString(s, spiralName, SYSTEM_OUT, displayRadians);
	printf("%s", output);
}

void geolib_idealab::displayBndry(Boundary b, char* boundaryName, int displayRadians)
{
	char* output;

	output = geolib_idealab::createBndryString(b, boundaryName, SYSTEM_OUT, displayRadians);
	printf("%s", output);
}

void geolib_idealab::displayComplexBndry(ComplexBoundary c, char* complexBoundaryName, int displayRadians)
{
	char* output;

	output = geolib_idealab::createComplexBndryString(c, complexBoundaryName, SYSTEM_OUT, displayRadians);
	printf("%s", output);
}

void geolib_idealab::displayMatlabPt(LLPoint p, char* pointName, int displayRadians)
{
	char* output;

	output = geolib_idealab::createPtString(p, pointName, MATLAB_OUT, displayRadians);
	printf("%s", output);
}

void geolib_idealab::displayMatlabGeo(Geodesic g, char* geoName, int displayRadians)
{
	char* output;

	output = geolib_idealab::createGeoString(g, geoName, MATLAB_OUT, displayRadians);
	printf("%s", output);
}

void geolib_idealab::displayMatlabArc(Arc a, char* arcName, int displayRadians)
{
	char* output;

	output = geolib_idealab::createArcString(a, arcName, MATLAB_OUT, displayRadians);
	printf("%s", output);
}

void geolib_idealab::displayMatlabLocus(Locus l, char* locusName, int displayRadians)
{
	char* output;

	output = geolib_idealab::createLocusString(l, locusName, MATLAB_OUT, displayRadians);
	printf("%s", output);
}

void geolib_idealab::displayMatlabSpiral(Spiral s, char* spiralName, int displayRadians)
{
	char* output;

	output = geolib_idealab::createSpiralString(s, spiralName, MATLAB_OUT, displayRadians);
	printf("%s", output);
}

void geolib_idealab::displayMatlabBndry(Boundary b, char* boundaryName, int displayRadians)
{
	char* output;

	output = geolib_idealab::createBndryString(b, boundaryName, MATLAB_OUT, displayRadians);
	printf("%s", output);
}

void geolib_idealab::displayMatlabComplexBndry(ComplexBoundary c, char* complexBoundaryName, int displayRadians)
{
	char* output;

	output = geolib_idealab::createComplexBndryString(c, complexBoundaryName, MATLAB_OUT, displayRadians);
	printf("%s", output);
}

void geolib_idealab::displayDataPt(LLPoint p, char* pointName, int displayRadians)
{
	char* output;

	output = geolib_idealab::createPtString(p, pointName, DATA_OUT, 1);
	printf("%s", output);
}

void geolib_idealab::displayDataGeo(Geodesic g, char* geoName, int displayRadians)
{
	char* output;

	output = geolib_idealab::createGeoString(g, geoName, DATA_OUT, 1);
	printf("%s", output);
}

void geolib_idealab::displayDataArc(Arc a, char* arcName, int displayRadians)
{
	char* output;

	output = geolib_idealab::createArcString(a, arcName, DATA_OUT, 1);
	printf("%s", output);
}

void geolib_idealab::displayDataLocus(Locus l, char* locusName, int displayRadians)
{
	char* output;

	output = geolib_idealab::createLocusString(l, locusName, DATA_OUT, 1);
	printf("%s", output);
}

void geolib_idealab::displayDataSpiral(Spiral s, char* spiralName, int displayRadians)
{
	char* output;

	output = geolib_idealab::createSpiralString(s, spiralName, DATA_OUT, 1);
	printf("%s", output);
}

void geolib_idealab::displayDataBndry(Boundary b, char* boundaryName, int displayRadians)
{
	char* output;

	output = geolib_idealab::createBndryString(b, boundaryName, DATA_OUT, 1);
	printf("%s", output);
}

void geolib_idealab::displayDataComplexBndry(ComplexBoundary c, char* complexBoundaryName, int displayRadians)
{
	char* output;

	output = geolib_idealab::createComplexBndryString(c, complexBoundaryName, DATA_OUT, 1);
	printf("%s", output);
}
