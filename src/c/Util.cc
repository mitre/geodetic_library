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
#include "Util.h"

using namespace geolib_idealab;

/*
 * NAME: Util.c
 *
 * DESCRIPTION:
 * 		This source file contains utility functions that are used for the
 * 		geolib project.  This source does not contain utility functions which are
 * 		strictly used for testing functions in the geolib project.
 *
 */

int geolib_idealab::sgn(double x)
{
    return x < 0.0 ? -1 : x > 0.0;
}

double geolib_idealab::reciprocal(double x)
{
    return modpos(x + M_PI, M_2PI);
}

double geolib_idealab::computeSubtendedAngle(double startCrs, double endCrs,
                             ArcDirection orient)
{

    double alpha;
    double temp;

    if (orient != COUNTERCLOCKWISE)
    {
        /* always use counter-clockwise orientaion */
        temp = startCrs;
        startCrs = endCrs;
        endCrs = temp;
    }

    if (startCrs > endCrs)
    {
        alpha = startCrs - endCrs;
    }
    else
    {
        alpha = M_2PI - (endCrs - startCrs);
    }

    alpha = ((double) orient) * alpha;

    return alpha;

}

/* Given two angles, representing courses, return the difference in the
 * smallest angle between them (i.e, in the range [0,PI]).
 */
ErrorSet geolib_idealab::minSubtendedAngle(double crs1, double crs2, double* angle)
{

    ErrorSet err = 0;

    *angle = crs1 - crs2;

    if (*angle > M_PI)
    {
        /* crs1 is to left of branch cut, crs2 is to right */
        *angle = M_2PI - *angle;
    }
    else if (*angle < -M_PI)
    {
        /* crs2 is to left of branch cut, crs1 is to right */
        *angle = M_2PI + *angle;
    }

    /* want only positive angles */
    *angle = fabs(*angle);

    return err;

}

double geolib_idealab::findN(LLPoint p)
{
    double f = FLATTENING;
    double ee = f * (2.0 - f);
    double lat = p.latitude;
    double sinLat = sin(lat);
    return SEMI_MAJOR_AXIS_NMI / sqrt(1.0 - ee * sinLat * sinLat);
}

double geolib_idealab::findM(LLPoint p)
{
    double f = FLATTENING;
    double ee = f * (2.0 - f);
    double lat = p.latitude;
    double sinLat = sin(lat);
    return SEMI_MAJOR_AXIS_NMI * (1.0 - ee) / pow(1.0 - ee * sinLat * sinLat,
            1.5);
}

double geolib_idealab::findRootSecantMethod(double* x, double* y, ErrorSet* err)
{

    if (x[0] == x[1])
    {
        return x[0];
    }
    if (y[0] == 0.0)
	{
		return x[0];
	}
    if (y[1] == 0.0)
	{
		return x[1];
	}
    if (y[0] == y[1])
    {
        return 0.5 * (x[0] + x[1]);
    }

    return -y[0] * (x[1] - x[0]) / (y[1] - y[0]) + x[0];

}

/*
 * NAME: geodeticLat
 *
 * DESCRIPTION:
 * 		Convert geocentric latitude to geodetic latitude
 *
 * INPUT(Type):
 * 		lat(double) = The geocentric latitude to be converted in radians.
 *
 * OUTPUT(Return Type):
 * 		geodeticLat(double) = The geodetic latitude in radians for the
 * 			given input geocentric latitude.
 *
 */
double geolib_idealab::geodeticLat(double lat)
{
    /* Angular eccentricity */
    double cosoe = 1.0 - FLATTENING;
    return atan(tan(lat) / cosoe / cosoe);
    return lat;
}

/* Convert geodetic latitude to geocentric latitude */
double geolib_idealab::geocentricLat(double lat)
{
    /* Angular eccentricity */
    double cosoe = 1.0 - FLATTENING;
    return atan(cosoe * cosoe * tan(lat));
    return lat;
}

LLPoint geolib_idealab::geodeticToGeocentric(LLPoint pt)
{
    LLPoint newPt = pt;
    newPt.latitude = geolib_idealab::geocentricLat(newPt.latitude);
    return newPt;

    return pt;
}

LLPoint geolib_idealab::geocentricToGeodetic(LLPoint pt)
{
    LLPoint newPt = pt;
    newPt.latitude = geolib_idealab::geodeticLat(newPt.latitude);
    return newPt;
    return pt;
}

/*
 * NAME: dms2rad
 *
 * DESCRIPTION:
 * 		Convert decimal degrees minutes seconds (DMS) to radians
 *
 * INPUT(Type):
 * 		d(int) = integer degrees
 * 		m(int) = integer minutes
 * 		s(double) = decimal seconds
 *
 * OUTPUT(Return Type):
 * 		dms2rad(double) = The decimal degrees minutes seconds in
 * 			radians.
 *
 */
double geolib_idealab::dms2rad(int d, int m, double s)
{

    if (d < 0 && m > 0)
        m = -m;
    if (d < 0 && s > 0)
        s = -s;

    return ((double) d + ( (double) m + s/60.0)/60.0)*M_PI/180.0;

}

double geolib_idealab::dms2radHemi(double d, double m, double s, char* hemi)
{
    double rad;

    rad = d + (m + s/60.0)/60.0;
    rad *= M_PI/180.0;

    if (strchr(hemi,'S') || strchr(hemi,'W') )
        rad *= -1.0;

    return rad;

}

/*
 * NAME: rad2dms
 *
 * DESCRIPTION:
 * 		Convert radians into decimal degrees minutes seconds (DMS)
 *
 * INPUT(Type):
 * 		crs(double) = course in radians
 * 		d(*int) = returned integer degrees
 * 		m(*int) = returned integer minutes
 * 		s(*double) = returned decimal seconds
 *
 * OUTPUT(Return Type):
 * 		d(*int) = returned integer degrees
 * 		m(*int) = returned integer minutes
 * 		s(*double) = returned decimal seconds
 *
 */
void geolib_idealab::rad2dms (double crs, int *d, int *m, double *s)
{

    double tempm;
    double tempd;

    crs = crs*180.0/M_PI;
    //  printf("crs = %f\n",crs);
    if (crs >= 0.0)
    {
        tempd = floor(crs);
    }
    else
    {
        tempd= ceil(crs);
    }

    tempm = fabs(crs - tempd)*60.0;

    *s= tempm - floor(tempm);

    *d = (int) tempd;
    *m = (int) floor(tempm);
    *s = *s*60;

    if (crs < 0.0)
    {
        *m *= -1;
        *s *= -1;
    }

}

//TODO Need to verify wherever this is called to confirm the test data will support
//using rad2dms in place of crsintrad2dms.  If so then rad2dms should be used
//and this function should go away - JAMEZCUA
void geolib_idealab::crsintrad2dms (double crs, int *d, int *m, double *s)
{

    double tempm;
    double tempd;

    crs = crs*180.0/M_PI;
    //  printf("crs = %f\n",crs);
    if (crs >= 0.0)
    {
        tempd = floor(crs);
    }
    else
    {
        tempd= ceil(crs);
    }

    //  printf("tempd = %f\n",tempd);

    tempm = fabs(crs - tempd)*60.0;

    //  printf("floor(tempm) = %f\n",floor(tempm));

    *s= tempm - floor(tempm);
    *d = (int) tempd;
    *m = (int) floor(tempm);
    *s = *s*60;

}

//TODO This actually looks like a more correct version of how the sgn function
//in libWGS84.c should operate since this version handles zero, need to verify
//the version in libWGS84.c handles zero properly since zero is neither
//positive nor negative - JAmezcua
double geolib_idealab::crsintsgn(double x)
{
    return (x > 0.0) ? 1 : ((x < 0) ? -1 : 0);
}

double geolib_idealab::minimum(double value1, double value2)
{
	if (value1 <= value2)
		return value1;
	else
		return value2;
}

double geolib_idealab::maximum(double value1, double value2)
{
	if (value1 >= value2)
		return value1;
	else
		return value2;
}

ErrorSet geolib_idealab::findSetMaxAndMin(double values[], int size, double* max, double* min)
{
	ErrorSet err = 0;
	double maxVal;
	double minVal;
	double temp, element;
	int i = 0;

    if (NULL == values)
    {
    	err |= NO_MEMORY_ALLOCATED_ERR;
        return err;
    }

	maxVal = values[0];
	minVal = values[0];

	for (i = 1; i < size; i++)
	{
		element = values[i];

		temp = geolib_idealab::maximum(maxVal, element);
		maxVal = temp;

		temp = geolib_idealab::minimum(minVal, element);
		minVal = temp;
	}

	if (NULL != max)
		(*max) = maxVal;

	if (NULL != min)
		(*min) = minVal;

	return err;
}

double geolib_idealab::modlat(double x)
{
    x = (x >= 0 ? x : x + M_2PI);
    return modpos(x + M_PI_2, M_2PI) - M_PI_2; /* M_PI_2 = PI/2, from math.h */
}

double geolib_idealab::modlon(double x)
{
    x = (x >= 0 ? x : x + M_2PI);
    return modpos(x + M_PI, M_2PI) - M_PI;
}

double geolib_idealab::modcrs(double crs)
{
    /* Map [-pi,pi] to [0,2*pi] */
    return modpos(crs + M_2PI, M_2PI);
}

double geolib_idealab::modpos(double x, double y)
{

    /* returns positive remainder of x/y. */

    x = fmod(x, y);
    if (x < 0.0)
        x = x + y;
    return x;

}

