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

/*  \file Vector.c
 *  \brief This file defines the arc functions used by geolib.
 *  \author Michael Mills, Richard Snow, Stuart Bowman, Juan Amezcua, John Landrigan
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
#include "Geolib.h"

namespace geolib_idealab {

Vector geodeticToECEF(LLPoint geo)
{

    double f = FLATTENING;
    double ee = f * (2.0 - f);
    double a = SEMI_MAJOR_AXIS_NMI;
    double lat = geo.latitude;
    double lon = geo.longitude;
    double sinLat = sin(lat);
    double N = a / sqrt(1.0 - ee * sinLat * sinLat);

    Vector v;

    /* Height is assumed to zero (on the ellipsoid) */

    v.x = N * cos(lat) * cos(lon);
    v.y = N * cos(lat) * sin(lon);
    v.z = N * (1.0 - ee) * sin(lat);

    return v;

}

Vector cross(Vector v1, Vector v2)
{

    Vector prod = { 0.0, 0.0, 0.0 };

    prod.x = v1.y * v2.z - v2.y * v1.z;
    prod.y = v1.z * v2.x - v1.x * v2.z;
    prod.z = v1.x * v2.y - v1.y * v2.x;

    return prod;

}

double dot(Vector v1, Vector v2)
{

    return (double) (v1.x * v2.x + v1.y * v2.y + v1.z * v2.z);

}

double norm(Vector v)
{

    return sqrt(dot(v, v));

}

Vector mapToUnitSphere(LLPoint p)
{

    Vector v;

    v.x = cos(p.latitude) * cos(p.longitude);
    v.y = cos(p.latitude) * sin(p.longitude);
    v.z = sin(p.latitude);

    return v;

}

void normalize(Vector* v)
{

    double length = norm(*v);

    v->x /= length;
    v->y /= length;
    v->z /= length;

}

void scalarMultiply(Vector* v, double scal)
{

    v->x *= scal;
    v->y *= scal;
    v->z *= scal;

}

Vector vectorAdd(Vector v1, Vector v2)
{

    Vector sum;

    sum.x = v1.x + v2.x;
    sum.y = v1.y + v2.y;
    sum.z = v1.z + v2.z;

    return sum;

}

Vector vectorSubtract(Vector v1, Vector v2)
{

    Vector sum;

    sum.x = v1.x - v2.x;
    sum.y = v1.y - v2.y;
    sum.z = v1.z - v2.z;

    return sum;

}

LLPoint mapVectorToSphere(Vector v)
{

    LLPoint p;
    normalize(&v);
    p.latitude = asin(v.z);
    p.longitude = atan2(v.y, v.x);

    return p;

}

} //namespace
