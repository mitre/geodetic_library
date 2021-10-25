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

/*  \file ErrorCodes.c
 *  \brief This file defines the error functions used by geolib.
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
#include "ErrorCodes.h"

#ifdef WIN32
#define strdup _strdup
#endif

using namespace geolib_idealab;

/*
 * Apply a mask to the error code and return the result
 * inputs:
 * 		err - An error value comprised of ERRCODES Enums
 *      mask - A mask value comprised of ERRCODES Enums
 *
 * */
ErrorSet geolib_idealab::getMaskedError(ErrorSet err, long mask)
{

    long maskComplement = 0;

    //Get the bitwise complement of the mask
    maskComplement = ~mask;

    //apply the mask complement to the err, & = bitwise AND
    return (maskComplement & err);

}

/*
 * Return a list of specified exit codes that are currently stored as error codes.
 * Passing in a non-zero value for an argument will add that exit code to the mask.
 *
 * */
long geolib_idealab::getMask(int err64, int err128, int err256, int err2048, int err4096,
             int err32768, int err65536)
{

    long mask = 0;

    if (err64)
        mask |= COLLINEAR_COURSE_ERR;
    if (err128)
        mask |= CONCENTRIC_CIRCLE_ERR;
    if (err256)
        mask |= NO_INTERSECTION_ERR;
    if (err2048)
        mask |= POINT_NOT_ON_LINE_ERR;
    if (err4096)
        mask |= NO_TANGENT_ARC_ERR;
    if (err32768)
        mask |= LINE_TOO_FAR_FROM_ARC_ERR;
    if (err65536)
        mask |= CIRCLE_INSIDE_CIRCLE_ERR;

    return mask;

}

/*
 * Return a list of all the exit codes that are currently stored as error codes.
 * The returned value should be 104896 which is the sum of the flagged error code values
 *
 * */
long geolib_idealab::getMaskAll()
{

    return geolib_idealab::getMask(1, 1, 1, 1, 1, 1, 1);
}

char* geolib_idealab::formatErrorMessage(ErrorSet err)
{
    char* buff = NULL;
    char* msgList[32]; /* Can be no more than 32 error messages */
    int msgCount = 0;
    int charCount = 0;
    int i = 0;

    if (err == 0)
    {
        msgList[msgCount] = strdup("No Error\n");
        charCount += (int) strlen(msgList[msgCount++]);
    }
    if (err & INVALID_CRS_ERR)
    {
        msgList[msgCount] = strdup("Invalid Course\n");
        charCount += (int) strlen(msgList[msgCount++]);
    }
    if (err & INVALID_LAT_ERR)
    {
        msgList[msgCount] = strdup("Invalid Latitude\n");
        charCount += (int) strlen(msgList[msgCount++]);
    }
    if (err & INVALID_LON_ERR)
    {
        msgList[msgCount] = strdup("Invalid Longitude\n");
        charCount += (int) strlen(msgList[msgCount++]);
    }
    if (err & INV_NOT_CONVERGED_ERR)
    {
        msgList[msgCount] = strdup("Vincenty Inverse did not converge\n");
        charCount += (int) strlen(msgList[msgCount++]);
    }
    if (err & NO_ROOT_ERR)
    {
        msgList[msgCount] = strdup("No root found in linear solver\n");
        charCount += (int) strlen(msgList[msgCount++]);
    }
    if (err & MALLOC_ERR)
    {
        msgList[msgCount] = strdup("Memory not allocated\n");
        charCount += (int) strlen(msgList[msgCount++]);
    }
    if (err & COLLINEAR_COURSE_ERR)
    {
        msgList[msgCount] = strdup(
                "Courses are colinear--infinite intersections\n");
        charCount += (int) strlen(msgList[msgCount++]);
    }
    if (err & CONCENTRIC_CIRCLE_ERR)
    {
        msgList[msgCount] = strdup("Circles are concentric\n");
        charCount += (int) strlen(msgList[msgCount++]);
    }
    if (err & NO_INTERSECTION_ERR)
    {
        msgList[msgCount] = strdup("No intersection found\n");
        charCount += (int) strlen(msgList[msgCount++]);
    }
    if (err & NO_MEMORY_ALLOCATED_ERR)
    {
        msgList[msgCount] = strdup("No memory allocated for return value\n");
        charCount += (int) strlen(msgList[msgCount++]);
    }
    if (err & NO_PROJECTED_POINT_ERR)
    {
        msgList[msgCount] = strdup("No projected point found\n");
        charCount += (int) strlen(msgList[msgCount++]);
    }
    if (err & POINT_NOT_ON_LINE_ERR)
    {
        msgList[msgCount] = strdup("Point not on geodesic\n");
        charCount += (int) strlen(msgList[msgCount++]);
    }
    if (err & NO_TANGENT_ARC_ERR)
    {
        msgList[msgCount] = strdup("No tangent arc exists\n");
        charCount += (int) strlen(msgList[msgCount++]);
    }
    if (err & NO_SPHERICAL_SOLUTION_ERR)
    {
        msgList[msgCount] = strdup("No spherical solution found\n");
        charCount += (int) strlen(msgList[msgCount++]);
    }
//    if (err & EPS_TOO_SMALL_ERR)
//    {
//        msgList[msgCount] = strdup("Epsilon parameter too small\n");
//        charCount += (int) strlen(msgList[msgCount++]);
//    }
    if (err & LINE_TOO_FAR_FROM_ARC_ERR)
    {
        msgList[msgCount] = strdup("Line too far from arc for intersection\n");
        charCount += (int) strlen(msgList[msgCount++]);
    }
    if (err & CIRCLE_INSIDE_CIRCLE_ERR)
    {
        msgList[msgCount] = strdup(
                "One circle contained in the other--no common tangents\n");
        charCount += (int) strlen(msgList[msgCount++]);
    }
    if (err & POINT_NOT_ON_ARC_ERR)
    {
        msgList[msgCount] = strdup("Point not on arc\n");
        charCount += (int) strlen(msgList[msgCount++]);
    }
    if (err & SUBTENDED_ANGLE_OUT_OF_RANGE_ERR)
    {
        msgList[msgCount]
                = strdup("Subtended angle greater than 360 degrees\n");
        charCount += (int) strlen(msgList[msgCount++]);
    }
    if (err & TOL_TOO_SMALL_ERR)
    {
        msgList[msgCount] = strdup("Tolerance parameter too small\n");
        charCount += (int) strlen(msgList[msgCount++]);
    }
    if (err & DUPLICATE_POINTER_ERR)
    {
        msgList[msgCount] = strdup("Duplicate pointers passed as arguments\n");
        charCount += (int) strlen(msgList[msgCount++]);
    }
    if (err & RADIUS_OUT_OF_RANGE_ERR)
    {
        msgList[msgCount] = strdup("Radius out of range\n");
        charCount += (int) strlen(msgList[msgCount++]);
    }
//    if (err & AZIMUTH_OUT_OF_RANGE_ERR)
//    {
//        msgList[msgCount] = strdup("Azimuth out of range\n");
//        charCount += (int) strlen(msgList[msgCount++]);
//    }
    if (err & INVALID_TYPE_ERR)
    {
        msgList[msgCount] = strdup("Invalid Type\n");
        charCount += (int) strlen(msgList[msgCount++]);
    }
    if (err & SHAPE_NOT_DEFINED_ERR)
    {
        msgList[msgCount] = strdup("Shape is not defined\n");
        charCount += (int) strlen(msgList[msgCount++]);
    }
    if (err & INVALID_SHAPE_ERR)
    {
        msgList[msgCount] = strdup("Invalid shape\n");
        charCount += (int) strlen(msgList[msgCount++]);
    }
    if (err & FWD_NOT_CONVERGED_ERR)
    {
        msgList[msgCount] = strdup("Vincenty Forward did not converge\n");
        charCount += (int) strlen(msgList[msgCount++]);
    }
    if (err & SEC_NOT_CONVERGED_ERR)
    {
        msgList[msgCount] = strdup("Secant Method did not converge\n");
        charCount += (int) strlen(msgList[msgCount++]);
    }
    if (err & ITERATION_MAX_REACHED_ERR)
    {
        msgList[msgCount] = strdup("Maximum iteration count reached\n");
        charCount += (int) strlen(msgList[msgCount++]);
    }
    if (err & ERROR_MAX_REACHED_ERR)
    {
        msgList[msgCount] = strdup(
                "Error exceeds bound; iteration not converged\n");
        charCount += (int) strlen(msgList[msgCount++]);
    }
    if (err & ANTIPODAL_POINTS_ERR)
    {
        msgList[msgCount] = strdup(
                "Nearly antipodal points; computation not accurate\n");
        charCount += (int) strlen(msgList[msgCount++]);
    }
    if (err & UNEXPECTED_ERR)
    {
        msgList[msgCount] = strdup("Unexpected Error\n");
        charCount += (int) strlen(msgList[msgCount++]);
    }

    buff = (char*) calloc(sizeof(char), charCount + 1);

    for (i = 0; i < msgCount; i++)
    {
        buff = strcat(buff, msgList[i]);
    }

    return buff;

}

