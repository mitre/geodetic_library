/*
 * Copyright 2007-2011 The MITRE Corporation.  All Rights reserved.
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

/** \file ErrorCodes.h
 *  \brief This file defines the error codes returned by Geolib functions.
 *  \author Michael Mills, Richard Snow, Stuart Bowman, Juan Amezcua, John Landrigan
 */

/*
 * A NOTE ABOUT PARAMETERS AND UNITS
 * Unless otherwise noted below, all distances are in nautical miles and
 * all angles, latitudes, longitudes, courses, and azimuths are in radians.
 * Standard conversions are:
 *   nautical miles to feet: multiply by 1852/0.3048 or use FEET_PER_NMI macro defined below
 *   radians to degrees: multiply by 180.0/pi, or use DEG_PER_RAD macro defined below
 */

#pragma once

namespace geolib_idealab {
#ifndef ERRORSET
/** @typedef ErrorSet
 * 
 * A long integer containing error code information.  Set to value of ErrorCode enum and returned by most geolib functions. Equivalent to type "long".
 */
typedef long ErrorSet;  /* Alias to make it easier to extend error handling later */
#define ERRORSET
#endif


/* Error codes */
/* These error codes are orthogonal: for any two error codes, A, B, A&B = 0. */
#ifndef ERRCODES
/** @enum ErrorCodes
 *
 * The ErrorCodes enumeration assigns a long integer value to
 * each error code that may be returned by a Geolib function.  The error
 * codes are orthogonal (i.e., for any two error codes, \f$e_1\f$ and
 * \f$e_2\f$, \f$e_1\ \&\ e_2 = 0\f$).  Therefore, multiple errors can
 * be stored in a single long integer value.
 */
enum ErrorCodes {
/** Indicates successful execution. */
  SUCCESS = 0x0,
/** Indicates that an azimuth value was out of range (e.g., desired course from pole other than \f$\pm\pi\f$). */
  INVALID_CRS_ERR = 0x1,     /*  = 1*1     = 1     = 0000 0000 0000 0000 0001 */
/** Indicates a latitude not in the range \f$[-\pi/2,\pi/2]\f$ */
  INVALID_LAT_ERR = 0x2,     /*  = 2*1     = 2     = 0000 0000 0000 0000 0010 */
/** Indicates a longitude not in the range \f$[-\pi,\pi]\f$ */
  INVALID_LON_ERR = 0x4,     /*  = 4*1     = 4     = 0000 0000 0000 0000 0100 */
/** Indicates that the inverse function did not converge to a solution MAX_ITERATION_COUNT iterations */
  INV_NOT_CONVERGED_ERR = 0x8,     /*  = 8*1     = 8     = 0000 0000 0000 0000 1000 */
/** Indicates that root could not be found to linear approximation of the error function */
  NO_ROOT_ERR = 0x10,    /*  = 1*16    = 16    = 0000 0000 0000 0001 0000 */
/** Indicates that memory could not be allocated.*/
  MALLOC_ERR = 0x20,    /*  = 2*16    = 32    = 0000 0000 0000 0010 0000 */
/** Indicates that two geodesics lie on top of each other and do not have discrete intersections.  Used as a status code. */
  COLLINEAR_COURSE_ERR =
  0x40,    /*  = 4*16    = 64    = 0000 0000 0000 0100 0000 */ //This will be a status indicator in future versions.
/** Indicates that two arcs or circles either do not intersect or are identical. Used as a status code. */
  CONCENTRIC_CIRCLE_ERR =
  0x80,    /*  = 8*16    = 128   = 0000 0000 0000 1000 0000 */ //This will be a status indicator in future versions.
/** Status code indicates that no intersection was found in the case that no intersection point gets returned. */
  NO_INTERSECTION_ERR =
  0x100,   /*  = 1*256   = 256   = 0000 0000 0001 0000 0000 */ //This will be a status indicator in future versions.
/** Indicates that a NULL pointer was passed for a required reference */
  NO_MEMORY_ALLOCATED_ERR = 0x200,   /*  = 2*256   = 512   = 0000 0000 0010 0000 0000 */
/** Indicates that an intermediate calculation point based on projecting to a geodesic could not be found */
  NO_PROJECTED_POINT_ERR = 0x400,   /*  = 4*256   = 1024  = 0000 0000 0100 0000 0000 */
/** Status code indicates that a point was not on the expected geodesic or locus */
  POINT_NOT_ON_LINE_ERR =
  0x800,   /*  = 8*256   = 2048  = 0000 0000 1000 0000 0000 */ //This will be a status indicator in future versions.
/** Status code indicates that no tangent arc could be found. */
  NO_TANGENT_ARC_ERR =
  0x1000,  /*  = 1*4096  = 4096  = 0000 0001 0000 0000 0000 */ //This will be a status indicator in future versions.
/** Indicates that no spherical solution was found so iteration could not begin. */
  NO_SPHERICAL_SOLUTION_ERR = 0x2000,  /*  = 2*4096  = 8192  = 0000 0010 0000 0000 0000 */
/** Status code indicates that no geodesic/locus-arc intersection could be found because the geodesic/locus is too far from the arc. */
  LINE_TOO_FAR_FROM_ARC_ERR =
  0x8000,  /*  = 8*4096  = 32768 = 0000 1000 0000 0000 0000 */ //This will be a status indicator in future versions.
/** Status code indicates that no arc-arc intersection was found because one arc lies entirely inside the other. */
  CIRCLE_INSIDE_CIRCLE_ERR =
  0x10000, /*  = 1*65536 = 65536 = 0001 0000 0000 0000 0000 */ //This will be a status indicator in future versions.
/** Indicates that start/end points used to define an arc are not equidistant from arc center. */
  POINT_NOT_ON_ARC_ERR = 0x20000,
/** Indicates that arc's subtended angle does not meet algorithm requirement. */
  SUBTENDED_ANGLE_OUT_OF_RANGE_ERR = 0x40000,
/** Indicates that the requested tolerance cannot be met due to large requested Vincenty algorithm precision. */
  TOL_TOO_SMALL_ERR = 0x80000,
/** Indicates that two or more passed reference share memory */
  DUPLICATE_POINTER_ERR = 0x100000,
/** Indicates that given radius does not meet algorithm requirement. */
  RADIUS_OUT_OF_RANGE_ERR = 0x200000,
/** Indicates that given azimuth does not meet algorithm requirement. */
  AZIMUTH_OUT_OF_RANGE_ERR = 0x400000,
/** Indicates that the incorrect variable type has been passed. */
  INVALID_TYPE_ERR = 0x800000,
/** Indicates that a shape object has not been defined. */
  SHAPE_NOT_DEFINED_ERR = 0x1000000,
/** Indicates that an invaled shape type has been passed. */
  INVALID_SHAPE_ERR = 0x2000000,
/** Indicates that the direct algorithm failed to converge.  */
  FWD_NOT_CONVERGED_ERR = 0x4000000,
/** Indicates that the secant method failed to converge. */
  SEC_NOT_CONVERGED_ERR = 0x8000000,
/** Indicates that a method has looped more than the allowed iteration count. */
  ITERATION_MAX_REACHED_ERR = 0x10000000,
/** Indicates that an error value has grown larger than allowed. */
  ERROR_MAX_REACHED_ERR = 0x20000000,
/** Indicates that two LLPoints are on opposite sides of the earth ellipsoid. */
  ANTIPODAL_POINTS_ERR = 0x40000000,
/** Indicates that an unknown error has occurred. */
  UNEXPECTED_ERR = 0x80000000
};
/* Cannot define larger err code in enum */
#define ERRCODES
#endif





/* The following functions are used to manipulate error codes */

/**
 * Apply a mask to remove certain error codes from return result
 *
 * @param err An error value comprised of ERRCODES enums, some of which shall be masked (i.e., removed)
 * @param mask A long integer value that is the result of applying AND operator to all 
 *             error code values to be masked.
 * @returns ErrorSet with masked ERRORCODES bits set to 0.
 *
 */
ErrorSet getMaskedError(ErrorSet err, long mask);

/**
 * Return a mask value to remove specific error codes from return value. These error codes are sometimes returned to 
 * indicate computation status that is not necessarily an error (i.e., no intersection may be legitimate return value). 
 * @param err64 Integer,set to 1 to mask COLLINEAR_COURSE_ERR;
 * @param err256 Integer, set to 1 to mask CONCENTRIC_CIRCLE_ERR;
 * @param err256 Integer, set to 1 to mask NO_INTERSECTION_ERR;
 * @param err204 Integer, set to 1 to mask POINT_NOT_ON_LINE_ERR;
 * @param err4096 Integer,set to 1 to mask NO_TANGENT_ARC_ERR;
 * @param err32768 Integer, set to 1 to mask LINE_TOO_FAR_FROM_ARC_ERR;
 * @param err65536 Integer, set to 1 to mask CIRCLE_INSIDE_CIRCLE_ERR;
 * @returns long integer with appropriate bits set to 1 to mask desired error codes. 
 */
long getMask(int err64, int err128, int err256, int err2048, int err4096, int err32768, int err65536);

/**
 * Return a mask value that will remove all error codes described in getMask() description from error value. Equivalent to 
 * calling getMask(1,1,1,1,1,1,1).
 *
 */
long getMaskAll();

/** Converts ErrorCode parameter into string Provides text output corresponding to all errors stored in err parameter
 * @param err ErrorCode (i.e., long integer) containing error code.
 * @returns variable length string containing error code description.
*/
char *formatErrorMessage(ErrorSet err);
} // namespace


