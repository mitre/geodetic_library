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

/*! \file Constants.h
 *  
 * This file defines the constants used by geolib.
 *
 */

/*
 * A NOTE ABOUT PARAMETERS AND UNITS
 * Unless otherwise noted below, all distances are in nautical miles and
 * all angles, latitudes, longitudes, courses, and azimuths are in radians.
 * Standard conversions are:
 *   nautical miles to feet: multiply by 1852/0.3048 or use FEET_PER_NMI macro defined below
 *   radians to degrees: multiply by 180.0/pi, or use DEG_PER_RAD macro defined below.
 *
 * \author Michael Mills, Richard Snow, Stuart Bowman, Juan Amezcua, John Landrigan
 *
 */

namespace geolib_idealab {
#ifndef WGS84_PARMS
/* Defining WGS-84 Parameters */
/** */
#define FLATTENING 0.00335281066474748
/** */
#define INVERSE_FLATTENING 298.257223563
/** */
#define SEMI_MAJOR_AXIS_FEET 20925646.3254593
/** */
#define SEMI_MAJOR_AXIS_METERS 6378137.0
/** */
#define SEMI_MAJOR_AXIS_NMI 3443.918466522678

#define SEMI_MINOR_AXIS_NMI 3432.3716599596  /*((1.0 - FLATTENING) * SEMI_MAJOR_AXIS_NMI) */

#define SEMI_MINOR_AXIS_METERS 6356752.314245179497564

//#define DEG2RAD 0.01745329251994329

#define ECCENTRICITY_SQ 0.0066943799901413169961

#define SECOND_ECCENTRICITY_SQ 0.0067394967422764349548

#define THIRD_FLATTENING 0.0016792203863837046951

/** */
#define FEET_PER_NMI 6076.115485564304           /* = (1852 m/nmi)/(0.3048 m/feet) (both exact) */
/** */
#define SPHERE_RADIUS_NMI 3438.140221487929
#define FLAT_TOL 5.0e-4                          /* Approximately 1 meter in nmi */
#define SPHERICAL_AZIMUTH_CUTOFF_DISTANCE 5.0e-4 /* Approximately 1 meter in nmi */
#define INTERNAL_ZERO 1e-15                      /* Not a scientific number, but something near zero is needed in a few cases */
#define ONE_NANOMETER_NMI 5.399568e-013          /* Sometimes we need a small value that can't be supplied by the library caller */
#define SMALL_DIST_THRESHOLD 1.0/1852.0          /* Calculations between points closer than this distance
                                                  * will use single-step projection algorithms */

#define ANTIPODAL_TOL 3.12413936106985 /* (179.0 * M_PI / 180.0)  Maximum allowable dihedral angle between two points */

/* Set sphere radius to default value = sqrt(a*b) (a macro in Geolib.h)*/
#define SPHERE_RADIUS SPHERE_RADIUS_NMI

#define MAX_ELLIPSOIDAL_ARC_RADIUS_NMI (M_PI_2 * (1.0 - FLATTENING) * SEMI_MAJOR_AXIS_NMI)

#define WGS84_PARMS
#endif /* WGS84_PARMS */

/** */
#define DEFAULT_EPS 1.0e-20
// Default is for iterateInverse and iterateDirect to use eps
// no larger than 1e-20
// To allow larger eps, pass -D KEEP_USER_EPS to compiler OR
// uncomment next line:
//#define KEEP_USER_EPS

/* Recursion limit */
#ifndef MAX_RECURSIONS
#define MAX_RECURSIONS 10
#endif

/* Iteration limit */
#ifndef MAX_ITERATIONS
/** */
#define MAX_ITERATIONS 100
#endif

/* Maximum distance error for refining algorithms */
#ifndef MAX_DISTANCE_ERROR
#define MAX_DISTANCE_ERROR 1000.0  //nautical miles
#endif

/* Vincenty iteration limit */
#ifndef MAX_V_ITERATIONS
#define MAX_V_ITERATIONS 100
#endif

#ifndef M_PI
/** */
#define M_PI 3.14159265358979323846
#endif

#ifndef M_PI_2
/** */
#define M_PI_2 1.57079632679489661923
#endif

#ifndef M_PI_4
/** */
#define M_PI_4 0.78539816339744830962
#endif

#ifndef M_2PI
/** */
#define M_2PI 6.28318530717958647692
#endif

#ifndef RAD_PER_DEG
/** */
#define RAD_PER_DEG 57.2957795130823
#endif
} // namespace