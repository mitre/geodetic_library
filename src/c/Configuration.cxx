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

#include<stdlib.h>  // needed for malloc() and free()
#include<stdio.h>
#include<math.h>
#if REPLACE_WITH_AMDLIBM
#include "amdlibm.h"
#endif

/* Header Code--Configuration.h
 *
 * TODO: Move to Configuration.h when implementing Configuration capability
 */

#include "Util.h"
#ifndef CONFIGURATION_H_
#define CONFIGURATION_H_

#endif /* CONFIGURATION_H_ */

typedef enum {
    WGS84        = 0,
    GRS80        = 1
    } Ellipsoid;

typedef struct _datum Datum;

/** @struct Configuration Configuration.h "include/Configuration.h"
 * @brief Defines numerical values used by most algorithms.
 *
 *  @param tol The tolerance (in NM) to which iterative solutions will be forced to converge
 *  @param eps (optional) The longitude error (in radians) allowed for the Vincently algorithms
 *  @param semiMajorAxis The ellipsoid's equatorial radius (in NM)
 *  @param flattening The ratio of the difference between the equatorial and polar radii to
 * the equatorial radius (dimensionless)
 *  @param lastErr Long parameter that stores the most recent error value
 */
#ifndef CONFIGURATION_STRUCT
typedef struct {
    double tol;
    double eps; // Leave this out?
    Datum* datum;
    double internalZero;
    double smallDistThreshold;
    int iterationCount;
    ErrorSet lastErr;
} Configuration;
#define CONFIGURATION_STRUCT
#endif

/* Functions to manipulate Configuration */
const Configuration* getDefaultConfigurationPtr();
ErrorSet newConfiguration(Configuration* cfg, double tol, double eps, Ellipsoid ellipsoid,
                                double smallDistThreshold);
ErrorSet removeConfiguration(Configuration** cfg);
ErrorSet setConfigurationTolerance(Configuration* cfg, double tol);
double getConfigurationTolerance(Configuration* cfg);
ErrorSet setConfigurationEpsilon(Configuration* cfg, double eps);
double getConfigurationEpsilon(Configuration* cfg);
ErrorSet setConfigurationSemiMajorAxis(Configuration* cfg, double sma);
double getConfigurationSemiMajorAxis(Configuration* cfg);
ErrorSet setConfigurationFlattening(Configuration* cfg, double f);
double getConfigurationFlattening(Configuration* cfg);
ErrorSet setConfigurationSphereRadius(Configuration* cfg, double r);
double getConfigurationSphereRadius(Configuration* cfg);
ErrorSet setConfigurationSmallDistThreshold(Configuration* cfg, double thr);
double getConfigurationSmallDistThreshold(Configuration* cfg, double thr);
ErrorSet setConfigurationIterationCount(Configuration* cfg, int n);
double getConfigurationIterationCount(Configuration* cfg);
ErrorSet setConfigurationLastError(Configuration* cfg, ErrorSet err);
ErrorSet getConfigurationLastError(Configuration* cfg);
ErrorSet displayConfiguration(Configuration* cfg);
ErrorSet setDatum(Configuration* cfg, Ellipsoid ellipsoidVal);

/* End Configuration.h */

#define DATUM_NAME_MAX_LENGTH 32
#define MAX_EPS_TOL_RATIO 1.0e-4
/*
 * NAME: Configuration.c
 *
 * DESCRIPTION:
 *      This source file defines a default configuration struct to hold
 *      basic constant values in one place.  The function getDefaultConfig()
 *      returns a pointer to this struct when called.
 *
 */

/** @struct Datum Configuration.h "include/Configuration.h"
 * @brief Stores semi-major axis and flattening to define ellipsoids
 * @param semiMajorAxis The semi-major axis of the ellipsoid in NM
 * @param flattening The flattening ratio of the ellipsoid, equal to \f$f=\frac{a-b}{a}\f$ where \f$a\f$ and
 * \f$b\f$ are the semi-major and semi-minor axes, respectively.
 * @param sphereRadius The radius (in NM) of the sphere that fits the ellipsoid the best, equal to
 * \f$\sqrt{ab} = a\sqrt{(1-f)}.\f$
 * @param label A string containing the common name of the ellipsoid (maximum 32 characters)
 */
#ifndef DATUM_STRUCT
struct _datum {
    char   label[DATUM_NAME_MAX_LENGTH];
    double semiMajorAxis;
    double flattening;
    double sphereRadius;
};
#define DATUM_STRUCT
#endif

/* Ellipsoid definitions.  These are static and will not change.  Additional elements will be added
 * as needed.  The Ellipsoid enum indexes this array and must be expanded as
 */
static Datum _datumList[2] = {
                               {"WGS84", SEMI_MAJOR_AXIS_NMI, FLATTENING, SEMI_MAJOR_AXIS_NMI }, // TODO fix non constant *sqrt(1.0-FLATTENING)},
                               {"GRS80", SEMI_MAJOR_AXIS_NMI, 0.00335281068118225, SEMI_MAJOR_AXIS_NMI}, // TODO *sqrt(1.0-0.00335281068118225)}
                             };


static Configuration DEFAULT_CONFIG = { 1.37e-9, 1.0e-20, &(_datumList[WGS84]), INTERNAL_ZERO,
                                       SMALL_DIST_THRESHOLD, 0, SUCCESS };

/* User won't be able to modify default configuration using this pointer */
const Configuration* getDefaultConfigurationPtr()
{
    return &DEFAULT_CONFIG;
}

ErrorSet newConfiguration(Configuration* cfg, double tol, double eps,
                          Ellipsoid ellipsoid, double smallDistThreshold)
{

    ErrorSet err = SUCCESS;

    if ((!cfg) || (sizeof(*cfg) != sizeof(Configuration)))
    {
        err |= MALLOC_ERR;
        return err;
    }

    setConfigurationTolerance(cfg, tol);
    setConfigurationEpsilon(cfg, eps);
    setConfigurationSmallDistThreshold(cfg, smallDistThreshold);
    setConfigurationLastError(cfg, SUCCESS);
    setConfigurationIterationCount(cfg, 0);
    setDatum(cfg,ellipsoid);

    return err;

}

ErrorSet setConfigurationTolerance(Configuration* cfg, double tol)
{
    if (cfg)
    {
        cfg->tol = tol;
        if (cfg->eps > MAX_EPS_TOL_RATIO*tol)
            cfg->eps = MAX_EPS_TOL_RATIO*tol;
        return SUCCESS;
    }
    else
        return NO_MEMORY_ALLOCATED_ERR;
}
double getConfigurationTolerance(Configuration* cfg)
{
    if (cfg)
        return cfg->tol;
    else
        return 0.0;
}

ErrorSet setConfigurationEpsilon(Configuration* cfg, double eps)
{
    if (cfg)
    {
        if (eps < MAX_EPS_TOL_RATIO*cfg->tol)
            cfg->eps = eps;
        else
            cfg->eps = MAX_EPS_TOL_RATIO*cfg->tol;
        return SUCCESS;
    }
    else
        return NO_MEMORY_ALLOCATED_ERR;

}
double getConfigurationEpsilon(Configuration* cfg)
{
    return cfg->eps;
}

//ErrorSet setConfigurationSemiMajorAxis(Configuration* cfg, double sma)
//{
//    if (cfg)
//    {
//        cfg->semiMajorAxis = sma;
//        return SUCCESS;
//    }
//    else
//        return NO_MEMORY_ALLOCATED_ERR;
//}
double getConfigurationSemiMajorAxis(Configuration* cfg)
{
    return cfg->datum->semiMajorAxis;
}

//ErrorSet setConfigurationFlattening(Configuration* cfg, double f)
//{
//    if (cfg)
//    {
//        cfg->flattening = f;
//        return SUCCESS;
//    }
//    else
//        return NO_MEMORY_ALLOCATED_ERR;
//}

double getConfigurationFlattening(Configuration* cfg)
{
    return cfg->datum->flattening;
}

//ErrorSet setConfigurationSphereRadius(Configuration* cfg, double r)
//{
//    if (cfg)
//    {
//        cfg->sphereRadius = r;
//        return SUCCESS;
//    }
//    else
//        return NO_MEMORY_ALLOCATED_ERR;
//
//}
double getConfigurationSphereRadius(Configuration* cfg)
{
    return cfg->datum->sphereRadius;
}

ErrorSet setConfigurationSmallDistThreshold(Configuration* cfg, double thr)
{
    if (cfg)
    {
        cfg->smallDistThreshold = thr;
        return SUCCESS;
    }
    else
        return NO_MEMORY_ALLOCATED_ERR;
}

double getConfigurationSmallDistThreshold(Configuration* cfg, double thr)
{
    return cfg->smallDistThreshold;
}

ErrorSet setConfigurationIterationCount(Configuration* cfg, int n)
{
    if (cfg)
    {
        cfg->iterationCount = n;
        return SUCCESS;
    }
    else
        return NO_MEMORY_ALLOCATED_ERR;
}
double getConfigurationIterationCount(Configuration* cfg)
{
    return cfg->iterationCount;
}

ErrorSet setConfigurationLastError(Configuration* cfg, ErrorSet err)
{
    if (cfg)
    {
        cfg->lastErr = err;
        return SUCCESS;
    }
    else
        return NO_MEMORY_ALLOCATED_ERR;
}

ErrorSet getConfigurationLastError(Configuration* cfg)
{
    return cfg->lastErr;
}

ErrorSet displayConfiguration(Configuration* cfg)
{
    if (cfg)
    {

        printf("Tolerance: %.15e NM\n", cfg->tol);
        printf("Epsilon:   %.15e radians\n", cfg->eps);
        printf("Datum name:     %s\n", cfg->datum->label);
        printf("Semimajor Axis: %.15e NM\n", cfg->datum->semiMajorAxis);
        printf("Flattening:     %.15e\n", cfg->datum->flattening);
        printf("Sphere Radius:  %.15e\n", cfg->datum->sphereRadius);
        printf("Small Distance Threshold: %.15e NM\n", cfg->smallDistThreshold);
        printf("Last Iteration Count: %03d\n", cfg->iterationCount);
        printf("Last Error: %#x = %s\n", (long) cfg->lastErr,
               formatErrorMessage(cfg->lastErr));
        return SUCCESS;
    }
    else
        return NO_MEMORY_ALLOCATED_ERR;

}

ErrorSet setDatum(Configuration* cfg, Ellipsoid ellipsoidVal)
{
    /* ellipsoid_val must be one of ellipsoid definitions enumerated in libWGS84.h */
    if (cfg)
    {
        cfg->datum = &(_datumList[ellipsoidVal]);
        return SUCCESS;
    }
    else
        return NO_MEMORY_ALLOCATED_ERR;

}

/* Get Datum pointer, but don't allow this pointer to be used to modify datum
 * parameters (due to const modifier) */
const Datum* lookupDatumPtr(Ellipsoid ellipsoidVal)
{

    return &(_datumList[ellipsoidVal]);

}
