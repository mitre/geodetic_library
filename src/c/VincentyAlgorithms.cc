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

/* \file Geodesic.c
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

using namespace geolib_idealab;

static void smallDistInverse(LLPoint p1, LLPoint p2, double* az, double* dist)
{
    double dlon = p2.longitude - p1.longitude;
    double dlat = p2.latitude - p1.latitude;
    double N = findN(p1);
    double M = findM(p1);
    double coslat = cos(p1.latitude);
    //    double tanAlpha = (N / M) * (cos(p1.latitude) * dlon / dlat);

    //    NdivM = coslat*coslat/BdivA/BdivA + sinlat*sinlat;
    //    tanAlpha = NdivM*(cosl((long double) p1.latitude)*dlambda/dphi);
    //    printf("tan(az') = %.16e\n",(double)tanAlpha);
    if (NULL != az)
    {

        *az = atan2(N * coslat * dlon, M * dlat);

        if (*az < 0)
        {
            *az += M_2PI;
        }
    }
    if (NULL != dist)
    {
        *dist = sqrt(pow(M * dlat, 2.0) + pow(N * coslat * dlon, 2.0));
    }

}

static LLPoint smallDistDirect(LLPoint p, double az, double ds)
{
    double N = findN(p);
    double M = findM(p);
    double dlat;
    double dlon;
	LLPoint q;

	/* Check for start at pole */
	if(fabs(p.latitude) == M_PI_2){
		dlat = (p.latitude > 0) ? -ds / M: ds / M;
		q.latitude = p.latitude + dlat;
		q.longitude = (p.latitude > 0) ? modlon((p.longitude + M_PI) - az): modlon((p.longitude) + az);
	} else {
		dlat = cos(az) * ds / M;
		dlon = sin(az) * ds / cos(p.latitude) / N;
		q.latitude = p.latitude + dlat;
		q.longitude = p.longitude + dlon;
	}

	/* Check for crossing over the pole */
	if(fabs(q.latitude) > M_PI_2){
		q.latitude = (q.latitude > 0) ? M_PI - q.latitude : -M_PI - q.latitude ;
		q.longitude = modlon(q.longitude + M_PI);
	}

    return q;
}

/******************************************************************************/
/*
 * Do basic iteration to compute inverse (distance/course from origin/dest)
 *
 * This function is only to be called by other functions in this library.
 *
 */
static ErrorSet iterateInverse(LLPoint origin, LLPoint dest, double eps,
                           double* tu1, double *tu2, double *sy, double* e,
                           double* cy, double* cz, double* y, double* c,
                           double* c2a, double* cu1, double* b1, double* su1,
                           double* cu2, double *x)
{
    /*
     EPS = 1e-9 gives results that exactly match Ed Williams' Javascript calculator
     See: <http://williams.best.vwh.net/gccalc.htm>

     In fact, this is a modified implementation of Ed Williams' code.

     The modifications made were in the manner in which the algorithm was coded.
     It was updated so as to minimize the impact of precision loss due to memory
     overflow and other rounding errors that may result when a mathematical
     operation is carried out using binary operations.

     Compare the general algorithm with Javascript at above site.

     ************* IMPORTANT NOTE FOR DEVELOPERS *************************
     ************* Please do not update this method without consulting
     ************* Dr. Michael Mills.  This method has been intentionally
     ************* coded in this manner for precision/memory reasons.
     ************* *******************************************************

     */
    ErrorSet err = 0;

    double ratio = 1 - FLATTENING; //ratio of semi-minor axis over semi-major axis
    double d = 0.0;
    double d2 = 0.0;
    double sa = 0.0;
    double s1 = 0.0, f1 = 0.0, sx = 0.0, cx = 0.0;
    double delta = 1.0; // Setting != 0.0 will force one iteration.
    int index = 0;
    double tanl = 0.0;
    double newtu1 = 0.0;
    SPECIAL_DOUBLE tandl = 0.0;
    SPECIAL_DOUBLE newtu2 = 0.0;
    SPECIAL_DOUBLE newtu1_2 = 0.0;
    SPECIAL_DOUBLE newtu1_2_1 = 0.0;
    SPECIAL_DOUBLE sqrt_newtu1_2_1 = 0.0;
    SPECIAL_DOUBLE newcu1 = 0.0;
    SPECIAL_DOUBLE newtu2_2 = 0.0;
    SPECIAL_DOUBLE newtu2_2_1 = 0.0;
    SPECIAL_DOUBLE sqrt_newtu2_2_1 = 0.0;
    SPECIAL_DOUBLE newcu2 = 0.0;
    SPECIAL_DOUBLE cu2Xcx = 0.0;
    SPECIAL_DOUBLE su1Xcu2Xcx = 0.0;
    SPECIAL_DOUBLE b1M_su1Xcu2Xcx = 0.0;
    SPECIAL_DOUBLE tu2Xtu2 = 0.0;
    SPECIAL_DOUBLE tu1Xtu1 = 0.0;
    SPECIAL_DOUBLE tu1Xtu1Ptu2Xtu2 = 0.0;
    SPECIAL_DOUBLE s1Xcx = 0.0;
    SPECIAL_DOUBLE s1Xsx = 0.0;
    SPECIAL_DOUBLE saXsa = 0.0;
    SPECIAL_DOUBLE czDc2a = 0.0;
    SPECIAL_DOUBLE czXcz = 0.0;
    SPECIAL_DOUBLE czXczX2 = 0.0;
    SPECIAL_DOUBLE N3Xc2a = 0.0;
    SPECIAL_DOUBLE N3Xc2aP4 = 0.0;
    SPECIAL_DOUBLE N3Xc2aP4_XFLATTENING = 0.0;
    SPECIAL_DOUBLE N3Xc2aP4_XFLATTENING_P4 = 0.0;
    SPECIAL_DOUBLE N3Xc2aP4_XFLATTENING_P4_Xc2a = 0.0;
    SPECIAL_DOUBLE N3Xc2aP4_XFLATTENING_P4_Xc2aXFLATTENING = 0.0;
    SPECIAL_DOUBLE eXcy = 0.0;
    SPECIAL_DOUBLE eXcyXc = 0.0;
    SPECIAL_DOUBLE eXcyXcPcz = 0.0;
    SPECIAL_DOUBLE eXcyXcPcz_Xsy = 0.0;
    SPECIAL_DOUBLE eXcyXcPcz_XsyXc = 0.0;
    SPECIAL_DOUBLE eXcyXcPcz_XsyXc_Py = 0.0;
    SPECIAL_DOUBLE eXcyXcPcz_XsyXc_Py_Xsa = 0.0;
    SPECIAL_DOUBLE oneMc = 0.0;
    SPECIAL_DOUBLE oneMcXFLATTENING = 0.0;
    SPECIAL_DOUBLE xXoneMcXFLATTENING = 0.0;
    SPECIAL_DOUBLE xXoneMcXFLATTENINGPdestlon = 0.0;
    SPECIAL_DOUBLE xXoneMcXFLATTENINGPdestlonMoriglon = 0.0;
    double dM10 = 0.0;
    double xM10 = 0.0;
    double dM10_PxM10 = 0.0;

    // Pass "-D KEEP_USER_EPS" to avoid overriding eps
#ifndef KEEP_USER_EPS
    // This will override any eps value larger than the default value passed from
    // calling functions
    if (eps> DEFAULT_EPS)
    eps = DEFAULT_EPS;
#endif

    // cu1,su1,cu2,b1
    cx = 0;
    sx = 0;
    *c = 0;
    *e = 0;
    *cy = 0;
    *cz = 0;
    *sy = 0;
    *y = 0;
    *c2a = 0;

    /* All angles are in radians */

    if ((origin.latitude + dest.latitude == 0.0) && (fabs(origin.longitude
            - dest.longitude) == M_PI))
    {

        origin.latitude = origin.latitude + 0.00001; //allow algorithm to complete
    }

    tanl = tan(origin.latitude);
    newtu1 = ratio * tanl;

    //    *tu1 = ratio * tan(origin.latitude);
    *tu1 = newtu1;
    tandl = tan(dest.latitude);
    newtu2 = ratio * tandl;
    //    *tu2 = ratio * tan(dest.latitude);
    *tu2 = newtu2;
    newtu1_2 = newtu1 * newtu1;
    newtu1_2_1 = 1.0 + newtu1_2;
    sqrt_newtu1_2_1 = sqrt(newtu1_2_1);
    newcu1 = 1.0 / sqrt_newtu1_2_1;
    //    *cu1 = 1.0 / sqrt(1.0 + *tu1 * *tu1);
    *cu1 = newcu1;
    *su1 = (*cu1) * (*tu1);
    newtu2_2 = newtu2 * newtu2;
    newtu2_2_1 = newtu2_2 + 1.0;
    sqrt_newtu2_2_1 = sqrt(newtu2_2_1);
    newcu2 = 1.0 / sqrt_newtu2_2_1;
    //    *cu2 = 1.0 / sqrt(1.0 + *tu2 * *tu2);
    *cu2 = newcu2;
    s1 = (*cu1) * (*cu2);
    *b1 = s1 * (*tu2);
    f1 = (*b1) * (*tu1);
    *x = dest.longitude - origin.longitude;
    d = (*x) + 1.0; //force 1 iteration


    while ((fabs(delta) > eps) && (index < MAX_V_ITERATIONS))
    {
        index = index + 1;

        sx = sin(*x);

        cx = cos(*x);

        *tu1 = (*cu2) * sx;

        cu2Xcx = (*cu2) * cx;
        su1Xcu2Xcx = (*su1) * cu2Xcx;
        b1M_su1Xcu2Xcx = (*b1) - su1Xcu2Xcx;
//        *tu2 = (*b1)-(*su1)*(*cu2)*cx;
        *tu2 = b1M_su1Xcu2Xcx;

        //        *tu2 = (*b1) - (*su1) * (*cu2) * cx;


        if (*tu1 == 0.0 && *tu2 == 0.0)
        {
            //Very close points with the same longitude will
            //produce tu1 and tu2 both equal to 0.0. This
            //was causing the next calculation to produce
            //a QNaN. If we set both values to a very small
            //number, the iteration completes successfully.
            *tu1 = eps;
            *tu2 = eps;
        }

        tu2Xtu2 = (*tu2) * (*tu2);
        tu1Xtu1 = (*tu1) * (*tu1);
        tu1Xtu1Ptu2Xtu2 = tu1Xtu1 + tu2Xtu2;
        *sy = sqrt(tu1Xtu1Ptu2Xtu2);
        //        *sy = sqrt((*tu1) * (*tu1) + (*tu2) * (*tu2));
        s1Xcx = s1 * cx;
        *cy = s1Xcx + f1;
        //        *cy = s1 * cx + f1;
        *y = atan2(*sy, *cy);
        s1Xsx = s1 * sx;
        sa = s1Xsx / (*sy);
        //        sa = s1 * sx / (*sy);
        saXsa = sa * sa;
        *c2a = 1 - saXsa;
        //        *c2a = 1 - sa * sa;

        *cz = f1 + f1;

        if (*c2a > 0.0)
        {

            czDc2a = (*cz) / (*c2a);
            *cz = *cy - czDc2a;
            //            *cz = *cy - *cz / (*c2a);

        }

        czXcz = (*cz) * (*cz);
        czXczX2 = czXcz * 2.0;
        *e = czXczX2 - 1.0;
        //        *e = (*cz) * (*cz) * 2.0 - 1.0;


        N3Xc2a = -3.0 * (*c2a);
        N3Xc2aP4 = N3Xc2a + 4.0;
        N3Xc2aP4_XFLATTENING = N3Xc2aP4 * FLATTENING;
        N3Xc2aP4_XFLATTENING_P4 = N3Xc2aP4_XFLATTENING + 4.0;
        N3Xc2aP4_XFLATTENING_P4_Xc2a = N3Xc2aP4_XFLATTENING_P4 * (*c2a);
        N3Xc2aP4_XFLATTENING_P4_Xc2aXFLATTENING = N3Xc2aP4_XFLATTENING_P4_Xc2a
                * FLATTENING;
        *c = N3Xc2aP4_XFLATTENING_P4_Xc2aXFLATTENING / 16.0;
        //        *c = ((-3.0 * (*c2a) + 4.0) * FLATTENING + 4.0) * (*c2a) * FLATTENING / 16.0;

        d2 = d;
        d = *x;

        eXcy = (*e) * (*cy);
        eXcyXc = eXcy * (*c);
        eXcyXcPcz = eXcyXc + (*cz);
        eXcyXcPcz_Xsy = eXcyXcPcz * (*sy);
        eXcyXcPcz_XsyXc = eXcyXcPcz_Xsy * (*c);
        eXcyXcPcz_XsyXc_Py = eXcyXcPcz_XsyXc + (*y);
        eXcyXcPcz_XsyXc_Py_Xsa = eXcyXcPcz_XsyXc_Py * sa;
        *x = eXcyXcPcz_XsyXc_Py_Xsa;
        //        *x = (((*e) * (*cy) * (*c) + (*cz)) * (*sy) * (*c) + (*y)) * sa;


        oneMc = 1.0 - (*c);
        oneMcXFLATTENING = oneMc * FLATTENING;
        xXoneMcXFLATTENING = eXcyXcPcz_XsyXc_Py_Xsa * oneMcXFLATTENING;
        xXoneMcXFLATTENINGPdestlon = xXoneMcXFLATTENING + dest.longitude;
        xXoneMcXFLATTENINGPdestlonMoriglon = xXoneMcXFLATTENINGPdestlon
                - origin.longitude;
        *x = xXoneMcXFLATTENINGPdestlonMoriglon;
        //        *x = (*x) * (1.0 - *c) * FLATTENING + dest.longitude - origin.longitude;


//        if (index > 20)
//        {
//            //            printf("Into damping scheme\n");
//            // After 20 iterations, if *x is oscillating in the last
//            // digit, then apply a damping operation to force convergence
//            // Multiplying by 10 appears to provide one additional bit of precision
//
//            dM10 = d * 10.0;
//            xM10 = (*x) * 10.0;
//            dM10_PxM10 = dM10 + xM10;
//            *x = (double) 0.5 * dM10_PxM10;
//            //            *x = (double) 0.5 * (10.0 * d + 10.0 * (*x));
//            *x *= (double) 0.1;
//            //            printf("Avgd x = %.16e\n",*x);
//        }
        if(*x == d2){
			dM10 = d * 10.0;
			xM10 = (*x) * 10.0;
			dM10_PxM10 = dM10 + xM10;
			*x = 0.5 * dM10_PxM10;
			*x *= 0.1;
        	delta = 0;
        	break;
        }
        delta = d - *x;
        //NOTE all of these printfs need to be removed before committing again
        //printf();
    }
    if (fabs(delta) > eps)
    {
        /* did not converge -- return an error */
        err |= INV_NOT_CONVERGED_ERR;
    }

    return err;

}

/*******************************************************************************/
/*
 */
static ErrorSet iterateDirect(LLPoint origin, double course, double distance,
                          double eps, double *c2a, double *cu, double *cy,
                          double *su, double *sy, double *cf, double *sf,
                          double *sa, double *e, double *cz, double *y)
{
    /*
     EPS = 1e-9 gives results that exactly match Ed Williams' Javascript calculator
     See: <http://williams.best.vwh.net/gccalc.htm>

     In fact, this is a modified implementation of Ed Williams' code.

     The modifications made were in the manner in which the algorithm was coded.
     It was updated so as to minimize the impact of precision loss due to memory
     overflow and other rounding errors that may result when a mathematical
     operation is carried out using binary operations.

     Compare the general algorithm with Javascript at above site.

     ************* IMPORTANT NOTE FOR DEVELOPERS *************************
     ************* Please do not update this method without consulting
     ************* Dr. Michael Mills.  This method has been intentionally
     ************* coded in this manner for precision/memory reasons.
     ************* *******************************************************

     */

    ErrorSet err = 0;

    double ratio = 1 - FLATTENING; //ratio of semi-minor axis over semi-major axis

    double b, x, c, d;
    double tu;
    //    double tu = ratio * tan(origin.latitude);
    double delta = 1.0; //Setting != 0.0 will force one iteration

    int index = 0;

    double tan_origlat;
    double atan_tucf;
    double tuXtu;
    double tuXtuP1;
    double sqrt_tuXtuP1;
    double saXsa;
    double ratioXratio;
    double oneDratioXratio;
    double oneDratioXratio_M1;
    double c2a__oneDratioXratio_M1;
    double oneP_c2a__oneDratioXratio_M1;
    double sqrt_oneP_c2a__oneDratioXratio_M1;
    double xM2;
    double xXx;
    double xXx_D4;
    double xXx_D4_M1;
    double pointthreesevenfiveXx;
    double pointthreesevenfiveXxXx;
    double pointthreesevenfiveXxXxM1;
    double cXratio;
    double cXratioXSEMI;
    double bPy;
    double czXcz;
    double twoXczXcz;
    double twoXe;
    double syXsy;
    double syXsyX4;
    double syXsyX4M3;
    double syXsyX4M3_Xy;
    double syXsyX4M3_XyXcz;
    double syXsyX4M3_XyXczXd;
    double syXsyX4M3_XyXczXd_D6;
    double syXsyX4Md_XyXczXd_D6_Px;
    double syXsyX4Md_XyXczXd_D6_Px_Xd;
    double syXsyX4Md_XyXczXd_D6_Px_Xd_D4;
    double syXsyX4Md_XyXczXd_D6_Px_Xd_D4_Mcz;
    double syXsyX4Md_XyXczXd_D6_Px_Xd_D4_Mcz_Msy;
    double syXsyX4Md_XyXczXd_D6_Px_Xd_D4_Mcz_MsyXd;
    double tenXc;
    double tenXy;
    double tenXc_P_tenXy;

    // Pass "-D KEEP_USER_EPS" to avoid overriding eps
#ifndef KEEP_USER_EPS
    // This will override any eps value larger than the default value passed from
    // calling functions
    if (eps> DEFAULT_EPS)
    eps = DEFAULT_EPS;
#endif

    *sy = 0.0;
    *cy = 0.0;

    tan_origlat = tan(origin.latitude);
    tu = ratio * tan_origlat;

    *sf = sin(course);
    *cf = cos(course);

    if (*cf == 0.0)
    {

        b = 0.0;
    }
    else
    {

        atan_tucf = atan2(tu, *cf);
        b = 2.0 * atan_tucf;
        //        b = 2.0 * atan2(tu, *cf);
    }

    tuXtu = tu * tu;
    tuXtuP1 = 1.0 + tuXtu;
    sqrt_tuXtuP1 = sqrt(tuXtuP1);
    *cu = 1.0 / sqrt_tuXtuP1;
    //    *cu = 1.0 / sqrt(1.0 + tu * tu);


    *su = tu * (*cu);
    *sa = (*cu) * (*sf);

    saXsa = (*sa) * (*sa);
    *c2a = 1.0 - saXsa;
    //    *c2a = 1.0 - (*sa) * (*sa);


    ratioXratio = ratio * ratio;
    oneDratioXratio = 1.0 / ratioXratio;
    oneDratioXratio_M1 = oneDratioXratio - 1.0;
    c2a__oneDratioXratio_M1 = (*c2a) * oneDratioXratio_M1;
    oneP_c2a__oneDratioXratio_M1 = 1.0 + c2a__oneDratioXratio_M1;
    sqrt_oneP_c2a__oneDratioXratio_M1 = sqrt(oneP_c2a__oneDratioXratio_M1);
    x = 1.0 + sqrt_oneP_c2a__oneDratioXratio_M1;
    //    x = 1.0 + sqrt(1.0 + *c2a * (1.0 / (ratio * ratio) - 1.0));


    xM2 = x - 2.0;
    x = xM2 / x;
    //    x = (x - 2.0) / x;


    c = 1.0 - x;

    xXx = x * x;
    xXx_D4 = xXx / 4.0;
    xXx_D4_M1 = xXx_D4 + 1.0;
    c = xXx_D4_M1 / c;
    //    c = (x * x / 4.0 + 1.0) / c;


    pointthreesevenfiveXx = 0.375 * x;
    pointthreesevenfiveXxXx = pointthreesevenfiveXx * x;
    pointthreesevenfiveXxXxM1 = pointthreesevenfiveXxXx - 1.0;
    d = pointthreesevenfiveXxXxM1 * x;
    //    d = (0.375 * x * x - 1.0) * x;


    cXratio = c * ratio;
    cXratioXSEMI = cXratio * SEMI_MAJOR_AXIS_NMI;
    tu = distance / cXratioXSEMI;
    //    tu = distance / (c * ratio * SEMI_MAJOR_AXIS_NMI); // SEMI_MAJOR_AXIS_FEET);


    *y = tu;
    c = *y + 1.0;

    /* NOTE: If starting latitude is at either pole (where azimuth is not uniquely
     * defined), then the point returned will have longitude equal to the input course */
    while ((fabs(delta) > eps) && (index < MAX_V_ITERATIONS))
    {
        index = index + 1;

        *sy = sin(*y);
        *cy = cos(*y);

        bPy = b + *y;
        *cz = cos(bPy);
        //        *cz = cos(b + *y);


        czXcz = (*cz) * (*cz);
        twoXczXcz = 2.0 * czXcz;
        *e = twoXczXcz - 1.0;
        //        *e = 2.0 * (*cz) * (*cz) - 1.0;


        c = *y;
        x = (*e) * (*cy);

        twoXe = 2.0 * (*e);
        *y = twoXe - 1.0;
        //        *y = 2.0 * (*e) - 1.0;


        syXsy = (*sy) * (*sy);
        syXsyX4 = syXsy * 4.0;
        syXsyX4M3 = syXsyX4 - 3.0;
        syXsyX4M3_Xy = syXsyX4M3 * (*y);
        syXsyX4M3_XyXcz = syXsyX4M3_Xy * (*cz);
        syXsyX4M3_XyXczXd = syXsyX4M3_XyXcz * d;
        syXsyX4M3_XyXczXd_D6 = syXsyX4M3_XyXczXd / 6.0;
        syXsyX4Md_XyXczXd_D6_Px = syXsyX4M3_XyXczXd_D6 + x;
        syXsyX4Md_XyXczXd_D6_Px_Xd = syXsyX4Md_XyXczXd_D6_Px * d;
        syXsyX4Md_XyXczXd_D6_Px_Xd_D4 = syXsyX4Md_XyXczXd_D6_Px_Xd / 4.0;
        syXsyX4Md_XyXczXd_D6_Px_Xd_D4_Mcz = syXsyX4Md_XyXczXd_D6_Px_Xd_D4
                - (*cz);
        syXsyX4Md_XyXczXd_D6_Px_Xd_D4_Mcz_Msy
                = syXsyX4Md_XyXczXd_D6_Px_Xd_D4_Mcz * (*sy);
        syXsyX4Md_XyXczXd_D6_Px_Xd_D4_Mcz_MsyXd
                = syXsyX4Md_XyXczXd_D6_Px_Xd_D4_Mcz_Msy * d;
        *y = syXsyX4Md_XyXczXd_D6_Px_Xd_D4_Mcz_MsyXd + tu;
        //        *y = ( ( ( (*sy) * (*sy) * 4.0 - 3.0) * (*y) * (*cz) * d / 6.0 + x) * d
        //                / 4.0 - (*cz)) * (*sy) * d + tu;

        if (index > 20)
        {
            // After 20 iterations, if *y is oscillating in the last
            // digit, then apply a damping operation to force convergence
            // Multiplying by 10 appears to provide one additional bit of precision
            //TODO The proof is left as an exercise for the reader.

            tenXc = 10.0 * c;
            tenXy = 10.0 * (*y);
            tenXc_P_tenXy = tenXc + tenXy;
            *y = (double) 0.5 * tenXc_P_tenXy;
            //            *y = (double) 0.5 * (10.0 * c + 10.0 * (*y));


            *y *= (double) 0.1;
        }
        delta = *y - c;
    }

    if (fabs(delta) > eps)
    {
        /* did not converge -- return an error */
        err |= FWD_NOT_CONVERGED_ERR;
    }

    return err;

}

/*******************************************************************************/
ErrorSet geolib_idealab::direct(LLPoint origin, double course, double distance, LLPoint* dest,
               double eps)
{

    ErrorSet err = 0;

    double ratio = 1 - FLATTENING; //ratio of semi-minor axis over semi-major axis
    double c2a, cu, cy, su, sy, cf, sf, sa, e, cz, y;
    double b, c, d, x;
    LLPoint npDest;

    /* Assign local storage if optional pointers are not provided */
	if (NULL == dest) dest = &npDest;

    /* Negative distance implies movement along reciprocal course*/
    if (distance < 0.0)
    {
        course = modcrs(course + M_PI);
        distance = fabs(distance);
    }
    else if (fabs(distance) < ONE_NANOMETER_NMI) //TODO This really should use a TOL check
    {
        *dest = origin;
        return err;
    }

    err |= iterateDirect(origin, course, distance, eps, &c2a, &cu, &cy, &su,
            &sy, &cf, &sf, &sa, &e, &cz, &y);

    /* Compute new latitude */
    b = cu * cy * cf - su * sy;
    c = ratio * sqrt(sa * sa + b * b);
    d = su * cy + cu * sy * cf;
    dest->latitude = modlat(atan2(d, c)); /* output radians */


	/* Compute new longitude */
	c = cu * cy - su * sy * cf;
	x = atan2(sy * sf, c);
	c = ((-3.0 * c2a + 4.0) * FLATTENING + 4.0) * c2a * FLATTENING / 16.0;
	d = ((e * cy * c + cz) * sy * c + y) * sa;

	dest->longitude = modlon(origin.longitude + x - (1.0 - c) * d
			* FLATTENING);

    return err;

}

/*******************************************************************************/
/*
 * Compute destination latitude
 *
 */
ErrorSet geolib_idealab::directLat(LLPoint origin, double course, double distance, double* lat,
                  double eps)
{

    ErrorSet err = 0;
    LLPoint dest = { 0.0, 0.0 };

    double ratio = 1 - FLATTENING; //ratio of semi-minor axis over semi-major axis
    double c2a, cu, cy, su, sy, cf, sf, sa, e, cz, y;
    double b, c, d;
    double npLat;

    /* Assign local storage if optional pointers are not provided */
	if (NULL == lat) lat = &npLat;

    /* Negative distance implies movement along reciprocal course*/
    if (distance < 0)
    {
        course = modcrs(course + M_PI);
        distance = fabs(distance);
    }
    else if (fabs(distance) < ONE_NANOMETER_NMI)
    {
        *lat = origin.latitude;
        return err;
    }

    err |= iterateDirect(origin, course, distance, eps, &c2a, &cu, &cy, &su,
            &sy, &cf, &sf, &sa, &e, &cz, &y);

    b = cu * cy * cf - su * sy;
    c = ratio * sqrt(sa * sa + b * b);
    d = su * cy + cu * sy * cf;
    *lat = modlat(atan2(d, c)); /* output radians */

    return err;

}

/*******************************************************************************/
/*
 * Compute destination longitude
 *
 */
ErrorSet geolib_idealab::directLon(LLPoint origin, double course, double distance, double* lon,
                  double eps)
{

    LLPoint dest = { 0.0, 0.0 };
    ErrorSet err = 0;

    double ratio = 1 - FLATTENING; //ratio of semi-minor axis over semi-major axis
    double c2a, cu, cy, su, sy, cf, sf, sa, e, cz, y;
    double b, c, d, x;
    double npLon;

    /* Assign local storage if optional pointers are not provided */
	if (NULL == lon) lon = &npLon;

    /* Negative distance implies movement along reciprocal course*/
    if (distance < 0)
    {
        course = modcrs(course + M_PI);
        distance = fabs(distance);
    }
    else if (fabs(distance) < ONE_NANOMETER_NMI)
    {
        *lon = origin.longitude;
        return err;
    }

        err |= iterateDirect(origin, course, distance, eps, &c2a, &cu, &cy,
                &su, &sy, &cf, &sf, &sa, &e, &cz, &y);

        b = cu * cy * cf - su * sy;
        c = ratio * sqrt(sa * sa + b * b);
        d = su * cy + cu * sy * cf;
        c = cu * cy - su * sy * cf;
        x = atan2(sy * sf, c);
        c = ((-3.0 * c2a + 4.0) * FLATTENING + 4.0) * c2a * FLATTENING / 16.0;
        d = ((e * cy * c + cz) * sy * c + y) * sa;
        *lon = modlon(origin.longitude + x - (1.0 - c) * d * FLATTENING);

    return err;

}

/*******************************************************************************/
/*
 * Carries out inverse computation, returning course and distance
 *
 */

ErrorSet geolib_idealab::inverse(LLPoint origin, LLPoint dest, double* crs,
		double *bcrs, double* dist, double eps)
{

    // origin = lat/lon of starting point
    // dest   = lat/lon of destination point
    // crs    = forward azimuth at origin point
    // bcrs   = back azimuth at dest point
    // dist   = length of path in NM
    // eps    = round-off tolerance

    ErrorSet err = 0;
    double ratio = 1 - FLATTENING;
    double tu1, tu2, sy, e, cy, cz, y, c, c2a;
    double x1, d1, cu1, b1, su1, cu2, x;
    double npCrs = 0.0, npBcrs = 0.0, npDist = 0.0;
    int computeCrs = 1, computeBcrs = 1, computeDist = 1;
//    int invCount = 0;
//    double invDelta = 0.0;

    /* Assign local storage if optional pointers are not provided */
    if (NULL == crs)
    {
        crs = &npCrs;
        computeCrs = 0;
    }
    if (NULL == bcrs)
    {
        bcrs = &npBcrs;
        computeBcrs = 0;
    }
    if (NULL == dist)
    {
        dist = &npDist;
        computeDist = 0;
    }

    if ((origin.latitude == dest.latitude) && (origin.longitude
            == dest.longitude))
    {
        *dist = 0.0;
        *crs = 0.0;
        *bcrs = M_PI;
    }
    else
    {

        /* Set status code if points are nearly antipodal.  Computation can continue,
         * but user will be informed via non-zero return code that results are suspect.
         */
        err |= ptsAreAntipodal(origin, dest);

        //    if ((fabs(origin.latitude - dest.latitude) > eps) || (fabs(origin.longitude
        //            - dest.longitude) > eps))
        //    {

        err |= iterateInverse(origin, dest, eps, &tu1, &tu2, &sy, &e, &cy, &cz,
                &y, &c, &c2a, &cu1, &b1, &su1, &cu2, &x);

        /* Set course output (it's in radians) */
        if (computeCrs)
            *crs = modcrs(atan2(tu1, tu2));
        if (computeBcrs)
            *bcrs = modcrs(atan2(cu1 * sin(x), b1 * cos(x) - su1 * cu2) + M_PI);
        if (computeDist)
        {
            x1 = sqrt((1 / (ratio * ratio) - 1.0) * c2a + 1.0);
            x1 += 1.0;
            x1 = (x1 - 2.0) / x1;

            c = 1.0 - x1;
            c = (x1 * x1 / 4.0 + 1.0) / c;
            d1 = (0.375 * x1 * x1 - 1.0) * x1;
            x1 = e * cy;
            /* Set distance output */
            *dist = ((((sy * sy * 4.0 - 3.0) * (1.0 - e - e) * cz * d1 / 6.0
                    - x1) * d1 / 4.0 + cz) * sy * d1 + y) * c
                    * SEMI_MAJOR_AXIS_NMI * ratio;
        }
    }

    return err;

}

/*******************************************************************************/
/*
 * Carries out inverse computation but returns only the distance.  Slightly
 * faster than inverse if you don't need the course.
 *
 */

ErrorSet geolib_idealab::invDist(LLPoint origin, LLPoint dest, double* dist, double eps)
{

    /* Use optional arguments to inverse to simplify this function */
    return geolib_idealab::inverse(origin, dest, NULL, NULL, dist, eps);

}

/*******************************************************************************/
/*
 * Carries out inverse computation but returns only the course.  Slightly
 * faster than inverse if you don't need to know the distance.
 *
 */

ErrorSet geolib_idealab::invCrs(LLPoint origin, LLPoint dest, double* fcrs, double* bcrs,
                 double eps)
{

    /* Use optional arguments to inverse to simplify this function */
    return geolib_idealab::inverse(origin, dest, fcrs, bcrs, NULL, eps);

}

/*******************************************************************************
 * Test whether a point lies on line segment defined by startPt and endPt
 *
 * INPUT
 *     startPt: Starting point that defines the line segment
 *       endPt: End point that defines the line segment
 *      testPt: Point to test against line segment
 *  lengthCode: Integer that represents boundedness of line segment:
 *              Allowable values: 0 = finite = macro SEGMENT
 *                                1 = semi-infinite (extends beyond endPt)
 *                                  = macro SEMIINFINITE
 *                                2 = infinite (extends beyond startPt and endPt)
 *                                  = macro INFINITE
 *         tol: Minimum distance from line segment in order for test point
 *              to be considered on it.
 *         eps: Accuracy tolerance for direct/inverse algorithms
 * OUTPUT
 *   1: if the point satisfies all conditions
 *   0: otherwise
 */

int geolib_idealab::ptIsOnGeo(LLPoint startPt, LLPoint endPt, LLPoint testPt,
                    LineType lengthCode, ErrorSet* err, double tol, double eps)
{

    LLPoint newStart, newEnd;

    double dist12, crs12, crs21;
    double crsTest1, crsTest2;
    double dist1Test, dist2Test;

    int returnVal = 0;
    int onCrs = 0;
    int betweenEnds = 0;

    ErrorSet newErr = 0;

    if (ptsAreSame(startPt, testPt, tol))
    {
        /* Point coincides with start point */
        returnVal = 1;
    }
    else if (ptsAreSame(endPt, testPt, tol))
    {
        /* Point coincides with end point */
        returnVal = 1;
    }
    else
    {

        returnVal = 0;
        newErr |= geolib_idealab::inverse(startPt, endPt, &crs12, &crs21, &dist12, eps);
        if (ptIsOnCrs(startPt, crs12, testPt, &crsTest1, &dist1Test,
                &newErr, tol, eps))
        {
            onCrs = 1;
            betweenEnds = ((dist1Test > 0.0) && (dist1Test < dist12));
        }
        else if ((dist1Test > 0.0) && (dist1Test < 10.0 / 1852.0))
        {
            /* Test point is extremely close to startPt.  Courses are not
             * accurate enough in this case.  Move startPt back 1.0 nm and test again */
            newErr |= geolib_idealab::direct(startPt, crs12 + M_PI, 1.0, &newStart, eps);
            newErr |= geolib_idealab::invCrs(newStart, endPt, &crs12, &crs21, eps);

            if (ptIsOnCrs(newStart, crs12, testPt, &crsTest1, &dist1Test,
                    &newErr, tol, eps))
            {
                onCrs = 1;
                betweenEnds
                        = ((dist1Test >= 1.0) && (dist1Test - 1.0 < dist12));
            }
        }

        if (onCrs && ((lengthCode != SEGMENT) || (betweenEnds)))
        {
            returnVal = 1;
        }
        else if (lengthCode == INFINITE)
        {

            if (ptIsOnCrs(endPt, crs21, testPt, &crsTest2, &dist2Test,
                    &newErr, tol, eps))
            {
                returnVal = 1;
            }
            else if ((dist2Test > 0.0) && (dist2Test < 10.0 / 1852.0))
            {

                newErr |= geolib_idealab::direct(startPt, crs12 + M_PI, dist12 + 1.0,
                        &newEnd, eps);
                newErr |= geolib_idealab::invCrs(startPt, newEnd, &crs12, &crs21, eps);

                if (ptIsOnCrs(newEnd, crs21, testPt, &crsTest2,
                        &dist2Test, &newErr, tol, eps))
                {
                    returnVal = 1;
                }
            }
            else
            {
                returnVal = 0;
            }
        }
        else
        {
            returnVal = 0;
        }

    }

    if (newErr != 0)
    {
        *err |= newErr;
        returnVal = 0;
    }

    return returnVal;

}

/*
 * Starting at startPt and following crs12, will one eventually reach testPt?
 */
int geolib_idealab::ptIsOnCrs(LLPoint startPt, double crs12, LLPoint testPt,
                   double* crsTest1, double* dist1Test, ErrorSet* err,
                   double tol, double eps)
{

    LLPoint comparePt;
    double crs1Test;
    double npCrsTest1;
    double npDist1Test;
    ErrorSet newErr = 0;

    /* Assign local storage if optional pointers are not provided */
    if (NULL == crsTest1) crsTest1 = &npCrsTest1;
    if (NULL == dist1Test) dist1Test = &npDist1Test;

    if (ptsAreSame(startPt, testPt, tol))
    {
        /* Point coincides with start point */
        /* crsTest1 undefined */
        *dist1Test = 0.0;
        return 1;
    }

    newErr |= geolib_idealab::inverse(startPt, testPt, &crs1Test, crsTest1, dist1Test,eps);
    newErr |= geolib_idealab::direct(startPt, crs12, *dist1Test, &comparePt, eps);
    if (fabs(modlon(crs12 - crs1Test)) > M_PI_2)
    {
        /* testPt is behind startPt, useful for calling function to know this */
        *dist1Test = -fabs(*dist1Test);
    }

    if ((newErr == 0) && ptsAreSame(testPt, comparePt, tol))
    {
        return 1;
    }
    else
    {
        //TODO This code below should be uncommented to fix bug 17561.
        /*Last possibility: the start point is at a pole in which case
         the test point is definitely on course. It is better to put
         this check at the end of this function so that it is rarely called. /
         if (ptIsAtPole(startPt,&newErr,tol,eps) != 0)
         {
         //Return 1.
         *err |= newErr;
         return 1;
         }
         else
         {*/
        *err |= newErr;
        return 0;
        //}
    }

}

/*******************************************************************************
 * Initialize a new Geodesic structure from given input parameters
 */
ErrorSet geolib_idealab::createGeo(Geodesic* geo, LLPoint geoStart, LLPoint geoEnd,
                      LineType lineType, double eps)
{

    double sCrs, eCrs, length;

    ErrorSet err = 0;

    Geodesic npGeo;
    /* Assign local storage if optional pointers are not provided */
    if (NULL == geo) geo = &npGeo;


    geo->startPoint = geoStart;
    geo->endPoint = geoEnd;
    geo->lineType = lineType;
    err |= geolib_idealab::inverse(geoStart, geoEnd, &sCrs, &eCrs, &length, eps);
    geo->startAz = sCrs;
    geo->endAz = modpos(eCrs + M_PI, M_2PI);
    geo->length = length;

    return err;

}

/******************************************************************************
 * Find the course of geodesic at given point. Valid return values are in
 * range [0, 2*PI).  Invalid return value is -1.0.
 */
double geolib_idealab::geoCrs(Geodesic geo, LLPoint testPt, double* startCrs,
                               double* revCrs, double* distToPt, ErrorSet* err,
                               double tol, double eps)
{

    double crs;
    double geoCrs, geoRevCrs;
    double distToStart, distToEnd, crsToStart, crsToEnd;
    double crsFromStart, crsFromEnd;
    double npCrs, npRevCrs, npDist;
    ErrorSet newErr = 0;

    /* Assign local storage if optional pointers are not provided */
    if (NULL == startCrs) startCrs = &npCrs;
    if (NULL == revCrs) revCrs = &npRevCrs;
    if (NULL == distToPt) distToPt = &npDist;

    /* Compute course of geodesic at important points */
    newErr |= geolib_idealab::invCrs(geo.startPoint, geo.endPoint, &geoCrs, &geoRevCrs,
            eps); //start to end
    if (newErr)
    {
        *err |= newErr;
        return -1.0;
    }
    else
    {
        // Return these courses
        *startCrs = geoCrs;
        *revCrs = geoRevCrs;
    }
    newErr |= geolib_idealab::inverse(testPt, geo.startPoint, &crsToStart, &crsFromStart,
            distToPt, eps); //test to start, return distToPt
    if (newErr)
    {
        *err |= newErr;
        return -1.0;
    }
    newErr |= geolib_idealab::inverse(testPt, geo.endPoint, &crsToEnd, &crsFromEnd,
            &distToEnd, eps); //test to end
    if (newErr)
    {
        *err |= newErr;
        return -1.0;
    }

    //Now determine what should be returned based on the location of testPt to start and end point.
    if (*distToPt < tol)
    {
        /* testPt is indistinguishable from startPt */
        crs = geoCrs;
    }
    else if (distToEnd < tol)
    {
        /* testPt is indistinguishable from endPt */
        crs = geoRevCrs + M_PI;
    }
    else if (geolib_idealab::ptIsOnCrs(geo.startPoint, geoCrs, testPt, &crsToStart,
            &distToStart, &newErr, tol, eps))
    {
        // testPt lies on geodesic extending from startPt and along geoCrs. If distToStart is negative,
        // then testPt is behind startPt.
        if (distToStart < 0.0)
        {
            // Behind start point. Check lineType
            if (geo.lineType == INFINITE)
            {
                // Return the course from test to start
                crs = crsToEnd;
            }
            else
            {
                return -1.0;
            }
        }
        else if (distToStart > geo.length)
        {
            // testPt is beyond endPt. Check lineType.
            if (geo.lineType != SEGMENT)
            {
                // Beyond end point
                crs = crsToStart + M_PI;
            }
            else
            {
                return -1.0;
            }
        }
        else
        {
            // testPt extends from startPt along geoCrs and is btwn startPt and endPt. Return crs of geodetic
            crs = crsToEnd;
        }
    }
    else if (geolib_idealab::ptIsOnCrs(geo.endPoint, geoRevCrs, testPt, &crsToEnd,
            &distToEnd, &newErr, tol, eps))
    {
        // testPt extends from endPt and along geoRevCrs.
        if (distToEnd < 0.0)
        {
            // testPt is beyond endPt, so look back toward startPt, then
            // return reciprocal course
            if (geo.lineType != SEGMENT)
            {
                // Beyond end point
                crs = crsToStart + M_PI;
            }
            else
            {
                return -1.0;
            }
        }
        else if (distToEnd > geo.length)
        {
            // testPt is behind startPt. Check lineType
            if (geo.lineType == INFINITE)
            {
                // Beyond start point
                crs = crsToEnd;
            }
            else
            {
                return -1.0;
            }
        }
        else
        {
            // testPt is btwn startPt and endPt
            crs = crsToEnd;
        }
    }
    else
    {
        /* No solution; return invalid value */
        return -1.0;
        //debug if (ptIsOnGeo(geo.startPoint,geo.endPoint,testPt,geo.length,&err,tol,eps))
        //debug {
        //debug     newErr |= UNEXPECTED_ERR;
        //debug     printf("DOH!");
        //debug }
    }

    /* Make sure return value is in valid range [0, 2*PI) */
    return modcrs(crs);
}
