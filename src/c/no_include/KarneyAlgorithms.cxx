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
#include<float.h>
#if REPLACE_WITH_AMDLIBM
#include "amdlibm.h"
#endif
#include "Geolib.h"

double sq(double x){
	return x*x;
}

void swap(double *x, double *y)
{
    double t;
    t = *x;
    *x = *y;
    *y = t;
}

static inline double AngNormalize(double x){
	// Place angle in [-180, 180).  Assumes x is in [-540, 540).

	return x >= 180 ? x - 360 : (x < -180 ? x + 360 : x);
}

static inline double AngRound(double x){
	// The makes the smallest gap in x = 1/16 - nextafter(1/16, 0) = 1/2^57
	// for doubles = 0.7 pm on the earth if x is an angle in degrees.  (This
	// is about 1000 times more resolution than we get with angles around 90
	// degrees.)  We use this to avoid having to deal with near singular
	// cases when x is non-zero but tiny (e.g., 1.0e-200).

	const double z = (double)(0.0625); // 1/16
	volatile double y = fabs(x);

	// The compiler mustn't "simplify" z - (z - y) to y
	y = y < z ? z - (z - y) : y;

	return x < 0 ? -y : y;
}

static inline void SinCosNorm(double* sinx, double* cosx){

	double r = hypot(*sinx, *cosx);
	*sinx /= r;
	*cosx /= r;

}

double SinCosSeries(int sinp, double sinx, double cosx, const double c[], int n){
	// Evaluate
	// y = sinp ? sum(c[i] * sin( 2*i    * x), i, 1, n) :
	//            sum(c[i] * cos((2*i+1) * x), i, 0, n-1) :
	// using Clenshaw summation.  N.B. c[0] is unused for sin series
	// Approx operation count = (n + 5) mult and (2 * n + 2) add
	c += (n + sinp);            // Point to one beyond last element
	double
	  ar = 2 * (cosx - sinx) * (cosx + sinx), // 2 * cos(2 * x)
	  y0 = n & 1 ? *--c : 0, y1 = 0;          // accumulators for sum
	// Now n is even
	n /= 2;
	while (n--) {
		// Unroll loop x 2, so accumulators return to their original role
		y1 = ar * y0 - y1 + *--c;
		y0 = ar * y1 - y0 + *--c;
	}
	return sinp
	  ? 2 * sinx * cosx * y0    // sin(2 * x) * y0
	  : cosx * (y0 - y1);       // cos(x) * (y0 - y1)
}


double A1m1f(double eps){
	// The scale factor A1-1 = mean value of I1-1
	double eps2 = sq(eps), t;

	t = eps2*(eps2*(eps2+4)+64)/256;
	return (t + eps) / (1 - eps);
}


double A2m1f(double eps){
	// The scale factor A2-1 = mean value of I2-1
	double eps2 = sq(eps), t;

	t = eps2*(eps2*(25*eps2+36)+64)/256;
	return t * (1 - eps) - eps;
}

void C1f(double eps, double c[]){
	// The coefficients C1[l] in the Fourier expansion of B1
	double eps2 = sq(eps), d = eps;

	c[1] = d*((6-eps2)*eps2-16)/32;
	d *= eps;
	c[2] = d*((64-9*eps2)*eps2-128)/2048;
	d *= eps;
	c[3] = d*(9*eps2-16)/768;
	d *= eps;
	c[4] = d*(3*eps2-5)/512;
	d *= eps;
	c[5] = -7*d/1280;
	d *= eps;
	c[6] = -7*d/2048;
}


void C2f(double eps, double c[]){
	// The coefficients C2[l] in the Fourier expansion of B2
	double eps2 = sq(eps), d = eps;

	c[1] = d*(eps2*(eps2+2)+16)/32;
	d *= eps;
	c[2] = d*(eps2*(35*eps2+64)+384)/2048;
	d *= eps;
	c[3] = d*(15*eps2+80)/768;
	d *= eps;
	c[4] = d*(7*eps2+35)/512;
	d *= eps;
	c[5] = 63*d/1280;
	d *= eps;
	c[6] = 77*d/2048;
}

void C1pf(double eps, double c[]){
	// The coefficients C1p[l] in the Fourier expansion of B1p
	double eps2 = sq(eps), d = eps;

	c[1] = d*(eps2*(205*eps2-432)+768)/1536;
	d *= eps;
	c[2] = d*(eps2*(4005*eps2-4736)+3840)/12288;
	d *= eps;
	c[3] = d*(116-225*eps2)/384;
	d *= eps;
	c[4] = d*(2695-7173*eps2)/7680;
	d *= eps;
	c[5] = 3467*d/7680;
	d *= eps;
	c[6] = 38081*d/61440;
}

void C3f(double eps, double c[]){
  // Evaluation C3 coeffs by Horner's method
  // Elements c[1] thru c[nC3_ - 1] are set
	int nC3x = (6 * (6 - 1)) / 2;
	  double n = THIRD_FLATTENING;
	double C3x[nC3x];
	//  C3coeff
	    C3x[0] = (1-n)/4;
	    C3x[1] = (1-n*n)/8;
	    C3x[2] = ((3-n)*n+3)/64;
	    C3x[3] = (2*n+5)/128;
	    C3x[4] = 3/(double)(128);
	    C3x[5] = ((n-3)*n+2)/32;
	    C3x[6] = ((-3*n-2)*n+3)/64;
	    C3x[7] = (n+3)/128;
	    C3x[8] = 5/(double)(256);
	    C3x[9] = (n*(5*n-9)+5)/192;
	    C3x[10] = (9-10*n)/384;
	    C3x[11] = 7/(double)(512);
	    C3x[12] = (7-14*n)/512;
	    C3x[13] = 7/(double)(512);
	    C3x[14] = 21/(double)(2560);

	int i, j, k;

  for (j = nC3x, k = 6 - 1; k; ) {
    double t = 0;
    for (i = 6 - k; i; --i)
      t = eps * t + C3x[--j];
    c[k--] = t;
  }

  double mult = 1;
  for (k = 1; k < 6; ) {
    mult *= eps;
    c[k++] *= mult;
  }
}

void C4f(double k2, double c[]){
  // Evaluation C4 coeffs by Horner's method
  // Elements c[0] thru c[nC4_ - 1] are set
    static const int nC4x = (6 * (6 + 1)) / 2;
//	const double f = FLATTENING;
//	const double f1 = 1 - f;
//	const double e2 = f*(2 - f);
//	const double ep2 = e2/sq(f1);
//	const double e2 = ECCENTRICITY_SQ;
	const double ep2 = SECOND_ECCENTRICITY_SQ;
    //  C4coeff
    double C4x[nC4x];
        C4x[0] = (ep2*(ep2*(ep2*((832-640*ep2)*ep2-1144)+1716)-3003)+
                  30030)/45045;
        C4x[1] = (ep2*(ep2*((832-640*ep2)*ep2-1144)+1716)-3003)/60060;
        C4x[2] = (ep2*((208-160*ep2)*ep2-286)+429)/18018;
        C4x[3] = ((104-80*ep2)*ep2-143)/10296;
        C4x[4] = (13-10*ep2)/1430;
        C4x[5] = -1/(double)(156);
        C4x[6] = (ep2*(ep2*(ep2*(640*ep2-832)+1144)-1716)+3003)/540540;
        C4x[7] = (ep2*(ep2*(160*ep2-208)+286)-429)/108108;
        C4x[8] = (ep2*(80*ep2-104)+143)/51480;
        C4x[9] = (10*ep2-13)/6435;
        C4x[10] = 5/(double)(3276);
        C4x[11] = (ep2*((208-160*ep2)*ep2-286)+429)/900900;
        C4x[12] = ((104-80*ep2)*ep2-143)/257400;
        C4x[13] = (13-10*ep2)/25025;
        C4x[14] = -1/(double)(2184);
        C4x[15] = (ep2*(80*ep2-104)+143)/2522520;
        C4x[16] = (10*ep2-13)/140140;
        C4x[17] = 5/(double)(45864);
        C4x[18] = (13-10*ep2)/1621620;
        C4x[19] = -1/(double)(58968);
        C4x[20] = 1/(double)(792792);

	int i, j, k;

  for (j = nC4x, k = 6; k; ) {
    double t = 0;
    for (i = 6 - k + 1; i; --i)
      t = k2 * t + C4x[--j];
    c[--k] = t;
  }

  double mult = 1;
  for (k = 1; k < 6; ) {
    mult *= k2;
    c[k++] *= mult;
  }
}

double A3f(double eps){
	// Evaluation sum(_A3c[k] * eps^k, k, 0, nA3x_-1) by Horner's method
	double v = 0;
	int i;
	int nA3x = 6;
//	double f = FLATTENING;
//	double n = f/(2-f);
	double n = THIRD_FLATTENING;
	double A3x[nA3x];
	//  A3coeff
	A3x[0] = 1;
	A3x[1] = (n-1)/2;
	A3x[2] = (n*(3*n-1)-2)/8;
	A3x[3] = ((-n-3)*n-1)/16;
	A3x[4] = (-2*n-3)/64;
	A3x[5] = -3/(double)(128);

	for (i = nA3x; i; )
	v = eps * v + A3x[--i];
	return v;
}

long DirectKarney(double lat1, double lon1, double azi1, int arcmode, double s12_a12, unsigned outmask,
                                 double* lat2, double* lon2, double* azi2,
                                 double* s12, double* m12, double* M12, double* M21,
                                 double* S12){
	long err = 0;

	// Underflow guard.  We require
	//   tiny_ * epsilon() > 0
	//   tiny_ + epsilon() == epsilon()
	const double tiny = sqrt(DBL_MIN);
//	const double tol0 = DBL_EPSILON;
	// Increase multiplier in defn of tol1_ from 100 to 200 to fix inverse case
	// 52.784459512564 0 -52.784459512563990912 179.634407464943777557
	// which otherwise failed for Visual Studio 10 (Release and Debug)
//	const double tol1 = 200 * tol0;
//	const double tol2 = sqrt(DBL_EPSILON);
//	const double xthresh = 1000 * tol2;

    double DEG2RAD = M_PI / 180.0;
	const double a = SEMI_MAJOR_AXIS_METERS;
	const double f = FLATTENING;
	const double f1 = 1 - f;
//	const double e2 = f*(2 - f);
//	const double ep2 = e2/sq(f1);
	const double e2 = ECCENTRICITY_SQ;
	const double ep2 = SECOND_ECCENTRICITY_SQ;
//	const double n = f/(2-f);
//	const double b = a*f1;
	const double b = SEMI_MINOR_AXIS_METERS;
//	const double c2 = (sq(a) + sq(b) *
//	           (e2 == 0 ? 1 :
//	            (e2 > 0 ? atanh(sqrt(e2)) : atan(sqrt(-e2))) /
//	            sqrt(fabs(e2))))/2; // authalic radius squared
	      // The sig12 threshold for "(double)ly short"
//	const double etol2 = 10 * tol2 / fmax(0.1, sqrt(fabs(e2)));

    double salp0, calp0, k2,
      salp1, calp1, ssig1, csig1, stau1, ctau1, somg1, comg1,
      A1m1, A2m1, A3c, B11, B21, B31, A4, B41;

    azi1 = AngNormalize(azi1);
    // Guard against underflow in salp0
    azi1 = AngRound(azi1);
    lon1 = AngNormalize(lon1);

    // alp1 is in [0, pi]
    double alp1 = azi1 * DEG2RAD;
    // Enforce sin(pi) == 0 and cos(pi/2) == 0.  Better to face the ensuing
    // problems directly than to skirt them.
    salp1 =      azi1  == -180 ? 0 : sin(alp1);
    calp1 = fabs(azi1) ==   90 ? 0 : cos(alp1);

    double cbet1, sbet1, phi;
    phi = lat1 * DEG2RAD;
    // Ensure cbet1 = +epsilon at poles
    sbet1 = f1 * sin(phi);
    cbet1 = fabs(lat1) == 90 ? tiny : cos(phi);
    SinCosNorm(&sbet1, &cbet1);

    // Evaluate alp0 from sin(alp1) * cos(bet1) = sin(alp0),
    salp0 = salp1 * cbet1; // alp0 in [0, pi/2 - |bet1|]
    // Alt: calp0 = hypot(sbet1, calp1 * cbet1).  The following
    // is slightly better (consider the case salp1 = 0).
    calp0 = hypot(calp1, salp1 * sbet1);
    // Evaluate sig with tan(bet1) = tan(sig1) * cos(alp1).
    // sig = 0 is nearest northward crossing of equator.
    // With bet1 = 0, alp1 = pi/2, we have sig1 = 0 (equatorial line).
    // With bet1 =  pi/2, alp1 = -pi, sig1 =  pi/2
    // With bet1 = -pi/2, alp1 =  0 , sig1 = -pi/2
    // Evaluate omg1 with tan(omg1) = sin(alp0) * tan(sig1).
    // With alp0 in (0, pi/2], quadrants for sig and omg coincide.
    // No atan2(0,0) ambiguity at poles since cbet1 = +epsilon.
    // With alp0 = 0, omg1 = 0 for alp1 = 0, omg1 = pi for alp1 = pi.
    ssig1 = sbet1; somg1 = salp0 * sbet1;
    csig1 = comg1 = sbet1 != 0 || calp1 != 0 ? cbet1 * calp1 : 1;
    SinCosNorm(&ssig1, &csig1); // sig1 in (-pi, pi]
    SinCosNorm(&somg1, &comg1);

    k2 = sq(calp0) * ep2;
    double eps = k2 / (2 * (1 + sqrt(1 + k2)) + k2);

    // index zero elements of _C1a, _C1pa, _C2a, _C3a are unused
    double C1a[6 + 1], C1pa[6 + 1], C2a[6 + 1], C3a[6];


	A1m1 = A1m1f(eps);
	C1f(eps, C1a);
	B11 = SinCosSeries(1, ssig1, csig1, C1a, 6);
	double s = sin(B11), c = cos(B11);
	// tau1 = sig1 + B11
	stau1 = ssig1 * c + csig1 * s;
	ctau1 = csig1 * c - ssig1 * s;
	// Not necessary because C1pa reverts C1a
	//    _B11 = -SinCosSeries(true, _stau1, _ctau1, _C1pa, nC1p_);

	C1pf(eps, C1pa);

//	A2m1 = A2m1f(eps);
//	C2f(eps, C2a);
//	B21 = SinCosSeries(1, ssig1, csig1, C2a, 6);

	C3f(eps, C3a);
	A3c = -f * salp0 * A3f(eps);
	B31 = SinCosSeries(1, ssig1, csig1, C3a, 6-1);

//	double C4a[6];    // all the elements of _C4a are used
//	C4f(k2, C4a);
//	// Multiplier = a^2 * e^2 * cos(alpha0) * sin(alpha0)
//	A4 = sq(a) * calp0 * salp0 * e2;
//	B41 = SinCosSeries(0, ssig1, csig1, C4a, 6);

//GenPosition
	double sig12, ssig12, csig12, B12 = 0, AB1 = 0;

	// Interpret s12_a12 as distance
	double tau12 = s12_a12 / (b * (1 + A1m1));

	s = sin(tau12);
	c = cos(tau12);
	// tau2 = tau1 + tau12
	B12 = - SinCosSeries(1, stau1 * c + ctau1 * s,
								 ctau1 * c - stau1 * s,
								 C1pa, 6);
	sig12 = tau12 - (B12 - B11);
	ssig12 = sin(sig12);
	csig12 = cos(sig12);

	double omg12, lam12, lon12;
	double ssig2, csig2, sbet2, cbet2, somg2, comg2, salp2, calp2;
	// sig2 = sig1 + sig12
	ssig2 = ssig1 * csig12 + csig1 * ssig12;
	csig2 = csig1 * csig12 - ssig1 * ssig12;

	AB1 = (1 + A1m1) * (B12 - B11);


	// sin(bet2) = cos(alp0) * sin(sig2)
	sbet2 = calp0 * ssig2;
	// Alt: cbet2 = hypot(csig2, salp0 * ssig2);
	cbet2 = hypot(salp0, calp0 * csig2);
	if (cbet2 == 0)
	// I.e., salp0 = 0, csig2 = 0.  Break the degeneracy in this case
	cbet2 = csig2 = tiny;
	// tan(omg2) = sin(alp0) * tan(sig2)
	somg2 = salp0 * ssig2; comg2 = csig2;  // No need to normalize
	// tan(alp0) = cos(sig2)*tan(alp2)
	salp2 = salp0; calp2 = calp0 * csig2; // No need to normalize
	// omg12 = omg2 - omg1
	omg12 = atan2(somg2 * comg1 - comg2 * somg1,
				comg2 * comg1 + somg2 * somg1);

	//Longitude
	lam12 = omg12 + A3c *
	( sig12 + (SinCosSeries(1, ssig2, csig2, C3a, 6-1)
			   - B31));
	lon12 = lam12 / DEG2RAD;
	// Can't use AngNormalize because longitude might have wrapped multiple
	// times.
	lon12 = lon12 - 360 * floor(lon12/360 + (double)(0.5));
	*lon2 = AngNormalize(lon1 + lon12);

	//Latitude
	*lat2 = atan2(sbet2, f1 * cbet2) / DEG2RAD;

	//Azimuth
	// minus signs give range [-180, 180). 0- converts -0 to +0.
	*azi2 = 0 - atan2(-salp2, calp2) / DEG2RAD;

//	double
//	  ssig1sq = sq( ssig1),
//	  ssig2sq = sq( ssig2),
//	  w1 = sqrt(1 + k2 * ssig1sq),
//	  w2 = sqrt(1 + k2 * ssig2sq),
//	  B22 = SinCosSeries(1, ssig2, csig2, C2a, 6),
//	  AB2 = (1 + A2m1) * (B22 - B21),
//	  J12 = (A1m1 - A2m1) * sig12 + (AB1 - AB2);
//	//REDUCEDLENGTH
//	// Add parens around (_csig1 * ssig2) and (_ssig1 * csig2) to ensure
//	// accurate cancellation in the case of coincident points.
//	*m12 = b * ((w2 * (csig1 * ssig2) - w1 * (ssig1 * csig2))
//			- csig1 * csig2 * J12);
//	//GEODESICSCALE
//	*M12 = csig12 + (k2 * (ssig2sq - ssig1sq) *  ssig2 / (w1 + w2)
//				  - csig2 * J12) * ssig1 / w1;
//	*M21 = csig12 - (k2 * (ssig2sq - ssig1sq) * ssig1 / (w1 + w2)
//				  - csig1 * J12) * ssig2 / w2;


	return err;

}

void Lengths(double eps, double sig12,
                       double ssig1, double csig1, double ssig2, double csig2,
                       double cbet1, double cbet2,
                       double* s12b, double* m12a, double* m0,
                       int scalep, double* M12, double* M21,
                       // Scratch areas of the right size
                       double C1a[], double C2a[]){

	const double f = FLATTENING;
	const double f1 = 1 - f;
//	const double e2 = f*(2 - f);
	const double e2 = ECCENTRICITY_SQ;

  // Return m12a = (reduced length)/_a; also calculate s12b = distance/_b,
  // and m0 = coefficient of secular term in expression for reduced length.
  C1f(eps, C1a);
  C2f(eps, C2a);
  double
    A1m1 = A1m1f(eps),
    AB1 = (1 + A1m1) * (SinCosSeries(1, ssig2, csig2, C1a, 6) -
                        SinCosSeries(1, ssig1, csig1, C1a, 6)),
    A2m1 = A2m1f(eps),
    AB2 = (1 + A2m1) * (SinCosSeries(1, ssig2, csig2, C2a, 6) -
                        SinCosSeries(1, ssig1, csig1, C2a, 6)),
    cbet1sq = sq(cbet1),
    cbet2sq = sq(cbet2),
    w1 = sqrt(1 - e2 * cbet1sq),
    w2 = sqrt(1 - e2 * cbet2sq),
    // Make sure it's OK to have repeated dummy arguments
    m0x = A1m1 - A2m1,
    J12 = m0x * sig12 + (AB1 - AB2);
  *m0 = m0x;
  // Missing a factor of _a.
  // Add parens around (csig1 * ssig2) and (ssig1 * csig2) to ensure accurate
  // cancellation in the case of coincident points.
  *m12a = (w2 * (csig1 * ssig2) - w1 * (ssig1 * csig2))
    - f1 * csig1 * csig2 * J12;
  // Missing a factor of _b
  *s12b = (1 + A1m1) * sig12 + AB1;
  if (scalep) {
    double csig12 = csig1 * csig2 + ssig1 * ssig2;
    J12 *= f1;
    *M12 = csig12 + (e2 * (cbet1sq - cbet2sq) * ssig2 / (w1 + w2)
                    - csig2 * J12) * ssig1 / w1;
    *M21 = csig12 - (e2 * (cbet1sq - cbet2sq) * ssig1 / (w1 + w2)
                    - csig1 * J12) * ssig2 / w2;
  }
}

double Astroid(double x, double y){
  // Solve k^4+2*k^3-(x^2+y^2-1)*k^2-2*y^2*k-y^2 = 0 for positive root k.
  // This solution is adapted from Geocentric::Reverse.
  double k;
  double
    p = sq(x),
    q = sq(y),
    r = (p + q - 1) / 6;
  if ( !(q == 0 && r <= 0) ) {
    double
      // Avoid possible division by zero when r = 0 by multiplying equations
      // for s and t by r^3 and r, resp.
      S = p * q / 4,            // S = r^3 * s
      r2 = sq(r),
      r3 = r * r2,
      // The discrimant of the quadratic equation for T3.  This is zero on
      // the evolute curve p^(1/3)+q^(1/3) = 1
      disc = S * (S + 2 * r3);
    double u = r;
    if (disc >= 0) {
      double T3 = S + r3;
      // Pick the sign on the sqrt to maximize abs(T3).  This minimizes loss
      // of precision due to cancellation.  The result is unchanged because
      // of the way the T is used in definition of u.
      T3 += T3 < 0 ? -sqrt(disc) : sqrt(disc); // T3 = (r * t)^3
      // N.B. cbrt always returns the double root.  cbrt(-8) = -2.
      double T = cbrt(T3); // T = r * t
      // T can be zero; but then r2 / T -> 0.
      u += T + (T != 0 ? r2 / T : 0);
    } else {
      // T is complex, but the way u is defined the result is double.
      double ang = atan2(sqrt(-disc), -(S + r3));
      // There are three possible cube roots.  We choose the root which
      // avoids cancellation.  Note that disc < 0 implies that r < 0.
      u += 2 * r * cos(ang / 3);
    }
    double
      v = sqrt(sq(u) + q),    // guaranteed positive
      // Avoid loss of accuracy when u < 0.
      uv = u < 0 ? q / (v - u) : u + v, // u+v, guaranteed positive
      w = (uv - q) / (2 * v);           // positive?
    // Rearrange expression for k to avoid loss of accuracy due to
    // subtraction.  Division by 0 not possible because uv > 0, w >= 0.
    k = uv / (sqrt(uv + sq(w)) + w);   // guaranteed positive
  } else {               // q == 0 && r <= 0
    // y = 0 with |x| <= 1.  Handle this case directly.
    // for y small, positive root is k = abs(y)/sqrt(1-x^2)
    k = 0;
  }
  return k;
}

double InverseStart(double sbet1, double cbet1,
                                  double sbet2, double cbet2,
                                  double lam12,
                                  double* salp1, double* calp1,
                                  // Only updated if return val >= 0
                                  double* salp2, double* calp2,
                                  // Scratch areas of the right size
                                  double C1a[], double C2a[]){

	// Underflow guard.  We require
	//   tiny_ * epsilon() > 0
	//   tiny_ + epsilon() == epsilon()
//	const double tiny = sqrt(DBL_MIN);
	const double tol0 = DBL_EPSILON;
	// Increase multiplier in defn of tol1_ from 100 to 200 to fix inverse case
	// 52.784459512564 0 -52.784459512563990912 179.634407464943777557
	// which otherwise failed for Visual Studio 10 (Release and Debug)
	const double tol1 = 200 * tol0;
	const double tol2 = sqrt(DBL_EPSILON);
	const double xthresh = 1000 * tol2;


//	const double a = SEMI_MAJOR_AXIS_METERS;
	const double f = FLATTENING;
	const double f1 = 1 - f;
//	const double e2 = f*(2 - f);
	const double e2 = ECCENTRICITY_SQ;
//	const double ep2 = e2/sq(f1);
	const double ep2 = SECOND_ECCENTRICITY_SQ;
//	const double n = f/(2-f);
	const double n = THIRD_FLATTENING;
//	const double b = a*f1;
//	const double c2 = (sq(a) + sq(b) *
//	           (e2 == 0 ? 1 :
//	            (e2 > 0 ? atanh(sqrt(e2)) : atan(sqrt(-e2))) /
//	            sqrt(fabs(e2))))/2; // authalic radius squared
	      // The sig12 threshold for "(double)ly short"
	const double etol2 = 10 * tol2 / fmax(0.1, sqrt(fabs(e2)));

  // Return a starting point for Newton's method in salp1 and calp1 (function
  // value is -1).  If Newton's method doesn't need to be used, return also
  // salp2 and calp2 and function value is sig12.
  double
    sig12 = -1,               // Return value
    // bet12 = bet2 - bet1 in [0, pi); bet12a = bet2 + bet1 in (-pi, 0]
    sbet12 = sbet2 * cbet1 - cbet2 * sbet1,
    cbet12 = cbet2 * cbet1 + sbet2 * sbet1;
#if defined(__GNUC__) && __GNUC__ == 4 && \
(__GNUC_MINOR__ < 6 || defined(__MINGW32__))
  // Volatile declaration needed to fix inverse cases
  // 88.202499451857 0 -88.202499451857 179.981022032992859592
  // 89.262080389218 0 -89.262080389218 179.992207982775375662
  // 89.333123580033 0 -89.333123580032997687 179.99295812360148422
  // which otherwise fail with g++ 4.4.4 x86 -O3 (Linux)
  // and g++ 4.4.0 (mingw) and g++ 4.6.1 (tdm mingw).
  double sbet12a;
  {
    volatile double xx1 = sbet2 * cbet1;
    volatile double xx2 = cbet2 * sbet1;
    sbet12a = xx1 + xx2;
  }
#else
  double sbet12a = sbet2 * cbet1 + cbet2 * sbet1;
#endif
  int shortline = cbet12 >= 0 && sbet12 < (double)(0.5) &&
    lam12 <= M_PI / 6;
  double
    omg12 = (!shortline ? lam12 :
             lam12 / sqrt(1 - e2 * sq((cbet1 + cbet2) / 2))),
    somg12 = sin(omg12), comg12 = cos(omg12);

  *salp1 = cbet2 * somg12;
  *calp1 = comg12 >= 0 ?
    sbet12 + cbet2 * sbet1 * sq(somg12) / (1 + comg12) :
    sbet12a - cbet2 * sbet1 * sq(somg12) / (1 - comg12);

  double
    ssig12 = hypot(*salp1, *calp1),
    csig12 = sbet1 * sbet2 + cbet1 * cbet2 * comg12;

  if (shortline && ssig12 < etol2) {
    // doublely short lines
    *salp2 = cbet1 * somg12;
    *calp2 = sbet12 - cbet1 * sbet2 * sq(somg12) / (1 + comg12);
    SinCosNorm(salp2, calp2);
    // Set return value
    sig12 = atan2(ssig12, csig12);
  } else if (csig12 >= 0 ||
             ssig12 >= 3 * fabs(f) * M_PI * sq(cbet1)) {
    // Nothing to do, zeroth order spherical approximation is OK
  } else {
    // Scale lam12 and bet2 to x, y coordinate system where antipodal point
    // is at origin and singular point is at y = 0, x = -1.
    double y, lamscale, betscale;
    // Volatile declaration needed to fix inverse case
    // 56.320923501171 0 -56.320923501171 179.664747671772880215
    // which otherwise fails with g++ 4.4.4 x86 -O3
    volatile double x;
    if (f >= 0) {            // In fact f == 0 does not get here
      // x = dlong, y = dlat
      {
        double
          k2 = sq(sbet1) * ep2,
          eps = k2 / (2 * (1 + sqrt(1 + k2)) + k2);
        lamscale = f * cbet1 * A3f(eps) * M_PI;
      }
      betscale = lamscale * cbet1;

      x = (lam12 - M_PI) / lamscale;
      y = sbet12a / betscale;
    } else {                  // _f < 0
      // x = dlat, y = dlong
      double
        cbet12a = cbet2 * cbet1 - sbet2 * sbet1,
        bet12a = atan2(sbet12a, cbet12a);
      double m12a, m0, dummy;
      // In the case of lon12 = 180, this repeats a calculation made in
      // Inverse.
      Lengths(n, M_PI + bet12a, sbet1, -cbet1, sbet2, cbet2,
              cbet1, cbet2, &dummy, &m12a, &m0, 0,
              &dummy, &dummy, C1a, C2a);
      x = -1 + m12a/(f1 * cbet1 * cbet2 * m0 * M_PI);
      betscale = x < -(double)(0.01) ? sbet12a / x :
        -f * sq(cbet1) * M_PI;
      lamscale = betscale / cbet1;
      y = (lam12 - M_PI) / lamscale;
    }

    if (y > -tol1 && x > -1 - xthresh) {
      // strip near cut
      if (f >= 0) {
        *salp1 = fmin((double)(1), -(double)(x));
        *calp1 = - sqrt(1 - sq(*salp1));
      } else {
        *calp1 = fmax((double)(x > -tol1 ? 0 : -1), (double)(x));
        *salp1 = sqrt(1 - sq(*calp1));
      }
    } else {
      // Estimate alp1, by solving the astroid problem.
      //
      // Could estimate alpha1 = theta + pi/2, directly, i.e.,
      //   calp1 = y/k; salp1 = -x/(1+k);  for _f >= 0
      //   calp1 = x/(1+k); salp1 = -y/k;  for _f < 0 (need to check)
      //
      // However, it's better to estimate omg12 from astroid and use
      // spherical formula to compute alp1.  This reduces the mean number of
      // Newton iterations for astroid cases from 2.24 (min 0, max 6) to 2.12
      // (min 0 max 5).  The changes in the number of iterations are as
      // follows:
      //
      // change percent
      //    1       5
      //    0      78
      //   -1      16
      //   -2       0.6
      //   -3       0.04
      //   -4       0.002
      //
      // The histogram of iterations is (m = number of iterations estimating
      // alp1 directly, n = number of iterations estimating via omg12, total
      // number of trials = 148605):
      //
      //  iter    m      n
      //    0   148    186
      //    1 13046  13845
      //    2 93315 102225
      //    3 36189  32341
      //    4  5396      7
      //    5   455      1
      //    6    56      0
      //
      // Because omg12 is near pi, estimate work with omg12a = pi - omg12
      double k = Astroid(x, y);
      double
        omg12a = lamscale * ( f >= 0 ? -x * k/(1 + k) : -y * (1 + k)/k ),
        somg12 = sin(omg12a), comg12 = -cos(omg12a);
      // Update spherical estimate of alp1 using omg12 instead of lam12
      *salp1 = cbet2 * somg12;
      *calp1 = sbet12a - cbet2 * sbet1 * sq(somg12) / (1 - comg12);
    }
  }
  SinCosNorm(salp1, calp1);
  return sig12;
}

double Lambda12(double sbet1, double cbet1, double sbet2, double cbet2,
                              double salp1, double calp1,
                              double* salp2, double* calp2,
                              double* sig12,
                              double* ssig1, double* csig1,
                              double* ssig2, double* csig2,
                              double* eps, double* domg12,
                              int diffp, double* dlam12,
                              // Scratch areas of the right size
                              double C1a[], double C2a[], double C3a[]){

	// Underflow guard.  We require
	//   tiny_ * epsilon() > 0
	//   tiny_ + epsilon() == epsilon()
	const double tiny = sqrt(DBL_MIN);
//	const double tol0 = DBL_EPSILON;
	// Increase multiplier in defn of tol1_ from 100 to 200 to fix inverse case
	// 52.784459512564 0 -52.784459512563990912 179.634407464943777557
	// which otherwise failed for Visual Studio 10 (Release and Debug)
//	const double tol1 = 200 * tol0;
//	const double tol2 = sqrt(DBL_EPSILON);
//	const double xthresh = 1000 * tol2;


//	const double a = SEMI_MAJOR_AXIS_METERS;
	const double f = FLATTENING;
//	const double f1 = 1 - f;
//	const double e2 = f*(2 - f);
//	const double ep2 = e2/sq(f1);
	const double e2 = ECCENTRICITY_SQ;
	const double ep2 = SECOND_ECCENTRICITY_SQ;
//	const double n = f/(2-f);
//	const double b = a*f1;
//	const double c2 = (sq(a) + sq(b) *
//	           (e2 == 0 ? 1 :
//	            (e2 > 0 ? atanh(sqrt(e2)) : atan(sqrt(-e2))) /
//	            sqrt(fabs(e2))))/2; // authalic radius squared
	      // The sig12 threshold for "(double)ly short"
//	const double etol2 = 10 * tol2 / fmax(0.1, sqrt(fabs(e2)));

  if (sbet1 == 0 && calp1 == 0)
    // Break degeneracy of equatorial line.  This case has already been
    // handled.
    calp1 = -tiny;

  double
    // sin(alp1) * cos(bet1) = sin(alp0)
    salp0 = salp1 * cbet1,
    calp0 = hypot(calp1, salp1 * sbet1); // calp0 > 0

  double somg1, comg1, somg2, comg2, omg12, lam12;
  // tan(bet1) = tan(sig1) * cos(alp1)
  // tan(omg1) = sin(alp0) * tan(sig1) = tan(omg1)=tan(alp1)*sin(bet1)
  *ssig1 = sbet1;
  somg1 = salp0 * sbet1;
  *csig1 = comg1 = calp1 * cbet1;
  SinCosNorm(ssig1, csig1);
  // SinCosNorm(somg1, comg1); -- don't need to normalize!

  // Enforce symmetries in the case abs(bet2) = -bet1.  Need to be careful
  // about this case, since this can yield singularities in the Newton
  // iteration.
  // sin(alp2) * cos(bet2) = sin(alp0)
  *salp2 = cbet2 != cbet1 ? salp0 / cbet2 : salp1;
  // calp2 = sqrt(1 - sq(salp2))
  //       = sqrt(sq(calp0) - sq(sbet2)) / cbet2
  // and subst for calp0 and rearrange to give (choose positive sqrt
  // to give alp2 in [0, pi/2]).
  *calp2 = cbet2 != cbet1 || fabs(sbet2) != -sbet1 ?
    sqrt(sq(calp1 * cbet1) +
         (cbet1 < -sbet1 ?
          (cbet2 - cbet1) * (cbet1 + cbet2) :
          (sbet1 - sbet2) * (sbet1 + sbet2))) / cbet2 :
    fabs(calp1);
  // tan(bet2) = tan(sig2) * cos(alp2)
  // tan(omg2) = sin(alp0) * tan(sig2).
  *ssig2 = sbet2;
  somg2 = salp0 * sbet2;
  *csig2 = comg2 = *calp2 * cbet2;
  SinCosNorm(ssig2, csig2);
  // SinCosNorm(somg2, comg2); -- don't need to normalize!

  // sig12 = sig2 - sig1, limit to [0, pi]
  *sig12 = atan2(fmax(*csig1 * *ssig2 - *ssig1 * *csig2, (double)(0)),
                *csig1 * *csig2 + *ssig1 * *ssig2);

  // omg12 = omg2 - omg1, limit to [0, pi]
  omg12 = atan2(fmax(comg1 * somg2 - somg1 * comg2, (double)(0)),
                comg1 * comg2 + somg1 * somg2);
  double B312, h0;
  double k2 = sq(calp0) * ep2;
  *eps = k2 / (2 * (1 + sqrt(1 + k2)) + k2);
  C3f(*eps, C3a);
  B312 = (SinCosSeries(1, *ssig2, *csig2, C3a, 6-1) -
          SinCosSeries(1, *ssig1, *csig1, C3a, 6-1));
  h0 = -f * A3f(*eps);
  *domg12 = salp0 * h0 * (*sig12 + B312);
  lam12 = omg12 + *domg12;

  if (diffp) {
    if (*calp2 == 0)
      *dlam12 = - 2 * sqrt(1 - e2 * sq(cbet1)) / sbet1;
    else {
      double dummy;
      Lengths(*eps, *sig12, *ssig1, *csig1, *ssig2, *csig2,
              cbet1, cbet2, &dummy, dlam12, &dummy,
              0, &dummy, &dummy, C1a, C2a);
      *dlam12 /= *calp2 * cbet2;
    }
  }

  return lam12;
}

double InverseKarney(double lat1, double lon1, double lat2, double lon2,
                                 unsigned outmask,
                                 double* s12, double* azi1, double* azi2,
                                 double* m12, double* M12, double* M21, double* S12){

	// Underflow guard.  We require
	//   tiny_ * epsilon() > 0
	//   tiny_ + epsilon() == epsilon()
	const double tiny = sqrt(DBL_MIN);
	const double tol0 = DBL_EPSILON;
	// Increase multiplier in defn of tol1_ from 100 to 200 to fix inverse case
	// 52.784459512564 0 -52.784459512563990912 179.634407464943777557
	// which otherwise failed for Visual Studio 10 (Release and Debug)
	const double tol1 = 200 * tol0;
//	const double tol2 = sqrt(DBL_EPSILON);
//	const double xthresh = 1000 * tol2;

    double DEG2RAD = M_PI / 180.0;
	const double a = SEMI_MAJOR_AXIS_METERS;
	const double f = FLATTENING;
	const double f1 = 1 - f;
//	const double e2 = f*(2 - f);
	const double e2 = ECCENTRICITY_SQ;
//	const double ep2 = e2/sq(f1);
//	const double n = f/(2-f);
	const double n = THIRD_FLATTENING;
//	const double b = a*f1;
	const double b = SEMI_MINOR_AXIS_METERS;
//	const double c2 = (sq(a) + sq(b) *
//	           (e2 == 0 ? 1 :
//	            (e2 > 0 ? atanh(sqrt(e2)) : atan(sqrt(-e2))) /
//	            sqrt(fabs(e2))))/2; // authalic radius squared
	      // The sig12 threshold for "(double)ly short"
//	const double etol2 = 10 * tol2 / fmax(0.1, sqrt(fabs(e2)));


   lon1 = AngNormalize(lon1);
   double lon12 = AngNormalize(AngNormalize(lon2) - lon1);
   // If very close to being on the same meridian, then make it so.
   // Not sure this is necessary...
   lon12 = AngRound(lon12);
   // Make longitude difference positive.
   int lonsign = lon12 >= 0 ? 1 : -1;
   lon12 *= lonsign;
   if (lon12 == 180)
     lonsign = 1;
   // If doubly close to the equator, treat as on equator.
   lat1 = AngRound(lat1);
   lat2 = AngRound(lat2);
   // Swap points so that point with higher (abs) latitude is point 1
   int swapp = fabs(lat1) >= fabs(lat2) ? 1 : -1;
   if (swapp < 0) {
     lonsign *= -1;
     swap(&lat1, &lat2);
   }
   // Make lat1 <= 0
   int latsign = lat1 < 0 ? 1 : -1;
   lat1 *= latsign;
   lat2 *= latsign;
   // Now we have
   //
   //     0 <= lon12 <= 180
   //     -90 <= lat1 <= 0
   //     lat1 <= lat2 <= -lat1
   //
   // longsign, swapp, latsign register the transformation to bring the
   // coordinates to this canonical form.  In all cases, 1 means no change was
   // made.  We make these transformations so that there are few cases to
   // check, e.g., on verifying quadrants in atan2.  In addition, this
   // enforces some symmetries in the results returned.

   double phi, sbet1, cbet1, sbet2, cbet2, s12x, m12x;

   phi = lat1 * DEG2RAD;
   // Ensure cbet1 = +epsilon at poles
   sbet1 = f1 * sin(phi);
   cbet1 = lat1 == -90 ? tiny : cos(phi);
   SinCosNorm(&sbet1, &cbet1);

   phi = lat2 * DEG2RAD;
   // Ensure cbet2 = +epsilon at poles
   sbet2 = f1 * sin(phi);
   cbet2 = fabs(lat2) == 90 ? tiny : cos(phi);
   SinCosNorm(&sbet2, &cbet2);

   // If cbet1 < -sbet1, then cbet2 - cbet1 is a sensitive measure of the
   // |bet1| - |bet2|.  Alternatively (cbet1 >= -sbet1), abs(sbet2) + sbet1 is
   // a better measure.  This logic is used in assigning calp2 in Lambda12.
   // Sometimes these quantities vanish and in that case we force bet2 = +/-
   // bet1 exactly.  An example where is is necessary is the inverse problem
   // 48.522876735459 0 -48.52287673545898293 179.599720456223079643
   // which failed with Visual Studio 10 (Release and Debug)

   if (cbet1 < -sbet1) {
     if (cbet2 == cbet1)
       sbet2 = sbet2 < 0 ? sbet1 : -sbet1;
   } else {
     if (fabs(sbet2) == -sbet1)
       cbet2 = cbet1;
   }

   double
     lam12 = lon12 * DEG2RAD,
     slam12 = lon12 == 180 ? 0 : sin(lam12),
     clam12 = cos(lam12);      // lon12 == 90 isn't interesting

   double a12, sig12, calp1, salp1, calp2, salp2;
   // index zero elements of these arrays are unused
   double C1a[6 + 1], C2a[6 + 1], C3a[6];

   int meridian = lat1 == -90 || slam12 == 0;

   if (meridian) {

     // Endpoints are on a single full meridian, so the geodesic might lie on
     // a meridian.

     calp1 = clam12; salp1 = slam12; // Head to the target longitude
     calp2 = 1; salp2 = 0;           // At the target we're heading north

     double
       // tan(bet) = tan(sig) * cos(alp)
       ssig1 = sbet1, csig1 = calp1 * cbet1,
       ssig2 = sbet2, csig2 = calp2 * cbet2;

     // sig12 = sig2 - sig1
     sig12 = atan2(fmax(csig1 * ssig2 - ssig1 * csig2, (double)(0)),
                   csig1 * csig2 + ssig1 * ssig2);
     {
       double dummy;
       Lengths(n, sig12, ssig1, csig1, ssig2, csig2,
               cbet1, cbet2, &s12x, &m12x, &dummy,
               0, M12, M21, C1a, C2a);
     }
     // Add the check for sig12 since zero length geodesics might yield m12 <
     // 0.  Test case was
     //
     //    echo 20.001 0 20.001 0 | Geod -i
     //
     // In fact, we will have sig12 > pi/2 for meridional geodesic which is
     // not a shortest path.
     if (sig12 < 1 || m12x >= 0) {
       m12x *= a;
       s12x *= b;
       a12 = sig12 / DEG2RAD;
     } else
       // m12 < 0, i.e., prolate and too close to anti-podal
       meridian = 0;
   }

   double omg12;
   if (!meridian &&
       sbet1 == 0 &&   // and sbet2 == 0
       // Mimic the way Lambda12 works with calp1 = 0
       (f <= 0 || lam12 <= M_PI - f * M_PI)) {

     // Geodesic runs along equator
     calp1 = calp2 = 0; salp1 = salp2 = 1;
     s12x = a * lam12;
     m12x = b * sin(lam12 / f1);
     if (0)
       *M12 = *M21 = cos(lam12 / f1);
     a12 = lon12 / f1;
     sig12 = omg12 = lam12 / f1;

   } else if (!meridian) {

     // Now point1 and point2 belong within a hemisphere bounded by a
     // meridian and geodesic is neither meridional or equatorial.

     // Figure a starting point for Newton's method
     sig12 = InverseStart(sbet1, cbet1, sbet2, cbet2,
                          lam12,
                          &salp1, &calp1, &salp2, &calp2,
                          C1a, C2a);

     if (sig12 >= 0) {
       // Short lines (InverseStart sets salp2, calp2)
       double wm = sqrt(1 - e2 * sq((cbet1 + cbet2) / 2));
       s12x = sig12 * a * wm;
       m12x = sq(wm) * a / f1 * sin(sig12 * f1 / wm);
       if (0)
         *M12 = *M21 = cos(sig12 * f1 / wm);
       a12 = sig12 / DEG2RAD;
       omg12 = lam12 / wm;
     } else {

       // Newton's method
       double ssig1, csig1, ssig2, csig2, eps;
       double ov = 0;
       unsigned numit = 0;
       unsigned trip = 0;
       unsigned maxit = 50;
       for (trip = 0; numit < maxit; ++numit) {
         double dv;
         double v = Lambda12(sbet1, cbet1, sbet2, cbet2, salp1, calp1,
                           &salp2, &calp2, &sig12, &ssig1, &csig1, &ssig2, &csig2,
                           &eps, &omg12, trip < 1, &dv, C1a, C2a, C3a) - lam12;
         if (!(fabs(v) > tiny) || !(trip < 1)) {
           if (!(fabs(v) <= fmax(tol1, ov)))
             numit = maxit;
           break;
         }
         double
           dalp1 = -v/dv;
         double
           sdalp1 = sin(dalp1), cdalp1 = cos(dalp1),
           nsalp1 = salp1 * cdalp1 + calp1 * sdalp1;
         calp1 = calp1 * cdalp1 - salp1 * sdalp1;
         salp1 = fmax((double)(0), nsalp1);
         SinCosNorm(&salp1, &calp1);
         // In some regimes we don't get quadratic convergence because slope
         // -> 0.  So use convergence conditions based on epsilon instead of
         // sqrt(epsilon).  The first criterion is a test on abs(v) against
         // 100 * epsilon.  The second takes credit for an anticipated
         // reduction in abs(v) by v/ov (due to the latest update in alp1) and
         // checks this against epsilon.
         if (!(fabs(v) >= tol1 && sq(v) >= ov * tol0)) ++trip;
         ov = fabs(v);
       }

//       if (numit >= maxit) {
//         // Signal failure.
//         if (1)
//           s12 = nan;
//         if (1)
//           azi1 = azi2 = nan;
//         if (1)
//           m12 = nan;
//         if (0)
//           M12 = M21 = nan;
//         if (0)
//           S12 = nan;
//         return -1;
//       }

       {
         double dummy;
         Lengths(eps, sig12, ssig1, csig1, ssig2, csig2,
                 cbet1, cbet2, &s12x, &m12x, &dummy,
                 0, M12, M21, C1a, C2a);
       }
       m12x *= a;
       s12x *= b;
       a12 = sig12 / DEG2RAD;
       omg12 = lam12 - omg12;
     }
   }

   if (1)
     *s12 = 0 + s12x;           // Convert -0 to 0

   if (1)
     *m12 = 0 + m12x;           // Convert -0 to 0

//   if (outmask & AREA) {
//     double
//       // From Lambda12: sin(alp1) * cos(bet1) = sin(alp0)
//       salp0 = salp1 * cbet1,
//       calp0 = hypot(calp1, salp1 * sbet1); // calp0 > 0
//     double alp12;
//     if (calp0 != 0 && salp0 != 0) {
//       double
//         // From Lambda12: tan(bet) = tan(sig) * cos(alp)
//         ssig1 = sbet1, csig1 = calp1 * cbet1,
//         ssig2 = sbet2, csig2 = calp2 * cbet2,
//         k2 = sq(calp0) * _ep2,
//         // Multiplier = a^2 * e^2 * cos(alpha0) * sin(alpha0).
//         A4 = sq(_a) * calp0 * salp0 * _e2;
//       SinCosNorm(ssig1, csig1);
//       SinCosNorm(ssig2, csig2);
//       double C4a[nC4_];
//       C4f(k2, C4a);
//       double
//         B41 = SinCosSeries(false, ssig1, csig1, C4a, nC4_),
//         B42 = SinCosSeries(false, ssig2, csig2, C4a, nC4_);
//       S12 = A4 * (B42 - B41);
//     } else
//       // Avoid problems with indeterminate sig1, sig2 on equator
//       S12 = 0;
//
//     if (!meridian &&
//         omg12 < double(0.75) * M_PI && // Long difference too big
//         sbet2 - sbet1 < double(1.75)) {            // Lat difference too big
//       // Use tan(Gamma/2) = tan(omg12/2)
//       // * (tan(bet1/2)+tan(bet2/2))/(1+tan(bet1/2)*tan(bet2/2))
//       // with tan(x/2) = sin(x)/(1+cos(x))
//       double
//         somg12 = sin(omg12), domg12 = 1 + cos(omg12),
//         dbet1 = 1 + cbet1, dbet2 = 1 + cbet2;
//       alp12 = 2 * atan2( somg12 * ( sbet1 * dbet2 + sbet2 * dbet1 ),
//                          domg12 * ( sbet1 * sbet2 + dbet1 * dbet2 ) );
//     } else {
//       // alp12 = alp2 - alp1, used in atan2 so no need to normalize
//       double
//         salp12 = salp2 * calp1 - calp2 * salp1,
//         calp12 = calp2 * calp1 + salp2 * salp1;
//       // The right thing appears to happen if alp1 = +/-180 and alp2 = 0, viz
//       // salp12 = -0 and alp12 = -180.  However this depends on the sign
//       // being attached to 0 correctly.  The following ensures the correct
//       // behavior.
//       if (salp12 == 0 && calp12 < 0) {
//         salp12 = tiny_ * calp1;
//         calp12 = -1;
//       }
//       alp12 = atan2(salp12, calp12);
//     }
//     S12 += _c2 * alp12;
//     S12 *= swapp * lonsign * latsign;
//     // Convert -0 to 0
//     S12 += 0;
//   }

   // Convert calp, salp to azimuth accounting for lonsign, swapp, latsign.
   if (swapp < 0) {
     swap(&salp1, &salp2);
     swap(&calp1, &calp2);
     if (0)
       swap(M12, M21);
   }

   salp1 *= swapp * lonsign; calp1 *= swapp * latsign;
   salp2 *= swapp * lonsign; calp2 *= swapp * latsign;

   if (1) {
     // minus signs give range [-180, 180). 0- converts -0 to +0.
     *azi1 = 0 - atan2(-salp1, calp1) / DEG2RAD;
     *azi2 = 0 - atan2(-salp2, calp2) / DEG2RAD;
   }

//    Returned value in [0, 180]
   return a12;
}


/*******************************************************************************/
ErrorSet direct(LLPoint origin, double course, double distance, LLPoint* dest,
               double eps)
{
	const double DEG2RAD = M_PI / 180.0;
	double lat1 = origin.latitude/DEG2RAD;
	double lon1 = origin.longitude/DEG2RAD;
	double azi1 = course/DEG2RAD;
	int arcmode = 0;
	double s12_a12 = distance*1852;
	unsigned outmask = 0;
	double lat2;
	double lon2;
	double azi2;
	double s12;
	double m12;
	double M12;
	double M21;
	double S12;
	ErrorSet err = 0;

	DirectKarney(lat1, lon1, azi1, arcmode, s12_a12, outmask,
	                                 &lat2, &lon2, &azi2,
	                                 &s12, &m12, &M12, &M21,
	                                 &S12);

		dest->latitude = lat2*DEG2RAD;
		dest->longitude = lon2*DEG2RAD;

	return err;
}

/*******************************************************************************/
/*
 * Compute destination latitude
 *
 */
ErrorSet directLat(LLPoint origin, double course, double distance, double* lat,
                  double eps)
{
	const double DEG2RAD = M_PI / 180.0;
	double lat1 = origin.latitude/DEG2RAD;
	double lon1 = origin.longitude/DEG2RAD;
	double azi1 = course/DEG2RAD;
	int arcmode = 0;
	double s12_a12 = distance*1852;
	unsigned outmask = 0;
	double lat2;
	double lon2;
	double azi2;
	double s12;
	double m12;
	double M12;
	double M21;
	double S12;
	ErrorSet err = 0;

	DirectKarney(lat1, lon1, azi1, arcmode, s12_a12, outmask,
	                                 &lat2, &lon2, &azi2,
	                                 &s12, &m12, &M12, &M21,
	                                 &S12);

	*lat = lat2*DEG2RAD;

	return err;

}

/*******************************************************************************/
/*
 * Compute destination longitude
 *
 */
ErrorSet directLon(LLPoint origin, double course, double distance, double* lon,
                  double eps)
{
	const double DEG2RAD = M_PI / 180.0;
	double lat1 = origin.latitude/DEG2RAD;
	double lon1 = origin.longitude/DEG2RAD;
	double azi1 = course/DEG2RAD;
	int arcmode = 0;
	double s12_a12 = distance*1852;
	unsigned outmask = 0;
	double lat2;
	double lon2;
	double azi2;
	double s12;
	double m12;
	double M12;
	double M21;
	double S12;
	ErrorSet err = 0;

	DirectKarney(lat1, lon1, azi1, arcmode, s12_a12, outmask,
	                                 &lat2, &lon2, &azi2,
	                                 &s12, &m12, &M12, &M21,
	                                 &S12);

	*lon = lon2*DEG2RAD;

	return err;

}

/*******************************************************************************/
/*
 * Carries out inverse computation, returning course and distance
 *
 */

ErrorSet inverse(LLPoint origin, LLPoint dest, double* crs, double *bcrs,
                  double* dist, double eps)
{
	const double DEG2RAD = M_PI / 180.0;
	double lat1 = origin.latitude/DEG2RAD;
	double lon1 = origin.longitude/DEG2RAD;
	double lat2 = dest.latitude/DEG2RAD;
	double lon2 = dest.longitude/DEG2RAD;
	unsigned outmask = 0;
	double s12, azi1, azi2, m12, M12, M21, S12;
	ErrorSet err = 0;

    if ((origin.latitude == dest.latitude) && (origin.longitude
            == dest.longitude))
    {
    	if( dist != NULL)
    		*dist = 0.0;
        if( crs != NULL)
        	*crs = 0.0;
        if( bcrs != NULL)
        	*bcrs = M_PI;
    }
    else
    {

		InverseKarney(lat1, lon1, lat2, lon2, outmask,
										 &s12, &azi1, &azi2,
										 &m12, &M12, &M21, &S12);
		if( crs != NULL)
			*crs = modcrs(azi1*DEG2RAD);
		if( bcrs != NULL)
			*bcrs = modcrs((azi2 + 180.0)*DEG2RAD);
		if( dist != NULL)
			*dist = s12/1852;

    }
	return err;

}

/*******************************************************************************/
/*
 * Carries out inverse computation but returns only the distance.  Slightly
 * faster than inverse if you don't need the course.
 *
 */

ErrorSet invDist(LLPoint origin, LLPoint dest, double* dist, double eps)
{

    /* Use optional arguments to inverse to simplify this function */
    return inverse(origin, dest, NULL, NULL, dist, eps);

}

/*******************************************************************************/
/*
 * Carries out inverse computation but returns only the course.  Slightly
 * faster than inverse if you don't need to know the distance.
 *
 */

ErrorSet invCrs(LLPoint origin, LLPoint dest, double* fcrs, double* bcrs,
                 double eps)
{

    /* Use optional arguments to inverse to simplify this function */
    return inverse(origin, dest, fcrs, bcrs, NULL, eps);

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

int ptIsOnGeo(LLPoint startPt, LLPoint endPt, LLPoint testPt,
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
        newErr |= inverse(startPt, endPt, &crs12, &crs21, &dist12, eps);
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
            newErr |= direct(startPt, crs12 + M_PI, 1.0, &newStart, eps);
            newErr |= invCrs(newStart, endPt, &crs12, &crs21, eps);

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

                newErr |= direct(startPt, crs12 + M_PI, dist12 + 1.0,
                        &newEnd, eps);
                newErr |= invCrs(startPt, newEnd, &crs12, &crs21, eps);

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
int ptIsOnCrs(LLPoint startPt, double crs12, LLPoint testPt,
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

    newErr |= inverse(startPt, testPt, &crs1Test, crsTest1, dist1Test,eps);
    newErr |= direct(startPt, crs12, *dist1Test, &comparePt, eps);
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
ErrorSet createGeo(Geodesic* geo, LLPoint geoStart, LLPoint geoEnd,
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
    err |= inverse(geoStart, geoEnd, &sCrs, &eCrs, &length, eps);
    geo->startAz = sCrs;
    geo->endAz = modpos(eCrs + M_PI, M_2PI);
    geo->length = length;

    return err;

}

/******************************************************************************
 * Find the course of geodesic at given point. Valid return values are in
 * range [0, 2*PI).  Invalid return value is -1.0.
 */
double geoCrs(Geodesic geo, LLPoint testPt, double* startCrs,
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
    newErr |= invCrs(geo.startPoint, geo.endPoint, &geoCrs, &geoRevCrs,
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
    newErr |= inverse(testPt, geo.startPoint, &crsToStart, &crsFromStart,
            distToPt, eps); //test to start, return distToPt
    if (newErr)
    {
        *err |= newErr;
        return -1.0;
    }
    newErr |= inverse(testPt, geo.endPoint, &crsToEnd, &crsFromEnd,
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
    else if (ptIsOnCrs(geo.startPoint, geoCrs, testPt, &crsToStart,
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
    else if (ptIsOnCrs(geo.endPoint, geoRevCrs, testPt, &crsToEnd,
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



