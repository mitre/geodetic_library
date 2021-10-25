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

/*
 * This is a direct translation of the NGS Forward and Inverse algorithms from FORTRAN into C.
 *
 * FORTRAN code can be downloaded from NGS at http://www.ngs.noaa.gov/TOOLS/Inv_Fwd/Inv_Fwd.html
 *
 * gpnarc, gpnloa, and gpnhri are from version 2
 * invers is from version 3
 *
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

void DIRCT1(double GLAT1, double GLON1, double* GLAT2, double* GLON2, double FAZ, double* BAZ, double S){
//C
//C *** SOLUTION OF THE GEODETIC DIRECT PROBLEM AFTER T.VINCENTY
//C *** MODIFIED RAINSFORD'S METHOD WITH HELMERT'S ELLIPTICAL TERMS
//C *** EFFECTIVE IN ANY AZIMUTH AND AT ANY DISTANCE SHORT OF ANTIPODAL
//C
//C *** A IS THE SEMI-MAJOR AXIS OF THE REFERENCE ELLIPSOID
//C *** F IS THE FLATTENING OF THE REFERENCE ELLIPSOID
//C *** LATITUDES AND LONGITUDES IN RADIANS POSITIVE NORTH AND EAST
//C *** AZIMUTHS IN RADIANS CLOCKWISE FROM NORTH
//C *** GEODESIC DISTANCE S ASSUMED IN UNITS OF SEMI-MAJOR AXIS A
//C
//C *** PROGRAMMED FOR CDC-6600 BY LCDR L.PFEIFER NGS ROCKVILLE MD 20FEB75
//C *** MODIFIED FOR SYSTEM 360 BY JOHN G GERGEN NGS ROCKVILLE MD 750608
//C

//    IMPLICIT REAL*8 (A-H,O-Z)
	double R, TU, SF, CF;
	double CU, SU, SA, C2A;
	double X, C, D, Y;
	double SY, CY, CZ, E;

//    COMMON/CONST/PI,RAD
	double PI = 4.0*atan(1.0); //	pi=4.d0*datan(1.d0)
//	double RAD = 180.0/PI; //	rad=180.d0/pi

//    COMMON/ELIPSOID/A,F
	double A = 6378137.0; //	a=6378137.d0
//	double F = 1.0/298.25722210088; //	f=1.d0/298.25722210088d0 this is the GRS80 flattening
	double F = FLATTENING; //	WGS84 flattening

//    DATA EPS/0.5D-13/
	double EPS = 0.5e-13;

//    R=1.-F
//    TU=R*DSIN(GLAT1)/DCOS(GLAT1)
//    SF=DSIN(FAZ)
//    CF=DCOS(FAZ)
	R = 1.0 - F;
	TU = R*sin(GLAT1)/cos(GLAT1);
	SF = sin(FAZ);
	CF = cos(FAZ);

//    BAZ=0.
//    IF(CF.NE.0.) BAZ=DATAN2(TU,CF)*2.
	*BAZ = 0.0;
	if(CF != 0.0) *BAZ = atan2(TU,CF)*2.0;

//    CU=1./DSQRT(TU*TU+1.)
//    SU=TU*CU
//    SA=CU*SF
//    C2A=-SA*SA+1.
	CU = 1.0/sqrt(TU*TU + 1);
	SU = TU*CU;
	SA = CU*SF;
	C2A = -SA*SA + 1.0;

//    X=DSQRT((1./R/R-1.)*C2A+1.)+1.
//    X=(X-2.)/X
//    C=1.-X
//    C=(X*X/4.+1)/C
//    D=(0.375D0*X*X-1.)*X
//    TU=S/R/A/C
//    Y=TU
	X = sqrt((1.0/R/R - 1.0)*C2A + 1.0) + 1.0;
	X = (X - 2.0)/X;
	C = 1.0-X;
	C = (X*X/4.0 + 1)/C;
	D = (0.375*X*X - 1.0)*X;
	TU = S/R/A/C;
	Y = TU;

	//100 SY=DSIN(Y)
	//    CY=DCOS(Y)
	//    CZ=DCOS(BAZ+Y)
	//    E=CZ*CZ*2.-1.
	//    C=Y
	//    X=E*CY
	//    Y=E+E-1.
	//    Y=(((SY*SY*4.-3.)*Y*CZ*D/6.+X)*D/4.-CZ)*SY*D+TU
	//    IF(DABS(Y-C).GT.EPS)GO TO 100
	while(fabs(Y-C) > EPS){
		SY = sin(Y);
		CY = cos(Y);
		CZ = cos(*BAZ+Y);
		E = CZ*CZ*2.0 - 1.0;
		C = Y;
		X = E*CY;
		Y = E + E - 1.0;
		Y = (((SY*SY*4.0 - 3.0)*Y*CZ*D/6.0 + X)*D/4.0 - CZ)*SY*D + TU;
	}

//    BAZ=CU*CY*CF-SU*SY
//    C=R*DSQRT(SA*SA+BAZ*BAZ)
//    D=SU*CY+CU*SY*CF
//    GLAT2=DATAN2(D,C)
	*BAZ = CU*CY*CF - SU*SY;
	C = R*sqrt(SA*SA + (*BAZ)*(*BAZ));
	D = SU*CY + CU*SY*CF;
	*GLAT2 = atan2(D,C);

//    C=CU*CY-SU*SY*CF
//    X=DATAN2(SY*SF,C)
//    C=((-3.*C2A+4.)*F+4.)*C2A*F/16.
//    D=((E*CY*C+CZ)*SY*C+Y)*SA
//    GLON2=GLON1+X-(1.-C)*D*F
	C = CU*CY - SU*SY*CF;
	X = atan2(SY*SF,C);
	C = ((-3.0*C2A + 4.0)*F + 4.0)*C2A*F/16.0;
	D = ((E*CY*C + CZ)*SY*C + Y)*SA;
	*GLON2 = GLON1 + X - (1.0 - C)*D*F;

//    BAZ=DATAN2(SA,BAZ)+PI
	*BAZ = atan2(SA,*BAZ) + PI;
}


/*******************************************************************************/
ErrorSet direct(LLPoint origin, double course, double distance, LLPoint* dest,
               double eps)
{

    ErrorSet err = 0;
    double GLAT1, GLON1, GLAT2, GLON2, FAZ, BAZ, S;

    GLAT1 = origin.latitude;
    GLON1 = origin.longitude;
    FAZ = course;
    S = distance*1852;

    DIRCT1(GLAT1, GLON1, &GLAT2, &GLON2, FAZ, &BAZ, S);

    dest -> latitude = GLAT2;
    dest -> longitude = GLON2;

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

    ErrorSet err = 0;
    LLPoint dest;

    err = direct(origin, course, distance, &dest ,eps);

    *lat = dest.latitude;

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

    ErrorSet err = 0;
    LLPoint dest;

    err = direct(origin, course, distance, &dest ,eps);

    *lon = dest.longitude;

    return err;

}


void gpnarc(double AMAX, double FLAT, double ESQ, double PI, double P1, double P2, double* ARC){
//C
//C********1*********2*********3*********4*********5*********6*********7*
//C
//C NAME:        GPNARC
//C VERSION:     200005.26
//C WRITTEN BY:  ROBERT (Sid) SAFFORD
//C PURPOSE:     SUBROUTINE TO COMPUTE THE LENGTH OF A MERIDIONAL ARC
//C              BETWEEN TWO LATITUDES
//C
//C INPUT PARAMETERS:
//C -----------------
//C AMAX         SEMI-MAJOR AXIS OF REFERENCE ELLIPSOID
//C FLAT         FLATTENING (0.0033528 ... )
//C ESQ          ECCENTRICITY SQUARED FOR REFERENCE ELLIPSOID
//C PI           3.14159...
//C P1           LAT STATION 1
//C P2           LAT STATION 2
//C
//C OUTPUT PARAMETERS:
//C ------------------
//C ARC          GEODETIC DISTANCE
//C
//C LOCAL VARIABLES AND CONSTANTS:
//C ------------------------------
//C GLOBAL VARIABLES AND CONSTANTS:
//C -------------------------------
//C
//C    MODULE CALLED BY:    GENERAL
//C
//C    THIS MODULE CALLS:
//C       LLIBFORE/ OPEN,   CLOSE,  READ,   WRITE,  INQUIRE
//C                 DABS,   DBLE,   FLOAT,  IABS,   CHAR,   ICHAR
//C
//C    INCLUDE FILES USED:
//C    COMMON BLOCKS USED:
//C
//C    REFERENCES: Microsoft FORTRAN 4.10 Optimizing Compiler, 1988
//C                MS-DOS Operating System
//C    COMMENTS:
//C********1*********2*********3*********4*********5*********6*********7*
//C::MODIFICATION HISTORY
//C::197507.05, RWS, VER 00 TENCOL RELEASED FOR FIELD USE
//C::198311.20, RWS, VER 01 MTEN   RELEASED TO FIELD
//C::198411.26, RWS, VER 07 MTEN2  RELEASED TO FIELD
//C::1985xx.xx, RWS, CODE   CREATED
//C::198506.10, RWS, WRK    ENHANCEMENTS RELEASED TO FIELD
//C::198509.01, RWS, VER 11 MTEN3  RELEASED TO FIELD
//C::198512.18, RWS, CODE   MODIFIED FOR MTEN3
//C::198708.10, RWS, CODE   MODIFIED TO USE NEW MTEN4 GPN RECORD FORMAT
//C::199112.31, RWS, VER 20 MTEN4 RELEASED TO FIELD
//C::200001.13, RWS, VER 21 MTEN4 RELEASED TO FIELD
//C::200005.26, RWS, CODE   RESTRUCTURED & DOCUMENTATION ADDED
//C::200012.31, RWS, VER 23 MTEN5 RELEASED
//C********1*********2*********3*********4*********5*********6*********7*
//CE::GPNARC
//C ---------------------------
//C     M T E N  (VERSION 3)
//C     M T E N  (VERSION 5.23)
//C ---------------------------

//      IMPLICIT REAL*8 (A-H,O-Z)

//      LOGICAL  FLAG
	int FLAG;
//      DATA TT/5.0D-15/
	double TT = 5.0e-15;
//C     CHECK FOR A 90 DEGREE LOOKUP

//      FLAG = .FALSE.
	FLAG = 0;
//      S1 = DABS(P1)
//      S2 = DABS(P2)
	double S1 = fabs(P1);
	double S2 = fabs(P2);

//      IF( (PI/2.0D0-TT).LT.S2 .AND. S2.LT.(PI/2.0D0+TT) )THEN
	if(((PI/2.0 - TT) < S2) && (S2 < (PI/2.0 + TT))){
//        FLAG = .TRUE.
		FLAG = 1;
	}//      ENDIF

//      IF( S1.GT.TT )THEN
	if(S1 > TT){
//        FLAG = .FALSE.
		FLAG = 0;
	}//      ENDIF

//      DA = (P2-P1)
//      S1 = 0.0D0
//      S2 = 0.0D0
	double DA = (P2 - P1);
	S1 = 0.0;
	S2 = 0.0;

//C     COMPUTE THE LENGTH OF A MERIDIONAL ARC BETWEEN TWO LATITUDES

//      E2 = ESQ
//      E4 = E2*E2
//      E6 = E4*E2
//      E8 = E6*E2
//      EX = E8*E2
	double E2 = ESQ;
	double E4 = E2*E2;
	double E6 = E4*E2;
	double E8 = E6*E2;
	double EX = E8*E2;

//      T1 = E2*(003.0D0/4.0D0)
//      T2 = E4*(015.0D0/64.0D0)
//      T3 = E6*(035.0D0/512.0D0)
//      T4 = E8*(315.0D0/16384.0D0)
//      T5 = EX*(693.0D0/131072.0D0)
	double T1 = E2*(3.0/4.0);
	double T2 = E4*(15.0/64.0);
	double T3 = E6*(35.0/512.0);
	double T4 = E8*(315.0/16384.0);
	double T5 = EX*(693.0/131072.0);

//      A  = 1.0D0+T1+3.0D0*T2+10.0D0*T3+35.0D0*T4+126.0D0*T5
	double A  = 1.0 + T1 + 3.0*T2 + 10.0*T3 + 35.0*T4 + 126.0*T5;

	double B, C, D, E, F;
	double DB, DC, DD, DE, DF;
//      IF( FLAG )THEN
	if( FLAG ){
//        GOTO 1
	}//      ENDIF
	else {
//      B  = T1+4.0D0*T2+15.0D0*T3+56.0D0*T4+210.0D0*T5
//      C  = T2+06.0D0*T3+28.0D0*T4+120.0D0*T5
//      D  = T3+08.0D0*T4+045.0D0*T5
//      E  = T4+010.0D0*T5
//      F  = T5
	B  = T1 + 4.0*T2 + 15.0*T3 + 56.0*T4 + 210.0*T5;
	C  = T2 + 6.0*T3 + 28.0*T4 + 120.0*T5;
	D  = T3 + 8.0*T4 + 45.0*T5;
	E  = T4 + 10.0*T5;
	F  = T5;

//      DB = DSIN(P2*2.0D0)-DSIN(P1*2.0D0)
//      DC = DSIN(P2*4.0D0)-DSIN(P1*4.0D0)
//      DD = DSIN(P2*6.0D0)-DSIN(P1*6.0D0)
//      DE = DSIN(P2*8.0D0)-DSIN(P1*8.0D0)
//      DF = DSIN(P2*10.0D0)-DSIN(P1*10.0D0)
	DB = sin(P2*2.0) - sin(P1*2.0);
	DC = sin(P2*4.0) - sin(P1*4.0);
	DD = sin(P2*6.0) - sin(P1*6.0);
	DE = sin(P2*8.0) - sin(P1*8.0);
	DF = sin(P2*10.0) - sin(P1*10.0);

//C     COMPUTE THE S2 PART OF THE SERIES EXPANSION

//      S2 = -DB*B/2.0D0+DC*C/4.0D0-DD*D/6.0D0+DE*E/8.0D0-DF*F/10.0D0
	S2 = -DB*B/2.0 + DC*C/4.0 - DD*D/6.0 + DE*E/8.0 - DF*F/10.0;

	}
//C     COMPUTE THE S1 PART OF THE SERIES EXPANSION
//    1 S1 = DA*A
	S1 = DA*A;

//C     COMPUTE THE ARC LENGTH

//      ARC = AMAX*(1.0D0-ESQ)*(S1+S2)
	*ARC = AMAX*(1.0 - ESQ)*(S1 + S2);

//      RETURN
}


void gpnloa(double AMAX, double FLAT, double ESQ, double PI, double DL, double* AZ1, double* AZ2, double AO, double BO, double* SMS){
//C
//C********1*********2*********3*********4*********5*********6*********7*
//C
//C NAME:        GPNLOA
//C VERSION:     200005.26
//C WRITTEN BY:  ROBERT (Sid) SAFFORD
//C PURPOSE:     SUBROUTINE TO COMPUTE THE LIFF-OFF-AZIMUTH CONSTANTS
//C
//C INPUT PARAMETERS:
//C -----------------
//C AMAX         SEMI-MAJOR AXIS OF REFERENCE ELLIPSOID
//C FLAT         FLATTENING (0.0033528 ... )
//C ESQ          ECCENTRICITY SQUARED FOR REFERENCE ELLIPSOID
//C PI           3.14159...
//C DL           LON DIFFERENCE
//C AZ1          AZI AT STA 1 -> STA 2
//C
//C OUTPUT PARAMETERS:
//C ------------------
//C AZ2          AZ2 AT STA 2 -> STA 1
//C AO           CONST
//C BO           CONST
//C SMS          DISTANCE ... EQUATORIAL - GEODESIC  (S - s)   "SMS"
//C
//C LOCAL VARIABLES AND CONSTANTS:
//C ------------------------------
//C GLOBAL VARIABLES AND CONSTANTS:
//C -------------------------------
//C
//C    MODULE CALLED BY:    GENERAL
//C
//C    THIS MODULE CALLS:
//C       LLIBFORE/ DSIN,   DCOS,   DABS,   DASIN
//C
//C    INCLUDE FILES USED:
//C    COMMON BLOCKS USED:
//C
//C    REFERENCES: Microsoft FORTRAN 4.10 Optimizing Compiler, 1988
//C                MS-DOS Operating System
//C    COMMENTS:
//C********1*********2*********3*********4*********5*********6*********7*
//C::MODIFICATION HISTORY
//C::1985xx.xx, RWS, CODE   CREATED
//C::198506.10, RWS, WRK    ENHANCEMENTS RELEASED TO FIELD
//C::198509.01, RWS, VER 11 MTEN3  RELEASED TO FIELD
//C::198512.18, RWS, CODE   MODIFIED FOR MTEN3
//C::198708.10, RWS, CODE   MODIFIED TO USE NEW MTEN4 GPN RECORD FORMAT
//C::199112.31, RWS, VER 20 MTEN4 RELEASED TO FIELD
//C::200001.13, RWS, VER 21 MTEN4 RELEASED TO FIELD
//C::200005.26, RWS, CODE   RESTRUCTURED & DOCUMENTATION ADDED
//C::200012.31, RWS, VER 23 MTEN5 RELEASED
//C********1*********2*********3*********4*********5*********6*********7*
//CE::GPNLOA
//C ---------------------------
//C     M T E N  (VERSION 3)
//C              (VERSION 4.22)
//C              (VERSION 5.23)
//C ---------------------------

//      IMPLICIT REAL*8 (A-H,O-Z)

//      DATA TT/5.0D-13/
	double TT = 5.0e-13;
//      DLON = DABS(DL)
	double DLON = fabs(DL);
//      CONS = (PI-DLON)/(PI*FLAT)
	double CONS = (PI - DLON)/(PI*FLAT);
//      F    = FLAT
	double F = FLAT;

//C     COMPUTE AN APPROXIMATE AZ

//      AZ   = DASIN(CONS)
	double AZ = asin(CONS);

//      T1   =    1.0D0
//      T2   =  (-1.0D0/4.0D0)*F*(1.0D0+F+F*F)
//      T4   =    3.0D0/16.0D0*F*F*(1.0D0+(9.0D0/4.0D0)*F)
//      T6   = (-25.0D0/128.0D0)*F*F*F
	double T1   =    1.0;
	double T2   =  (-1.0/4.0)*F*(1.0 + F + F*F);
	double T4   =    3.0/16.0*F*F*(1.0 + (9.0/4.0)*F);
	double T6   = (-25.0/128.0)*F*F*F;

//      ITER = 0
	double ITER = 0;

	double S, C2, CS;
	while((fabs(S-AZ) > TT) && ITER < 6){
//    1 ITER = ITER+1
//      S    = DCOS(AZ)
		S = cos(AZ);
//      C2   = S*S
		C2 = S*S;

//C     COMPUTE NEW AO

//      AO   = T1 + T2*C2 + T4*C2*C2 + T6*C2*C2*C2
		AO   = T1 + T2*C2 + T4*C2*C2 + T6*C2*C2*C2;
//      CS   = CONS/AO
		CS   = CONS/AO;
//      S    = DASIN(CS)
		S    = asin(CS);

//      IF( DABS(S-AZ).LT.TT )THEN
//        GOTO 2
//      ENDIF

//      AZ   = S
		AZ = S;
//      IF( ITER.LE.6 )THEN
//        GOTO 1
//      ENDIF
	}
//    2 AZ1  = S
	*AZ1 = S;

//      IF( DL.LT.0.0D0 )THEN
	if(DL < 0.0){
//        AZ1 = 2.0D0*PI-AZ1
		*AZ1 = 2.0*PI - *AZ1;
	}//      ENDIF

//      AZ2  = 2.0D0*PI-AZ1
	*AZ2  = 2.0*PI - *AZ1;

//C     EQUATORIAL - GEODESIC  (S - s)   "SMS"

//      ESQP = ESQ/(1.0D0-ESQ)
	double ESQP = ESQ/(1.0 - ESQ);
//      S    = DCOS(AZ1)
	S = cos(*AZ1);

//      U2   = ESQP*S*S
//      U4   = U2*U2
//      U6   = U4*U2
//      U8   = U6*U2
	double U2   = ESQP*S*S;
	double U4   = U2*U2;
	double U6   = U4*U2;
	double U8   = U6*U2;

//      T1   =     1.0D0
//      T2   =    (1.0D0/4.0D0)*U2
//      T4   =   (-3.0D0/64.0D0)*U4
//      T6   =    (5.0D0/256.0D0)*U6
//      T8   = (-175.0D0/16384.0D0)*U8
	T1   =     1.0;
	T2   =    (1.0/4.0)*U2;
	T4   =   (-3.0/64.0)*U4;
	T6   =    (5.0/256.0)*U6;
	double T8   = (-175.0/16384.0)*U8;


//      BO   = T1 + T2 + T4 + T6 + T8
	BO   = T1 + T2 + T4 + T6 + T8;
//      S    = DSIN(AZ1)
	S    = sin(*AZ1);
//      SMS  = AMAX*PI*(1.0D0 - FLAT*DABS(S)*AO - BO*(1.0D0-FLAT))
	*SMS  = AMAX*PI*(1.0 - FLAT*fabs(S)*AO - BO*(1.0 - FLAT));

//      RETURN
}

void gpnhri (double a, double f, double esq, double pi, double p1, double e1, double p2, double e2, double* az1, double* az2, double* s){
//c
//c********1*********2*********3*********4*********5*********6*********7*
//c
//c name:        gpnhri
//c version:     200208.09
//c written by:  robert (sid) safford
//c purpose:     subroutine to compute helmert rainsford inverse problem
//c
//c     solution of the geodetic inverse problem after t. vincenty
//c     modified rainsford's method with helmert's elliptical terms
//c     effective in any azimuth and at any distance short of antipocal
//c     from/to stations must not be the geographic pole.
//c     parameter a is the semi-major axis of the reference ellipsoid
//c     finv=1/f is the inverse flattening of the reference ellipsoid
//c     latitudes and longitudes in radians positive north and west
//c     forward and back azimuths returned in radians clockwise from south
//c     geodesic distance s returned in units of semi-major axis a
//c     programmed for ibm 360-195   09/23/75
//c
//c     note - note - note -
//c     1. do not use for meridional arcs and be careful on the equator.
//c     2. azimuths are from north(+) clockwise and
//c     3. longitudes are positive east(+)
//c
//c input parameters:
//c -----------------
//c a            semi-major axis of reference ellipsoid      meters
//c f            flattening (0.0033528...)
//c esq          eccentricity squared
//c pi           3.14159...
//c p1           lat station 1                               radians
//c e1           lon station 1                               radians
//c p2           lat station 2                               radians
//c e2           lon station 2                               radians
//c
//c output parameters:
//c ------------------
//c az1          azi at sta 1 -> sta 2                       radians
//c az2          azi at sta 2 -> sta 1                       radians
//c s            geodetic dist between sta(s) 1 & 2          meters
//c
//c local variables and constants:
//c ------------------------------
//c aa               constant from subroutine gpnloa
//c alimit           equatorial arc distance along the equator   (radians)
//c arc              meridional arc distance latitude p1 to p2 (in meters)
//c az1              azimuth forward                          (in radians)
//c az2              azimuth back                             (in radians)
//c bb               constant from subroutine gpnloa
//c dlon             temporary value for difference in longitude (radians)
//c equ              equatorial distance                       (in meters)
//c r1,r2            temporary variables
//c s                ellipsoid distance                        (in meters)
//c sms              equatorial - geodesic distance (S - s) "Sms"
//c ss               temporary variable
//c tol0             tolerance for checking computation value
//c tol1             tolerance for checking a real zero value
//c tol2             tolerance for close to zero value
//c twopi            two times constant pi
//c
//c global variables and constants:
//c -------------------------------
//c
//c    module called by:    general
//c
//c    this module calls:   gpnarc, gpnloa
//c       llibfore/ dsin,   dcos,   dsqrt,  dabs,  datan2, write
//c
//c    include files used:
//c    common blocks used:
//c
//c    references: microsoft fortran 4.10 optimizing compiler, 1988
//c                ms-dos operating system
//c    comments:
//c********1*********2*********3*********4*********5*********6*********7*
//c::modification history
//c::197507.05, rws, ver 00 tencol released for field use
//c::198311.20, rws, ver 01 mten   released to field
//c::198411.26, rws, ver 07 mten2  released to field
//c::198506.10, rws, wrk    enhancements released to field
//c::198507.22, rws, code   modified for mten3
//c::198509.01, rws, ver 11 mten3  released to field
//c::198708.10, rws, code   modified to use new mten4 gpn record format
//c::199112.31, rws, ver 20 mten4 released to field
//c::200001.13, rws, ver 21 mten4 released to field
//c::200005.26, rws, code   restructured & documentation added
//c::200012.31, rws, ver 23 mten5 released
//c::200104.09, rws, code   added to calblin program
//c::200208.09, rws, code   added subroutines gpnarc & gpnloa
//c********1*********2*********3*********4*********5*********6*********7*
//ce::gpnhri
//c  -------------------------------
//c     m t e n  (version 3)
//c              (version 4.22)
//c              (version 5.23)
//c  -------------------------------

// implicit real*8 (a-h,o-z)


	// data tol0 /5.0d-15/
	// data tol1 /5.0d-14/
	// data tol2 /7.0d-03/
	double tol0 = 5.0e-15;
	double tol1 = 5.0e-14;
	double tol2 = 7.0e-03;

	// twopi = 2.0d0*pi
	double twopi = 2.0*pi;

//c     test the longitude difference with tol1
//c     tol1 is approximately 0.000000001 arc seconds

	// ss = e2-e1
	double ss = e2 - e1;

	double r2, r1;
	double arc;

	// if( dabs(ss).lt.tol1 )then
	if(fabs(ss) < tol1){
		// e2 = e2+tol1
		e2 = e2 + tol1;

		// write(*,*) ' longitudal difference is near zero '

		// r2 = p2
		// r1 = p1
		r2 = p2;
		r1 = p1;

		// call gpnarc ( a, f, esq, pi, r1, r2, arc )
		gpnarc(a, f, esq, pi, r1, r2, &arc);

		// s  = dabs( arc )
		*s  = fabs(arc);

		// if( p2.gt.p1 )then
		if(p2 > p1){
			// az1 = 0.0d0
			// az2 = pi
			*az1 = 0.0;
			*az2 = pi;
		} else {
			// az1 = pi
			// az2 = 0.0d0
			*az1 = pi;
			*az2 = 0.0;
		}		// endif

		return;	// return
	}	// endif

//c     test for longitude over 180 degrees

	// dlon = e2-e1
	double dlon = e2 - e1;

	// if( dlon.ge.0.0d0 )then
	if(dlon >= 0.0){

		// if( pi.le.dlon .and. dlon.lt.twopi )then
		if((pi <= dlon) && (dlon < twopi)){
			// dlon = dlon-twopi
			dlon = dlon - twopi;

		}	// endif

	} else {
		// ss = dabs(dlon)
		ss = fabs(dlon);

		// if( pi.le.ss .and. ss.lt.twopi )then
		if((pi <= ss) && (ss < twopi)){
			// dlon = dlon+twopi
			dlon = dlon + twopi;
		}	// endif

	}	// endif


	// ss = dabs( dlon )
	ss = fabs(dlon);

	// if( ss.gt.pi )then
	if(ss > pi ){
//c::     write(*,*) '  '
//c::     write(*,*) ' Longitude difference over 180 degrees  '
//c::     write(*,*) ' Turn it around '

		// ss = twopi-ss
		ss = twopi - ss;
	}//      endif

//c     compute the limit in longitude (alimit), it is equal
//c     to twice the distance from the equator to the pole,
//c     as measured along the equator (east/ewst)

	// alimit = pi*(1.0d0-f)
	double alimit = pi*(1.0 - f);

//c     test for anti-nodal difference

	// if( ss.ge.alimit )then
	if( ss >= alimit ){
		// r1 = dabs(p1)
		// r2 = dabs(p2)
		r1 = fabs(p1);
		r2 = fabs(p2);

//c       latitudes r1 & r2 are not near the equator

		// if( r1.gt.tol2 .and. r2.gt.tol2 )then
		if((r1 > tol2) && (r2 > tol2)){
			//  goto 60

		}//        endif
		else {

//c       longitude difference is greater than lift-off point
//c       now check to see if  "both"  r1 & r2 are on equator

			// if( r1.lt.tol1 .and. r2.gt.tol2 )then
			if((r1 < tol1) && (r2 > tol2)){
				// goto 60
			}//        endif
			else {

				// if( r2.lt.tol1 .and. r1.gt.tol2 )then
				if((r2 < tol1) && (r1 > tol2)){
					// goto 60
				}//        endif
				else {

//c       check for either r1 or r2 just off the equator but < tol2

					// if( r1.gt.tol1. or. r2.gt.tol1 )then
					if((r1 > tol1) || (r2 > tol1)){
						// az1 = 0.0d0
						// az2 = 0.0d0
						// s   = 0.0d0
				          *az1 = 0.0;
				          *az2 = 0.0;
				          *s   = 0.0;
				          // return
						return;
					}//        endif

//c       compute the azimuth to anti-nodal point

//c::     write(*,*) '  '
//c::     write(*,*) ' Longitude difference beyond lift-off point '
//c::     write(*,*) '  '

					// call gpnloa (a,f,esq,pi,dlon,az1,az2,aa,bb,sms)
					double aa, bb, sms;
					gpnloa(a, f, esq, pi, dlon, az1, az2, aa, bb, &sms);
//c       compute the equatorial distance & geodetic

					// equ = a*dabs(dlon)
					double equ = a*fabs(dlon);
					// s   = equ-sms
					*s = equ-sms;

					// return
					return;
				}
			}
		}
	}//      endif
	// 60 continue

	// f0   = (1.0d0-f)
	double f0 = (1.0 - f);
	// b    = a*f0
	double b = a*f0;
	// epsq = esq/(1.0d0-esq)
	double epsq = esq/(1.0 - esq);
	// f2   = f*f
	double f2 = f*f;
	// f3   = f*f2
	double f3 = f*f2;
	// f4   = f*f3
	double f4 = f*f3;

//c     the longitude difference

	// dlon  = e2-e1
	dlon = e2 - e1;
	// ab    = dlon
	double ab = dlon;
	// kount = 0
	int kount = 0;

//c     the reduced latitudes

	// u1    = f0*dsin(p1)/dcos(p1)
	// u2    = f0*dsin(p2)/dcos(p2)
	double u1 = f0*sin(p1)/cos(p1);
	double u2 = f0*sin(p2)/cos(p2);

	// u1    = datan(u1)
	// u2    = datan(u2)
	u1 = atan(u1);
	u2 = atan(u2);

	// su1   = dsin(u1)
	// cu1   = dcos(u1)
	double su1   = sin(u1);
	double cu1   = cos(u1);

	// su2   = dsin(u2)
	// cu2   = dcos(u2)
	double su2   = sin(u2);
	double cu2   = cos(u2);

	double xy = 1;
	double clon, slon;
	double csig, ssig;
	double sig, sinalf;
	double w, t4, t6;
	double ao, a2, a4, a6;
	double qo;
	double q2, q4, q6;
	double r3;
	double xz;

	while((kount <= 7) && (xy > 0.5e-13)){
//c     counter for the iteration operation
		// 1 kount = kount+1
		kount = kount + 1;

		// clon  = dcos(ab)
		// slon  = dsin(ab)
		clon = cos(ab);
		slon = sin(ab);

		// csig  = su1*su2+cu1*cu2*clon
		// ssig  = dsqrt((slon*cu2)**2+(su2*cu1-su1*cu2*clon)**2)
		csig = su1*su2+cu1*cu2*clon;
		ssig = sqrt(pow((slon*cu2),2) + pow((su2*cu1-su1*cu2*clon),2));

		// sig   = datan2(ssig,csig)
		// sinalf=cu1*cu2*slon/ssig
		sig = atan2(ssig, csig);
		sinalf = cu1*cu2*slon/ssig;

		// w   = (1.0d0-sinalf*sinalf)
		// t4  = w*w
		// t6  = w*t4
		w = (1.0-sinalf*sinalf);
		t4 = w*w;
		t6 = w*t4;

//c     the coefficients of type a

		// ao  = f-f2*(1.0d0+f+f2)*w/4.0d0+3.0d0*f3*(1.0d0+
		// 1        9.0d0*f/4.0d0)*t4/16.0d0-25.0d0*f4*t6/128.0d0
		ao = f - f2*(1.0 + f + f2)*w/4.0 + 3.0*f3*(1.0 +
				9.0*f/4.0)*t4/16.0 - 25.0*f4*t6/128.0;

		// a2  = f2*(1.0d0+f+f2)*w/4.0d0-f3*(1.0d0+9.0d0*f/4.0d0)*t4/4.0d0+
		// 1        75.0d0*f4*t6/256.0d0
	    a2 = f2*(1.0 + f + f2)*w/4.0 - f3*(1.0 + 9.0*f/4.0)*t4/4.0 +
	    		75.0*f4*t6/256.0;

	    // a4  = f3*(1.0d0+9.0d0*f/4.0d0)*t4/32.0d0-15.0d0*f4*t6/256.0d0
	    a4 = f3*(1.0 + 9.0*f/4.0)*t4/32.0 - 15.0*f4*t6/256.0;

	    // a6  = 5.0d0*f4*t6/768.0d0
	    a6 = 5.0*f4*t6/768.0;

//c     the multiple angle functions

	    // qo  = 0.0d0
	    qo  = 0.0;
	    // if( w.gt.tol0 )then
	    if(w > tol0){
	    	// qo = -2.0d0*su1*su2/w
	    	qo = -2.0*su1*su2/w;
	    }// endif

	    // q2  = csig+qo
	    q2  = csig + qo;
	    // q4  = 2.0d0*q2*q2-1.0d0
	    q4  = 2.0*q2*q2 - 1.0;
	    // q6  = q2*(4.0d0*q2*q2-3.0d0)
	    q6 = q2*(4.0*q2*q2 - 3.0);
	    // r2  = 2.0d0*ssig*csig
	    r2 = 2.0*ssig*csig;
	    // r3  = ssig*(3.0d0-4.0d0*ssig*ssig)
	    r3 = ssig*(3.0 - 4.0*ssig*ssig);

//c     the longitude difference

	    // s   = sinalf*(ao*sig+a2*ssig*q2+a4*r2*q4+a6*r3*q6)
	    *s = sinalf*(ao*sig + a2*ssig*q2 + a4*r2*q4 + a6*r3*q6);

	    // xz  = dlon+s
	    xz  = dlon + *s;

	    // xy  = dabs(xz-ab)
	    xy  = fabs(xz-ab);
	    // ab  = dlon+s
	    ab  = dlon + *s;

	    // if( xy.lt.0.5d-13 )then
	    	// goto 4
	    // endif

	    // if( kount.le.7 )then
	    	// goto 1
	    // endif
	}

//c     the coefficients of type b

	// 4 z   = epsq*w
	double z = epsq*w;

	// bo  = 1.0d0+z*(1.0d0/4.0d0+z*(-3.0d0/64.0d0+z*(5.0d0/256.0d0-
	// 1         z*175.0d0/16384.0d0)))
	double bo  = 1.0 + z*(1.0/4.0 + z*(-3.0/64.0 + z*(5.0/256.0 -
		z*175.0/16384.0)));

	// b2  = z*(-1.0d0/4.0d0+z*(1.0d0/16.0d0+z*(-15.0d0/512.0d0+
	// 1         z*35.0d0/2048.0d0)))
	double b2 = z*(-1.0/4.0 + z*(1.0/16.0 + z*(-15.0/512.0 +
		z*35.0/2048.0)));

	// b4  = z*z*(-1.0d0/128.0d0+z*(3.0d0/512.0d0-z*35.0d0/8192.0d0))
	double b4 = z*z*(-1.0/128.0 + z*(3.0/512.0 - z*35.0/8192.0));
	// b6  = z*z*z*(-1.0d0/1536.0d0+z*5.0d0/6144.0d0)
	double b6 = z*z*z*(-1.0/1536.0 + z*5.0/6144.0);

//c     the distance in meters

	// s   = b*(bo*sig+b2*ssig*q2+b4*r2*q4+b6*r3*q6)
	*s = b*(bo*sig + b2*ssig*q2 + b4*r2*q4 + b6*r3*q6);

//c     first compute the az1 & az2 for along the equator

	// if( dlon.gt.pi )then
	if(dlon > pi){
		// dlon = (dlon-2.0d0*pi)
		dlon = (dlon - 2.0*pi);
	}// endif

	// if( dabs(dlon).gt.pi )then
	if(fabs(dlon) > pi){
		// dlon = (dlon+2.0d0*pi)
		dlon = (dlon + 2.0*pi);
	}// endif

	// az1 = pi/2.0d0
	*az1 = pi/2.0;
	// if( dlon.lt.0.0d0 )then
	if(dlon < 0.0){
		// az1 = 3.0d0*az1
        *az1 = 3.0*(*az1);
	}// endif

	// az2 = az1+pi
	*az2 = *az1 + pi;
	// if( az2.gt.2.0d0*pi )then
	if(*az2 > 2.0*pi ){
		// az2 = az2-2.0d0*pi
		*az2 = *az2 - 2.0*pi;
	}// endif

//c     now compute the az1 & az2 for latitudes not on the equator

	// if( .not.(dabs(su1).lt.tol0 .and. dabs(su2).lt.tol0) )then
	if( !((fabs(su1) < tol0) && (fabs(su2) < tol0)) ){
		// tana1 =  slon*cu2/(su2*cu1-clon*su1*cu2)
		double tana1 =  slon*cu2/(su2*cu1 - clon*su1*cu2);
		// tana2 =  slon*cu1/(su1*cu2-clon*su2*cu1)
		double tana2 =  slon*cu1/(su1*cu2 - clon*su2*cu1);
		// sina1 =  sinalf/cu1
		double sina1 =  sinalf/cu1;
		// sina2 = -sinalf/cu2
		double sina2 = -sinalf/cu2;

//c       azimuths from north,longitudes positive east

		// az1   = datan2(sina1,sina1/tana1)
		*az1 = atan2(sina1,sina1/tana1);
		// az2   = pi-datan2(sina2,sina2/tana2)
		*az2 = pi - atan2(sina2,sina2/tana2);
	}// endif

	// if( az1.lt.0.0d0 )then
	if(*az1 < 0.0){
		// az1 = az1+2.0d0*pi
		*az1 = *az1 + 2.0*pi;
	}// endif

	// if( az2.lt.0.0d0 )then
	if(*az2 < 0.0){
		// az2 = az2+2.0d0*pi
		*az2 = *az2 + 2.0*pi;
	}// endif

	// return
}

double d_sign(x,y)
double x;
double y;
{
    if (x > 0.0) {
	if (y > 0.0) {
	    return (x);
	} else {
	    return (-x);
	}
    } else {
	if (y < 0.0) {
	    return (x);
	} else {
	    return (-x);
	}
    }
}



int invers(double a, double rf, double b1,
	double l1, double b2, double l2, double *faz,
	double *baz, double *s)
{
    /* System generated locals */
    double d__1, d__2;

int kind, it;
double lam, sig;

double c_b3 = 1e-15;

    /* Local variables */
    double c__, d__, l, z__, ep2, boa, tem1, tem2, biga, bigb,
	    bige, bigf, dsig, temp, prev, test, beta1, beta2, cosu1, cosu2,
	    sinu1, sinu2, sinal, costm, cosal2, costm2, coslam, sinlam,
	    cossig, sinsig;

/* ** inverse for long-line and antipodal cases. */
/* ** latitudes may be 90 degrees exactly. */
/* ** latitude positive north, longitude positive east, radians. */
/* ** azimuth clockwise from north, radians. */
/* ** original programmed by thaddeus vincenty, 1975, 1976 */
/* ** removed back side solution option, debugged, revised -- 2011may01 -- dgm */
/* ** this version of code is interim -- antipodal boundary needs work */
/* ** output (besides faz,baz, and s): */
/* ** it,   iteration count */
/* ** sig,  spherical distance on auxiliary sphere */
/* ** lam,  longitude difference on auxiliary sphere */
/* ** kind, solution flag:  kind=1, long-line;  kind=2, antipodal */
    boa = 1. - 1. / rf;
/* **** sinu1 = boa*dsin(b1)/dsqrt((boa*dsin(b1))**2+dcos(b1)**2) */
/* **** cosu1 = dsqrt(-sinu1**2+1.d0)                               !*** roundoff */
/* **** sinu2 = boa*dsin(b2)/dsqrt((boa*dsin(b2))**2+dcos(b2)**2) */
/* **** cosu2 = dsqrt(-sinu2**2+1.d0)                               !*** roundoff */
    beta1 = atan(boa * tan(b1));
/* *** better reduced latitu */
    sinu1 = sin(beta1);
    cosu1 = cos(beta1);
    beta2 = atan(boa * tan(b2));
/* *** better reduced latitu */
    sinu2 = sin(beta2);
    cosu2 = cos(beta2);
    l = l2 - l1;
/* *** longitude difference */
    if (l > 3.1415926535897932384626433832795) {
	l = l - 3.1415926535897932384626433832795 -
		3.1415926535897932384626433832795;
    }
    if (l < -3.1415926535897932384626433832795) {
	l = l + 3.1415926535897932384626433832795 +
		3.1415926535897932384626433832795;
    }
    prev = l;
    test = l;
    it = 0;
    kind = 1;
    lam = l;
/* ** top of the long-line loop (kind=1) */
/* *** v13   (rap */
L2:
    sinlam = sin(lam);
/* **** if(dabs(pi-dabs(l)).lt.0.2d-10) sinlam=0.d0        !*** no--troublesome */
    coslam = cos(lam);
    temp = cosu1 * sinu2 - sinu1 * cosu2 * coslam;
/* Computing 2nd power */
    d__1 = cosu2 * sinlam;
/* Computing 2nd power */
    d__2 = temp;
    sinsig = sqrt(d__1 * d__1 + d__2 * d__2);
/* *** v14   (rap */
    cossig = sinu1 * sinu2 + cosu1 * cosu2 * coslam;
/* *** v15   (rap */
    sig = atan2(sinsig, cossig);
/* **** sinal = cosu1*cosu2*sinlam/sinsig                  !*** v17   (rapp part II) */
    if (fabs(sinsig) < 1e-15) {
	sinal = cosu1 * cosu2 * sinlam / d_sign(&c_b3, &sinsig);
/* *** avoid div */
    } else {
	sinal = cosu1 * cosu2 * sinlam / sinsig;
    }
/* Computing 2nd power */
    d__1 = sinal;
    cosal2 = -(d__1 * d__1) + 1.;
/* **** costm = -2.d0*sinu1*sinu2/(cosal2+tol)+cossig      !*** v18   (rapp part II) */
    if (fabs(cosal2) < 1e-15) {
	costm = sinu1 * sinu2 / d_sign(&c_b3, &cosal2) * -2. + cossig;
/* *** avoi */
    } else {
	costm = sinu1 * sinu2 / cosal2 * -2. + cossig;
    }
    costm2 = costm * costm;
    c__ = ((cosal2 * -3. + 4.) / rf + 4.) * cosal2 / rf / 16.;
/* ** entry point of the antipodal loop (kind=2) */
/* *** v10   (rap */
L6:
if(it > 1e3) goto L100;
    ++(it);
    d__ = (((costm2 * 2. - 1.) * cossig * c__ + costm) * sinsig * c__ + sig)
	    * (1. - c__) / rf;
    if (kind == 1) {
	lam = l + d__ * sinal;
	if ((d__1 = lam - test, fabs(d__1)) < 1e-14) {
	    goto L100;
	}
	if (fabs(lam) > 3.1415926535897932384626433832795) {
	    kind = 2;
	    lam = 3.1415926535897932384626433832795;
	    if (l < 0.) {
		lam = -(lam);
	    }
	    sinal = 0.;
	    cosal2 = 1.;
	    test = 2.;
	    prev = test;
	    sig = 3.1415926535897932384626433832795 - (d__1 = atan(sinu1 /
		    cosu1) + atan(sinu2 / cosu2), fabs(d__1));
	    sinsig = sin(sig);
	    cossig = cos(sig);
	    c__ = ((cosal2 * -3. + 4.) / rf + 4.) * cosal2 / rf / 16.;
	    if ((d__1 = sinal - prev, fabs(d__1)) < 1e-14) {
		goto L100;
	    }
/* ****     costm = -2.d0*sinu1*sinu2/(cosal2+tol)+cossig */
	    if (fabs(cosal2) < 1e-15) {
		costm = sinu1 * sinu2 / d_sign(&c_b3, &cosal2) * -2. + cossig;
/* *** */
	    } else {
		costm = sinu1 * sinu2 / cosal2 * -2. + cossig;
	    }
	    costm2 = costm * costm;
	    goto L6;
	}
	if ((lam - test) * (test - prev) < 0. && it > 5) {
	    lam = (lam * 2. + test * 3. + prev) / 6.;
	}
/* *** refined */
	prev = test;
	test = lam;
	goto L2;
    } else {
	sinal = (lam - l) / d__;
	if ((sinal - test) * (test - prev) < 0. && it > 5) {
	    sinal = (sinal * 2. + test * 3. + prev) / 6.;
	}
/* *** refined */
	prev = test;
	test = sinal;
/* Computing 2nd power */
	d__1 = sinal;
	cosal2 = -(d__1 * d__1) + 1.;
	sinlam = sinal * sinsig / (cosu1 * cosu2);
/* Computing 2nd power */
	d__2 = sinlam;
	coslam = -sqrt((d__1 = -(d__2 * d__2) + 1., fabs(d__1)));
	lam = atan2(sinlam, coslam);
	temp = cosu1 * sinu2 - sinu1 * cosu2 * coslam;
/* Computing 2nd power */
	d__1 = cosu2 * sinlam;
/* Computing 2nd power */
	d__2 = temp;
	sinsig = sqrt(d__1 * d__1 + d__2 * d__2);
	cossig = sinu1 * sinu2 + cosu1 * cosu2 * coslam;
	sig = atan2(sinsig, cossig);
	c__ = ((cosal2 * -3. + 4.) / rf + 4.) * cosal2 / rf / 16.;
	if ((d__1 = sinal - prev, fabs(d__1)) < 1e-14) {
	    goto L100;
	}
/* ****   costm = -2.d0*sinu1*sinu2/(cosal2+tol)+cossig */
	if (fabs(cosal2) < 1e-15) {
	    costm = sinu1 * sinu2 / d_sign(&c_b3, &cosal2) * -2. + cossig;
/* *** av */
	} else {
	    costm = sinu1 * sinu2 / cosal2 * -2. + cossig;
	}
	costm2 = costm * costm;
	goto L6;
    }
/* ** convergence */
L100:
    if (kind == 2) {
/* *** antipodal */
	*faz = sinal / cosu1;
/* Computing 2nd power */
	d__1 = *faz;
	*baz = sqrt(-(d__1 * d__1) + 1.);
	if (temp < 0.) {
	    *baz = -(*baz);
	}
	*faz = atan2(*faz, *baz);
	tem1 = -sinal;
	tem2 = sinu1 * sinsig - cosu1 * cossig * *baz;
	*baz = atan2(tem1, tem2);
    } else {
/* *** long-line */
	tem1 = cosu2 * sinlam;
	tem2 = cosu1 * sinu2 - sinu1 * cosu2 * coslam;
	*faz = atan2(tem1, tem2);
	tem1 = -cosu1 * sinlam;
	tem2 = sinu1 * cosu2 - cosu1 * sinu2 * coslam;
	*baz = atan2(tem1, tem2);
    }
    if (*faz < 0.) {
	*faz = *faz + 3.1415926535897932384626433832795 +
		3.1415926535897932384626433832795;
    }
    if (*baz < 0.) {
	*baz = *baz + 3.1415926535897932384626433832795 +
		3.1415926535897932384626433832795;
    }
/* ** Helmert 1880 from Vincenty "Geodetic inverse solution between antipodal points" */
    ep2 = 1. / (boa * boa) - 1.;
    bige = sqrt(ep2 * cosal2 + 1.);
    bigf = (bige - 1.) / (bige + 1.);
    biga = (bigf * bigf / 4. + 1.) / (1. - bigf);
    bigb = bigf * (1. - bigf * .375 * bigf);
/* Computing 2nd power */
    d__1 = sinsig;
    z__ = bigb / 6. * costm * (d__1 * d__1 * 4. - 3.) * (costm2 * 4. - 3.);
    dsig = bigb * sinsig * (costm + bigb / 4. * (cossig * (costm2 * 2. - 1.)
	    - z__));
    *s = boa * a * biga * (sig - dsig);
    return 0;
} /* invers_ */

ErrorSet inverse(LLPoint origin, LLPoint dest, double* crs, double *bcrs,
                  double* dist, double eps){

	ErrorSet err = 0;
	//    COMMON/CONST/PI,RAD
		double PI = 4.0*atan(1.0); //	pi=4.d0*datan(1.d0)

	//    COMMON/ELIPSOID/A,F
		double A = 6378137.0; //	a=6378137.d0
	//	double F = 1.0/298.25722210088; //	f=1.d0/298.25722210088d0 this is the GRS80 flattening
		double F = FLATTENING; //	WGS84 flattening
		double esq = F*(2.0-F);


	double npCrs = 0.0, npBcrs = 0.0, npDist = 0.0;
    /* Assign local storage if optional pointers are not provided */
    if (NULL == crs)
    {
        crs = &npCrs;
    }
    if (NULL == bcrs)
    {
        bcrs = &npBcrs;
    }
    if (NULL == dist)
    {
        dist = &npDist;
    }

//	gpnhri(A, F, esq, PI, origin.latitude, origin.longitude, dest.latitude, dest.longitude, crs, bcrs, dist);
	invers(A, 1/F, origin.latitude, origin.longitude, dest.latitude, dest.longitude, crs, bcrs, dist);
	*dist = *dist/1852.0;

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

/* Table of constant values */



