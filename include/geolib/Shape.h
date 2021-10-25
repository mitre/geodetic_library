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

/*! \file Shape.h
 *  \brief This file defines the shape structures used by geolib.
 */

/* A NOTE ABOUT PARAMETERS AND UNITS
 * Unless otherwise noted below, all distances are in nautical miles and
 * all angles, latitudes, longitudes, courses, and azimuths are in radians.
 * Standard conversions are:
 *   nautical miles to feet: multiply by 1852/0.3048 or use FEET_PER_NMI macro defined below
 *   radians to degrees: multiply by 180.0/pi, or use DEG_PER_RAD macro defined below.
 *
 * \author Michael Mills, Richard Snow, Stuart Bowman, Juan Amezcua, John Landrigan
 *
 */

#pragma once

namespace geolib_idealab {
#ifdef TARGETS
#define SPECIAL_DOUBLE double
#else
#define SPECIAL_DOUBLE long double
#endif

#undef INFINITE

#ifndef LINETYPE
/** @enum LineType
 *
 * Applies to and describes the exent of a Geodesic or Locus struct 
 */
typedef enum {
  /* Values given explicitly for clarity */
  SEGMENT = 0, /**< The geodesic or locus exists only between its start and end points */
  SEMIINFINITE = 1, /**< The geodesic or locus extends beyond the end point */
  INFINITE =
  2 /**< The geodesic or locus extends before the start point and beyond the end point (unbounded would be a better term) */
} LineType;
#define LINETYPE
#endif

#ifndef ARCDIRECTION
/** @enum ArcDirection
 *
 * Defines the orientation of and Arc struct from the perspective of viewing the Arc from space.
 */
typedef enum {
  CLOCKWISE = 1, /**< Clockwise means the Arc turns to the right from the startPoint */
  COUNTERCLOCKWISE = -1  /**< Counterclockwise means the Arc turns to the left from the startPoint */
} ArcDirection;
#define ARCDIRECTION
#endif

#ifndef LLPOINT_STRUCT
/** @struct LLPoint

 * The LLPoint struct represents a geodetic position defined by a latitude
 * and longitude pair. The latitude and longitude are stored in radians.
 *
 **/
typedef struct {
  /* Coordinates must be in radians. */
  double latitude;   /**< Number representing the latitude. Valid range is \f$[-\pi/2,\pi/2]\f$.
			  A positive latitude indicates a position in the northern hemisphere*/
  double longitude;  /**< Number representing the longitude.
			* Valid range is \f$[-\pi,\pi]\f$.  A positive longitude indicates a position in the eastern hemisphere */
} LLPoint;
#define LLPOINT_STRUCT
#endif

#ifndef LLPOINTPAIR
/** @typedef LLPointPair
 *
 * LLPointPair is an alias for a fixed-length two element array of LLPoints
 *
 */
typedef LLPoint LLPointPair[2];
#define LLPOINTPAIR
#endif

#ifndef LLPOINTSET
/** @struct LLPointSet
 * 
 * An array of pointers to LLPoint structs
 */
typedef struct {
  int length;      /**< Number of elements currently in array */
  LLPoint **elements; /**< Array of LLPoint pointers                    */
  int maxLength;   /**< max number of pointers that can be held */
}
    LLPointSet;
#define LLPOINTSET
#endif

#define LLPOINTSET_ARRAY_INCREMENT 20

/** Initialize and return an LLPointSet object.
 * @return Returns an LLPointSet with length = 0, elements pointer initialized and
 * pointing to unitialized memory sufficiently large to hold LLPOINTSET_ARRAY_INCREMENT Shape pointers.
 */
LLPointSet createPtSet();

/** Deallocate memory reserved for LLPointSet struct. Memory pointed to by given
 * pointer is freed.
 * @param set LLPointSet object to be cleared. (LLPointSet*)
 * @return Nothing
 */
void clearPtSet(LLPointSet *set);

/** Adds a point to an LLPointSet object
 * @param set LLPointSet to add a point to (LLPointSet*)
 * @param point LLPoint to be added to the LLPointSet (LLPoint*)
 * @return Updates the LLPointSet with the additional point.
 */

ErrorSet addPtToPtSet(LLPointSet *set, LLPoint *point);

#ifndef VECTOR_STRUCT
/** @struct Vector
 * 
 * A vector in Cartesian 3-space.
 * @param x Double 
 * @param y Double
 * @param z Double
 *
 */
typedef struct {
  SPECIAL_DOUBLE x;
  SPECIAL_DOUBLE y;
  SPECIAL_DOUBLE z;
} Vector;
#define VECTOR_STRUCT
#endif

/** @struct Locus
 * The Locus struct represents a geodetic curve that is defined by a central
 * geodesic and distance function.  Given any point \f$p\f$ on the defining geodesic
 * and a distance function \f$d = f(p),\f$ then the corresponding point on the locus
 * is found by projecting a point perpendicular to the geodesic at \f$p\f$ and moving a
 * distance \f$f(p)\f$.  Negative distances indicate that the locus point is
 * to the left of the geodesic (from the point of view of someone facing in the
 * direction of the geodesic; positive distances indicate that the locus is to
 * the right. (Locus)
 *
 * FUTURE ENHANCEMENT: Replace the geo* parameters with a Geodesic struct.
 */
#ifndef LOCUS_STRUCT
typedef struct {
  LLPoint geoStart; /**< Start point of defining geodesic                   */
  LLPoint geoEnd;   /**< End point of defining geodesic                     */
  double geoLength; /**< Length of defining geodesic [NMI] */
  double geoAz;     /**< Azimuth from geoStart to geoEnd [radians] */
  double geoRevAz;  /**< Azimuth from geoEnd to geoStart [radians] */
  LLPoint locusStart; /**< start point of locus (computed at initialization) */
  LLPoint locusEnd; /**< end point of locus (computed at initialization) */
  double startDist; /**< distance to locus at geoStart [nmi]       */
  double endDist;   /**< distance to locus at geoEnd [nmi]         */
  LineType lineType; /* enum that defines the extent of the locus */
  double slope;      /* slope of the locus                */
} Locus;
#define LOCUS_STRUCT
#endif

/** @struct Locus2ndOrder
 * The Locus2ndOrder struct represents a geodetic curve that is defined by a central
 * locus and distance function.  Given any point \f$p\f$ on the defining locus
 * and a distance function \f$d = f(p),\f$ then the corresponding point on the 2nd order locus
 * is found by projecting a point perpendicular to the locus at \f$p\f$ and moving a
 * distance \f$f(p)\f$.  Negative distances indicate that the 2nd order locus point is
 * to the left of the locus (from the point of view of someone facing in the
 * direction of the locus; positive distances indicate that the 2nd order locus is to
 * the right.
 */
#ifndef LOCUS2NDORDER_STRUCT
typedef struct {
  Locus locus; /**< Defining locus */
  double secondDist; /**< distance to 2nd order locus */
} Locus2ndOrder;
#define LOCUS2NDORDER_STRUCT
#endif

#ifndef SPIRAL_STRUCT
/** @struct Spiral
 * The Spiral struct describes an archimedean (or arithmetic) spiral on the surface of the ellipsoid.  The spiral is defined by
 * a center LLPoint, a start LLPoint, an end LLPoint, and an orientation.  The orientation convention is identical to the one for
 * Arcs.  The radius for any point between the start and end points is determined by a linear function of azimuth.  The slope of this
 * function is the spiral's "growth rate".
 *
 */
typedef struct {
  LLPoint centerPoint; /**< The center of the Spiral */
  LLPoint startPoint; /**< The start of the Spiral */
  LLPoint endPoint; /**< The end of the Spiral */
  double startRadius; /**< The distance from centerPoint to startPoint (computed at initialization)   */
  double startAz; /**< azimuth from centerPoint to startPoint [radians] (computed at initialization)   */
  double endRadius; /**< The distance from centerPoint to endPoint (computed at initialization)   */
  double endAz; /**< azimuth from center to endPoint [radians] (computed at initialization)   */
  double growthRate; /**< rate of change of radius with azimuth angle [NMI/radian] */
  double subtendedAngle; /**< Difference betweeen endAz and startAz */
  ArcDirection dir; /* -1 => counterclockwise, +1 => clockwise     */
} Spiral;
#define SPIRAL_STRUCT
#endif

#ifndef GEODESIC_STRUCT
/** @struct Geodesic
 *
 * A Geodesic represents the shortest path between two LLPoints on the surface of the ellipsoid. It is defined by its start and end LLPoints.
 *
 */
typedef struct {
  LLPoint startPoint; /**< Start of geodesic */
  LLPoint endPoint; /**< End of geodesic */
  LineType lineType; /**< enumerator indicating the extend of the geodesic */
  double
      startAz;  /**< Azimuth of geodesic at start point facing toward end point (computed at time of initialization) */
  double endAz; /**< Azimuth of geodesic at end point, facing away from start point (computed at time of definition) */
  double length; /** Distance from startPoint to endPoint in nautical miles (computed at time of definition) */
} Geodesic;
#define GEODESIC_STRUCT
#endif

/** @struct Arc
 * Defines all or part of a small circle of the surface of the ellipsoid. Consists of all points 
 * a distance equal to the radius from the centerPoint, and between the startPoint and endPoint. What is
 * meant by "between" depends upon the value of the ArcDirection parameter. 
 * 
 */
#ifndef ARC_STRUCT
typedef struct {
  LLPoint centerPoint; /**< The center of the arc */
  LLPoint startPoint;  /**< The start of the arc */
  double startAz;      /**< Azimuth at centerPoint of geodesic from centerPoint to startPoint */
  LLPoint endPoint;    /**< The end of the arc */
  double endAz;        /**< Azimuth at centerPoint of geodesic from centerPoint to endPoint */
  double radius;       /**< The distance from the centerPoint to any point on the arc */
  ArcDirection dir;    /**< -1 => counterclockwise, +1 => clockwise */
  double
      subtendedAngle; /**< Normalized difference betweeen startAz and endAz. If dir = CLOCKWISE, then subtendedAngle > 0 */
} Arc;
#define ARC_STRUCT
#endif

/** @enum ShapeType
 *
 * Label which indicates the type of element in a Boundary's shape array.
 */
#ifndef SHAPETYPE
typedef enum {
  ARC = 0,
  GEODESIC = 1,
  LOCUS = 2,
  SPIRAL = 3,
  LLPOINT = 4,
  BOUNDARY = 5,
  LOCUS2NDORDER = 6
}
    ShapeType;
#define SHAPETYPE
#endif

#ifndef SHAPE
/** @struct Shape
 *
 * Facilitates storage of mutable shapes in Boundary structs. Consists of a void pointer and label.
 *
 */
typedef struct {
  void *this_shape;  /**< Pointer to shape struct */
  ShapeType type; /**< Enumerator that indicates the type of shape */
}
    Shape;
#define SHAPE
#endif

#ifndef BOUNDARY
/** @struct Boundary
 *
 * A Boundary is a collection of shapes. Boundaries must be closed (every startPoint coincides with an endPoint) and simply connected.
 *
 */
typedef struct {
  int length;      /**< Number of elements currently in array */
  Shape *elements; /**< Array of shape                    */
  int
      maxLength;   /**< Current max number of Shape pointers that can be held. Incremented if maxLength+1 elements are added */
}
    Boundary;
#define BOUNDARY
#endif

#ifndef COMPLEX_BOUNDARY
typedef struct {
  int length;      /* Number of boundaries currently in array */
  Boundary **elements; /* Array of boundary pointers                 */
  int maxLength;   /* max number of Boundary pointers that can be held */
}
    ComplexBoundary;
#define COMPLEX_BOUNDARY
#endif

/** Initialize an Arc structure from defining parameters
 * @param arc Pointer to Arc structure that is to be initialized. (Arc*)
 * @param center Center point of geodesic (LLPoint)
 * @param startPoint Start point of arc (LLPoint)
 * @param endPoint End point of arc (LLPoint)
 * @param dir ArcDirection parameter (ArcDirection)
 * @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
 * @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
 * @return Updates the input Arc struct with the constructed Arc
 * @retval SUCCESS Indicates successful execution. */
ErrorSet createArc(Arc *arc, LLPoint center, LLPoint startPoint,
                   LLPoint endPoint, ArcDirection dir, double tol, double eps);

/** Initialize a Geodesic structure from defining parameters
 * @param geo Pointer to Geodesic structure that is to be initialized. (Geodesic*)
 * @param geoStart Start point of geodesic (LLPoint)
 * @param geoEnd End point of geodesic (LLPoint)
 * @param lineType LineType enum indicating extent of locus (LineType)
 *                   SEGMENT => locus exists only between start and end
 *                   SEMIINFINITE => locus extends beyond end
 *                   INFINITE => locus is unbounded in both directions
 * @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
 * @return Updates the input Geodesic struct with the constructed Geodesic
 * @retval SUCCESS Indicates successful execution. */
ErrorSet createGeo(Geodesic *geo, LLPoint geoStart, LLPoint geoEnd,
                   LineType lineType, double eps);

/** Initialize a Locus structure from defining parameters
 * @param loc Pointer to locus structure that will be initialized (Locus)
 * @param geoStart Start point of defining geodesic (LLPoint)
 * @param geoEnd End point of defining geodesic (LLPoint)
 * @param startDist Distance from geoStart to locus in nmi
 *                    startDist > 0 if locus is to right of geoStart
 *                    startDist < 0 if locus is to left of geoStart (double)
 * @param endDist Distance from geodesic to locus at endPoint in nmi
 *                  endDist > 0 if locus is to right of geoEnd
 *                  endDist < 0 if locus is to left of geoEnd (double)
 * @param lineType LineType enum indicating extent of locus
 *                   SEGMENT => locus exists only between start and end
 *                   SEMIINFINITE => locus extends beyond end
 *                   INFINITE => locus is unbounded in both directions (LineType)
 * @param tol Accuracy tolerance in nmi (max distance from found solution to true solution) (double)
 * @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
 * @return Updates the input Locus struct with the constructed Locus
 * @retval SUCCESS Indicates successful execution. */
ErrorSet createLocus(Locus *loc, LLPoint geoStart, LLPoint geoEnd,
                     double startDist, double endDist, LineType lineType,
                     double tol, double eps);

/** Initialize a new Spiral structure.
 * @param sp Pointer to new Spiral (Spiral)
 * @param center Center point of the spiral (LLPoint)
 * @param startRad Initial radius of the spiral (double)
 * @param endRad Final radius of the spiral (double)
 * @param startAz Azimuth from the center point to the start point of the spiral (double)
 * @param endAz Azimuth from the center point to the end point of the spiral (double)
 * @param dir Turn direction of the spiral (ArcDirection)
 * @param eps Convergence tolerance for Vincenty forward/inverse algorithms (double)
 * @return Updates the input Spiral struct with the constructed Spiral
 * @retval SUCCESS Indicates successful execution. */

ErrorSet createSpiral(Spiral *sp,
                      LLPoint center,
                      double startRad,
                      double endRad,
                      double startAz,
                      double endAz,
                      ArcDirection dir,
                      double eps);

/** Initialize a LLPoint structure from defining parameters
 * @param llpoint Pointer to LLPoint structure that will be initialized (LLPoint)
 * @param lat Latitude (in radians) (double)
 * @param lon Longitude (in radians) (double)
 * @return Updates the input LLPoint struct with the constructed LLPoint
 * @retval SUCCESS Indicates successful execution.*/
ErrorSet createPt(LLPoint *llpoint, double lat, double lon);

} // namespace