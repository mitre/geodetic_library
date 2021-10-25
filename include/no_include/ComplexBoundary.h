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

#ifndef COMPLEX_BOUNDARY_H_
#define COMPLEX_BOUNDARY_H_

#endif /*COMPLEX_BOUNDARY_H_*/

#include "Geolib.h"

#define COMPLEX_BOUNDARY_ARRAY_INCREMENT 20



/** Initialize and return a Complex Boundary object.
 * @return A Complex Boundary with length = 0, elements pointer initialized and
 * pointing to unitialized memory sufficiently large to hold COMPLEX_BOUNDARY_ARRAY_INCREMENT Boundary ???.
 */
ComplexBoundary newComplexBoundary();

/** Deallocate memory reserved for ComplexBoundary struct. Memory pointed to by given
 * pointer is freed.
 * @return Nothing
 */
void clearComplexBoundary(ComplexBoundary* b);

ErrorSet addElementToComplexBoundary(ComplexBoundary* c, Boundary* element);

ErrorSet boundaryCircleListIntersections(Boundary boundary, Arc circles[], int circlesSize,
                                    ComplexBoundary* commonAreas, double tol, double eps);

ErrorSet longitudinallyPartitionBoundary(Boundary orderedBoundary, ComplexBoundary* complexBoundary, double tol, double eps);

ErrorSet complexBoundaryCircleIntersectionExists(Boundary boundary, Arc circles[], int circlesSize, int checkSurface,
                                    Boundary* intersectionList, Boundary* noIntersectionList, double tol, double eps);

