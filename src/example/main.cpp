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

//
// Example code for using the geolib c-library.
//

#include <string>
#include <cstring>
#include <iostream>
#include <cfloat>
#include "geolib/Geolib.h"
#include "geolib/Util.h"
#include "geolib/Shape.h"
#include "geolib/Constants.h"
#include "geolib/ErrorCodes.h"
#include "geolib/cppmanifest.h"

const std::string buildinfoFlag("--buildinfo");

using namespace geolib_idealab;

int main(int argc,
         char *argv[]) {

  // handle command line flag --buildinfo
  if (argc == 2) {
    std::string arg1(argv[1]);
    if (arg1 == buildinfoFlag) {
      std::cout << "\n\n**************" << std::endl;
      std::cout << "geolib-idealab build info:" << std::endl;
      std::cout << "Build version: " << cppmanifest::getVersion() << std::endl;
      std::cout << "Created by: " << cppmanifest::getUserName() << std::endl;
      std::cout << "Created date-time: " << cppmanifest::getBuildTimeStamp() << std::endl;
      std::cout << "Built with GCC version: " << cppmanifest::getBuildCompilerVersion() << std::endl;
      std::cout << "Built on system name: " << cppmanifest::getBuildSystemName() << std::endl;
      std::cout << "Built on system processor: " << cppmanifest::getBuildSystemProcessor() << std::endl;
      std::cout << "Built with system ver: " << cppmanifest::getBuildSystemVersion() << std::endl;
      std::cout << "Built on system host name: " << cppmanifest::getBuildHostName() << std::endl;
      std::cout << "Built from git branch: " << cppmanifest::getGitBranch() << std::endl;
      if (cppmanifest::getGitIsClean()) {
        std::cout << "Built from git hash: " << cppmanifest::getGitHash() << std::endl;
      } else {
        std::cout << "Built from git hash: " << cppmanifest::getGitHash() << "-DIRTY" << std::endl;
      }
      std::cout << "**************\n\n" << std::endl;
    }
  }

  // Write some output
  std::string geolib_idealab_version = cppmanifest::getVersion();
  std::cout << "Example calls using the MITRE geolib library: " << geolib_idealab_version.c_str() << std::endl;

  std::cout << "Start at (0, 0) and project north to (1.0, 0)..." << std::endl;

  // Prepare some variables declared by geolib
  LLPoint origin;
  origin.latitude = 0.0;
  origin.longitude = 0.0;
  std::string pt_name = "origin";
  int n = pt_name.length();
  char char_array[n + 1];
  strcpy(char_array, pt_name.c_str());
  displayPt(origin, char_array, false);

  LLPoint destination;
  destination.latitude = 1.0*M_PI/180.0;
  destination.longitude = 0.0;
  pt_name = "destination";
  n = pt_name.length();
  char char_array2[n + 1];
  strcpy(char_array2, pt_name.c_str());
  displayPt(destination, char_array2, false);

  // Call inverse operation
  double epsilon = DEFAULT_EPS;
  double crs_radians_calculated = DBL_MIN;
  double reverse_course_radians_calculated = DBL_MIN;
  double distance_nm_calculated = DBL_MIN;
  ErrorSet error_set_inverse = inverse(
      origin,
      destination,
      &crs_radians_calculated,
      &reverse_course_radians_calculated,
      &distance_nm_calculated,
      epsilon);
  if (error_set_inverse == SUCCESS) {
    std::cout << "\ndistance_nm_calculated: " << distance_nm_calculated << std::endl;
    std::cout << "course: " << crs_radians_calculated << std::endl;
  } else {
    std::cout << "BLAH. Something went wrong with the inverse operation! " << formatErrorMessage(error_set_inverse) << std::endl;
  }

  // Call direct operation, using the above calculated values to see if we get the
  // originally defined destination.
  std::cout << "\n\n**************" << std::endl;
  std::cout << "Now reverse the operation to see if we can calulate the destination..." << std::endl;
  LLPoint destination_calculated;
  ErrorSet error_set_direct = direct(
      origin,
      crs_radians_calculated,
      distance_nm_calculated,
      &destination_calculated,
      epsilon);
  if (error_set_direct == SUCCESS) {
    std::string pt_name = "destination_calculated";
    int n = pt_name.length();
    char char_array[n + 1];
    strcpy(char_array, pt_name.c_str());
    displayPt(destination_calculated,char_array,false);
  } else {
    std::cout << "BLAH. Something went wrong with the direct() operation! " << formatErrorMessage(error_set_direct) << std::endl;
  }

  // return
  return 0;
}