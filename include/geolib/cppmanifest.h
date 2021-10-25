// ****************************************************************************
// NOTICE
//
// This is the copyright work of The MITRE Corporation, and was produced
// for the U. S. Government under Contract Number DTFAWA-10-C-00080, and
// is subject to Federal Aviation Administration Acquisition Management
// System Clause 3.5-13, Rights In Data-General, Alt. III and Alt. IV
// (Oct. 1996).  No other use other than that granted to the U. S.
// Government, or to those acting on behalf of the U. S. Government,
// under that Clause is authorized without the express written
// permission of The MITRE Corporation. For further information, please
// contact The MITRE Corporation, Contracts Office, 7515 Colshire Drive,
// McLean, VA  22102-7539, (703) 983-6000. 
//
// Copyright 2020 The MITRE Corporation. All Rights Reserved.
// ****************************************************************************

/*
 * This is just an aggregator header file. 
 */

#if !defined (GEOLIB_IDEALAB_CPPMANIFEST_H)
#define GEOLIB_IDEALAB_CPPMANIFEST_H

#if defined (GEOLIB_IDEALAB_CPPMANIFEST_HAVE_PRAGMA_ONCE)
#pragma once
#endif

#include "geolib_build_info.h"

namespace geolib_idealab {
  namespace cppmanifest {

  static std::string getUserName() {
    return GEOLIB_IDEALAB_CPPMANIFEST_USERNAME_STR;
  }

  static std::string getBuildTimeStamp() {
    return GEOLIB_IDEALAB_CPPMANIFEST_BUILDTIMESTAMP_STR;
  }

  static std::string getBuildCompilerVersion() {
    return GEOLIB_IDEALAB_CPPMANIFEST_BUILDGCCVERSION_STR;
  }

  static std::string getBuildSystemVersion() {
    return GEOLIB_IDEALAB_CPPMANIFEST_BUILDSYSTEMVERSION_STR;
  }

  static std::string getBuildSystemName() {
    return GEOLIB_IDEALAB_CPPMANIFEST_BUILDSYSTEMNAME_STR;
  }

  static std::string getBuildSystemProcessor() {
    return GEOLIB_IDEALAB_CPPMANIFEST_BUILDSYSTEMPROCESSOR_STR;
  }

  static std::string getBuildHostName() {
    return GEOLIB_IDEALAB_CPPMANIFEST_BUILDHOSTNAME_STR;
  }

  static std::string getGitBranch() {
    return GEOLIB_IDEALAB_CPPMANIFEST_GIT_BRANCH;
  }

  static std::string getGitHash() {
    return GEOLIB_IDEALAB_CPPMANIFEST_GIT_HASH;
  }

  static bool getGitIsClean() {
    if (GEOLIB_IDEALAB_CPPMANIFEST_GIT_LOCAL_CHANGES == std::string("CLEAN")) {
      return 1;
    } else {
      return 0;
    }
  }

   static std::string getVersion() {
       std::string verStr(GEOLIB_VERSION_STR);
       bool stripLastChar = verStr.find("-") == verStr.length() - 1;
       if (stripLastChar) {
           // dump the last char
           verStr.resize(verStr.find("-"));
       }
       return verStr;
   }
  }
}

#endif
/*
 * Do not put content below here.
 */
