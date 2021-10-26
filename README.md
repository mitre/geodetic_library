# MITRE Geodetic Library

`Geodetic library` (or `geolib`) is a library for performing WGS-84 calculations with high precision. We think it's very handy and use it regularly in our internal aviation software. Enjoy!

## Algorithm Documentation

This work was originally performed by MITRE for the FAA. Technical documentation can be found in [FAA Order 8260.58A](https://www.faa.gov/regulations_policies/orders_notices/index.cfm/go/document.information/documentid/1029267), Appendix E. It describes important algorithms used for aviation purposes. This library is our software implementation of the algorithms discussed there.

### Calling Algorithms In This Library

There are two basic mathematical operations on the ellipsoid: `direct` and `inverse`. (Read more about the [mathematical background](https://en.wikipedia.org/wiki/Vincenty%27s_formulae) on Wikipedia.) We solve these here with high precision using the `Vincenty` formulation (default implementation, see compile options for other mathematical implementations). Once you've linked the library, the algorithms are available in the `geolib_idealab` namespace. You can call these two algorithms like this:

```c
using namespace geolib_idealab;

LLPoint origin;
origin.latitude = 0.0;
origin.longitude = 0.0;
std::string pt_name = "origin";
int n = pt_name.length();
char char_array[n + 1];
strcpy(char_array, pt_name.c_str());

LLPoint destination;
destination.latitude = 1.0*M_PI/180.0;
destination.longitude = 0.0;
pt_name = "destination";
n = pt_name.length();
char char_array2[n + 1];
strcpy(char_array2, pt_name.c_str());

// direct operation
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
        std::cout << "Something went wrong with the direct() operation! " << formatErrorMessage(error_set_direct) << std::endl;
}

// inverse operation
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
        std::cout << "Something went wrong with the inverse operation! " << formatErrorMessage(error_set_inverse) << std::endl;
}
```

### Example Binary

See also the `examples` binary built by this repository in /src/example. It uses [CPM](https://github.com/cpm-cmake/CPM.cmake) to retrieve, build, and link against `geolib`. Then it makes some sample calls to show `direct()` and `inverse()` operations on the ellipsoid.

To build, perform this `cmake` operation:

```bash
mkdir build
cd build
cmake ../src/example
make
```

This will produce a binary in /src/example/bin called `examples` that links against `libgeolib.a` and calls it in the `main` routine. 

## Build

Use `cmake` in the normal way.

```bash
mkdir build
cd build
cmake ..
make
```

You'll get two artifacts:

* the geodetic library `libgeolib.a` in /lib. 
* a test binary `testGeolib` in /bin.

This library & sample binary builds and runs on **CentOs7** and **MacOsX**.

**NOTE 1**: Several c-code source files have been intentionally ignored in the build of the library. The goal in removing unneeded code from the library was to keep callers from using code that has not been recently tested. The ignored source code may be added back into the library build as needed, but it should be tested as it is added back in. See [CMakeLists.txt](/geolib/src/main/c/CMakeLists.txt) for more information.

**NOTE 2**: The `geolib` library can be compiled to run one of three solvers: Vincenty, Karney, or NGSVincenty. The solver decision is required to occur at compile-time. Only one solver may be compiled into the library. See [CMakeLists.txt](/CMakeLists.txt) for more details.

## Automatic Builds

We've been using CI internally for this library for years. We'll try to get a public CI working soon.

## Issues

We want to hear from you. Please use standard GitHub Issue and Discussion tools for interaction with the developers.

## Tests

Yup, we have 'em. See /test.

We currently use a custom test harness for our testing. The goal is to switch over to running with `gtest` at some point in the future.

## CPM Snippet

`geolib` is a library that you'll want to link against, not an application to run. A recommended way to do this with your software is through the [C++ Package Manager](https://github.com/cpm-cmake/CPM.cmake) (or CPM). Add this into your CMake profile and you should be good-to-go.

Drop this in the appropriate place for your application:
```cmake
CPMAddPackage(
        NAME geolib_library
        GITHUB_REPOSITORY mitre/geodetic_library
        VERSION 3.2.7-SNAPSHOT
        GIT_TAG main
)
```
When using this library, we strongly recommend using a tagged, versioned build of geolib whenever possible.

Then, where you do your linking commands in cmake, also add this:
```cmake
target_include_directories(my-application PUBLIC ${geolib_library_SOURCE_DIR}/include)
target_link_libraries(my-application geolib)
```

Once these are done, your application (`my-application`) will be linked against `geolib_library` and you can call it from your code. Like this:

```c++
#include "geolib/cppmanifest.h"
int main(int argc, char *argv[]) {
  std::cout << "geolib-library build version: " << geolib_idealab::cppmanifest::getVersion() << std::endl;
}
```

## Open Source Notice

This content has been approved for public release.

This is the copyright work of The MITRE Corporation, and was produced for the U. S. Government under Contract Number DTFAWA-10-C-00080, and is subject to Federal Aviation Administration Acquisition Management System Clause 3.5-13, Rights In Data-General, Alt. III and Alt. IV (Oct. 1996). No other use other than that granted to the U. S. Government, or to those acting on behalf of the U. S. Government, under that Clause is authorized without the express written permission of The MITRE Corporation. For further information, please contact The MITRE Corporation, Contracts Office, 7515 Colshire Drive, McLean, VA 22102-7539, (703) 983-6000.

Copyright 2021 The MITRE Corporation. All Rights Reserved. Approved for Public Release; Distribution Unlimited. 15-1482

This project contains content developed by The MITRE Corporation. If this code is used in a deployment or embedded within another project, it is requested that you send an email to opensource@mitre.org in order to let us know where this software is being used.

## License

This library is made available via the [Apache 2.0 license](https://www.apache.org/licenses/LICENSE-2.0). Please send an email to opensource@mitre.org to get answers to licensing questions.