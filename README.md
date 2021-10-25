# MITRE Geodetic Library

`Geodetic library` is a library for performing WGS-84 calculations with high precision. It's very handy. Enjoy!

## Build

Use `cmake` in the normal way. You'll get two artifacts:

* the geodetic library `libgeolib.a` in project_root/lib. 
* an example binary called `examples` that links against `libgeolib.a` and calls it in the `main`. See [main.cpp](/geolib/src/main/examples/main.cpp) and project_root/bin. 

This library & sample binary builds and runs on **CentOs7** and **MacOsX**.

**NOTE**: Several c-code source files have been intentionally ignored in the build of the `geolib` library. These implementations appeared unneeded from the perspective of the AAESim project. The goal in removing unneeded code from the library was to keep callers from using code that has not been recently tested. The ignored source code may be added back into the library build as needed, but it should be tested as it is added back in. See [CMakeLists.txt](/geolib/src/main/c/CMakeLists.txt) for more information.  

**NOTE**: The `geolib` library can be compiled to run one of three solvers: Vincenty, Karney, or NGSVincenty. For AAESim purposes, Vincenty is selected. See [CMakeLists.txt](/CMakeLists.txt) for more details. 

## Container

You can develop and build this for the idealab. See .docker/Dockerfile and associated README.

## Automatic Builds

![Build Status](https://pandafood.mitre.org/plugins/servlet/wittified/build-status/AAES-GEOLIB)

https://pandafood.mitre.org/browse/AAES-GEOLIB
 
This project does not deploy its library artifact anywhere. But, there is a builder on Bamboo. See the link above. It provides both build artifacts as built from the main branch on **CentOs7** and **MacOsX**. 

## Issues

Please post issues to JIRA at the link below. Use the project AAES and the component "geolib".

https://anthill.mitre.org/

## Tests

Yup, we have 'em. Thanks to the parent project for making a tested project. See the Bamboo job for details. 

## CPM Snippet

`geolib-idealab` is really a library that you'll want to link against. A recommended way to do this with your IDEA Lab software is through the [C++ Package Manager](https://github.com/cpm-cmake/CPM.cmake) (or CPM). Do this and you should good-to-go.

Drop this in the appropriate place for your application:
```cmake
CPMAddPackage(
        NAME geolib_idealab
        GIT_REPOSITORY https://mustache.mitre.org/scm/aaes/geolib.git
        VERSION 3.2.7-SNAPSHOT
        GIT_TAG main
)
if (geolib_idealab_ADDED)
    message("found and retrieved geolib_idealab dependency")
    message("geolib_idealab_SOURCE_DIR: ${geolib_idealab_SOURCE_DIR}")
    message("geolib_idealab_BINARY_DIR: ${geolib_idealab_BINARY_DIR}")

    # set a variable for the includes. Use the variable later inside target_include_directories()
    set(geolib_idealab_INCLUDE_DIRS ${geolib_idealab_SOURCE_DIR}/geolib/src/main/include)
endif ()
```
In the above, please use a tagged, versioned build of geolib whenever possible.

Then, where you do your linking commands in cmake, also add this:
```cmake
target_include_directories(my-application PUBLIC
        ${geolib_idealab_INCLUDE_DIRS}
        )
target_link_libraries(my-application geolib-idealab)
```

Once these are done, your application (`my-application`) will be linked against `geolib-idealab` and you can call it from your code. Like this:

```c++
#include "geolib/cppmanifest.h"
int main(int argc, char *argv[]) {
  std::cout << "geolib-idealab build version: " << geolib_idealab::cppmanifest::getVersion() << std::endl;
}
```

## Examples

See also the `examples` binary built by this repository. 