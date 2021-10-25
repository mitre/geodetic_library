# MITRE Geodetic Library

`Geodetic library` (or `geolib`) is a library for performing WGS-84 calculations with high precision. We think it's very handy and use it regularly in our internal software. Enjoy!

## Build

Use `cmake` in the normal way. You'll get two artifacts:

* the geodetic library `libgeolib.a` in /lib. 
* an example binary called `examples` that links against `libgeolib.a` and calls it in the `main`. See [main.cpp](/src/example/main.cpp) and /bin. 

This library & sample binary builds and runs on **CentOs7** and **MacOsX**.

**NOTE 1**: Several c-code source files have been intentionally ignored in the build of the library. The goal in removing unneeded code from the library was to keep callers from using code that has not been recently tested. The ignored source code may be added back into the library build as needed, but it should be tested as it is added back in. See [CMakeLists.txt](/geolib/src/main/c/CMakeLists.txt) for more information.

**NOTE 2**: The `geolib` library can be compiled to run one of three solvers: Vincenty, Karney, or NGSVincenty. The solver decision is required to occur at compile-time. Only one solver may be compiled into the library. See [CMakeLists.txt](/CMakeLists.txt) for more details.

## Automatic Builds

We've been using CI internally for this library for years. We'll try to get a public builder working soon.

## Issues

Please use standard GitHub Issue and Discussion tools for interaction with the developers.

## Tests

Yup, we have 'em. See /test.

We do use a custom test harness for these tests. We'd prefer to use `gtest` at this point, but it will take some time to make the switch.

## CPM Snippet

`geolib` is a library that you'll want to link against, not an application to run. A recommended way to do this with your software is through the [C++ Package Manager](https://github.com/cpm-cmake/CPM.cmake) (or CPM). Add this into your CMake profile and you should be good-to-go.

Drop this in the appropriate place for your application:
```cmake
CPMAddPackage(
        NAME geolib_library
        GIT_REPOSITORY https://github.com/mitre/geodetic_library.git
        VERSION 3.2.7
        GIT_TAG 3.2.7
)
```
In the above, please use a tagged, versioned build of geolib whenever possible.

Then, where you do your linking commands in cmake, also add this:
```cmake
target_include_directories(my-application PUBLIC
        ${geolib_library_INCLUDE_DIRS}
        )
target_link_libraries(my-application geolib_library)
```

Once these are done, your application (`my-application`) will be linked against `geolib_library` and you can call it from your code. Like this:

```c++
#include "geolib/cppmanifest.h"
int main(int argc, char *argv[]) {
  std::cout << "geolib-library build version: " << geolib_idealab::cppmanifest::getVersion() << std::endl;
}
```

## Examples

See also the `examples` binary built by this repository. 