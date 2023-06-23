# VT C++ Linear Algebra & Kalman Filter Library

This is the repository for C++ implementation of numeric vectors, numeric matrices and linear algebra operations.
Additional standard library implementations are included for C++11 embedded applications without `std` STLs.

This library does not depend on any STLs and built on C++11 standard, optimized for speed and safety, minimized
runtime with a lot of static checking and safety checking during compile-time.

The data structure is built on top of templated C-style static array, so this will avoid heap usages and comply
with small and safety-critical embedded systems. Most of the computations are evaluated at compile-time.

Useful functions to create vectors and matrices, including block matrices, are included.
These data structures can be created directly from array. They are also `constexpr`-compliant.

## Installation

### General Installation

1. Copy header files from *include* folder to your include folder.
2. `#include "vt_linalg"` for linear algebra data structures and tools.
3. `#include "vt_kalman"` for Kalman filter library.
4. See *examples* folder for further usages.
5. See *test* folder for testing and usages.

### PlatformIO

See [vt_linalg registry](https://registry.platformio.org/libraries/vtneil/vt-linalg) on PlatformIO for further
instructions.

## Documentation

To read documentation, visit [https://vt.in.th/docs/vt-linalg/index.html](https://vt.in.th/docs/vt-linalg/index.html).
