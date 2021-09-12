# IncompleteBeta

Simple C++11 implementation of:
1. Beta function `beta(a, b)`
2. Incomplete beta function `incbeta(x, a, b)`
3. Regularized incomplete beta function `incbeta_reg(x, a, b) = incbeta(x, a, b) / beta(a, b)`

The definitions of below special functions can be found in [wiki](https://en.wikipedia.org/wiki/Beta_function). 

The implementation is based on C99 library [incbeta](https://github.com/codeplea/incbeta) which is utilized Lentz's algorithm to solve the continued fraction (continued fraction representation of regularized incomplete beta function can be found [here](https://dlmf.nist.gov/8.17#SS5.p1)). The original solution has been modified to provide better robustness and automatic type promotion of different input parameters using C++11 features.

## Example of usage
```c++
#include "JustMath/Beta.hpp"

...
// define x, a, b
...

/* Beta function */
double r1 = JustMath::beta(a, b);

/* Incomplete Beta function */
double r2 = JustMath::incbeta(x, a, b);

/* Regularized Incomplete Beta function */
double r2 = JustMath::incbeta_reg(x, a, b);
```

## Integration
There are two avaliable options:
1. Just copy `Beta.hpp` into your project
    ```c++
    #include "JustMath/Beta.hpp"
    ```
    and set C++11 support, for example:
    - Compiler flag: `-std=c++11` for GCC and Clang, `/std:c++11` for MSVC
    - Using CMake:
        - Specify compile features for specific target: `target_compile_features(<target> <PRIVATE|PUBLIC|INTERFACE> cxx_std_11)`
        - Set global property: `set(CMAKE_CXX_STANDARD 11 CACHE STRING "The C++ standard to use")`

2. Using CMake:
    ```cmake
    add_subdirectory(IncompleteBeta)
    ...
    target_link_libraries(<YOUR_TARGET_NAME> PRIVATE IncompleteBeta)
    ```