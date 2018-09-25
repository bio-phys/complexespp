# Inastemp - Intrinsics as template

Version : 0.1 (29/06/2016)

[TOC]

    
## Introduction

- What is Inastemp?
Inastemp provides a set of `C++` classes to make vectorization with intrinsics easier.
It aims at developing numerical kernels by separating the algorithm from the hardware target.
Inastemp comes with several examples and patterns related to widespread use-cases.

- License
The code is is published under the MIT license.
See the `LICENSE` file or [https://opensource.org/licenses/MIT] for details.

- Target audience and use cases
The library requires basic C++ knowledge about classes, templates, etc.
It is not mandatory to really know what is vectorization all about, but it certainly helps.
For example, the increment of loops is usually tied to the underlying size of the vector which means that one should understand why and where the increment is not equal to 1.
A good way to use Inastemp would be to have a software engineer managing the inclusion and hiding the small complexities of the template process,
and to have the scientists work with basic simple functions templatized against a vector type.

- Overview
The library itself purely consists of headers. However, it comes with CMAKE files.
CMAKE will detect the hardware capacities and turn on some intrinsics by default.
One can turn them on/off using `ccmake` or other configuration tool if needed.
Unit tests and examples can be compiled and the headers installed with the usual CMAKE stage.


### Templates --- Breaking the Genericity/Optimization Opposition

The main point of Inastemp is to factorize the intrinsics code to let the developers implement kernels without knowing the final destination of the code.
It has been common practice to write a kernel for each architecture (partly because, similarly, it has been common practice to develop in C/Fortran).  At this point, Inastemp fills the gap being responsible of the architecture.
When developing a kernel for a specific target we usually make some specific optimizations.
But any optimization usually means being less generic: Optimization by targeting an hardware, by using special instructions, by considering that the input data respect a given pattern, etc.
Using C++ templates, the abstraction allows to have a generic algorithm and to choose only at compile time the optimizations that should be applied.
It offers a huge benefit in term of factorization and maintainability of the code.
```
                 Generic <-----------------------> Optimized
       All architectures                           One architecture
               All input                           Special input
        
Template C++ source-code                           Compiled Template C++ code
```

## Features of Inastemp

- The following x86 SIMD types are currently supported:
    - SSE3, SSSE3, SSE4.1, SSE4.2, AVX, AVX2, AVX512-KNL, AVX512-SKL
- The following Powere PC SIMD types are currently supported:
    - Power-8 Altivec/VMX
- arithmetic operators `*/+-` are provided
- CPU capacities are detected automatically during the CMake stage
- The compiler capacities are detected automatically during the CMake stage
- The library purely contains of headers, no linkage is necessary.
- CPU detection may use Intel®-SDE
- Unit-tests may use Intel®-SDE
- Fast intrinsic `exp()` function (if not supported natively by the compiler)
- Explicit branches vectorization
- several patterns which represent many applications are demonstrated


## Compilation

### Basic compilation

```bash
# Create a Build directory
mkdir Build
# Go into the Build directory
cd Build
# Run cmake - considering Inastemp is the upper directory
cmake ..
# OR with more output messages
VERBOSE=1 cmake ..
# Then make will build the unit test and examples
make
```

### CMake variables - Hardware detection (X86)

There are two hardware detections:
- the instructions supported by the compiler
- the instructions supported by the current CPU

For each existing vector type in Inastemp, the cmake stage runs both tests.
The results are print out if `VERBOSE=1`.
If a vector type can be compiled, it creates the cmake variable `INASTEMP_USE_X` where `X` could be `SSE3`, `AVX2`, etc.
Then, all the vector types that are supported by the CPU are turned to `ON` by default.
So it is possible to turn on specific options even if the CPU does not support them (to execute over intel SDE for example).
But if an options does not exist (even if the hardware seems appropriate) it means that the compiler does not support it.

For example, here is a part of the output of the `ccmake ..` command on a AVX2 CPU and using GCC 6.1:
```
 INASTEMP_USE_AVX                 ON
 INASTEMP_USE_AVX2                ON
 INASTEMP_USE_AVX512KNL           OFF
 INASTEMP_USE_AVX512SKL           OFF
 INASTEMP_USE_SSE3                ON
 INASTEMP_USE_SSE41               ON
 INASTEMP_USE_SSE42               ON
 INASTEMP_USE_SSSE3               ON
```

`AVX512KNL` and `AVX512SKL` are supported by the compiler but not by the hardware, so they are turned `OFF` but could be turn to `ON` if needed.

By turning the cmake variable `INASTEMP_ISDE_CPU` to `ON` the hardware detection is done over intel SDE.
In this case, one can ask Inastemp to check any hardware (passing the appropriate options to isde).

### CMake variables - Hardware detection (Power PC)

Inastemp perform the following test to enable the VMX classes (and to disable all the X86 classes):
```
if(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "ppc64le")
```
If the detection looks incorrect, please check the value returned by `${CMAKE_SYSTEM_PROCESSOR}` and open an issue on our gitlab.

### Using Inastemp as a sub-project with CMake

Considering that the folder `inastemp exist in the root directory of the project, one could use:
```
## inastemp inclusion in CMakeLists.txt of the other project
# tell inastemp we do not want examples/tests
# (will be detected automatically otherwise because it is included as a subproject)
set(INASTEMP_JUST_LIB TRUE)
# Set the C++ standard required by the master project (default inastamp is 11)
set(CMAKE_CXX_STANDARD 11) # Could be 11 or 14 ...
set(CMAKE_CXX_STANDARD_REQUIRED ON) # Could be OFF in specific cases
# add the cmakelist directory
add_subdirectory(inastemp)
# use the filled variables from inastemp
INCLUDE_DIRECTORIES(
         ${INASTEMP_BINARY_DIR}/Src    
         ${INASTEMP_SOURCE_DIR}/Src 
         ${INASTEMP_INCLUDE_DIR}   
    )
# propagate the flags to be able to compile 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${INASTEMP_CXX_FLAGS}")
```

The statement `set(INASTEMP_JUST_LIB TRUE)` asks Inastemp not to compile tests/examples and to limit inclusion of extra warning flags in `INASTEMP_CXX_FLAGS`.
The inclusion of Inastemp as a subproject does not modify the variables of the upper project.
This is why we use the variables `${INASTEMP_CXX_FLAGS}` explicitly to transfer the hardware specific flags.

When Inastemp detects it is included as a subproject, it turns the variable `INASTEMP_AS_SUBPROJECT` to `ON` and minimizes the additional flags.
As a result, the variable `INASTEMP_CXX_FLAGS` will contain a flag for `C++11` and `-funroll-loops` (from `INASTEMP_EXTRA_CXX_FLAGS`), plus the flags needed for the target instructions.
Such that, if one wants to compile only some files with these specific flags, it must use `INASTEMP_CXX_FLAGS` carefully.

### Compilers support

Inastemp was developed and tested using the following compilers on the x86_64 architecture.
- Gcc 6.1 (earlier versions if AVX512/KNL/SKL are not used, like 4.9)
- Clang 3.5
- Intel 16.0
Earlier versions may work as well.

It was also tested on IBM POWER 8 NVL/4023GHz (Openpower) using the following compilers:
- Gcc 6.2.0

- Intel special flags
We pass `-diag-disable 2304 -diag-disable 10121 -diag-disable 10120` to intel compiler to remove the implicit conversion warnings between Inastemp classes.
In fact, we voluntarily want implicit conversion because it helps us to factorize and reduce the number of lines drastically.

User can add the following flags to push for more inlining `-inline-forceinline -no-inline-max-total-size` for intel compiler.
However, it looks like such options do not work for several compiler versions, and thus we do not want to add them by default.

### Multiple hardwares compilation

Inastemp is not designed to compile a binary targeting multiple hardwares (like having an execution path for the different possibilities).
If such need appears, let us know and we might enable it in a next release.

## Best practices

To ensure good performance, developers should obey some rules:
- Never try to use an abstract data type above the vector-type (especially never use `InaVecInterface`).
In fact, the Inastemp functions must be inlined and known at compile time to be efficient.
While in some cases the compiler may know what is behind an abstract type this is not always true and may lead to very poor performance.
- Do not use the SCALAR classes for real development. SSE is available almost everywhere so SCALAR classes should be used only when really needed, like for example at the end of the loop if it is not possible to allocate a loop with a size equal to a multiple of the length of the vector.
- Do not silently convert scalar float/double to vector type in a loop, as shown in the following example:

```cpp
for(){
    // ...
    VecType a = ....
    VecType v = 1.0 * a; // Silent conversion of 1.0
    // Or 
    VecType v = VecType(1.0) * a; // Explicit conversion 
    // In the previous lines 1.0 is converted to VecType
    // It might be optimized by the compiler but please manually move it before the loop
}
// Do
const VecType VecOne = 1;
for(){
    // ...
    VecType a = ....
    VecType v = VecOne * a;
```

- Never update a variable in a lambda function in a if-else structure (see If-else description)

```cpp
a = 1
b = VecType::If(cond).Then([&](){ // take 'a' by ref
        a += 1
        return a;
    }).Else([&](){
        return a; // here 'a' is 1+1
    });

```

## Exp Function (SIMD)

Inastemp provides an implementation of the `exponential` function using intrinsics.
It follows the nice paper _Fast Exponential Computation on SIMD Architectures_
available at
[https://www.researchgate.net/publication/272178514_Fast_Exponential_Computation_on_SIMD_Architectures].
The Remez polynomial is used, and the Scilab file to find out the coefficient is provided.


## Intel® Software Development Emulator (iSDE/SDE)

When using SDE, Inastemp looks for `sde64`, therefore you must update your path or put a symbolic link:
```bash
// Most system:
export PATH=/PATH-TO-SDE-BIN/:$PATH
// On ubuntu if sde is not installed in the system folder, I used
sudo update-alternatives --install /usr/bin/sde64 sde64 /PATH-TO-SDE-BIN/sde64 1
```


## Simple Examples

### Patterns/sumArrays1L
Taken from the pattern examples, we simply sum two arrays.

```cpp
template < class VecType >
void SumArrays(double* __restrict__ dest, const double* __restrict__ src1,
               const double* __restrict__ src2, const int nbToProceed) {

    for (int idx = 0; idx < nbToProceed; idx += VecType::VecLength) {
        const VecType v1(&src1[idx]);
        const VecType v2(&src2[idx]);
        const VecType res = v1 + v2;
        res.storeInArray(&dest[idx]);
        // In one line:
        //VecType(VecType(&src1[idx]) + VecType(&src2[idx]).storeInArray(&dest[idx]);
    }
}
```

### Examples/ParticlesWithBranch

Here we compute the interactions between particles.  The scalar sequential code reads:

```cpp
void ScalarFunction(const int nbParticles, const double* __restrict__ positionsX,
                        const double* __restrict__ positionsY,const double* __restrict__ positionsZ,
                        const double* __restrict__ physicalValues, double* __restrict__ potentials,
                        const double cutDistance, const double constantIfCut){
    for(int idxTarget = 0 ; idxTarget < nbParticles ; ++idxTarget){
        for(int idxSource = idxTarget+1 ; idxSource < nbParticles ; ++idxSource){
            const double dx = positionsX[idxTarget] - positionsX[idxSource];
            const double dy = positionsY[idxTarget] - positionsY[idxSource];
            const double dz = positionsZ[idxTarget] - positionsZ[idxSource];

            const double distance = std::sqrt(dx*dx + dy*dy + dz*dz);
            const double inv_distance = 1/distance;

            if(distance < cutDistance){
                potentials[idxTarget]  += ( inv_distance * physicalValues[idxSource] );
                potentials[idxSource] += ( inv_distance * physicalValues[idxTarget] );
            }
            else{
                potentials[idxTarget]  += ( inv_distance * (physicalValues[idxTarget]-constantIfCut) );
                potentials[idxSource] += ( inv_distance * (physicalValues[idxTarget]-constantIfCut) );
            }
        }
    }
}
```

A part of it is rewritten :

```cpp
template <class VecType>
void VectorizedFunction(const int nbParticles, const double* __restrict__ positionsX,
                        const double* __restrict__ positionsY,const double* __restrict__ positionsZ,
                        const double* __restrict__ physicalValues, double* __restrict__ potentials,
                        const double cutDistance, const double constantIfCut){

    const VecType VecOne = 1;
    const VecType VecConstantIfCut = constantIfCut;
    const VecType VecCutDistance = cutDistance;

    for(int idxTarget = 0 ; idxTarget < nbParticles ; ++idxTarget){

        const VecType targetX = positionsX[idxTarget];
        const VecType targetY = positionsY[idxTarget];
        const VecType targetZ = positionsZ[idxTarget];

        const VecType targetPhysicalValue = physicalValues[idxTarget];
        VecType targetPotential = VecType::GetZero();

        const int lastToCompute = ((nbParticles-(idxTarget+1))/VecType::VecLength)*VecType::VecLength+(idxTarget+1);

        for(int idxSource = idxTarget+1 ; idxSource < lastToCompute ; idxSource += VecType::VecLength){
            const VecType dx = targetX - VecType(&positionsX[idxSource]);
            const VecType dy = targetY - VecType(&positionsY[idxSource]);
            const VecType dz = targetZ - VecType(&positionsZ[idxSource]);

            const VecType distance = VecType(dx*dx + dy*dy + dz*dz).sqrt();
            const VecType inv_distance = VecOne/distance;

            const typename VecType::MaskType testRes = (distance < VecCutDistance);

            const VecType sourcesPhysicalValue = VecType(&physicalValues[idxSource]);

            targetPotential += inv_distance * VecType::IfElse(testRes, sourcesPhysicalValue,
                                                                      sourcesPhysicalValue-VecConstantIfCut);
            const VecType resSource = inv_distance * VecType::IfElse(testRes, targetPhysicalValue,
                                                                       targetPhysicalValue-VecConstantIfCut);
            const VecType currentSource = VecType(&potentials[idxSource]);
            (resSource+currentSource).storeInArray(&potentials[idxSource]);
        }

        potentials[idxTarget] += targetPotential.horizontalSum();
        // .....
```


## Removing Branches

In many cases applications have branches inside their kernels which makes them very difficult to vectorize.
Inastemp provides tools to explicitly manage the branches.  Note that the approach computes all branches!
Therefore, large branches may lead to large execution time.


### Principle

```cpp
// Scalar code
if( a < b )
   res += c * d // Case A()
else
   res -= d // Case B()
```

The code example is equivalent to a multiplication by zero for the `false` case:
```cpp
// scalar code
conditionA = ( a < b ) ? 1.0 : 0
conditionB = !( a < b ) ? 1.0 : 0
res += conditionA * c * d // Case A()
res -= conditionB * d // Case B()
```

It is however faster to use bit masks as follows:
```cpp
// scalar code
conditionATrue = ( a < b ) ? 0xFF..FF : 0
res += conditionA AND (c * d) // Case A()
res -= (NOT conditionATrue) AND d // Case B()
```

Using Inastemp one may explicitely use bit-masks.  Alternatively, the Inastemp if-else function can be used.


### Explicit use of bit masks

- Getting a mask:
```cpp
// Get mask by calling functions
VecType::MaskType mask = VecType().Is{Lower,Greater}[OrEqual](a,b)
VecType::MaskType mask = VecType().Is[Not]Equal(a,b)
// and much more: IsPositive, IsZero,....
// Or based on C++ operator overloading
VecType::MaskType mask = a < b; // a (<|>)[=] b
// and much more: ==, !=, etc..
```

- Using a mask:
The VecType classes provide methods such as AND (&), OR (|), XOR (\^), NOT_AND.


### If-else functions

The VecType classes provide three methods to deal with conditions.
Note that Inastemp computes all sub-results and performs all computation,
ie. even for the false condition.

- IfTrue(condition, value)

```cpp
res = VecType::IfTrue(a < b , c);
// res is set to c if a < b, else to zero
```

The IfTrue method is mainly useful when a part of a formula is conditional:
```cpp
res += a * b + .IfTrue(a < b , c);
```

- IfElse(condition, valueIfTrue, valueIfFalse)

```cpp
res = VecType::IfElse(a < b , c, d);
// res is set to c if a < b, else to d
```

- If-ifelse-else
The vectorizer classes provide an `if` method which allows to return a value or to call a lambda function.

```cpp
// Faster
c = VecType::If(test).Then(a).Else(b);
// Should not be slower using functions
c = VecType::If(test).Then([&](){return a;}).Else([&](){return b;});
```

Also, the lambda/anonymous functions must pass varaibles per references using `[&](){}` and not `[=](){}`.

Here are some examples.

```cpp
z = something()
a = VecType::If(cond1).Then( z )
    .ElseIf(cond2).Then([&](){
        VecType x = func(z)
        return x;
    })
    .Else([&](){
        return z*4;
    });
```

The structure can be expressed as follows:
```cpp
VecType::If( a condition )
                .Then(a value or a lambda function) // Then must follow If
          .ElseIf( a condition ) // Optional else if
                .Then(a value or a lambda function) // Then must follow ElseIf
          .Else(a value or a lambda function) // Optional else

```


### Branch optimization

Knowing that there are branches to be computed one can think about numerical optimizations.
In the next example we compute two times `b * c` because both the branches are executed:
```cpp
a = VecType::If(cond1)
    .Then([&](){
        return b * c + 1;
    })
    .Else([&](){
        return b * c + 2;
    });
```

The developer should simplify the formulas as best as possible to reduce the number of instructions.
```cpp
bc = b * c
a = VecType::If(cond1)
    .Then([&](){
        return bc + 1;
    })
    .Else([&](){
        return bc + 2;
    });
```

### Avoid the computation of all branches when not necessary

In some cases, some branches are much more probable than others such that it could make sense to test if
at least one value in the vector requieres to compute the exeptional path.
It is possible to test if a MaskType is "all true" or "all false" and to act consequently.
But this might reduce the performance while doing less work because of CPU missprediction and un-pipelining.

```cpp
// Scalar code
if( a < b ){ // this is true 99% of the time
   c = exp(a);
}
else {
   c = exp(b);
}

// Vectorized code with all branches computed (2 exp call per scalar value)
c = VecType::IfElse(a < b).Then( a.Exp() ).Else( b.Exp() );

// Improved one
VecType::MaskType mask_ab = (a < b);
if( mask_ab.isAllTrue() ){ // This is higly probable
    c = exp(a)
}
else{ // isAllFalse is not probable here
    c = VecType::IfElse(a < b).Then(exp(a)).Else(exp(b));
}
```


## How to add a new vector type

To add a new intrisics type the following steps need to be taken.
- Ensure the detection of the CPU if you want to automatically enable the SIMD types supported by your hardware.
    See `CMakeModules/GetCpuInfos.cmake` and `CMakeModules/getCpuInfos.cpp`
- Ensure that the compiler supports the seleced SIMD types (mandatory!).
    See `CMakeModules/GetCompilerInfos.cmake`
    - Add a folder `TYPE` in `CMakeModules` with the same test as the others :
    	- `compileTestTYPE.cpp` to test if the compiler can compile the code.
    	- `checkTYPEpe.cpp` to test if the compiler provide operators for this type.
    - Modify the `GetCompilerInfos.cmake` file where the key `(ADD-NEW-HERE)` is located
- Add single/double classes following the same syntax as the other:
	- `InaTYPEOperators.cpp` for the operators (a simple modification of an existing file should do the job)
	- `InaVecTYPE<{double,float>.hpp` to put the class for the new type in double and single precision
- Add your type in the main `CMakeLists.txt` file


## Common Compilation Errors

### `Inlining failed` compilation error using GCC

In case an error message similar to
```cpp
inlining failed in call to always_inline A-SIMD-FUNCTION : target specific option mismatch
 ```
is thrown it means that the flags passed to GCC are not enough to enable the
function `A-SIMD-FUNCTION`. One should carefuly check the version of
SSE/AVX/AVX512 related to the given function. A `VERBOSE=1 make` should show the
passed flags.


### `Ignoring attributes` GCC compilation warning
```cpp
warning: ignoring attributes on template argument ‘__mX {aka __vector(Y) float/double}’ [-Wignored-attributes]
```
This warning may safely be ignored because for the Inastemp implementation this
is not a problem. In fact, Gcc tells us that the type passed by template has its
attributes removed, because we do not use the template type to instantiate
variables.

### `Segmentation fault (core dumped)` while using If/ElseIf/Then anonymous/lambda functions

Please ensure to pass the values per references, that is in the lambda function declarations using `[&](){}` and not `[=](){}`.
In fact, it seems that the compiler is not propagate the memory alignement correctly when copying variable for anonymous functions.
Therefore, you can also tryied your binary with iSDE and may obtain something like:
```
SDE ERROR:  TID: 0 executed instruction with an unaligned memory reference to address 0x1d2e490 INSTR: 0x00042495e: IFORM: VMOVAPS_MEMqq_YMMqq :: vmovaps ymmword ptr [rbx], ymm0
	FUNCTION: _ZNSt14_Function_base13_Base_managerIZ24VectorizedFunctionLambdaI10InaVecAVX2IfEfEviPKT0_S6_S6_S6_PS4_S4_S4_EUlvE_E15_M_init_functorERSt9_Any_dataOS8_St17integral_constantIbLb0EE
	FUNCTION ADDR: 0x000424923
```

# Authors
- Berenger Bramas (berenger.bramas@mpcdf.mpg.de)
- Logo Credit to D.M.
