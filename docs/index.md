# GAL

[GAL](https://github.com/jeremyong/gal) is a realtime suitable C++17 library designed to be simultaneously competitive with and complementary to traditional methods for computing geometry (e.g. linear algebra, vector spaces, quaternions, dual quaternions). As the majority of its work is done at compile-time, the library is naturally "header-only" due to the implicit inlining requirements.

For computing in 3D, Geometric Algebra (henceforth, just "GA") promotes the familiar 3-coordinate space to an 8 dimensional one (or 16 if projectivized) in order to directly encode the various geometric entities naturally and in a way that promotes uniform expressiveness. For example, formulae involving transformations (rotations, translations, rigid body motion) and incidence (projections from one entity to another, metric measurements) are all expressed in a *uniform* way between points, lines, planes, etc. (as opposed to needing a bespoke formula for each situation).

## Getting Started

* `git clone git@github.com:jeremyong/gal.git` - Clone the [project](https://github.com/jeremyong/gal) (directly, as a submodule in your project, via CMake's [FetchContent](https://cmake.org/cmake/help/latest/module/FetchContent.html) module, etc)
* `gal` - Link the `gal` interface target if incorporating the library into your own code (will disable test and sample compilation)
* `mkdir build && cd build && cmake .. -G Ninja && ninja` - Build the tests
* `./bin/gal_test` - From the build folder, run the tests

When building as a standalone project, the following CMake options are available:

| Option | Default | Description
--- | --- | ---
`GAL_TESTS_ENABLED` | `ON` | Compiles the tests
`GAL_SAMPLES_ENABLED` | `ON` | Compiles the samples (none as of yet, stay tuned!)
`GAL_PROFILE_COMPILATION_ENABLED` | `OFF` | Enables timing data generation (traces if using clang, reports if using gcc)

If using CMake to integrate GAL into your project, here's a quick snippet you can use (requires CMake 3.14 or above):

```cmake
include(FetchContent)
FetchContent_Declare(
    gal
    GIT_REPOSITORY https://github.com/jeremyong/gal.git
    GIT_TAG master # or specific commit hash
)
FetchContent_MakeAvailable(gal)

target_link_libraries(your_target PRIVATE gal)
```

## Requirements

* A working C++17 compiler (currently gcc and clang are tested but MSVC will come soon)
* [CMake](https://cmake.org) 3.14+ if building the tests
* (optional) [Ninja](https://ninja-build.org/) Backend build system

## Usage

### Pick A Model

To use GAL, you need to minimally include a single header corresponding to the model you wish to work with. The model defines the [metric signature](https://en.wikipedia.org/wiki/Metric_signature) as well as the dimensionality of the space you wish to represent. Typically, this will either be the Euclidean plane or the fundamental Euclidean 3-space. If you installed or linked to the library through conventional means, you would include things with the `gal` path prefix (as in, `#include <gal/pga.hpp>` for example).

| Header | Space | Namespace | Description
--- | --- | --- | ---
`gal/ega.hpp` | \(\mathbb{R}_{3, 0, 0}\) | `gal::ega` | Canonical 3D Euclidean space
`gal/pga2.hpp` | \(\mathbf{P}(\mathbb{R}^*_{2, 0, 1})\) | `gal::pga2` | Projectivized 2D Euclidean plane
`gal/pga.hpp` | \(\mathbf{P}(\mathbb{R}^*_{3, 0, 1})\) | `gal::pga` | Projectivized 3D Euclidean space
`gal/cga2.hpp` | \(\mathbf{C}(\mathbb{R}_{3, 1, 0})\) | `gal::cga2` | Conformalized 2D Euclidean plane
`gal/cga.hpp` | \(\mathbf{C}(\mathbb{R}_{4, 1, 0})\) | `gal::cga` | Conformalized 3D Euclidean plane

All functions and entities specific to the model will be defined in the namespace according to the table above. General functions and operations will be defined in the `gal` namespace.

### Create some entities

Regardless of the space you choose, you will generally have access to the geometric entities you would expect, each templatized based on the precision you wish to work with.

| Entity | Class
--- | ---
Points | `point<T>`
Lines | `line<T>`
Planes | `plane<T>`
Rotors | `rotor<T>`
Translator | `translator<T>`
Motor | `motor<T>`

As these entities have representations that are model-dependent, they reside in the model namespace from the table above. If a class is missing, it is likely not implemented as of yet (most work will happen in the \(\mathbf{P}(\mathbb{R}^*_{3, 0, 1})\)) so please file an issue if there's something you need (or even better, submit a PR!).

All the entities above are representations of *Euclidean* entities and take on a different shape internally (at compile-time) depending on the space and model you're working in.

All entities support initialization the way you'd expect. For example `point<double> p{0, 1, 2}` defines the point at \((0, 1, 2)\) when using any of the models that represent Euclidean 3-space. They also support indexed access (i.e. `p[1] == 1` will evaluate true) and in some cases, accessors are aliased using nested unions for convience. The point coordinates, for example, may be referred to via `u`, `v`, and `w`, or even `s` and `t` depending on the model dimensionality.

Where appropriate, entities may support the `normalize` method, which mutates the entity as appropriate for its type. A divide-by-zero check is not performed. In all cases in fact, the library adopts the philosophy that the user should not pay for what the user does not use.

Entities are stored sparsely regardless of the model they are in with *no additional overhead* for multivector basis-index bookkeeping. For example, `sizeof(point<float>)` is 12 (3 floats) regardless of whether you are working with CGA (where the point is represented as \(\mathbf{o} + p_x\mathbf{e_x} + p_y\mathbf{e_y} + p_z\mathbf{e_z} + \frac{1}{2}p^2\boldsymbol\infty\)) or in EGA (where the point is simply \(p_x\mathbf{e_x} + p_y\mathbf{e_y} + p_z\mathbf{e_z}\)).

### Evaluating expressions

Expressions are evaluated lazily at compile-time in a manner that maximally promotes [CSE](https://en.wikipedia.org/wiki/Common_subexpression_elimination). Thus, we sacrifice a little bit of code-authorship convenience to gain a significant speedup at runtime.

A computation will take the form `compute([](...) { ... }, entities...)` where the first argument is a lambda and the second argument are the entities you created in the prior step.

A simple calculation might look like:

```c++
using namespace gal::pga;
using gal::compute;

point<> p1{1, 0, 0}; // The default value type is `float`
point<> p2{0, 1, 1};

// Take the regressive product of the two points to construct the line passing through them
line<> l = compute([](auto p1, auto p2) { return p1 & p2; }, p1, p2);

// You can verify that l defines Plucker coordinates that pass through (1, 0, 0) and (0, 1, 1)
```

In the lambda (the first argument to compute), we do not capture variables, and we do not produce side-effects. The lambda itself is never actually evaluated (instead, `decltype` is used to compute the final result based solely on the type the lambda returns). In the body of the lambda, we shadow the outside variable names for convenience. The actual types that are passed to the lambda are building blocks of a template expression tree. Thus, operations that are supported in the body of a compute lambda are *not* supported outside the lambda.

```c++
point<> p1{1, 0, 0};
point<> p2{0, 1, 1};

line<> l = p1 & p2; // Fails to compile!
```

This restriction may be relaxed in the future so that such operations create an immediate computation context for prototyping convenience.

For convience, computations may return multiple results that will be reified to the same number of results.

```c++
point<> p1{1, 0, 0};
point<> p2{0, 1, 1};
point<> p3{1, 0, 1};

// Make a line from p1 to p2, and a plane containing p1, p2, and p3
auto&& [p12, p123] = compute(
    [](auto p1, auto p2, auto p3) {
        auto p12 = p1 & p2;
        auto p123 = p12 & p3;
        return std::make_tuple(p12, p123);
    }, p1, p2, p3);
```

The binary operations supported in the compute lambda are as follows:

| Symbol | Operation | Description
--- | --- | ---
`*` | \(ab\) | The fundamental geometric product
`^` | \(a \wedge b \) | The exterior (wedge or outer) product
`&` | \(a \vee b = (a^\star \wedge b^\star)^\star\) | The regressive product
`>>` | \(a \rfloor b \) | The left contraction product
`|` | \(a \cdot b\) | The symmetric inner product
`%` | \(ba \tilde b\) | The sandwich operator
`+` | \(a + b\) | Multivector sum
`-` | \(a + b\) | Multivector difference
`scalar_product` | \(\langle ab\rangle_0\) | Scalar product

And the following table describes the unary operations

| Symbol | Operation | Description
--- | --- | ---
`~` | \(\tilde a\) | Reversion
`!` | \(a^*\) | The Poincare dual map
`-` | \(a + b\) | Multivector negation
`extract<uint8_t... E>(a)` | \(\Sigma_{\{i \in E\}} a_i\) | Extract a specified set of components into a new multivector

Generally, the operations above work with multivectors, The main exception is the use of `+`, `-`, `*`, and `/` in order to shift or scale a multivector by a compile-time constant. For this, one (and at most one) of the operands must be of type `frac<int, int>`. For example `frac<1, 2> * a` would divide the multivector `a` by 2. Equivalently, this could be done with `a / frac<2>` as you would expect (the denominator defaults to 1).

### Reading results

As we saw in the above example, we can either return from a computation a single result, or just as easily a variadic number of results. What is happening behind the scenes here is an implicit cast whenever you specify the type. It starts
as a generic `entity` type which is templatized on the data it contains. For example, if we did a computation which produced a result of the form \(a + be_{01}\), the computation would produce an entity of type `entity<A, T, 0, 0b11>` where `A` is the algebraic model we are using, `T` is the value type, and `0` and `0b11` (3 expressed as binary) are bitfields encoding the coordinates of the multivector. The entity type can be worked with directly using the index operator `[size_t]`. Alternatively, it can be cast at any point to a concrete type like so:

```c++
point<float> p{an_entity};
```

There is no runtime cost associated with the cast unless the source entity is aliased in a manner the compiler must do a copy. No RTTI or inheritance is used. You cannot cast the type to a different model or precision (at this time).

Entities, do *not* need to be casted into concrete types before being used in a computation. For example:

```c++
point<> p1{1, 0, 0};
point<> p2{0, 1, 1};
point<> p3{1, 0, 1};

// This returns a line, but lets say we don't know that
auto p12 = compute(
    [](auto p1, auto p2) {
        return p1 & p2;
    }, p1, p2);

static_assert(p12.size() == 6, "p12 is a 6 dimensional bivector");

// We use p12 without caring about its shape (same efficiency as above) before finally
// storing it as a plane. Note that we could have continued to proceed without storing
// it as a concrete type.
plane<> p123 = compute(
    [](auto p12, auto p3)
        return p12 & p3;
    }, p12, p3);
```

If needed, a custom multivector can be easily created manually in the geometric model of your choosing. For example:

```c++
using algebra = gal::pga;
entity<algebra, double, 0b11, 0b101> my_multivector{3.2, 1.2};
```

This creates the quantity \(3.2e_{01} + 1.2_{02})\ and can be used in a compute context like any other entity (concrete or otherwise). The basis elements are expressed as a bitfield with the lower indices corresponding to the least significant bits. It is important that they be specified *in ascending order* as this is not currently checked. Internally, all multivectors, polynomials, and indeterminates are kept sorted to achieve optimal compiler throughput and many algorithms may break if this total ordering is not respected.

## Roadmap

(not ordered)

- Support MSVC
- Support transcendental function usage in lazily evaluated contexts
- Provide samples and visualization framework
- Provide examples using dual and complex numbers as the field type
- Provide more instructional documentation for newer practitioners of GA
- Support much longer evaluation contexts (compile-time acceleration)
- Integrate CI/automated testing

## Project layout

    docs/
        index.md  # This documentation homepage.
        ...       # Other markdown pages, images and other files.
    public/                     # Added to your include path (public headers)
        gal/
            algebra.hpp         # Routines for manipulating multivectors (product, sum, negation, etc)
            algorithm.hpp       # Compile-time routines (i.e. sorting, rearrangement)
            cga.hpp             # Provides conformal geometric algebra
            cga2.hpp            # Provides 2D conformal geometric algebra (aka compass ruler algebra)
            ega.hpp             # Provides 3D geometric algebra
            engine.hpp          # Defines various mechanisms for evaluating expressions at runtime
            entity.hpp          # Describes the statically-typed representation of runtime multivectors
            expression_debug.hpp    # Debug facilities
            expression.hpp      # Expression template interface
            format.hpp          # Various string-conversion routines
            geometric_algebra.hpp   # Implements the various products and operations defined in GA
            null_algebra.hpp    # Routines for converting to and from the null-basis
            numeric.hpp         # Compile time numeric facilities (rational numbers, fast pow, etc)
            pga.hpp             # Provides the 3D projective geometric algebra P(R3*)
            pga2.hpp            # Provides the 2D projective geometric algebra P(R2*)
    samples/
        main.cpp    # Primary entrypoint (coming soon!)
    test/
        ...         # Various files to test different functionality
