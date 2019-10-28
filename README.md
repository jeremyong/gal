# Geometric Algebra Library

[![Build Status](https://travis-ci.org/jeremyong/gal.svg?branch=master)](https://travis-ci.org/jeremyong/gal)

[**GAL**](https://www.jeremyong.com/gal/) (project page link) is a C++17 expression compiler and engine for computing with geometric algebra. It focuses primarily on speed and customizability with the ambition for being suitable for use in production environments where real-time speed is a factor.

Some things that make GAL unique are:

- GAL deals with the intrinsically higher dimensionality of GA by reducing computation to a sparse subset at compile time.
- Often, when dealing with GA constructions, many terms may eliminate but in ways that are not optimizable by the compiler. GAL uses compile-time techniques to ensure that term cancellation which is guaranteed to be correct can be elided.
- GAL separates computation into three layers: data, algebraic expression evaluation, and computation. Each layer is customizable, and examples of why you might want to do so include supporting exotic metrics or number systems or generating shader ISA code (as opposed to CPU evaluation).
- While GAL is template-heavy, fast compile times are a top priority for GAL, which leverages numerous modern C++17 techniques (C++2a soon) to improve compile times.

For the most complete and up-to-date documentation, please refer to the project page [here](https://www.jeremyong.com/gal/). If you wish to learn more about the internals of the library, feel free to continue reading.

## Runtime and Compile-time performance

Early results are showing that GAL is >2x faster than versor but with slower compile times which depend on
the expression complexity.

Runtime and compile time performance is achieved by using the following approach to expression evaluation.

1. Expressions are encoded using standard expression tree templates.
2. Expression inputs are entities that are encoded with flat representations (e.g. an R3 point is just 3 floats and not tied to a multivector representation).
3. Prior to evaluation, expression inputs are expressed in indeterminate multivector form (parameters of each input are expressed via integral tags, not floating-point values). This is encoded using 3 flat compile-time arrays per multivector (storing intederminates, monomials, and polynomials) for fast simplification and evaluation.
4. Expressions are expanded in indeterminate form using polynomial coefficients so that term arithmetic can happen exactly over the field of rationals. A number of techniques are used to put strict upper bounds on compile time memory and CPU usage whever possible.
5. The final indeterminate form is evaluated coefficient by coefficient. For CSE, the compiler is relied on at this time, although future work in doing compile time multivariate polynomial reduction is possible. For maximal throughput, this does not rely on features like `std::tuple` which are known to have poor compilation performance where possible.
6. The results are optionally cast back into the flat entity form which extracts multivector components and applies scaling as appropriate. This operation is also a compile time operation which may cause additional computation to drop out trivially.

## Usage

Being a template-library, GAL is header-only and can be installed by either linking the `gal` interface target via cmake or by copying the files in `src` to a known include path.

To build the tests, you need a C++ compiler (currently untested via MSVC) that supports C++17 or greater.

```sh
# From the root directory of this project
mkdir build
cd build

# Optionally supply release type, flags, etc and pick your favorite generator
cmake .. -G Ninja

# ... or whichever (ideally multicore-friendly) build system you choose
ninja

# Run the tests
./bin/gal_test
```

It is recommend when using clang that `-fno-math-errno` be passed to your compiler's build settings as this was found to be an obstacle for clang to generate optimal code in many circumstances. Many platforms set this as the default, but depending on the math library that is linked on your system, your mileage may vary.

A secondary recommendation is that usage should comfortably enable both `-ffast-math` and `-mfma`. The former
is usually not used to retain finer control over numerical stability. However, using GA improves stability
considerably over traditional linear algebra approaches. The second flag enables fused-multiply-add
instructions which both improves precision and is available on most hardware.

## Motivation

Geometric Algebra promises (and fulfills) a unified algebraic system for manipulating geometric objects
(points, lines, planes, spheres, etc) in a way that is coordinate-free, consistent, and logical. A
typical computer graphics library is often the union of a number of disjoint algebras and algorithms.
Examples include dual quaternions for skinning, quaternions for rotation, 1-up projective space for linear
rotations and translations when interpolation is not needed, and special handling for objects that do not
transform covariantly with the metric (e.g. normal vectors).

A (very) loose way of understanding how GA does this is by encoding operations in a way that loses less
relevant information. For example, consider the cross-product. Once computed, the information that it
originated from two vectors with some orientation relative to each other is lost. Indeed, given nothing
but the coordinates of a normal vector, it is impossible to describe the orientation of the vectors that
produced it. In GA, such information is encoded by higher ordered basis elements beyond the standard ones
proferred in the typical R3 vector space. This is what makes GA a *graded* algebra. The second important
(and necessary) piece that gives GA its power are the operations between elements of this algebra. Without
going into much detail, the geometric product encodes both projections and rejections nicely within the
framework supplies us with the *conjugate* (aka "sandwich") operator to elegantly supply interpolatable (read.
differentiable) rotations and translations of any geometric object. All of this is to say, the state
space of the *code* needed to express a vast range of computation is compressed significantly.

There is a "problem" though, which is that, being a graded algebra, geometric algebras are necessarily
higher-dimensional than the spaces they represent. A 3D space typically requires 3 coordinates but in GA,
this would require 8 (2^3) coordinates to describe a fully general multivector. The 5D conformal space
(conformal meaning that homomorphisms are angle-preserving) requires 32 coordinates! On top of that,
there are often a number of term cancellations as expressions are evaluated (as quantities contract one
one another in ways that are degenerate for example). This results in a higher operation count, all else
being equal. Generally, actual computation is done in smaller embeddings within the full tensor algebra,
and runtime compression of the data is unacceptable.

To combat this, GAL provides a fully compile time expression evaluation system and computational engine
to fully simplify expressions. Perhaps the most egregious example is a CGA (Conformal Geometric Algebra)
point being contracted onto itself (the contraction operator resembles the dot product but is generalized
to act on higher-order elements than just vectors). Under CGA, the contraction of a point to itself is
*exactly* zero. However, if you implemented a contraction operator using standard floating point math,
the compiler would be unable to optimize this as such in general. GAL makes the following code possible:

```c++
#include <gal/cga.hpp>

using point  = gal::cga::point<float>;

float point_norm(point p)
{
    return compute([](auto p)
    {
        // Contract a cga point back onto itself
        // Note that p here contains no actual data! It is just the type that represents a CGA point
        // with internal tags that refer to the values contained in the outer scope p (we locally
        // shadow the variable name for brevity)
        return p >> p;
    }, p);
}
```

The assembly for the routine above looks like the following (compiled with -O1, not even -O2):

```assembly
point_norm(gal::cga::point<float>):
        pxor    xmm0, xmm0
        ret
```

and it should go without saying that `xor` of a register to itself is assembly-shorthand for just zeroing
the register. In other words, all the multiplications of a 5-coordinate vector to itself got optimized away
completely!

## API

The library is still under flux, but for now, the snippet below should give you a decent idea of how to use GAL.
For more complete documentation, please refer to the [project page](https://jeremyong.com/gal).

Example usage:

```c++
#include <gal/cga.hpp>

// Let's work with the projectivized dual space of R3
using point = gal::pga::point<float>;

// Let's construct a few random points.
// Each of these points occupies no more than 12 bytes (+ alignment padding) but
// carry with them compile time representations of the PGA
point p1{2.4f, 3.6f, 1.3f};
point p2{-1.1f, 2.7f, 5.0f};
point p3{-1.8f, -2.7f, -4.3f};

// We issue a computation using the `compute` method which accepts a lambda
plane<float> p = compute<gal::pga::plane<float>>([](auto p1, auto p2, auto p3)
{
    // The p1, p2, and p3 variables here are shadow types of the points residing
    // "in the engine" and we operate with them using any of the operations:
    // operator*        := Geometric product
    // operator^        := Exterior (aka wedge) product
    // operator~        := Reversion
    // operator!        := (Poincare) Dual
    // operator&        := Regressive product (point meet, plane join)
    // operator+        := Vector space addition
    // operator-        := Vector space subtraction
    // operator|        := Symmetric inner product
    // operator>>       := Left contraction
    // operator%        := Sandwich operator (rhs ^ lhs ^ ~rhs)

    // Operations that are permitted are chosen because they respect associativity
    // in the way you would expect.

    // Here we just use the regressive product to construct a plane which passes through
    // the three points.
    return p1 & p2 & p3;
}, p1, p2, p3);

// The results have now been computed and placed into the constructed plane which
// is parametered by the equation ax + by + cz + d = 0
// Here, we check to ensure that the three points we used in the construction do
// in fact lie on the plane.
auto epsilon = /* a small float value */;
CHECK_EQ(p1.x * p.x + p1.y * p.y + p1.z * p.z + p.d, epsilon);
CHECK_EQ(p2.x * p.x + p2.y * p.y + p2.z * p.z + p.d, epsilon);
CHECK_EQ(p3.x * p.x + p3.y * p.y + p3.z * p.z + p.d, epsilon);
```