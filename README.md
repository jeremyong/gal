
<p align="left">
    <img src="./doc/logo.png" alt="GAL" width="300">
</p>

# Geometric Algebra Library

GAL is a C++17 expression compiler and engine for computing with geometric algebra. It focuses primarily on speed and customizability with the ambition for being suitable for use in production environments where real-time speed is a factor.

Some things that make GAL unique are:

- GAL deals with the intrinsically higher dimensionality of GA by reducing computation to a sparse subset at compile time.
- Often, when dealing with GA constructions, many terms may eliminate but in ways that are not optimizable by the compiler. GAL uses compile-time techniques to ensure that term cancellation which is guaranteed to be correct can be elided.
- GAL separates computation into three layers: data, algebraic expression evaluation, and computation. Each layer is customizable, and examples of why you might want to do so include supporting exotic metrics or number systems or generating shader ISA code (as opposed to CPU evaluation).
- While GAL is template-heavy, fast compile times are a top priority for GAL, which leverages numerous modern C++17 techniques (C++2a soon) to improve compile times.

GAL is in the early stages of development, so please stay tuned for more!

## Usage

Being a template-library, GAL is header-only and can be installed by either linking the `gal` interface target via cmake or by copying the files in `src` to a known include path.

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
framework, allow the *conjugate* (aka "sandwich") operator to elegantly supply interpolatable (read.
differentiable) rotations and translations of any geometric object. All of this is to say, the state
space of the *code* needed to express a vast range of computation is compressed significantly.

There is a "problem" though, which is that, being a graded algebra, geometric algebras are necessarily
higher-dimensional than the spaces they represent. A 3D space typically requires 3 coordinates but in GA,
this would require 8 (2^3) coordinates to describe a fully general multivector. The 5D conformal space
(conformal meaning that homomorphisms are angle-preserving) requires 32 coordinates! On top of that,
there are often a number of term cancellations as expressions are evaluated (as quantities contract one
one another in ways that are degenerate for example). This results in a higher operation count, all else
being equal.

To combat this, GAL provides a fully compile time expression evaluation system and computational engine
to fully simplify expressions. Perhaps the most egregious example is a CGA (Conformal Geometric Algebra)
point being contracted onto itself (the contraction operator resembles the dot product but is generalized
to act on higher-order elements than just vectors). Under CGA, the contraction of a point to itself is
*exactly* zero. However, if you implemented a contraction operator using standard floating point math,
the compiler would be unable to optimize this as such in general. GAL makes the following code possible:

```c++
using gal::scalar;
using gal::cga::point;

float construct_plane(point<float> p)
{
    gal::engine engine{p};

    return engine.compute<scalar<float>>([](auto p)
    {
        // Contract a cga point back onto itself
        return p >> p;
    });
}
```

The assembly for the routine above looks like the following (compiled with -O1, not even -O2):

```assembly
construct_plane(gal::cga::point<float>):
        pxor    xmm0, xmm0
        ret
```

and it should go without saying that `xor` of a register to itself is assembly-shorthand for just zeroing
the register. In other words, all the multiplications of a 5-coordinate vector to itself got optimized away
completely!

## API

The library is still under flux, but for now, the file `test/test_engine.cpp` should give you a decent idea of how to use GAL.

Example usage:

```c++
// Let's work with the projectivized dual space of R3
using point = gal::pga::point<float>;

// Let's construct a few random points.
// Each of these points occupies no more than 12 bytes (+ alignment padding) but
// carry with them compile time representations of the PGA
point p1{2.4f, 3.6f, 1.3f};
point p2{-1.1f, 2.7f, 5.0f};
point p3{-1.8f, -2.7f, -4.3f};

// Any time we wish to evaluate an expression, we create an instance
// of an engine and pass all objects we wish to compute with to the constructor.
// Like other aspects of the library, the engine is a compile time construct
// and occupies no space (unless you choose to store it).
gal::engine engine{p1, p2, p3};

// We issue a computation using the `compute` method which takes a type parameter
// to indicate the return type and a lambda which will perform the computation.
auto plane = engine.compute<gal::pga::plane<float>>([](auto p1, auto p2, auto p3)
{
    // The p1, p2, and p3 variables here are shadow types of the points residing
    // "in the engine" and we operate with them using any of the operations:
    // operator*        := Geometric product
    // operator^        := Exterior (aka wedge) product
    // operator~        := Reversion
    // operator!        := (Poincare) Dual
    // operator|        := Join
    // operator+        := Vector space addition
    // operator-        := Vector space subtraction
    // operator>>       := Left Contraction
    // conjugate(a, b)  := a ^ b ^ ~a

    // Operations that are permitted are chosen because they respect associativity
    // in the way you would expect.

    // Here we just use the join operator to construct a plane which passes through
    // the three points.
    return p1 | p2 | p3;
});

// The results have now been computed and placed into the constructed plane which
// is parametered by the equation ax + by + cz + d = 0
// Here, we check to ensure that the three points we used in the construction do
// in fact lie on the plane.
auto epsilon = /* a small float value */;
CHECK_EQ(p1.x * p.x + p1.y * p.y + p1.z * p.z + p.d, epsilon);
CHECK_EQ(p2.x * p.x + p2.y * p.y + p2.z * p.z + p.d, epsilon);
CHECK_EQ(p3.x * p.x + p3.y * p.y + p3.z * p.z + p.d, epsilon);
```


## TODO

- Implement additional geometric objects (even sub-algebra versors, etc)
- Implement additional operators (exponential maps and logarithms)
- Implement the Conformal Geometric Algebra
- Add sample applications and provide benchmark
- Add an additional engine that compiles expressions to SPIR-V or shader code
- Add additional documentation
- Support arithmetic constraints between finite field elements