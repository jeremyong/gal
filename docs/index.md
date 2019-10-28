# GAL

[![Build Status](https://travis-ci.org/jeremyong/gal.svg?branch=master)](https://travis-ci.org/jeremyong/gal)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

[GAL](https://github.com/jeremyong/gal) is a realtime suitable C++17 library designed to be simultaneously competitive with and complementary to traditional methods for computing geometry (e.g. linear algebra, vector spaces, quaternions, dual quaternions). As the majority of its work is done at compile-time, the library is naturally "header-only" due to the implicit inlining requirements.

For computing in 3D, Geometric Algebra (henceforth, just *GA*) promotes the familiar 3-coordinate space to an 8 dimensional one (or 16 if projectivized) in order to directly encode the various geometric entities naturally and in a way that promotes uniform expressiveness. For example, formulae involving transformations (rotations, translations, rigid body motion) and incidence (projections from one entity to another, metric measurements) are all expressed in a *uniform* way between points, lines, planes, etc. (as opposed to needing a bespoke formula for each situation).

!!! Tip "GA? Is it fast though?"
    GAL is pretty fast! Check out the [benchmarks](benchmarks.md) page.