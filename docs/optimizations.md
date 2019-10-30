# Optimizations

This is a non-exhaustive set of notes to document the various optimizations that are incorporated in GAL.

## Compilation time

- Multivector compile-time representations are templated based on the algebra they reside in (metric + basis) and the size of its constituent compile-time expressions. Representing multivectors as typechains (the traditional approach) is completely untenable for non-trivial work (your compiler will quickly run out of both time and memory).
- Internally, multivectors contain three flat (compile-time constant expression) arrays refering to indeterminates, monomials, and polynomials. Monomials have a count and offset into the indeterminates array, and polynomials (multivector terms) have a count and offset into the monomials array. All arrays are sorted according to well-ordering defined on each (indterminates, monomials, and terms) so that operations like summation, multiplication, etc occur with minimal algorithmic complexity.
- Function template instantiation is kept to a minimum, and use of SFINAE, no matter how convenient, is forbidden.
- The final reification does NOT rely on `constexpr` expansion because this would introduce a function call that the compiler cannot reasonably inline (affects final performance by 10x!). This is done in a type expansion instead at the cost of some compile-time performance.
- Compile-time rational addition and multiplication is guarded against overflows carefully and if an overflow is about to occur, mediant approximation is used.
- Exponentiation is unrolled out compile time using binary right-to-left expansion since the exponent itself is known at compile time. When fractional exponents are present, a single call to `std::pow` is used instead.
- GAL is extremely careful about how its headers are organized so that the parser has minimal work to do. In particular, GAL does whatever it can to avoid the inclusion of system headers unless absolutely necessary. The substitutions made for leaner classes do *not* show up in the final build (they are only used for intermediate computation) and are measurably quicker to compile.

## Runtime

- Because expressions are reduced first over the field of rational coefficients in the polynomial ring of finitely generated indeterminates (the expression inputs), term cancellation occurs exactly.
- When terms would not appear in the final computed result, no instructions are generated.
- If terms drop out in the final result, the type the computation is reified to does not include that term (reflected in `sizeof(result)`).
- When converting an entity result into a concrete entity that has a greater size, a zero is written to the unoccupied terms as cheaply as possible (this is just zero-initialization).
- All expression computation is zero-copy, meaning it is up to the compiler if it wishes to rearrange the data (e.g. in SIMD registers).
- While needless sprinkling of `inline` and inline-forcing attributes is generally ill-advised, it's absolutely critical in GAL. GCC in particular has trouble inlining certain code that is actually completely elidable at compile time due to the expression complexity. As a result, using the force-inline attribute guarantees inlines in places where it was deemed important (and guaranteed not to introduce unnecessary code bloat). In places in GAL where `GAL_FORCE_INLINE` is used, the code will actually both shrink in size *and* accelerate in speed.

## Recommended compiler flags

The philosophy adopted by the author of GAL is that algorithms should be numerically robust such that floating point math can be regarded as associative. While this is sometimes more difficult to achieve, the gains are often well worth it as associativity enables the compiler to more easily vectorize and fuse operations. GA affords a good degree of numerical stability already. For example, a rigid body motion or rotation in GA is minimally parameterized. In contrast, the multiplication of several rotation matrices can quickly get out of hand, as the correlated off-diagonal elements drift breaking normalization. As such, users of GAL are recommended to consider the following compiler flags, at the very least, on a translation unit by translation unit basis.

!!! info "MSVC Usage"
    Unlike GCC and Clang, MSVC provides a single control flag to modify the floating point environment. For MSVC, `/fp:fast` is recommended which enables all recommendations described in more detail below.

The flags mentioned below apply to both GCC and Clang:

| Flag | Purpose | Rationale |
--- | --- | ---
`-fno-signed-zeros` | Disables the signed representation of zeros. | Signed zeros are traditionally used to ensure that a divide by zero will create the correct infinity (there is a big difference between \(+\infty\) and \(-\infty\)). Under projectivization, many operations that would have historically produced an infinity map instead to an ideal. The other case where division by zero may arise is during normalization, but in such a case, the result will be incorrect regardless of the sign of the infinity produced.
`-ffinite-math-only` | The compiler will no longer honor infinities. | As mentioned, infinities under projectivization map to a finite quantity in all cases, and honoring infinities is wasteful.
`-funsafe-math-optimizations` | Permit associative operation reordering. | Floating point expressions like \((x + y) + z\) and \(x + (y + z)\) are not equivalent in general due to rounding errors. It is the opinion of the author that computation should generally be done in ideal precision ranges where possible. GA naturally avoids a number of pathological cases such as adding/multiplying a large number followed by subtracting/dividing by a large number.
`-fno-trapping-math` | Disables trapping floating point exceptions such as divide by zero and overflow. | Normalizing entities such as motors, points, etc. are the only instance where a divide by zero can occur and overflow is conditioned on a reasable coordinate system (object space vs world space, etc). It is preferable to just check for what would have been a possible exception in the isolated cases where it is needed than to enable exception trapping everywhere.

For convenience, the flag `-ffast-math` enables all the flags above and should be usable with GAL.

Additionally, the hardware-dependent flags below are recommended if permitted by your deployment requirements:

| GCC/Clang Flag | MSVC Flag | AMD Support | Intel Support | Rationale |
--- | --- | --- | --- | ---
`-mfma` | `/arch:AVX2` and `/fp:fast` | Piledriver and later (2013) | Haswell and later (2014) | Compilers have gotten exceptionally good at identifying patterns and expressions that benefit from the fused-multiply-add instruction. This instruction not only increases speed but it also increases precision because the entire compound operation is completed at full precision before rounding. These flags correspond to `FMA3`[^1]
`-march=native` | | Any | Any | If you control the final target, using `-march=native` will provide the most aggressive optimizations possible for the target hardware being compiled on.

It goes without saying that all other standard practices apply (prefer 64-bit due to the additional register count, etc).

[^1]: The `FMA4` set of instructions can help in some circumstances, but personally, I haven't found them to be worth the loss in compatibility afforded by using the slightly older (but very good) `FMA3` instructions.

*[SFINAE]: Substitution Failure is Not an Error
*[SIMD]: Single Instruction Multiple Data