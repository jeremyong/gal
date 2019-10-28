# Optimizations

This is a non-exhaustive set of notes to document the various optimizations that are incorporated in GAL.

## Compilation time

- Multivector compile-time representations are templated based on the algebra they reside in (metric + basis) and the size of its constituent compile-time expressions. Representing multivectors as typechains (the traditional approach) is completely untenable for non-trivial work (your compiler will quickly run out of both time and memory).
- Internally, multivectors contain 3 flat (compile-time constant expression) arrays refering to indeterminates, monomials, and polynomials. Monomials have a count and offset into the indeterminates array, and polynomials (multivector terms) have a count and offset into the monomials array. All arrays are sorted according to well-ordering defined on each (indterminates, monomials, and terms) so that operations like summation, multiplication, etc occur with minimal algorithmic complexity.
- Function template instantiation is kept to a minimum, and use of SFINAE, no matter how convenient, is forbidden.
- The final reification does NOT rely on `constexpr` expansion because this would introduce a function call that the compiler cannot reasonably inline (affects final performance by 10x!). This is done in a type expansion instead at the cost of some compile-time performance.
- Compile-time rational addition and multiplication is guarded against overflows carefully and if an overflow is about to occur, mediant approximation is used.
- Exponentiation is unrolled out compile time using binary right-to-left expansion since the exponent itself is known at compile time. When fractional exponents are present, a single call to `std::pow` is used instead.

## Runtime

- Because expressions are reduced first over the field of rational coefficients in the polynomial ring of finitely generated indeterminates (the expression inputs), term cancellation occurs exactly.
- When terms would not appear in the final computed result, no instructions are generated.
- If terms drop out in the final result, the type the computation is reified to does not include that term (reflected in `sizeof(result)`).
- When converting an entity result into a concrete entity that has a greater size, a zero is written to the unoccupied terms as cheaply as possible (this is just zero-initialization).
- All expression computation is zero-copy, meaning it is up to the compiler if it wishes to rearrange the data (e.g. in SIMD registers).