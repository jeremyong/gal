# Coding Conventions

This document describes a number of patterns you should expect to encounter as either a user or a contributor of the library.

## General Guidelines

- Like the STL, the canonical casing is `snake_case` for classes, structs, unions, and variables.
- GAL never allocates on the heap.
- No virtual dispatch is used.
- No inheritance is used.
- Don't-pay-for-what-you-don't-use is a manifest property and no *runtime* overhead is imposed unintentionally.
    - Zero-initialization is not provided for runtime types
    - Division-by-zero is not checked (i.e. during normalization)
- Every feature must have compile times taken into consideration.
- User facing APIs are in the `gal` namespace or in one of the model namespaces.
- Hidden functionality or implementation detail is contained in a nested `detail` namespace.
- Methods are decorated appropriately with attributes and qualifiers (`[[nodiscard]]`, `noexcept`, `const`).

## Notes to Developers

- The largely header-only aspect of the library is due to the heavy reliance on `constexpr` which implicitly inlines the function.
- Type names are often kept short and abbreviated to reduce parsing time. This isn't done unless profiling shows it to be necessary.
- The codegen of the final expression reification is extremely sensitive to changes in the code, so if operating in this area, a benchmark is necessary to verify no regression took place.
- When juggling templates, we endeavor (greatly!) to avoid template instantiations wherever possible. By the same token, we avoid relying on `std::tuple`, SFINAE, recursion, and other techniques known to have adverse effects on compile time wherever possible.
- GAL tiptoes around undefined behavior as much as possible, preferring to be standards conforming. The notable exception is the reliance of well-behaved alignment for structs and data members colocated within a union. This is a reasonable expectation for all modern compilers today.

*[STL]: Standard Template Library
*[SFINAE]: Subtitution Failure Is Not An Error