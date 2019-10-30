# Benchmarks

!!! danger "Obligatory disclaimer"
    The author isn't honestly a believer in benchmarks as a proxy for how fast *your* code will run. Use these timings only as a heuristic on which to base your judgement and please, please, benchmark your own code. If you find something isn't fast, please file an issue! Thanks :)

These results were produced on October 28, 2019 for the 3D CGA inverse kinematics benchmark provided by [gabenchmark](https://github.com/ga-developers/ga-benchmark/blob/master/AlgorithmInverseKinematics.cpp) (all benchmarks compiled with `-O2 -mfma` using GCC 10). Prior to each benchmark, CPU scaling was turned off for consistent results.

For possible explanations to explain both the compile-time and runtime behavior of GAL, please refer to the [optimizations](optimizations.md) section of this page.

## Runtime

| Library | Commit Hash | Time | Iterations | Iterations/ns (higher is better) | Speedup (higher is better) |
--- | --- | --- | --- | --- | ---
[Versor](https://github.com/wolftype/versor) | [2c5b455](https://github.com/wolftype/versor/commit/2c5b455187a162bd429a3acf6ffd2f06fb4732ab) | 1138 ns | 608337 | 535/ns | 1x
[GATL](https://github.com/laffernandes/gatl) | [c47ad5](https://github.com/laffernandes/gatl/commit/c47ad5d13c18f0116797d1ea614557cca5e99bb3) | 469 ns | 1505207 | 3209/ns | 5.998x
==[GAL](https://github.com/jeremyong/gal)== | [e3f49ae](https://github.com/jeremyong/gal/commit/e3f49ae28680a89d328b2d37320fe61773c41b56) | 376 ns | 1873236 | **4982/ns** | **9.312x**

This was tested on an Intel i9-9900K. The variable runtime is due to the way [google/benchmark](https://github.com/google/benchmark) works. It will keep running new iterations until the statistical significance of the measurement increases above a predetermined amount. The important metric is the rate in the fourth column. The speedup factor in the final column is computed as the ratio to the throughput of Versor (first row).

## Compilation time

The compilation times for the same example are shown below:

| Compiler | Library | Compile Time (lower is better) |
--- | --- | ---
g++ (GCC) 10.0.0 20191012 (experimental) | Versor | **5.50 s**
g++ (GCC) 10.0.0 20191012 (experimental) | GATL | 20.38 s
g++ (GCC) 10.0.0 20191012 (experimental) | ==GAL== | 12.30 s
clang++ 9.0.0 | Versor | **3.47 s**
clang++ 9.0.0 | GATL | *could not compile*[^1]
clang++ 9.0.0 | ==GAL== | 7.19 s

As with the benchmark, CPU scaling was turned off for the test and all compilations were run from a clean build folder. As you can see, GAL doesn't compile the fastest, but does enjoy an edge over GATL due to the availability of `constexpr` at the time of authorship to avoid slow-to-compile techniques like SFINAE (`std::enable_if`).[^2] Improving compilation time (or at the very least) not regressing is a constant goal for GAL.

## File size

The file size of the compiled test can be a useful indicator of future bloat if used in your library or application. Obviously, you will need to measure this in your own code since the final size may not scale in any predictable way.

The sizes of various segments (`.strtab`, `.text`, and `.symtab`) are included as an indication of what bloat may be attributed to. Smaller sections are omitted for brevity. Also, note that the raw sizes are obviously more important than the percentages. Statistics were taken by using the [google/bloaty](https://github.com/google/bloaty) ELF parser.

!!! tip ""
    The file sizes below include the benchmarking code, so your file savings may be even more dramatic (just remember that no guarantees can be made).

| Library | File size (lower is better) | Percent reduction (higher is better) | `.strtab` | `.text` | `.symtab`
--- | --- | --- | --- | --- | ---
Versor | 282.7 kB | 0% | 92.2 kB (32.6%) | 157 kB (55.8%) | 7.88 kB (2.8%)
GATL | 72.2 kB | 74.461% | 42.2 kB (58.5%) | 12.9 kB (17.8%) | 3.52 kB (4.9%)
==GAL== | **26.5 kB** | **90.626%**| **1.85 kB** (7.0%) | **9.9 kB** (37.4%) | **2.18 kB** (8.2%)

GAL tries very hard not to include headers unnecessarily which has an impact both on compile times and executable bloat.

For reference, the `.strtab` segment refers to the string table (stores null-terminated strings). The `.text` segment refers to the code segment, and the `.symtab` refers to the symbol table (needed to relocate symbols and references).

## Notes for reproducing the benchmark

* GAL provides a mechanism for specifying arbitrary multivectors as types with exact rational coefficients. This means that the conformal point located at \((1, 0, 0)\), represented by \(\mathbf{o} + e_2 + \frac{1}{2}\boldsymbol\infty\) can be represented exactly as opposed to the canonical point representation with \(p_x = p_y = 0.0\) being stored and computed with floating point zeros.
* When this optimization was made, all results were completely accurate up till the final reference value for `Jg_f[4]` for which the published result should be `3.31122e6`. Instead, GAL computes `1.831977e6` for this entry. Reintroducing floating point zeros into the expression immediately recovers the reference result, indicating that this lone result is likely incorrect as stated.
* A check of all 16 other reference results though matches identically (all of which are in more reasonable floating point ranges).
* To compile the GATL implementation of the algorithm, two changes needed to be introduced due to a breaking change with more modern compilers in handling `std::enable_if_t`. Please refer to the following snippets:

!!! warning "Old (doesn't compile under g++10)"
    ```c++
        // GATL
        // lazy_context.hpp:652
        // (whitespace added here for readability on web)
        template<
            typename CoefficientType,
            typename Expression,
            typename = std::enable_if_t<super::stored_inputs_count() != 0>>
        constexpr decltype(auto)
        eval(clifford_expression<CoefficientType, Expression> const &expression) const {
            return detail::eval<
                base_id + 1,
                base_id + (detail::tag_t)super::stored_inputs_count()>(
                    expression, super::stored_inputs_tuple());
        }

        template<
            typename CoefficientType,
            typename Expression,
            typename = std::enable_if_t<super::stored_inputs_count() == 0>>
        constexpr decltype(auto)
        eval(clifford_expression<CoefficientType, Expression> &&) const noexcept {
            return clifford_expression<CoefficientType, Expression>();
        }
    ```

was changed to:

!!! success "New (compiles successfully)"
    ```c++
        // GATL
        // lazy_context.hpp:652
        template<typename CoefficientType, typename Expression>
        constexpr decltype(auto)
        eval(clifford_expression<CoefficientType, Expression> const &expression) const {
            if constexpr (super::stored_inputs_count() > 0)
            {
                return detail::eval<
                    base_id + 1,
                    base_id + (detail::tag_t)super::stored_inputs_count()>(
                        expression, super::stored_inputs_tuple());
            }
            else
            {
                return clifford_expression<CoefficientType, Expression>();
            }
        }
    ```

* As this was an entirely compile-time transformation, it is not believed that either the runtime nor compilation time would be effected much. It would likely improve compilation time if anything due to the removal of an extra template instantiation in certain situations.

## Roadmap

Benchmarking is a fine-art, and while I do not believe they should be misused, they are still helpful for identifying future areas of improvement and, perhaps more importantly, detecting regression. Because of the heavy compile-time nature of the library, small perturbations in the code can (in specific sections), cause dramatic effects downstream. A primary concern of GAL is to quickly work on getting demos working so that timings can be collected for a variety of useful workloads.

***

## Note to other library authors

Please note that nothing on this site is intended to disparage other implementers of GA. The world needs more of them! There are certainly disadvantages of GAL at this time, namely that it is a bit less ergonomic to use at this time, and it's less battle-tested than other libraries. GAL had the benefit of being developed with newer compilers which allowed newer techniques (although there was still difficulty to be sure). In any case, I feel obligated to state a profound appreciation for the implementers and practitioners that paved a way before me.

If I've neglected to include a library or if you think I have made any errors, please don't hesitate to contact me and I'll make sure things get sorted. If there are situations where you think your library or code might perform better, let me know as well. Thanks!

[^1]: Because of the newer clang compiler, GATL couldn't be compiled, possibly due to changes in the C++ standard or because older versions of clang were inadvertedly more permissive. Note that for g++, I opted to test all libraries against the current trunk version, as GAL uses various C++17 techniques handled better in GCC10 and improvements to the compiler generally benefit all libraries.

[^2]: SFINAE refers to a C++ programming technique to restrict resolution of a particular overload to types that match a set of specified constraints. While powerful, SFINAE is known to be one of the worst offenders in compile time cost.

*[MSVC]: Microsoft Visual C++
*[GCC]: GNU C Compiler
*[CGA]: Conformal Geometric Algebra
*[GATL]: Geometric Algebra Template Library
*[ELF]: Executable and Linker Format
*[SFINAE]: Substitution Failure is Not an Error