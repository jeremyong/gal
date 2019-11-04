#pragma once

#include "algebra.hpp"
#include "numeric.hpp"

// Templatized routines and operations parameterized by metric signature

namespace gal
{
// Abstractions for compile-time computation of algebraic operations
// Here, we restrict ourselves to associative (but not necessarily commutative) algebras over a
// field

// Template parameters encode the [metric signature](https://en.wikipedia.org/wiki/Metric_signature)
// of the metric tensor underlying the algebra. P: # of basis elements with positive norm V: # of
// basis elements that have negative norm R: # of basis elements that have zero norm Note that
// degenerate metric tensors are not permitted. The metric tensor encoded by this type is
// diagonalized and normalized. Computations can always be expressed using non-orthonormal metrics
// via change-of-basis.
//
// Examples:
//
// metric<3, 0, 0> := 3D Euclidean
// metric<3, 0, 1> := 3D Projective (alt. 4D euclidean)
// metric<3, 1, 0> := Minkowski Spacetime
// metric<4, 1, 0> := Conformal Geometric Algebra w/o change of basis
template <size_t P, size_t V, size_t R>
struct metric
{
    constexpr static size_t p = P;
    constexpr static size_t v = V;
    constexpr static size_t r = R;

    // The dimension of the metric corresponds to the total number of basis elements represented.
    // Note that the metric can induce a larger dimension in an algebra if the multivector space
    // associated with it is [graded](https://en.wikipedia.org/wiki/Graded_multivector_space) (e.g.
    // Clifford Algebra).
    constexpr static size_t dimension = P + V + R;

    // The dot product between two 1-grade basis elements fully defines the Cayley table of an
    // algebra. The ordering convention here is degenerate elements are followed by positive norm
    // elements which are finally followed by negative norm elements.
    [[nodiscard]] constexpr static int dot(size_t lhs, size_t rhs) noexcept
    {
        // The default metric features a diagonal Cayley table
        if (lhs != rhs || lhs < r)
        {
            return 0;
        }
        else if (lhs >= r + p)
        {
            return -1;
        }
        else
        {
            return 1;
        }
    }

    // This is a helper function derived from the Cayley table. Given a basis element, it returns
    // whether the blade supplied contains a vector with a non-zero dot product as the index in the
    // first position, and the dot product in the second position. Returning {-1, .} indicates that
    // the blade and element are orthogonal
    [[nodiscard]] constexpr static std::pair<int, int> intercept(size_t e, size_t blade) noexcept
    {
        // This is the implementation for a diagonal Cayley table
        if (blade & (1 << e))
        {
            return {e, dot(e, e)};
        }
        else
        {
            return {-1, 0};
        }
    }
};

// The specialization with a metric signature as defined above fully specifies a tensor algebra
template <typename Metric>
struct algebra
{
    using metric_t = Metric;
    constexpr static mv<algebra<Metric>, 0, 1, 1> pseudoscalar{
        mv_size{0, 1, 1},
        {},
        {mon{one, zero, 0, 0}},
        {term{1, 0, (1 << metric_t::dimension) - 1}}};
    constexpr static mv<algebra<Metric>, 0, 1, 1> pseudoscalar_inv{
        mv_size{0, 1, 1},
        {},
        {mon{((metric_t::dimension * (metric_t::dimension - 1) / 2 + metric_t::v) % 2 == 0 ? one
                                                                                           : minus_one),
             zero,
             0,
             0}},
        {term{1, 0, (1 << metric_t::dimension) - 1}}};

    // Return a multivector fully spanning the even subalgebra
    constexpr static auto even_mv(uint32_t id)
    {
        constexpr width_t even_dim = 1 << (metric_t::dimension - 1);
        mv<algebra<Metric>, even_dim, even_dim, even_dim> out;
        out.size = mv_size{even_dim, even_dim, even_dim};

        uint32_t e = 0;
        for (width_t i = 0; i != even_dim; ++i)
        {
            out.inds[i]  = ind{id++, one};
            out.mons[i]  = mon{one, one, 1, i};
            out.terms[i] = term{1, i, e};
            e            = next_even(e);
        }
        return out;
    }

    // Return a bivector as an indeterminate multivector
    constexpr static auto bivector_mv(uint32_t id)
    {
        constexpr width_t bivector_dim = metric_t::dimension * (metric_t::dimension - 1) / 2;
        mv<algebra<Metric>, bivector_dim, bivector_dim, bivector_dim> out;
        out.size = mv_size{bivector_dim, bivector_dim, bivector_dim};

        uint32_t e = 0b11;
        for (width_t i = 0; i != bivector_dim; ++i)
        {
            out.inds[i]  = ind{id++, one};
            out.mons[i]  = mon{one, one, 1, i};
            out.terms[i] = term{1, i, e};
            e            = next_even(e);
        }
        return out;
    }

    // For each operation, the static product function returns a generator id and multiplier given
    // two generators. At this point, non-diagonal metric tensors are not supported.

    struct geometric
    {
        [[nodiscard]] constexpr static std::pair<elem_t, int> product(elem_t g1, elem_t g2) noexcept
        {
            if (g1 == 0)
            {
                return {g2, 1};
            }
            else if (g2 == 0)
            {
                return {g1, 1};
            }
            else
            {
                elem_t g     = g1 ^ g2;
                elem_t swaps = 0;

                // The geometric product contracts incident generators based on the metric signature
                // and produces higher-grade tensor products for non-incident generators.
                while (g1 > 0)
                {
                    auto lhs_g        = leading_set_index(g1);
                    auto [index, dot] = metric_t::intercept(lhs_g, g2);
                    if (index == -1)
                    {
                        // Exterior
                        swaps += pop_count(g2 & ((1 << lhs_g) - 1));
                        g2 |= 1 << lhs_g;
                    }
                    else
                    {
                        // Contract
                        swaps += pop_count(g2 & ((1 << index) - 1));
                        if (dot == 0)
                        {
                            return {0, 0};
                        }
                        else if (dot == -1)
                        {
                            ++swaps;
                        }
                        g2 &= ~(1 << index);
                    }

                    g1 &= ~(1 << lhs_g);
                }
                return {g, swaps % 2 == 0 ? 1 : -1};
            }
        }
    };

    struct exterior
    {
        [[nodiscard]] constexpr static std::pair<elem_t, int> product(elem_t g1, elem_t g2) noexcept
        {
            if (g1 == 0)
            {
                return {g2, 1};
            }
            else if (g2 == 0)
            {
                return {g1, 1};
            }
            else
            {
                elem_t intersection = g1 & g2;
                if (intersection != 0)
                {
                    return {0, 0};
                }
                else
                {
                    elem_t g     = g1 | g2;
                    elem_t swaps = 0;

                    while (g1 > 0)
                    {
                        auto lhs_g = leading_set_index(g1);
                        swaps += pop_count(g2 & ((1 << lhs_g) - 1));
                        g1 &= ~(1 << lhs_g);
                        g2 |= 1 << lhs_g;
                    }
                    return {g, swaps % 2 == 0 ? 1 : -1};
                }
            }
        }
    };

    struct contract
    {
        [[nodiscard]] constexpr static std::pair<elem_t, int> product(elem_t g1, elem_t g2) noexcept
        {
            if (g1 == 0)
            {
                return {g2, 1};
            }
            else if (pop_count(g1) > pop_count(g2))
            {
                return {0, 0};
            }
            else
            {
                elem_t swaps = 0;

                while (g1 > 0)
                {
                    auto lhs_g        = gal::leading_set_index(g1);
                    auto [index, dot] = metric_t::intercept(lhs_g, g2);
                    if (index == -1 || dot == 0)
                    {
                        return {0, 0};
                    }
                    else
                    {
                        swaps += pop_count(g2 & ((1 << index) - 1));
                        g1 &= ~(1 << lhs_g);
                        g2 &= ~(1 << index);
                        if (dot == -1)
                        {
                            ++swaps;
                        }
                    }
                }

                return {g2, swaps % 2 == 0 ? 1 : -1};
            }
        }
    };

    struct symmetric_inner
    {
        [[nodiscard]] constexpr static std::pair<elem_t, int> product(elem_t g1, elem_t g2) noexcept
        {
            if (g1 == 0 || g2 == 0)
            {
                return {0, 0};
            }
            else
            {
                auto [g, multiplier] = geometric::product(g1, g2);
                if (multiplier == 0)
                {
                    return {0, 0};
                }

                int grade1       = static_cast<int>(pop_count(g1));
                int grade2       = static_cast<int>(pop_count(g2));
                int target_grade = ::gal::detail::abs(grade1 - grade2);
                int grade        = static_cast<int>(pop_count(g));

                if (grade == target_grade)
                {
                    return {g, multiplier};
                }
                else
                {
                    return {0, 0};
                }
            }
        }
    };
};
} // namespace gal
