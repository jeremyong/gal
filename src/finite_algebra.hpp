#pragma once

#include "ring_generator.hpp"

#include <type_traits>

namespace gal
{
// Abstractions for compile-time computation of algebraic operations
// Here, we restrict ourselves to associative (but not necessarily commutative) algebras over a field

// Template parameters encode the [metric signature](https://en.wikipedia.org/wiki/Metric_signature) of
// the metric tensor underlying the algebra.
// P: # of basis elements with positive norm
// V: # of basis elements that have negative norm
// R: # of basis elements that have zero norm
// Note that degenerate metric tensors are not permitted.
// The metric tensor encoded by this type is diagonalized and normalized. Computations can always be
// expressed using non-orthonormal metrics via change-of-basis.
//
// Examples:
//
// metric<3, 0, 0> := 3D Euclidean
// metric<3, 0, 1> := 3D Projective (alt. 4D euclidean)
// metric<3, 1, 0> := Minkowski Spacetime
// metric<4, 1, 0> := Conformal Space w/o change of basis
template <size_t P, size_t V, size_t R>
struct metric
{
    constexpr static bool is_diagonal = true;
    constexpr static size_t p = P;
    constexpr static size_t v = V;
    constexpr static size_t r = R;

    // The dimension of the metric corresponds to the total number of basis elements represented.
    // Note that the metric can induce a larger dimension in an algebra if the multivector space
    // associated with it is [graded](https://en.wikipedia.org/wiki/Graded_multivector_space) (e.g. Clifford Algebra).
    constexpr static size_t dimension = P + V + R;

    // The dot product between two 1-grade basis elements fully defines the Cayley table of an algebra. This metric may
    // be specialized to support non-diagonal 1-grade Cayley tables (e.g. Conformal Geometric Algebra)
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

    // This is a helper function derived from the Cayley table. Given a basis element, it returns whether the blade
    // supplied contains a vector with a non-zero dot product as the index in the first position, and the dot product in
    // the second position. Returning {-1, .} indicates that the blade and element are orthogonal
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

// The coefficient of a term is a linear combination of factors which are kept separate for the purposes of determining
// an optimal computational strategy. An `monomial` is a product of `factor`s.
// Invariant: the factors are sorted in order
template <typename Q, typename... Generators>
struct monomial
{
    constexpr static size_t size  = sizeof...(Generators) + 1;
    constexpr static int degree   = (Generators::degree + ...);
    constexpr static bool is_zero = false;
    using rational_t              = Q;
};

template <typename Q>
struct monomial<Q>
{
    constexpr static size_t size = 1;
    constexpr static int degree = 0;
    constexpr static bool is_zero = false;
    using rational_t              = Q;
};

template <typename Q, typename G>
struct monomial<Q, G>
{
    constexpr static size_t size  = 1;
    constexpr static int degree   = G::degree;
    constexpr static bool is_zero = false;
    using tag_t                   = typename G::tag_t;
    using rational_t              = Q;
};

template <>
struct monomial<zero>
{
    constexpr static size_t size  = 0;
    constexpr static int degree   = 0;
    constexpr static bool is_zero = true;
    using rational_t              = zero;
};

template <typename M>
[[nodiscard]] constexpr bool is_monomial_identity() noexcept
{
    return std::is_same<M, monomial<one>>::value;
}

template <typename T, size_t O, typename D1, typename D2>
[[nodiscard]] constexpr auto operator*(generator<T, D1, O> lhs, generator<T, D2, O> rhs)
{
    using lhs_t      = decltype(lhs);
    using rhs_t      = decltype(rhs);
    constexpr auto d = D1::value + D2::value;
    if constexpr (d == O)
    {
        return zero_generator{};
    }
    else
    {
        return generator<T, degree<d>, O>{};
    }
}

namespace detail
{
    template <typename M, typename... I, typename G1, typename... J, typename G2, typename... K>
    [[nodiscard]] constexpr auto
    product(monomial<M, I...> accum, monomial<one, G1, J...> lhs, monomial<one, G2, K...> rhs)
    {
        static_assert(!G1::tag_t::untagged && !G2::tag_t::untagged);
        if constexpr (G1::tag_t::untagged && G2::tag_t::untagged)
        {
            // reduce an extra recursion if both factors are untagged
            if constexpr (sizeof...(J) == 0)
            {
                return monomial<M, I..., G1, G2, K...>{};
            }
            else if constexpr (sizeof...(K) == 0)
            {
                return monomial<M, I..., G1, G2, J...>{};
            }
            else
            {
                return product(monomial<M, I..., G1, G2>{}, monomial<one, J...>{}, monomial<one, K...>{});
            }
        }
        else if constexpr (G1::tag_t::untagged || less<typename G1::tag_t, typename G2::tag_t>::value)
        {
            if constexpr (sizeof...(J) == 0)
            {
                return monomial<M, I..., G1, G2, K...>{};
            }
            else
            {
                return product(monomial<M, I..., G1>{}, monomial<one, J...>{}, rhs);
            }
        }
        else if constexpr (G2::tag_t::untagged || less<typename G2::tag_t, typename G1::tag_t>::value)
        {
            if constexpr (sizeof...(K) == 0)
            {
                return monomial<M, I..., G2, G1, J...>{};
            }
            else
            {
                return product(monomial<M, I..., G2>{}, lhs, monomial<one, K...>{});
            }
        }
        else
        {
            // Both G1 and G2 are tagged and represent the same factor
            using g = decltype(G1{} * G2{});
            if constexpr (sizeof...(J) == 0 || sizeof...(K) == 0)
            {
                if constexpr (g::is_zero)
                {
                    return monomial<zero>{};
                }
                else
                {
                    return monomial<M, I..., g, J..., K...>{};
                }
            }
            else
            {
                if constexpr (g::is_zero)
                {
                    return monomial<zero>{};
                }
                else
                {
                    return product(monomial<M, I..., g>{}, monomial<one, J...>{}, monomial<one, K...>{});
                }
            }
        }
    }
} // namespace detail

// Order preserving monomial multiplication
template <typename Q1, typename Q2, typename... G1, typename... G2>
[[nodiscard]] constexpr auto operator*(monomial<Q1, G1...>, monomial<Q2, G2...>)
{
    if constexpr (Q1::is_zero || Q2::is_zero)
    {
        return monomial<zero>{};
    }
    else if constexpr (sizeof...(G1) == 0 || sizeof...(G2) == 0)
    {
        using Q = decltype(Q1{} * Q2{});
        return monomial<Q, G1..., G2...>{};
    }
    else
    {
        using Q = decltype(Q1{} * Q2{});
        return detail::product(monomial<Q>{}, monomial<one, G1...>{}, monomial<one, G2...>{});
    }
}

template <typename Q, typename... G>
[[nodiscard]] constexpr auto operator-(monomial<Q, G...>)
{
    return monomial<typename Q::neg_t, G...>{};
}

template <typename Q1, typename Q2, typename... G>
[[nodiscard]] constexpr auto operator+(monomial<Q1, G...> lhs, monomial<Q2, G...> rhs) noexcept
{
    using Q = decltype(Q1{} + Q2{});
    if constexpr (Q::is_zero)
    {
        return monomial<zero>{};
    }
    else if constexpr (Q1::is_zero)
    {
        return rhs;
    }
    else if constexpr (Q2::is_zero)
    {
        return lhs;
    }
    else
    {
        return monomial<Q, G...>{};
    }
};

template <typename Q1, typename Q2, typename... G>
[[nodiscard]] constexpr auto operator-(monomial<Q1, G...>, monomial<Q2, G...>) noexcept
{
    using Q = decltype(Q1{} - Q2{});
    if constexpr (Q::is_zero)
    {
        return monomial<zero>{};
    }
    else
    {
        return monomial<Q, G...>{};
    }
};

template <int N, int D, typename Q, typename... G>
[[nodiscard]] constexpr auto operator*(rational<N, D> lhs, monomial<Q, G...>) noexcept
{
    if constexpr (N == 0 || Q::is_zero)
    {
        return monomial<zero>{};
    }
    else
    {
        return monomial<decltype(lhs * Q{}), G...>{};
    }
}

// Unlike other functions in this module, `simplify` does not assume its input is ordered with respect to
// any monomial ordering. It relies on the addition operator to coalesce like-terms and the output is
// respect to be well-ordered.
template <template <typename, typename...> typename T, typename C, typename... I>
[[nodiscard]] constexpr auto simplify(T<C, I...>) noexcept
{
    if constexpr (sizeof...(I) == 0)
    {
        return T<C>{};
    }
    else
    {
        return (T<C, I>{} + ...);
    }
}

namespace detail
{
    // Fully generic addition operator for ring or multivector elements
    // C := common template parameter (overloaded to mean different things)
    template <template <typename, typename...> typename T, typename C, typename... I, typename... J, typename... K>
    [[nodiscard]] constexpr auto sum(T<C, I...> accum, T<C, J...> lhs, T<C, K...> rhs) noexcept
    {
        using lhs_t = decltype(lhs);
        using rhs_t = decltype(rhs);
        if constexpr (lhs.is_zero)
        {
            return T<C, I..., K...>{};
        }
        else if constexpr (rhs.is_zero)
        {
            return T<C, I..., J...>{};
        }
        else
        {
            if constexpr (equals<typename lhs_t::first_t, typename rhs_t::first_t>::value())
            {
                using next_t = decltype(typename lhs_t::first_t{} + typename rhs_t::first_t{});
                if constexpr (next_t::is_zero)
                {
                    return sum(accum, typename lhs_t::subsequent_t{}, typename rhs_t::subsequent_t{});
                }
                else
                {
                    return sum(T<C, I..., next_t>{}, typename lhs_t::subsequent_t{}, typename rhs_t::subsequent_t{});
                }
            }
            else if constexpr (less<typename lhs_t::first_t, typename rhs_t::first_t>::value())
            {
                return sum(T<C, I..., typename lhs_t::first_t>{}, typename lhs_t::subsequent_t{}, rhs);
            }
            else
            {
                return sum(T<C, I..., typename rhs_t::first_t>{}, lhs, typename rhs_t::subsequent_t{});
            }
        }
    }

    // Using the product operator defined by Algebra, distributes the term lhs into the monomials of rhs
    template <typename Algebra, template <typename, typename...> typename T, typename C, typename I, typename... J>
    [[nodiscard]] constexpr auto distribute(I lhs, T<C, J...> rhs) noexcept
    {
        return (Algebra::product(lhs, J{}) + ...);
    }

    // Fully generic multiplication operator for ring and multivector elements
    // Algebra := algebra which defines the product operator between terms
    template <typename Algebra, template <typename, typename...> typename T, typename C, typename... I, typename... J>
    [[nodiscard]] constexpr auto product(T<C, I...> lhs, T<C, J...> rhs) noexcept
    {
        if constexpr (lhs.size == 0 || rhs.size == 0)
        {
            // No generators to distribute
            return T<C>{};
        }
        else
        {
            if constexpr (Algebra::has_order_preserving_product)
            {
                // We have a choice of distributing items on the left to items on the right or vice versa. We'll choose
                // the former.
                return (distribute<Algebra>(I{}, rhs) + ...);
            }
            else
            {
                // In the event that the algebraic product is not order preserving, we need an additional pass to group
                // like elements together.
                return (simplify(distribute<Algebra>(I{}, rhs)) + ...);
            }
        }
    }
} // namespace detail

template <typename Q1, typename Q2, typename... G1, typename... G2>
struct equals<monomial<Q1, G1...>, monomial<Q2, G2...>>
{
    [[nodiscard]] constexpr static bool value() noexcept
    {
        if constexpr (sizeof...(G1) == sizeof...(G2))
        {
            return (std::is_same<G1, G2>::value && ...);
        }
        else
        {
            return false;
        }
    }
};

namespace detail
{
    template <typename T1, typename... T1s, typename T2, typename... T2s>
    [[nodiscard]] constexpr bool compare_lex(std::tuple<T1, T1s...>, std::tuple<T2, T2s...>)
    {
        if constexpr (equals<T1, T2>::value)
        {
            if constexpr (T1::degree < T2::degree)
            {
                return true;
            }
            else if constexpr (T1::degree > T2::degree)
            {
                return false;
            }
            else
            {
                if constexpr (sizeof...(T2s) == 0)
                {
                    return false;
                }
                else if constexpr (sizeof...(T1s) == 0)
                {
                    return true;
                }
                else
                {
                    return compare_lex(std::tuple<T1s...>{}, std::tuple<T2s...>{});
                }
            }
        }
        else if constexpr (less<T1, T2>::value)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
} // namespace detail

// Monomials are compared using the graded lexicographic ordering which is preserved under monomial multiplication
template <typename M1, typename M2, typename... G1, typename... G2>
struct less<monomial<M1, G1...>, monomial<M2, G2...>>
{
    [[nodiscard]] constexpr static bool value() noexcept
    {
        if constexpr (monomial<M1, G1...>::degree < monomial<M2, G2...>::degree)
        {
            return true;
        }
        else if constexpr (monomial<M1, G1...>::degree > monomial<M2, G2...>::degree)
        {
            return false;
        }
        else if constexpr (sizeof...(G1) == 0)
        {
            // Both monomials are empty so neither is stricly less than the other
            return false;
        }
        else if constexpr (sizeof...(G1) == sizeof...(G2))
        {
            // Monomial sizes are equal so we can apply lexicographic comparison with a fold
            return (less<G1, G2>::value || ...);
        }
        else
        {
            return detail::compare_lex(std::tuple<G1...>{}, std::tuple<G2...>{});
        }
    }
};

template <size_t E>
struct element
{
    constexpr static size_t value = E;
};

// A term encodes a single linearly independent term in a multivector expression.
// Operations that would contribute or subtract in a manner that cancels out the term
// will annihilate it during the final reduction.
// E := the basis element associated with this term
// As := the variadic addends that comprise the polynomial coefficient of the term.
// The basis element index is used to impose a partial ordering on terms
template <typename E, typename... As>
struct term
{
    using element_t = E;
    using first_t                 = void;
    constexpr static size_t size  = 0;
    constexpr static bool is_zero = true;
};

template <typename E, typename A, typename... As>
struct term<E, A, As...>
{
    using element_t = E;
    using first_t                         = A;
    using subsequent_t                    = term<E, As...>;
    constexpr static size_t size          = sizeof...(As) + 1;
    constexpr static size_t basis_element = E::value;
    constexpr static bool is_zero         = A::is_zero && sizeof...(As) == 0;
};

template <typename E1, typename E2, typename... A1, typename... A2>
struct equals<term<E1, A1...>, term<E2, A2...>>
{
    [[nodiscard]] constexpr static bool value() noexcept
    {
        return E1::value == E2::value;
    }
};

template <typename E1, typename E2, typename... A1, typename... A2>
struct less<term<E1, A1...>, term<E2, A2...>>
{
    [[nodiscard]] constexpr static bool value() noexcept
    {
        return E1::value < E2::value;
    }
};

template <int N, int D, typename E, typename... A>
[[nodiscard]] constexpr auto operator*(rational<N, D> lhs, term<E, A...>) noexcept
{
    if constexpr (N == 0)
    {
        return term<element<0>>{};
    }
    else if constexpr (sizeof...(A) == 0)
    {
        return term<element<0>>{};
    }
    else
    {
        return term<E, decltype(lhs * A{})...>{};
    }
}

template <typename E, typename... I, typename... J>
[[nodiscard]] constexpr auto operator+(term<E, I...> lhs, term<E, J...> rhs) noexcept
{
    return ::gal::detail::sum(term<E>{}, lhs, rhs);
}

template <typename E, typename... I>
[[nodiscard]] constexpr auto operator-(term<E, I...> in) noexcept
{
    return minus_one{} * in;
}

template <typename E, typename... I, typename... J>
[[nodiscard]] constexpr auto operator-(term<E, I...> lhs, term<E, J...> rhs) noexcept
{
    return ::gal::detail::sum(term<E>{}, lhs, -rhs);
}

// The ring algebra wraps the traditional product between monomials
template <typename E>
struct ring_algebra
{
    constexpr static bool has_order_preserving_product = true;
    template <typename M1, typename M2>
    [[nodiscard]] constexpr static auto product(M1 lhs, M2 rhs)
    {
        return term<E, decltype(lhs * rhs)>{};
    }
};

template <typename E, typename... M1, typename... M2>
[[nodiscard]] constexpr auto operator*(term<E, M1...> lhs, term<E, M2...> rhs) noexcept
{
    return ::gal::detail::product<ring_algebra<E>>(lhs, rhs);
}

// With respect to the code here, a multivector is the additive union of terms
// Vectors are tagged with an identifier used to coalesce annihilating terms during expression evaluation.
// An tag of 0 indicates the multivector is untagged (common for intermediate multivector results)
template <typename V, typename... Terms>
struct multivector
{
    using first_t                 = void;
    using subsequent_t            = void;
    constexpr static size_t size  = 0;
    constexpr static bool is_zero = true;
};

template <typename T, typename... Ts>
struct multivector<void, T, Ts...>
{
    using first_t                 = T;
    using subsequent_t            = multivector<void, Ts...>;
    constexpr static size_t size  = sizeof...(Ts) + 1;
    constexpr static bool is_zero = T::is_zero && sizeof...(Ts) == 0;
};

template <int N, int D, typename... T>
[[nodiscard]] constexpr auto operator*(rational<N, D> lhs, multivector<void, T...>) noexcept
{
    if constexpr (N == 0)
    {
        return multivector<void>{};
    }
    else if constexpr (sizeof...(T) == 0)
    {
        return multivector<void>{};
    }
    else
    {
        return multivector<void, decltype(lhs * T{})...>{};
    }
}

template <typename... I, typename... J>
[[nodiscard]] constexpr auto operator+(multivector<void, I...> lhs, multivector<void, J...> rhs) noexcept
{
    return ::gal::detail::sum(multivector<void>{}, lhs, rhs);
}

template <int N, int D, typename... I>
[[nodiscard]] constexpr auto operator+(rational<N, D>, multivector<void, I...> rhs) noexcept
{
    return multivector<void, term<element<0>, monomial<rational<N, D>>>>{} +  rhs;
}

template <int N, int D, typename... I>
[[nodiscard]] constexpr auto operator+(multivector<void, I...> lhs, rational<N, D>) noexcept
{
    return lhs + multivector<void, term<element<0>, monomial<rational<N, D>>>>{};
}

template <typename... I>
[[nodiscard]] constexpr auto operator-(multivector<void, I...> in) noexcept
{
    return minus_one{} * in;
}

template <typename... I, typename... J>
[[nodiscard]] constexpr auto operator-(multivector<void, I...> lhs, multivector<void, J...> rhs) noexcept
{
    return ::gal::detail::sum(multivector<void>{}, lhs, -rhs);
}

// A vector space is not equipped with any sort of product, however, it is useful to provide a product
// operator to use as the basis of any algebraic products to be defined with algebras later.
// Because products are not necessarily commutative in general, we do not assume commutativity here.
// If the algebra attached to the vector space below is not a graded algebra, then the multivectors are
// just normal vectors.
// Algebra := defines invokable product between basis elements (size_t, size_t) -> size_t
template <typename Algebra, typename... I, typename... J>
[[nodiscard]] constexpr auto operator*(multivector<void, I...> lhs, multivector<void, J...> rhs) noexcept
{
    return ::gal::detail::product<Algebra>(lhs, rhs);
}

template <typename E, typename... G>
[[nodiscard]] constexpr auto fiber(term<E, G...>) noexcept
{
    return term<E, monomial<one>>{};
}

// Returns the multivector with all polynomial generators stripped
template <typename... I>
[[nodiscard]] constexpr auto fiber(multivector<void, I...>) noexcept
{
    return multivector<void, decltype(fiber(I{}))...>{};
}

// An algebra over a field is a multivector space equipped with a bilinear product.
// The bilinear product is assumed to be associative but not commutative, and should take two elements and produce
// a third with a sign depending on the metric tensor element.
// Here, the `algebra` class template is not directly linked to the product as many algebras feature multiple
// bilinear products (left/right contraction, exterior, geometric, etc).
// This library assumes that the algebra is a *finite* algebra (read. not finitely generated or infinite)
// When the multivector space of the algebra is graded, the gradation is assumed to be the number of basis elements.
template <typename Field, typename Metric, bool Graded = true, bool HasScalar = true>
struct algebra
{
    using field_t      = Field;
    using metric_t     = Metric;
    using graded_t     = std::bool_constant<Graded>;
    using has_scalar_t = std::bool_constant<HasScalar>;

    [[nodiscard]] constexpr static bool is_graded() noexcept
    {
        return graded_t::value;
    }

    [[nodiscard]] constexpr static bool has_scalar() noexcept
    {
        return has_scalar_t::value;
    }

    // The "dimension" of the algebra is not the dimension of the space the algebra may be intended to model.
    // For example, the Euclidean space (3D) may be embedded in the Conformal Geometric Algebra (32D).
    [[nodiscard]] constexpr static size_t dimension() noexcept
    {
        constexpr auto metric_dimension = metric_t::dimension;
        if constexpr (is_graded())
        {
            // TODO(C++23): Change literal suffix to match size_t
            return (1ull << metric_dimension) - (has_scalar() ? 0 : 1);
        }
        else
        {
            return metric_dimension - (has_scalar() ? 0 : 1);
        }
    }

    static_assert(dimension() <= std::numeric_limits<size_t>::max(),
                  "The size of the graded ring associated with this metric tensor cannot be encoded in a std::array");
};
} // namespace gal