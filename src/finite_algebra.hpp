#pragma once

#include <array>
#include <cstdint>
#include <type_traits>
#include <tuple>
#include <utility>

namespace gal
{
// Abstractions for compile-time computation of algebraic operations
// Here, we restrict ourselves to associative (but not necessarily commutative) algebras over a field

// Template parameters encode the [metric signature](https://en.wikipedia.org/wiki/Metric_signature) of
// the metric tensor underlying the algebra.
// V: # of basis elements that have negative norm
// P: # of basis elements with positive norm
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
template <size_t V, size_t P, size_t R> struct metric
{
    constexpr static size_t v = V;
    constexpr static size_t p = P;
    constexpr static size_t r = R;

    // The dimension of the metric corresponds to the total number of basis elements represented.
    // Note that the metric can induce a larger dimension in an algebra if the multivector space
    // associated with it is [graded](https://en.wikipedia.org/wiki/Graded_multivector_space) (e.g. Clifford Algebra).
    constexpr static size_t dimension = V + P + R;

    [[nodiscard]] constexpr static int element_norm(size_t e) noexcept
    {
        if (e < v)
        {
            return -1;
        }
        else if (e >= v + p)
        {
            return 0;
        }
        else
        {
            return 1;
        }
    }
};

template <int D>
struct degree
{
    constexpr static int value = D;
};

// Although the multivector space will ultimately be defined over a field, we decompose the field into the
// product of scalars (essentially factoring out a free module). The free module over the ring of
// integers has the nice property that we can condense computation by performing arithmetic exactly
// at compiole time.
// The free module we factor out is a bimodule (i.e. there is no preference for left or right
// multiplication by the scalar).
// A factor encodes the factor weight, the factor degree, as well as its identifier (if available).
// NOTE: "Degree" here is meant in the sense of a polynimal (e.g. x^2 has degree 2).
// We permit negative degrees to allow expressing linear combinations of nth-roots as well.
template <typename Degree = degree<1>, size_t ID = 0> struct factor
{
    // The precedence for the partial ordering of the factor is input_id < factor_id < degree
    constexpr static int degree = Degree::value;
    constexpr static int id = ID;
};

template <typename T1, typename T2>
struct equals
{
    constexpr static bool value = false;
};

// The factors that make up a monomial are weakly ordered based on the source identifiers.
// If all factors are identified (ID != 0), the ordering becomes a total order.
// NOTE: the comparison operators here do not consider the actual quantities held by the factors.
template <typename D1, size_t ID1, typename D2, size_t ID2>
struct equals<factor<D1, ID1>, factor<D2, ID2>>
{
    constexpr static bool value = ID1 != 0 && ID2 != 0 && ID1 == ID2 && D1::value == D2::value;
};

template <typename T1, typename T2>
struct less
{
    constexpr static bool value = false;
};

template <typename D1, size_t ID1, typename D2, size_t ID2>
struct less<factor<D1, ID1>, factor<D2, ID2>>
{
    // NOTE: if the factors are untagged (id == 0), it is irrelevant which factor gets ordered first.
    constexpr static bool value = ID1 < ID2 || D1::value < D2::value;
};

template <int M>
struct multiplier
{
    constexpr static int value = M;
};

// The coefficient of a term is a linear combination of factors which are kept separate for the purposes of determining
// an optimal computational strategy. An `monomial` is a product of `factor`s.
// Invariant: the factors are sorted in order
template <typename M, typename... Factors> struct monomial
{
    constexpr static size_t size = sizeof...(Factors) + 1;
    constexpr static int multiplier = M::value;
    constexpr static int degree = (Factors::degree + ...);
};

template <typename... Factors> struct monomial<multiplier<0>, Factors...>
{
    constexpr static size_t size = 0;
    constexpr static int multiplier = 0;
};

template <typename D1, typename D2, size_t ID1, size_t ID2>
[[nodiscard]] constexpr auto operator*(factor<D1, ID1>, factor<D2, ID2>)
{
    using lhs_t = factor<D1, ID1>;
    using rhs_t = factor<D2, ID2>;
    if constexpr (ID1 == ID2)
    {
        return monomial<multiplier<1>, factor<degree<D1::value + D2::value>, ID1>>{};
    }
    else if constexpr (less<lhs_t, rhs_t>::value)
    {
        return monomial<multiplier<1>, lhs_t, rhs_t>{};
    }
    else
    {
        return monomial<multiplier<1>, rhs_t, lhs_t>{};
    }
}

namespace impl
{
template <typename M, typename... I, typename F1, typename... J, typename F2, typename... K>
[[nodiscard]] constexpr auto
product(monomial<M, I...> accum, monomial<multiplier<1>, F1, J...> lhs, monomial<multiplier<1>, F2, K...> rhs)
{
    if constexpr (F1::id == 0 && F2::id == 0)
    {
        // reduce an extra recursion if both factors are untagged
        if constexpr (sizeof...(J) == 0)
        {
            return monomial<M, I..., F1, F2, K...>{};
        }
        else if constexpr (sizeof...(K) == 0)
        {
            return monomial<M, I..., F1, F2, J...>{};
        }
        else
        {
            return product(monomial<M, I..., F1, F2>{}, monomial<multiplier<1>, J...>{}, monomial<multiplier<1>, K...>{});
        }
    }
    if constexpr (F1::id == 0 || F1::id < F2::id)
    {
        if constexpr (sizeof...(J) == 0)
        {
            return monomial<M, I..., F1, F2, K...>{};
        }
        else
        {
            return product(monomial<M, I..., F1>{}, monomial<multiplier<1>, J...>{}, rhs);
        }
    }
    else if constexpr (F2::id == 0 || F2::id < F1::id)
    {
        if constexpr (sizeof...(K) == 0)
        {
            return monomial<M, I..., F2, F1, J...>{};
        }
        else
        {
            return product(monomial<M, I..., F2>{}, lhs, monomial<multiplier<1>, K...>{});
        }
    }
    else
    {
        // Both F1 and F2 are tagged and represent the same factor
        using next_factor_t = factor<degree<F1::degree + F2::degree>, F1::id>;
        if constexpr (sizeof...(J) == 0)
        {
            return monomial<M, I..., next_factor_t, K...>{};
        }
        else if constexpr (sizeof...(K) == 0)
        {
            return monomial<M, I..., next_factor_t, J...>{};
        }
        else
        {
            return product(monomial<M, I..., next_factor_t>{}, monomial<multiplier<1>, J...>{}, monomial<multiplier<1>, K...>{});
        }
    }
}
}

// Order preserving monomial multiplication
template <typename M1, typename M2, typename... F1, typename... F2>
[[nodiscard]] constexpr auto operator*(monomial<M1, F1...>, monomial<M2, F2...>)
{
    using m = multiplier<M1::value * M2::value>;
    if constexpr (M1::value == 0 || M2::value == 0)
    {
        return monomial<multiplier<0>>{};
    }
    else if constexpr (sizeof...(F1) == 0)
    {
        return monomial<m, F2...>{};
    }
    else if constexpr (sizeof...(F2) == 0)
    {
        return monomial<m, F1...>{};
    }
    else
    {
        return ::gal::impl::product(monomial<m>{}, monomial<multiplier<1>, F1...>{}, monomial<multiplier<1>, F2...>{});
    }
}

template <typename M, typename...F>
[[nodiscard]] constexpr auto operator-(monomial<M, F...>)
{
    return monomial<multiplier<-M::value>, F...>{};
}

template <typename M1, typename M2, typename... F>
[[nodiscard]] constexpr auto operator+(monomial<M1, F...> lhs, monomial<M2, F...> rhs) noexcept
{
    if constexpr (M1::value == 0)
    {
        return rhs;
    }
    else if constexpr (M2::value == 0)
    {
        return lhs;
    }
    {
        return monomial<multiplier<M1::value + M2::value>, F...>{};
    }
};

template <typename M1, typename M2, typename... F>
[[nodiscard]] constexpr auto operator-(monomial<M1, F...>, monomial<M2, F...>) noexcept
{
    if constexpr (M1::value == M2::value)
    {
        return monomial<multiplier<0>>{};
    }
    else
    {
        return monomial<multiplier<M1::value - M2::value>, F...>{};
    }
};

namespace impl
{
// Fully generic addition operator for ring or multivector elements
// C := common template parameter (overloaded to mean different things)
template <template <typename, typename...> typename T, typename C, typename...I, typename...J, typename... K>
[[nodiscard]] constexpr auto
sum(T<C, I...> accum, T<C, J...> lhs, T<C, K...> rhs) noexcept
{
    using lhs_t = decltype(lhs);
    using rhs_t = decltype(rhs);
    if constexpr (lhs.size == 0)
    {
        return T<C, I..., K...>{};
    }
    else if constexpr (rhs.size == 0)
    {
        return T<C, I..., J...>{};
    }
    else
    {
        if constexpr (equals<typename lhs_t::first_t, typename rhs_t::first_t>::value())
        {
            using next_t = decltype(typename lhs_t::first_t{} + typename rhs_t::first_t{});
            if constexpr (next_t::size == 0)
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
[[nodiscard]] constexpr auto
distribute(I lhs, T<C, J...> rhs) noexcept
{
    return (Algebra::product(lhs, J{}) + ...);
}

// Unlike other functions in this module, `simplify` does not assume its input is ordered with respect to
// any monomial ordering. It relies on the addition operator to coalesce like-terms and the output is
// respect to be well-ordered.
template <template <typename, typename...> typename T, typename C, typename... I>
[[nodiscard]] constexpr auto
simplify(T<C, I...>) noexcept
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

// Fully generic multiplication operator for ring and multivector elements
// Algebra := algebra which defines the product operator between terms
template <typename Algebra, template <typename, typename...> typename T, typename C, typename...I, typename... J>
[[nodiscard]] constexpr auto
product(T<C, I...> lhs, T<C, J...> rhs) noexcept
{
    if constexpr (lhs.size == 0 || rhs.size == 0)
    {
        // No factors to distribute
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
}

// Lack of constexpr arguments prevent this from being implementable as an operator
template <int S, typename T>
struct scale
{
};

template <int S, typename M, typename... F>
struct scale<S, monomial<M, F...>>
{
    using type = monomial<multiplier<S * M::value>, F...>;
};

template <typename M1, typename M2, typename... F1, typename... F2>
struct equals<monomial<M1, F1...>, monomial<M2, F2...>>
{
    [[nodiscard]] constexpr static bool value() noexcept
    {
        if constexpr (sizeof...(F1) == sizeof...(F2))
        {
            return (std::is_same<F1, F2>::value && ...);
        }
        else
        {
            return false;
        }
    }
};

namespace impl
{
template <typename T1, typename... T1s, typename T2, typename... T2s>
[[nodiscard]] constexpr bool compare_lex(std::tuple<T1, T1s...>, std::tuple<T2, T2s...>)
{
    if constexpr (T1::id < T2::id)
    {
        return true;
    }
    else if (T1::id > T2::id)
    {
        return false;
    }
    else if (T1::id == T2::id && T1::degree < T2::degree)
    {
        return true;
    }
    else if (T1::id == T2::id && T1::degree > T2::degree)
    {
        return false;
    }
    else if constexpr (sizeof...(T2s) == 0)
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
} // namespace impl

// Monomials are compared using the graded lexicographic ordering which is preserved under monomial multiplication
template <typename M1, typename M2, typename... F1, typename... F2>
struct less<monomial<M1, F1...>, monomial<M2, F2...>>
{
    [[nodiscard]] constexpr static bool value() noexcept
    {
        if constexpr (monomial<M1, F1...>::degree < monomial<M2, F2...>::degree)
        {
            return true;
        }
        else if constexpr (monomial<M1, F1...>::degree > monomial<M2, F2...>::degree)
        {
            return false;
        }
        else if constexpr (sizeof...(F1) == 0)
        {
            // Both monomials are empty so neither is stricly less than the other
            return false;
        }
        else if constexpr (sizeof...(F1) == sizeof...(F2))
        {
            // Monomial sizes are equal so we can apply lexicographic comparison with a fold
            return ((F1::id < F2::id) || ...);
        }
        else
        {
            return impl::compare_lex(std::tuple<F1...>{}, std::tuple<F2...>{});
        }
    }
};

template <size_t E>
struct element
{
    constexpr static size_t value = E;
};

// A term encodes all the originating factors that multiplicatively comprise the term.
// Operations that would contribute or subtract in a manner that cancels out the term
// will annihilate it during the final reduction.
// E := the basis element associated with this term
// The basis element index is used to impose a partial ordering on terms
template <typename E, typename... As> struct term
{
    using first_t = void;
    constexpr static size_t size = 0;
};

template <typename E, typename A, typename... As>
struct term<E, A, As...>
{
    using first_t = A;
    using subsequent_t = term<E, As...>;
    constexpr static size_t size = sizeof...(As) + 1;
    constexpr static size_t basis_element = E::value;
};

template <int S, typename E, typename... A>
struct scale<S, term<E, A...>>
{
    // Scale the monomials of the term by the lhs
    using type = term<E, typename scale<S, A>::type...>;
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

template <typename E, typename... I, typename... J>
[[nodiscard]] constexpr auto operator+(term<E, I...> lhs, term<E, J...> rhs) noexcept
{
    return ::gal::impl::sum(term<E>{}, lhs, rhs);
}

template <typename E, typename... I>
[[nodiscard]] constexpr auto operator-(term<E, I...>) noexcept
{
    return typename scale<-1, term<E, I...>>::type{};
}

template <typename E, typename... I, typename... J>
[[nodiscard]] constexpr auto operator-(term<E, I...> lhs, term<E, J...> rhs) noexcept
{
    return ::gal::impl::sum(term<E>{}, lhs, -rhs);
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
    return ::gal::impl::product<ring_algebra<E>>(lhs, rhs);
}

// With respect to the code here, a multivector is the additive union of terms
// Vectors are tagged with an identifier used to coalesce annihilating terms during expression evaluation.
// An tag of 0 indicates the multivector is untagged (common for intermediate multivector results)
template <typename V, typename... Terms> struct multivector
{
    using first_t = void;
    using subsequent_t = void;
    constexpr static size_t size = 0;
};

template <typename T, typename... Ts> struct multivector<void, T, Ts...>
{
    using first_t = T;
    using subsequent_t = multivector<void, Ts...>;
    constexpr static size_t size = sizeof...(Ts) + 1;
};

template <int S, typename... T>
struct scale<S, multivector<void, T...>>
{
    using type = multivector<void, typename scale<S, T>::type...>;
};

template <typename... I, typename... J>
[[nodiscard]] constexpr auto operator+(multivector<void, I...> lhs, multivector<void, J...> rhs)
{
    return ::gal::impl::sum(multivector<void>{}, lhs, rhs);
}

template <typename... I>
[[nodiscard]] constexpr auto operator-(multivector<void, I...>)
{
    return typename scale<-1, multivector<void, I...>>::type{};
}

template <typename... I, typename... J>
[[nodiscard]] constexpr auto operator-(multivector<void, I...> lhs, multivector<void, J...> rhs)
{
    return ::gal::impl::sum(multivector<void>{}, lhs, -rhs);
}

// A vector space is not equipped with any sort of product, however, it is useful to provide a product
// operator to use as the basis of any algebraic products to be defined with algebras later.
// Because products are not necessarily commutative in general, we do not assume commutativity here.
// If the algebra attached to the vector space below is not a graded algebra, then the multivectors are
// just normal vectors.
// Algebra := defines invokable product between basis elements (size_t, size_t) -> size_t
template <typename Algebra, typename...I, typename... J>
[[nodiscard]] constexpr auto operator*(multivector<void, I...> lhs, multivector<void, J...> rhs)
{
    return ::gal::impl::product<Algebra>(lhs, rhs);
}

// An algebra over a field is a multivector space equipped with a bilinear product.
// The bilinear product is assumed to be associative but not commutative, and should take two elements and produce
// a third with a sign depending on the metric tensor element.
// Here, the `algebra` class template is not directly linked to the product as many algebras feature multiple
// bilinear products (left/right contraction, exterior, geometric, etc).
// This library assumes that the algebra is a *finite* algebra (read. not finitely generated or infinite)
// When the multivector space of the algebra is graded, the gradation is assumed to be the number of basis elements.
template <typename Field, typename Metric, bool Graded = true, bool HasScalar = true> struct algebra
{
    using field_t = Field;
    using metric_t = Metric;
    using graded_t = std::bool_constant<Graded>;
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