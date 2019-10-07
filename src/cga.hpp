#pragma once

#include "ga.hpp"

namespace gal
{
namespace cga
{
    // The metric is defined here as the standard Minkowski spacetime. To extract the conformal representations,
    // a change of basis is required where o = 1/2 * (e3 + e4) and inf = e4 - e3.
    // The element e3 here is the added unit norm basis element and the elements e0, e1, and e2 correspond to
    // the canonical Euclidean R3 basis representation.
    using metric = ::gal::metric<4, 1, 0>;

    // The CGA is a graded algebra with 32 basis elements
    using algebra = ga::algebra<metric>;

    using e_o = multivector<void, term<element<0b100>, monomial<one_half>>, term<element<0b1000>, monomial<one_half>>>;
    using e_inf = multivector<void, term<element<0b100>, monomial<minus_one_half>>, term<element<0b1000>, monomial<one_half>>>;

    GAL_OPERATORS(algebra);

    template <int X, int Y, int Z>
    using point_t = multivector<void,
        term<element<0b1>, monomial<rational<X>>>,
        term<element<0b10>, monomial<rational<Y>>>,
        term<element<0b100>, monomial<rational<Z>>>,
        term<element<0b1000>, monomial<one_half>, monomial<rational<-(X * X + Y * Y + Z * Z), 2>>>,
        term<element<0b10000>, monomial<one_half>, monomial<rational<X * X + Y * Y + Z * Z, 2>>>>;

    template <typename T = float>
    struct alignas(16) point
    {
        using value_t = T;
        constexpr static size_t size = 3;

        template <typename Q, size_t ID>
        using p2
            = monomial<Q, generator<tag<ID, 0>, degree<2>>, generator<tag<ID, 1>, degree<2>>, generator<tag<ID, 2>, degree<2>>>;

        template <size_t ID>
        using type = multivector<void,
                                 term<element<0b1>, monomial<one, generator<tag<ID, 0>>>>,
                                 term<element<0b10>, monomial<one, generator<tag<ID, 1>>>>,
                                 term<element<0b100>, monomial<one, generator<tag<ID, 2>>>>,
                                 term<element<0b1000>,
                                      monomial<one_half>,
                                      monomial<rational<-1, 2>, generator<tag<ID, 0>, degree<2>>>,
                                      monomial<rational<-1, 2>, generator<tag<ID, 1>, degree<2>>>,
                                      monomial<rational<-1, 2>, generator<tag<ID, 2>, degree<2>>>>,
                                 term<element<0b10000>,
                                      monomial<one_half>,
                                      monomial<one_half, generator<tag<ID, 0>, degree<2>>>,
                                      monomial<one_half, generator<tag<ID, 1>, degree<2>>>,
                                      monomial<one_half, generator<tag<ID, 2>, degree<2>>>>>;

        union
        {
            T x;
            T u;
        };

        union
        {
            T y;
            T v;
        };

        union
        {
            T z;
            T w;
        };

        GAL_ACCESSORS

        template <typename Engine, typename... I>
        [[nodiscard]] constexpr static point<T> convert(const Engine& engine, multivector<void, I...> mv) noexcept
        {
            auto x_e = extract<0b1>(mv);
            auto y_e = extract<0b10>(mv);
            auto z_e = extract<0b100>(mv);

            auto&& [x, y, z] = engine.template evaluate_terms<T>(x_e, y_e, z_e);
            return {x, y, z};
        }
    };

    // TODO: provide representations for planes, spheres, flats, etc.
} // namespace cga
} // namespace gal