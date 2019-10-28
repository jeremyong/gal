#pragma once

#include <gal/cga.hpp>

// SANITY TEST
//
// When the input angles are:
//
//   ang1 = deg2rad(14.0)
//   ang2 = deg2rad(-25.0)
//   ang3 = deg2rad(32.6)
//   ang4 = deg2rad(66.9)
//   ang5 = deg2rad(-42.0)
//
// the result must be:
//
//   R1 = 0.992546 + 0.121869 * e1^e2
//   R2 = 0.976296 + 0.210006 * e1^e3 - 0.0523604 * e2^e3 + 142.804 * e1^ni - 35.6051 * e2^ni - 43.2871 * e3^ni
//   R3 = 0.959806 - 0.272314 * e1^e3 + 0.0678954 * e2^e3 - 404.827 * e1^ni + 100.935 * e2^ni + 161.69 * e3^ni
//   T2 = 1 - 182.475 * e1^ni + 45.4961 * e2^ni + 41.6926 * e3^ni
//   R4 = 0.834423 + 0.296658 * e1^e2 + 0.112228 * e1^e3 + 0.450123 * e2^e3 + 145.475 * e1^ni + 583.469 * e2^ni
//   Rg = 0.933654 + 0.277405 * e1^e2 + 0.0937376 * e1^e3 - 0.206198 * e2^e3 + 112.644 * e1^ni - 763.223 * e2^ni
//        - 174.171 * e3^ni
//   Jg_f = 1351.52 * e1 - 498.052 * e2 + 2132.49 * e3 + 0.99996 * no + 3.31122e+06 * ni
//
// where 'no' is the null point at the origin and 'ni' is the null point at infinity.

namespace gabenchmark
{
using namespace gal::cga;
using point  = point<real_t>;
using scalar = gal::scalar<cga_algebra, real_t>;
using namespace gal;

template <typename T = float>
struct point_z
{
    using algebra_t = cga_algebra;
    using value_t   = T;

    [[nodiscard]] constexpr static mv<algebra_t, 0, 3, 3> ie(uint32_t id) noexcept
    {
        // A CGA point is represented as no + p + 1/2 p^2 ni
        return {mv_size{0, 3, 3},
                {
                },
                {
                    mon{one, one, 0, 0},         // p_z
                    mon{one, zero, 0, 0},        // no
                    mon{one_half, rat{2}, 0, 0}, // 1/2 p_z^2
                },
                {
                    term{1, 0, 0b100},  // p_z
                    term{1, 1, 0b1000}, // no
                    term{1, 2, 0b10000} // 1/2 p^2 ni
                }};
    }

    [[nodiscard]] constexpr static size_t size() noexcept
    {
        return 0;
    }
};

template <typename T = float>
union point_xz
{
    using algebra_t = cga_algebra;
    using value_t   = T;

    std::array<T, 2> data;
    struct
    {
        union
        {
            T x;
            T u;
        };

        union
        {
            T z;
            T w;
        };
    };

    constexpr point_xz(T a, T c) noexcept
        : x{a}
        , z{c}
    {}

    template <uint8_t... E>
    constexpr point_xz(entity<algebra_t, T, E...> in) noexcept
        : data{in.template select<0b1, 0b100>()}
    {}

    [[nodiscard]] constexpr static mv<algebra_t, 4, 5, 4> ie(uint32_t id) noexcept
    {
        // A CGA point is represented as no + p + 1/2 p^2 ni
        return {mv_size{4, 5, 4},
                {
                    ind{id, rat{1}},     // ind0 = p_x
                    ind{id + 1, rat{1}}, // ind2 = p_z
                    ind{id, rat{2}},     // ind3 = p_x^2
                    ind{id + 1, rat{2}}, // ind5 = p_z^2
                },
                {
                    mon{one, one, 1, 0},         // p_x
                    mon{one, one, 1, 1},         // p_z
                    mon{one, zero, 0, 0},        // no
                    mon{one_half, rat{2}, 1, 1}, // 1/2 p_x^2
                    mon{one_half, rat{2}, 1, 2}, // 1/2 p_z^2
                },
                {
                    term{1, 0, 0b1},    // p_x
                    term{1, 1, 0b100},  // p_z
                    term{1, 2, 0b1000}, // no
                    term{2, 3, 0b10000} // 1/2 p^2 ni
                }};
    }

    [[nodiscard]] constexpr static size_t size() noexcept
    {
        return 2;
    }

    [[nodiscard]] constexpr T const& operator[](size_t index) const noexcept
    {
        return data[index];
    }

    [[nodiscard]] constexpr T& operator[](size_t index) noexcept
    {
        return data[index];
    }
};

// Fourth order exp expansion
template <typename T>
inline auto expp(const T& arg)
{
    auto arg2 = compute([](auto arg) { return arg * arg; }, arg);
    return compute(
        [](auto arg, auto arg2) {
            auto arg3 = arg * arg2;
            auto arg4 = arg2 * arg2;
            return frac<1> + arg + arg2 / frac<2> + arg3 / frac<6> + arg4 / frac<24>;
        },
        arg,
        arg2);
}

template <typename Scalar>
auto InverseKinematics(const Scalar& ang1, const Scalar& ang2, const Scalar& ang3, const Scalar& ang4, const Scalar& ang5)
{
    real_t d1 = 200.0, d2 = 680.0, d3 = 150.0, d4 = 140.0, d5 = 114.2;
    real_t l12 = 890.0, l23 = 880.0;

    real_t J1_x = d1;
    real_t J1_y = 0.0;
    real_t J1_z = d2;
    real_t J2_x = d1;
    real_t J2_y = 0.0;
    real_t J2_z = d2 + l12;
    real_t J3_x = d1 + l23;
    real_t J3_y = 0.0;
    real_t J3_z = d2 + l12 + d3;
    real_t Jg_x = d1 + l23 + d4 + d5;
    real_t Jg_y = 0.0;
    real_t Jg_z = d2 + l12 + d3;

    point_xz<real_t> J1{J1_x, J1_z};
    point_xz<real_t> J2{J2_x, J2_z};
    point_xz<real_t> J3{J3_x, J3_z};
    point_xz<real_t> Jg{Jg_x, Jg_z};
    point_z<real_t> Pz;

    auto Lz = compute(
        [](auto Pz, auto ang1) { return frac<1, 2> * ang1 * ((n_o<real_t> ^ Pz ^ n_i<real_t>) >> ips<real_t>); },
        Pz,
        scalar{ang1});
    auto R1 = expp(Lz);

    point P2_help{J1_x, J1_y + 1.0, J1_z};

    auto L2 = compute(
        [](auto R1, auto J1, auto P2_help, auto ang2) {
            auto L2init = (J1 ^ P2_help ^ n_i<real_t>) >> ips<real_t>;
            return ang2 * (L2init % R1) / frac<2>;
        },
        R1,
        J1,
        P2_help,
        scalar{ang2});
    auto R2 = expp(L2);

    point P3_help{J2_x, J2_y + 1.0, J2_z};

    auto R21 = compute([](auto R1, auto R2) { return R2 * R1; }, R1, R2);

    auto [J2_f, L3] = compute(
        [](auto R21, auto J2, auto P3_help, auto ang3) {
            auto L3init = (J2 ^ P3_help ^ n_i<real_t>) >> ips<real_t>;
            auto J2_f   = J2 % R21;
            return std::make_tuple(J2_f, frac<1, 2> * ang3 * (L3init % R21));
        },
        R21,
        J2,
        P3_help,
        scalar{ang3});

    auto R3 = expp(L3);

    auto [J2_rot1, t2_help] = compute(
        [](auto R1, auto J2, auto J2_f) {
            auto J2_rot1 = J2 % R1;
            auto t2      = extract<0b1, 0b10, 0b100>{}(J2_f)-extract<0b1, 0b10, 0b100>{}(J2_rot1);
            return std::make_tuple(J2_rot1, frac<-1, 2> * t2 ^ n_i<real_t>);
        },
        R1,
        J2,
        J2_f);

    auto T2 = expp(t2_help);

    auto [L4init, L4weight, R3T2R1] = compute(
        [](auto J3, auto Jg, auto R3, auto T2, auto R1) {
            auto L4init   = (J3 ^ Jg ^ n_i<real_t>) >> ips<real_t>;
            auto L4weight = L4init >> ~L4init;
            return std::make_tuple(L4init, L4weight, R3 * T2 * R1);
        },
        J3,
        Jg,
        R3,
        T2,
        R1);

    auto norm = L4weight[0] < 0 ? -std::sqrt(-L4weight[0]) : std::sqrt(L4weight[0]);
    for (auto& component : L4init)
    {
        component /= norm;
    }

    auto L4 = compute([](auto L4init, auto R3T2R1, auto ang4) { return frac<1, 2> * ang4 * (L4init % R3T2R1); },
                      L4init,
                      R3T2R1,
                      scalar{ang4});
    auto R4 = expp(L4);

    point Pg_help{J3_x, J3_y + 1.0, J3_z};
    auto [Lginit, R4R3T2R1] = compute(
        [](auto R4, auto R3T2R1, auto J3, auto Pg_help) {
            auto Lginit   = (J3 ^ Pg_help ^ n_i<real_t>) >> ips<real_t>;
            auto R4R3T2R1 = R4 * R3T2R1;
            return std::make_tuple(Lginit, R4R3T2R1);
        },
        R4,
        R3T2R1,
        J3,
        Pg_help);

    auto Lg = compute([](auto Lginit, auto R4R3T2R1, auto ang5) { return frac<1, 2> * ang5 * (Lginit % R4R3T2R1); },
                      Lginit,
                      R4R3T2R1,
                      scalar{ang5});

    auto Rg = expp(Lg);

    auto Rfinal = compute([](auto Rg, auto R4R3T2R1) { return Rg * R4R3T2R1; }, Rg, R4R3T2R1);

    auto Jg_f = compute([](auto Rfinal, auto Jg) { return Jg % Rfinal; }, Rfinal, Jg);

    return std::make_tuple(R1, R2, R3, T2, R4, Rg, Jg_f);
}
} // namespace gabenchmark