#pragma once

#include <cga.hpp>
#include <engine.hpp>

using gal::rational;
using gal::multivector;
using gal::extract;

namespace gabenchmark
{
template <typename... T>
inline auto expp(multivector<void, T...> arg)
{
    // return rational<1>{} + arg
           // + arg * ((rational<1, 2>{} * arg) + arg * ((rational<1, 6>{} * arg) + arg * (rational<1, 24>{} * arg)));
    // return rational<1>{} + arg
           // + arg * ((rational<1, 2>{} * arg) + arg * ((rational<1, 6>{} * arg)));
    // TODO reduce *exactness* in the rational polynomial expansion to kill terms to approach epsilon
    return rational<1>{} + arg;
}

struct foo
{
    real_t a;
    real_t b;

    template <typename Engine, typename... I>
    [[nodiscard]] constexpr static foo convert(const Engine& engine, multivector<void, I...> mv) noexcept
    {
        auto a_e = extract<0>(mv);
        auto b_e = extract<0b1100>(mv);
        auto&& [a, b] = engine.template evaluate_terms<real_t>(a_e, b_e);
        return {a, b};
    }
};

using point = gal::cga::point<real_t>;
using gal::cga::e_o;
using gal::cga::e_inf;
template <typename Scalar>
constexpr auto
InverseKinematics(const Scalar& ang1, const Scalar& ang2, const Scalar& ang3, const Scalar& ang4, const Scalar& ang5)
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

    point J1{J1_x, J1_y, J1_z};
    point J2{J2_x, J2_y, J2_z};
    point J3{J3_x, J3_y, J3_z};
    point Jg{Jg_x, Jg_y, Jg_z};

    point Pz{0.0, 0.0, 1.0};

    point P2_help{J1_x, J1_y + 1.0, J1_z};
    point P3_help{J2_x, J2_y + 1.0, J2_z};
    point Pg_help{J3_x, J3_y + 1.0, J3_z};

    gal::engine engine{J1, J2, J3, Jg, Pz, P2_help, P3_help, Pg_help, ang1, ang2, ang3, ang4, ang5};

    // TODO variadic return
    auto f = engine.template compute<foo>([](auto J1,
                                             auto J2,
                                             auto J3,
                                             auto Jg,
                                             auto Pz,
                                             auto P2_help,
                                             auto P3_help,
                                             auto Pg_help,
                                             auto ang1,
                                             auto ang2,
                                             auto ang3,
                                             auto ang4,
                                             auto ang5) {
        using namespace gal::cga;
        auto Lz = ~(e_o{} ^ Pz ^ e_inf{});
        auto R1 = expp(Lz * (rational<1, 2>{} * ang1));
        return R1;
    });
    printf("ang1: %f\n", ang1);
    printf("%f + %f * e1^e2\n", f.a, f.b);
}
} // namespace gabenchmark