#include "test_util.hpp"

#include <gal/cga.hpp>
#include <gal/format.hpp>
#include <iostream>

#include <doctest/doctest.h>

using real_t = double;
#include "../benchmark/ga-benchmark/SpecializedAlgorithmInverseKinematics.hpp"

real_t deg_to_rad(real_t deg)
{
    return M_PI * deg / real_t{180};
}

using namespace gal::cga;
using namespace gal;

TEST_SUITE_BEGIN("inverse-kinematics");

TEST_CASE("cga-ik-rt")
{
    // using point  = point<real_t>;
    // using scalar = gal::scalar<cga_algebra, real_t>;
    //point Pz{0, 0, 1};
    // auto no = decltype(cga::e_o<real_t>)::lhs;
    // auto ni = decltype(cga::e_inf<real_t>)::lhs;
    // auto ip = decltype(cga::ips<real_t>)::lhs;
    // auto pz1 = gal::detail::product(cga_algebra::exterior{}, no, point::ie(0));
    // auto pz2 = gal::detail::product(cga_algebra::exterior{}, pz1, ni);
    // auto lz = gal::detail::product(cga_algebra::contract{}, pz2, ip);
    // auto lz2 = gal::detail::product(cga_algebra::geometric{}, lz, lz);
    // auto lz3 = gal::detail::product(cga_algebra::geometric{}, lz2, lz);
    // auto lz4 = gal::detail::product(cga_algebra::geometric{}, lz2, lz2);
    // auto ex = gal::detail::scalar_sum(rat{1}, lz2);
    // auto ex1 = gal::detail::sum(ex, lz2);
    // auto ex2 = gal::detail::sum(ex1, lz3);
    // auto ex3 = gal::detail::sum(ex2, lz4);
}

TEST_CASE("cga-ik")
{
    auto ang1{deg_to_rad(14.0)};
    auto ang2{deg_to_rad(-25.0)};
    auto ang3{deg_to_rad(32.6)};
    auto ang4{deg_to_rad(66.9)};
    auto ang5{deg_to_rad(-42.0)};
    auto&& [R1, R2, R3, T2, R4, Rg, Jg_f] = gabenchmark::InverseKinematics(ang1, ang2, ang3, ang4, ang5);
    std::cout << "CGA IK benchmark results\n"
              << "R1: " << to_string(R1) << '\n'
              << "R2: " << to_string(R2) << '\n'
              << "R3: " << to_string(R3) << '\n'
              << "T2: " << to_string(T2) << '\n'
              << "R4: " << to_string(R4) << '\n'
              << "Rg: " << to_string(Rg) << '\n'
              << "Jg_f: " << to_string(Jg_f) << std::endl;

    CHECK_EQ(R1[0], doctest::Approx(0.992546).epsilon(0.01));
    CHECK_EQ(R1[1], doctest::Approx(0.121869).epsilon(0.01));
    CHECK_EQ(T2[0], doctest::Approx(1).epsilon(0.01));
    CHECK_EQ(T2[1], doctest::Approx(-182.475).epsilon(0.01));
    CHECK_EQ(T2[2], doctest::Approx(45.4961).epsilon(0.01));
    CHECK_EQ(Rg[0], doctest::Approx(0.933654).epsilon(0.01));
    CHECK_EQ(Rg[1], doctest::Approx(0.277405).epsilon(0.01));
    CHECK_EQ(Rg[2], doctest::Approx(0.0937376).epsilon(0.01));
    CHECK_EQ(Rg[3], doctest::Approx(-0.206198).epsilon(0.01));
    CHECK_EQ(Rg[4], doctest::Approx(112.648).epsilon(0.01));
    CHECK_EQ(Rg[5], doctest::Approx(-763.223).epsilon(0.01));
    CHECK_EQ(Rg[6], doctest::Approx(-174.171).epsilon(0.01));
    CHECK_EQ(Jg_f[0], doctest::Approx(1350.82).epsilon(0.01));
    CHECK_EQ(Jg_f[1], doctest::Approx(-498.052).epsilon(0.01));
    CHECK_EQ(Jg_f[2], doctest::Approx(2132.49).epsilon(0.01));
    CHECK_EQ(Jg_f[3], doctest::Approx(0.99996).epsilon(0.01));
    CHECK_EQ(Jg_f[4], doctest::Approx(3.31122e6).epsilon(0.01));
}

TEST_SUITE_END();