#include "test_util.hpp"

#include <doctest/doctest.h>
#include <gal/cga.hpp>

using namespace gal;
using namespace gal::cga;

// Fourth order exp expansion
template <typename T>
inline auto expp(const T& arg)
{
    auto arg2 = compute([](auto arg) { return arg * arg; }, arg);
    return compute(
        [](auto arg, auto arg2) {
            auto arg3 = arg * arg2;
            auto arg4 = arg2 * arg2;
            return 1 + arg + arg2 / 2 + arg3 / 6 + arg4 / 24;
        },
        arg,
        arg2);
}

TEST_SUITE_BEGIN("conformal-geometric-algebra");

TEST_CASE("null-basis-conversion")
{
    SUBCASE("to-natural-basis")
    {
        auto p  = point<>::ie(0);
        auto pn = gal::detail::to_natural_basis(p);
        CHECK_EQ(gal::detail::uses_null_basis<gal::cga::cga_algebra>, true);
        // TODO: verify
    }
}

TEST_CASE("point-norm")
{
    point<float> p{3.9f, 1.2f, -29.f};
    auto rpn = evaluate<point<float>>::rpnf([](auto p) { return p >> p; });
    std::printf("point norm: %s\n", gal::to_string(rpn).c_str());
    auto d      = evaluate<point<float>>::ie([](auto p) { return p >> p; });
    auto p_norm = compute([](auto p) { return p >> p; }, p);
    // Notice that we can compute that the norm of a point is exactly zero at compile time.
    // static_assert(p_norm.size() == 0);
    CHECK_EQ(p_norm.size(), 0);
}

TEST_CASE("line-construction")
{
    scalar<cga_algebra, float> ang{5.5};
    auto rpn = evaluate<scalar<cga_algebra, float>>::rpnf_reshaped(
        [](auto ang) { return ((1_no ^ (1_no + 1_e3 + 1_ni / 2) ^ 1_ni) >> 1_ips) * ang / 2; });
    std::printf("line construction: %s\n", gal::to_string(rpn).c_str());
    auto ie = evaluate<scalar<cga_algebra, float>>::ie_reshaped(
        [](auto ang) { return ((1_no ^ (1_no + 1_e3 + 1_ni / 2) ^ 1_ni) >> 1_ips) * ang / 2; });
    auto L1 = compute(
        [](auto ang) { return ((1_no ^ (1_no + 1_e3 + 1_ni / 2) ^ 1_ni) >> 1_ips) * ang / 2; }, ang);
    auto R1 = expp(L1);
    point<> P2_help{200.0f, 1.0, 680.0f};
    point<> J1{200.0f, 0.0f, 680.0f};
    scalar<cga_algebra, float> ang2{2.1f};
    auto L2 = compute(
        [](auto R1, auto J1, auto P2_help, auto ang2) {
            auto L2init = (J1 ^ P2_help ^ 1_ni) >> 1_ips;
            return (L2init % R1) * ang2 / 2;
        },
        R1,
        J1,
        P2_help,
        scalar{ang2});
    auto R2 = expp(L2);
}

TEST_SUITE_END();
