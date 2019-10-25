#include "test_util.hpp"

#include <doctest/doctest.h>
#include <gal/ega.hpp>
#include <gal/engine.hpp>
#include <gal/format.hpp>

#include <cstdio>

using namespace gal::ega;
using gal::compute;

constexpr inline ega_algebra::geometric gp{};

TEST_SUITE_BEGIN("euclidean-geometric-algebra");

TEST_CASE("rotors")
{
    SUBCASE("vector-reflection-through-vector")
    {
        vector<float> v1{1, 1, 0};
        vector<float> v2{1, 0, 0};

        auto ie1 = vector<float>::ie(0);
        auto ie2 = vector<float>::ie(3);
        auto result = gal::detail::product(ega_algebra::geometric{}, ie1, ie2);
        auto reverse = gal::detail::reverse(ie1);

        auto reflect = compute([](auto v1, auto v2) { return v1 % v2; }, v1, v2);
        CHECK_EQ(reflect[0], doctest::Approx(1));
        CHECK_EQ(reflect[1], doctest::Approx(-1));
        CHECK_EQ(reflect[2], doctest::Approx(0));
    }

    SUBCASE("rotate-component")
    {
        using namespace gal;

        // Simplified single coordinate rotation about a single axis

        // Rotor oriented in the +z direction
        mv<void, 2, 2, 2> r{mv_size{2, 2, 2},
                            {
                                ind{0, 1}, // cos(t/2)
                                ind{1, 1}, // sin(t/2)
                            },
                            {
                                mon{one, 1, 0, 1}, // cos(t/2)
                                mon{one, 1, 1, 1}, // sin(t/2)
                            },
                            {
                                term{1, 0, 0},   // cos(t/2)
                                term{1, 1, 0b11} // sin(t/2) * e0e1
                            }};

        // Vector oriented in the +x direction
        mv<void, 0, 1, 1> v{mv_size{0, 1, 1}, {}, {mon{one, 0, 0, 1}}, {term{1, 0, 0b1}}};

        auto rr   = gal::detail::reverse(r);
        auto vrr  = gal::detail::product(gp, v, rr);
        auto rvrr = gal::detail::product(gp, r, vrr);

        // Expected: (cos^2(t/2) - sin^2(t/2))e0 - 2sin(t/2)cos(t/2)e1

        CHECK_EQ(rvrr.size.ind, 4);
        CHECK_EQ(rvrr.size.mon, 3);
        CHECK_EQ(rvrr.size.term, 2);
        CHECK_EQ(rvrr.inds[0].id, 0);
        CHECK_EQ(rvrr.inds[0].degree, 2);
        CHECK_EQ(rvrr.inds[1].id, 1);
        CHECK_EQ(rvrr.inds[1].degree, 2);
        CHECK_EQ(rvrr.inds[2].id, 0);
        CHECK_EQ(rvrr.inds[2].degree, 1);
        CHECK_EQ(rvrr.inds[3].id, 1);
        CHECK_EQ(rvrr.inds[3].degree, 1);
        CHECK_EQ(rvrr.mons[0].q.num, 1);
        CHECK_EQ(rvrr.mons[0].count, 1);
        CHECK_EQ(rvrr.mons[0].ind_offset, 0);
        CHECK_EQ(rvrr.mons[0].degree, 2);
        CHECK_EQ(rvrr.mons[1].q.num, -1);
        CHECK_EQ(rvrr.mons[1].count, 1);
        CHECK_EQ(rvrr.mons[1].ind_offset, 1);
        CHECK_EQ(rvrr.mons[1].degree, 2);
        CHECK_EQ(rvrr.mons[2].q.num, -2);
        CHECK_EQ(rvrr.mons[2].count, 2);
        CHECK_EQ(rvrr.mons[2].ind_offset, 2);
        CHECK_EQ(rvrr.mons[2].degree, 2);
        CHECK_EQ(rvrr.terms[0].count, 2);
        CHECK_EQ(rvrr.terms[0].mon_offset, 0);
        CHECK_EQ(rvrr.terms[0].element, 1);
        CHECK_EQ(rvrr.terms[1].count, 1);
        CHECK_EQ(rvrr.terms[1].mon_offset, 2);
        CHECK_EQ(rvrr.terms[1].element, 0b10);
    }

    SUBCASE("rotate-vector")
    {
        rotor<float> r1{M_PI / 2, 0, 0, 1};
        vector<float> v1{1, 0, 0};

        vector<float> rotated = compute([](auto r1, auto v1) { return v1 % r1; }, r1, v1);
        CHECK_EQ(rotated.x, doctest::Approx(0.0));
        CHECK_EQ(rotated.y, doctest::Approx(1.0));
        CHECK_EQ(rotated.z, doctest::Approx(0.0));
    }

    SUBCASE("rotor-composition")
    {
        // Use 2 pi/4 rotations to generate the same pi/2 rotation from above
        rotor<float> r1{M_PI / 4, 0, 0, 1};
        rotor<float> r2{M_PI / 4, 0, 0, 1};
        vector<float> v1{1, 0, 0};

        vector<float> rotated = compute(
            [](auto r1, auto r2, auto v1) {
                auto r = r1 * r2;
                return v1 % r;
            },
            r1,
            r2,
            v1);
        auto r = gal::evaluate<rotor<float>, rotor<float>, vector<float>>{}.debug([](auto r1, auto r2, auto v1) {
            auto r = r1 * r2;
            return v1 % r;
        });
        CHECK_EQ(rotated.x, doctest::Approx(0.0));
        CHECK_EQ(rotated.y, doctest::Approx(1.0));
        CHECK_EQ(rotated.z, doctest::Approx(0.0));
    }
}

TEST_SUITE_END();