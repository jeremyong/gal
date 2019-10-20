#include <gal/algorithm.hpp>
#include <doctest/doctest.h>

using gal::detail::sort;

TEST_SUITE_BEGIN("algorithm");

TEST_CASE("sort")
{
    SUBCASE("single-element")
    {
        std::array<int, 1> a = {1};
        sort(a.begin(), a.end());

        CHECK_EQ(a[0], 1);
    }

    SUBCASE("two-element")
    {
        std::array<int, 2> a = {1, 2};
        sort(a.begin(), a.end());

        CHECK_EQ(a[0], 1);
        CHECK_EQ(a[1], 2);
    }

    SUBCASE("two-element-2")
    {
        std::array<int, 2> a = {2, 1};
        sort(a.begin(), a.end());

        CHECK_EQ(a[0], 1);
        CHECK_EQ(a[1], 2);
    }

    SUBCASE("three-element")
    {
        std::array<int, 3> a = {1, 2, 3};
        sort(a.begin(), a.end());

        CHECK_EQ(a[0], 1);
        CHECK_EQ(a[1], 2);
        CHECK_EQ(a[2], 3);
    }

    SUBCASE("three-element-2")
    {
        std::array<int, 3> a = {2, 1, 3};
        sort(a.begin(), a.end());

        CHECK_EQ(a[0], 1);
        CHECK_EQ(a[1], 2);
        CHECK_EQ(a[2], 3);
    }

    SUBCASE("three-element-3")
    {
        std::array<int, 3> a = {1, 3, 2};
        sort(a.begin(), a.end());

        CHECK_EQ(a[0], 1);
        CHECK_EQ(a[1], 2);
        CHECK_EQ(a[2], 3);
    }

    SUBCASE("three-element-4")
    {
        std::array<int, 3> a = {3, 1, 2};
        sort(a.begin(), a.end());

        CHECK_EQ(a[0], 1);
        CHECK_EQ(a[1], 2);
        CHECK_EQ(a[2], 3);
    }

    SUBCASE("three-element-5")
    {
        std::array<int, 3> a = {3, 2, 1};
        sort(a.begin(), a.end());

        CHECK_EQ(a[0], 1);
        CHECK_EQ(a[1], 2);
        CHECK_EQ(a[2], 3);
    }

    SUBCASE("three-element-6")
    {
        std::array<int, 3> a = {2, 3, 1};
        sort(a.begin(), a.end());

        CHECK_EQ(a[0], 1);
        CHECK_EQ(a[1], 2);
        CHECK_EQ(a[2], 3);
    }

    SUBCASE("four-element")
    {
        std::array<int, 4> a = {2, 4, 3, 1};
        sort(a.begin(), a.end());

        CHECK_EQ(a[0], 1);
        CHECK_EQ(a[1], 2);
        CHECK_EQ(a[2], 3);
        CHECK_EQ(a[3], 4);
    }

    SUBCASE("longer-sequence")
    {
        std::array<int, 9> a = {9, 2, 6, 4, 8, 5, 3, 1, 7};
        sort(a.begin(), a.end());

        CHECK_EQ(a[0], 1);
        CHECK_EQ(a[1], 2);
        CHECK_EQ(a[2], 3);
        CHECK_EQ(a[3], 4);
        CHECK_EQ(a[4], 5);
        CHECK_EQ(a[5], 6);
        CHECK_EQ(a[6], 7);
        CHECK_EQ(a[7], 8);
        CHECK_EQ(a[8], 9);
    }
}

TEST_SUITE_END();