#include <array>
#include <doctest/doctest.h>
#include <gal/algorithm.hpp>
#include <gal/crc.hpp>
#include <gal/tuple.hpp>

using gal::detail::sort;

template <typename T>
struct add_n
{
    int incr;
    constexpr auto operator()(T in)
    {
        return in + incr;
    }
};

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

TEST_CASE("crc32")
{
    static_assert(gal::detail::crc32(0xDEADBEEF) == 0x1A5A601F);
}

TEST_CASE("tuple")
{
    using gal::tuple;
    tuple<int, char> t{1039, 'b'};
    auto i = t.template get<0>();
    CHECK_EQ(i, 1039);
    auto c = t.template get<1>();
    CHECK_EQ(c, 'b');

    SUBCASE("tuple-find-type")
    {
        auto index = t.find_t<char>();
        CHECK_EQ(index, 1);
    }

    SUBCASE("tuple-split")
    {
        auto sp = t.template split_at_type<char>();
        CHECK_EQ(decltype(sp.first)::size(), 1);
        CHECK_EQ(decltype(sp.second)::size(), 0);
        CHECK_EQ(sp.first.template get<0>(), 1039);
    }

    SUBCASE("tuple-push")
    {
        auto p = t.push(3.5f);
        CHECK_EQ(decltype(p)::size(), 3);
        CHECK_EQ(p.template get<0>(), 3.5f);
        CHECK_EQ(p.template get<1>(), 1039);
        CHECK_EQ(p.template get<2>(), 'b');
    }

    SUBCASE("tuple-set")
    {
        t.template set<0>(42);
        CHECK_EQ(t.template get<0>(), 42);
    }

    SUBCASE("tuple-map")
    {
        t.template mutate<add_n>(2);
        CHECK_EQ(t.template get<0>(), 1041);
        CHECK_EQ(t.template get<1>(), 'd');
    }

    SUBCASE("tuple-pop")
    {
        auto pop = t.pop();
        CHECK_EQ(decltype(pop.second)::size(), 1);
        CHECK_EQ(pop.first, 1039);
        CHECK_EQ(pop.second.template get<0>(), 'b');
    }
}
TEST_SUITE_END();
