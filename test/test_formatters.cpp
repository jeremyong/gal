#include "test_util.hpp"

#include <doctest/doctest.h>
#include <gal/formatters.hpp>
#include <gal/pga.hpp>
#include <string>

using namespace gal;
using namespace std::string_literals;

TEST_SUITE_BEGIN("formatters");

TEST_CASE("generator-format")
{
    SUBCASE("empty-multivector")
    {
        CHECK_EQ(fmt::format("{}", ::gal::multivector<void>{}), "0"s);
    }
}

TEST_SUITE_END();