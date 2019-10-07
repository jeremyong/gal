#include "test_util.hpp"

#include <doctest/doctest.h>
#include <formatters.hpp>

using real_t = double;
#include "../benchmark/ga-benchmark/SpecializedAlgorithmInverseKinematics.hpp"

real_t deg_to_rad(real_t deg)
{
    return M_PI * deg / real_t{180};
}

TEST_SUITE_BEGIN("inverse-kinematics");

TEST_CASE("cga-ik")
{
    using scalar = gal::scalar<real_t>;
    
    scalar ang1{deg_to_rad(14.0)};
    scalar ang2{deg_to_rad(-25.0)};
    scalar ang3{deg_to_rad(32.6)};
    scalar ang4{deg_to_rad(66.9)};
    scalar ang5{deg_to_rad(-42.0)};
    gabenchmark::InverseKinematics(ang1, ang2, ang3, ang4, ang5);
}

TEST_SUITE_END();