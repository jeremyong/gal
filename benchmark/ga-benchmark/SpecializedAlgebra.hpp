#pragma once

#define GA_DEFAULT_FLOATING_POINT_TYPE gabenchmark::real_t

#if GABENCHMARK_CHECK_MODEL(ConformalModel)
#if GABENCHMARK_D_DIMENSIONS == 2
#include <gal/cga2.hpp>
namespace gabenchmark
{
using namespace gal;
using namespace gal::cga2;
}
#elif GABENCHMARK_D_DIMENSIONS == 3
#include <gal/cga.hpp>
namespace gabenchmark
{
using namespace gal;
using namespace gal::cga;
}
#endif
#else
// Only support the IK benchmark for now
#define GABENCHMARK_DOES_NOT_IMPLEMENT_THE_MODEL
#endif