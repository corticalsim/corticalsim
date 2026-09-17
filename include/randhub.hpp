#ifndef RANDHUB_H
#define RANDHUB_H

#include <limits>
#include <chrono>
#include <numeric>
#include <random>

template <class real_t = double, class int_t = unsigned long>
class RandHub
{
    // This class contains an RNG engine and some convenience attributes
    // and methods for drawing random numbers from various distributions.
    std::mt19937_64 mt;
    std::uniform_real_distribution<real_t> real_dist{ 0.0, 1.0 };
    std::uniform_int_distribution<int_t> int_dist{ std::numeric_limits<int_t>::min(),
                                                   std::numeric_limits<int_t>::max() };
    std::normal_distribution<real_t> norm_dist{ 0.0, 1.0 };

  public:

    RandHub():
        mt(std::random_device{}())
    {
    }

    RandHub(unsigned s):
        mt(s)
    {
    }

    // Real number in [0,1].
    //
    // NOTE
    // Technically, the standard says that `uniform_real_distribution` should
    // return a number in the *half-open* interval [0, ubound).
    // However, there seems to be a bug in most compiler implementations
    // whereby on rare occasions it could also return the upper bound,
    // thus inadvertently matching the desired behaviour for this method.
    // However, we cannot rely on that in case it is 'fixed' later on.
    real_t rand(const real_t lbound, const real_t ubound)
    {
        // Set the parameters of the distribution.
        real_dist.param(typename std::uniform_real_distribution<real_t>::param_type(
        lbound, std::nexttoward(ubound, std::numeric_limits<long double>::max())));

        return std::min(real_dist(mt), ubound);
    }

    real_t rand(const real_t ubound = 1.0) { return rand(0.0, ubound); }

    // Real number in [0,ubound).
    real_t randExc(const real_t lbound, const real_t ubound)
    {
        return std::min(rand(lbound, ubound), std::nexttoward(ubound, static_cast<real_t>(0.0)));
    }

    real_t randExc(const real_t ubound = 1.0) { return randExc(0.0, ubound); }

    // Real number in (0,ubound).
    real_t randDblExc(const real_t lbound, const real_t ubound)
    {
        return std::max(std::nexttoward(static_cast<real_t>(0.0), 1.0), randExc(lbound, ubound));
    }

    real_t randDblExc(const real_t ubound = 1.0) { return randDblExc(0.0, ubound); }

    // Integer in [0,ubound] for ubound < int_t::max.
    int_t randInt(const int_t lbound, const int_t ubound)
    {
        // Set the parameters of the distribution.
        int_dist.param(typename std::uniform_int_distribution<int_t>::param_type(lbound, ubound));
        return int_dist(mt);
    }

    int_t randInt(const int_t ubound = std::numeric_limits<int_t>::max()) { return randInt(0, ubound); }

    // Normal distribution.
    real_t randNorm(const real_t mean = 0.0, const real_t stddev = 1.0)
    {
        // Set the parameters of the distribution.
        norm_dist.param(typename std::normal_distribution<real_t>::param_type(mean, stddev));
        return norm_dist(mt);
    }

    // Same as rand().
    real_t operator()() { return rand(); }

    // Seed the RNG engine.
    void seed(const unsigned s) { mt.seed(s); }
};

#endif // RANDHUB_H
