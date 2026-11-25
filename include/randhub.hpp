#ifndef RANDHUB_H
#define RANDHUB_H

#include <limits>
#include <chrono>
#include <numeric>
#include <random>

template <class RealT = double, class IntT = unsigned long>
class RandHub
{
    // This class contains an RNG engine and some convenience attributes
    // and methods for drawing random numbers from various distributions.
    std::mt19937_64 mt;
    std::uniform_real_distribution<RealT> real_dist{ 0.0, 1.0 };
    std::uniform_int_distribution<IntT> int_dist{ std::numeric_limits<IntT>::min(),
                                                   std::numeric_limits<IntT>::max() };
    std::normal_distribution<RealT> norm_dist{ 0.0, 1.0 };

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
    RealT rand(const RealT lbound, const RealT ubound)
    {
        // Set the parameters of the distribution.
        real_dist.param(typename std::uniform_real_distribution<RealT>::param_type(
        lbound, std::nexttoward(ubound, std::numeric_limits<long double>::max())));

        return std::min(real_dist(mt), ubound);
    }

    RealT rand(const RealT ubound = 1.0) { return rand(0.0, ubound); }

    // Real number in [0,ubound).
    RealT rand_exc(const RealT lbound, const RealT ubound)
    {
        return std::min(rand(lbound, ubound), std::nexttoward(ubound, static_cast<RealT>(0.0)));
    }

    RealT rand_exc(const RealT ubound = 1.0) { return rand_exc(0.0, ubound); }

    // Real number in (0,ubound).
    RealT rand_dbl_exc(const RealT lbound, const RealT ubound)
    {
        return std::max(std::nexttoward(static_cast<RealT>(0.0), 1.0), rand_exc(lbound, ubound));
    }

    RealT rand_dbl_exc(const RealT ubound = 1.0) { return rand_dbl_exc(0.0, ubound); }

    // Integer in [0,ubound] for ubound < int_t::max.
    IntT rand_int(const IntT lbound, const IntT ubound)
    {
        // Set the parameters of the distribution.
        int_dist.param(typename std::uniform_int_distribution<IntT>::param_type(lbound, ubound));
        return int_dist(mt);
    }

    IntT rand_int(const IntT ubound = std::numeric_limits<IntT>::max()) { return rand_int(0, ubound); }

    // Normal distribution.
    RealT rand_norm(const RealT mean = 0.0, const RealT stddev = 1.0)
    {
        // Set the parameters of the distribution.
        norm_dist.param(typename std::normal_distribution<RealT>::param_type(mean, stddev));
        return norm_dist(mt);
    }

    // Same as rand().
    RealT operator()() { return rand(); }

    // Seed the RNG engine.
    void seed(const unsigned s) { mt.seed(s); }
};

#endif // RANDHUB_H
