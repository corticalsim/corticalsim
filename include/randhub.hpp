#ifndef RANDHUB_H
#define RANDHUB_H

#include <limits>
#include <chrono>
#include <numeric>
#include <random>

template <class Real = double, class Int = unsigned long>
class RandHub
{
    // This class contains an RNG engine and some convenience attributes
    // and methods for drawing random numbers from various distributions.
    std::mt19937_64 mt;
    std::uniform_real_distribution<Real> real_dist{ 0.0, 1.0 };
    std::uniform_int_distribution<Int> int_dist{ std::numeric_limits<Int>::min(),
                                                   std::numeric_limits<Int>::max() };
    std::normal_distribution<Real> norm_dist{ 0.0, 1.0 };

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
    Real rand(const Real lbound, const Real ubound)
    {
        // Set the parameters of the distribution.
        real_dist.param(typename std::uniform_real_distribution<Real>::param_type(
        lbound, std::nexttoward(ubound, std::numeric_limits<long double>::max())));

        return std::min(real_dist(mt), ubound);
    }

    Real rand(const Real ubound = 1.0) { return rand(0.0, ubound); }

    // Real number in [0,ubound).
    Real rand_exc(const Real lbound, const Real ubound)
    {
        return std::min(rand(lbound, ubound), std::nexttoward(ubound, static_cast<Real>(0.0)));
    }

    Real rand_exc(const Real ubound = 1.0) { return rand_exc(0.0, ubound); }

    // Real number in (0,ubound).
    Real rand_dbl_exc(const Real lbound, const Real ubound)
    {
        return std::max(std::nexttoward(static_cast<Real>(0.0), 1.0), rand_exc(lbound, ubound));
    }

    Real rand_dbl_exc(const Real ubound = 1.0) { return rand_dbl_exc(0.0, ubound); }

    // Integer in [0,ubound] for ubound < int_t::max.
    Int rand_int(const Int lbound, const Int ubound)
    {
        // Set the parameters of the distribution.
        int_dist.param(typename std::uniform_int_distribution<Int>::param_type(lbound, ubound));
        return int_dist(mt);
    }

    Int rand_int(const Int ubound = std::numeric_limits<Int>::max()) { return rand_int(0, ubound); }

    // Normal distribution.
    Real rand_norm(const Real mean = 0.0, const Real stddev = 1.0)
    {
        // Set the parameters of the distribution.
        norm_dist.param(typename std::normal_distribution<Real>::param_type(mean, stddev));
        return norm_dist(mt);
    }

    // Same as rand().
    Real operator()() { return rand(); }

    // Seed the RNG engine.
    void seed(const unsigned s) { mt.seed(s); }
};

#endif // RANDHUB_H
