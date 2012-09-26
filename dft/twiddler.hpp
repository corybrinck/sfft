#pragma once

#include <complex>
#include <vector>
#include <cmath>
#include <limits>

namespace dft { namespace detail
{
  template<typename T>
  std::vector<std::complex<T> > getTwiddleFactors(size_t N, bool fwd,
    size_t Lm1Qm1 = std::numeric_limits<size_t>::max)
  {
    std::vector<std::complex<T> > factors(std::min(N, Lm1Qm1));
    if(factors.size() > 0)
    {
      factors[0] = 1;
      if(factors.size() > 1)
      {
        // Use high-precision multiplication of phasors.
        typedef std::complex<long double> HighPrecComplex;
        long double phase = 2*3.14159265358979323846/N;
        HighPrecComplex phasor(cos(phase), fwd ? -sin(phase) : sin(phase));
        HighPrecComplex tmp(phasor);
        factors[1] = tmp;
        for(size_t i = 2; i < factors.size(); ++i)
          factors[i] = tmp = tmp * phasor;
      }
    }
    return factors;
  }
} // namespace detail

template<typename T>
struct Twiddler
{
  Twiddler(size_t N_, bool fwd,
    size_t Lm1Qm1 = std::numeric_limits<size_t>::max())
    :factors(detail::getTwiddleFactors<T>(N_, fwd, Lm1Qm1)),
    N(N_) {}

  const std::vector<std::complex<T> > factors;
  const size_t N;
};

} //namespace dft

