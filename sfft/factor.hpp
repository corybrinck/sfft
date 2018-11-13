/* Copyright 2012, Cory Brinck
 * Distributed under the terms version 3 of the GNU General Public License
 * without any warranty. */

#pragma once

#include <vector>
#include <cmath>
namespace sfft
{
  namespace math
  {
    inline std::vector<std::size_t> getPossiblePrimeFactors(std::size_t n)
    {
      std::vector<std::size_t> primes = {2,3,5,7,11,13,17,19,23};

      std::size_t sqrt_n = static_cast<std::size_t>(std::sqrt<size_t>(n));
      while (!primes.empty() && primes.back() > sqrt_n)
      {
        if (primes.back() == n)
          return std::vector<std::size_t>(1, n);
        primes.pop_back();
      }
      return primes;
    }

    inline std::size_t numFactors2(std::size_t n)
    {
      std::size_t count = 0;
      if (n > 0)
      {
        while ((n & 1) == 0)
        {
          n = n >> 1;
          ++count;
        }
      }
      return count;
    }

    inline std::vector<std::size_t> factor(std::size_t n)
    {
      std::vector<std::size_t> factors(numFactors2(n), 2);
      n = n >> factors.size();
      std::vector<std::size_t> primes = getPossiblePrimeFactors(n);

      for (std::size_t i = 1; n > 1 && i < primes.size(); ++i)
      {
        std::size_t factored = n / primes[i];
        while (n > 1 && factored*primes[i] == n)
        {
          factors.push_back(primes[i]);
          n = factored;
          factored /= primes[i];
        }
      }

      if (n > 1) factors.push_back(n);

      return factors;
    }

  } //namespace math
}

