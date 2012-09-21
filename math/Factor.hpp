#pragma once

#include <vector>

namespace math
{
  std::vector<size_t> getPossiblePrimeFactors(size_t n)
  {
    std::vector<size_t> primes;
    primes.push_back(2);
    primes.push_back(3);
    primes.push_back(5);
    primes.push_back(7);
    primes.push_back(11);
    primes.push_back(13);
    primes.push_back(17);
    primes.push_back(19);
    primes.push_back(23);

    size_t sqrt_n = sqrt(n);
    while(!primes.empty() && primes.back() > sqrt_n)
    {
      if(primes.back() == n)
        return std::vector<size_t>(1,n);
      primes.pop_back();
    }
    return primes;
  }

  size_t numFactors2(size_t n)
  {
    size_t count = 0;
    if(n > 0)
    {
      while((n & 1) == 0)
      {
        n = n >> 1;
        ++count;
      }
    }
    return count;
  }

  std::vector<size_t> factor(size_t n)
  {
    std::vector<size_t> factors(numFactors2(n), 2);
    n = n >> factors.size();
    std::vector<size_t> primes = getPossiblePrimeFactors(n);
	
    for(size_t i = 1; n > 1 && i < primes.size(); ++i)
    {
      size_t factored = n/primes[i];
      while(n > 1 && factored*primes[i] == n)
      {
        factors.push_back(primes[i]);
        n = factored;
        factored /= primes[i];
      }
    }

    if(n > 1) factors.push_back(n);

    return factors;
  }

} //namespace math

