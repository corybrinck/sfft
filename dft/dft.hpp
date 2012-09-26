#pragma once

#include "twiddler.hpp"

namespace dft
{
  template<typename CmplxIter, typename SrcIter>
  void dft(CmplxIter dstBegin, SrcIter srcBegin,
    const Twiddler<typename CmplxIter::value_type::value_type>& twiddler,
    size_t N, size_t dstStride = 1, size_t srcStride = 1)
  {
    size_t twiddleStride = twiddler.N/N;
    //k == 0
    dstBegin[0] = srcBegin[0];
    for(size_t n = 1; n < N; ++n)
      dstBegin[0] += srcBegin[n*srcStride];

    for(size_t k = 1; k < N; ++k)
    {
      dstBegin[k*dstStride] = srcBegin[0];
      for(size_t n = 1; n < N; ++n)
      {
        dstBegin[k*dstStride] += srcBegin[n*srcStride]*
          twiddler.factors[twiddleStride*((k*n)%N)];
      }
    }
  }
} // namespace dft






