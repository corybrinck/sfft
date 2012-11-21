/* Copyright 2012, Cory Brinck
 * Distributed under the terms version 3 of the GNU General Public License
 * without any warranty. */

#pragma once

#include "twiddler.hpp"

namespace dft
{
  template<typename CmplxIter, typename SrcIter>
  void fftOdd(CmplxIter dstBegin, SrcIter srcBegin,
    const Twiddler<typename CmplxIter::value_type::value_type>& twiddler,
    size_t N, size_t dstStride = 1, size_t srcStride = 1)
  {
    typedef typename CmplxIter::value_type Complex_t;
    typedef typename SrcIter::value_type Source_t;
    size_t L = N/2 + 1;

    size_t twiddleStride = twiddler.N/N;

    for(size_t l = 0; l < L; ++l)
      dstBegin[l*dstStride] = *srcBegin;
    for(size_t n = L; n < N; ++n)
      dstBegin[n*dstStride] = Complex_t(0,0);

    for(size_t l = 1; l < L; ++l)
    {
      Complex_t sum = srcBegin[l*srcStride] + srcBegin[(N-l)*srcStride];
      Complex_t dif = srcBegin[l*srcStride] - srcBegin[(N-l)*srcStride];
      Complex_t jdif(-dif.imag(),dif.real());
      *dstBegin += sum;

      for(size_t m = 1, n = l; m < L; n -= N)
      {
        for(; n < L && m < L; n += l, ++m)
        {
          const Complex_t& factor = twiddler.factors[twiddleStride*n];
          dstBegin[m*dstStride] += factor.real()*sum;
          dstBegin[(N-m)*dstStride] -= factor.imag()*jdif;
        }
        for(; n < N && m < L; n += l, ++m)
        {
          const Complex_t& factor = twiddler.factors[twiddleStride*(N - n)];
          dstBegin[m*dstStride] += factor.real()*sum;
          dstBegin[(N-m)*dstStride] += factor.imag()*jdif;
        }
      }
    }

    for(size_t l = 1; l < L; ++l)
    {
      size_t idx1 = l*dstStride, idx2 = (N-l) * dstStride;
      Complex_t tmp = dstBegin[idx1];
      dstBegin[idx1] -= dstBegin[idx2];
      dstBegin[idx2] += tmp;
    }
  }

} // namespace dft








