/* Copyright 2012, Cory Brinck
 * Distributed under the terms version 3 of the GNU General Public License
 * without any warranty. */

#pragma once

#include "dft.hpp"
#include "fft_radix2.hpp"
#include "fft_odd.hpp"
#include <assert.h>

namespace dft
{
  struct DFT
  {
    template<typename CmplxIter, typename SrcIter>
    void operator()(CmplxIter dstBegin, SrcIter srcBegin,
      const Twiddler<typename CmplxIter::value_type::value_type>& twiddler,
      size_t N, size_t dstStride, size_t srcStride)
    {
      dft(dstBegin, srcBegin, twiddler, N, dstStride, srcStride);
    }
  };

  struct FFTRadix2
  {
    template<typename CmplxIter, typename SrcIter>
    void operator()(CmplxIter dstBegin, SrcIter srcBegin,
      const Twiddler<typename CmplxIter::value_type::value_type>& twiddler,
      size_t N, size_t dstStride, size_t srcStride)
    {
      //TODO: Do an in-place Radix 2 FFT to avoid needing this tmp
      size_t power = math::numFactors2(N);
      std::vector<typename CmplxIter::value_type> tmp(N);
      detail::fftRadix2(dstBegin, srcBegin, tmp.begin(), twiddler, power, dstStride, srcStride, 1);
    }
  };

  struct FFTOdd
  {
    template<typename CmplxIter, typename SrcIter>
    void operator()(CmplxIter dstBegin, SrcIter srcBegin,
      const Twiddler<typename CmplxIter::value_type::value_type>& twiddler,
      size_t N, size_t dstStride, size_t srcStride)
    {
      fftOdd(dstBegin, srcBegin, twiddler, N, dstStride, srcStride);
    }
  };

  template<typename Alg>
  struct FFTCols
  {
    template<typename CmplxIter, typename SrcIter>
    void operator()(CmplxIter dstBegin, SrcIter srcBegin,
      const Twiddler<typename CmplxIter::value_type::value_type>& twiddler,
      size_t rows, size_t cols, size_t dstStride = 1, size_t srcStride = 1)
    {
      for(size_t c = 0; c < cols; ++c)
      {
        Alg()(dstBegin + c*dstStride, srcBegin + c*srcStride,
          twiddler, rows, cols*dstStride, cols*srcStride);
      }
    }
  };

} // namespace dft
