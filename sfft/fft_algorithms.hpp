/* Copyright 2012, Cory Brinck
 * Distributed under the terms version 3 of the GNU General Public License
 * without any warranty. */

#pragma once

#include "dft.hpp"
#include "fft_radix2.hpp"
#include "fft_odd.hpp"

namespace sfft
{
  struct DFT
  {
    template<typename CmplxIter, typename SrcIter>
    void operator()(CmplxIter dstBegin, SrcIter srcBegin,
      const Twiddler<typename DereferencedType<CmplxIter>::value_type>& twiddler,
      size_t N, size_t dstStride, size_t srcStride)
    {
      dft(dstBegin, srcBegin, twiddler, N, dstStride, srcStride);
    }
  };

  struct FFTRadix2
  {
    template<typename CmplxIter1, typename SrcIter, typename CmplxIter2>
    void operator()(CmplxIter1 dstBegin, SrcIter srcBegin, CmplxIter2 tmpBegin,
      const Twiddler<typename DereferencedType<CmplxIter1>::value_type>& twiddler,
      size_t N, size_t dstStride, size_t srcStride, size_t tmpStride)
    {
      size_t power = math::numFactors2(N);
      detail::fftRadix2(dstBegin, srcBegin, tmpBegin, twiddler, power, dstStride, srcStride, tmpStride);
    }
  };

  struct FFTOdd
  {
    template<typename CmplxIter, typename SrcIter>
    void operator()(CmplxIter dstBegin, SrcIter srcBegin,
      const Twiddler<typename DereferencedType<CmplxIter>::value_type>& twiddler,
      size_t N, size_t dstStride, size_t srcStride)
    {
      fftOdd(dstBegin, srcBegin, twiddler, N, dstStride, srcStride);
    }

    template<typename CmplxIter1, typename SrcIter, typename CmplxIter2>
    void operator()(CmplxIter1 dstBegin, SrcIter srcBegin, CmplxIter2,
      const Twiddler<typename DereferencedType<CmplxIter1>::value_type>& twiddler,
      size_t N, size_t dstStride, size_t srcStride, size_t)
    {
      operator()(dstBegin, srcBegin, twiddler, N, dstStride, srcStride);
    }
  };

  template<typename Alg>
  struct FFTCols
  {
    template<typename CmplxIter1, typename SrcIter, typename CmplxIter2>
    void operator()(CmplxIter1 dstBegin, SrcIter srcBegin, CmplxIter2 tmpBegin,
      const Twiddler<typename DereferencedType<CmplxIter1>::value_type>& twiddler,
      size_t rows, size_t cols, size_t dstStride, size_t srcStride, size_t tmpStride)
    {
      for(size_t c = 0; c < cols; ++c)
      {
        Alg()(dstBegin + c*dstStride, srcBegin + c*srcStride, tmpBegin,
          twiddler, rows, cols*dstStride, cols*srcStride, tmpStride);
      }
    }
  };

} // namespace sfft
