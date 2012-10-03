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
