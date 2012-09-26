#pragma once

#include "dft.hpp"
#include "fft_radix2.hpp"
#include "fft_odd.hpp"
#include <assert.h>

namespace dft
{
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

  struct FFT2
  {
    template<typename CmplxIter, typename SrcIter>
    void operator()(CmplxIter dstBegin, SrcIter srcBegin,
      size_t dstStride, size_t srcStride)
    {
      dftRadix2(dstBegin, srcBegin, dstStride, srcStride);
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

  template<typename Alg>
  struct FFTCooleyTukey
  {
    template<typename CmplxIter1, typename SrcIter, typename CmplxIter2>
    void fftDecomposed(CmplxIter1 dstBegin, SrcIter srcBegin, CmplxIter2 tmpBegin,
      const Twiddler<typename CmplxIter1::value_type::value_type>& twiddler, size_t N,
      FactorIt factorsBegin, FactorIt factorsEnd,
      size_t dstStride, size_t srcStride, size_t tmpStride)
    {
      auto factors = std::distance(factorsBegin, factorsEnd);

      if(factors < 2)
      {
        if(factors != 0)
          Alg()(dstBegin, srcBegin, twiddler, N, dstStride, srcStride);
        return;
      }

      size_t cols = (*factorsBegin);
      size_t rows = N/cols;

      if(factors == 2)
        Alg()(tmpBegin, srcBegin, twiddler, rows, cols, tmpStride, srcStride);
      else
        fftColsDecomposed<Alg>(tmpBegin, srcBegin, dstBegin,
          twiddler, rows, cols, ++factorsBegin, factorsEnd,
          tmpStride, srcStride, dstStride);

      size_t twiddleStride = twiddler.N/N;
      for(size_t r = 1; r < rows; ++r)
      {
        for(size_t c = 1; c < cols; ++c)
        {
          tmpBegin[(r*cols + c)*tmpStride] *=
            twiddler.factors[r*c*twiddleStride];
        }
      }

      transposeDftCols<Alg>(dstBegin, tmpBegin, twiddler, rows, cols, dstStride, tmpStride);
    }
  };
} // namespace dft
