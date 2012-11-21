/* Copyright 2012, Cory Brinck
 * Distributed under the terms version 3 of the GNU General Public License
 * without any warranty. */

#pragma once

#include "twiddler.hpp"
#include "fft_algorithms.hpp"

namespace sfft
{
  template<typename Alg>
  struct FFTCooleyTukey
  {
    template<typename CmplxIter, typename SrcIter>
    void transposeDftCols(CmplxIter dstBegin, SrcIter srcBegin,
      const Twiddler<typename CmplxIter::value_type::value_type>& twiddler,
      size_t rows, size_t cols, size_t dstStride = 1, size_t srcStride = 1)
    {
      for(size_t r = 0; r < rows; ++r)
      {
        FFTOdd()(dstBegin + r*dstStride, srcBegin + r*cols*srcStride,
          twiddler, cols, rows*dstStride, srcStride);
      }
    }

    typedef std::vector<size_t>::iterator FactorIt;

    template<typename CmplxIter1, typename SrcIter, typename CmplxIter2>
    void fftColsDecomposed(CmplxIter1 dstBegin, SrcIter srcBegin, CmplxIter2 tmpBegin,
      const Twiddler<typename CmplxIter1::value_type::value_type>& twiddler,
      size_t rows,  size_t cols, FactorIt factorsBegin, FactorIt factorsEnd,
      size_t dstStride, size_t srcStride, size_t tmpStride)
    {
      for(size_t c = 0; c < cols; ++c)
      {
        fftDecomposed(dstBegin + c*dstStride, srcBegin + c*srcStride, tmpBegin, twiddler,
          rows, factorsBegin, factorsEnd, dstStride*cols, srcStride*cols, tmpStride);
      }
    }

    template<typename CmplxIter1, typename SrcIter, typename CmplxIter2>
    void fftDecomposed(CmplxIter1 dstBegin, SrcIter srcBegin, CmplxIter2 tmpBegin,
      const Twiddler<typename CmplxIter1::value_type::value_type>& twiddler, size_t N,
      FactorIt factorsBegin, FactorIt factorsEnd,
      size_t dstStride, size_t srcStride, size_t tmpStride)
    {
      size_t factors = std::distance(factorsBegin, factorsEnd);

      if(factors < 2)
      {
        if(factors != 0) Alg()(dstBegin, srcBegin, twiddler, N, dstStride, srcStride);
        return;
      }

      size_t cols = (*factorsBegin);
      size_t rows = N/cols;

      if(factors == 2)
        FFTCols<Alg>()(tmpBegin, srcBegin, twiddler, rows, cols, tmpStride, srcStride);
      else
        fftColsDecomposed(tmpBegin, srcBegin, dstBegin, twiddler, rows, cols,
          ++factorsBegin, factorsEnd, tmpStride, srcStride, dstStride);

      size_t twiddleStride = twiddler.N/N;

      for(size_t r = 1; r < rows; ++r)
        for(size_t c = 1; c < cols; ++c)
          tmpBegin[(r*cols + c)*tmpStride] *= twiddler.factors[r*c*twiddleStride];

      transposeDftCols(dstBegin, tmpBegin, twiddler, rows, cols, dstStride, tmpStride);
    }
  }; // struct FFTCooleyTukey
}
