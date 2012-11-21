/* Copyright 2012, Cory Brinck
 * Distributed under the terms version 3 of the GNU General Public License
 * without any warranty. */

#pragma once

#include "twiddler.hpp"
#include "../math/factor.hpp"

namespace dft { namespace detail
{
  template<typename CmplxIter, typename SrcIter>
  void dftRadix2(CmplxIter dstBegin, SrcIter srcBegin,
    size_t dstStride = 1, size_t srcStride = 1)
  {
    (*dstBegin) = (*srcBegin) + srcBegin[srcStride];
    dstBegin[dstStride] = (*srcBegin) - srcBegin[srcStride];
  }

  template<typename CmplxIter, typename SrcIter>
  void dftColsRadix2(CmplxIter dstBegin, SrcIter srcBegin, size_t cols,
    size_t dstStride = 1, size_t srcStride = 1)
  {
    for(size_t c = 0; c < cols; ++c)
    {
      dftRadix2(dstBegin + c*dstStride, srcBegin + c*srcStride,
        cols*dstStride, cols*srcStride);
    }
  }

  template<typename CmplxIter, typename SrcIter>
  void transposeDftColsRadix2(CmplxIter dstBegin, SrcIter srcBegin, size_t rows,
    size_t dstStride = 1, size_t srcStride = 1)
  {
    for(size_t r = 0; r < rows; ++r)
    {
      dftRadix2(dstBegin + r*dstStride, srcBegin + 2*r*srcStride,
        rows*dstStride, srcStride);
    }
  }

  template<typename CmplxIter1, typename SrcIter, typename CmplxIter2>
  void fftRadix2Cols(CmplxIter1 dstBegin, SrcIter srcBegin, CmplxIter2 tmpBegin,
    const Twiddler<typename CmplxIter1::value_type::value_type>& twiddler,
    size_t power, size_t cols, size_t dstStride, size_t srcStride, size_t tmpStride)
  {
    for(size_t c = 0; c < cols; ++c)
    {
      fftRadix2(dstBegin + c*dstStride, srcBegin + c*srcStride, tmpBegin,
        twiddler, power, dstStride*cols, srcStride*cols, tmpStride);
    }
  }

  template<typename CmplxIter1, typename SrcIter, typename CmplxIter2>
  void fftRadix2(CmplxIter1 dstBegin, SrcIter srcBegin, CmplxIter2 tmpBegin,
    const Twiddler<typename CmplxIter1::value_type::value_type>& twiddler,
    size_t power, size_t dstStride, size_t srcStride, size_t tmpStride)
  {
    if(power < 2)
    {
      if(power == 1)
        dftRadix2(dstBegin, srcBegin, dstStride, srcStride);
      return;
    }

    if(power == 2)
      dftColsRadix2(tmpBegin, srcBegin, 2, tmpStride, srcStride);
    else
      fftRadix2Cols(tmpBegin, srcBegin, dstBegin, twiddler,
        power-1, 2, tmpStride, srcStride, dstStride);


    // Multiply
    size_t twiddleStride = twiddler.N >> power;
    size_t rows = 1 << (power-1);
    for(size_t r = 1; r < rows; ++r)
      tmpBegin[tmpStride*(2*r + 1)] *= twiddler.factors[r*twiddleStride];

    transposeDftColsRadix2(dstBegin, tmpBegin, rows, dstStride, tmpStride);
  }

  template<typename CmplxIter, typename SrcIter>
  void fftMultipleRadix2(CmplxIter dstBegin, size_t dstSampleStride, size_t dstTransformStride,
    SrcIter srcBegin, size_t srcSampleStride, size_t srcTransformStride,
    size_t transforms, size_t power, bool fwd)
  {
    typedef typename CmplxIter::value_type Complex_t;
    typedef typename Complex_t::value_type Float_t;
    size_t N = 1 << power;
    Twiddler<Float_t> twiddler(N, fwd, (N >> 1) - 1);
    std::vector<Complex_t> tmp(N);

    fftRadix2Cols(dstBegin, srcBegin, tmp.begin(), twiddler, power, transforms,
      dstSampleStride, srcSampleStride, 1);
  }
} //namespace detail

  template<typename CmplxIter, typename SrcIter>
  void fftRadix2(CmplxIter dstBegin, SrcIter srcBegin, size_t power, bool fwd,
    size_t dstStride, size_t srcStride)
  {
    size_t N = 1 << power;
    detail::fftMultipleRadix2(dstBegin, dstStride, N, srcBegin, srcStride, N, 1, power, fwd);
  }

} //namespace dft

