/* Copyright 2012, Cory Brinck
 * Distributed under the terms version 3 of the GNU General Public License
 * without any warranty. */

#pragma once

#include "factor.hpp"
#include "fft_cooley_tukey.hpp"

namespace sfft
{
  namespace detail
  {
    template<typename CmplxIter, typename SrcIter>
    void fftMultiple(CmplxIter dstBegin, size_t dstSampleStride, size_t dstTransformStride,
      SrcIter srcBegin, size_t srcSampleStride, size_t srcTransformStride,
      size_t transforms, size_t N, bool fwd)
    {
      typedef DereferencedType<CmplxIter> Complex_t;
      typedef typename Complex_t::value_type Float_t;
      Twiddler<Float_t> twiddler(N, fwd, N);
      size_t pow2 = math::numFactors2(N);
      std::vector<size_t> factors = math::factor(N >> pow2);
      std::vector<Complex_t> tmp(N);

      if (pow2 == 0)
      {
        if (factors.empty())
          factors.push_back(1);

        for (size_t transform = 0; transform < transforms; ++transform)
        {
          FFTCooleyTukey<FFTOdd, FFTOdd>()(
            dstBegin + transform*dstTransformStride, srcBegin + transform*srcTransformStride,
            tmp.data(), twiddler, N, factors.begin(), factors.end(),
            dstSampleStride, srcSampleStride, 1);
        }
      }
      else
      {
        factors.push_back(size_t(1) << pow2);
        for (size_t transform = 0; transform < transforms; ++transform)
        {
          FFTCooleyTukey<FFTOdd, FFTRadix2>()(
            dstBegin + transform*dstTransformStride, srcBegin + transform*srcTransformStride,
            tmp.data(), twiddler, N, factors.begin(), factors.end(),
            dstSampleStride, srcSampleStride, 1);
        }
      }
    }
  } // namespace detail

  template<typename CmplxIter, typename SrcIter>
  void fft(CmplxIter dstBegin, SrcIter srcBegin, size_t N, bool fwd,
    size_t dstStride = 1, size_t srcStride = 1)
  {
    detail::fftMultiple(dstBegin, dstStride, N, srcBegin, srcStride, N, 1, N, fwd);
  }

  template<typename CmplxIter, typename SrcIter>
  void fftRows(CmplxIter dstBegin, SrcIter srcBegin, size_t rows, size_t cols, bool fwd)
  {
    detail::fftMultiple(dstBegin, 1, cols, srcBegin, 1, cols, rows, cols, fwd);
  }

  template<typename CmplxIter, typename SrcIter>
  void fftCols(CmplxIter dstBegin, SrcIter srcBegin, size_t rows, size_t cols, bool fwd)
  {
    detail::fftMultiple(dstBegin, cols, 1, srcBegin, cols, 1, cols, rows, fwd);
  }

  template<typename CmplxIter, typename SrcIter>
  void fft2D(CmplxIter dstBegin, SrcIter srcBegin, size_t rows, size_t cols, bool fwd)
  {
    std::vector<DereferencedType<CmplxIter>> tmp(rows*cols);
    fftRows(tmp.data(), srcBegin, rows, cols, fwd);
    fftCols(dstBegin, tmp.data(), rows, cols, fwd);
  }
} // namespace sfft







