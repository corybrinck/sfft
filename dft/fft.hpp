#pragma once

#include "twiddler.hpp"
#include "fft_algorithms.hpp"
#include "../math/factor.hpp"

namespace dft
{
  namespace detail
  {
    template<typename CmplxIter, typename SrcIter>
    void fftMultiple(CmplxIter dstBegin, size_t dstSampleStride, size_t dstTransformStride,
      SrcIter srcBegin, size_t srcSampleStride, size_t srcTransformStride,
      size_t transforms, size_t N, bool fwd)
    {
      typedef typename CmplxIter::value_type Complex_t;
      typedef typename Complex_t::value_type Float_t;
      Twiddler<Float_t> twiddler(N, fwd, N);
      std::vector<size_t> factors = math::factor(N);
      std::vector<Complex_t> tmp(N);

      // TODO: FFTOdd should only be used after removing all factors of 2
      FFTCooleyTukey<FFTOdd>::fftColsDecomposed(dstBegin, srcBegin, tmp.begin(), twiddler, N, transforms,
        factors.begin(), factors.end(), dstSampleStride, srcSampleStride, 1);
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
} // namespace dft







