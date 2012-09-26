#pragma once

#include "twiddler.hpp"
#include "fft_algorithms.hpp"
#include "../math/factor.hpp"

namespace dft
{
  namespace detail
  {
    template<typename Alg, typename CmplxIter, typename SrcIter>
    void dftCols(CmplxIter dstBegin, SrcIter srcBegin,
      const Twiddler<typename CmplxIter::value_type::value_type>& twiddler,
      size_t rows, size_t cols, size_t dstStride = 1, size_t srcStride = 1)
    {
      for(size_t c = 0; c < cols; ++c)
      {
        FFT<Alg>()(dstBegin + c*dstStride, srcBegin + c*srcStride,
          twiddler, rows, cols*dstStride, cols*srcStride);
      }
    }


    template<typename Alg, typename CmplxIter, typename SrcIter>
    void transposeDftCols(CmplxIter dstBegin, SrcIter srcBegin,
      const Twiddler<typename CmplxIter::value_type::value_type>& twiddler,
      size_t rows, size_t cols, size_t dstStride = 1, size_t srcStride = 1)
    {
      for(size_t r = 0; r < rows; ++r)
      {
        FFT<Alg>()(dstBegin + r*dstStride, srcBegin + r*cols*srcStride,
          twiddler, cols, rows*dstStride, srcStride);
      }
    }

    typedef std::vector<size_t>::iterator FactorIt;

    template<typename Alg, typename CmplxIter1, typename SrcIter, typename CmplxIter2>
    void fftColsDecomposed(CmplxIter1 dstBegin, SrcIter srcBegin, CmplxIter2 tmpBegin,
      const Twiddler<typename CmplxIter1::value_type::value_type>& twiddler,
      size_t rows,  size_t cols, FactorIt factorsBegin, FactorIt factorsEnd,
      size_t dstStride, size_t srcStride, size_t tmpStride);

    template<typename Alg, typename CmplxIter1, typename SrcIter, typename CmplxIter2>
    void fftDecomposed(CmplxIter1 dstBegin, SrcIter srcBegin, CmplxIter2 tmpBegin,
      const Twiddler<typename CmplxIter1::value_type::value_type>& twiddler, size_t N,
      FactorIt factorsBegin, FactorIt factorsEnd,
      size_t dstStride, size_t srcStride, size_t tmpStride)
    {
      auto factors = std::distance(factorsBegin, factorsEnd);

      if(factors < 2)
      {
        if(factors != 0)
          FFT<Alg>()(dstBegin, srcBegin, twiddler, N, dstStride, srcStride);
        return;
      }

      size_t cols = (*factorsBegin);
      size_t rows = N/cols;

      if(factors == 2)
        dftCols<Alg>(tmpBegin, srcBegin, twiddler, rows, cols, tmpStride, srcStride);
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

    template<typename Alg, typename CmplxIter1, typename SrcIter, typename CmplxIter2>
    void fftColsDecomposed(CmplxIter1 dstBegin, SrcIter srcBegin, CmplxIter2 tmpBegin,
      const Twiddler<typename CmplxIter1::value_type::value_type>& twiddler,
      size_t rows,  size_t cols, FactorIt factorsBegin, FactorIt factorsEnd,
      size_t dstStride, size_t srcStride, size_t tmpStride)
    {
      for(size_t c = 0; c < cols; ++c)
      {
        fftDecomposed<Alg>(dstBegin + c*dstStride, srcBegin + c*srcStride,
          tmpBegin, twiddler, rows, factorsBegin, factorsEnd,
          dstStride*cols, srcStride*cols, tmpStride);
      }
    }

    template<typename CmplxIter, typename SrcIter>
    void fftMultiple(CmplxIter dstBegin, size_t dstSampleStride, size_t dstTransformStride,
      SrcIter srcBegin, size_t srcSampleStride, size_t srcTransformStride,
      size_t transforms, size_t N, bool fwd)
    {
      typedef typename CmplxIter::value_type Complex_t;
      typedef typename Complex_t::value_type Float_t;
      Twiddler<Float_t> twiddler(N, fwd, N);
      std::vector<size_t> factors = factor(N);
      std::vector<Complex_t> tmp(N);

      fftColsDecomposed<Odd>(dstBegin, srcBegin, tmp.begin(), twiddler, N, transforms,
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







