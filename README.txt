sfft

A simple templated Fast Fourier Transform (FFT) library written in C++. Performs FFTs of arbitrary sizes on arbitrary data-types using a templated interface with random-access iterators, such as pointers.  The output sequence is expected to use std::complex or a type with a similar public interface.

EXAMPLE:
---------------------------------------------------------------------------------------

#include "sfft/fft.hpp"
#include <deque>
#include <iostream>
#include <vector>

// Forward and inverse transform a delta function
void test_sfft()
{
  size_t N = 2*3*4*5*6*7;
  std::deque<int> delta_function(N, 0);
  delta_function[N/2] = 1.0;

  std::cout << "Delta function length: " << N << std::endl;

  std::deque<std::complex<float> > forward_output(N);
  std::vector<std::complex<double> > inverse_output(N);

  sfft::fft(forward_output.begin(), delta_function.begin(), N, true);
  sfft::fft(inverse_output.begin(), forward_output.begin(), N, false);

  std::cout << "Delta peak before transforms: "  << delta_function[N/2] << std::endl;
  std::cout << "Delta peak after transforms: "  << inverse_output[N/2] << std::endl;
}

---------------------------------------------------------------------------------------

Both in-place and out-of-place trasforms are supported for std::complex inputs. Currently a buffer of memory is allocated and freed internally with most calls regardless of in-place or out-of-place use.  Memory is also allocated internally for Twiddle factors.

Cooley-Tukey and Radix-2 algorithms have been implemented, as well as an O(N^2/4) algorithm for prime numbers. No O(N*log(N)) algorithm has yet been implemented for prime DFT sizes.

Distributed under the GPL Public Use License Version 3.
