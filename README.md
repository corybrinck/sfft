sfft
====

A simple templated Fast Fourier Transform (FFT) library written in C++. Performs FFTs of arbitrary sizes on arbitrary data-types using a templated interface with random-access iterators, such as pointers.  The output sequence is expected to use std::complex or a type with a similar public interface.

Both in-place and out-of-place trasforms are supported. Currently a buffer of memory is allocated and freed internally with most calls regardless of in-place or out-of-place use.  Memory is also allocated internally for Twiddle factors.

Cooley-Tukey and Radix-2 algorithms have been implemented, as well as an O(N^2/4) algorithm for prime numbers. No O(N*log(N)) algorithm has yet been implemented for prime DFT sizes.  

Distributed under the GPL Public Use License Version 3.