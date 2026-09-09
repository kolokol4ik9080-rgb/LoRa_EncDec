// fft.h - FFT of arbitrary length (matching MATLAB `fft(x, N)`).
//
// The LoRa dechirp needs a zero-padded FFT of length fft_len = 2^(sf+1) * zpr.
// With the default zero-padding ratio of 10 this is not a power of two
// (e.g. sf=8 -> 5120 = 2^10 * 5), so we provide a Bluestein transform that
// works for any length, wrapped around a fast radix-2 core.
#pragma once
#include "lora_types.h"
#include <cstddef>

namespace lora {


void fft_pow2(cvec& a, bool invert);


std::size_t next_pow2(std::size_t n);


class DftPlan {
public:
    DftPlan() : N_(0), pow2_(true), M_(0) {}
    explicit DftPlan(std::size_t N);

    cvec forward(const cd* x, std::size_t len) const;
    cvec forward(const cvec& x) const { return forward(x.data(), x.size()); }

    std::size_t size() const { return N_; }

private:
    std::size_t N_;
    bool        pow2_;
    std::size_t M_;         // Bluestein convolution length (power of two)
    cvec        w_;         // chirp weights exp(-i*pi*n^2/N), n = 0..N-1
    cvec        kernelF_;   // FFT of the Bluestein kernel, length M_
};

} // namespace lora
