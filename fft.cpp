// fft.cpp - see fft.h
#include "fft.h"
#include <cmath>
#include <algorithm>

namespace lora {

static const double PI = 3.14159265358979323846;

std::size_t next_pow2(std::size_t n) {
    std::size_t p = 1;
    while (p < n) p <<= 1;
    return p;
}

void fft_pow2(cvec& a, bool invert) {
    const std::size_t n = a.size();
    if (n <= 1) return;

    for (std::size_t i = 1, j = 0; i < n; ++i) {
        std::size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }

    for (std::size_t len = 2; len <= n; len <<= 1) {
        double ang = 2.0 * PI / (double)len * (invert ? 1.0 : -1.0);
        cd wlen(std::cos(ang), std::sin(ang));
        for (std::size_t i = 0; i < n; i += len) {
            cd w(1.0, 0.0);
            for (std::size_t k = 0; k < len / 2; ++k) {
                cd u = a[i + k];
                cd v = a[i + k + len / 2] * w;
                a[i + k]           = u + v;
                a[i + k + len / 2] = u - v;
                w *= wlen;
            }
        }
    }

    if (invert) {
        double inv = 1.0 / (double)n;
        for (auto& x : a) x *= inv;
    }
}

DftPlan::DftPlan(std::size_t N) : N_(N), pow2_(true), M_(0) {
    if (N_ == 0) return;
    pow2_ = (N_ & (N_ - 1)) == 0;
    if (pow2_) return;

    M_ = next_pow2(2 * N_ - 1);
    w_.resize(N_);
    for (std::size_t n = 0; n < N_; ++n) {
        unsigned long long nn = ((unsigned long long)n * (unsigned long long)n) % (2ULL * (unsigned long long)N_);
        double ang = -PI * (double)nn / (double)N_;
        w_[n] = cd(std::cos(ang), std::sin(ang));
    }
    cvec b(M_, cd(0.0, 0.0));
    b[0] = cd(1.0, 0.0);
    for (std::size_t m = 1; m < N_; ++m) {
        cd v = std::conj(w_[m]);
        b[m]      = v;
        b[M_ - m] = v;
    }
    fft_pow2(b, false);
    kernelF_ = std::move(b);
}

cvec DftPlan::forward(const cd* x, std::size_t len) const {
    if (N_ == 0) return cvec();
    const std::size_t L = std::min(len, N_);

    if (pow2_) {
        cvec a(N_, cd(0.0, 0.0));
        for (std::size_t i = 0; i < L; ++i) a[i] = x[i];
        fft_pow2(a, false);
        return a;
    }

    cvec A(M_, cd(0.0, 0.0));
    for (std::size_t n = 0; n < L; ++n) A[n] = x[n] * w_[n];
    fft_pow2(A, false);
    for (std::size_t i = 0; i < M_; ++i) A[i] *= kernelF_[i];
    fft_pow2(A, true);

    cvec out(N_);
    for (std::size_t k = 0; k < N_; ++k) out[k] = w_[k] * A[k];
    return out;
}

}
