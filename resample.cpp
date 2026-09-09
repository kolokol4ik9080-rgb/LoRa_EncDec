// resample.cpp - see resample.h
#include "resample.h"
#include <cmath>
#include <algorithm>

namespace lora {

static const double PI = 3.14159265358979323846;

long gcd_ll(long a, long b) {
    a = std::labs(a);
    b = std::labs(b);
    while (b) { long t = a % b; a = b; b = t; }
    return a;
}

static double sinc(double x) {              // normalized sinc: sin(pi x)/(pi x)
    if (std::fabs(x) < 1e-12) return 1.0;
    double px = PI * x;
    return std::sin(px) / px;
}

static double bessel_i0(double x) {         
    double sum = 1.0, term = 1.0, y = x * x / 4.0;
    for (int k = 1; k < 60; ++k) {
        term *= y / (double)(k * k);
        sum  += term;
        if (term < 1e-14 * sum) break;
    }
    return sum;
}

static std::vector<double> kaiser(int L, double beta) {
    std::vector<double> w(L);
    double denom = bessel_i0(beta);
    double M = (L - 1) / 2.0;
    for (int n = 0; n < L; ++n) {
        double r   = (M > 0) ? (n - M) / M : 0.0;
        double arg = beta * std::sqrt(std::max(0.0, 1.0 - r * r));
        w[n] = bessel_i0(arg) / denom;
    }
    return w;
}

static cvec upfirdn(const cvec& x, const std::vector<double>& h, long p, long q) {
    long Lx = (long)x.size();
    if (Lx == 0) return cvec();
    long Lu = (Lx - 1) * p + 1;
    long Lh = (long)h.size();
    long Ly = Lu + Lh - 1;
    cvec y((size_t)std::max<long>(Ly, 0), cd(0.0, 0.0));

    for (long i = 0; i < Lx; ++i) {
        cd xu = x[i];
        long base = i * p;
        for (long k = 0; k < Lh; ++k) y[base + k] += xu * h[k];
    }

    long Lo = (Ly > 0) ? (Ly + q - 1) / q : 0;
    cvec out((size_t)Lo);
    for (long i = 0; i < Lo; ++i) out[i] = y[i * q];
    return out;
}

cvec resample(const cvec& x, long p, long q) {
    long g = gcd_ll(p, q);
    p /= g; q /= g;
    if (p == 1 && q == 1) return x;

    const int    N    = 10;
    const double beta = 5.0;
    long   pqmax  = std::max(p, q);
    double fc     = 1.0 / 2.0 / (double)pqmax;
    int    L      = (int)(2 * N * pqmax + 1);
    double M      = (L - 1) / 2.0;
    double cutoff = 2.0 * fc;

    std::vector<double> kw = kaiser(L, beta);
    std::vector<double> h(L);
    for (int n = 0; n < L; ++n)
        h[n] = cutoff * sinc(cutoff * ((double)n - M)) * kw[n];

    double s = 0.0;
    for (int n = 0; n < L; n += p) s += h[n];
    for (int n = 0; n < L; ++n) h[n] = (double)p * h[n] / s;

    double Lhalf = (L - 1) / 2.0;
    long   Lx    = (long)x.size();
    long   nz    = (long)std::floor((double)q - std::fmod(Lhalf, (double)q));
    std::vector<double> hh;
    hh.reserve((size_t)(L + nz + 8));
    for (long i = 0; i < nz; ++i) hh.push_back(0.0);
    for (int  n = 0; n < L;  ++n) hh.push_back(h[n]);
    double Lhalf2 = Lhalf + nz;
    long   delay  = (long)std::floor(std::ceil(Lhalf2) / (double)q);

    long nz1 = 0;
    auto target = (long)std::ceil((double)Lx * p / (double)q);
    while ((long)std::ceil(((double)((Lx - 1) * p + (long)hh.size() + nz1)) / (double)q) - delay < target)
        ++nz1;
    for (long i = 0; i < nz1; ++i) hh.push_back(0.0);

    cvec y  = upfirdn(x, hh, p, q);
    long Ly = target;

    cvec out;
    if ((long)y.size() > delay) {
        long avail = (long)y.size() - delay;
        long take  = std::min(Ly, avail);
        out.assign(y.begin() + delay, y.begin() + delay + take);
    }
    out.resize((size_t)Ly, cd(0.0, 0.0));
    return out;
}

cvec lowpass(const cvec& x, double fpass, double fs) {
    double nyq = fs / 2.0;
    double fcn = fpass / nyq;
    if (fcn >= 1.0 || x.empty()) return x;

    const double steep = 0.85;
    double tw = (1.0 - steep) * (1.0 - fcn);
    if (tw < 1e-3) tw = 1e-3;

    const double A = 60.0;
    double beta = (A > 50.0) ? 0.1102 * (A - 8.7)
                : (A >= 21.0 ? 0.5842 * std::pow(A - 21.0, 0.4) + 0.07886 * (A - 21.0) : 0.0);

    int Nord = (int)std::ceil((A - 8.0) / (2.285 * tw * PI));
    if (Nord < 4) Nord = 4;
    int L = Nord + 1;
    if ((L % 2) == 0) ++L;
    int Mid = (L - 1) / 2;
    double cutoff = fcn + tw / 2.0;

    std::vector<double> kw = kaiser(L, beta);
    std::vector<double> h(L);
    double s = 0.0;
    for (int n = 0; n < L; ++n) { h[n] = cutoff * sinc(cutoff * (double)(n - Mid)) * kw[n]; s += h[n]; }
    for (int n = 0; n < L; ++n) h[n] /= s;

    long Lx = (long)x.size();
    cvec y((size_t)Lx, cd(0.0, 0.0));
    for (long n = 0; n < Lx; ++n) {
        cd acc(0.0, 0.0);
        long lo = std::max<long>(0, n + Mid - (L - 1));
        long hi = std::min<long>(Lx - 1, n + Mid);
        for (long idx = lo; idx <= hi; ++idx) acc += h[(int)(n + Mid - idx)] * x[idx];
        y[n] = acc;
    }
    return y;
}

}
