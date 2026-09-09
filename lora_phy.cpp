// lora_phy.cpp - see lora_phy.h
#ifdef _WIN32
#  define NOMINMAX
#  define WIN32_LEAN_AND_MEAN
#  include <windows.h>
#endif

#include "lora_phy.h"
#include "resample.h"
#include "gf2.h"
#include "hamming.h"
#include "crc16.h"
#include "whitening.h"

#include <cmath>
#include <algorithm>
#include <fstream>
#include <iostream>

namespace lora {

static const double PI = 3.14159265358979323846;

namespace {
#ifdef _WIN32
std::wstring to_wide(const std::string& s, UINT cp) {
    if (s.empty()) return std::wstring();
    int n = MultiByteToWideChar(cp, 0, s.data(), (int)s.size(), nullptr, 0);
    if (n <= 0) return std::wstring();
    std::wstring w((size_t)n, L'\0');
    MultiByteToWideChar(cp, 0, s.data(), (int)s.size(), &w[0], n);
    return w;
}
#endif

std::ifstream open_binary_in(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
#ifdef _WIN32
    if (!f) { f.clear(); f.open(to_wide(path, CP_UTF8).c_str(), std::ios::binary); }
    if (!f) { f.clear(); f.open(to_wide(path, CP_ACP ).c_str(), std::ios::binary); }
#endif
    return f;
}

std::ofstream open_binary_out(const std::string& path) {
    std::ofstream f(path, std::ios::binary);
#ifdef _WIN32
    if (!f) { f.clear(); f.open(to_wide(path, CP_UTF8).c_str(), std::ios::binary); }
    if (!f) { f.clear(); f.open(to_wide(path, CP_ACP ).c_str(), std::ios::binary); }
#endif
    return f;
}
}

LoRaPHY::LoRaPHY(double rf_freq_, int sf_, double bw_, double fs_)
    : rf_freq(rf_freq_), sf(sf_), bw(bw_), fs(fs_) {
    cfo_val = 0.0;
    init();
}

void LoRaPHY::init() {
    long N = (long)std::llround(std::pow(2.0, sf));
    bin_num    = N * zero_padding_ratio;
    sample_num = 2 * N;
    fft_len    = sample_num * zero_padding_ratio;

    downchirp = chirp(false, sf, bw, 2.0 * bw, 0.0, cfo_val, 0.0);
    upchirp   = chirp(true,  sf, bw, 2.0 * bw, 0.0, cfo_val, 0.0);

    ldr = (std::pow(2.0, sf) / bw > 16e-3) ? 1 : 0;

    plan = DftPlan((size_t)fft_len);
}

// chirp
cvec LoRaPHY::chirp(bool is_up, int sf, double bw, double fs,
                    double h, double cfo, double tdelta, double tscale) {
    double N = std::pow(2.0, sf);
    double T = N / bw;
    long   samp_per_sym = (long)std::llround(fs / bw * N);
    double h_orig = h;
    h = std::round(h);
    cfo = cfo + (h_orig - h) / N * bw;

    double k, f0;
    if (is_up) { k =  bw / T; f0 = -bw / 2.0 + cfo; }
    else       { k = -bw / T; f0 =  bw / 2.0 + cfo; }

    long snum = (long)std::floor(samp_per_sym * (N - h) / N) + 1;
    cvec c1((size_t)std::max<long>(snum, 0));
    for (long i = 0; i < snum; ++i) {
        double t = (double)i / fs * tscale + tdelta;
        double phase = 2.0 * PI * (t * (f0 + k * T * h / N + 0.5 * k * t));
        c1[(size_t)i] = std::exp(cd(0.0, phase));
    }
    double phi = (snum == 0) ? 0.0 : std::arg(c1[(size_t)(snum - 1)]);

    long snum2 = (long)std::floor(samp_per_sym * h / N - 1.0) + 1;
    if (snum2 < 0) snum2 = 0;
    cvec c2((size_t)snum2);
    for (long i = 0; i < snum2; ++i) {
        double t = (double)i / fs + tdelta;
        double phase = phi + 2.0 * PI * (t * (f0 + 0.5 * k * t));
        c2[(size_t)i] = std::exp(cd(0.0, phase));
    }

    cvec y;
    y.reserve((size_t)((snum > 0 ? snum - 1 : 0) + snum2));
    for (long i = 0; i < snum - 1; ++i) y.push_back(c1[(size_t)i]);
    for (long i = 0; i < snum2;    ++i) y.push_back(c2[(size_t)i]);
    return y;
}

// dechirp
LoRaPHY::Peak LoRaPHY::dechirp(long x, bool is_up) {
    const cvec& c = is_up ? downchirp : upchirp;
    if (x < 1 || (x - 1 + sample_num) > (long)sig_.size()) return { 0.0, 1 };

    cvec prod((size_t)sample_num);
    for (long i = 0; i < sample_num; ++i)
        prod[(size_t)i] = sig_[(size_t)(x - 1 + i)] * c[(size_t)i];

    cvec ft = plan.forward(prod.data(), (size_t)sample_num);

    double best = -1.0; long besti = 0;
    for (long i = 0; i < bin_num; ++i) {
        double v = std::abs(ft[(size_t)i]) + std::abs(ft[(size_t)(fft_len - bin_num + i)]);
        if (v > best) { best = v; besti = i; }
    }
    return { best, besti + 1 };
}

// detect
long LoRaPHY::detect(long start_idx) {
    long ii = start_idx;
    std::vector<long> pk_bin_list;
    while (ii < (long)sig_.size() - sample_num * preamble_len) {
        if ((long)pk_bin_list.size() == preamble_len - 1) {
            return ii - (long)std::llround((double)(pk_bin_list.back() - 1) / zero_padding_ratio * 2.0);
        }
        Peak pk0 = dechirp(ii);
        if (!pk_bin_list.empty()) {
            long bin_diff = imod(pk_bin_list.back() - pk0.idx, bin_num);
            if (bin_diff > bin_num / 2) bin_diff = bin_num - bin_diff;
            if (bin_diff <= zero_padding_ratio) pk_bin_list.push_back(pk0.idx);
            else { pk_bin_list.clear(); pk_bin_list.push_back(pk0.idx); }
        } else {
            pk_bin_list.push_back(pk0.idx);
        }
        ii += sample_num;
    }
    return -1;
}

// sync
long LoRaPHY::sync(long x) {
    bool found = false;
    while (x < (long)sig_.size() - sample_num) {
        Peak up_peak   = dechirp(x, true);
        Peak down_peak = dechirp(x, false);
        if (std::abs(down_peak.h) > std::abs(up_peak.h)) found = true;
        x += sample_num;
        if (found) break;
    }
    if (!found) return x;

    Peak pkd = dechirp(x, false);
    long to = (pkd.idx > bin_num / 2)
            ? (long)std::llround((double)(pkd.idx - 1 - bin_num) / zero_padding_ratio)
            : (long)std::llround((double)(pkd.idx - 1) / zero_padding_ratio);
    x = x + to;

    Peak pku = dechirp(x - 4 * sample_num, true);
    preamble_bin = pku.idx;
    cfo_val = (preamble_bin > bin_num / 2)
            ? (double)(preamble_bin - bin_num - 1) * bw / bin_num
            : (double)(preamble_bin - 1) * bw / bin_num;

    Peak pku2 = dechirp(x - sample_num, true);
    Peak pkd2 = dechirp(x - sample_num, false);
    if (std::abs(pku2.h) > std::abs(pkd2.h)) return x + (long)std::llround(2.25 * sample_num);
    else                                     return x + (long)std::llround(1.25 * sample_num);
}

// demodulate
DemodResult LoRaPHY::demodulate(const cvec& sig_in) {
    cfo_val = 0.0;
    init();

    cvec sig = sig_in;
    if (!fast_mode) sig = lowpass(sig, bw / 2.0, fs);

    long p = (long)std::llround(2.0 * bw);
    long q = (long)std::llround(fs);
    sig_ = resample(sig, p, q);

    const double N = std::pow(2.0, sf);
    DemodResult res;
    long x = 1;
    while (x < (long)sig_.size()) {
        x = detect(x);
        if (x < 0) break;
        x = sync(x);

        Peak pk_netid1 = dechirp((long)std::llround(x - 4.25 * sample_num));
        Peak pk_netid2 = dechirp((long)std::llround(x - 3.25 * sample_num));
        double nid1 = dmod(((double)pk_netid1.idx + bin_num - preamble_bin) / zero_padding_ratio, N);
        double nid2 = dmod(((double)pk_netid2.idx + bin_num - preamble_bin) / zero_padding_ratio, N);

        if (x > (long)sig_.size() - 8 * sample_num + 1) break;

        std::vector<double> symbols;
        symbols.reserve(64);
        for (int ii = 0; ii < 8; ++ii) {
            Peak pk = dechirp(x + (long)ii * sample_num);
            symbols.push_back(dmod(((double)pk.idx + bin_num - preamble_bin) / zero_padding_ratio, N));
        }
        if (has_header) {
            if (!parse_header(symbols)) { x += 7 * sample_num; continue; }
        }

        int sym_num = calc_sym_num(payload_len);
        if (x > (long)sig_.size() - (long)sym_num * sample_num + 1) break;
        for (int ii = 8; ii < sym_num; ++ii) {
            Peak pk = dechirp(x + (long)ii * sample_num);
            symbols.push_back(dmod(((double)pk.idx + bin_num - preamble_bin) / zero_padding_ratio, N));
        }
        x += (long)sym_num * sample_num;

        std::vector<double> comp = dynamic_compensation(symbols);
        std::vector<double> col(comp.size());
        for (size_t i = 0; i < comp.size(); ++i) col[i] = dmod(std::round(comp[i]), N);

        res.symbols.push_back(std::move(col));
        res.cfo.push_back(cfo_val);
        res.netid.push_back({ nid1, nid2 });
    }
    return res;
}

// parse_header
bool LoRaPHY::parse_header(const std::vector<double>& data) {
    std::vector<double> comp = dynamic_compensation(data);
    std::vector<int> symbols_g = gray_coding(comp);
    std::vector<int> first8(symbols_g.begin(), symbols_g.begin() + 8);
    std::vector<int> codewords = diag_deinterleave(first8, sf - 2);
    std::vector<int> nibbles = hamming_decode(codewords, 8, hamming_decoding_en);

    payload_len = nibbles[0] * 16 + nibbles[1];
    crc = nibbles[2] & 1;
    cr  = nibbles[2] >> 1;

    int recv0 = nibbles[3] & 1;
    int n5 = nibbles[4];
    int recv[5] = { recv0, (n5 >> 3) & 1, (n5 >> 2) & 1, (n5 >> 1) & 1, n5 & 1 };
    std::array<int, 5> calc = header_checksum(nibbles[0], nibbles[1], nibbles[2]);
    for (int i = 0; i < 5; ++i) if (calc[i] != recv[i]) return false;
    return true;
}

// calc_sym_num
int LoRaPHY::calc_sym_num(int plen) {
    double num = (double)(2 * plen - sf + 7 + 4 * crc - 5 * (1 - has_header));
    double val = std::ceil(num / (double)(sf - 2 * ldr));
    double m = std::max((4.0 + cr) * val, 0.0);
    return (int)(8 + m);
}

// dynamic_compensation
std::vector<double> LoRaPHY::dynamic_compensation(const std::vector<double>& data) {
    const double N = std::pow(2.0, sf);
    size_t n = data.size();
    std::vector<double> symbols(n);
    for (size_t i = 0; i < n; ++i) {
        double drift = (1.0 + (double)(i + 1)) * N * cfo_val / rf_freq;
        symbols[i] = dmod(data[i] - drift, N);
    }
    if (ldr) {
        double bin_offset = 0.0, v_last = 1.0;
        for (size_t i = 0; i < n; ++i) {
            double v = symbols[i];
            double bin_delta = dmod(v - v_last, 4.0);
            if (bin_delta < 2.0) bin_offset -= bin_delta;
            else                 bin_offset += -bin_delta + 4.0;
            v_last = v;
            symbols[i] = dmod(v + bin_offset, N);
        }
    }
    return symbols;
}

// gray_coding
std::vector<int> LoRaPHY::gray_coding(std::vector<double> din) {
    const double N = std::pow(2.0, sf);
    size_t n = din.size();
    for (size_t i = 0; i < n && i < 8; ++i) din[i] = std::floor(din[i] / 4.0);
    for (size_t i = 8; i < n; ++i) {
        if (ldr) din[i] = std::floor(din[i] / 4.0);
        else     din[i] = dmod(din[i] - 1.0, N);
    }
    std::vector<int> symbols(n);
    for (size_t i = 0; i < n; ++i) {
        unsigned s = (unsigned)std::llround(din[i]);
        symbols[i] = (int)(s ^ (s >> 1));
    }
    return symbols;
}

// diag_deinterleave
std::vector<int> LoRaPHY::diag_deinterleave(const std::vector<int>& symbols, int ppm) {
    int n = (int)symbols.size();
    std::vector<int> pre(ppm, 0);
    for (int j = 0; j < ppm; ++j) {
        int val = 0;
        for (int i = 0; i < n; ++i) {
            int jj  = (j + i) % ppm;                      // left-rotate row i by i
            int bit = (symbols[i] >> (ppm - 1 - jj)) & 1; // MSB-first bit
            val |= (bit << i);
        }
        pre[j] = val;
    }
    std::vector<int> codewords(ppm);
    for (int k = 0; k < ppm; ++k) codewords[k] = pre[ppm - 1 - k];   // flipud
    return codewords;
}

// decode
DecodeResult LoRaPHY::decode(const std::vector<std::vector<double>>& symbols_m) {
    DecodeResult res;
    for (const auto& col : symbols_m) {
        std::vector<int> symbols_g = gray_coding(col);
        std::vector<int> first8(symbols_g.begin(), symbols_g.begin() + 8);
        std::vector<int> codewords = diag_deinterleave(first8, sf - 2);

        std::vector<int> nibbles;
        if (!has_header) {
            nibbles = hamming_decode(codewords, 8, hamming_decoding_en);
        } else {
            std::vector<int> nb = hamming_decode(codewords, 8, hamming_decoding_en);
            payload_len = nb[0] * 16 + nb[1];
            crc = nb[2] & 1;
            cr  = nb[2] >> 1;
            int recv0 = nb[3] & 1;
            int n5 = nb[4];
            int recv[5] = { recv0, (n5 >> 3) & 1, (n5 >> 2) & 1, (n5 >> 1) & 1, n5 & 1 };
            std::array<int, 5> calc = header_checksum(nb[0], nb[1], nb[2]);
            for (int i = 0; i < 5; ++i) {
                if (calc[i] != recv[i]) { std::cerr << "[warn] invalid header checksum\n"; break; }
            }
            nibbles.assign(nb.begin() + 5, nb.end());
        }

        int rdd  = cr + 4;
        int ppm2 = sf - 2 * ldr;
        for (int ii = 8; ii + rdd <= (int)symbols_g.size(); ii += rdd) {
            std::vector<int> block(symbols_g.begin() + ii, symbols_g.begin() + ii + rdd);
            std::vector<int> cw = diag_deinterleave(block, ppm2);
            std::vector<int> nb2 = hamming_decode(cw, rdd, hamming_decoding_en);
            nibbles.insert(nibbles.end(), nb2.begin(), nb2.end());
        }

        int nbytes = std::min(255, (int)(nibbles.size() / 2));
        std::vector<uint8_t> bytes(nbytes);
        for (int i = 0; i < nbytes; ++i)
            bytes[i] = (uint8_t)((nibbles[2 * i] & 0xF) | ((nibbles[2 * i + 1] & 0xF) << 4));

        int len = payload_len;
        std::vector<uint8_t> data, checksum;
        if (crc) {
            int pl = std::min(len, (int)bytes.size());
            std::vector<uint8_t> pay(bytes.begin(), bytes.begin() + pl);
            data = whiten(pay);                                  // dewhiten
            if (len + 1 < (int)bytes.size()) {                   
                data.push_back(bytes[len]);
                data.push_back(bytes[len + 1]);
            }
            int cl = std::min(len, (int)data.size());
            std::vector<uint8_t> pc(data.begin(), data.begin() + cl);
            checksum = calc_crc(pc);
        } else {
            int pl = std::min(len, (int)bytes.size());
            std::vector<uint8_t> pay(bytes.begin(), bytes.begin() + pl);
            data = whiten(pay);                                  // dewhiten
        }
        res.data.push_back(std::move(data));
        res.checksum.push_back(std::move(checksum));
    }
    return res;
}

// transmit
std::vector<int> LoRaPHY::gen_header(int plen) {
    std::vector<int> h(5, 0);
    h[0] = plen >> 4;
    h[1] = plen & 15;
    h[2] = (2 * cr) | crc;
    std::array<int, 5> cs = header_checksum(h[0], h[1], h[2]);
    h[3] = cs[0];
    h[4] = 0;
    for (int i = 0; i < 4; ++i) h[4] |= cs[i + 1] << (3 - i);
    return h;
}

std::vector<int> LoRaPHY::diag_interleave(const std::vector<int>& codewords, int rdd) {
    int ncw = (int)codewords.size();
    std::vector<int> symbols_i(rdd);
    for (int r = 0; r < rdd; ++r) {
        int val = 0;
        for (int i0 = 0; i0 < ncw; ++i0) {
            int idx = (i0 + r) % ncw;
            int bit = (codewords[idx] >> r) & 1;
            val |= (bit << i0);
        }
        symbols_i[r] = val;
    }
    return symbols_i;
}

std::vector<int> LoRaPHY::gray_decoding(const std::vector<int>& symbols_i) {
    const double N = std::pow(2.0, sf);
    std::vector<int> symbols(symbols_i.size());
    for (size_t i = 0; i < symbols_i.size(); ++i) {
        unsigned num = (unsigned)symbols_i[i];
        unsigned mask = num >> 1;
        while (mask) { num ^= mask; mask >>= 1; }
        if ((int)i <= 7 || ldr) symbols[i] = (int)dmod((double)num * 4.0 + 1.0, N);
        else                    symbols[i] = (int)dmod((double)num + 1.0, N);
    }
    return symbols;
}

std::vector<int> LoRaPHY::encode(const std::vector<uint8_t>& payload) {
    std::vector<uint8_t> data = payload;
    if (crc) { auto cc = calc_crc(payload); data.push_back(cc[0]); data.push_back(cc[1]); }

    int plen = (int)payload.size();
    int sym_num = calc_sym_num(plen);
    int nibble_num = (sf - 2) + (sym_num - 8) / (cr + 4) * (sf - 2 * ldr);

    int pad = (int)std::ceil((nibble_num - 2 * (int)data.size()) / 2.0);
    if (pad < 0) pad = 0;
    std::vector<uint8_t> data_w = data;
    for (int i = 0; i < pad; ++i) data_w.push_back(255);

    {   // whiten only the first plen bytes
        std::vector<uint8_t> head(data_w.begin(), data_w.begin() + plen);
        std::vector<uint8_t> wh = whiten(head);
        for (int i = 0; i < plen; ++i) data_w[i] = wh[i];
    }

    std::vector<int> data_nibbles(nibble_num);
    for (int i = 1; i <= nibble_num; ++i) {
        int idx = (int)std::ceil(i / 2.0);
        if (i % 2 == 1) data_nibbles[i - 1] = data_w[idx - 1] & 0xF;
        else            data_nibbles[i - 1] = data_w[idx - 1] >> 4;
    }

    std::vector<int> all_nibbles;
    if (has_header) { auto hn = gen_header(plen); all_nibbles = hn; }
    all_nibbles.insert(all_nibbles.end(), data_nibbles.begin(), data_nibbles.end());

    std::vector<int> codewords = hamming_encode(all_nibbles, sf, cr);

    std::vector<int> first(codewords.begin(), codewords.begin() + (sf - 2));
    std::vector<int> symbols_i = diag_interleave(first, 8);
    int ppm = sf - 2 * ldr;
    int rdd = cr + 4;
    for (int i = sf - 1; i + ppm - 1 <= (int)codewords.size(); i += ppm) {
        std::vector<int> block(codewords.begin() + (i - 1), codewords.begin() + (i - 1) + ppm);
        std::vector<int> si = diag_interleave(block, rdd);
        symbols_i.insert(symbols_i.end(), si.begin(), si.end());
    }
    return gray_decoding(symbols_i);
}

cvec LoRaPHY::modulate(const std::vector<int>& symbols) {
    cvec uc = chirp(true,  sf, bw, fs, 0.0, cfo_val, 0.0);
    cvec dc = chirp(false, sf, bw, fs, 0.0, cfo_val, 0.0);
    long chirp_len = (long)uc.size();

    cvec out;
    for (int i = 0; i < preamble_len; ++i) out.insert(out.end(), uc.begin(), uc.end());
    cvec nid1 = chirp(true, sf, bw, fs, 24.0, cfo_val, 0.0);
    cvec nid2 = chirp(true, sf, bw, fs, 32.0, cfo_val, 0.0);
    out.insert(out.end(), nid1.begin(), nid1.end());
    out.insert(out.end(), nid2.begin(), nid2.end());

    out.insert(out.end(), dc.begin(), dc.end());
    out.insert(out.end(), dc.begin(), dc.end());
    long qd = (long)std::llround(chirp_len / 4.0);
    out.insert(out.end(), dc.begin(), dc.begin() + qd);

    for (int s : symbols) {
        cvec ch = chirp(true, sf, bw, fs, (double)s, cfo_val, 0.0);
        out.insert(out.end(), ch.begin(), ch.end());
    }
    return out;
}

// file I/Q
cvec LoRaPHY::read_pcm_int16(const std::string& filename) {
    std::ifstream f = open_binary_in(filename);
    if (!f) return cvec();
    f.seekg(0, std::ios::end);
    std::streamoff sz = f.tellg();
    f.seekg(0);
    size_t nsamp = (size_t)(sz / 2);
    std::vector<int16_t> raw(nsamp);
    if (nsamp) f.read((char*)raw.data(), (std::streamsize)(nsamp * 2));
    size_t npairs = nsamp / 2;
    cvec out(npairs);
    for (size_t i = 0; i < npairs; ++i)
        out[i] = cd((double)raw[2 * i] / 32768.0, (double)raw[2 * i + 1] / 32768.0);
    return out;
}

cvec LoRaPHY::read_cfile_float(const std::string& filename) {
    std::ifstream f = open_binary_in(filename);
    if (!f) return cvec();
    f.seekg(0, std::ios::end);
    std::streamoff sz = f.tellg();
    f.seekg(0);
    size_t nfloat = (size_t)(sz / 4);
    std::vector<float> raw(nfloat);
    if (nfloat) f.read((char*)raw.data(), (std::streamsize)(nfloat * 4));
    size_t npairs = nfloat / 2;
    cvec out(npairs);
    for (size_t i = 0; i < npairs; ++i)
        out[i] = cd((double)raw[2 * i], (double)raw[2 * i + 1]);
    return out;
}

bool LoRaPHY::write_cfile_float(const cvec& data, const std::string& filename) {
    std::ofstream f = open_binary_out(filename);
    if (!f) return false;
    std::vector<float> raw(data.size() * 2);
    for (size_t i = 0; i < data.size(); ++i) {
        raw[2 * i]     = (float)data[i].real();
        raw[2 * i + 1] = (float)data[i].imag();
    }
    f.write((char*)raw.data(), (std::streamsize)(raw.size() * 4));
    return true;
}

bool LoRaPHY::write_pcm_int16(const cvec& data, const std::string& filename) {
    std::ofstream f = open_binary_out(filename);
    if (!f) return false;
    std::vector<int16_t> raw(data.size() * 2);
    for (size_t i = 0; i < data.size(); ++i) {
        double re = std::max(-1.0, std::min(1.0, data[i].real()));
        double im = std::max(-1.0, std::min(1.0, data[i].imag()));
        raw[2 * i]     = (int16_t)std::llround(re * 32767.0);
        raw[2 * i + 1] = (int16_t)std::llround(im * 32767.0);
    }
    f.write((char*)raw.data(), (std::streamsize)(raw.size() * 2));
    return true;
}

double LoRaPHY::time_on_air(int plen) {
    int sym_num = calc_sym_num(plen);
    return (sym_num + 4.25 + preamble_len) * (std::pow(2.0, sf) / bw) * 1000.0;
}

}
