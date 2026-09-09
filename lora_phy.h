
#pragma once
#include "lora_types.h"
#include "fft.h"
#include <string>
#include <vector>
#include <array>
#include <cstdint>

namespace lora {

struct DemodResult {
    std::vector<std::vector<double>> symbols;
    std::vector<double>              cfo;
    std::vector<std::array<double,2>> netid;
};

// Decoded bytes per packet.
struct DecodeResult {
    std::vector<std::vector<uint8_t>> data;
    std::vector<std::vector<uint8_t>> checksum;
};

class LoRaPHY {
public:
    
    double rf_freq;              // carrier frequency (Hz), used for SFO drift correction
    int    sf;                  // spreading factor 7..12
    double bw;                  // bandwidth (Hz)
    double fs;                  // sampling frequency (Hz)
    int    cr = 4;              // code rate 1..4 -> 4/5..4/8 (from header in explicit mode)
    int    payload_len = 0;     // payload length in bytes (from header in explicit mode)
    int    has_header = 1;      // 1 = explicit header, 0 = implicit header
    int    crc = 1;             // 1 = payload CRC enabled
    int    ldr = 0;             // low data rate optimization (auto-set in init)
    int    preamble_len = 6;    // number of basic preamble up-chirps
    int    zero_padding_ratio = 10;  // FFT zero-padding ratio
    bool   fast_mode = false;   // true skips the pre-resample low-pass filter
    bool   hamming_decoding_en = true;
    bool   is_debug = false;

    LoRaPHY(double rf_freq, int sf, double bw, double fs);
    void init();                

    // receiver
    DemodResult  demodulate(const cvec& sig_in);
    DecodeResult decode(const std::vector<std::vector<double>>& symbols_m);

    // transmitter (for self-test / round-trip)
    std::vector<int> encode(const std::vector<uint8_t>& payload);
    cvec             modulate(const std::vector<int>& symbols);

    // signal file I/O
    static cvec read_pcm_int16(const std::string& filename);
    static cvec read_cfile_float(const std::string& filename);
    static bool write_cfile_float(const cvec& data, const std::string& filename);
    static bool write_pcm_int16(const cvec& data, const std::string& filename);

    double time_on_air(int plen);

    static cvec chirp(bool is_up, int sf, double bw, double fs,
                      double h, double cfo, double tdelta, double tscale = 1.0);

private:
    // derived state 
    long    bin_num = 0;        // FFT bins of interest (2^sf * zpr)
    long    sample_num = 0;     // samples per symbol at 2*bw (2 * 2^sf)
    long    fft_len = 0;        // zero-padded FFT length (2 * bin_num)
    cvec    downchirp, upchirp; // reference chirps at 2*bw sampling
    long    preamble_bin = 0;   // reference bin used to cancel CFO
    double  cfo_val = 0.0;      // current carrier frequency offset estimate
    DftPlan plan;               // FFT plan for length fft_len
    cvec    sig_;               // resampled working signal

    struct Peak { double h; long idx; };   // idx is 1-based (MATLAB convention)

    Peak dechirp(long x, bool is_up = true);
    long detect(long start_idx);
    long sync(long x);
    bool parse_header(const std::vector<double>& data);
    int  calc_sym_num(int plen);
    std::vector<double> dynamic_compensation(const std::vector<double>& data);
    std::vector<int>    gray_coding(std::vector<double> din);
    std::vector<int>    diag_deinterleave(const std::vector<int>& symbols, int ppm);

    std::vector<int> gen_header(int plen);
    std::vector<int> diag_interleave(const std::vector<int>& codewords, int rdd);
    std::vector<int> gray_decoding(const std::vector<int>& symbols_i);
};

}
