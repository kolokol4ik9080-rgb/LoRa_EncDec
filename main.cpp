#include "lora_phy.h"

#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <cmath>
#ifdef _WIN32
#include <windows.h>
#endif

using namespace lora;

// input helpers
static std::string trim(std::string s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

static std::string ask_string(const std::string& p, const std::string& def) {
    std::cout << p << " [" << def << "]: ";
    std::string s;
    if (!std::getline(std::cin, s)) return def;
    s = trim(s);
    return s.empty() ? def : s;
}

static double ask_double(const std::string& p, double def) {
    std::ostringstream os; os << def;
    std::string s = ask_string(p, os.str());
    try { return std::stod(s); } catch (...) { return def; }
}

static long ask_long(const std::string& p, long def) {
    std::ostringstream os; os << def;
    std::string s = ask_string(p, os.str());
    try { return std::stol(s); } catch (...) { return def; }
}

static bool ask_bool(const std::string& p, bool def) {
    std::string s = ask_string(p + (def ? " (Y/n)" : " (y/N)"), def ? "y" : "n");
    if (s.empty()) return def;
    char c = s[0];
    return c == 'y' || c == 'Y' || c == '1';
}

// pretty output
static void print_packet(LoRaPHY& phy, const DemodResult& dm, const DecodeResult& dc, size_t p) {
    std::cout << "\n=== Packet " << (p + 1) << " ===\n";
    std::cout << "CFO   : " << std::fixed << std::setprecision(2) << dm.cfo[p] << " Hz\n";
    std::cout << "NetID : [" << (long)std::llround(dm.netid[p][0]) << ", "
              << (long)std::llround(dm.netid[p][1]) << "]\n";
    std::cout << "Symbols: " << dm.symbols[p].size() << "\n";

    const auto& data = dc.data[p];
    std::vector<uint8_t> payload = data;
    if (phy.crc && data.size() >= 2) payload.assign(data.begin(), data.end() - 2);

    if (phy.crc && data.size() >= 2 && !dc.checksum[p].empty()) {
        int recv = data[data.size() - 2] * 256 + data[data.size() - 1];
        int calc = dc.checksum[p][0] * 256 + dc.checksum[p][1];
        std::cout << "CRC   : " << (calc == recv ? "OK" : "FAIL")
                  << " (calc=" << std::hex << std::uppercase << std::setw(4) << std::setfill('0') << calc
                  << ", recv=" << std::setw(4) << std::setfill('0') << recv
                  << std::dec << std::setfill(' ') << ")\n";
    }

    std::cout << "HEX   : ";
    for (uint8_t b : payload)
        std::cout << std::hex << std::uppercase << std::setw(2) << std::setfill('0') << (int)b << " ";
    std::cout << std::dec << std::setfill(' ') << "\n";

    std::cout << "ASCII : \"";
    for (uint8_t b : payload) std::cout << (char)((b >= 0x20 && b <= 0x7E) ? (char)b : '.');
    std::cout << "\"\n";

    std::cout << "Time on air: " << std::fixed << std::setprecision(2)
              << phy.time_on_air((int)payload.size()) << " ms\n";
}

// run Rx
static void run_rx() {
    std::cout << "\n-- Demodulate & decode from file --\n";
    std::string fname = ask_string("Signal file path", "14-28-27_436886418Hz_01_Fd1M-kv.pcm");
    std::string fmt   = ask_string("File format (pcm=int16 I/Q, cfile=float32 I/Q)", "pcm");
    double fs = ask_double("Sampling rate Fs, Hz", 1e6);
    double bw = ask_double("Bandwidth BW, Hz", 125e3);
    long   sf = ask_long  ("Spreading factor SF (7..12)", 8);
    double rf = ask_double("Center frequency RF, Hz", 438e6);
    bool   hdr  = ask_bool("Explicit header?", true);
    bool   fast = ask_bool("Fast mode (skip low-pass)?", false);
    long   pre  = ask_long("Preamble length (basic up-chirps)", 6);

    LoRaPHY phy(rf, (int)sf, bw, fs);
    phy.has_header   = hdr ? 1 : 0;
    phy.fast_mode    = fast;
    phy.preamble_len = (int)pre;
    if (!hdr) {
        phy.cr          = (int)ask_long("Code rate CR (1=4/5..4=4/8)", 4);
        phy.crc         = ask_bool("CRC enabled?", true) ? 1 : 0;
        phy.payload_len = (int)ask_long("Payload length (bytes)", 0);
    }

    cvec sig = (fmt == "cfile") ? LoRaPHY::read_cfile_float(fname)
                                : LoRaPHY::read_pcm_int16(fname);
    if (sig.empty()) { std::cout << "Error: cannot read or empty file: " << fname << "\n"; return; }
    std::cout << "Loaded samples: " << sig.size() << "  (" << (double)sig.size() / fs << " s)\n";

    std::cout << "Demodulating...\n";
    DemodResult dm = phy.demodulate(sig);
    if (dm.symbols.empty()) { std::cout << "No packets found.\n"; return; }
    std::cout << "Packets found: " << dm.symbols.size() << "\n";

    DecodeResult dc = phy.decode(dm.symbols);
    for (size_t p = 0; p < dc.data.size(); ++p) print_packet(phy, dm, dc, p);
}

// run self-test
static void run_selftest() {
    std::cout << "\n-- Self-test round-trip (encode -> modulate -> demodulate -> decode) --\n";
    double fs = ask_double("Fs, Hz", 1e6);
    double bw = ask_double("BW, Hz", 125e3);
    long   sf = ask_long  ("SF (7..12)", 8);
    double rf = ask_double("RF, Hz", 470e6);
    long   cr = ask_long  ("CR (1..4)", 4);
    bool   crc = ask_bool ("CRC?", true);
    std::string msg = ask_string("Message (payload text)", "Hello LoRa 123");

    LoRaPHY tx(rf, (int)sf, bw, fs);
    tx.has_header = 1; tx.cr = (int)cr; tx.crc = crc ? 1 : 0; tx.fast_mode = false; tx.preamble_len = 8;

    std::vector<uint8_t> payload(msg.begin(), msg.end());
    std::vector<int> syms = tx.encode(payload);
    cvec sig = tx.modulate(syms);
    std::cout << "Encoded symbols: " << syms.size() << ", signal samples: " << sig.size() << "\n";

    if (ask_bool("Save generated signal to selftest.pcm?", false))
        LoRaPHY::write_pcm_int16(sig, "selftest.pcm");

    LoRaPHY rx(rf, (int)sf, bw, fs);
    rx.has_header = 1; rx.fast_mode = false; rx.preamble_len = 8;
    DemodResult dm = rx.demodulate(sig);
    if (dm.symbols.empty()) { std::cout << "Demod: no packets found!\n"; return; }
    DecodeResult dc = rx.decode(dm.symbols);

    auto data = dc.data[0];
    std::vector<uint8_t> got = data;
    if (rx.crc && data.size() >= 2) got.assign(data.begin(), data.end() - 2);
    std::string gotstr(got.begin(), got.end());
    std::cout << "Demodulated symbols: " << dm.symbols[0].size() << "\n";
    std::cout << "Decoded payload    : \"" << gotstr << "\"\n";
    std::cout << (gotstr == msg ? ">> ROUND-TRIP OK\n" : ">> ROUND-TRIP FAILED\n");
}

// run CLI
static int run_cli(int argc, char** argv) {
    std::string cmd = argv[1];
    if (cmd != "rx" || argc < 11) {
        std::cerr << "usage: lora rx <file> <pcm|cfile> <fs> <bw> <sf> <rf> <hdr0|1> <fast0|1> <preamble>\n";
        return 2;
    }
    std::string fname = argv[2], fmt = argv[3];
    double fs = std::stod(argv[4]);
    double bw = std::stod(argv[5]);
    int    sf = std::stoi(argv[6]);
    double rf = std::stod(argv[7]);
    int    hdr = std::stoi(argv[8]);
    int    fast = std::stoi(argv[9]);
    int    pre = std::stoi(argv[10]);

    LoRaPHY phy(rf, sf, bw, fs);
    phy.has_header = hdr; phy.fast_mode = (fast != 0); phy.preamble_len = pre;

    cvec sig = (fmt == "cfile") ? LoRaPHY::read_cfile_float(fname)
                                : LoRaPHY::read_pcm_int16(fname);
    if (sig.empty()) { std::cerr << "cannot read file\n"; return 1; }

    DemodResult dm = phy.demodulate(sig);
    DecodeResult dc = phy.decode(dm.symbols);

    std::cout << "num_packets " << dm.symbols.size() << "\n";
    for (size_t p = 0; p < dm.symbols.size(); ++p) {
        std::cout << "packet " << (p + 1)
                  << " cfo " << std::fixed << std::setprecision(6) << dm.cfo[p]
                  << " netid " << (long)std::llround(dm.netid[p][0]) << " " << (long)std::llround(dm.netid[p][1]) << "\n";
        std::cout << "symbols";
        for (double v : dm.symbols[p]) std::cout << " " << (long)std::llround(v);
        std::cout << "\n";
        const auto& data = dc.data[p];
        std::cout << "nbytes " << data.size() << "\n";
        std::cout << "bytes";
        for (uint8_t b : data) std::cout << " " << (int)b;
        std::cout << "\n";
        std::vector<uint8_t> payload = data;
        if (phy.crc && data.size() >= 2) payload.assign(data.begin(), data.end() - 2);
        std::cout << "ascii ";
        for (uint8_t b : payload) std::cout << (char)((b >= 0x20 && b <= 0x7E) ? (char)b : '.');
        std::cout << "\n";
    }
    return 0;
}

// main
int main(int argc, char** argv) {
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
#endif
    if (argc >= 2) return run_cli(argc, argv);

    std::cout << "======================================\n"
                 "  LoRa PHY (C++ port of LoRaPHY)\n"
                 "======================================\n";
    while (true) {
        std::cout << "\nMenu:\n"
                     "  1) Demodulate & decode a signal file (Rx)\n"
                     "  2) Self-test round-trip (Tx -> Rx)\n"
                     "  0) Exit\n"
                     "Choice: ";
        std::string s;
        if (!std::getline(std::cin, s)) break;
        s = trim(s);
        if (s == "1") run_rx();
        else if (s == "2") run_selftest();
        else if (s == "0" || s.empty()) break;
        else std::cout << "Unknown choice.\n";
    }
    return 0;
}
