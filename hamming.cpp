// hamming.cpp - see hamming.h
#include "hamming.h"
#include "gf2.h"

namespace lora {

// XOR of selected (1-based) bits of w.
static int bit_reduce(unsigned w, const int* pos, int n) {
    int b = bitget(w, pos[0]);
    for (int i = 1; i < n; ++i) b ^= bitget(w, pos[i]);
    return b;
}

std::vector<int> hamming_decode(const std::vector<int>& codewords, int rdd, bool hamming_en) {
    std::vector<int> nibbles(codewords.size());

    auto parity_fix = [](int p) -> int {
        switch (p) {
            case 3: return 4;   // 011 -> b3 wrong
            case 5: return 8;   // 101 -> b4 wrong
            case 6: return 1;   // 110 -> b1 wrong
            case 7: return 2;   // 111 -> b2 wrong
            default: return 0;
        }
    };

    for (size_t i = 0; i < codewords.size(); ++i) {
        unsigned cw = (unsigned)codewords[i];
        if (!hamming_en) { nibbles[i] = (int)(cw & 0xF); continue; }

        switch (rdd) {
            case 5:
            case 6:
                nibbles[i] = (int)(cw & 0xF);
                break;
            case 7:
            case 8: {
                static const int P2[] = {7, 4, 2, 1};
                static const int P3[] = {5, 3, 2, 1};
                static const int P5[] = {6, 4, 3, 2};
                int p2 = bit_reduce(cw, P2, 4);
                int p3 = bit_reduce(cw, P3, 4);
                int p5 = bit_reduce(cw, P5, 4);
                int parity = p2 * 4 + p3 * 2 + p5;
                cw ^= (unsigned)parity_fix(parity);
                nibbles[i] = (int)(cw & 0xF);
                break;
            }
            default:
                nibbles[i] = (int)(cw & 0xF);
                break;
        }
    }
    return nibbles;
}

std::vector<int> hamming_encode(const std::vector<int>& nibbles, int sf, int cr) {
    std::vector<int> codewords(nibbles.size());
    for (size_t i = 0; i < nibbles.size(); ++i) {
        unsigned nib = (unsigned)nibbles[i];
        static const int A1[] = {1, 3, 4};
        static const int A2[] = {1, 2, 4};
        static const int A3[] = {1, 2, 3};
        static const int A4[] = {1, 2, 3, 4};
        static const int A5[] = {2, 3, 4};
        int p1 = bit_reduce(nib, A1, 3);
        int p2 = bit_reduce(nib, A2, 3);
        int p3 = bit_reduce(nib, A3, 3);
        int p4 = bit_reduce(nib, A4, 4);
        int p5 = bit_reduce(nib, A5, 3);

        int cr_now = ((int)i <= sf - 3) ? 4 : cr;

        unsigned cw = 0;
        switch (cr_now) {
            case 1:
                cw = ((unsigned)p4 << 4) | nib;
                break;
            case 2:
                cw = ((unsigned)p5 << 5) | ((unsigned)p3 << 4) | nib;
                break;
            case 3:
                cw = ((unsigned)p2 << 6) | ((unsigned)p5 << 5) | ((unsigned)p3 << 4) | nib;
                break;
            case 4:
                cw = ((unsigned)p1 << 7) | ((unsigned)p2 << 6) | ((unsigned)p5 << 5) | ((unsigned)p3 << 4) | nib;
                break;
            default:
                cw = nib;
                break;
        }
        codewords[i] = (int)cw;
    }
    return codewords;
}

}