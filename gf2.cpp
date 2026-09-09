// gf2.cpp - see gf2.h
#include "gf2.h"

namespace lora {

// header_checksum_matrix (5 rows x 12 columns).
static const int H[5][12] = {
    {1,1,1,1,0,0,0,0,0,0,0,0},
    {1,0,0,0,1,1,1,0,0,0,0,1},
    {0,1,0,0,1,0,0,1,1,0,1,0},
    {0,0,1,0,0,1,0,1,0,1,1,1},
    {0,0,0,1,0,0,1,0,1,1,1,1},
};

std::array<int, 5> header_checksum(int n1, int n2, int n3) {
    int bits[12];
    int nib[3] = { n1, n2, n3 };
    for (int k = 0; k < 3; ++k)
        for (int b = 0; b < 4; ++b)
            bits[k * 4 + b] = (nib[k] >> (3 - b)) & 1;

    std::array<int, 5> c{};
    for (int r = 0; r < 5; ++r) {
        int acc = 0;
        for (int j = 0; j < 12; ++j) acc ^= (H[r][j] & bits[j]);
        c[r] = acc;
    }
    return c;
}

}
