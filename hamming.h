// hamming.h - LoRa Hamming (4/5 .. 4/8) encoding and decoding.
#pragma once
#include <vector>
#include <cstdint>

namespace lora {

std::vector<int> hamming_decode(const std::vector<int>& codewords, int rdd, bool hamming_en = true);

std::vector<int> hamming_encode(const std::vector<int>& nibbles, int sf, int cr);

}
