// whitening.h - LoRa data whitening sequence (LFSR x^8+x^6+x^5+x^4+1).
#pragma once
#include <vector>
#include <cstdint>

namespace lora {

// The 255-byte whitening sequence used by LoRa (see LoRaPHY.m).
extern const uint8_t WHITENING_SEQ[255];

// XOR the first data.size() bytes with the whitening sequence.
// Whitening is its own inverse, so this serves both whiten and dewhiten.
std::vector<uint8_t> whiten(const std::vector<uint8_t>& data);

}
