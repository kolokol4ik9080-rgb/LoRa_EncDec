// crc16.h - LoRa payload CRC (CRC-16/CCITT, polynomial x^16+x^12+x^5+1).
#pragma once
#include <vector>
#include <cstdint>

namespace lora {

// Raw CRC-16/CCITT over a byte buffer, MSB-first, init 0x0000, no reflection,
// no final XOR (matches comm.CRCGenerator with 'X^16 + X^12 + X^5 + 1').
uint16_t crc16_ccitt(const uint8_t* data, size_t n);

// LoRa payload checksum (port of LoRaPHY.calc_crc). Returns the two checksum
// bytes {b1, b2} appended after the payload.
std::vector<uint8_t> calc_crc(const std::vector<uint8_t>& data);

} // namespace lora
