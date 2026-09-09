#include "crc16.h"

namespace lora {

uint16_t crc16_ccitt(const uint8_t* data, size_t n) {
    uint16_t crc = 0x0000;
    for (size_t i = 0; i < n; ++i) {
        crc ^= (uint16_t)data[i] << 8;
        for (int b = 0; b < 8; ++b) {
            if (crc & 0x8000) crc = (uint16_t)((crc << 1) ^ 0x1021);
            else              crc = (uint16_t)(crc << 1);
        }
    }
    return crc;
}

std::vector<uint8_t> calc_crc(const std::vector<uint8_t>& data) {
    size_t len = data.size();
    switch (len) {
        case 0:  return {0, 0};
        case 1:  return {data[0], 0};
        case 2:  return {data[1], data[0]};
        default: {
            uint16_t crc = crc16_ccitt(data.data(), len - 2);
            uint8_t b1 = (uint8_t)((crc & 0xFF)        ^ data[len - 1]);
            uint8_t b2 = (uint8_t)(((crc >> 8) & 0xFF) ^ data[len - 2]);
            return {b1, b2};
        }
    }
}

}
