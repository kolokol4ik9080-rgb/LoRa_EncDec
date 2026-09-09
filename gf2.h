#pragma once
#include <array>

namespace lora {


inline int bitget(unsigned v, int pos) { return (int)((v >> (pos - 1)) & 1u); }

std::array<int, 5> header_checksum(int n1, int n2, int n3);

}
