#pragma once
#include "lora_types.h"

namespace lora {

long gcd_ll(long a, long b);

cvec resample(const cvec& x, long p, long q);


cvec lowpass(const cvec& x, double fpass, double fs);

}
