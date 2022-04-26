#pragma once

#include "../../DSP_SHARED/simd_granodi_math.h"

namespace jon_dsp {

using namespace simd_granodi;

// For scaling attack / release times to be the same as many compressors /
// expanders, equal to 1/2.2
// Use on the time, NOT the alpha
static constexpr double conventional_time_scale_ = 0.45454545454545453;

// Utility function for exponential smoothing filters (single pole lowpass
// filters calc from time)
// Gives an alpha of 0.0 for true instant attack / release if time_ms = 0.0
template <typename VecType>
inline VecType exp_smoothing_alpha_p3(const double sample_rate,
    const VecType& time_ms)
{
    VecType denom = time_ms * (sample_rate * 0.001);
    auto dgtz = denom > 0.0;
    denom = dgtz.choose(denom, 1.0);
    VecType result = exp_p3(-1.0 / denom);
    return dgtz.choose_else_zero(result);
}

} // namespace jon_dsp
