#pragma once

#include <cmath>
#include "../simd_granodi_math/simd_granodi_math.h"

namespace jon_dsp {

static constexpr double VOLT_TO_DB_SMALLEST_DB = -120.0, // 10e-6 volts
v2db_loge_scale_post = 8.685889638065035,
db2v_expe_scale_pre = 1.1512925464970229e-1;

//
//
// STANDARD LIBRARY (slow but accurate)

template <typename VecType>
inline VecType volt_to_db_std(const VecType& x,
    const typename VecType::elem_t min_db = VOLT_TO_DB_SMALLEST_DB)
{
    // 20*log10(x)
    // std::log() of a negative number gives NaN
    // std::log(0.0) gives -inf
    // negative x as a voltage is a valid db, so we use abs(x)
    return VecType::max_fast(v2db_loge_scale_post * x.abs().std_log(), min_db);
}

template <typename VecType>
inline VecType db_to_volt_std(const VecType& x) {
    return (x * typename VecType::elem_t(db2v_expe_scale_pre)).std_exp();
}

//
//
// CEPHES

template <typename VecType>
inline VecType volt_to_db_cm(const VecType& x,
    const double min_db = VOLT_TO_DB_SMALLEST_DB)
{
    return VecType::max(v2db_loge_scale_post * log_cm(x.abs()), min_db);
}

template <typename VecType>
inline VecType db_to_volt_cm(const VecType& x) {
    return exp_cm(x * db2v_expe_scale_pre);
}

//
//
// CUBIC APPROX (smooth but inaccurate)

static constexpr double v2db_log2_scale_post = 6.020599913279623,
db2v_exp2_scale_pre = 1.660964047443681e-1;

template <typename VecType>
inline VecType volt_to_db_p3(const VecType& x,
    const double min_db = VOLT_TO_DB_SMALLEST_DB)
{
    return VecType::max(v2db_log2_scale_post * log2_p3(x), min_db);
}

template <typename VecType>
inline VecType db_to_volt_p3(const VecType& x) {
    exp2_p3(x * db2v_exp2_scale_pre);
}

} // namespace jon_dsp
