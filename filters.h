#pragma once

#include <cassert>

// These filters use cmath for setting coefficients,
// but SSE for evaluating the filter.
#include <cmath>

#include "db_conv.h"
#include "param.h"

namespace jon_dsp {

static constexpr double pi = 3.141592653589793;
static constexpr double two_pi = 6.283185307179586;
static constexpr double sqrt_0p5 = 0.7071067811865476;

static constexpr double k_hipass_hz_ = 38.0;
static constexpr double k_hipass_q_ = 0.5;
static constexpr double k_hipass_scale_ = 1.004997542994373; // 0.0433 dB

static constexpr double k_shelf_db_ = 4.0;
static constexpr double k_shelf_hz_ = 1500.0;
static constexpr double k_shelf_q_ = sqrt_0p5;

// 10^(gain_db/40) = sqrt(10^(gain_db/20)), for peak and shelf filters
template <typename VecType>
inline VecType sqrt_db_to_v_(const VecType& gain_db) {
    return (gain_db * 5.7564627324851146e-2).std_exp();
}

template <typename VecType>
struct Taps1p {
    struct ParamGroup : public ParamGroupManager<ParamGroup, VecType> {
        ParamValidated<VecType> xnm1, ynm1;
    } param_;

    void init() {
        param_.init();
        reset();
    }

    void reset() {
        param_.xnm1.set(0.0); param_.ynm1.set(0.0);
    }

    VecType xnm1() const { return param_.xnm1.get(); }
    VecType ynm1() const { return param_.ynm1.get(); }

    void set_xnm1(const VecType& a) { param_.xnm1.set(a); }
    void set_ynm1(const VecType& a) { param_.ynm1.set(a); }
};

template <typename VecType>
struct Taps2p {
    struct ParamGroup : public ParamGroupManager<ParamGroup, VecType> {
        ParamValidated<VecType> xnm1, xnm2, ynm1, ynm2;
    } param_;

    void init() { 
        param_.init();
        reset();
    }

    void reset() {
        param_.xnm1.set(0.0); param_.xnm2.set(0.0);
        param_.ynm1.set(0.0); param_.ynm2.set(0.0);
    }

    VecType xnm1() const { return param_.xnm1.get(); }
    VecType xnm2() const { return param_.xnm2.get(); }
    VecType ynm1() const { return param_.ynm1.get(); }
    VecType ynm2() const { return param_.ynm2.get(); }

    void set_xnm1(const VecType& a) { param_.xnm1.set(a); }
    void set_xnm2(const VecType& a) { param_.xnm2.set(a); }
    void set_ynm1(const VecType& a) { param_.ynm1.set(a); }
    void set_ynm2(const VecType& a) { param_.ynm2.set(a); }
};

template <typename VecType>
struct Coeff1p {
    struct ParamGroup : public ParamGroupManager<ParamGroup, VecType> {
        ParamValidated<VecType> b0, b1, a1, c, d;
    } param_;

    void init() { param_.init(); }

    void set_identity() {
        param_.d.set(1.0);
        param_.b0.set(0.0); param_.b1.set(0.0); param_.a1.set(0.0);
        param_.c.set(0.0);
    }

    void set_lpf(const typename VecType::elem_t sample_rate,
        const VecType& corner_freq)
    {
        VecType temp {(-two_pi * corner_freq) / sample_rate};
        temp = temp.std_exp();
        param_.b0.set(1.0 - temp);
        param_.b1.set(0.0);
        param_.a1.set(-temp);
        param_.c.set(1.0); param_.d.set(0.0);
    }

    void set_hpf(const typename VecType::elem_t sample_rate,
        const VecType& corner_freq)
    {
        VecType temp {(-two_pi * corner_freq) / sample_rate};
        temp = temp.std_exp();
        param_.b0.set(0.5 * (1.0 + temp));
        param_.b1.set(-param_.b0.get());
        param_.a1.set(-temp);
        param_.c.set(1.0); param_.d.set(0.0);
    }

    // Delay as a proportion of one sample, regardless of sample rate.
    // Old comment says this is identity filter at sample_delay = 0, but not
    // sure of this?
    // No output with sample_delay = 1
    void set_apf(const VecType& sample_delay) {
        assert(((sample_delay >= 0.0) && (sample_delay < 1.0))
            .debug_valid_eq(true));
        param_.b0.set((1.0 - sample_delay) / (1.0 + sample_delay));
        param_.b1.set(1.0);
        param_.a1.set(param_.b0.get());
        param_.c.set(1.0); param_.d.set(0.0);
    }

    void set_low_shelf(const typename VecType::elem_t sample_rate,
        const VecType& corner_freq,
        const VecType& gain_db)
    {
        const VecType theta_c = (two_pi * corner_freq) / sample_rate;
        const VecType mu = db_to_volt_std(gain_db);
        const VecType beta = 4 / (1.0 + mu);
        const VecType delta = beta * (theta_c * 0.5).std_tan();
        const VecType gamma = (1.0 - delta) / (1.0 + delta);
        param_.b0.set(0.5 * (1.0 - gamma));
        param_.b1.set(param_.b0.get());
        param_.a1.set(-gamma);
        param_.c.set(mu - 1.0);
        param_.d.set(1.0);
    }

    void set_hi_shelf(const typename VecType::elem_t sample_rate,
        const VecType& corner_freq,
        const VecType& gain_db)
    {
        const VecType theta_c = (two_pi * corner_freq) / sample_rate;
        const VecType mu = db_to_volt_std(gain_db);
        const VecType beta = (1.0 + mu) * 0.25;
        const VecType delta = beta * (theta_c * 0.5).std_tan();
        const VecType gamma = (1.0 - delta) / (1.0 + delta);
        param_.b0.set(0.5 * (1.0 + gamma));
        param_.b1.set(-0.5 * (1.0 + gamma));
        param_.a1.set(-gamma);
        param_.c.set(mu - 1.0);
        param_.d.set(1.0);
    }

    template <typename ArgType>
    ArgType iterate(Taps1p<ArgType>& taps, const ArgType& x) const {
        ArgType y = x * ArgType::from(param_.b0.get());
        y += taps.xnm1() * ArgType::from(param_.b1.get());
        y -= taps.ynm1() * ArgType::from(param_.a1.get());
        taps.set_xnm1(x);
        taps.set_ynm1(y);
        return ArgType::from(param_.c.get())*y +
            ArgType::from(param_.d.get())*x;
    }
};

template <typename VecType>
struct Coeff2p {
    struct ParamGroup : public ParamGroupManager<ParamGroup, VecType> {
        ParamValidated<VecType> b0, b1, b2, a1, a2, c, d;
    } param_;

    void init() { param_.init(); }

    void set_identity() {
        param_.d.set(1.0);
        param_.b0.set(0.0); param_.b1.set(0.0); param_.b2.set(0.0);
        param_.a1.set(0.0); param_.a2.set(0.0);
        param_.c.set(0.0);
    }

    // butterworth lpf when q is default value
    void set_lpf(const typename VecType::elem_t sample_rate,
        const VecType& corner_freq,
        const VecType& q = typename VecType::elem_t {sqrt_0p5})
    {
        const VecType w0 = two_pi * corner_freq / sample_rate;
        const VecType sinw0 = w0.std_sin();
        const VecType cosw0 = w0.std_cos();
        const VecType alpha = sinw0 / (2.0 * q);

        // We divide all coeffs here, instead of during evaluation
        const VecType a0 = 1.0 + alpha;
        param_.b0.set(((1.0 - cosw0) / 2.0) / a0);
        param_.b1.set((1.0 - cosw0) / a0);
        param_.b2.set(param_.b0.get()); // b0 already divided by a0
        param_.a1.set((-2.0 * cosw0) / a0);
        param_.a2.set((1.0 - alpha) / a0);
        param_.c.set(1.0);
        param_.d.set(0.0);
    }

    // butterworth hpf
    void set_hpf_cookbook(const typename VecType::elem_t sample_rate,
        const VecType& corner_freq,
        const VecType& q = typename VecType::elem_t {sqrt_0p5})
    {
        const VecType w0 = two_pi * corner_freq / sample_rate;
        const VecType sinw0 = w0.std_sin();
        const VecType cosw0 = w0.std_cos();
        const VecType alpha = sinw0 / (2.0 * q);

        const VecType a0 = 1.0 + alpha;
        param_.b0.set(((1.0 + cosw0) / 2.0)/a0);
        param_.b1.set(-(1.0 + cosw0)/a0);
        param_.b2.set(param_.b0.get()); // b0 already scaled
        param_.a1.set((-2.0 * cosw0)/a0);
        param_.a2.set((1.0 - alpha)/a0);
        param_.c.set(1.0);
        param_.d.set(0.0);
    }
    void set_hpf(const typename VecType::elem_t sample_rate,
        const VecType& corner_freq,
        const VecType& q = typename VecType::elem_t {sqrt_0p5})
    {
        const VecType theta = two_pi * corner_freq / sample_rate;
        const VecType q_inv = 1.0 / q;
        const VecType d2sintheta = q_inv * 0.5 * theta.std_sin();
        const VecType beta = 0.5 * ((1.0 - d2sintheta) / (1.0 + d2sintheta));
        const VecType gamma = (0.5 + beta) * theta.std_cos();
        const VecType sum_half_beta_gamma = 0.5 + beta + gamma;
        param_.b0.set(sum_half_beta_gamma * 0.5);
        param_.b1.set(-sum_half_beta_gamma);
        param_.b2.set(param_.b0.get());
        param_.a1.set(-2.0 * gamma);
        param_.a2.set(2.0 * beta);
        param_.c.set(1.0);
        param_.d.set(0.0);
    }

    // This does NOT null with the official coefficients given for
    // 48kHz sample rate, but is extremely close
    void set_k_hipass(const typename VecType::elem_t sample_rate) {
        set_hpf(sample_rate, k_hipass_hz_, k_hipass_q_);
        // Scale the output to get closer to official filters
        param_.c.set(k_hipass_scale_);
    }

    // This nulls with the official coefficients given for 48kHz sample rate
    void set_k_shelf(const typename VecType::elem_t sample_rate) {
        set_hi_shelf(sample_rate, k_shelf_db_, k_shelf_hz_, k_shelf_q_);
    }

    void set_hi_shelf(const typename VecType::elem_t sample_rate,
        const VecType& gain_db,
        const VecType& corner_freq,
        const VecType& q = typename VecType::elem_t {sqrt_0p5})
    {
        const VecType A = sqrt_db_to_v_<VecType>(gain_db);
        const VecType w0 = two_pi * corner_freq / sample_rate;
        const VecType cosw0 = w0.std_cos();
        const VecType alpha = w0.std_sin() / (2.0*q);
        const VecType two_sqrtA_alpha = 2.0 * A.std_sqrt() * alpha;

        const VecType a0 = (A+1) - (A-1)*cosw0 + two_sqrtA_alpha;
        param_.b0.set((A*((A+1) + (A-1)*cosw0 + two_sqrtA_alpha))  /a0);
        param_.b1.set((-2.0*A*((A-1) + (A+1)*cosw0))  /a0);
        param_.b2.set((A*((A+1) + (A-1)*cosw0 - two_sqrtA_alpha))  /a0);
        param_.a1.set((2.0*((A-1) - (A+1)*cosw0))  /a0);
        param_.a2.set(((A+1) - (A-1)*cosw0 - two_sqrtA_alpha)  /a0);

        param_.c.set(1.0); param_.d.set(0.0);
    }

    template <typename ArgType>
    ArgType iterate(Taps2p<ArgType>& taps, const ArgType x) const {
        ArgType y = x * ArgType::from(param_.b0.get());
        y += taps.xnm1() * ArgType::from(param_.b1.get());
        y += taps.xnm2() * ArgType::from(param_.b2.get());
        y -= taps.ynm1() * ArgType::from(param_.a1.get());
        y -= taps.ynm2() * ArgType::from(param_.a2.get());
        taps.set_xnm2(taps.xnm1());
        taps.set_xnm1(x);
        taps.set_ynm2(taps.ynm1());
        taps.set_ynm1(y);
        return ArgType::from(param_.c.get())*y +
            ArgType::from(param_.d.get())*x;
    }
};

} // namespace jon_dsp
