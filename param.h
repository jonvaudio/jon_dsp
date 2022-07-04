// Copyright 2020-2022 Jon Ville

#pragma once

#include <algorithm>
#include <atomic>
#include <cassert>
#include <cstring> // memcpy
#include <initializer_list>
#include "../simd_granodi/simd_granodi.h"
#include "denormal_disable.h"

#ifndef NDEBUG
#define JON_DSP_VALIDATE_PARAMS
#endif

namespace jon_dsp {

using namespace simd_granodi;

// Simplest possible wrapper for an array that only does bounds checking in a
// debug build. Absolutely none of the functionality of std::array.
template <typename T, std::size_t SIZE>
class RTArray {
    T data_[SIZE];

    void assert_i_(const std::size_t i) const {
        assert(i < SIZE);
        (void) i;
    }
public:
    RTArray() {}
    RTArray(const std::initializer_list<T>& init_vals) {
        if (init_vals.size() == 1) set_all(*(init_vals.begin()));
        else if (init_vals.size() == SIZE) {
            std::size_t i = 0;
            for (const T x : init_vals) {
                if (i < SIZE) data_[i++] = x;
            }
        } else assert(false); // "Initialization list is the wrong size"
    }

    const T* data() const { return data_; }
    T* data() { return data_; }

    void sg_vectorcall(set_all)(const T x) {
        for (std::size_t i = 0; i < SIZE; ++i) data_[i] = x;
    }

    T sg_vectorcall(operator[])(const std::size_t i) const {
        assert_i_(i); return data_[i];
    }

    T& operator[](const std::size_t i) {
        assert_i_(i); return data_[i];
    }

    static constexpr std::size_t size = SIZE;
};

template <typename VecType>
class ParamValidated {
    VecType current_;
    #ifdef JON_DSP_VALIDATE_PARAMS
    bool valid_;
    #endif
public:
    ParamValidated() { init_(); }

    void sg_vectorcall(init_)() {
        #ifdef JON_DSP_VALIDATE_PARAMS
        valid_ = false;
        #endif
    }

    void sg_vectorcall(set)(const VecType x) {
        #ifdef JON_DSP_VALIDATE_PARAMS
        valid_ = true;
        #endif
        current_ = x;
    }

    VecType sg_vectorcall(get)() const {
        #ifdef JON_DSP_VALIDATE_PARAMS
        assert(valid_);
        #endif
        return current_;
    }
};

template<typename VecType, int32_t SMOOTH_TIME_MS = 1>
class SmoothParamValidated {
    VecType current_, delta_, target_;
    static constexpr typename VecType::elem_t smooth_time_ms_ =
        typename VecType::elem_t {SMOOTH_TIME_MS};
    #ifdef JON_DSP_VALIDATE_PARAMS
    bool valid_current_, valid_target_, can_advance_;
    #endif
public:
    SmoothParamValidated() { init_(); }
    void sg_vectorcall(init_)() {
        #ifdef JON_DSP_VALIDATE_PARAMS
        valid_current_ = false; valid_target_ = false; can_advance_ = false;
        #endif
    }

    void sg_vectorcall(set)(const VecType new_target, const float sample_rate) {
        #ifdef JON_DSP_VALIDATE_PARAMS
        assert(sample_rate > 0);
        valid_target_ = true;
        // We don't assert(valid_current_) as we might not have reset at
        // the start yet (see comment from reset_() as to why we don't reset
        // when setting the parameter for the first time)
        #endif
        target_ = new_target;
        using elem_t = typename VecType::elem_t;
        delta_ = (target_ - current_) /
            (smooth_time_ms_ * elem_t(0.001) * elem_t(sample_rate));
    }

    // This is a weird function that was originally called from more than one
    // place, but needs to be cleaned up / moved / integrated into where it is
    // called from
    void sg_vectorcall(calc_next_current_delta_)(VecType& next_current,
        VecType& next_delta) const
    {
        #ifdef JON_DSP_VALIDATE_PARAMS
        assert(valid_current_);
        assert(valid_target_);
        #endif
        VecType next = current_ + delta_;
        // Snap if:
        // - Delta too small (next == current_),
        // - We've overshot ((next > target) == (current_ < target)
        // - We've hit the target exactly (next == target_)
        auto should_snap = (next == current_) ||
            ((next > target_) == (current_ < target_)) ||
            (next == target_);
        next_current = should_snap.choose(target_, next);
        next_delta = (!should_snap).choose_else_zero(delta_);
    }

    VecType sg_vectorcall(advance_then_get)() {
        #ifdef JON_DSP_VALIDATE_PARAMS
        assert(can_advance_);
        can_advance_ = false;
        #endif
        calc_next_current_delta_(current_, delta_);
        return current_;
    }

    VecType sg_vectorcall(get_current)() const {
        #ifdef JON_DSP_VALIDATE_PARAMS
        assert(valid_current_);
        #endif
        return current_;
    }

    VecType sg_vectorcall(get_target)() const {
        #ifdef JON_DSP_VALIDATE_PARAMS
        assert(valid_target_);
        #endif
        return target_;
    }

    // Reset is called just before we start using the parameter, not when we
    // set the parameter for the first time. As there might be several
    // "first times" where the parameter is loaded in from different places
    // (defaults, presets etc) before the audio starts playing.
    void sg_vectorcall(reset_)() {
        #ifdef JON_DSP_VALIDATE_PARAMS
        assert(valid_target_);
        valid_current_ = true;
        #endif
        current_ = target_;
        delta_ = 0;
    }

    void sg_vectorcall(unlock_advance)() {
        #ifdef JON_DSP_VALIDATE_PARAMS
        can_advance_ = true;
        #endif
    }
};

namespace hidden {

// Parameter method caller for a struct of type ParamGroupType whose only
// members are parameters, all of type ParamType
template <typename ParamGroupType, typename ParamType, typename MethodPolicy>
inline void param_group_method_caller_(char *param_group) {
    static constexpr std::size_t param_group_size = sizeof(ParamGroupType),
        param_size = sizeof(ParamType),
        count = param_group_size / param_size;
    static_assert(param_group_size >= param_size, "Wrong way round");
    static_assert((param_group_size % param_size) == 0,
        "Problem with parameter struct padding");
    for (std::size_t i = 0; i < count; ++i, param_group += param_size) {
        ParamType dummy;
        std::memcpy(&dummy, param_group, param_size);
        MethodPolicy::call_method(dummy);
        std::memcpy(param_group, &dummy, param_size);
    }
}

} // namespace hidden

template <typename ParamGroupType, typename VecType>
class ParamGroupManager {
    struct Init_ {
        static void call_method(ParamValidated<VecType>& p) { p.init_(); }
    };
public:
    void init() {
        hidden::param_group_method_caller_<ParamGroupType,
            ParamValidated<VecType>, Init_>(reinterpret_cast<char*>(this));
    }
};

template <typename ParamGroupType, typename VecType>
class SmoothParamGroupManager {
    struct Init_ {
        static void call_method(SmoothParamValidated<VecType>& p) { p.init_(); }
    };
    struct Reset_ {
        static void call_method(SmoothParamValidated<VecType>& p) { p.reset_(); }
    };
    struct UnlockAdvance_ {
        static void call_method(SmoothParamValidated<VecType>& p) {
            p.unlock_advance();
        }
    };
public:
    void init() {
        hidden::param_group_method_caller_<ParamGroupType,
            SmoothParamValidated<VecType>, Init_>(
                    reinterpret_cast<char*>(this));
    }
    void reset() {
        hidden::param_group_method_caller_<ParamGroupType,
            SmoothParamValidated<VecType>, Reset_>(
                reinterpret_cast<char*>(this));
    }
    void unlock_advance() {
        hidden::param_group_method_caller_<ParamGroupType,
            SmoothParamValidated<VecType>, UnlockAdvance_>(
                reinterpret_cast<char*>(this));
    }
};

template <typename T>
inline constexpr T max_constexpr(const T a, const T b) {
    return a > b ? a : b;
}

static constexpr int32_t STANDARD_SAMPLE_RATES[] = {
      8000,
     11025,
     16000,
     22050,
     44100,
     48000,
     88200,
     96000,
    176400,
    192000,
    352800,
    384000
};

static constexpr float STANDARD_SAMPLE_RATES_MIN = 8000.0f,
    STANDARD_SAMPLE_RATES_MAX = 384000.0f;

// Scale some kind of size by common sample rates. Values outside the range of
// common sample rates will be constrained.
inline constexpr int32_t sg_vectorcall(scale_size_by_sample_rate)(
    const int32_t sample_rate, const int32_t size_at_4448,
    const int32_t smallest)
{
    return max_constexpr(smallest,
        sample_rate <= 11025 ? size_at_4448 / 4 :       // 8, 11.025
        (sample_rate <= 22050 ? size_at_4448 / 2 :      // 16, 22.05
        (sample_rate <= 48000 ? size_at_4448 :          // 44.1, 48
        (sample_rate <= 96000 ? size_at_4448 * 2 :      // 88.2, 96
        (sample_rate <= 192000 ? size_at_4448 * 4 :     // 176.4, 192
        size_at_4448 * 8)))));                          // 352.8, 384
}

inline constexpr int32_t sg_vectorcall(max_size_needed)(
    const int32_t size_at_4448)
{
    return scale_size_by_sample_rate(384000, size_at_4448, 1);
}

// If you have some top level effect, then do:
// class MyTopLevelEffect : public TopLevelEffectManager<MyEffect> { ... }
// and implement all the methods it calls
template <typename TopLevelEffectType, int32_t BlockSizeAt4448 = 256>
class TopLevelEffectManager {
    // float can represent any conventional sample rate exactly, and be casted
    // to double without warning
    float sample_rate_{0.0f};
    int32_t timer_size_{0};
    bool should_reset_this_block_{true};
    TopLevelEffectType& effect() {
        return *static_cast<TopLevelEffectType*>(this);
    }

    void sg_vectorcall(assert_ready_)() const { assert(initialised()); }
public:
    static constexpr int32_t BlockSizeAt4448_ = BlockSizeAt4448;

    bool sg_vectorcall(initialised)() const { return sample_rate_ != 0.0f; }
    bool sg_vectorcall(atomic_params_have_been_read)() const {
        assert_ready_();
        return !should_reset_this_block_;
    }

    float sg_vectorcall(sample_rate)() const {
        assert_ready_(); return sample_rate_;
    }

    void sg_vectorcall(init)(const float sample_rate) {
        assert(STANDARD_SAMPLE_RATES_MIN <= sample_rate &&
            sample_rate <= STANDARD_SAMPLE_RATES_MAX);
        // For non-debug builds, we also want to constrain the sample rate
        // to avoid buffer overflows etc
        const float safe_sample_rate =
            std::min(std::max(STANDARD_SAMPLE_RATES_MIN, sample_rate),
                STANDARD_SAMPLE_RATES_MAX);
        sample_rate_ = safe_sample_rate;
        timer_size_ = scale_size_by_sample_rate(
            static_cast<int32_t>(safe_sample_rate), BlockSizeAt4448, 1);
        should_reset_this_block_ = true;

        ScopedDenormalDisable sdd;
        effect().init_();
    }

    void sg_vectorcall(reset)() {
        assert_ready_();
        ScopedDenormalDisable sdd;
        effect().reset_();
    }

    template <typename FloatType>
    void process_sc(FloatType *const left_io,
        FloatType *right_io,
        const FloatType *left_sc,
        const FloatType *right_sc,
        const int32_t block_size)
    {
        // Input checks
        assert_ready_();
        assert(left_io != nullptr);
        if (left_io == nullptr) return;
        right_io = right_io != nullptr ? right_io : left_io;
        left_sc = left_sc != nullptr ? left_sc : left_io;
        // We don't pick left_sc here as
        // left_sc == nullptr && right_sc == nullptr
        // would force a mono sidechain
        right_sc = right_sc != nullptr ? right_sc : right_io;

        // Actual processing
        ScopedDenormalDisable sdd;
        for (int32_t i = 0; i < block_size; i += timer_size_) {
            effect().read_atomic_params_();
            if (should_reset_this_block_) {
                effect().reset_();
                should_reset_this_block_ = false;
            }
            const int32_t limit = std::min(i+timer_size_, block_size);
            effect().iterate_block_io_(left_io+i, right_io+i,
                left_sc+i, right_sc+i, limit-i);
            effect().publish_meter_();
        }
    }

    template <typename FloatType>
    void process(FloatType *const left_io,
        FloatType *right_io,
        const int32_t block_size)
    {
        process_sc<FloatType>(left_io, right_io, nullptr, nullptr, block_size);
    }
};

//
// (Juce's implementation of audio parameter listeners uses a writer lock,
// which guarantees juce params will only be written to by one writer)

// One writer (GUI) does:
//     param.store(new_value);
// One reader (audio thread) does:
//     if(param.consume()) internal_param = param.load();
// All writes are guaranteed to be read, but may be redundantly read more than
// once under some instruction interleavings
template <typename FloatOrDouble>
class AtomicParam {
    std::atomic<FloatOrDouble> atomic_param_ {FloatOrDouble{}};
    // Init to false, so nothing is read until an actual store is received
    std::atomic<bool> flag_ {false};
    #ifdef JON_DSP_VALIDATE_PARAMS
    bool init_first_time_ {false};
    #endif
public:
    void sg_vectorcall(store)() {
        flag_.store(true);
    }

    void sg_vectorcall(store)(const FloatOrDouble new_value) {
        atomic_param_.store(new_value);
        #ifdef JON_DSP_VALIDATE_PARAMS
        init_first_time_ = true;
        #endif
        store();
    }

    bool sg_vectorcall(consume)() {
        bool expected = true;
        return flag_.compare_exchange_strong(expected, false);
    }

    // Has NO effect on flags
    FloatOrDouble sg_vectorcall(load)() const {
        #ifdef JON_DSP_VALIDATE_PARAMS
        assert(init_first_time_);
        #endif
        return atomic_param_.load();
    }
};

class AtomicParamGroup {
    std::atomic<bool> flag_{false};
public:
    void sg_vectorcall(store)() { flag_.store(true); }

    bool sg_vectorcall(consume)() {
        bool expected = true;
        return flag_.compare_exchange_strong(expected, false);
    }
};

// AtomicParam does not need validation, as it won't be read from until
// it's been written to for the first time.
// But AtomicParamGroupMemeber does, as store()ing to another member of the
// group will trigger a read from all members.
// The validation for AtomicParamGroupMember is only for very first call
// of init() due to complications, but this is good enough
template <typename FloatOrDouble>
class AtomicParamGroupMember {
    std::atomic<FloatOrDouble> atomic_param_{FloatOrDouble{}};
    AtomicParamGroup& group_;
    #ifdef JON_DSP_VALIDATE_PARAMS
    bool init_first_time_ {false};
    #endif
public:
    AtomicParamGroupMember(AtomicParamGroup& group)
        : group_{group} {}

    void sg_vectorcall(store)(const FloatOrDouble new_value) {
        atomic_param_.store(new_value);
        #ifdef JON_DSP_VALIDATE_PARAMS
        init_first_time_ = true;
        #endif
        group_.store();
    }

    FloatOrDouble sg_vectorcall(load)() const {
        #ifdef JON_DSP_VALIDATE_PARAMS
        assert(init_first_time_);
        #endif
        return atomic_param_.load();
    }
};

#ifdef JON_DSP_JUCE
// Utility class to use lambdas with juce audio parameters, instead of
// implementing listeners.
struct LambdaizedAPPL : public juce::AudioProcessorParameter::Listener {
    // Warning: changing this to a reference causes crashes
    const std::function<void()> f_;
    LambdaizedAPPL(const std::function<void()>& f) : f_{f} {}

    void parameterValueChanged(int i, float new_val_normalized) override {
        (void) i; (void) new_val_normalized; f_();
    }
    void parameterGestureChanged(int i, bool g) override {
        (void) i; (void) g;
    }
};
#endif

} // namespace jon_dsp

// Optional convenience macros. Note that you cannot declare two of the same
// type of PARAM_GROUP in the same scope.
// These can't be used if the effect uses a template for its parameter types.
#define PARAM_GROUP_SS(NAME,...) struct ParamGroupSS : \
public jon_dsp::ParamGroupManager<ParamGroupSS, jon_dsp::Vec_ss> \
{ jon_dsp::ParamValidated<jon_dsp::Vec_ss> __VA_ARGS__; } NAME

#define PARAM_GROUP_PS(NAME,...) struct ParamGroupPS : \
public jon_dsp::ParamGroupManager<ParamGroupPS, jon_dsp::Vec_ps> \
{ jon_dsp::ParamValidated<jon_dsp::Vec_ps> __VA_ARGS__; } NAME

#define PARAM_GROUP_SD(NAME,...) struct ParamGroupSD : \
public jon_dsp::ParamGroupManager<ParamGroupSD, jon_dsp::Vec_sd> \
{ jon_dsp::ParamValidated<jon_dsp::Vec_sd> __VA_ARGS__; } NAME

#define PARAM_GROUP_PD(NAME,...) struct ParamGroupPD : \
public jon_dsp::ParamGroupManager<ParamGroupPD, jon_dsp::Vec_pd> \
{ jon_dsp::ParamValidated<jon_dsp::Vec_pd> __VA_ARGS__; } NAME

#define SMOOTH_PARAM_GROUP_SS(NAME,...) struct SmoothParamGroupSS : \
public jon_dsp::SmoothParamGroupManager<SmoothParamGroupSS, jon_dsp::Vec_ss> \
{ jon_dsp::SmoothParamValidated<jon_dsp::Vec_ss> __VA_ARGS__; } NAME

#define SMOOTH_PARAM_GROUP_PS(NAME,...) struct SmoothParamGroupPS : \
public jon_dsp::SmoothParamGroupManager<SmoothParamGroupPS, jon_dsp::Vec_ps> \
{ jon_dsp::SmoothParamValidated<jon_dsp::Vec_ps> __VA_ARGS__; } NAME

#define SMOOTH_PARAM_GROUP_SD(NAME,...) struct SmoothParamGroupSD : \
public jon_dsp::SmoothParamGroupManager<SmoothParamGroupSD, jon_dsp::Vec_sd> \
{ jon_dsp::SmoothParamValidated<jon_dsp::Vec_sd> __VA_ARGS__; } NAME

#define SMOOTH_PARAM_GROUP_PD(NAME,...) struct SmoothParamGroupPD : \
public jon_dsp::SmoothParamGroupManager<SmoothParamGroupPD, jon_dsp::Vec_pd> \
{ jon_dsp::SmoothParamValidated<jon_dsp::Vec_pd> __VA_ARGS__; } NAME

#define PARAM_GROUP_SI32(NAME,...) struct ParamGroupSI32 : \
public jon_dsp::ParamGroupManager<ParamGroupSI32, jon_dsp::Vec_s32x1> \
{ jon_dsp::ParamValidated<jon_dsp::Vec_s32x1> __VA_ARGS__; } NAME

#define PARAM_GROUP_PI32(NAME,...) struct ParamGroupPI32 : \
public jon_dsp::ParamGroupManager<ParamGroupPI32, jon_dsp::Vec_pi32> \
{ jon_dsp::ParamValidated<jon_dsp::Vec_pi32> __VA_ARGS__; } NAME

#define PARAM_GROUP_SI64(NAME,...) struct ParamGroupSI64 : \
public jon_dsp::ParamGroupManager<ParamGroupSI64, jon_dsp::Vec_s64x1> \
{ jon_dsp::ParamValidated<jon_dsp::Vec_s64x1> __VA_ARGS__; } NAME

#define PARAM_GROUP_PI64(NAME,...) struct ParamGroupPI64 : \
public jon_dsp::ParamGroupManager<ParamGroupPI64, jon_dsp::Vec_pi64> \
{ jon_dsp::ParamValidated<jon_dsp::Vec_pi64> __VA_ARGS__; } NAME
