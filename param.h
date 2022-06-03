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
// However, it does zero itself in its ctor.
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
            for (const T& x : init_vals) {
                if (i < SIZE) data_[i++] = x;
            }
        } else assert(false); // "Initialization list is the wrong size"
    }

    void set_all(const T& x) {
        for (std::size_t i = 0; i < SIZE; ++i) data_[i] = x;
    }

    const T& operator[](const std::size_t i) const {
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

    void init_() {
        #ifdef JON_DSP_VALIDATE_PARAMS
        valid_ = false;
        #endif
    }

    void set(const VecType& x) {
        #ifdef JON_DSP_VALIDATE_PARAMS
        valid_ = true;
        #endif
        current_ = x;
    }

    VecType get() const {
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
    void init_() {
        #ifdef JON_DSP_VALIDATE_PARAMS
        valid_current_ = false; valid_target_ = false; can_advance_ = false;
        #endif
    }

    void set(const VecType& new_target, const int32_t sample_rate) {
        #ifdef JON_DSP_VALIDATE_PARAMS
        assert(sample_rate > 0);
        valid_target_ = true;
        // We don't assert(valid_current_) as we might not have snapped at
        // the start yet (see comment from snap_() as to why we don't snap
        // when setting the parameter for the first time)
        #endif
        target_ = new_target;
        using elem_t = typename VecType::elem_t;
        delta_ = (target_ - current_) /
            (smooth_time_ms_ * elem_t(0.001) * elem_t(sample_rate));
    }

    void calc_next_current_delta_(VecType& next_current,
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

    VecType advance_then_get() {
        #ifdef JON_DSP_VALIDATE_PARAMS
        assert(can_advance_);
        can_advance_ = false;
        #endif
        calc_next_current_delta_(current_, delta_);
        return current_;
    }

    VecType preview_next_get() const {
        VecType next_get, discard;
        calc_next_current_delta_(next_get, discard);
        return next_get;
    }

    VecType get_current() const {
        #ifdef JON_DSP_VALIDATE_PARAMS
        assert(valid_current_);
        #endif
        return current_;
    }

    VecType get_target() const {
        #ifdef JON_DSP_VALIDATE_PARAMS
        assert(valid_target_);
        #endif
        return target_;
    }

    // Snap is called just before we start using the parameter, not when we
    // set the parameter for the first time. As there might be several
    // "first times" where the parameter is loaded in from different places
    // (defaults, presets etc) before the audio starts playing.
    void snap_() {
        #ifdef JON_DSP_VALIDATE_PARAMS
        assert(valid_target_);
        valid_current_ = true;
        #endif
        current_ = target_;
        delta_ = 0;
    }

    void unlock_advance_() {
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
    struct Snap_ {
        static void call_method(SmoothParamValidated<VecType>& p) { p.snap_(); }
    };
    struct UnlockAdvance_ {
        static void call_method(SmoothParamValidated<VecType>& p) {
            p.unlock_advance_();
        }
    };
public:
    void init() {
        hidden::param_group_method_caller_<ParamGroupType,
            SmoothParamValidated<VecType>, Init_>(
                    reinterpret_cast<char*>(this));
    }
    void snap() {
        hidden::param_group_method_caller_<ParamGroupType,
            SmoothParamValidated<VecType>, Snap_>(
                reinterpret_cast<char*>(this));
    }
    void unlock_advance() {
        hidden::param_group_method_caller_<ParamGroupType,
            SmoothParamValidated<VecType>, UnlockAdvance_>(
                reinterpret_cast<char*>(this));
    }
};

template <typename T>
inline constexpr T max_constexpr(const T& a, const T& b) {
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

static constexpr int32_t STANDARD_SAMPLE_RATES_MIN = STANDARD_SAMPLE_RATES[0],
    STANDARD_SAMPLE_RATES_MAX = 384000;

// Scale some kind of size by common sample rates. Values outside the range of
// common sample rates will be constrained.
inline constexpr int32_t scale_size_by_sample_rate(const int32_t sample_rate,
    const int32_t size_at_4448, const int32_t smallest)
{
    return max_constexpr(smallest,
        sample_rate <= 11025 ? size_at_4448 / 4 :       // 8, 11.025
        (sample_rate <= 22050 ? size_at_4448 / 2 :      // 16, 22.05
        (sample_rate <= 48000 ? size_at_4448 :          // 44.1, 48
        (sample_rate <= 96000 ? size_at_4448 * 2 :      // 88.2, 96
        (sample_rate <= 192000 ? size_at_4448 * 4 :     // 176.4, 192
        size_at_4448 * 8)))));                          // 352.8, 384
}

inline constexpr int32_t max_size_needed(const int32_t size_at_4448) {
    return scale_size_by_sample_rate(384000, size_at_4448, 1);
}

namespace param_impl {

// Block transformer for two arrays of IOFloatType (stereo), and an
// effect that processes VecType
template <typename IOFloatType, typename VecType>
struct BlockTransformer_IO {
    IOFloatType *const left_io_, *const right_io_;
    BlockTransformer_IO(IOFloatType* const left_io,
        IOFloatType* const right_io)
        : left_io_ {left_io}, right_io_ {right_io} {}

    bool is_nullptr() const {
        return !(left_io_ && right_io_);
    }

    template <typename EffectType>
    void iterate_(const int32_t i, EffectType& effect) const {
        using elem_t = typename VecType::elem_t;
        VecType sample = VecType::set_duo(static_cast<elem_t>(right_io_[i]),
            static_cast<elem_t>(left_io_[i]));
        sample = effect.iterate_(sample);
        left_io_[i] = static_cast<IOFloatType>(sample.template get<0>());
        right_io_[i] = static_cast<IOFloatType>(sample.template get<1>());
    }
};

// BlockTransformer_IO with a sidechain input. It's fine for the sidechain to
// be nullptr
template <typename IOFloatType, typename VecType>
struct BlockTransformer_IO_SC {
    IOFloatType *const left_io_, *const right_io_;
    const IOFloatType *const left_sc_, *const right_sc_;
    BlockTransformer_IO_SC(IOFloatType* const left_io,
        IOFloatType* const right_io,
        const IOFloatType *const left_sc,
        const IOFloatType *const right_sc)
        : left_io_ {left_io}, right_io_ {right_io},
        left_sc_ {left_sc == nullptr ? left_io : left_sc},
        right_sc_ {right_sc == nullptr ? right_io : right_sc} {}

    bool is_nullptr() const {
        return !(left_io_ && right_io_);
    }

    template <typename EffectType>
    void iterate_(const int32_t i, EffectType& effect) const {
        using elem_t = typename VecType::elem_t;
        VecType sample = VecType::set_duo(static_cast<elem_t>(right_io_[i]),
            static_cast<elem_t>(left_io_[i]));
        const VecType sidechain = VecType::set_duo(
            static_cast<elem_t>(right_sc_[i]),
            static_cast<elem_t>(left_sc_[i]));
        sample = effect.iterate_(sample, sidechain);
        left_io_[i] = static_cast<IOFloatType>(sample.template get<0>());
        right_io_[i] = static_cast<IOFloatType>(sample.template get<1>());
    }
};

} // namespace param_impl

// Todo: use std::conditional and a type function to make this object more
// generic (accept PD or PS, and with / without a SC)
// If you have some top level effect, then do:
// class MyTopLevelEffect : public TopLevelEffectManager<MyEffect> { ... }
// and implement all the methods it calls
template <typename TopLevelEffectType, int32_t BlockSizeAt4448 = 256>
class TopLevelEffectManager {
    int32_t sample_rate_, timer_size_;
    bool snap_smooth_params_;
    TopLevelEffectType& effect_;

    void assert_ready_() const { assert(initialised()); }
public:
    TopLevelEffectManager() : sample_rate_{0}, timer_size_{0},
        snap_smooth_params_{true},
        effect_{*static_cast<TopLevelEffectType*>(this)} {}

    bool initialised() const { return sample_rate_ != 0; }
    bool atomic_params_have_been_read() const {
        assert_ready_();
        return !snap_smooth_params_;
    }

    int32_t sample_rate() const { assert_ready_(); return sample_rate_; }

    void init(const int32_t sample_rate) {
        assert(STANDARD_SAMPLE_RATES_MIN <= sample_rate &&
            sample_rate <= STANDARD_SAMPLE_RATES_MAX);
        // For non-debug builds, we also want to constrain the sample rate
        // to avoid buffer overflows etc
        const int32_t safe_sample_rate =
            std::min(std::max(STANDARD_SAMPLE_RATES_MIN, sample_rate),
                STANDARD_SAMPLE_RATES_MAX);
        sample_rate_ = safe_sample_rate;
        timer_size_ = scale_size_by_sample_rate(
            safe_sample_rate, BlockSizeAt4448, 1);
        snap_smooth_params_ = true;

        ScopedDenormalDisable sdd;
        effect_.init_();
    }

    template <typename BlockTransformer>
    void process_from_transformer_(BlockTransformer& block,
        const int32_t block_size)
    {
        assert_ready_();
        if (block.is_nullptr()) return;
        ScopedDenormalDisable sdd;
        for (int32_t j = 0; j < block_size; j += timer_size_) {
            effect_.read_atomic_params_();
            if (snap_smooth_params_) {
                effect_.snap_();
                snap_smooth_params_ = false;
            }
            int32_t limit = std::min(j+timer_size_, block_size);
            for (int32_t i = j; i < limit; ++i) {
                effect_.unlock_advance_();
                block.iterate_(i, effect_);
            }
            effect_.publish_meter_();
        }
    }

    template<typename FloatType>
    void process_io(FloatType* const left_io, FloatType* const right_io,
        const int32_t block_size)
    {
        using VecType = typename TopLevelEffectType::VecType;
        process_from_transformer_(
            param_impl::BlockTransformer_IO<FloatType, VecType>(
                left_io, right_io), block_size);
    }

    template<typename FloatType>
    void process_io_sc(FloatType* const left_io,
        FloatType* const right_io,
        const FloatType* const left_sc,
        const FloatType* const right_sc,
        const int32_t block_size)
    {
        using VecType = typename TopLevelEffectType::VecType;
        process_from_transformer_(
               param_impl::BlockTransformer_IO_SC<FloatType, VecType>(
                   left_io, right_io, left_sc, right_sc), block_size);
    }

    /*template <typename FloatOrDouble>
    void process(FloatOrDouble* const left_io,
        FloatOrDouble* const right_io,
        const FloatOrDouble* const left_sc,
        const FloatOrDouble* const right_sc,
        const int32_t block_size)
    {
        assert_ready_();
        if (!(left_io && right_io)) return;
        const bool valid_sc = left_sc && right_sc;
        ScopedDenormalDisable sdd;
        for (int32_t j = 0; j < block_size; j += timer_size_) {
            effect_.read_atomic_params_();
            if (snap_smooth_params_) {
                effect_.snap_();
                snap_smooth_params_ = false;
            }
            int32_t limit = std::min(j+timer_size_, block_size);
            for (int32_t i = j; i < limit; ++i) {
                effect_.unlock_advance_();
                Vec_pd sample = {static_cast<double>(right_io[i]),
                        static_cast<double>(left_io[i])},
                    sidechain; // default ctor init with value 0
                if (valid_sc) {
                    sidechain = {static_cast<double>(right_sc[i]),
                        static_cast<double>(left_sc[i])};
                }
                sample = effect_.iterate_(sample, sidechain);
                left_io[i] = static_cast<FloatOrDouble>(sample.d0());
                right_io[i] = static_cast<FloatOrDouble>(sample.d1());
            }
            effect_.publish_meter_();
        }
    }*/
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
public:
    void store(const FloatOrDouble& new_value) {
        atomic_param_.store(new_value);
        flag_.store(true);
    }

    bool consume() {
        bool expected = true;
        return flag_.compare_exchange_strong(expected, false);
    }

    // Has NO effect on flags
    FloatOrDouble load() const {
        return atomic_param_.load();
    }
};

class AtomicParamGroup {
    std::atomic<bool> flag_{false};
public:
    void store() { flag_.store(true); }

    bool consume() {
        bool expected = true;
        return flag_.compare_exchange_strong(expected, false);
    }
};

template <typename FloatOrDouble>
class AtomicParamGroupMember {
    std::atomic<FloatOrDouble> atomic_param_{FloatOrDouble{}};
    AtomicParamGroup& group_;
public:
    AtomicParamGroupMember(AtomicParamGroup& group)
        : group_{group} {}

    void store(const FloatOrDouble& new_value) {
        atomic_param_.store(new_value);
        group_.store();
    }

    FloatOrDouble load() const {
        return atomic_param_.load();
    }
};

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
