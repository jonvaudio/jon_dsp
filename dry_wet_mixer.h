// Copyright 2020-2022 Jon Ville

#pragma once

#include <cassert>

#include "db_conv.h"
#include "param.h"

namespace jon_dsp {

template<typename VecType>
inline VecType sg_vectorcall(dry_wet_mix)(const VecType dry, const VecType wet,
    const VecType wet_lin)
{
    assert(((VecType {0.0} <= wet_lin) && (wet_lin <= VecType {1.0}))
        .debug_valid_eq(true));
    return wet_lin*wet + (1.0-wet_lin)*dry;
}

// MIN_DB must be negative
// If wet is set to MIN_DB, then it snaps to -inf dB (mul by zero)
template <typename ParamType, int32_t MIN_DB>
struct DryWetMixer_Instant {
    struct ParamGroup : public ParamGroupManager<ParamGroup, ParamType> {
        ParamValidated<ParamType> wet_lin, wet_lin_prev;
    } param_;

    static constexpr double min_db_ = static_cast<double>(MIN_DB);

    void init() {
        static_assert(MIN_DB < 0, "");
        param_.init();
    }

    void reset() {
        param_.wet_lin_prev.set(param_.wet_lin.get());
    }

    void set_wet(const ParamType wet) {
        assert(((ParamType {0.0} <= wet) && (wet <= ParamType {1.0}))
            .debug_valid_eq(true));
        ParamType wet_lin = ((wet == 0.0) || (wet == 1.0))
            .choose(wet, db_to_volt_std(min_db_ + (-min_db_ * wet)));
        param_.wet_lin.set(wet_lin);
    }

    void set_wet_pc(const ParamType wet_pc) { set_wet(wet_pc * 0.01); }

    bool any_dry_signal() {
        static_assert(ParamType::elem_count == 1, "");
        return param_.wet_lin.get().data() != 1.0;
    }

    bool any_wet_signal() {
        static_assert(ParamType::elem_count == 1, "");
        return param_.wet_lin.get().data() != 0.0;
    }

    // Call BEFORE iterate()
    bool sg_vectorcall(should_reset_wet_before_iterate)() const {
        static_assert(ParamType::elem_count == 1, "");
        return param_.wet_lin_prev.get().data() == 0.0 &&
            param_.wet_lin.get().data() != 0.0;
    }

    // Call BEFORE iterate()
    bool sg_vectorcall(should_reset_dry_before_iterate)() const {
        static_assert(ParamType::elem_count == 1, "");
        return param_.wet_lin_prev.get().data() == 1.0 &&
            param_.wet_lin.get().data() != 1.0;
    }

    template <typename VecType>
    VecType iterate(const VecType dry, const VecType wet) {
        reset();
        const VecType wet_lin = param_.wet_lin.get()
            .template to<VecType>();
        return dry_wet_mix(dry, wet, wet_lin);
    }
};

// This one is smoothed
template <typename ParamType, int32_t MIN_DB>
struct DryWetMixer {
    struct SmoothParamGroup :
        public SmoothParamGroupManager<SmoothParamGroup, ParamType>
    {
        SmoothParamValidated<ParamType> wet_lin;
    } param_;

    void sg_vectorcall(init)() {
        static_assert(MIN_DB < 0, "");
        param_.init();
    }

    void sg_vectorcall(reset)() {
        param_.reset();
    }

    static constexpr double min_db_ = static_cast<double>(MIN_DB);

    void sg_vectorcall(set_wet)(const ParamType wet, const float sample_rate) {
        assert(((ParamType {0.0} <= wet) && (wet <= ParamType {1.0}))
            .debug_valid_eq(true));
        ParamType wet_lin = ((wet == 0.0) || (wet == 1.0))
            .choose(wet, db_to_volt_std(min_db_ + (-min_db_ * wet)));
        param_.wet_lin.set(wet_lin, sample_rate);
    }
    void sg_vectorcall(set_wet_pc)(const ParamType wet_pc,
        const float sample_rate)
    {
        set_wet(wet_pc * 0.01, sample_rate);
    }

    // Call AFTER advance()
    bool sg_vectorcall(any_dry_signal)() {
        static_assert(ParamType::elem_count == 1, "");
        return get_current().data() != 1.0;
    }

    // Call AFTER advance()
    bool sg_vectorcall(any_wet_signal)() {
        static_assert(ParamType::elem_count == 1, "");
        return get_current().data() != 0.0;
    }

    ParamType sg_vectorcall(get_current)() const {
        return param_.wet_lin.get_current();
    }
    ParamType sg_vectorcall(get_target)() const {
        return param_.wet_lin.get_target();
    }

    // Call BEFORE advance()
    bool sg_vectorcall(should_reset_wet_before_advance)() const {
        static_assert(ParamType::elem_count == 1, "");
        return get_current().data() == 0.0 && get_target().data() != 0.0;
    }

    // Call BEFORE advance()
    bool sg_vectorcall(should_reset_dry_before_advance)() const {
        static_assert(ParamType::elem_count == 1, "");
        return get_current().data() == 1.0 && get_target().data() != 1.0;
    }

    void sg_vectorcall(advance)() {
        param_.unlock_advance();
        (void) param_.wet_lin.advance_then_get();
    }

    template <typename ArgType>
    ArgType sg_vectorcall(advance_then_iterate)(const ArgType dry,
        const ArgType wet)
    {
        param_.unlock_advance();
        const ArgType wet_lin = param_.wet_lin.advance_then_get()
            .template to<ArgType>();
        return dry_wet_mix(dry, wet, wet_lin);
    }

    // This iterate is stateless and can be called multiple times per input
    // sample
    template <typename ArgType>
    ArgType sg_vectorcall(iterate_stateless)(const ArgType dry,
        const ArgType wet)
    {
        const ArgType wet_lin = param_.wet_lin.get_current()
            .template to<ArgType>();
        return dry_wet_mix(dry, wet, wet_lin);
    }
};

} // namespace jon_dsp
