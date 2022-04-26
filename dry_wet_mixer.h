#pragma once

#include <cassert>

#include "db_conv.h"
#include "param.h"

namespace jon_dsp {

template<typename VecType>
inline VecType dry_wet_mix(const VecType& dry, const VecType& wet,
    const VecType& wet_lin)
{
    assert(((VecType {0.0} <= wet_lin) && (wet_lin <= VecType {1.0}))
        .debug_valid_eq(true));
    return wet_lin*wet + (1.0-wet_lin)*dry;
}

// MIN_DB must be negative
// If wet is set to MIN_DB, then it snaps to -inf dB (mul by zero)
template <typename ParamType, int32_t MIN_DB>
struct DryWetMixer {
    struct SmoothParamGroup :
        public SmoothParamGroupManager<SmoothParamGroup, ParamType>
    {
        SmoothParamValidated<ParamType> wet_lin;
    } param_;

    void init() {
        assert(MIN_DB < 0);
        param_.init();
    }

    void snap() {
        param_.snap();
    }

    void unlock_advance() {
        param_.unlock_advance();
        (void) param_.wet_lin.advance_then_get();
    }

    void unlock_advance(bool& reset_dry, bool& reset_wet) {
        reset_dry = false;
        reset_wet = false;
        auto prev_wet = param_.wet_lin.get_current().data();
        param_.unlock_advance();
        auto current_wet = param_.wet_lin.advance_then_get().data();
        if (prev_wet == 1.0 && current_wet != 1.0) reset_dry = true;
        if (prev_wet == 0.0 && current_wet != 0.0) reset_wet = true;
    }

    static constexpr double min_db_ = static_cast<double>(MIN_DB);

    void set_wet(const ParamType& wet, const int32_t sample_rate) {
        assert(((ParamType {0.0} <= wet) && (wet <= ParamType {1.0}))
            .debug_valid_eq(true));
        ParamType wet_lin = ((wet == 0.0) || (wet == 1.0))
            .choose(wet, db_to_volt_std(min_db_ + (-min_db_ * wet)));
        param_.wet_lin.set(wet_lin, sample_rate);
    }
    void set_wet_pc(const ParamType wet_pc, const int32_t sample_rate) {
        set_wet(wet_pc * 0.01, sample_rate);
    }
    bool this_sample_100pc_wet_and_on_target() {
        return param_.wet_lin.get_current().data() == 1.0
            && param_.wet_lin.get_target().data() == 1.0;
    }
    bool this_sample_has_any_dry_signal() {
        return !this_sample_100pc_wet_and_on_target();
    }
    bool this_sample_100pc_dry_and_on_target() {
        return param_.wet_lin.get_current().data() == 0.0
            && param_.wet_lin.get_target().data() == 0.0;
    }
    bool this_sample_has_any_wet_signal() {
        return !this_sample_100pc_dry_and_on_target();
    }
    ParamType get_target() const {
        return param_.wet_lin.get_target();
    }

    // Stateless, can be called multiple times per sample
    template <typename ArgType>
    ArgType iterate(const ArgType& dry, const ArgType& wet) {
        const ArgType wet_lin = ArgType::from(param_.wet_lin.get_current());
        return dry_wet_mix(dry, wet, wet_lin);
    }
};

} // namespace jon_dsp
