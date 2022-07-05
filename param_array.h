// Copyright 2020-2022 Jon Ville

// Array / buffer based parameter tools

#pragma once

#include "param.h"

namespace jon_dsp {

template <typename VecType>
struct BufferWrapperIOMono {
    typedef typename VecType::elem_t elem_t;
    typedef VecType vec_t;
    
    elem_t *const io_;

    BufferWrapperIOMono(elem_t *const io) : io_{io} {}

    VecType get(const std::size_t i) const {
        return VecType{io_[i]};
    }
    void set(const std::size_t i, const VecType x) {
        io_[i] = x.template get<0>();
    }
};

template <typename VecType>
struct BufferWrapperIOStereo {
    typedef typename VecType::elem_t elem_t;
    typedef VecType vec_t;

    elem_t *const left_io_, *right_io_;

    BufferWrapperIOStereo(elem_t *const left_io, elem_t *const right_io) :
        left_io_{left_io}, right_io_{right_io} {}

    VecType get(const std::size_t i) const {
        return VecType::set_duo(right_io_[i], left_io_[i]);
    }
    void set(const std::size_t i, const VecType x) {
        left_io_[i] = x.template get<0>();
        right_io_[i] = x.template get<1>();
    }
};

template <typename VecType>
struct BufferWrapperVec {
    VecType *const data_;
    typedef typename VecType::elem_t elem_t;
    typedef VecType vec_t;

    BufferWrapperVec(VecType *const data) : data_{data} {}

    VecType get(const std::size_t i) const { return data_[i]; }
    void set(const std::size_t i, const VecType x) { data_[i] = x; }
};

template <typename VecType, typename StoreType>
struct BufferWrapperVecConvert {
    StoreType *const data_;
    typedef typename VecType::elem_t elem_t;
    typedef VecType vec_t;

    BufferWrapperVecConvert(StoreType *const data) : data_{data} {}

    VecType get(const std::size_t i) const {
        return data_[i].template to<VecType>();
    }
    void set(const std::size_t i, const VecType x) {
        data_[i] = x.template to<StoreType>();
    }
};

struct EnabledSwitch {
    struct ParamGroup : public ParamGroupManager<ParamGroup, bool> {
        ParamValidated<bool> current, target;
    } param_;

    void init() { param_.init(); }
    void reset() { param_.current.set(param_.target.get()); }

    void set_enabled(const bool enabled) { param_.target.set(enabled); }

    bool should_reset_wet() const {
        return !param_.current.get() && param_.target.get();
    }
    bool should_reset_dry() const {
        return param_.current.get() && !param_.target.get();
    }

    bool enabled() const {
        return param_.target.get();
    }
};

template <typename VecType, int32_t MIN_DB>
struct DryWetMixBuf {
    struct ParamGroup : public ParamGroupManager<ParamGroup, VecType> {
        ParamValidated<VecType> wet_lin, wet_lin_prev;
    } param_;

    typedef typename VecType::elem_t elem_t;

    void init() {
        static_assert(MIN_DB < 0, "");
        param_.init();
    }

    void reset() {
        param_.wet_lin_prev.set(param_.wet_lin.get());
    }

    void set_wet(const VecType wet) {
        assert(((VecType {0.0} <= wet) && (wet <= VecType {1.0}))
            .debug_valid_eq(true));
        VecType wet_lin = ((wet == 0.0) || (wet == 1.0))
            .choose(wet, db_to_volt_std(elem_t{MIN_DB} + (-elem_t{MIN_DB} * wet)));
        param_.wet_lin.set(wet_lin);
    }

    void set_wet_pc(const VecType wet_pc) { set_wet(wet_pc * 0.01); }

    bool any_dry_signal() {
        static_assert(VecType::elem_count == 1, "");
        return param_.wet_lin.get().data() != 1.0;
    }
    bool is_100pc_dry() {
        static_assert(VecType::elem_count == 1, "");
        return param_.wet_lin.get().data() == 0.0;
    }

    bool any_wet_signal() {
        static_assert(VecType::elem_count == 1, "");
        return param_.wet_lin.get().data() != 0.0;
    }
    bool is_100pc_wet() {
        static_assert(VecType::elem_count == 1, "");
        return param_.wet_lin.get().data() == 1.0;
    }

    bool should_reset_wet() const {
        static_assert(VecType::elem_count == 1, "");
        return param_.wet_lin_prev.get().data() == 0.0 &&
            param_.wet_lin.get().data() != 0.0;
    }
    bool should_reset_dry() const {
        static_assert(VecType::elem_count == 1, "");
        return param_.wet_lin_prev.get().data() == 1.0 &&
            param_.wet_lin.get().data() != 1.0;
    }
    bool on_target() const {
        return (param_.wet_lin.get() == param_.wet_lin_prev.get())
            .debug_valid_eq(true);
    }

    // Put the result in the wet array, the dry array is untouched
    template <typename ArgType>
    void process_instant(const ArgType *const dry, ArgType *const wet,
        const int32_t n) const
    {
        reset();
        typedef typename ArgType::elem_t elem_t;
        const ArgType wet_lin = param_.wet_lin.get().to<ArgType>();
        for (int32_t i = 0; i < n; ++i) {
            wet[i] *= wet_lin;
            wet[i] += dry[i] * (elem_t{1} - wet_lin);
        }
    }

    // wet is untouched, result goes in dry
    template <typename DryBufType, typename WetBufType>
    void process_fade_linear(DryBufType dry_buf, const WetBufType wet_buf, const int32_t n)
    {
        assert(!on_target());
        typedef typename DryBufType::elem_t elem_t;
        typedef typename DryBufType::vec_t vec_t;
        const vec_t prev_wet = param_.wet_lin_prev.get().to<vec_t>(),
            scale = (param_.wet_lin.get().to<vec_t>() - prev_wet) / static_cast<elem_t>(n);
        for (int32_t i = 0; i < n-1; ++i) {
            // i+1 means we start moving on first sample
            const vec_t wet_lin = (prev_wet + (static_cast<elem_t>(i+1) * scale));
            vec_t result = wet_lin * wet_buf.get(i);
            result += (elem_t{1} - wet_lin) * dry_buf.get(i);
            dry_buf.set(i, result);
        }
        // Set final sample to target
        const vec_t wet_lin = param_.wet_lin.get().to<vec_t>();
        vec_t result = wet_lin * wet_buf.get(n-1);
        result += (elem_t{1} - wet_lin) * dry_buf.get(n-1);
        dry_buf.set(n-1, result);
        reset();
    }
};

// Untested
template <typename VecType>
struct SmoothParamBuf {
    struct ParamGroup : public ParamGroupManager<ParamGroup, VecType> {
        ParamValidated<VecType> current, target;
    } param_;

    void init() { param_.init(); }
    void reset() { param_.current.set(param_.target.get()); }

    void set_target(const VecType target) { param_.target.set(target); }
    VecType get_target() const { return param_.target.get(); }
    VecType get_current() const { return param_.current.get(); }
    bool on_target() const {
        return (param_.target.get() == param_.current.get())
            .debug_valid_eq(true);
    }

    template <typename BufType>
    void process_linear(BufType buf, const int32_t n) {
        assert(!on_target());
        typedef typename BufType::elem_t elem_t;
        typedef typename BufType::vec_t vec_t;
        const vec_t current = param_.current.get().to<vec_t>(),
            scale = (param_.target.get().to<vec_t>() - current) / static_cast<elem_t>(n);
        for (int32_t i = 0; i < n-1; ++i) {
            // i+1 means we start moving on first sample
            buf.set(i, (current + (static_cast<elem_t>(i+1) * scale)));
        }
        // Guarantee on target by final sample in buffer
        buf.set(n-1, param_.target.get().to<vec_t>());
        // Set current = target
        reset();
    }
};

} // namespace jon_dsp
