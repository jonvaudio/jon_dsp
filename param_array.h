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

    BufferWrapperIOMono() : io_{nullptr} {}
    BufferWrapperIOMono(elem_t *const io) : io_{io} {}

    VecType get(const std::size_t i) const {
        return VecType{io_[i]};
    }
    void set(const std::size_t i, const VecType x) const {
        io_[i] = x.template get<0>();
    }
};

template <typename VecType>
struct BufferWrapperIOStereo {
    typedef typename VecType::elem_t elem_t;
    typedef VecType vec_t;

    elem_t *const left_io_, *right_io_;

    BufferWrapperIOStereo() : left_io_{nullptr}, right_io_{nullptr} {}
    BufferWrapperIOStereo(elem_t *const left_io, elem_t *const right_io) :
        left_io_{left_io}, right_io_{right_io} {}

    VecType sg_vectorcall(get)(const std::size_t i) const {
        return VecType::set_duo(right_io_[i], left_io_[i]);
    }
    template <typename GetType>
    GetType sg_vectorcall(get_type)(std::size_t i) const {
        return GetType::set_duo(right_io_[i], left_io_[i]);
    }
    void sg_vectorcall(set)(const std::size_t i, const VecType x) const {
        left_io_[i] = x.template get<0>();
        right_io_[i] = x.template get<1>();
    }
    template <typename SetType>
    void sg_vectorcall(set_type)(const std::size_t i, const SetType x) const {
        left_io_[i] = x.template get<0>();
        right_io_[i] = x.template get<1>();
    }
};

template <typename VecType>
struct BufferWrapperVec {
    VecType *const data_;
    typedef typename VecType::elem_t elem_t;
    typedef VecType vec_t;

    BufferWrapperVec() : data_{nullptr} {}
    BufferWrapperVec(VecType *const data) : data_{data} {}

    VecType get(const std::size_t i) const { return data_[i]; }
    void set(const std::size_t i, const VecType x) const { data_[i] = x; }
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
    bool any_wet_signal() const {
        // current = false, target = false -> NO
        // current = false, target = true -> YES
        // current = true, target = false -> YES
        // current = true, target = true -> YES
        return param_.current.get() || param_.target.get();
    }
    bool any_dry_signal() const {
        // current = false, target = false -> YES
        // current = false, target = true -> YES
        // current = true, target = false -> YES
        // current = true, target = true -> NO
        return !(param_.current.get() && param_.target.get());
    }

    bool on_target() const {
        return param_.target.get() == param_.current.get();
    }

    bool enabled() const {
        return param_.target.get();
    }
};

// Stateless, and for one-time use processing a block
// Not to be used as a member
template <typename VecType>
struct LinearFade {
    const VecType start_, scale_;
    LinearFade() = delete;
    LinearFade(const VecType start, const VecType end, const int32_t n) :
        start_{start}, 
        scale_{(end-start) / static_cast<typename VecType::elem_t>(n)} {}
    LinearFade(const EnabledSwitch& es) {
        static_assert(VecType::elem_count == 1);
        assert(!es.on_target());
        start_ = es.enabled() ? 0.0 : 1.0;
        scale_ = 1.0 / static_cast<typename VecType::elem_t>(n);
        if (!es.enabled()) scale_ = -scale_;
    }

    const VecType start() const { return start_; }
    const VecType scale() const { return scale_; }

    static VecType sg_vectorcall(calculate_i)(const VecType start,
        const VecType scale, const int32_t i)
    {
        // Use i+1 so that we start moving on the first sample, and hit the
        // target (although with a small amount of fp error) on the final sample
        return start + (static_cast<typename VecType::elem_t>(i+1) * scale);
    }

    template <typename BufType>
    void process(BufType buf, const int32_t n) {
        const VecType start = start_, scale = scale_;
        for (int32_t i = 0; i < n; ++i) {
            buf.set(i, calculate_i(start, scale, i));
        }
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

    /*template <typename ArgType>
    void process_instant(const ArgType *const dry, ArgType *const wet,
        const int32_t n) const
    {
        reset();
        typedef typename ArgType::elem_t elem_t;
        const ArgType wet_lin = param_.wet_lin.get().template to<ArgType>();
        for (int32_t i = 0; i < n; ++i) {
            wet[i] *= wet_lin;
            wet[i] += dry[i] * (elem_t{1} - wet_lin);
        }
    }*/

    // result goes in dry, wet is untouched
    template <typename DryBufType, typename WetBufType>
    void process_instant(DryBufType dry_buf, const WetBufType wet_buf, const int32_t n)
    {
        assert(on_target());
        typedef typename DryBufType::elem_t elem_t;
        typedef typename DryBufType::vec_t vec_t;
        const vec_t wet_lin = param_.wet_lin.get().template to<vec_t>();
        for (int32_t i = 0; i < n; ++i) {
            vec_t result = wet_lin * wet_buf.get(i).template to<vec_t>();
            result += (elem_t{1} - wet_lin) * dry_buf.get(i);
            dry_buf.set(i, result);
        }
    }

    // result goes in dry, wet is untouched
    template <typename DryBufType, typename WetBufType>
    void process_fade_linear(DryBufType dry_buf, const WetBufType wet_buf, const int32_t n, bool reset_after_fade=true)
    {
        assert(!on_target());
        typedef typename DryBufType::elem_t elem_t;
        typedef typename DryBufType::vec_t vec_t;
        const vec_t prev_wet = param_.wet_lin_prev.get().template to<vec_t>(),
            scale = (param_.wet_lin.get().template to<vec_t>() - prev_wet) / static_cast<elem_t>(n);
        for (int32_t i = 0; i < n-1; ++i) {
            // i+1 means we start moving on first sample
            const vec_t wet_lin = (prev_wet + (static_cast<elem_t>(i+1) * scale));
            vec_t result = wet_lin * wet_buf.get(i).template to<vec_t>();
            result += (elem_t{1} - wet_lin) * dry_buf.get(i);
            dry_buf.set(i, result);
        }
        // Set final sample to target
        const vec_t wet_lin = param_.wet_lin.get().template to<vec_t>();
        vec_t result = wet_lin * wet_buf.get(n-1).template to<vec_t>();
        result += (elem_t{1} - wet_lin) * dry_buf.get(n-1);
        dry_buf.set(n-1, result);
        if (reset_after_fade) reset();
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
    void process_linear(BufType buf, const int32_t n, bool reset_after_fade=true) {
        assert(!on_target());
        typedef typename BufType::elem_t elem_t;
        typedef typename BufType::vec_t vec_t;
        const vec_t current = param_.current.get().template to<vec_t>(),
            scale = (param_.target.get().template to<vec_t>() - current) / static_cast<elem_t>(n);
        for (int32_t i = 0; i < n-1; ++i) {
            // i+1 means we start moving on first sample
            buf.set(i, (current + (static_cast<elem_t>(i+1) * scale)));
        }
        // Guarantee on target by final sample in buffer
        buf.set(n-1, param_.target.get().template to<vec_t>());
        // Set current = target
        if (reset_after_fade) reset();
    }
};

} // namespace jon_dsp
