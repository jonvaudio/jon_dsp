// Copyright 2020-2022 Jon Ville

#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>

#ifdef JON_DSP_JUCE
#include <JuceHeader.h>
#endif

namespace jon_dsp {

template <typename NumType>
inline NumType constrain(const NumType min, const NumType max, const NumType x)
{
    return std::min(max, std::max(min, x));
}

template <typename NumType>
inline NumType reflect(const NumType mirror, const NumType x) {
    return mirror + mirror - x;
}

// Private class used for implementation, use RangeLin instead.
template <typename FloatType>
struct RangeLin_Inc_ {
    FloatType lowerb_, upperb_;
    RangeLin_Inc_() : lowerb_(0.0), upperb_(1.0) {}
    RangeLin_Inc_(const FloatType start, const FloatType end) : lowerb_(start),
        upperb_(end)
    {
        if (end < start) std::swap(lowerb_, upperb_);
        assert(lowerb_ < upperb_);
    }

    FloatType lowerb() const { return lowerb_; }
    FloatType upperb() const { return upperb_; }
    FloatType start() const { return lowerb_; }
    FloatType end() const { return upperb_; }

    FloatType normalise(const FloatType unnorm_val) const {
        // This assertion caused a problem with a range where:
        // lowerb = double{1.01}
        // because float{double{1.01}} < double{1.01}
        // So we now constrain instead of assert
        // assert(lowerb_ <= unnorm_val && unnorm_val <= upperb_);
        constrain(lowerb_, upperb_, unnorm_val);
        FloatType norm = (unnorm_val-lowerb_) / (upperb_-lowerb_);
        // in case of rounding error
        return constrain<FloatType>(0.0, 1.0, norm);
    }

    FloatType unnormalise(const FloatType norm_val) const {
        assert(0.0 <= norm_val && norm_val <= 1.0);
        FloatType unnorm = lowerb_ + norm_val*(upperb_-lowerb_);
        return constrain(lowerb_, upperb_, unnorm);
    }
};

// Remembers whether or not to flip a value in the normal domain, for
// ranges where start > end
template <typename FloatType>
struct FlipNorm {
    bool flip_;
    FlipNorm() : flip_(false) {}
    FlipNorm(const FloatType start, const FloatType end) : flip_(start > end) {}
    FlipNorm(const bool flip) : flip_(flip) {}

    bool flipped() const { return flip_; }

    FloatType flipn(const FloatType norm_val) const {
        assert(0 <= norm_val && norm_val <= 1.0);
        return flip_ ? static_cast<FloatType>(1.0) - norm_val : norm_val;
    }
};

// Linear range, increasing or decreasing
// RangeLin_Inc_ does much of the assertions & constraining
template <typename FloatType>
struct RangeLin {
    RangeLin_Inc_<FloatType> range_;
    FlipNorm<FloatType> fn_;
    RangeLin() {}
    RangeLin(const FloatType start, const FloatType end)
        : range_(start, end), fn_(start, end) {}

    FloatType normalise(const FloatType unnorm_val) const {
        return fn_.flipn(range_.normalise(unnorm_val));
    }

    FloatType unnormalise(const FloatType norm_val) const {
        return range_.unnormalise(fn_.flipn(norm_val));
    }

    FloatType lowerb() const { return range_.lowerb(); }
    FloatType upperb() const { return range_.upperb(); }
    FloatType start() const {
        return fn_.flipped() ? range_.upperb() : range_.lowerb();
    }
    FloatType end() const {
        return fn_.flipped() ? range_.lowerb() : range_.upperb();
    }
};

// Skewed range, increasing or decreasing
// Aims to be compatible with juce skew, but with mirroring OFF
// If you need juce mirroring, use RangSymmSkew with a mirror point in the
// center
template <typename FloatType>
struct RangeSkewed {
    RangeLin_Inc_<FloatType> range_;
    FloatType skew_exp_;
    FlipNorm<FloatType> fn_;

    RangeSkewed() : skew_exp_(1.0) {}
    RangeSkewed(const FloatType start, const FloatType end)
        : range_(start, end), skew_exp_(1.0), fn_(start, end) {}
    RangeSkewed(const FloatType start, const FloatType end,
        const FloatType skew)
        : range_(start, end), skew_exp_(skew), fn_(start, end) {}

    // Tested
    RangeSkewed(const FloatType start, const FloatType end,
        const FloatType map_unnorm_val, const FloatType to_norm_val)
        : range_(start, end), fn_(start, end)
    {
        assert(0.0 <= to_norm_val && to_norm_val <= 1.0);
        assert(range_.lowerb() <= map_unnorm_val &&
            map_unnorm_val <= range_.upperb());
        skew_exp_ = std::log(to_norm_val) /
            std::log(range_.normalise(map_unnorm_val));
    }

    // Skew a value in the range [0, 1] into the range [0, 1]
    FloatType skew_(const FloatType unsk_val) const {
        assert(0.0 <= unsk_val && unsk_val <= 1.0);
        // We assume that a std lib implementation of pow() is always the
        // identity for args of 0.0 or 1.0
        return std::pow(unsk_val, skew_exp_);
    }

    // Unskew a value in the range [0, 1] into the range [0, 1]
    FloatType unskew_(const FloatType skewed_val) const {
        assert(0.0 <= skewed_val && skewed_val <= 1.0);
        // less code method that gives almost-same result as juce with some
        // rounding error
        // return std::pow(skewed_val, 1.0 / skew_exp_);
        // juce method that should give identical result
        return skewed_val > 0.0 ?
            std::exp(std::log(skewed_val) / skew_exp_) : skewed_val;
    }

    FloatType normalise(const FloatType unnorm_val) const {
        return fn_.flipn(skew_(range_.normalise(unnorm_val)));
    }

    FloatType unnormalise(const FloatType norm_val) const {
        return range_.unnormalise(unskew_(fn_.flipn(norm_val)));
    }

    FloatType lowerb() const { return range_.lowerb(); }
    FloatType upperb() const { return range_.upperb(); }
    FloatType start() const {
        return fn_.flipped() ? range_.upperb() : range_.lowerb();
    }
    FloatType end() const {
        return fn_.flipped() ? range_.lowerb() : range_.upperb();
    }
};

// Note that to_partial_norm_val is not normalised over the whole range,
// as this would require numerical methods to solve. Instead, it represents
// the proportion between mirror_point and the upper bound. Eg:
// RangeSymmSkew(-12.0, 12.0, 6.0, 0.5, 0.0) which maps 6.0 to 0.5 about 0.0,
// gives a straight line
template <typename FloatType>
struct RangeSymmSkew {
    RangeSkewed<FloatType> lower_range_;
    RangeLin_Inc_<FloatType> whole_range_;
    FloatType lowerb_, upperb_, mirror_point_;
    FlipNorm<FloatType> fn_;

    RangeSymmSkew(const FloatType start, const FloatType end,
        const FloatType map_unnorm_val, const FloatType to_partial_norm_val,
        const FloatType mirror_point) : fn_(start, end)
    {
        lowerb_ = start, upperb_ = end;
        if (fn_.flipped()) std::swap(lowerb_, upperb_);
        assert(lowerb_ <= mirror_point && mirror_point <= upperb_);
        assert(mirror_point != map_unnorm_val);
        mirror_point_ = mirror_point;
        FloatType map_un = map_unnorm_val;
        FloatType to_pn = to_partial_norm_val;
        if (map_unnorm_val > mirror_point_) {
            map_un = reflect(mirror_point_, map_un);
            to_pn = 1.0 - to_pn;
        }
        FloatType extended_lower_range = std::max(mirror_point_ - lowerb_,
            upperb_ - mirror_point_);
        lower_range_ = RangeSkewed<FloatType>(
            mirror_point_-extended_lower_range, mirror_point_, map_un, to_pn);
        whole_range_ = RangeLin_Inc_<FloatType>(partial_normalise_(lowerb_),
            partial_normalise_(upperb_));
    }

    FloatType partial_normalise_(const FloatType unnorm_val) const {
        assert(lowerb_ <= unnorm_val && unnorm_val <= upperb_);
        bool upper = unnorm_val > mirror_point_;
        FloatType un = unnorm_val;
        if (upper) un = reflect<FloatType>(mirror_point_, un);
        FloatType n = lower_range_.normalise(un);
        if (upper) n = reflect<FloatType>(1.0, n);
        return n;
        // Don't need to check for small fp out-of-bounds error here as this
        // function is used to define whole_range_ bounds in the first place
    }

    FloatType normalise(const FloatType unnorm_val) const {
        FloatType pn = partial_normalise_(unnorm_val);
        return fn_.flipn(whole_range_.normalise(pn));
    }

    FloatType unnormalise(const FloatType norm_val) const {
        FloatType pu = whole_range_.unnormalise(fn_.flipn(norm_val));
        bool upper = pu > 1.0;
        if (upper) pu = reflect<FloatType>(1.0, pu);
        FloatType un = lower_range_.unnormalise(pu);
        if (upper) un = reflect<FloatType>(mirror_point_, un);
        // small fp error mid-way through a cut-off skewed curve
        return constrain(lowerb_, upperb_, un);
    }

    FloatType lowerb() const { return lowerb_; }
    FloatType upperb() const { return upperb_; }
    FloatType start() const { return fn_.flipped() ? upperb_ : lowerb_; }
    FloatType end() const { return fn_.flipped() ? lowerb_ : upperb_; }
};

// Eg, base 10, start 30, end 20000 will give nice Hz plot
// For any base B, all powers of B will be equidistant from each other
// in the normalised domain
// Start and end must both be > 0, but end < start is allowed
template <typename FloatType>
struct RangeLog {
    FloatType base_, loge_base_;
    // range_ is to remeber our upper and lower bounds, as recalc
    // from log_range_ gives fp errors that fail assertions
    RangeLin<FloatType> log_range_, range_;

    RangeLog(const FloatType base, const FloatType start, const FloatType end) {
        assert(base > 0.0);
        assert(start > 0.0 && end > 0.0);
        base_ = base;
        loge_base_ = std::log(base);
        log_range_ = RangeLin<FloatType>(std::log(start) /
            loge_base_, std::log(end) / loge_base_);
        range_ = RangeLin<FloatType>(start, end);
    }
    RangeLog(const FloatType start, const FloatType end)
        : RangeLog(10.0, start, end) {}

    FloatType normalise(const FloatType unnorm_val) const {
        // This doesn't assert that unnorm_val is within range, but
        // prevents log of an invalid value. The range_.normalise() will then
        // have a valid range check assertion in the log domain
        assert(0.0 < unnorm_val);
        FloatType logb = std::log(unnorm_val) / loge_base_;
        return log_range_.normalise(logb);
    }

    FloatType unnormalise(const FloatType norm_val) const {
        // lin_range_.unnormalise() does the norm_val in range check before
        // the pow()
        const FloatType result = std::pow(base_,
            log_range_.unnormalise(norm_val));
        // Constrain needed due to FP error
        return constrain(lowerb(), upperb(), result);
    }

    FloatType lowerb() const { return range_.lowerb(); }
    FloatType upperb() const { return range_.upperb(); }
    FloatType start() const { return range_.start(); }
    FloatType end() const { return range_.end(); }
};

//
//
// SNAPPERS
// Snapping notes:
// - Snappers are expected by juce to constrain their output to legal values!
// - snap_unnorm(unnorm_val) should NOT assert that unnorm_val is in range,
//   because the juce gui code often puts 0.0 through snappers.
// - Snappers interpret all ranges as increasing (but can be used with
//   decreasing ranges, just set the snap as if it is increasing)

// Doesn't snap, but still constrains as juce expects it to
template <typename FloatType>
struct SnapperNone {
    FloatType lowerb_, upperb_;

    SnapperNone(const FloatType lower_bound, const FloatType upper_bound)
        : lowerb_(lower_bound), upperb_(upper_bound)
    {
        assert(lower_bound < upper_bound);
    }

    FloatType snap_unnorm(const FloatType unnorm_val) const {
        return constrain(lowerb_, upperb_, unnorm_val);
    }

    FloatType lowerb() const { return lowerb_; }
    FloatType upperb() const { return upperb_; }
};

// eg 3 sig places = 23.2 or 232000
// WARNING: If the lower bound is 0, this will behave strangely!
template <typename FloatType>
struct SnapperSigFigures {
    FloatType lowerb_, upperb_; int sig_places_;

    SnapperSigFigures(const FloatType lower_bound, const FloatType upper_bound,
        const int sig_places)
        : lowerb_(lower_bound), upperb_(upper_bound), sig_places_(sig_places)
    {
        assert(lowerb_ < upperb_);
    }

    FloatType snap_unnorm(const FloatType unnorm_val) const {
        FloatType un = unnorm_val;
        FloatType snapped;
        if (un != 0.0) {
            FloatType sign = 1;
            if (un < 0.0) {
                un = -un;
                sign = -1.0;
            }
            FloatType lg10 = floor(std::log10(un));
            FloatType scale = std::pow(10,
                lg10+1.0-static_cast<FloatType>(sig_places_));
            un /= scale;
            un = std::round(un);
            snapped = sign * un * scale;
        } else {
            snapped = un;
        }
        return constrain(lowerb_, upperb_, snapped);
    }
    FloatType lowerb() const { return lowerb_; }
    FloatType upperb() const { return upperb_; }
};

// Same as juce
template <typename FloatType>
struct SnapperInterval {
    FloatType lowerb_, upperb_, snap_;
    // A default is ctor is needed to have an array of these, it invariants
    // intentionally to save processing
    SnapperInterval() : lowerb_(0.0), upperb_(0.0), snap_(0.0) {}
    SnapperInterval(const FloatType lowerb, const FloatType upperb,
        const FloatType snap)
        : lowerb_(lowerb), upperb_(upperb), snap_(snap)
    {
        assert(lowerb_ < upperb_);
        // we don't care if a snap value is larger than the interval: that just
        // means it always snaps to the lower bound of that interval
        assert(snap >= 0.0);
    }

    FloatType snap_unnorm(const FloatType unnorm_val) const {
        FloatType snapped = unnorm_val;
        FloatType unnv_offset = snapped - lowerb_;
        snapped = lowerb_ + std::round(unnv_offset / snap_) * snap_;
        // small chance to snap outside interval
        return constrain(lowerb_, upperb_, snapped);
    }
    FloatType lowerb() const { return lowerb_; }
    FloatType upperb() const { return upperb_; }
};

// Multiple intervals across different parts of the range.
// It is acceptable to have areas of the range not covered by any interval.
// Intervals may not overlap, but the upperbound of one may equal the lower
// bound of another.
// They must be passed to the ctor in ascending order for assertions to work.
//
// To construct:
// SnapperIntervals(range.lowerb(), range.upperb(),
//     {{{lb,ub,snap}, {lb,ub,snap}, {...}}})
template <typename FloatType, int NumIntervals>
struct SnapperIntervals {
    FloatType lowerb_, upperb_;
    std::array<SnapperInterval<FloatType>, NumIntervals> snaps_;

    SnapperIntervals() { assert(NumIntervals >= 0); }
    SnapperIntervals(const FloatType lower_bound, const FloatType upper_bound,
        const std::array<SnapperInterval<FloatType>, NumIntervals>& snaps)
        : lowerb_(lower_bound), upperb_(upper_bound), snaps_(snaps)
    {
        assert(NumIntervals >= 0);
        // intervals are not allowed to overlap, and they must increase (to help
        // find errors, we don't care about efficiency here)
        // but they DON'T have to cover the whole range (flexible)
        FloatType lastupperb = lowerb_;
        for (int i = 0; i < NumIntervals; ++i) {
            assert(snaps[i].lowerb() >= lastupperb);
            lastupperb = snaps[i].upperb();
            assert(lastupperb <= upperb_);
        }
    }

    FloatType snap_unnorm(const FloatType unnorm_val) const {
        FloatType snapped = unnorm_val;
        for (int i = 0; i < NumIntervals; ++i) {
            // Use of <= and < is intentional to avoid snapping unnorm val in
            // next interval
            if (snaps_[i].lowerb() <= unnorm_val &&
                unnorm_val < snaps_[i].upperb())
            {
                snapped = snaps_[i].snap_unnorm(unnorm_val);
                // Break for efficiency. Program should behave identically
                // without it:
                break;
            }
        }
        // We still have to constrain here, because there is a chance the
        // unnorm_val was not in any snapping interval (intervals don't have
        // to cover the whole range).
        return constrain(lowerb_, upperb_, snapped);
    }

    FloatType lowerb() const { return lowerb_; }
    FloatType upperb() const { return upperb_; }
};

#ifdef JON_DSP_JUCE

// OUTPUT TO JUCE. The resulting range can be used in
// juce::AudioParameterFloat (pass float type via constructor)
// or with GUI sliders (pass double type via ::setNormalisableRange(),
// then use ::setNumDecimalPlacesToDisplay() to set
// sane display digits)
template <typename FloatType, typename RangeType, typename SnapperType>
inline juce::NormalisableRange<FloatType> get_juce_norm_range(
    const RangeType& range, const SnapperType& snapper)
{
    assert(range.lowerb() == snapper.lowerb() &&
        range.upperb() == snapper.upperb());
    return juce::NormalisableRange<FloatType>(
        static_cast<FloatType>(range.lowerb()),
        static_cast<FloatType>(range.upperb()),
    [range](FloatType start, FloatType end, FloatType norm_val) -> FloatType
    {
        (void) start; (void) end;
        return static_cast<FloatType>(
            range.unnormalise(static_cast<double>(norm_val)));
    },
    [range](FloatType start, FloatType end, FloatType unnorm_val) -> FloatType {
        (void) start; (void) end;
        return static_cast<FloatType>(
            range.normalise(static_cast<double>(unnorm_val)));
    },
    [snapper](FloatType start, FloatType end, FloatType unnorm_val)
        -> FloatType {
         (void) start; (void) end;
         return static_cast<FloatType>(
             snapper.snap_unnorm(static_cast<double>(unnorm_val))); });
}

// Convenience function
template<typename RangeType, typename SnapperType>
inline void get_juce_norm_ranges(const RangeType& range,
    const SnapperType& snapper,
    juce::NormalisableRange<double>& doubleDest,
    juce::NormalisableRange<float>& floatDest)
{
    doubleDest = get_juce_norm_range<double>(range, snapper);
    floatDest = get_juce_norm_range<float>(range, snapper);
}

#endif

} // namespace jon_dsp
