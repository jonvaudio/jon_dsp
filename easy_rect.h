// Copyright 2020-2022 Jon Ville

#pragma once

#include <cassert>

#ifdef JON_DSP_JUCE
#include <JuceHeader.h>
#endif

namespace jon_dsp {

// Immutable constexpr 2D rectangle implementation designed to make grid-based
// layouts as easy as possible. Can export a juce rectangle if
// JON_DSP_JUCE is defined.
// (0.0, 0.0) is top left corner
template <typename NumType>
struct EasyRect {
    // These are const, but making them so deletes copy and move constructors
    NumType x_, y_, w_, h_;
    constexpr EasyRect() : x_(0), y_(0), w_(0), h_(0) {}
    constexpr EasyRect(const NumType x, const NumType y, const NumType w,
        const NumType h)
        : x_(x), y_(y), w_(w), h_(h) {}
    #ifdef JON_DSP_JUCE
    constexpr EasyRect(const juce::Rectangle<NumType>& jrect)
        : x_(jrect.getX()), y_(jrect.getY()), w_(jrect.getWidth()),
        h_(jrect.getHeight()) {}
    #endif

    constexpr NumType x() const { return x_; }
    constexpr NumType y() const { return y_; }
    constexpr NumType w() const { return w_; }
    constexpr NumType h() const { return h_; }

    #ifdef JON_DSP_JUCE
    template <typename DestNumType>
    juce::Rectangle<DestNumType> to_juce() const {
        return juce::Rectangle<DestNumType>(static_cast<DestNumType>(x_),
            static_cast<DestNumType>(y_),
            static_cast<DestNumType>(w_), static_cast<DestNumType>(h_));
    }
    template <typename SrcNumType>
    static EasyRect from_juce(const juce::Rectangle<SrcNumType>& jr) {
        return EasyRect(static_cast<NumType>(jr.getX()),
            static_cast<NumType>(jr.getY()),
            static_cast<NumType>(jr.getWidth()),
            static_cast<NumType>(jr.getHeight()));
    }
    #endif

    template <typename DestNumType>
    constexpr EasyRect to() const {
        EasyRect<DestNumType> result(static_cast<DestNumType>(x_),
            static_cast<DestNumType>(y_),
            static_cast<DestNumType>(w_), static_cast<DestNumType>(h_));
        return result;
    }

    // Read as "increase top by...", "decrease bottom by..."
    // Only the incr/decr edge moves
    constexpr EasyRect inc_top(const NumType units) const {
        return EasyRect(x_, y_-units, w_, h_+units);
    }
    constexpr EasyRect dec_top(const NumType units) const {
        return inc_top(-units);
    }

    constexpr EasyRect inc_left(const NumType units) const {
        return EasyRect(x_-units, y_, w_+units, h_);
    }
    constexpr EasyRect dec_left(const NumType units) const {
        return inc_left(-units);
    }

    constexpr EasyRect inc_right(const NumType units) const {
        return EasyRect(x_, y_, w_+units, h_);
    }
    constexpr EasyRect dec_right(const NumType units) const {
        return inc_right(-units);
    }

    constexpr EasyRect inc_bottom(const NumType units) const {
        return EasyRect(x_, y_, w_, h_+units);
    }
    constexpr EasyRect dec_bottom(const NumType units) const {
        return inc_bottom(-units);
    }

    // Grow and shrink each edge keeping centre origin
    constexpr EasyRect grow(const NumType units) const {
        return EasyRect(x_-units, y_-units, w_+2*units, h_+2*units);
    }
    constexpr EasyRect shrink(const NumType units) const {
        return grow(-units);
    }

    // take a new rectangle consisting of ... units from the ...
    // the ctor asserts units > 0.0
    constexpr EasyRect take_top(const NumType units) const {
        return EasyRect(x_, y_, w_, units);
    }
    constexpr EasyRect take_left(const NumType units) const {
        return EasyRect(x_, y_, units, h_);
    }
    constexpr EasyRect take_right(const NumType units) const {
        return EasyRect(x_+w_-units, y_, units, h_);
    }
    constexpr EasyRect take_bottom(const NumType units) const {
        return EasyRect(x_, y_+h_-units, w_, units);
    }

    constexpr void split_horiz(const NumType units, EasyRect& left,
        EasyRect& right) const
    {
        EasyRect safe_copy = *this; // In case user passes a reference to *this
        left = safe_copy.take_left(units); right = safe_copy.dec_left(units);
    }
    constexpr void split_from_left(const NumType units, EasyRect& left,
        EasyRect& right) const
    {
        return split_horiz(units, left, right);
    }
    constexpr void split_from_right(const NumType units, EasyRect& left,
        EasyRect& right) const
    {
        EasyRect safe_copy = *this;
        left = safe_copy.dec_right(units); right = safe_copy.take_right(units);
    }

    constexpr void split_vert(const NumType units, EasyRect& top,
        EasyRect& bottom) const
    {
        EasyRect safe_copy = *this;
        top = safe_copy.take_top(units); bottom = safe_copy.dec_top(units);
    }
    constexpr void split_from_top(const NumType units, EasyRect& top,
        EasyRect& bottom) const
    {
        return split_vert(units, top, bottom);
    }
    constexpr void split_from_bottom(const NumType units, EasyRect& top,
        EasyRect& bottom) const
    {
        EasyRect safe_copy = *this;
        top = safe_copy.dec_bottom(units); bottom = safe_copy.take_bottom(units);
    }

    // Easily move around within grids
    constexpr EasyRect translate(const NumType x, const NumType y) const {
        return EasyRect(x_+x, y_+y, w_, h_);
    }
    constexpr EasyRect translate_right(const NumType units) const {
        return translate(units, 0);
    }
    constexpr EasyRect translate_left(const NumType units) const {
        return translate(-units, 0);
    }
    constexpr EasyRect translate_up(const NumType units) const {
        return translate(0, -units);
    }
    constexpr EasyRect translate_down(const NumType units) const {
        return translate(0, units);
    }
    constexpr EasyRect grid_right(const NumType pad=0) const {
        return translate(w_+pad, 0);
    }
    constexpr EasyRect grid_left(const NumType pad=0) const {
        return translate(-w_-pad, 0);
    }
    constexpr EasyRect grid_up(const NumType pad=0) const {
        return translate(0, -h_-pad);
    }
    constexpr EasyRect grid_down(const NumType pad=0) const {
        return translate(0, h_+pad);
    }

    constexpr EasyRect place_right(const EasyRect& rect,
        const NumType pad=0) const
    {
        return EasyRect(rect.x()+rect.w()+pad, rect.y(), w_, h_);
    }
    constexpr EasyRect place_left(const EasyRect& rect,
        const NumType pad=0) const
    {
        return EasyRect(rect.x()-w_-pad, rect.y(), w_, h_);
    }
    constexpr EasyRect place_above(const EasyRect& rect,
        const NumType pad=0) const
    {
        return EasyRect(rect.x(), rect.y()-y_-pad, w_, h_);
    }
    constexpr EasyRect place_below(const EasyRect& rect,
        const NumType pad=0) const
    {
        return EasyRect(rect.x(), rect.y()+rect.h()+pad, w_, h_);
    }


    // Build grids
    constexpr EasyRect grid_container(const EasyRect& cell,
        const NumType num_cols, const NumType num_rows)
    {
        return EasyRect(cell.x(), cell.y(), cell.x()*num_cols,
            cell.y().num_rows);
    }
    constexpr EasyRect row_container(const EasyRect& cell,
        const NumType num_cols)
    {
        return grid_container(cell, num_cols, 1);
    }
    constexpr EasyRect col_container(const EasyRect& cell,
        const NumType num_rows)
    {
        return grid_container(cell, 1, num_rows);
    }

    // contain_pad(): given a content rectangle C and some padding, generate a
    // container rectangle with the same x and y coordinates as C
    constexpr EasyRect contain_pad(const NumType pad_left,
        const NumType pad_right, const NumType pad_top,
        const NumType pad_bottom)
    {
        return EasyRect(x_, y_, pad_left+pad_right, pad_top+pad_bottom);
    }
    // content_pad(): given a content rectangle C and some padding, move C so
    // that it can be padded. Its size will not change
    constexpr EasyRect content_pad(const NumType pad_left,
        const NumType pad_top)
    {
        return translate(pad_left, pad_top);
    }
};

} // namespace jon_dsp
