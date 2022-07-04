// Copyright 2020-2022 Jon Ville

#pragma once

#include "../simd_granodi/simd_granodi.h"

namespace jon_dsp {

template <typename CoeffType>
struct Poly3 {
    CoeffType a_, b_, c_, d_;

    Poly3() : a_ {0.0}, b_ {0.0}, c_ {1.0}, d_ {0.0} {}
    Poly3(const CoeffType a, const CoeffType b,
        const CoeffType c, const CoeffType d)
        : a_ {a}, b_ {b}, c_ {c}, d_ {d} {}

    // Calc coefficients for cube that starts at (x1, y1) with gradient g1 and
    // ends at (x2, y2) with gradient g2
    static Poly3<CoeffType> sg_vectorcall(calc_2p_grad)(const CoeffType x1,
        const CoeffType y1, const CoeffType g1, const CoeffType x2,
        const CoeffType y2, const CoeffType g2)
    {
        CoeffType gscale = x2 - x1;
        CoeffType scale = 1.0 / gscale;
        CoeffType sg1 = gscale * g1;
        CoeffType sg2 = gscale * g2;

        CoeffType a = 2.0*y1 - 2.0*y2 + sg1 + sg2;
        CoeffType b = -3.0*y1 + 3.0*y2 - 2.0*sg1 - sg2;
        CoeffType c = sg1;
        CoeffType d = y1;

        CoeffType offset = -x1;
        CoeffType scale2 = scale * scale;
        CoeffType scale3 = scale2 * scale;
        CoeffType offset2 = offset * offset;
        CoeffType offset3 = offset2 * offset;

        Poly3<CoeffType> result;
        result.a_ = a * scale3;
        result.b_ = 3.0*a*scale3*offset + b*scale2;
        result.c_ = 3.0*a*scale3*offset2 + 2.0*b*scale2*offset + c*scale;
        result.d_ = a*scale3*offset3 + b*scale2*offset2 + c*scale*offset + d;

        return result;
    }

    template <typename VecType>
    VecType sg_vectorcall(eval)(const VecType x) const {
        return x.mul_add(a_.template to<VecType>(), b_.template to<VecType>())
            .mul_add(x, c_.template to<VecType>())
            .mul_add(x, d_.template to<VecType>());
    }
};

template <typename CoeffType>
struct Poly2 {
    CoeffType a_, b_, c_;

    Poly2() : a_ {0.0}, b_ {1.0}, c_ {0.0} {}
    Poly2(const CoeffType a, const CoeffType b,
        const CoeffType c) : a_ {a}, b_ {b}, c_ {c} {}

    // Calc coefficients for quadratic that starts at (x1, y1) with gradient g1
    // and ends at (x2, y2) with gradient g2.
    // Might not work as well as Poly3 for some inputs
    static Poly2<CoeffType> sg_vectorcall(calc_2p_grad)(const CoeffType x1,
        const CoeffType y1, const CoeffType g1,
        const CoeffType x2, const CoeffType y2, const CoeffType g2)
    {
        CoeffType gscale = x2 - x1;
        CoeffType scale = 1.0 / gscale;
        CoeffType sg1 = gscale * g1;
        CoeffType sg2 = gscale * g2;

        CoeffType a = -3.0*y1 + 3.0*y2 - 2.0*sg1 - sg2;
        CoeffType b = sg1;
        CoeffType c = y1;

        CoeffType offset = -x1;
        CoeffType scale2 = scale * scale;
        CoeffType offset2 = offset * offset;

        Poly2<CoeffType> result;
        result.a_ = a*scale2;
        result.b_ = 2.0*a*scale2*offset + b*scale;
        result.c_ = a*scale2*offset2 + b*scale*offset + c;

        return result;
    }

    static Poly2<CoeffType> sg_vectorcall(calc_3p)(const CoeffType x1,
        const CoeffType y1,
        const CoeffType x2, const CoeffType y2,
        const CoeffType x3, const CoeffType y3)
    {
        CoeffType x2_s = x2 * x2;
        CoeffType x3_s = x3 * x3;

        Poly2<CoeffType> result;

        result.b_ = (y2 - y1)*x3_s + (y1 - y2)*x2_s + y2 - y3;
        result.b_ /= (x3_s - x2_s)*(x2 - x1) + x2 - x3;
        result.a_ = y3 - y2 - result.b_*(x3 - x2);
        result.a_ /= x3_s - x2_s;
        result.c_ = y2 - result.a_*x2_s - result.b_*x2;

        return result;
    }

    template <typename VecType>
    VecType sg_vectorcall(eval)(const VecType x) const {
        //return (a_*x + b_)*x + c_;
        return x.mul_add(a_.template to<VecType>(), b_.template to<VecType>())
            .mul_add(x, c_.template to<VecType>());
    }

    template <typename VecType>
    VecType sg_vectorcall(deriv)(const VecType x) const {
        //return 2.0*a_*x + b_;
        return x.mul_add(2.0*a_.template to<VecType>(),
            b_.template to<VecType>());
    }
};

} // namespace jon_dsp
