// Copyright 2020-2022 Jon Ville

#pragma once

#include "../simd_granodi/simd_granodi.h"

#define JON_DSP_DENORMAL_ERROR "Cannot modify floating point environment"

// Note: I have disabled FENV_ACCESS controls in both clang and MSVC++. It was
// causing a huge slowdown in MSVC++. I haven't tested whether that is also
// the case in clang, but I am avoiding using it for now.

#ifdef __clang__
    // We need clang 12 or newer on sse to use this pragma
    //#if __clang_major__ >= 12 && defined SIMD_GRANODI_ARCH_SSE
        //#pragma STDC FENV_ACCESS ON
    //#else
        #define JON_DSP_DENORMAL_NOINLINE __attribute__((noinline, optnone))
    //#endif
#elif defined __GNUC__
// GCC treats our intrinsics / asm as a signal fence, so this is probably
// not needed. But using it in case of future changes
#define JON_DSP_DENORMAL_NOINLINE __attribute__((noinline, optimize("-O0")))
#elif defined _MSC_VER
//#pragma fenv_access (on)
#define JON_DSP_DENORMAL_NOINLINE __declspec(noinline)
#else
#error JON_DSP_DENORMAL_ERROR
#endif

namespace jon_dsp {

#ifdef SIMD_GRANODI_ARCH_ARM32
typedef uint32_t fp_status;
#elif defined SIMD_GRANODI_ARCH_ARM64
typedef uint64_t fp_status;
#elif defined SIMD_GRANODI_ARCH_SSE
// Same size on x86 and x64
typedef uint32_t fp_status;
#else
#error JON_DSP_DENORMAL_ERROR
#endif

class ScopedDenormalDisable {
    fp_status fp_status_;

    #ifdef JON_DSP_DENORMAL_NOINLINE
    JON_DSP_DENORMAL_NOINLINE
    #endif
    static fp_status get_fp_status() {
        fp_status fs = 0;
        #ifdef SIMD_GRANODI_ARCH_ARM32
        __asm__ volatile("vmrs &0, fpscr" : "=r"(fs));
        #elif defined SIMD_GRANODI_ARCH_ARM64
        __asm__ volatile("mrs %0, fpcr" : "=r"(fs));
        #elif defined SIMD_GRANODI_ARCH_SSE
        fs = _mm_getcsr();
        #endif
        return fs;
    }

    #ifdef JON_DSP_DENORMAL_NOINLINE
    JON_DSP_DENORMAL_NOINLINE
    #endif
    static void set_fp_status(const fp_status fs) {
        #ifdef SIMD_GRANODI_ARCH_ARM32
        __asm__ volatile("vmsr fpscr, %0" : : "ri"(fs));
        #elif defined SIMD_GRANODI_ARCH_ARM64
        __asm__ volatile("msr fpcr, %0" : : "ri"(fs));
        #elif defined SIMD_GRANODI_ARCH_SSE
        _mm_setcsr(fs);
        #endif
    }

    static fp_status disable_denormals() {
        fp_status previous_fs = get_fp_status();
        #if defined SIMD_GRANODI_ARCH_ARM64 || SIMD_GRANODI_ARCH_ARM32
        set_fp_status(previous_fs | 0x1000000);
        #elif defined SIMD_GRANODI_ARCH_SSE
        // flush_to_zero:
        //     "sets denormal results from floating-point calculations to zero"
        // denormals_are_zero:
        //     "treats denormal values used as input to floating-point
        //     instructions as zero"
        // flush_to_zero = 0x8000, denormals_are_zero = 0x40
        set_fp_status(previous_fs | 0x8040);
        #endif
        return previous_fs;
    }

public:
    ScopedDenormalDisable() {
        fp_status_ = disable_denormals();
    }
    ~ScopedDenormalDisable() {
        set_fp_status(fp_status_);
    }
};

} // namespace jon_dsp
