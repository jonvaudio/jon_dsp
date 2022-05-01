#include <cstdlib>
#include <cstdio>

#include "../../denormal_disable.h"

#define test_assert(cond) do { \
    if (!(cond)) { \
        printf("test failed on line %d.\n", __LINE__); \
        exit(1); \
    } } while(0)

using namespace jon_dsp;

int main() {
    #ifdef JON_DSP_DENORMAL_NOINLINE
    typedef volatile float denormal_test_float;
    typedef volatile double denormal_test_double;
    #else
    typedef float denormal_test_float;
    typedef double denormal_test_double;
    #endif

    denormal_test_float denormal_f = sg_bitcast_u32x1_f32x1(0x007FFFFF);
    denormal_test_double denormal_d =
        sg_bitcast_u64x1_f64x1(0x000FFFFFFFFFFFFF);

    { ScopedDenormalDisable sdd;
        test_assert(denormal_f == 0.0f);
        test_assert(denormal_d == 0.0);
    }

    test_assert(denormal_f != 0.0f);
    test_assert(denormal_d != 0.0);

    printf("test success\n");

    return 0;
}
