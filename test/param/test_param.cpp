#include <cstdio>
#include <string>

#include "../../param.h"

struct TopLevelEffect : jon_dsp::TopLevelEffectManager<TopLevelEffect> {
    // In practice, these DSP parameters will be inside sub-effects
    PARAM_GROUP_SD(param_, instant);
    SMOOTH_PARAM_GROUP_SD(smooth_param_, smooth);

    jon_dsp::AtomicParam<double> ap_instant, ap_smooth;

    /*jon_dsp::AtomicParamGroup related_group;
    jon_dsp::AtomicParamGroupMember<double> ap_a{related_group},
        ap_b{related_group}, apc{related_group};*/

    // Uncomment below to test all macros
    /*PARAM_GROUP_SS(param_ss_, a, b, c);
    PARAM_GROUP_PS(param_ps_, a, b, c);
    PARAM_GROUP_SD(param_sd_, a, b, c);
    PARAM_GROUP_PD(param_pd_, a, b, c);
    SMOOTH_PARAM_GROUP_SS(smooth_param_ss_, a, b, c);
    SMOOTH_PARAM_GROUP_PS(smooth_param_ps_, a, b, c);
    SMOOTH_PARAM_GROUP_SD(smooth_param_sd_, a, b, c);
    SMOOTH_PARAM_GROUP_PD(smooth_param_pd_, a, b, c);
    PARAM_GROUP_SI32(param_si32_, a, b, c);
    PARAM_GROUP_PI32(param_pi32_, a, b, c);
    PARAM_GROUP_SI64(param_si64_, a, b, c);
    PARAM_GROUP_PI64(param_pi64_, a, b, c);*/

    void set_instant_from_gui(const double instant) {
        ap_instant.store(instant);
    }

    void set_smooth_from_gui(const double smooth) {
        ap_smooth.store(smooth);
    }

    void init_() {
        param_.init(); smooth_param_.init();

        /*param_ss_.init(); param_ps_.init();
        param_sd_.init(); param_pd_.init();
        smooth_param_ss_.init(); smooth_param_ps_.init();
        smooth_param_sd_.init(); smooth_param_pd_.init();
        param_si32_.init(); param_pi32_.init();
        param_si64_.init(); param_pi64_.init();*/
    }

    jon_dsp::Vec_pd iterate_(const jon_dsp::Vec_pd& sample) {
        // todo: uncomment this when can construct vec from scalar type
        return sample + param_.instant.get().set1() +
            smooth_param_.smooth.advance_then_get().set1();
    }
    jon_dsp::Vec_pd iterate_(const jon_dsp::Vec_pd& sample,
        const jon_dsp::Vec_pd& sidechain)
    {
        (void) sidechain;
        return iterate_(sample);
    }

    void snap_() {
        smooth_param_.snap();

        /*smooth_param_ss_.snap(); smooth_param_ps_.snap();
        smooth_param_sd_.snap(); smooth_param_pd_.snap();*/
    }

    void unlock_advance_() {
        smooth_param_.unlock_advance();

        /*smooth_param_ss_.unlock_advance(); smooth_param_ps_.unlock_advance();
        smooth_param_sd_.unlock_advance(); smooth_param_pd_.unlock_advance();*/
    }

    void read_atomic_params_() {
        if (ap_instant.consume()) {
            param_.instant.set(ap_instant.load());
        }

        if (ap_smooth.consume()) {
            smooth_param_.smooth.set(ap_smooth.load(), sample_rate());
        }
    }

    void publish_meter_() {}
};

int main() {
    TopLevelEffect effect;
    static constexpr int32_t blocksize = 103;
    static constexpr int32_t process_block_calls = 3;
    double left_io[blocksize*process_block_calls],
        right_io[blocksize*process_block_calls];

    #ifdef NDEBUG
    printf("NDEBUG is defined. This is an optimized build.\n");
    #else
    printf("NDEBUG is NOT defined. This is a DEBUG build.\n");
    #endif

    #ifdef JON_DSP_VALIDATE_PARAMS
    printf("JON_DSP_VALIDATE_PARAMS is defined.\n");
    #else
    printf("JON_DSP_VALIDATE_PARAMS is NOT defined.\n");
    #endif

    for (int32_t i = 0; i < blocksize*process_block_calls; ++i) {
        left_io[i] = right_io[i] = 0.0;
    }

    effect.init(48000);

    // the set from gui methods can happen before or after a top level effect
    // is initialised, they are only read when the effect starts to process
    // data
    effect.set_instant_from_gui(0.0);
    effect.set_smooth_from_gui(0.0);

    effect.process<double>(left_io, right_io, nullptr, nullptr, blocksize);

    effect.set_instant_from_gui(5.0);
    effect.set_smooth_from_gui(6.0);

    effect.process<double>(&left_io[blocksize], &right_io[blocksize],
        nullptr, nullptr, blocksize);

    effect.set_smooth_from_gui(12.0);
    effect.process<double>(&left_io[blocksize*2], &right_io[blocksize*2],
        nullptr, nullptr, blocksize);


    #ifdef JON_DSP_VALIDATE_PARAMS
    std::string filename = "test_validated_output";
    #else
    std::string filename = "test_opt_output";
    #endif
    printf("The output filename is %s\n\n", filename.data());

    FILE *output = fopen((filename + ".csv").data(), "w");
    for (int32_t i = 0; i < blocksize*process_block_calls; ++i) {
        fprintf(output, "%.9f,%.9f\n", left_io[i], right_io[i]);
    }
    fclose(output);
}
