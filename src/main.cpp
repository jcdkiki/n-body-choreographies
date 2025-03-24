#include "simulation.hpp"
#include "render.hpp"
#include <cstdio>

double speed = 1.0;
AdaptiveParams params = {
    .min_dt = 1e-6,
    .max_dt = 1e-2,
    .tolerance = 1e-8,
    .safety = 0.8,
    .dt_scale = 1.2,
    .current_dt = 1e-3,
    .sim_time = 0.0
};

unsigned int callback_base(unsigned int interval, void *name)
{
    // для явного и неявного метода
    double dt = (double)interval / 1000.f * speed;

    S_Step_Explicit(dt);
    R_Step(dt);

    return interval;
}

unsigned int callback_predictor_corrector(unsigned int interval, void* name) {
    double dt = (double)interval / 1000.0 * speed;
    S_Step_AdamsPredictorCorrector(dt);
    R_Step(dt);
    
    return interval;
}

unsigned int callback_adaptive(unsigned int interval, void* name) {
    // для метода с адитивным шагом
    double real_dt = (double)interval / 1000.0 * speed;
    
    double target_time = params.sim_time + real_dt;
    while (params.sim_time < target_time) {
        double remaining = target_time - params.sim_time;
        if (remaining < params.current_dt) {
            params.current_dt = remaining;
        }
        
        if (!S_Step_Adaptive(&params)) {
            continue;
        }
    }
    R_Step(real_dt);
    
    return interval;
}

int main(int argc, char **argv)
{
    if (argc == 3) {
        sscanf(argv[2], "%lf", &speed);
    }

    S_Init(argv[1]);
    R_Init(callback_adaptive);

    while (1) {
        if (!R_GetInput()) break;
        R_Draw();
    }

    R_Deinit();
    //S_Deinit();
}
