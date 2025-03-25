#include "logging.hpp"
#include "simulation.hpp"
#include "render.hpp"
#include <cstdio>
#include <string.h>
typedef enum {
    METHOD_EXPLICIT,
    METHOD_IMPLICIT,
    METHOD_ADAMS,
    METHOD_ADAPTIVE
} MethodType;

MethodType current_method = METHOD_ADAPTIVE;
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

void set_method(const char* method_name) {
    if (strcmp(method_name, "explicit") == 0) {
        current_method = METHOD_EXPLICIT;
    } else if (strcmp(method_name, "implicit") == 0) {
        current_method = METHOD_IMPLICIT;
    } else if (strcmp(method_name, "adams") == 0) {
        current_method = METHOD_ADAMS;
    } else if (strcmp(method_name, "adaptive") == 0) {
        current_method = METHOD_ADAPTIVE;
    } else {
        printf("Unknown method: %s. Using adaptive.\n", method_name);
    }
}

unsigned int callback_base(unsigned int interval, void *name) {
    static double t = 0;
    double dt = (double)interval / 1000.f * speed;
    
    if (current_method == METHOD_EXPLICIT) {
        S_Step_Explicit(dt);
        log_position(t);
    } else if (current_method == METHOD_IMPLICIT) {
        S_Step_Implicit(dt);
        log_position(t);
    }
    
    R_Step(dt);
    t += dt;
    return interval;
}

unsigned int callback_adams(unsigned int interval, void* name) {
    static double t = 0;
    double dt = (double)interval / 1000.0 * speed;
    
    S_Step_AdamsPredictorCorrector(dt);
    log_position(t);
    
    R_Step(dt);
    t += dt;
    return interval;
}

unsigned int callback_adaptive(unsigned int interval, void* name) {
    static double t = 0;
    double real_dt = (double)interval / 1000.0 * speed;
    double target_time = params.sim_time + real_dt;
    
    while (params.sim_time < target_time) {
        double remaining = target_time - params.sim_time;
        if (remaining < params.current_dt) params.current_dt = remaining;
        S_Step_Adaptive(&params);
    }
    
    log_position(params.sim_time);
    R_Step(real_dt);
    t += real_dt;
    return interval;
}

int main(int argc, char **argv) {
    if (argc < 3) {
        printf("Usage: %s <method_file> <speed> [method]\n", argv[0]);
        printf("Available methods: explicit, implicit, adams, adaptive\n");
        return 1;
    }

    init_logs();

    if (argc >= 4) {
        set_method(argv[3]);
    }

    S_Init(argv[1]);
    sscanf(argv[2], "%lf", &speed);

    if (current_method == METHOD_ADAMS) {
        R_Init(callback_adams);
    } else if (current_method == METHOD_ADAPTIVE) {
        R_Init(callback_adaptive);
    } else {
        R_Init(callback_base);
    }
    
    while (1) {
        if (!R_GetInput()) break;
        R_Draw();
    }

    close_logs();
    
    R_Deinit();
    return 0;
}
