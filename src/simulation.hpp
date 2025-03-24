#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include "common.hpp"

constexpr int TRAIL_SIZE = 1024;
constexpr int STEP_MS = 1;

struct BodyData {
    Vec3 trail[TRAIL_SIZE];
    double mass;
};

struct BodyState {
    Vec3 position;
    Vec3 velocity;
};

extern int N;
extern BodyState *state;
extern BodyData  *data;

void S_Step_AdamsPredictorCorrector(double dt);
void S_Step(double dt);
void S_Init(const char *method_filename);

#endif
