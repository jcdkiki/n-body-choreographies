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

typedef struct {
    double min_dt;     // Минимальный шаг
    double max_dt;     // Максимальный шаг
    double tolerance;  // Допустимая ошибка
    double safety;     // Коэффициент безопасности
    double dt_scale;   // Коэффициент увеличения шага
    double current_dt;  // Текущий шаг интегрирования
    double sim_time;    // Накопленное время симуляции
} AdaptiveParams;

extern int N;
extern BodyState *state;
extern BodyData  *data;

void S_Step_AdamsPredictorCorrector(double dt);
void S_Step_Explicit(BodyState* current_state, BodyState* next_state, double dt);
int S_Step_Adaptive(AdaptiveParams* params);
void S_Step_Implicit(double dt);
void S_Step_Explicit(double dt);
void S_Init(const char *method_filename);

#endif
