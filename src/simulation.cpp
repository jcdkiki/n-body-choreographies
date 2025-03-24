#include "simulation.hpp"
#include "common.hpp"
#include <cassert>
#include <cfloat>
#include <string.h>
#include <stdio.h>

struct Method {
    int n;
    double **a;
    double *k;
};

int N;
BodyState *state;
BodyData  *data;

static Method method;

void S_CalculateDelta(BodyState *cur_state, BodyState *delta)
{
    for (int i = 0; i < N; i++) {
        delta[i].position = cur_state[i].velocity;
        delta[i].velocity = { 0, 0, 0 };
        
        for (int j = 0; j < N; j++) {
            if (i == j) continue;

            Vec3 r = cur_state[j].position - cur_state[i].position;
            double r2 = length2(r);
            delta[i].velocity += (r / sqrt(r2)) * (GRAVITY_CONSTANT * data[j].mass / r2);  
        }
    }
}

void S_CalculateJacobianBlock(BodyState *state, int i, int j, double *out, int stride)
{
    if (i == j) {
        out[stride*0 + 3] = 1.0;
        out[stride*1 + 4] = 1.0;
        out[stride*2 + 5] = 1.0;
        
        double res_xx = 0.0, res_xy = 0.0, res_xz = 0.0, res_yy = 0.0, res_yz = 0.0, res_zz = 0.0;
        for (int k = 0; k < N; k++) {
            if (k == i) continue;
            Vec3 R = state[k].position - state[i].position;
            double r2 = length2(R);
            double r3 = r2*sqrt(r2);
            double r5 = r3*r2;
            double m = data[k].mass;

            res_xx += m * (1.0 / r3 - 3*R.x*R.x / r5);
            res_xy += m * R.x * R.y / r5;
            res_xz += m * R.x * R.z / r5;
            res_yy += m * (1.0 / r3 - 3*R.y*R.y / r5);
            res_yz += m * R.y * R.z / r5;
            res_zz += m * (1.0 / r3 - 3*R.z*R.z / r5);
        }

        out[stride*3 + 0] = -GRAVITY_CONSTANT * res_xx;
        out[stride*3 + 1] = 3*GRAVITY_CONSTANT * res_xy;
        out[stride*3 + 2] = 3*GRAVITY_CONSTANT * res_xz;
        out[stride*4 + 0] = 3*GRAVITY_CONSTANT * res_xy;
        out[stride*4 + 1] = -GRAVITY_CONSTANT * res_yy;
        out[stride*4 + 2] = 3*GRAVITY_CONSTANT * res_yz;
        out[stride*5 + 0] = 3*GRAVITY_CONSTANT * res_xz;
        out[stride*5 + 1] = 3*GRAVITY_CONSTANT * res_yz;
        out[stride*5 + 2] = -GRAVITY_CONSTANT * res_zz;
    }
    else {
        Vec3 R = state[j].position - state[i].position;
        double r2 = length2(R);
        double r3 = r2*sqrt(r2);
        double r5 = r3*r2;
        double m = data[j].mass;

        out[stride*3 + 0] = GRAVITY_CONSTANT * m * (1.0 / r3 - 3.0*R.x*R.x / r5);
        out[stride*3 + 1] = 3 * GRAVITY_CONSTANT * m * R.x * R.y / r5;
        out[stride*3 + 2] = 3 * GRAVITY_CONSTANT * m * R.x * R.z / r5;

        out[stride*4 + 0] = out[stride*3 + 1];
        out[stride*4 + 1] = GRAVITY_CONSTANT * m * (1.0 / r3 - 3.0*R.y*R.y / r5);
        out[stride*4 + 2] = 3 * GRAVITY_CONSTANT * m * R.y * R.z / r5;

        out[stride*5 + 0] = out[stride*3 + 2];
        out[stride*5 + 1] = out[stride*4 + 2];
        out[stride*5 + 2] = GRAVITY_CONSTANT * m * (1.0 / r3 - 3.0*R.z*R.z / r5);
    }
}

void S_BuildJacobian(BodyState *state, double *jacobian)
{
    int row = N * 6;
    for (int i = 0; i < row*row; i++) {
        jacobian[i] = 0.0;
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double *loc = &jacobian[i*6*row + j*6];
            S_CalculateJacobianBlock(state, i, j, loc, row);
        }
    }
}

void S_AddState(BodyState *dst, BodyState *src, double k)
{
    for (int i = 0; i < N; i++) {
        dst[i].position += src[i].position * k;
        dst[i].velocity += src[i].velocity * k;
    }
}

double S_CalculateStateLength(BodyState *state)
{
    double *arr = (double*)state;
    double res = 0.0;

    for (int i = 0; i < 6*N; i++) {
        res += arr[i]*arr[i];
    }

    return sqrt(res);
}

void S_MulMatVec(double *mat, double *vec, int n, double *res)
{
    for (int i = 0; i < n; i++) {
        res[i] = 0.0;
        for (int j = 0; j < n; j++) {
            res[i] += mat[i*n + j] * vec[j];
        }
    }
}

void S_Step_Implicit(double dt) {
    const double tolerance = 1e-10;
    const int max_iter = 10;
    
    double *J = (double*)malloc(N * N * 36 * sizeof(double));
    double *invJ = (double*)malloc(N * N * 36 * sizeof(double));
    BodyState *F = (BodyState*)malloc(N * sizeof(BodyState));
    BodyState *tmp = (BodyState*)malloc(N * sizeof(BodyState));
    BodyState *delta = (BodyState*)malloc(N * sizeof(BodyState));
    BodyState *guess = (BodyState*)malloc(N * sizeof(BodyState));

    // euler
    // guess = state + delta(state) * dt
    memcpy(guess, state, sizeof(BodyState) * N);
    S_CalculateDelta(state, delta);
    S_AddState(guess, delta, dt);

    int iter = 0;
    double length;
    for (; iter < max_iter; iter++)
    {
        // F = guess - state - delta(guess) * dt;
        memcpy(F, guess, sizeof(BodyState) * N);
        S_CalculateDelta(guess, delta);
        S_AddState(F, state, -1);
        S_AddState(F, delta, -dt);
        
        // if (F.length() < tolerance) break;
        length = S_CalculateStateLength(F);
        if (length < tolerance) {
            break;
        }

        S_BuildJacobian(state, J);
        M_InverseMatrix(J, N * 6, invJ);

        // guess = guess - invJ * F;
        S_MulMatVec(invJ, (double*)F, N * 6, (double*)tmp);
        S_AddState(guess, tmp, -1);
    }

    printf("iter = %d\n, length=%lf\n", iter, length);
    memcpy(state, guess, N * sizeof(BodyState));

    free(J);
    free(invJ);
    free(delta);
    free(guess);
    free(tmp);
    free(F);
}

void S_Step_AdamsPredictorCorrector(double dt) {
    static BodyState *prev_k[4];
    static bool initialized = false;

    if (!initialized) {
        for (int i = 0; i < 4; i++) {
            prev_k[i] = (BodyState*)malloc(N * sizeof(BodyState));
            S_CalculateDelta(state, prev_k[i]);
        }
        initialized = true;
    }

    // Предиктор (Адамс-Бэшфорт 4-го порядка)
    BodyState *predictor = (BodyState*)malloc(N * sizeof(BodyState));
    for (int i = 0; i < N; i++) {
        predictor[i].position = state[i].position + (
            prev_k[0][i].position * 55.0 - prev_k[1][i].position * 59.0 +
            prev_k[2][i].position * 37.0 - prev_k[3][i].position * 9.0
        ) * dt / 24.0;

        predictor[i].velocity = state[i].velocity + (
            prev_k[0][i].velocity * 55.0 - prev_k[1][i].velocity * 59.0 +
            prev_k[2][i].velocity * 37.0 - prev_k[3][i].velocity * 9.0
        ) * dt / 24.0;
    }

    // Корректор (Адамс-Моултон 4-го порядка)
    BodyState *k_new = (BodyState*)malloc(N * sizeof(BodyState));
    S_CalculateDelta(predictor, k_new);

    for (int i = 0; i < N; i++) {
        state[i].position += (
            k_new[i].position * 9.0 + prev_k[0][i].position * 19.0 -
            prev_k[1][i].position * 5.0 + prev_k[2][i].position
        ) * dt / 24.0;

        state[i].velocity += (
            k_new[i].velocity * 9.0 + prev_k[0][i].velocity * 19.0 -
            prev_k[1][i].velocity * 5.0 + prev_k[2][i].velocity
        ) * dt / 24.0;
    }
    
    for (int i = 3; i > 0; i--) {
        memcpy(prev_k[i], prev_k[i-1], N * sizeof(BodyState));
    }
    memcpy(prev_k[0], k_new, N * sizeof(BodyState));

    free(predictor);
    free(k_new);
}

void S_Step_Explicit(BodyState* current_state, BodyState* next_state, double dt)
{
    BodyState *tmp_state = (BodyState*)malloc(sizeof(BodyState) * N);
    BodyState **deltas = (BodyState**)malloc(sizeof(BodyState*) * method.n);
    
    for (int i = 0; i < method.n; i++) {
        deltas[i] = (BodyState*)malloc(sizeof(BodyState) * N);
    }

    if (next_state != current_state) {
        memcpy(next_state, current_state, sizeof(BodyState) * N);
    }

    for (int i = 0; i < method.n; i++) {
        memcpy(tmp_state, current_state, sizeof(BodyState) * N);
        for (int j = 0; j < i; j++) {
            for (int k = 0; k < N; k++) {
                tmp_state[k].position += deltas[j][k].position * method.a[i][j] * dt;
                tmp_state[k].velocity += deltas[j][k].velocity * method.a[i][j] * dt;
            }
        }
        
        S_CalculateDelta(tmp_state, deltas[i]);
    }

    for (int i = 0; i < N; i++) {
        BodyState sum = {};
        for (int j = 0; j < method.n; j++) {
            sum.position += deltas[j][i].position * method.k[j] * dt;
            sum.velocity += deltas[j][i].velocity * method.k[j] * dt;
        }

        next_state[i].position += sum.position;
        next_state[i].velocity += sum.velocity;
    }

    for (int i = 0; i < method.n; i++) {
        free(deltas[i]);
    }
    free(deltas);
    free(tmp_state);
}

void S_Step_Explicit(double dt) {
    S_Step_Explicit(state, state, dt);
}

int S_Step_Adaptive(AdaptiveParams* params) {
    static BodyState *state1 = NULL;
    static BodyState *state2 = NULL;
    static BodyState *temp_state = NULL;

    if (state1 == NULL) {
        state1 = (BodyState*)malloc(sizeof(BodyState) * N);
        state2 = (BodyState*)malloc(sizeof(BodyState) * N);
        temp_state = (BodyState*)malloc(sizeof(BodyState) * N);
    }
    
    // Первый шаг с current_dt
    S_Step_Explicit(state, state1, params->current_dt);
    
    // Два шага с current_dt/2
    memcpy(temp_state, state, sizeof(BodyState) * N);
    S_Step_Explicit(temp_state, temp_state, params->current_dt/2);
    S_Step_Explicit(temp_state, state2, params->current_dt/2);
    
    // Оценка ошибки
    double error = 0.0;
    for (int i = 0; i < N; i++) {
        Vec3 pos_diff = state1[i].position - state2[i].position;
        Vec3 vel_diff = state1[i].velocity - state2[i].velocity;
        error += length2(pos_diff) + length2(vel_diff);
    }
    error = sqrt(error / (6*N));
    
    // Адаптация шага
    if (error < params->tolerance) {
        memcpy(state, state2, sizeof(BodyState) * N);
        params->sim_time += params->current_dt;

        if (error < params->tolerance / 10) {
            params->current_dt = fmin(params->current_dt * params->dt_scale, params->max_dt);
        }
        return 1;
    } else {
        params->current_dt = fmax(params->current_dt * params->safety * pow(params->tolerance/error, 0.2), 
                  params->min_dt);
        return 0; 
    }
}

void S_ReadState()
{
    scanf("%d", &N);

    state = (BodyState*)malloc(sizeof(BodyState) * N);
    data = (BodyData*)malloc(sizeof(BodyData) * N);

    for (int i = 0; i < N; i++) {
        scanf("%lf %lf %lf %lf %lf %lf %lf",
            &data[i].mass,
            &state[i].position.x, &state[i].position.y, &state[i].position.z,
            &state[i].velocity.x, &state[i].velocity.y, &state[i].velocity.z
        );

        for (int j = 0; j < TRAIL_SIZE; j++) {
            data[i].trail[j] = state[i].position;
        }
    }

    printf("%d\n", N);
}

void S_ReadMethod(const char *filename)
{
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        printf("Error: Could not open method file %s\n", filename);
        exit(1);
    }

    fscanf(fp, "%d", &method.n);

    method.k = (double*)malloc(sizeof(double) * method.n);
    method.a = (double**)malloc(sizeof(double*) * method.n);
    
    for (int i = 1; i < method.n; i++) {
        method.a[i] = (double*)malloc(sizeof(double) * i);
        
        for (int j = 0; j < i; j++) {
            char buffer[32];
            fscanf(fp, "%s", buffer);

            if (strchr(buffer, '/')) {
                double a, b;
                sscanf(buffer, "%lf/%lf", &a, &b);
                method.a[i][j] = a / b;
            }
            else {
                double a;
                sscanf(buffer, "%lf", &a);
                method.a[i][j] = a;
            }
        }
    }

    for (int i = 0; i < method.n; i++) {
        char buffer[32];
        fscanf(fp, "%s", buffer);
        
        if (strchr(buffer, '/')) {
            double a, b;
            sscanf(buffer, "%lf/%lf", &a, &b);
            method.k[i] = a / b;
        } else {
            double a;
            sscanf(buffer, "%lf", &a);
            method.k[i] = a;
        }
    }

    fclose(fp);
}

void S_Init(const char *method_filename)
{
    S_ReadMethod(method_filename);
    S_ReadState();
}
