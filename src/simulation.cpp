#include "simulation.hpp"
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

void S_Step(double dt)
{
    BodyState* tmp_state = (BodyState*)malloc(sizeof(BodyState) * N);
    BodyState* deltas = (BodyState*)malloc(sizeof(BodyState) * N * method.n);

    for (int i = 0; i < method.n; i++) {
        memcpy(tmp_state, state, sizeof(BodyState) * N);

        for (int j = 0; j < i; j++) {
            for (int k = 0; k < N; k++) {
                tmp_state[k].position += deltas[j * N + k].position * method.a[i][j] * dt;
                tmp_state[k].velocity += deltas[j * N + k].velocity * method.a[i][j] * dt;
            }
        }

        S_CalculateDelta(tmp_state, &deltas[i * N]);
    }

    for (int i = 0; i < N; i++) {
        BodyState sum = {};
        for (int j = 0; j < method.n; j++) {
            sum.position += deltas[j * N + i].position * method.k[j] * dt;
            sum.velocity += deltas[j * N + i].velocity * method.k[j] * dt;
        }
        state[i].position += sum.position;
        state[i].velocity += sum.velocity;
    }

    free(deltas);
    free(tmp_state);
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
