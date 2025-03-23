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

/*
void S_CalculateJacobianCij(double mj, Vec3 *ri, Vec3 *rj, double *out, int stride)
{
    double r = length(*rj - *ri);
    double r3 = r*r*r;
    double r5 = r3*r*r;
    double gmj = GRAVITY_CONSTANT * mj;

    double Rx = (rj->x - ri->x);
    double Ry = (rj->y - ri->y);
    double Rz = (rj->z - ri->z);

    out[stride*0 + 0] = gmj * (1.0 / r3 - 3.0*Rx*Rx / r5);  // C_ij_xx
    out[stride*0 + 1] = -3.0 * gmj * Rx * Ry / r5;          // C_ij_xy   
    out[stride*0 + 2] = -3.0 * gmj * Rx * Rz / r5;          // C_ij_xz 

    out[stride*1 + 0] = out[stride*0 + 1];                  // C_ij_yx
    out[stride*1 + 1] = gmj * (1.0 / r3 - 3.0*Ry*Ry / r5);  // C_ij_yy
    out[stride*1 + 2] = -3.0 * gmj * Ry * Rz / r5;          // C_ij_yz
    
    out[stride*2 + 0] = out[stride*0 + 2];                  // C_ij_zx
    out[stride*2 + 1] = out[stride*1 + 2];                  // C_ij_zy
    out[stride*2 + 2] = gmj * (1.0 / r3 - 3.0*Rz*Rz / r5);  // C_ij_zz
}

void S_CalculateJacobianABi(BodyState *state, int i, int n_bodies, double *out, int stride)
{
    double res_xx = 0.0, res_xy = 0.0, res_xz = 0.0,
                         res_yy = 0.0, res_yz = 0.0,
                                       res_zz = 0.0;

    for (int j = 0; j < n_bodies; j++) {
        if (i == j) continue;

        Vec3 R = state[j].position - state[i].position;
        double r = length(R);
        double r3 = r*r*r;
        double r5 = r3*r*r;
        double mj = data[j].mass;

        res_xx += mj * (1.0 / r3 - 3.0*R.x*R.x / r5);
        res_yy += mj * (1.0 / r3 - 3.0*R.y*R.y / r5);
        res_zz += mj * (1.0 / r3 - 3.0*R.z*R.z / r5);

        res_xy += mj * R.x * R.y / r5;
        res_xz += mj * R.x * R.z / r5;
        res_yz += mj * R.y * R.z / r5;
    }

    out[stride*0 + 0] = -GRAVITY_CONSTANT * res_xx;
    out[stride*1 + 1] = -GRAVITY_CONSTANT * res_yy;
    out[stride*2 + 2] = -GRAVITY_CONSTANT * res_zz;

    out[stride*0 + 1] = 3*GRAVITY_CONSTANT * res_xy;
    out[stride*0 + 2] = 3*GRAVITY_CONSTANT * res_xz;
    out[stride*1 + 2] = 3*GRAVITY_CONSTANT * res_yz;

    out[stride*1 + 0] = out[stride*0 + 1];
    out[stride*2 + 0] = out[stride*0 + 2];
    out[stride*2 + 1] = out[stride*1 + 2];
}

void S_BuildJacobian(BodyState *state, int n_bodies, double *jacobian)
{
    int row = n_bodies * 6;
    for (int i = 0; i < row * row; i++) {
        jacobian[i] = 0.0;
    }
    
    for (int i = 0; i < n_bodies * 3; i++) {
        jacobian[i*row + i + n_bodies*3] = 1.0;
    }

    for (int i = 0; i < n_bodies; i++) {
        double *loc = &jacobian[(i*3 + n_bodies*3) * row + i*3];
        S_CalculateJacobianABi(state, i, n_bodies, loc, row);

        for (int j = 0; j < n_bodies; j++) {
            if (i == j) continue;

            loc = &jacobian[(i*3 + n_bodies*3) * row + j*3];
            S_CalculateJacobianCij(data[j].mass, &state[i].position, &state[j].position, loc, row);
        }
    }
}*/

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

/*
void S_Step(double dt) {
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
*/


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
