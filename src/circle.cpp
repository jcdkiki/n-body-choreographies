#include <stdio.h>
#include <math.h>
#include "common.hpp"

int main(int argc, char **argv)
{
    int N;
    sscanf(argv[1], "%d", &N);

    printf("%d\n", N);
    for (int i = 0; i < N; i++) {
        double mass = 1.0;
        double required_velocity = sqrt(GRAVITY_CONSTANT * mass / sqrt(N));
        Vec3 position = { cos(i*2*M_PI/N), 0, sin(i*2*M_PI/N) };
        Vec3 velocity = { -sin(i*2*M_PI/N) * required_velocity, 0, cos(i*2*M_PI/N) * required_velocity };
        
        printf("%lf ", mass);
        printf("%lf %lf %lf ", position.x, position.y, position.z);
        printf("%lf %lf %lf\n", velocity.x, velocity.y, velocity.z);
    }
}
