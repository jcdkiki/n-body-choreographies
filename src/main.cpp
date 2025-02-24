#include "simulation.hpp"
#include "render.hpp"
#include <cstdio>

unsigned int callback(unsigned int interval, void *name)
{
    double dt = (double)interval / 1000.f;

    S_Step(dt);
    R_Step(dt);

    return interval;
}

int main(int argc, char **argv)
{
    S_Init(argv[1]);
    R_Init(callback);

    while (1) {
        if (!R_GetInput()) break;
        R_Draw();
    }

    R_Deinit();
    //S_Deinit();
}
