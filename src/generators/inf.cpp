#include <stdio.h>
#include <math.h>
#include "../simulation.hpp"

int main()
{
    printf("%d\n", 3);
    
    BodyState bodies[3] = {
        BodyState { {-1, 0, 0}, {0, 0, -1} },
        BodyState { {0, 0, 0}, {sqrt(2)/2, 0, sqrt(2)/2} },
        BodyState { {1, 0, 0}, {0, 0, 1} }
    };

    for (int i = 0; i < 3; i++) {
        //bodies[i].velocity *= 0.5;

        printf("1 %lf %lf %lf %lf %lf %lf\n",
            bodies[i].position.x, bodies[i].position.y, bodies[i].position.z,
            bodies[i].velocity.x, bodies[i].velocity.y, bodies[i].velocity.z
        );
    }

    return 0;
}
