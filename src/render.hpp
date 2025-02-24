#ifndef RENDER_HPP
#define RENDER_HPP

void R_Init(unsigned int (*callback)(unsigned int interval, void *name));
void R_Deinit();
void R_Draw();
int  R_GetInput();
void R_Step(double dt);

#endif
