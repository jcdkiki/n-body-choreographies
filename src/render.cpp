#include "render.hpp"
#include "SDL_timer.h"
#include "SDL_video.h"
#include "simulation.hpp"
#include "common.hpp"

#include <GL/glut.h>
#include <SDL2/SDL.h>

constexpr double ROLL_FACTOR = 0.05;
constexpr double TRAIL_UPDATE_TIME = 1.0 / 30.0;
constexpr int WINDOW_WIDTH = 800;
constexpr int WINDOW_HEIGHT = 600;
constexpr double CAMERA_SPEED = 10.0;

static SDL_Window *window;
static SDL_Renderer *renderer;
static SDL_TimerID timer_id;
static SDL_GLContext context;

static int cur_trail = 0;
static Vec3 camera_pos = { 0, 5, 0 };
static Vec3 camera_rotation = { 0, 0, -M_PI/2 };
static Vec3 camera_direction;
static Vec3 camera_direction_normal;
static Vec3 camera_velocity;
static double trail_time = 0.0;

static const Uint8 *keys;

void R_Draw()
{
    glClearColor(0.05f, 0.05f, 0.1f, 1.f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glViewport(0, 0, WINDOW_WIDTH, WINDOW_HEIGHT);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    double aspect = (double)WINDOW_WIDTH / WINDOW_HEIGHT;
    gluPerspective(60.0, aspect, 0.01, 100.0);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    double sideways = (camera_direction_normal.x * camera_velocity.x + camera_direction_normal.y * camera_velocity.y) * ROLL_FACTOR;

    Vec3 center_pos = camera_pos + camera_direction;
    gluLookAt(
        camera_pos.x, camera_pos.y, camera_pos.z, 
        center_pos.x, center_pos.y, center_pos.z,
        -sin(sideways)*camera_direction.y, sin(sideways)*camera_direction.x, cos(sideways)
    );

    glDisable(GL_LIGHTING);
    glBegin(GL_LINES);
    glColor3f(1.f, 0.f, 0.f);
    glVertex3f(0, 0, 0);
    glVertex3f(10, 0, 0);
    
    glColor3f(0.f, 1.f, 0.f);
    glVertex3f(0, 0, 0);
    glVertex3f(0, 10, 0);
    
    glColor3f(0.f, 0.f, 1.f);
    glVertex3f(0, 0, 0);
    glVertex3f(0, 0, 10);
    glEnd();
    
    float ambient[3] = { 0.05, 0.05, 0.1 };
    float light0_pos[3] = { 0.f, 0.f, 0.f };
    float light1_pos[3] = { (float)camera_pos.x, (float)camera_pos.y, (float)camera_pos.z };
    float specular[4] = { 0.9, 0.5, 1, 1 };

    glLightModelfv(GL_AMBIENT, ambient);
    glShadeModel(GL_SMOOTH);

    glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_AMBIENT_AND_DIFFUSE);


	glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
	glLightfv(GL_LIGHT0, GL_POSITION, light0_pos);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
	glMaterialf (GL_FRONT_AND_BACK, GL_SHININESS, 30);

    for (int i = 0; i < N; i++) {
        glEnable(GL_LIGHTING);
        float color[3] = { 1.f * (i + 1) / (float)N, 1.f * (N - i) / N, 0.7f }; 
        glColor3fv(color);
        glPushMatrix();
        glTranslatef(state[i].position.x, state[i].position.y, state[i].position.z);
        
        GLUquadric *quadric = gluNewQuadric();
        gluSphere(quadric, 0.1, 32, 32);
        gluDeleteQuadric(quadric);

        glPopMatrix();
        
        glDisable(GL_LIGHTING);
        glBegin(GL_LINE_STRIP);
        for (int j = 0; j < TRAIL_SIZE; j++) {
            int idx = (cur_trail - j + TRAIL_SIZE) % TRAIL_SIZE;
            float alpha = 1 - j / (float)TRAIL_SIZE;
            glColor4f(color[0], color[1], color[2], alpha*alpha);
            glVertex3f(data[i].trail[idx].x, data[i].trail[idx].y, data[i].trail[idx].z);
        }
        glEnd();
    }

    glColor4f(1, 1, 1, 0.25);
    glBegin(GL_LINES);
    for (int i = -10; i <= 10; i++) {
        glVertex3f(i, 0, -10);
        glVertex3f(i, 0, 10);
        glVertex3f(-10, 0, i);
        glVertex3f(10, 0, i);
    }
    glEnd();
    
    SDL_GL_SwapWindow(window);
}

void R_Step(double dt)
{
    camera_direction = {
        cos(camera_rotation.z)*cos(camera_rotation.x),
        sin(camera_rotation.z)*cos(camera_rotation.x),
        sin(camera_rotation.x)
    };

    camera_direction_normal = {
        -sin(camera_rotation.z),
        cos(camera_rotation.z),
        0
    };

    SDL_PumpEvents();
    if (keys[SDL_SCANCODE_A]) camera_velocity += camera_direction_normal * CAMERA_SPEED * dt;
    if (keys[SDL_SCANCODE_D]) camera_velocity -= camera_direction_normal * CAMERA_SPEED * dt;
    if (keys[SDL_SCANCODE_W]) camera_velocity += camera_direction * CAMERA_SPEED * dt;
    if (keys[SDL_SCANCODE_S]) camera_velocity -= camera_direction * CAMERA_SPEED * dt;
    camera_velocity *= pow(0.1, dt);
    camera_pos += camera_velocity * dt;

    trail_time += dt;
    while (trail_time >= TRAIL_UPDATE_TIME) {
        cur_trail = (cur_trail + 1) % TRAIL_SIZE;
        for (int i = 0; i < N; i++) {
            data[i].trail[cur_trail] = state[i].position;
        }

        trail_time -= TRAIL_UPDATE_TIME;
    }
}

void R_Init(unsigned int (*callback)(unsigned int interval, void *name))
{
    SDL_Init(SDL_INIT_EVERYTHING | SDL_INIT_NOPARACHUTE);

    window = SDL_CreateWindow("main",
        SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
        WINDOW_WIDTH, WINDOW_HEIGHT,
        SDL_WINDOW_OPENGL
    );

    context = SDL_GL_CreateContext(window);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    SDL_CaptureMouse(SDL_TRUE);
    SDL_SetRelativeMouseMode(SDL_TRUE);
    SDL_TimerID timer_id = SDL_AddTimer(STEP_MS, callback, NULL);
    keys = SDL_GetKeyboardState(NULL);
}

void R_Deinit()
{
    SDL_RemoveTimer(timer_id);
    SDL_GL_DeleteContext(context);
    SDL_DestroyWindow(window);
    SDL_Quit();
}

int R_GetInput(void) {
    SDL_Event event;
    while (SDL_PollEvent(&event)) {
        switch (event.type) {
            case SDL_QUIT: return 0;
            case SDL_KEYDOWN:
                switch (event.key.keysym.sym) {
                    case SDLK_ESCAPE:
                        return 0;
                }
                break;
        
            case SDL_MOUSEMOTION:
                camera_rotation.z -= event.motion.xrel * 0.01;
                camera_rotation.x -= event.motion.yrel * 0.01;
                if (camera_rotation.x < -19*M_PI/40) camera_rotation.x = -19*M_PI/40;
                if (camera_rotation.x > 19*M_PI/40) camera_rotation.x = 19*M_PI/40;
                break;
        }
    }

    return 1;
}
