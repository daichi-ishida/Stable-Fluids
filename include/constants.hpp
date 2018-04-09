#pragma once
#include <iostream>
#include <string>

constexpr float INTERACTION = 80.0f;

/* Scene Constants */
static const char *TITLE = "Stable Fluids";
constexpr int WIDTH = 720;
constexpr int HEIGHT = 720;
constexpr int FPS_CAP = 120;
constexpr float FOV = 70;
constexpr float NEAR_PLANE = 0.1f;
constexpr float FAR_PLANE = 1000;

/* Simulator Constants */
constexpr int DIM = 2;
constexpr int LENGTH = 1.0;
constexpr int N = 32;
constexpr int PARTICLE_NUM = 100;
constexpr bool USE_VORTEX_PARTICLES = false;
constexpr float TEMPERATURE_MAX = 200.0f;

constexpr float VISCOSITY = 0.001f;
constexpr float VORT_EPS = 0.01f;
constexpr float GRAVITY_Y = -0.98f;
constexpr float AMBIENT_TEMP = 30.0f;
constexpr float DT = 0.01f;
constexpr float FINISH_TIME = 3.0f;

constexpr int MATSIZE()
{
    switch (DIM)
    {
    case 2:
        return N * N;
        break;
    case 3:
        return N * N * N;
        break;
    default:
        static_assert(true, "Invalid dimension assignment.");
        exit(EXIT_FAILURE);
    }
}

constexpr int SIZE = MATSIZE();

constexpr int POS(int i, int j) { return i + N * j; };
constexpr int POS(int i, int j, int k) { return i + N * j + N * N * k; };