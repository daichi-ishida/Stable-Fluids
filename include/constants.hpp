#pragma once
#include <iostream>

constexpr int DIM = 2;
constexpr int LENGTH = 1.0f;
constexpr int N = 64;
constexpr bool USE_VORTEX_PARTICLES = true;

constexpr float VISCOSITY = 0.001f;
constexpr float DIFFUSION = 0.001f;
constexpr float VORT_EPS = 1000.0;
constexpr float GRAVITY_Y = -9.8f;
constexpr float AMBIENT_TEMP = 30.0f;
constexpr float DT = 0.01f;

enum
{
    SHOW_FORCE,
    SHOW_SPEED,
    SHOW_DENSITY,
    SHOW_TEXTURE
};

int showItem = SHOW_DENSITY;

constexpr int MATSIZE()
{
    switch (DIM)
    {
    case 2:
        return (N + 2) * (N + 2);
        break;
    case 3:
        return (N + 2) * (N + 2) * (N + 2);
        break;
    default:
        static_assert(true, "Invalid dimension assignment.");
        exit(EXIT_FAILURE);
    }
}

constexpr int SIZE = MATSIZE();

constexpr int POS(int i, int j) { return i + (N + 2) * j; };
constexpr int POS(int i, int j, int k) { return i + (N + 2) * j + (N + 2) * (N + 2) * k; };