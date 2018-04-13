#pragma once
#include <iostream>
#include <string>

/* Scene Constants */
static const char *TITLE = "Stable Fluids";
constexpr int WIDTH = 720;
constexpr int HEIGHT = 720;

/* Simulator Constants */
constexpr int DIM = 2;
constexpr int LENGTH = 1.0;
constexpr int N = 64;
constexpr int SOURCE_SIZE = N / 3;

constexpr float INTERACTION = 1000000.0f;

constexpr float VISCOSITY = 0.001f;
constexpr float GRAVITY_Y = 9.8f;
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

enum EMode
{
    E_Continuous = 0,
    E_Once = 1
};