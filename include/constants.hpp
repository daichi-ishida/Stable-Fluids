#pragma once
#include <iostream>

constexpr int DIM = 2;
constexpr int LENGTH = 300;
constexpr int N = 100;

constexpr float VISCOSITY = 1.0E-6f;
constexpr float DIFFUSION = 1.0E-6f;
constexpr float GRAVITY_Y = -9.8f;
constexpr float DT = 0.1f;

constexpr int MATSIZE()
{
    switch (DIM)
    {
    case 2:
        return (LENGTH + 2) * (LENGTH + 2);
        break;
    case 3:
        return (LENGTH + 2) * (LENGTH + 2) * (LENGTH + 2);
        break;
    default:
        static_assert(true, "Invalid dimension assignment.");
        exit(EXIT_FAILURE);
    }
}

constexpr int SIZE = MATSIZE();

constexpr int POS(int i, int j) { return i + (LENGTH + 2) * j; };
constexpr int POS(int i, int j, int k) { return i + (LENGTH + 2) * j + (LENGTH + 2) * (LENGTH + 2) * k; };