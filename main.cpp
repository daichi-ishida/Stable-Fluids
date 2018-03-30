#include <iostream>
#include <cmath>
#include <Eigen/Core>

constexpr int DIM = 2;
constexpr int LENGTH = 300;
constexpr int SIZE = (LENGTH + 2) * (LENGTH + 2);
constexpr int POS(int i, int j) { return i + (LENGTH + 2) * j; };

static double u[SIZE], v[SIZE], u_prev[SIZE], v_prev[SIZE];
