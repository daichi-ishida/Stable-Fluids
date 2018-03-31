#pragma once

constexpr int DIM = 2;
constexpr int LENGTH = 300;
constexpr int SIZE = (LENGTH + 2) * (LENGTH + 2);
constexpr int POS(int i, int j) { return i + (LENGTH + 2) * j; };