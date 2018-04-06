#pragma once
#include "constants.hpp"

class GridCells2D
{
public:
  GridCells2D();
  ~GridCells2D();

  const float getDensColor(int i, int j) { return dens[POS(i, j)]; };
  const float getTempColor(int i, int j) { return temp[POS(i, j)] / TEMPERATURE_MAX; };

  float u[SIZE], v[SIZE];
  float u0[SIZE], v0[SIZE];
  float dens[SIZE];
  float temp[SIZE];
};