#pragma once
#include "constants.hpp"

class GridCells2D
{
public:
  GridCells2D();
  ~GridCells2D();

  float v_x[SIZE], v_y[SIZE];
  float v_x_prev[SIZE], v_y_prev[SIZE];
  float dens[SIZE], dens_prev[SIZE];
};