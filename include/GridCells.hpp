#pragma once
#include "constants.hpp"

class GridCells
{
public:
    GridCells();
    ~GridCells();

    float u[SIZE], v[SIZE];
    float u_prev[SIZE], v_prev[SIZE];
    float dens[LENGTH], dens_prev[LENGTH];
}