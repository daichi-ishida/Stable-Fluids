#pragma once
#include "GridCells.hpp"

class Simulator
{
public:
    Simulator(GridCells *grid_cells);
    ~Simulator();

    void initSet();
    void update();

private:
    void vStep(float u[], float u_prev[], float visc, float F[], float dt);
    void sStep(float dens[], float dens_prev[], float kS, float aS, float u[], float source, float dt);

    void swap(float u[], float u_prev[]);
    void addForce(float u_prev_i, float &F_i, float dt);
    void transport(float u_i, float u_prev_i, float u_prev[], float dt);
    void diffuse(float u_prev_i, float u_i, float visc, float dt);
    void project(float u[], float u_prev[], float dt);
    void addSource(int LENGTH, float &x, float &s, float dt);

    GridCells* m_grid_cells;
}