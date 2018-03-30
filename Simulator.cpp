#include "Simulator.hpp"

Simulator::Simulator(GridCells grid_cells) : m_grid_cells(grid_cells)
{
}

Simulator::~Simulator()
{
}

void Simulator::update()
{
}

/* private */
void Simulator::vStep(float u[], float u_prev[], float visc, float F[], float dt)
{
    for (int i = 0; i < LENGTH; ++i)
    {
        addForce(u_prev[i], F[i], dt);
    }
    for (int i = 0; i < LENGTH; ++i)
    {
        transport(u[i], u_prev[i], u_prev, dt);
    }
    for (int i = 0; i < LENGTH; ++i)
    {
        diffuse(u_prev[i], u[i], visc, dt);
    }
}
void Simulator::sStep(float dens[], float dens_prev[], float kS, float aS, float u[], float source, float dt)
{
    addForce(dens_prev, source, dt);
    transport(s, s_prev, u, dt);
}

void Simulator::swap(float arr[], float arr_prev[], int size)
{
    for (int i = 0; i < LENGTH; ++i)
    {
        float tmp = arr_prev[i];
        arr_prev[i] = arr[i];
        arr[i] = arr_prev[i];
    }
}

void Simulator::addSource(int LENGTH, float x[], float s[], float dt)
{
    for (int i = 0; i < LENGTH; ++i)
    {
        x[i] += dt * s[i];
    }
}