#include <cmath>
#include "Simulator2D.hpp"
#include "interpolation.hpp"

Simulator2D::Simulator2D(GridCells2D &grid_cells) : m_grid_cells(grid_cells)
{
}

Simulator2D::~Simulator2D()
{
}

void Simulator2D::update()
{
}

/* private */

void Simulator2D::velocityStep()
{
    addForce();
    advect();
    FFT();
    diffuse();
    project();
    IFFT();
}

// void Simulator2D::densityStep()
// {
//     addSource(m_grid_cells.dens, m_grid_cells.dens_prev);
//     diffuse(m_grid_cells.dens, m_grid_cells.dens_prev, 0);
//     advect(m_grid_cells.dens, m_grid_cells.dens_prev, 0);
// }

void Simulator2D::addForce()
{
    for (int i = 0; i < SIZE; ++i)
    {
        m_grid_cells.v_x[i] += DT * m_grid_cells.v_x_prev[i];
        m_grid_cells.v_y[i] += DT * m_grid_cells.v_y_prev[i];

        m_grid_cells.v_x_prev[i] = m_grid_cells.v_x[i];
        m_grid_cells.v_y_prev[i] = m_grid_cells.v_y[i];
    }
}

void Simulator2D::advect()
{
    for (int i = 1; i < N + 1; ++i)
    {
        for (int j = 1; j < N + 1; ++j)
        {
            float x = constrainValue(i - N * DT * m_grid_cells.v_x_prev[POS(i, j)]);
            float y = constrainValue(j - N * DT * m_grid_cells.v_y_prev[POS(i, j)]);

            int i0 = (int)std::floor(x);
            int j0 = (int)std::floor(y);
            float s = x - (float)i0;
            float t = y - (float)j0;

            int i1 = i0 + 1;
            int j1 = j0 + 1;

            // grid interpolation - linear
            float intrpl_vx_vert = lerp::linear(m_grid_cells.v_x_prev[POS(i0, j0)], m_grid_cells.v_x_prev[POS(i0, j1)], t);
            float intrpl_vx_horiz = lerp::linear(m_grid_cells.v_x_prev[POS(i1, j0)], m_grid_cells.v_x_prev[POS(i1, j1)], t);

            float intrpl_vy_vert = lerp::linear(m_grid_cells.v_y_prev[POS(i0, j0)], m_grid_cells.v_y_prev[POS(i0, j1)], t);
            float intrpl_vy_horiz = lerp::linear(m_grid_cells.v_y_prev[POS(i1, j0)], m_grid_cells.v_y_prev[POS(i1, j1)], t);

            m_grid_cells.v_x[POS(i, j)] = lerp::linear(intrpl_vx_vert, intrpl_vx_horiz, s);
            m_grid_cells.v_y[POS(i, j)] = lerp::linear(intrpl_vy_vert, intrpl_vy_horiz, s);
        }
    }
}

void Simulator2D::FFT()
{

}

// void Simulator2D::diffuse(float current[], float previous[], int b)
// {
//     float diff_rate = DT * DIFFUSION * LENGTH * LENGTH;
//     for (int k = 0; k < ITERATIONS; ++k)
//     {
//         for (int i = 1; i < LENGTH; ++i)
//         {
//             for (int j = i; j < LENGTH; ++j)
//             {
//                 current[POS(i, j)] = (previous[POS(i, j)] + diff_rate * calNeighborSum(current, i, j)) / (1 + 4 * diff_rate);
//             }
//         }
//         setBounds(current, b);
//     }
// }

float Simulator2D::constrainValue(float value)
{
    if (value < 0.5)
    {
        return 0.5;
    }
    else if (value > N + 0.5)
    {
        return N + 0.5;
    }
    else
    {
        return value;
    }
}