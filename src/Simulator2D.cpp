#include <cmath>
#include "Simulator2D.hpp"

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
    difuse();
    project();
    IFFT();
}

void Simulator2D::densityStep()
{
    addSource(m_grid_cells.dens, m_grid_cells.dens_prev);
    diffuse(m_grid_cells.dens, m_grid_cells.dens_prev, 0);
    advect(m_grid_cells.dens, m_grid_cells.dens_prev, 0);
}

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
    float x = LENGTH / (2 * N);
    float y = LENGTH / (2 * N);
    for (int i = 0; i < N; ++i, x += LENGTH / N)
    {
        for (int j = 0; j < N; ++j, y += LENGTH / N)
        {
            float x_prev = N * (x - DT * m_grid_cells.v_x_prev[POS(i, j)]) - 0.5;
            float y_prev = N * (y - DT * m_grid_cells.v_y_prev[POS(i, j)]) - 0.5;
            int grid_x0 = (int)std::floor(x_prev);
            int grid_y0 = (int)std::floor(y_prev);
            float s = x_prev - (float)grid_x0;
            float t = y_prev - (float)grid_y0;

            grid_x0 = (N + (grid_x0 % N)) % N;
            grid_y0 = (N + (grid_y0 % N)) % N;

            int grid_x1 = (grid_x0 + LENGTH) % N;
            int grid_y1 = (grid_y0 + LENGTH) % N;

            // grid interpolation - linear
            v_x[POS(i, j)] = (1 - s) * () +
                             s * ((1 - t) * u0[i1 + N * j0] + t * u0[i1 + N * j1]);
            v[i + N * j] = (1 - s) * ((1 - t) * v0[i0 + N * j0] + t * v0[i0 + N * j1]) +
                           s * ((1 - t) * v0[i1 + N * j0] + t * v0[i1 + N * j1]);
        }
    }
}

void Simulator2D::diffuse(float current[], float previous[], int b)
{
    float diff_rate = DT * DIFFUSION * LENGTH * LENGTH;
    for (int k = 0; k < ITERATIONS; ++k)
    {
        for (int i = 1; i < LENGTH; ++i)
        {
            for (int j = i; j < LENGTH; ++j)
            {
                current[POS(i, j)] = (previous[POS(i, j)] + diff_rate * calNeighborSum(current, i, j)) / (1 + 4 * diff_rate);
            }
        }
        setBounds(current, b);
    }
}
