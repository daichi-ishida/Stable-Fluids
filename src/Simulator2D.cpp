#include <cmath>
#include <random>
#include "Simulator2D.hpp"
#include "interpolation.hpp"

Simulator2D::Simulator2D(GridCells2D &grid_cells) : m_grid_cells(grid_cells),
                                                    m_is_pause(false),
                                                    m_use_vor_particles(USE_VORTEX_PARTICLES)
{
    m_fft_U = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * (N+2) * (N+2));
    m_fft_V = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * (N+2) * (N+2));
}

Simulator2D::~Simulator2D()
{
}

void Simulator2D::update()
{
    if (m_is_pause)
    {
        return;
    }

    velocityStep();
    densityStep();
}

/* private */

void Simulator2D::velocityStep()
{
    addForce();
    advect();
    FFT();
    diffuse();
    IFFT();
}

void Simulator2D::densityStep()
{
    addSource();
    resetForce();
    calVorticity();
    advectDensity();
}

void Simulator2D::addForce()
{
    for (int i = 0; i < SIZE; ++i)
    {
        m_grid_cells.u[i] += DT * m_grid_cells.u0[i];
        m_grid_cells.v[i] += DT * m_grid_cells.v0[i];
    }
}

void Simulator2D::advect()
{
    for (int i = 1; i < N + 1; ++i)
    {
        for (int j = 1; j < N + 1; ++j)
        {
            float x = constrainValue(i - N * DT * m_grid_cells.u0[POS(i, j)]);
            float y = constrainValue(j - N * DT * m_grid_cells.v0[POS(i, j)]);

            int i0 = (int)std::floor(x);
            int j0 = (int)std::floor(y);
            float s = x - (float)i0;
            float t = y - (float)j0;

            int i1 = i0 + 1;
            int j1 = j0 + 1;

            // grid interpolation - linear
            float intrpl_u_vert = lerp::linear(m_grid_cells.u0[POS(i0, j0)], m_grid_cells.u0[POS(i0, j1)], t);
            float intrpl_u_horiz = lerp::linear(m_grid_cells.u0[POS(i1, j0)], m_grid_cells.u0[POS(i1, j1)], t);

            float intrpl_v_vert = lerp::linear(m_grid_cells.v0[POS(i0, j0)], m_grid_cells.v0[POS(i0, j1)], t);
            float intrpl_v_horiz = lerp::linear(m_grid_cells.v0[POS(i1, j0)], m_grid_cells.v0[POS(i1, j1)], t);

            m_grid_cells.u[POS(i, j)] = lerp::linear(intrpl_u_vert, intrpl_u_horiz, s);
            m_grid_cells.v[POS(i, j)] = lerp::linear(intrpl_v_vert, intrpl_v_horiz, s);
        }
    }
}

void Simulator2D::FFT()
{
    // prepare arrays for FFTW
    for (int i = 0; i < N + 2; ++i)
    {
        for (int j = 0; j < N + 2; ++j)
        {
            m_grid_cells.u0[i] = m_grid_cells.u[i];
            m_grid_cells.v0[i] = m_grid_cells.v[i];
        }
    }

    m_plan_rc = fftwf_plan_dft_r2c_2d(N, N, m_grid_cells.u0, m_fft_U, FFTW_ESTIMATE);
    fftwf_execute(m_plan_rc);
    m_plan_rc = fftwf_plan_dft_r2c_2d(N, N, m_grid_cells.v0, m_fft_V, FFTW_ESTIMATE);
    fftwf_execute(m_plan_rc);
}

void Simulator2D::diffuse()
{
    // damp viscosity and conserve mass
    // in fourier space
    for (int i = 0; i < N + 2; ++i)
    {
        float x = 0.5 * i;
        for (int j = 0; j < N + 2; ++j)
        {
            float y = (j <= N / 2) ? j : j - N;
            float r = x * x + y * y;
            if (r == 0.0)
            {
                continue;
            }
            float f = std::exp(-r * DT * VISCOSITY);
            float U0 = m_fft_U[POS(i, j)][0];
            float V0 = m_fft_V[POS(i, j)][0];

            float U1 = m_fft_U[POS(i, j)][1];
            float V1 = m_fft_V[POS(i, j)][1];

            m_fft_U[POS(i, j)][0] = f * ((1 - x * x / r) * U0 - x * y / r * V0);
            m_fft_U[POS(i, j)][1] = f * ((1 - x * x / r) * U1 - x * y / r * V1);
            m_fft_V[POS(i, j)][0] = f * (-y * x / r * U0 + (1 - y * y / r) * V0);
            m_fft_V[POS(i, j)][1] = f * (-y * x / r * U1 + (1 - y * y / r) * V1);
        }
    }
}

void Simulator2D::IFFT()
{
    m_plan_cr = fftwf_plan_dft_c2r_2d(N, N, m_fft_U, m_grid_cells.u0, FFTW_ESTIMATE);
    fftwf_execute(m_plan_cr);
    m_plan_cr = fftwf_plan_dft_c2r_2d(N, N, m_fft_V, m_grid_cells.v0, FFTW_ESTIMATE);
    fftwf_execute(m_plan_cr);

    // normalize
    float f = 1.0 / (N * N);
    for (int i = 0; i < N + 2; ++i)
    {
        for (int j = 0; j < N + 2; ++j)
        {
            m_grid_cells.u[POS(i, j)] = f * m_grid_cells.u0[POS(i, j)];
            m_grid_cells.v[POS(i, j)] = f * m_grid_cells.v0[POS(i, j)];
        }
    }
}

void Simulator2D::addSource()
{
    for (int j = 0; j < N / 30 + 1; ++j)
    {
        // initialize smoke
        for (int i = N / 2 - 7; i < N / 2 + 7; ++i)
        {
            m_grid_cells.dens[POS(i, j)];
        }
        // initialize temp
        for (int i = N / 2 - 5; i < N / 2 + 5; ++i)
        {
            
            static std::mt19937 mt64(0);

            std::uniform_real_distribution<float> f_uni_rand(0.0, 1.0);

            m_grid_cells.temp[POS(i, j)] = std::sqrt(f_uni_rand(mt64)) * 200.0;
        }
    }
}

void Simulator2D::resetForce()
{
    for (int i = 0; i < N + 2; ++i)
    {
        for (int j = 0; j < N + 2; ++j)
        {
            const float weight = GRAVITY_Y;
            const float bouy = 1.0;

            float f_x = 0.0f;
            float f_y = weight * m_grid_cells.dens[POS(i, j)] + bouy * (m_grid_cells.temp[POS(i, j)] - AMBIENT_TEMP);

            m_grid_cells.u0[POS(i, j)] = DT * f_x;
            m_grid_cells.v0[POS(i, j)] = DT * f_y;
        }
    }
}

void Simulator2D::calVorticity()
{
    Eigen::Vector3f eta;

    for (int i = 0; i < N + 2; ++i)
    {
        for (int j = 0; j < N + 2; ++j)
        {
            int i0 = (i - 1 + N) % (N + 2);
            int j0 = (j - 1 + N) % (N + 2);
            int i1 = (i + 1) % (N + 2);
            int j1 = (j + 1) % (N + 2);

            vortg[POS(i, j)][0] = 0.0;
            vortg[POS(i, j)][1] = 0.0;
            vortg[POS(i, j)][2] = (m_grid_cells.u[POS(i1, j)] - m_grid_cells.u[POS(i0, j)] - m_grid_cells.v[POS(i, j1)] + m_grid_cells.v[POS(i, j0)]) * 0.5 * N / LENGTH;
        }
    }
    for (int i = 0; i < N + 2; ++i)
    {
        for (int j = 0; j < N + 2; ++j)
        {
            int i0 = (i - 1 + N) % (N + 2);
            int j0 = (j - 1 + N) % (N + 2);
            int i1 = (i + 1) % (N + 2);
            int j1 = (j + 1) % (N + 2);

            eta[0] = (vortg[POS(i1, j)].norm() - vortg[POS(i0, j)].norm()) * 0.5 * N / LENGTH;
            eta[1] = (vortg[POS(i, j1)].norm() - vortg[POS(i, j0)].norm()) * 0.5 * N / LENGTH;
            eta[2] = 0.0f;
            eta.normalize();

            Eigen::Vector3f f = VORT_EPS * (LENGTH / N) * (eta.cross(vortg[POS(i, j)]));

            m_grid_cells.u0[POS(i, j)] += DT * f[0];
            m_grid_cells.v0[POS(i, j)] += DT * f[1];
        }
    }

    // else
    // {
    //     for (int i = 0; i < N + 2; ++i)
    //     {
    //         for (int j = 0; j < N + 2; ++j)
    //         {
    //             int i0 = (i - 1 + N) % (N + 2);
    //             int j0 = (j - 1 + N) % (N + 2);
    //             int i1 = (i + 1) % (N + 2);
    //             int j1 = (j + 1) % (N + 2);

    //             du[POS(i, j)][0] = (m_grid_cells.u[POS(i1, j)] - m_grid_cells.u[POS(i0, j)]) * 0.5 * N / LENGTH;
    //             du[POS(i, j)][1] = (m_grid_cells.u[POS(i, j1)] - m_grid_cells.u[POS(i, j0)]) * 0.5 * N / LENGTH;
    //             du[POS(i, j)][2] = 0.0; // for 2D

    //             dv[POS(i, j)][0] = (m_grid_cells.v[POS(i1, j)] - m_grid_cells.v[POS(i0, j)]) * 0.5 * N / LENGTH;
    //             dv[POS(i, j)][1] = (m_grid_cells.v[POS(i, j1)] - m_grid_cells.v[POS(i, j0)]) * 0.5 * N / LENGTH;
    //             dv[POS(i, j)][2] = 0.0; // for 2D
    //         }
    //     }
    // }
}

void Simulator2D::advectDensity()
{
    for (int i = 1; i < N + 1; ++i)
    {
        for (int j = 1; j < N + 1; ++j)
        {
            float x = constrainValue(i - N * DT * m_grid_cells.u0[POS(i, j)]);
            float y = constrainValue(j - N * DT * m_grid_cells.v0[POS(i, j)]);

            int i0 = (int)std::floor(x);
            int j0 = (int)std::floor(y);
            float s = x - (float)i0;
            float t = y - (float)j0;

            int i1 = i0 + 1;
            int j1 = j0 + 1;

            // grid interpolation - linear
            float intrpl_dens_vert = lerp::linear(m_grid_cells.dens[POS(i0, j0)], m_grid_cells.dens[POS(i0, j1)], t);
            float intrpl_dens_horiz = lerp::linear(m_grid_cells.dens[POS(i1, j0)], m_grid_cells.dens[POS(i1, j1)], t);

            float intrpl_temp_vert = lerp::linear(m_grid_cells.temp[POS(i0, j0)], m_grid_cells.temp[POS(i0, j1)], t);
            float intrpl_temp_horiz = lerp::linear(m_grid_cells.temp[POS(i1, j0)], m_grid_cells.temp[POS(i1, j1)], t);

            m_grid_cells.dens[POS(i, j)] = lerp::linear(intrpl_dens_vert, intrpl_dens_horiz, s);
            m_grid_cells.temp[POS(i, j)] = lerp::linear(intrpl_temp_vert, intrpl_temp_horiz, s);
        }
    }
}

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