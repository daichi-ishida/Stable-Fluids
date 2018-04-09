#include <cmath>
#include <random>
#include <iostream>
#include "Simulator2D.hpp"
#include "interpolation.hpp"

GridCells2D *Simulator2D::m_grid_cells;
bool Simulator2D::m_is_dragging;
glm::ivec2 Simulator2D::m_new_pos;
glm::ivec2 Simulator2D::m_old_pos;

void Simulator2D::mouseEvent(GLFWwindow *window, int button, int action, int mods)
{
    double px, py;
    glfwGetCursorPos(window, &px, &py);

    if (action == GLFW_PRESS)
    {
        // std::cout << "Dragging Start" << std::endl;

        if (!m_is_dragging)
        {
            m_is_dragging = true;
            m_old_pos = glm::ivec2(px, py);
            m_new_pos = glm::ivec2(px, py);
        }
    }
    else
    {
        m_is_dragging = false;
        m_old_pos = glm::ivec2(0, 0);
        m_new_pos = glm::ivec2(0, 0);
    }
}

void Simulator2D::mouseMoveEvent(GLFWwindow *window, double xpos, double ypos)
{
    if (m_is_dragging)
    {
        // std::cout << "move" << std::endl;
        // マウスの現在位置を更新
        m_new_pos = glm::ivec2(xpos, ypos);

        // マウスがあまり動いていない時は処理をしない
        const double dx = m_new_pos.x - m_old_pos.x;
        const double dy = m_new_pos.y - m_old_pos.y;
        const double length = dx * dx + dy * dy;
        if (length < 2.0f * 2.0f)
        {
            return;
        }
        else
        {
            double fx, fy;
            unsigned int i, j;

            // 加える力を概算
            fx = N * (m_new_pos.x - m_old_pos.x) / LENGTH;
            fy = N * (m_new_pos.y - m_old_pos.y) / LENGTH;

            // 力を加える
            i = std::fmax(0.0, std::fmin(N - 1 - 1e-6, N * m_new_pos.x / LENGTH));
            j = std::fmax(0.0, std::fmin(N - 1 - 1e-6, N * m_new_pos.y / LENGTH));
            if (i > 0 && j > 0 && i < N - 1 && j < N - 1)
            { // 端のグリッドは避ける
                // 重み関数を計算（概算なので，なくてもいい）
                double wx = N * m_new_pos.x / LENGTH - i;
                double wy = N * m_new_pos.y / LENGTH - j;

                // 直接流速に力を与える
                m_grid_cells->u0[POS(i, j)] += (1.0 - wx) * fx;
                m_grid_cells->u0[POS(i + 1, j)] += wx * fx;
                m_grid_cells->v0[POS(i, j)] += (1.0 - wy) * fy;
                m_grid_cells->v0[POS(i, j + 1)] += wy * fy;
            }
            m_old_pos = m_new_pos;
        }
    }
}

Simulator2D::Simulator2D(GridCells2D *grid_cells) : m_is_pause(false)
{
    m_grid_cells = grid_cells;
    m_fft_U = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * N * N);
    m_fft_V = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * N * N);
}

Simulator2D::~Simulator2D()
{
    fftwf_destroy_plan(m_plan_rc);
    fftwf_destroy_plan(m_plan_cr);
    fftwf_free(m_fft_U);
    fftwf_free(m_fft_V);
}

void Simulator2D::update()
{
    if (m_is_pause || !m_grid_cells)
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
    for (int i = 0; i < N * N; ++i)
    {
        m_grid_cells->u[i] += DT * m_grid_cells->u0[i];
        m_grid_cells->v[i] += DT * m_grid_cells->v0[i];
    }
}

void Simulator2D::advect()
{
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            float x = constrainValue(i - N * DT * m_grid_cells->u0[POS(i, j)]);
            float y = constrainValue(j - N * DT * m_grid_cells->v0[POS(i, j)]);

            int i0 = (int)std::floor(x);
            int j0 = (int)std::floor(y);
            float s = x - (float)i0;
            float t = y - (float)j0;

            int i1 = i0 + 1;
            int j1 = j0 + 1;

            // grid interpolation - linear
            float intrpl_u_vert = lerp::linear(m_grid_cells->u0[POS(i0, j0)], m_grid_cells->u0[POS(i0, j1)], t);
            float intrpl_u_horiz = lerp::linear(m_grid_cells->u0[POS(i1, j0)], m_grid_cells->u0[POS(i1, j1)], t);

            float intrpl_v_vert = lerp::linear(m_grid_cells->v0[POS(i0, j0)], m_grid_cells->v0[POS(i0, j1)], t);
            float intrpl_v_horiz = lerp::linear(m_grid_cells->v0[POS(i1, j0)], m_grid_cells->v0[POS(i1, j1)], t);

            m_grid_cells->u[POS(i, j)] = lerp::linear(intrpl_u_vert, intrpl_u_horiz, s);
            m_grid_cells->v[POS(i, j)] = lerp::linear(intrpl_v_vert, intrpl_v_horiz, s);
        }
    }
}

void Simulator2D::FFT()
{
    // prepare arrays for FFTW
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            m_grid_cells->u0[POS(i, j)] = m_grid_cells->u[POS(i, j)];
            m_grid_cells->v0[POS(i, j)] = m_grid_cells->v[POS(i, j)];
        }
    }

    m_plan_rc = fftwf_plan_dft_r2c_2d(N, N, m_grid_cells->u0, m_fft_U, FFTW_ESTIMATE);
    fftwf_execute(m_plan_rc);
    m_plan_rc = fftwf_plan_dft_r2c_2d(N, N, m_grid_cells->v0, m_fft_V, FFTW_ESTIMATE);
    fftwf_execute(m_plan_rc);
}

void Simulator2D::diffuse()
{
    // damp viscosity and conserve mass
    // in fourier space
    for (int i = 0; i < N; ++i)
    {
        float x = 0.5 * i;
        for (int j = 0; j < N; ++j)
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
    m_plan_cr = fftwf_plan_dft_c2r_2d(N, N, m_fft_U, m_grid_cells->u0, FFTW_ESTIMATE);
    fftwf_execute(m_plan_cr);
    m_plan_cr = fftwf_plan_dft_c2r_2d(N, N, m_fft_V, m_grid_cells->v0, FFTW_ESTIMATE);
    fftwf_execute(m_plan_cr);

    // normalize
    float f = 1.0 / (N * N);
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            m_grid_cells->u[POS(i, j)] = f * m_grid_cells->u0[POS(i, j)];
            m_grid_cells->v[POS(i, j)] = f * m_grid_cells->v0[POS(i, j)];
        }
    }
}

void Simulator2D::addSource()
{
    //for (int j = 0; j < N / 30 + 1; ++j)
    for (int j = N / 2 - 7; j < N / 2 + 7; ++j)
    {
        // initialize smoke
        for (int i = N / 2 - 7; i < N / 2 + 7; ++i)
        {
            m_grid_cells->dens[POS(i, j)] = 1.0f;
        }
        // initialize temp
        for (int i = N / 2 - 5; i < N / 2 + 5; ++i)
        {

            static std::mt19937 mt64(0);

            std::uniform_real_distribution<float> f_uni_rand(0.0f, 1.0f);

            m_grid_cells->temp[POS(i, j)] = std::sqrt(f_uni_rand(mt64)) * TEMPERATURE_MAX;
        }
    }
}

void Simulator2D::resetForce()
{
    for (int j = 0; j < N; ++j)
    {
        for (int i = 0; i < N; ++i)
        {
            const float weight = GRAVITY_Y;
            const float bouy = 1.0f;

            float f_x = 0.0f;
            float f_y = weight * m_grid_cells->dens[POS(i, j)] + bouy * (m_grid_cells->temp[POS(i, j)] - AMBIENT_TEMP);

            m_grid_cells->u0[POS(i, j)] = DT * f_x;
            m_grid_cells->v0[POS(i, j)] = DT * f_y;
        }
    }
}

void Simulator2D::calVorticity()
{
    Eigen::Vector3f eta;

    for (int j = 0; j < N; ++j)
    {
        for (int i = 0; i < N; ++i)
        {
            int i0 = (i - 1 + N) % N;
            int j0 = (j - 1 + N) % N;
            int i1 = (i + 1) % N;
            int j1 = (j + 1) % N;

            vortg[POS(i, j)][0] = 0.0;
            vortg[POS(i, j)][1] = 0.0;
            vortg[POS(i, j)][2] = (m_grid_cells->v[POS(i1, j)] - m_grid_cells->v[POS(i0, j)] - m_grid_cells->u[POS(i, j1)] + m_grid_cells->u[POS(i, j0)]) * 0.5 * N / LENGTH;
        }
    }
    for (int j = 0; j < N; ++j)
    {
        for (int i = 0; i < N; ++i)
        {
            int i0 = (i - 1 + N) % N;
            int j0 = (j - 1 + N) % N;
            int i1 = (i + 1) % N;
            int j1 = (j + 1) % N;

            eta[0] = (vortg[POS(i1, j)].norm() - vortg[POS(i0, j)].norm()) * 0.5 * N / LENGTH;
            eta[1] = (vortg[POS(i, j1)].norm() - vortg[POS(i, j0)].norm()) * 0.5 * N / LENGTH;
            eta[2] = 0.0f;
            if (eta.norm() != 0.0)
            {
                eta.normalize();
            }

            Eigen::Vector3f f = VORT_EPS * (LENGTH / N) * (eta.cross(vortg[POS(i, j)]));

            m_grid_cells->u0[POS(i, j)] += DT * f[0];
            m_grid_cells->v0[POS(i, j)] += DT * f[1];
        }
    }
}

void Simulator2D::advectDensity()
{
    for (int j = 0; j < N; ++j)
    {
        for (int i = 0; i < N; ++i)
        {
            float x = constrainValue(i - N * DT * m_grid_cells->u0[POS(i, j)]);
            float y = constrainValue(j - N * DT * m_grid_cells->v0[POS(i, j)]);

            int i0 = (int)std::floor(x);
            int j0 = (int)std::floor(y);
            float s = x - (float)i0;
            float t = y - (float)j0;

            int i1 = i0 + 1;
            int j1 = j0 + 1;

            // grid interpolation - linear
            float intrpl_dens_vert = lerp::linear(m_grid_cells->dens[POS(i0, j0)], m_grid_cells->dens[POS(i0, j1)], t);
            float intrpl_dens_horiz = lerp::linear(m_grid_cells->dens[POS(i1, j0)], m_grid_cells->dens[POS(i1, j1)], t);

            float intrpl_temp_vert = lerp::linear(m_grid_cells->temp[POS(i0, j0)], m_grid_cells->temp[POS(i0, j1)], t);
            float intrpl_temp_horiz = lerp::linear(m_grid_cells->temp[POS(i1, j0)], m_grid_cells->temp[POS(i1, j1)], t);

            m_grid_cells->dens[POS(i, j)] = lerp::linear(intrpl_dens_vert, intrpl_dens_horiz, s);
            m_grid_cells->temp[POS(i, j)] = lerp::linear(intrpl_temp_vert, intrpl_temp_horiz, s);
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