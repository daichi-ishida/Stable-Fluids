#include <cmath>
#include <iostream>
#include "Simulator2D.hpp"

GridCells2D *Simulator2D::m_grid_cells;
bool Simulator2D::m_is_dragging;
glm::ivec2 Simulator2D::m_new_pos;
glm::ivec2 Simulator2D::m_old_pos;

Simulator2D::Simulator2D(GridCells2D *grid_cells, EMode mode) : m_is_pause(false), m_mode(mode)
{
    m_grid_cells = grid_cells;
    m_fft_U = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * N * N);
    m_fft_V = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * N * N);
    m_plan_u_rc = fftwf_plan_dft_r2c_2d(N, N, m_grid_cells->u0, m_fft_U, FFTW_MEASURE);
    m_plan_v_rc = fftwf_plan_dft_r2c_2d(N, N, m_grid_cells->v0, m_fft_V, FFTW_MEASURE);

    m_plan_u_cr = fftwf_plan_dft_c2r_2d(N, N, m_fft_U, m_grid_cells->u0, FFTW_MEASURE);
    m_plan_v_cr = fftwf_plan_dft_c2r_2d(N, N, m_fft_V, m_grid_cells->v0, FFTW_MEASURE);

    if (m_mode == E_Once)
    {
        addSource();
    }
}

Simulator2D::~Simulator2D()
{
    fftwf_destroy_plan(m_plan_u_rc);
    fftwf_destroy_plan(m_plan_v_rc);
    fftwf_destroy_plan(m_plan_u_cr);
    fftwf_destroy_plan(m_plan_v_cr);
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

/* callback function */
void Simulator2D::mouseEvent(GLFWwindow *window, int button, int action, int mods)
{
    if (button == GLFW_MOUSE_BUTTON_LEFT)
    {
        double px, py;
        glfwGetCursorPos(window, &px, &py);

        switch (action)
        {
        case GLFW_PRESS:
            m_is_dragging = true;
            m_old_pos = glm::ivec2(px, py);
            m_new_pos = glm::ivec2(px, py);
            break;
        default:
            m_is_dragging = false;
            m_old_pos = glm::ivec2(0, 0);
            m_new_pos = glm::ivec2(0, 0);
            break;
        }
    }
}

void Simulator2D::mouseMoveEvent(GLFWwindow *window, double xpos, double ypos)
{
    if (m_is_dragging)
    {
        // update mouse position
        m_new_pos = glm::ivec2(xpos, ypos);

        // ignore slight movement
        float dx = m_new_pos.x - m_old_pos.x;
        float dy = m_new_pos.y - m_old_pos.y;
        float length = dx * dx + dy * dy;
        if (length < 2.0f)
        {
            return;
        }
        else
        {
            float tmp_fx, tmp_fy;
            unsigned int i, j;

            // calculate force
            tmp_fx = INTERACTION * N * (m_new_pos.x - m_old_pos.x) / (float)WIDTH;
            tmp_fy = INTERACTION * N * (m_new_pos.y - m_old_pos.y) / (float)HEIGHT;

            // specify the index to add force
            i = std::fmax(0.0, std::fmin(N - 1, N * m_new_pos.x / (float)WIDTH));
            j = std::fmax(0.0, std::fmin(N - 1, N * m_new_pos.y / (float)HEIGHT));
            if (i > 0 && j > 0 && i < N - 1 && j < N - 1)
            { // avoid edge of grid
                // calculate weight
                float wx = N * m_new_pos.x / (float)WIDTH - i;
                float wy = N * m_new_pos.y / (float)HEIGHT - j;

                // add force
                m_grid_cells->fx[POS(i, j)] = (1.0 - wx) * tmp_fx;
                m_grid_cells->fx[POS(i + 1, j)] = wx * tmp_fx;
                m_grid_cells->fy[POS(i, j)] = (1.0 - wy) * tmp_fy;
                m_grid_cells->fy[POS(i, j + 1)] = wy * tmp_fy;
            }
            m_old_pos = m_new_pos;
        }
    }
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
    if (m_mode == E_Continuous)
    {
        addSource();
    }
    resetForce();
    advectDensity();
}

void Simulator2D::addForce()
{
    for (int i = 0; i < SIZE; ++i)
    {
        m_grid_cells->u[i] += DT * m_grid_cells->fx[i];
        m_grid_cells->v[i] += DT * m_grid_cells->fy[i];

        m_grid_cells->u0[i] = m_grid_cells->u[i];
        m_grid_cells->v0[i] = m_grid_cells->v[i];
    }
}

void Simulator2D::advect()
{
    for (unsigned int j = 0; j < N; ++j)
    {
        for (unsigned int i = 1; i < N; ++i)
        {
            float x = i * WIDTH / (float)N;
            float y = (j + 0.5) * HEIGHT / (float)N;

            x = x - DT * interp(x, y - 0.5 * HEIGHT / (float)N, m_grid_cells->u0, N + 1, N);
            y = y - DT * interp(x - 0.5 * WIDTH / (float)N, y, m_grid_cells->v0, N, N + 1);

            m_grid_cells->u[POS(i, j)] = interp(x, y - 0.5 * HEIGHT / (float)N, m_grid_cells->u0, N + 1, N);
        }
    }

    for (unsigned int j = 1; j < N; ++j)
    {
        for (unsigned int i = 0; i < N; ++i)
        {
            float x = (i + 0.5) * WIDTH / (float)N;
            float y = j * HEIGHT / (float)N;

            x = x - DT * interp(x, y - 0.5 * HEIGHT / (float)N, m_grid_cells->u0, N + 1, N);
            y = y - DT * interp(x - 0.5 * WIDTH / (float)N, y, m_grid_cells->v0, N, N + 1);

            m_grid_cells->v[POS(i, j)] = interp(x - 0.5 * WIDTH / (float)N, y, m_grid_cells->v0, N, N + 1);
        }
    }
}

void Simulator2D::FFT()
{
    // prepare arrays for FFTW
    for (int j = 0; j < N; ++j)
    {
        for (int i = 0; i < N; ++i)
        {
            m_grid_cells->u0[POS(i, j)] = m_grid_cells->u[POS(i, j)];
            m_grid_cells->v0[POS(i, j)] = m_grid_cells->v[POS(i, j)];
        }
    }

    fftwf_execute(m_plan_u_rc);
    fftwf_execute(m_plan_v_rc);
}

void Simulator2D::diffuse()
{
    // damp viscosity and conserve mass
    // in fourier space
    for (int j = 0; j < N; ++j)
    {
        float ky = (j <= N / 2) ? j : j - N;
        for (int i = 0; i <= N / 2; ++i)
        {
            float kx = i;
            float kk = kx * kx + ky * ky;

            if (kk < 0.001)
            {
                continue;
            }
            float f = std::exp(-kk * DT * VISCOSITY);
            int idx = i + j * (N / 2 + 1);
            float U0 = m_fft_U[idx][0];
            float V0 = m_fft_V[idx][0];

            float U1 = m_fft_U[idx][1];
            float V1 = m_fft_V[idx][1];

            m_fft_U[idx][0] = f * ((1 - kx * kx / kk) * U0 - (kx * ky / kk) * V0);
            m_fft_U[idx][1] = f * ((1 - kx * kx / kk) * U1 - (kx * ky / kk) * V1);
            m_fft_V[idx][0] = f * ((-kx * ky / kk) * U0 + (1 - ky * ky / kk) * V0);
            m_fft_V[idx][1] = f * ((-kx * ky / kk) * U1 + (1 - ky * ky / kk) * V1);
        }
    }
}

void Simulator2D::IFFT()
{
    fftwf_execute(m_plan_u_cr);
    fftwf_execute(m_plan_v_cr);

    // normalize
    float f = 1.0 / (float)(N * N);
    for (int j = 0; j < N; ++j)
    {
        for (int i = 0; i < N; ++i)
        {
            m_grid_cells->u[POS(i, j)] = f * m_grid_cells->u0[POS(i, j)];
            m_grid_cells->v[POS(i, j)] = f * m_grid_cells->v0[POS(i, j)];
        }
    }
}

void Simulator2D::addSource()
{
    for (int j = N / 2 - SOURCE_SIZE / 2; j < N / 2 + SOURCE_SIZE / 2; ++j)
    {
        // initialize smoke
        for (int i = N / 2 - SOURCE_SIZE / 2; i < N / 2 + SOURCE_SIZE / 2; ++i)
        {
            m_grid_cells->dens[POS(i, j)] = 1.0f;
        }
    }
}

void Simulator2D::resetForce()
{
    for (int j = 0; j < N; ++j)
    {
        for (int i = 0; i < N; ++i)
        {
            m_grid_cells->fx[POS(i, j)] = 0.0f;
            m_grid_cells->fy[POS(i, j)] = GRAVITY_Y;
        }
    }
}

void Simulator2D::advectDensity()
{
    for (unsigned int j = 0; j < N; ++j)
    {
        for (unsigned int i = 0; i < N; ++i)
        {
            float x = i * WIDTH / (float)N;
            float y = j * HEIGHT / (float)N;

            x = x - DT * interp(x, y, m_grid_cells->u, N, N);
            y = y - DT * interp(x, y, m_grid_cells->v, N, N);

            m_grid_cells->dens[POS(i, j)] = interp(x, y, m_grid_cells->dens, N, N);
        }
    }
}

float Simulator2D::interp(float x, float y, float q[], unsigned int Nx, unsigned int Ny)
{
    x = std::fmax(0.0, std::fmin(Nx - 1 - 1e-6, N * x / (float)WIDTH));
    y = std::fmax(0.0, std::fmin(Ny - 1 - 1e-6, N * y / (float)HEIGHT));

    unsigned int i = x;
    unsigned int j = y;

    float f[4] = {q[POS(i, j)], q[POS(i, j + 1)], q[POS(i + 1, j)], q[POS(i + 1, j + 1)]};

    x = x - i;
    y = y - j;

    float c[4] = {(1.0f - x) * (1.0f - y), (1.0f - x) * y, x * (1.0f - y), x * y};

    return c[0] * f[0] + c[1] * f[1] + c[2] * f[2] + c[3] * f[3];
}