#pragma once
#include <fftw3.h>

#define GLFW_INCLUDE_GLU
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>

#include "GridCells2D.hpp"

class Simulator2D
{
public:
  Simulator2D(GridCells2D *grid_cells, EMode mode);
  ~Simulator2D();

  void update();
  // void pause() { m_is_pause = true; };
  // void restart() { m_is_pause = false; };

  static void mouseEvent(GLFWwindow *window, int button, int action, int mods);
  static void mouseMoveEvent(GLFWwindow *window, double xpos, double ypos);

private:
  void velocityStep();
  void densityStep();

  void addForce();
  void advect();
  void FFT();
  void diffuse();
  void IFFT();

  void addSource();
  void resetForce();
  void advectDensity();

  float interp(float x, float y, float q[], unsigned int Nx, unsigned int Ny);

  static GridCells2D *m_grid_cells;

  fftwf_plan m_plan_u_rc;
  fftwf_plan m_plan_u_cr;
  fftwf_plan m_plan_v_rc;
  fftwf_plan m_plan_v_cr;

  fftwf_complex *m_fft_U;
  fftwf_complex *m_fft_V;

  EMode m_mode;
  bool m_is_pause;
  static bool m_is_dragging;
  static glm::ivec2 m_new_pos;
  static glm::ivec2 m_old_pos;
};