#pragma once
#include <fftw3.h>
#include <Eigen/Core>
#include "GridCells2D.hpp"

class Simulator2D
{
public:
  Simulator2D(GridCells2D &grid_cells);
  ~Simulator2D();

  void update();
  void pause() { m_is_pause = true; };
  void restart() { m_is_pause = false; };
  float getTime() { return m_time; };

private:
  void swap();

  void velocityStep();
  void densityStep();

  void addForce();
  void advect();
  void FFT();
  void diffuse();
  void IFFT();

  void addSource();
  void resetForce();
  void calVorticity();

  float constrainValue(float value);

  GridCells2D &m_grid_cells;
  fftwf_plan m_plan_rc;
  fftwf_plan m_plan_cr;
  fftwf_complex *m_fft_Vx;
  fftwf_complex *m_fft_Vy;

  Eigen::Vector3f vortg[SIZE];

  bool m_is_pause;
  bool m_use_vor_particles;
  float m_time;
};