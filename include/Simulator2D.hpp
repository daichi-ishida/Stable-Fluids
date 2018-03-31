#pragma once
#include "GridCells2D.hpp"

class Simulator2D
{
public:
  Simulator2D(GridCells2D &grid_cells);
  ~Simulator2D();

  void update();

private:
  void swap();

  void velocityStep();
  void densityStep();

  void addForce();
  void advect();
  void FFT();
  void diffuse();
  void project();
  void IFFT();

  GridCells2D &m_grid_cells;
};