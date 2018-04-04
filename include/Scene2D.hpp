#pragma once
#include "GridCells2D.hpp"

class Scene2D
{
public:
  Scene2D(GridCells2D &grid_cells);
  ~Scene2D();

  void update();
  void draw();

private:
  GridCells2D &m_grid_cells;
};