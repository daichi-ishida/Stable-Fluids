#pragma once
#include "GridCells2D.hpp"

class Scene2D
{
public:
  Scene2D(GridCells2D *grid_cells, const float &time);
  ~Scene2D();

  void update();
  void draw();

  void writeData();

private:
  void drawGrid();
  void drawVelocity();
  void drawDensity();
  void writeData_inVtuFormat();

  int m_file_num;

  const float &m_time;

  GridCells2D *m_grid_cells;
};