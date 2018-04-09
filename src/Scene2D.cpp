#define GLM_ENABLE_EXPERIMENTAL
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <GL/gl.h>
#include <glm/glm.hpp>

#include "Scene2D.hpp"
#include "constants.hpp"

Scene2D::Scene2D(GridCells2D *grid_cells, const float &time) : m_grid_cells(nullptr), m_time(time), m_file_num(0)
{
    m_grid_cells = grid_cells;
}

Scene2D::~Scene2D()
{
}

void Scene2D::update()
{
}

void Scene2D::draw()
{
    // 背景色と深度値のクリア
    glClear(GL_COLOR_BUFFER_BIT);
    // 見る範囲の指定
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glOrtho(-0.02, LENGTH + 0.02, -0.02, LENGTH + 0.02, -1.0, 1.0);

    if (!m_grid_cells)
    {
        return;
    }
    drawDensity();

    glColor4f(0.0f, 1.0f, 0.0f, 0.5f);
    drawGrid();

    glColor4f(0.0f, 0.0f, 1.0f, 0.5f);
    drawVelocity();
}

/* private */
void Scene2D::drawGrid()
{
    glBegin(GL_LINES);
    for (unsigned int i = 0; i < N + 1; i++)
    {
        glVertex2d(LENGTH / N * i, 0.0);
        glVertex2d(LENGTH / N * i, LENGTH);
    }
    for (unsigned int j = 0; j < N + 1; j++)
    {
        glVertex2d(0.0, LENGTH / N * j);
        glVertex2d(LENGTH, LENGTH / N * j);
    }
    glEnd();
}

void Scene2D::drawVelocity()
{
    glBegin(GL_LINES);
    for (unsigned int i = 0; i < N; i++)
        for (unsigned int j = 0; j < N; j++)
        {
            // 中心の流速を可視化する
            glm::vec2 p = {(i + 0.5) * LENGTH / N, (j + 0.5) * LENGTH / N};
            glm::vec2 vel = {0.5 * (m_grid_cells->u[POS(i, j)] + m_grid_cells->u[POS(i + 1, j)]), 0.5 * (m_grid_cells->v[POS(i, j)] + m_grid_cells->v[POS(i, j + 1)])};
            double ks = 5.0;
            glVertex2d(p.x, p.y);
            glVertex2d(p.x + ks * LENGTH / N * vel.x, p.y + ks * LENGTH / N * vel.y);
        }
    glEnd();
}

void Scene2D::drawDensity()
{
    double maxv = -1e18;
    double minv = 1e18;
    for (unsigned int i = 0; i < N; i++)
    {
        for (unsigned int j = 0; j < N; j++)
        {
            maxv = std::fmax(maxv, m_grid_cells->dens[POS(i, j)]);
            minv = std::fmin(minv, m_grid_cells->dens[POS(i, j)]);
        }
    }

    double det = maxv - minv;
    if (det)
    {
        for (unsigned int i = 0; i < N; i++)
        {
            for (unsigned int j = 0; j < N; j++)
            {
                // 中心の流速を計算する
                double normp = 2.0 * (m_grid_cells->dens[POS(i, j)] - minv) / det - 1.0;
                glColor4d(normp > 0, 0.3, normp <= 0, std::fabs(normp));
                glRectf(i * LENGTH / N, j * LENGTH / N, (i + 1) * LENGTH / N, (j + 1) * LENGTH / N);
            }
        }
    }
}

void Scene2D::writeData()
{
    // _mkdir("output");
    writeData_inVtuFormat();
    ++m_file_num;
}

void Scene2D::writeData_inVtuFormat()
{
    std::ostringstream sout;
    sout << std::setfill('0') << std::setw(3) << std::right << m_file_num;

    std::string file_name = "output/stable_fluids_" + sout.str() + ".vti";
    std::ofstream ofs;
    ofs.open(file_name);
    if (!ofs)
    {
        std::cout << "ERROR : file open error at writing data in .vti format\n"
                  << file_name << " cannot open" << std::endl;
        exit(EXIT_FAILURE);
    }

    /* header */
    ofs << "<?xml version='1.0' encoding='UTF-8'?>" << std::endl;
    ofs << "<VTKFile xmlns='VTK' byte_order='LittleEndian' version='0.1' type='ImageData'>" << std::endl;
    ofs << "<ImageData WholeExtent='0 " << N - 1 << " 0 " << N - 1 << " 0 0' Origin='0 0 0' Spacing='1.0 1.0 1.0'>" << std::endl;
    ofs << "<Piece Extent='0 " << N - 1 << " 0 " << N - 1 << " 0 0'>" << std::endl;
    ofs << "<PointData Scalars='density' Vectors='velocity'>" << std::endl;

    ofs << "<DataArray type='Float32' Name='density' NumberOfComponents='1' format='ascii'>" << std::endl;
    for (int j = 0; j < N; ++j)
    {
        for (int i = 0; i < N; ++i)
        {
            ofs << m_grid_cells->dens[POS(i, j)] << " ";
        }
        ofs << std::endl;
    }
    ofs << "</DataArray>" << std::endl;

    ofs << "<DataArray type='Float32' Name='velocity' NumberOfComponents='3' format='ascii'>" << std::endl;

    for (int j = 0; j < N; ++j)
    {
        for (int i = 0; i < N; ++i)
        {
            ofs << m_grid_cells->u[POS(i, j)] << " " << m_grid_cells->v[POS(i, j)] << " 0" << std::endl;
        }
    }
    ofs << "</DataArray>" << std::endl;

    ofs << "</PointData>" << std::endl;
    ofs << "<CellData></CellData>" << std::endl;
    ofs << "</Piece>" << std::endl;
    ofs << "</ImageData>" << std::endl;
    ofs << "</VTKFile>" << std::endl;

    ofs.close();
}