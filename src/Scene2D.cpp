#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>

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
    glClear(GL_COLOR_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glViewport(0, 0, WIDTH, HEIGHT);
    glLoadIdentity();
    glOrtho(0, WIDTH, HEIGHT, 0, -1.0, 1.0);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    if (!m_grid_cells)
    {
        return;
    }

    drawDensity();

    glColor4f(1.0f, 1.0f, 1.0f, 0.5f);
    drawGrid();

    glColor4f(1.0f, 0.0f, 0.0f, 0.5f);
    drawVelocity();
}

/* private */
void Scene2D::drawGrid()
{
    glBegin(GL_LINES);
    for (unsigned int i = 0; i < N + 1; ++i)
    {
        glVertex2d((WIDTH / (float)N) * i, 0.0);
        glVertex2d((WIDTH / (float)N) * i, HEIGHT);
    }
    for (unsigned int j = 0; j < N + 1; ++j)
    {
        glVertex2d(0.0, (HEIGHT / (float)N) * j);
        glVertex2d(WIDTH, (HEIGHT / (float)N) * j);
    }
    glEnd();
}

void Scene2D::drawVelocity()
{
    glBegin(GL_LINES);
    for (unsigned int y = 0; y < N; ++y)
    {
        for (unsigned int x = 0; x < N; ++x)
        {
            glm::vec2 p = {(x + 0.5) * WIDTH / (float)N, (y + 0.5) * HEIGHT / (float)N};
            glm::vec2 vel = {m_grid_cells->u[POS(x, y)], m_grid_cells->v[POS(x, y)]};
            vel = glm::normalize(vel);
            float ks = 1.0f;
            glVertex2d(p.x, p.y);
            glVertex2d(p.x + ks * (WIDTH / (float)N) * vel.x, p.y + ks * (HEIGHT / (float)N) * vel.y);
        }
    }
    glEnd();
}

void Scene2D::drawDensity()
{
    for (unsigned int y = 0; y < N; ++y)
    {
        for (unsigned int x = 0; x < N; ++x)
        {
            glColor4d(1.0f, 1.0f, 1.0f, m_grid_cells->dens[POS(x, y)]);
            glRectf(x * WIDTH / (float)N, y * HEIGHT / (float)N, (x + 1) * WIDTH / (float)N, (y + 1) * HEIGHT / (float)N);
        }
    }
}

void Scene2D::writeData()
{
    writeData_inVtiFormat();
    ++m_file_num;
}

void Scene2D::writeData_inVtiFormat()
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