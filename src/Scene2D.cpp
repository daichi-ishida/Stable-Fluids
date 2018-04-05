#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip> 
#include <string>
#include <direct.h>
#include "Scene2D.hpp"
#include "constants.hpp"

Scene2D::Scene2D(GridCells2D &grid_cells, const float &time) : m_grid_cells(grid_cells), m_time(time),m_file_num(0)
{
}

Scene2D::~Scene2D()
{
}

void Scene2D::writeData()
{
    _mkdir("output");
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
        std::cout << "ERROR : file open error at writing data in .vti format\n" << file_name << " cannot open" << std::endl;
        exit(EXIT_FAILURE);
    }

    /* header */
    ofs << "<?xml version='1.0' encoding='UTF-8'?>" << std::endl;
    ofs << "<VTKFile xmlns='VTK' byte_order='LittleEndian' version='0.1' type='ImageData'>" << std::endl;
    ofs << "<ImageData WholeExtent='1 " <<  N+2 << " 1 " << N+2 << " 0 1' Origin='0 0 0” Spacing='1.0 1.0 1.0'>" << std::endl;
    ofs << "<Piece Extent='1 " <<  N+2 << " 1 " << N+2 << " 0 1' Origin='0 0 0' Spacing='1.0 1.0 1.0'>" << std::endl;
    ofs << "<PointData Scalars='density'>" << std::endl;
    ofs << "<DataArray Name='density' type='Float32' format='ascii'>" << std::endl;

    for (int i = 0; i < N + 2; ++i)
    {
        for (int j = 0; j < N + 2;++j)
        {
            ofs << m_grid_cells.dens[POS(i, j)] << " ";
        }
        ofs << std::endl;
    }

    ofs << "</DataArray>" << std::endl;
    ofs << "</PointData>" << std::endl;
    ofs << "</Piece>" << std::endl;
    ofs << "</ImageData>" << std::endl;
    ofs << "</VTKFile>" << std::endl;

    ofs.close();
}