#include <iostream>
#include "constants.hpp"
#include "GridCells2D.hpp"
#include "Simulator2D.hpp"

int main()
{
    GridCells2D grid_cells;
    Simulator2D simulator(grid_cells);

    while(1)
    {
        simulator.update();
    }
    return 0;
}

