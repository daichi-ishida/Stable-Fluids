#include <iostream>
#include "constants.hpp"
#include "GridCells.hpp"
#include "Simulator.hpp"

int main()
{
    GridCells grid_cells;
    Simulator simulator(&grid_cells);

    while(1)
    {
        simulator.update();
    }
    return 0;
}

