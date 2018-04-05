#include <iostream>
#include "constants.hpp"
#include "GridCells2D.hpp"
#include "Scene2D.hpp"
#include "Simulator2D.hpp"

int main()
{
    float time = 0.0f;
    GridCells2D grid_cells;
    Scene2D scene(grid_cells, time);
    Simulator2D simulator(grid_cells);

    std::cout << "\n*** START SIMULATION ***\n";
    scene.writeData();
    while (1)
    {
        time += DT;
        simulator.update();
        scene.writeData();
        if (time >= FINISH_TIME)
        {
            break;
        }
    }
    std::cout << "\n*** END SIMULATION ***\n";

    return 0;
}
