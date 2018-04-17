#define GLFW_INCLUDE_GLU
#include <GLFW/glfw3.h>

#include "constants.hpp"
#include "GridCells2D.hpp"
#include "Scene2D.hpp"
#include "Simulator2D.hpp"

int main(int argc, char *argv[])
{
    // set default density mode
    EMode mode = E_Continuous;

    if (argc > 2)
    {
        fprintf(stderr, "too much arguments\n");
        exit(EXIT_FAILURE);
    }
    else if (argc == 2)
    {
        char *p = argv[1];
        if (*p == '-')
        {
            p++;
            // set density mode
            switch (*p)
            {
            case 'o':
                mode = E_Once;
                break;
            case 'c':
                mode = E_Continuous;
                break;
            default:
                fprintf(stderr, "exceptional argument\n");
                exit(EXIT_FAILURE);
                break;
            }
        }
        else
        {
            fprintf(stderr, "exceptional argument\n");
            exit(EXIT_FAILURE);
        }
    }

    float time = 0.0f;
    GridCells2D *grid_cell = new GridCells2D();
    Scene2D scene(grid_cell, time);
    Simulator2D *simulator = new Simulator2D(grid_cell, mode);

    // initialize OpenGL
    if (!glfwInit())
    {
        fprintf(stderr, "Initialization failed!\n");
        exit(EXIT_FAILURE);
    }

    // Create Window
    GLFWwindow *window = glfwCreateWindow(WIDTH, HEIGHT, TITLE, nullptr, nullptr);

    if (!window)
    {
        fprintf(stderr, "Window creation failed!");
        glfwTerminate();
        exit(EXIT_FAILURE);
    }

    glfwMakeContextCurrent(window);

    // register event callback function
    glfwSetMouseButtonCallback(window, simulator->mouseEvent);
    glfwSetCursorPosCallback(window, simulator->mouseMoveEvent);

    // initialize scene
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

    std::cout << "\n*** START SIMULATION ***\n";

    while (!glfwWindowShouldClose(window))
    {
        time += DT;

        simulator->update();
        scene.draw();

        // swap draw buffer
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    std::cout << "*** END ***\n\n";

    if (simulator)
    {
        delete simulator;
    }
    if (grid_cell)
    {
        delete grid_cell;
    }

    glfwTerminate();
    return 0;
}