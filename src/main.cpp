#define GLFW_INCLUDE_GLU
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "constants.hpp"
#include "GridCells2D.hpp"
#include "Scene2D.hpp"
#include "Simulator2D.hpp"

//void resizeGL(GLFWwindow *window, int width, int height);

int main()
{
    float time = 0.0f;
    GridCells2D *grid_cell = new GridCells2D();
    Scene2D scene(grid_cell, time);
    Simulator2D *simulator = new Simulator2D(grid_cell);

    // initialize OpenGL
    if (!glfwInit())
    {
        fprintf(stderr, "Initialization failed!\n");
        exit(EXIT_FAILURE);
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

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
    //glfwSetScrollCallback(window, wheelEvent);

    //initialize GLEW
    glewExperimental = true;
    if (glewInit() != GLEW_OK)
    {
        fprintf(stderr, "GLEW initialization failed!\n");
        exit(EXIT_FAILURE);
    }

    // register resize callback function
    // glfwSetWindowSizeCallback(window, resizeGL);

    // initialize scene
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

    std::cout << "\n*** START PARTICLE-SIMULATION ***\n";

    scene.writeData();

    while (!glfwWindowShouldClose(window))
    {
        time += DT;

        simulator->update();
        scene.draw();
        scene.writeData();

        // swap draw buffer
        glfwSwapBuffers(window);
        glfwPollEvents();
        if (time >= FINISH_TIME)
        {
            break;
        }
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

// void resizeGL(GLFWwindow *window, int width, int height)
// {

//     // GLFW管理のウィンドウサイズを変更
//     glfwSetWindowSize(window, WIDTH, HEIGHT);

//     // 実際のウィンドウサイズ (ピクセル数) を取得
//     int renderBufferWidth, renderBufferHeight;
//     glfwGetFramebufferSize(window, &renderBufferWidth, &renderBufferHeight);

//     // ビューポート変換の更新
//     glViewport(0, 0, renderBufferWidth, renderBufferHeight);
//     glLoadIdentity();
//     glOrtho(-0.02, LENGTH + 0.02, -0.02, LENGTH + 0.02, -1.0, 1.0);
// }