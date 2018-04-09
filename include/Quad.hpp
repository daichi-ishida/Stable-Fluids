#pragma once
#include <memory>
#define GLFW_INCLUDE_GLU

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>

class Quad
{
public:
  typedef std::shared_ptr<Quad> Ptr;
  Quad(GLfloat size);
  ~Quad();

  void Draw();

private:
  const int ELEMENT_SIZE = 6;
  GLuint vaoID;
  GLuint vboID;
  GLuint index_bufferID;
  GLuint color_bufferID;
};