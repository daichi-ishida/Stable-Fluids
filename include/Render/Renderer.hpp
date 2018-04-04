#pragma once
#include <Eigen/Core>

class Renderer
{
  public:
    Renderer();
    ~Renderer();

    void prepare();
    void render();
    void creatProjectionMatrix();

  private:
    const float m_fov;
    const float m_near_plane;
    const float m_far_plane;

    Eigen::Matrix4f m_projection_matrix;
};