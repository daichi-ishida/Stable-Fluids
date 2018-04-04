#include "Render/Renderer.hpp"
#include "constants.hpp"

Renderer::Renderer() : m_fov(FOV),
                       m_near_plane(NEAR_PLANE),
                       m_far_plane(FAR_PLANE)
{
    creatProjectionMatrix();
}

Renderer::~Renderer()
{
}
