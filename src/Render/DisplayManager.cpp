#include "Render/DisplayManager.hpp"
#include "constants.hpp"

DisplayManager::DisplayManager() : m_width(WIDTH),
                                   m_height(HEIGHT),
                                   m_fps_cap(FPS_CAP)
{
}

DisplayManager::~DisplayManager()
{
}
