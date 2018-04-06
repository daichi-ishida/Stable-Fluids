#pragma once

class DisplayManager
{
  public:
    DisplayManager();
    ~DisplayManager();

  private:
    const int m_width;
    const int m_height;
    const int m_fps_cap;
};