#pragma once

class Loader
{
  public:
    Loader();
    ~Loader();

    void loadToVAO(float positions[], int indices[]);
    void cleanUp();

  private:
    int createVAO();
    void storeDataInAttributeList(int attribute_number, int coordinate_size, float data[]);
    void unbindVAO();
    void bindIndicesBuffer(int indices[]);

    const int m_vaos;
    const int m_vbos;
};