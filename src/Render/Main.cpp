#include <iostream>

class Triangle {
public:
  typedef std::shared_ptr<triangle> Ptr;
  Triangle(GLfloat size) {
    GLfloat low = -size/2.0f;
    GLfloat high = size/2.0f;
    GLfloat vertices[3*3] = {
      low, low, low,
      low, high, low,
      high, high, low,
    };
    GLuint indecise[ELEMENT_SIZE] = {
      0, 1, 2
    };
    //VAO（バッファデータ・設定群)の生成
    glGenVertexArrays(1, &vertex_array_object_);
    //現在のVAOに設定
    glBindVertexArray(vertex_array_object_);
    //VBO（バッファ）の生成
    glGenBuffers(1, &vertex_buffer_object_);
    //現在のバッファに設定
    glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer_object_);
    //現在のバッファにデータを送信
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
    //データ型の設定(float * 3)
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, nullptr);
    //シェーダで0番地を使用
    glEnableVertexAttribArray(0);
    //インデックスの設定
    glGenBuffers(1, &index_buffer_object_);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer_object_);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indecise), indecise, GL_DYNAMIC_DRAW);
    //VAO設定を終了する
    glBindVertexArray(0);
  }
 
  ~Triangle() {
    glDeleteBuffers(1, &index_buffer_object_);
    glDeleteBuffers(1, &vertex_buffer_object_);
    glDeleteVertexArrays(1, &vertex_array_object_);
  }
   
  void Draw() {
    //VAOの有効化（VAOに割り当てた設定とバッファが復元される)
    glBindVertexArray(vertex_array_object_);
    //インデックスを用いて描画する
    glDrawElements(GL_TRIANGLES, ELEMENT_SIZE, GL_UNSIGNED_INT, nullptr);
    //VAOの無効化
    glBindVertexArray(0);
  }
private:
  enum { ELEMENT_SIZE = 3};
  GLuint vertex_array_object_;
  GLuint vertex_buffer_object_;
  GLuint index_buffer_object_;
};

