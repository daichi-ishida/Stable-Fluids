#include <cmath>
#include "interpolation.hpp"

/*
 * Linear Interpolation and Easing Functions
 *  start and end and parameter t is required.
 */

namespace lerp
{

float linear(float start, float end, float t)
{
}

return float t * float end + (1 - float t) * float start

                                 float
                                 InSine(float start, float end, float t) :
    ''' 開始側が緩やかなSineのイージング * 第一引数 : 始点 * 第二引数 : 終点 * 第三引数 : 0 -
           1の時刻パラメータ
    ''' th = math.pi * float t / 2 s = math.cos(th) return s * float start + (1 - s) * float end

                                                                                  float
                                                                                  OutSine(float start, float end, float t) :
    ''' 終了側が緩やかなSineのイージング * 第一引数 : 始点 * 第二引数 : 終点 * 第三引数 : 0 -
                                        1の時刻パラメータ
    ''' th = math.pi * float t / 2 s = math.sin(th) return s * float end + (1 - s) * float start

                                                                                float
                                                                                InOutSine(float start, float end, float t) :
    ''' 開始・終了が緩やかなSineのイージング * 第一引数 : 始点 * 第二引数 : 終点 * 第三引数 : 0 -
                                        1の時刻パラメータ
    ''' th = math.pi * (float t - 0.5)
                            s = 0.5 * (math.sin(th) + 1) return s * float end + (1 - s) * float start

                                                                                    float
                                                                                    InQuad(float start, float end, float t) :
    ''' 開始側が緩やかな二次関数イージング * 第一引数 : 始点 * 第二引数 : 終点 * 第三引数 : 0 -
                                1の時刻パラメータ
    ''' s = float t * float t return s * float end + (1 - s) * float start

                                                          float
                                                          OutQuad(float start, float end, float t) :
    ''' 終了側が緩やかな二次関数のイージング * 第一引数 : 始点 * 第二引数 : 終点 * 第三引数 : 0 -
             1の時刻パラメータ
    ''' u = float t - 1 s = u * u return (1 - s) * float end + s * float start

                                                                    float
                                                                    InOutQuad(float start, float end, float t) :
    ''' 開始・終了が緩やかな二次関数のイージング * 第一引数 : 始点 * 第二引数 : 終点 * 第三引数 : 0 -
                                 1の時刻パラメータ
    ''' if float t <= 0.5 : return 2 * float t * float t * float end + (1 - 2 * float t *float t) * float start else : s = float t - 1 return (1 - 2 * s * s) * float end + 2 * s * s * float start

                                                                                                                                                                                 float
                                                                                                                                                                                 InExpo(float start, float end, float t) :
    ''' 開始側が緩やかな指数関数のイージング * 第一引数 : 始点 * 第二引数 : 終点 * 第三引数 : 0 - 1の時刻パラメータ
    ''' s = math.exp(10 * (float t - 1)) return s * float end + (1 - s) * float start

                                                                     float
                                                                     OutExpo(float start, float end, float t) :
    ''' 終了側が緩やかな指数関数のイージング * 第一引数 : 始点 * 第二引数 : 終点 * 第三引数 : 0 -
             1の時刻パラメータ
    ''' s = math.exp(-10 * float t) return (1 - s) * float end + s * float start

                                                                      float
                                                                      InOutExpo(float start, float end, float t) :
    ''' 開始・終了が緩やかな指数関数のイージング * 第一引数 : 始点 * 第二引数 : 終点 * 第三引数 : 0 -
                 1の時刻パラメータ
    ''' if float t <=
             0.5 : s = 0.5 * math.exp(10 * (float t - 0.5)) return s * float end + (1 - s) * float start else : s = -0.5 * math.exp(10 * (0.5 - float t)) + 1 return s * float end + (1 - s) * float start

                                                                                                                                                                                         float
                                                                                                                                                                                         InCirc(float start, float end, float t) :
    ''' 開始側が緩やかな円のイージング * 第一引数 : 始点 * 第二引数 : 終点 * 第三引数 : 0 - 1の時刻パラメータ
    ''' s = math.sqrt(1 - float t *float t) return (1 - s) * float end + s * float start

                                                                              float
                                                                              OutCirc(float start, float end, float t) :
    ''' 終了側が緩やかな円のイージング * 第一引数 : 始点 * 第二引数 : 終点 * 第三引数 : 0 -
             1の時刻パラメータ
    ''' u = float t - 1 s = math.sqrt(1 - u * u) return s * float end + (1 - s) * float start

                                                                             float
                                                                             InOutCirc(float start, float end, float t) :
    ''' 開始・終了が緩やかな円のイージング * 第一引数 : 始点 * 第二引数 : 終点 * 第三引数 : 0 -
                                 1の時刻パラメータ
    ''' if float t <= 0.5 : s = math.sqrt(0.25 - float t *float t) return (0.5 + s) * float start + (0.5 - s) * float end else : u = float t - 1 s = math.sqrt(0.25 - u * u) return (0.5 - s) * float start + (0.5 + s) * float end

                                                                                                                                                                                                                   float
                                                                                                                                                                                                                   InBack(float start, float end, float t, b = 1.7) :
    ''' 開始側が緩やかな三次関数のイージング * 第一引数 : 始点 * 第二引数 : 終点 * 第三引数 : 0 - 1の時刻パラメータ 第四引数 : 戻り具合（極値に関係）,
         デフォルト : 1.7
    ''' s = float t * float t * ((b + 1) * float t - b) return s * float end + (1 - s) * float start

                                                                                    float
                                                                                    OutBack(float start, float end, float t, b = 1.7) :
    ''' 終了側が緩やかな三次関数のイージング * 第一引数 : 始点 * 第二引数 : 終点 * 第三引数 : 0 -
             1の時刻パラメータ 第四引数 : 追い越し具合（極値に関係）,
         デフォルト : 1.7
    ''' u = 1 - float t
                     s = -u * u * ((b + 1) * u - b) + 1 return s * float end + (1 - s) * float start

                                                                                   float
                                                                                   InOutBack(float start, float end, float t, b = 1.7) :
    ''' 開始・終了が緩やかな三次関数のイージング * 第一引数 : 始点 * 第二引数 : 終点 * 第三引数 : 0 -
                         1の時刻パラメータ 第四引数 : 戻り・追い越し具合（極値に関係）,
         デフォルト : 1.7
    ''' if float t <= 0.5 : s = 2 * float t * float t * (2 * (b + 1) * float t - b) return s * float end + (1 - s) * float start else : u = 1 - float t
                                                                                                                                                     s = -2 * u * u * (2 * (b + 1) * u - b) + 1 return s * float end + (1 - s) * float start

                                                                                                                                                                                                                           float
                                                                                                                                                                                                                           InElastic(float start, float end, float t, a = 8) :
    ''' 開始側が緩やかな減衰曲線のイージング * 第一引数 : 始点 * 第二引数 : 終点 * 第三引数 : 0 - 1の時刻パラメータ 第四引数 : 減衰具合,
         デフォルト : 8
    ''' omg = 6.5 * math.pi
                         s = math.exp(-a * (1 - float t)) * math.cos(omg * (1 - float t)) return s * float end + (1 - s) * float start

                                                                                                                     float
                                                                                                                     OutElastic(float start, float end, float t, a = 8) :
    ''' 終了側が緩やかな減衰曲線のイージング * 第一引数 : 始点 * 第二引数 : 終点 * 第三引数 : 0 -
                             1の時刻パラメータ 第四引数 : 減衰具合,
         デフォルト : 8
    ''' omg = 6.5 * math.pi
                         s = -math.exp(-a * float t) * math.cos(omg * float t) + 1 return s * float end + (1 - s) * float start

                                                                                                              float
                                                                                                              InOutElastic(float start, float end, float t, a = 8) :
    ''' 開始・終了が緩やかな減衰曲線のイージング * 第一引数 : 始点 * 第二引数 : 終点 * 第三引数 : 0 -
                             1の時刻パラメータ 第四引数 : 減衰具合,
         デフォルト : 8
    ''' omg = 7 * math.pi
                       u = float t - 0.5 if float t <= 0.5 : s = 0.5 * math.exp(-a * (-u)) * math.cos(omg * (-u)) return s * float end + (1 - s) * float start else : s = -0.5 * math.exp(-a * u) * math.cos(omg * u) + 1 return s * float end + (1 - s) * float start

                                                                                                                                                                                                                                                     float
                                                                                                                                                                                                                                                     InBounce(float start, float end, float t) :
    ''' 開始側が緩やかなバウンドのイージング * 第一引数 : 始点 * 第二引数 : 終点 * 第三引数 : 0 - 1の時刻パラメータ
    ''' omg = 3.5 * math.pi
                         k = math.floor((omg * (1 - float t) + math.pi / 2) / math.pi)
                                 u = math.exp(-k)
                                         s = math.fabs(u * math.cos(omg * (1 - float t))) return s * float end + (1 - s) * float start

                                                                                                                     float
                                                                                                                     OutBounce(float start, float end, float t) :
    ''' 終了側が緩やかなバウンドのイージング * 第一引数 : 始点 * 第二引数 : 終点 * 第三引数 : 0 -
                                             1の時刻パラメータ
    ''' omg = 3.5 * math.pi
                         k = math.floor((omg * float t + math.pi / 2) / math.pi)
                                 u = math.exp(-k)
                                         s = 1 - math.fabs(u * math.cos(omg * float t)) return s * float end + (1 - s) * float start
}
