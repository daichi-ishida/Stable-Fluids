#pragma once

/*
 * Linear interpolation and Easing Functions
 *  start and end and parameter t is required.
 */

namespace lerp
{
    inline float linear(float start, float end, float t);
    inline float inSine(float start, float end, float t);
    inline float outSine(float start, float end, float t);
    inline float inOutSine(float start, float end, float t);
    inline float inQuad(float start, float end, float t);
    inline float outQuad(float start, float end, float t);
    inline float inOutQuad(float start, float end, float t);
    inline float inExp(float start, float end, float t);
    inline float outExp(float start, float end, float t);
    inline float inOutExp(float start, float end, float t);
    inline float inCirc(float start, float end, float t);
    inline float outCirc(float start, float end, float t);
    inline float inOutCirc(float start, float end, float t);
    inline float inBack(float start, float end, float t, float b = 1.7);
    inline float outBack(float start, float end, float t, float b = 1.7);
    inline float inOutBack(float start, float end, float t, float b = 1.7);
    inline float inElastic(float start, float end, float t, float a = 8);
    inline float outElastic(float start, float end, float t, float a = 8);
    inline float inOutElastic(float start, float end, float t, float a = 8);
    inline float inBounce(float start, float end, float t);
    inline float outBounce(float start, float end, float t);;
}