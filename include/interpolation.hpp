#pragma once

/*
 * Linear interpolation and Easing Functions
 *  start and end and parameter t is required.
 */

namespace lerp
{
    float linear(float start, float end, float t);
    float inSine(float start, float end, float t);
    float outSine(float start, float end, float t);
    float inOutSine(float start, float end, float t);
    float inQuad(float start, float end, float t);
    float outQuad(float start, float end, float t);
    float inOutQuad(float start, float end, float t);
    float inExp(float start, float end, float t);
    float outExp(float start, float end, float t);
    float inOutExp(float start, float end, float t);
    float inCirc(float start, float end, float t);
    float outCirc(float start, float end, float t);
    float inOutCirc(float start, float end, float t);
    float inBack(float start, float end, float t, float b = 1.7);
    float outBack(float start, float end, float t, float b = 1.7);
    float inOutBack(float start, float end, float t, float b = 1.7);
    float inElastic(float start, float end, float t, float a = 8);
    float outElastic(float start, float end, float t, float a = 8);
    float inOutElastic(float start, float end, float t, float a = 8);
    float inBounce(float start, float end, float t);
    float outBounce(float start, float end, float t);
;
}