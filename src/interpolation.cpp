#define _USE_MATH_DEFINES
#include <cmath>
#include "interpolation.hpp"

/*
 * Linear interpolation and Easing Functions
 *  start and end and parameter t is required.
 */

namespace lerp
{
    float linear(float start, float end, float t)
    {
        return t * end + (1 - t) * start;
    }

    float inSine(float start, float end, float t)
    {
        float th = M_PI * t / 2;
        float s = std::cos(th);
        return s * start + (1 - s) * end;
    }

    float outSine(float start, float end, float t)
    {
        float th = M_PI * t / 2;
        float s = std::sin(th);
        return s * end + (1 - s) * start;
    }

    float inOutSine(float start, float end, float t)
    {
        float th = M_PI * (t - 0.5);
        float s = 0.5 * (std::sin(th) + 1);
        return s * end + (1 - s) * start;
    }

    float inQuad(float start, float end, float t)
    {
        float s = t * t;
        return s * end + (1 - s) * start;
    }

    float outQuad(float start, float end, float t)
    {
        float u = t - 1;
        float s = u * u;
        return (1 - s) * end + s * start;
    }

    float inOutQuad(float start, float end, float t)
    {
        if (t <= 0.5)
        {
            return 2 * t * t * end + (1 - 2 * t * t) * start;
        }
        else
        {
            float s = t - 1;
            return (1 - 2 * s * s) * end + 2 * s * s * start;
        }
    }

    float inExp(float start, float end, float t)
    {
        float s = std::exp(10 * (t - 1));
        return s * end + (1 - s) * start;
    }

    float outExp(float start, float end, float t)
    {
        float s = std::exp(-10 * t);
        return (1 - s) * end + s * start;
    }

    float inOutExp(float start, float end, float t)
    {
        if (t <= 0.5)
        {
            float s = 0.5 * std::exp(10 * (t - 0.5));
            return s * end + (1 - s) * start;
        }
        else
        {
            float s = -0.5 * std::exp(10 * (0.5 - t)) + 1;
            return s * end + (1 - s) * start;
        }
    }

    float inCirc(float start, float end, float t)
    {
        float s = std::sqrt(1 - t * t);
        return (1 - s) * end + s * start;
    }

    float outCirc(float start, float end, float t)
    {
        float u = t - 1;
        float s = std::sqrt(1 - u * u);
        return s * end + (1 - s) * start;
    }

    float inOutCirc(float start, float end, float t)
    {
        if (t <= 0.5)
        {
            float s = std::sqrt(0.25 - t * t);
            return (0.5 + s) * start + (0.5 - s) * end;
        }
        else
        {
            float u = t - 1;
            float s = std::sqrt(0.25 - u * u);
            return (0.5 - s) * start + (0.5 + s) * end;
        }
    }

    float inBack(float start, float end, float t, float b)
    {
        float s = t * t * ((b + 1) * t - b);
        return s * end + (1 - s) * start;
    }

    float outBack(float start, float end, float t, float b)
    {
        float u = 1 - t;
        float s = -u * u * ((b + 1) * u - b) + 1;
        return s * end + (1 - s) * start;
    }

    float inOutBack(float start, float end, float t, float b)
    {
        if (t <= 0.5)
        {
            float s = 2 * t * t * (2 * (b + 1) * t - b);
            return s * end + (1 - s) * start;
        }
        else
        {
            float u = 1 - t;
            float s = -2 * u * u * (2 * (b + 1) * u - b) + 1;
            return s * end + (1 - s) * start;
        }
    }

    float inElastic(float start, float end, float t, float a)
    {
        float omg = 6.5 * M_PI;
        float s = std::exp(-a * (1 - t)) * std::cos(omg * (1 - t));
        return s * end + (1 - s) * start;
    }

    float outElastic(float start, float end, float t, float a)
    {
        float omg = 6.5 * M_PI;
        float s = -std::exp(-a * t) * std::cos(omg * t) + 1;
        return s * end + (1 - s) * start;
    }

    float inOutElastic(float start, float end, float t, float a)
    {
        float omg = 7 * M_PI;
        float u = t - 0.5;
        if (t <= 0.5)
        {
            float s = 0.5 * std::exp(-a * (-u)) * std::cos(omg * (-u));
            return s * end + (1 - s) * start;
        }
        else
        {
            float s = -0.5 * std::exp(-a * u) * std::cos(omg * u) + 1;
            return s * end + (1 - s) * start;
        }
    }

    float inBounce(float start, float end, float t)
    {
        float omg = 3.5 * M_PI;
        float k = std::floor((omg * (1 - t) + M_PI / 2) / M_PI);
        float u = std::exp(-k);
        float s = std::fabs(u * std::cos(omg * (1 - t)));
        return s * end + (1 - s) * start;
    }

    float outBounce(float start, float end, float t)
    {
        float omg = 3.5 * M_PI;
        float k = std::floor((omg * t + M_PI / 2) / M_PI);
        float u = std::exp(-k);
        float s = 1 - std::fabs(u * std::cos(omg * t));
        return s * end + (1 - s) * start;
    }
}
