#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <malloc.h>
#include "skyfill_tif.h"
#include "colorspace_conversions.h"

#define min_f(a, b, c)  (fminf(a, fminf(b, c)))
#define max_f(a, b, c)  (fmaxf(a, fmaxf(b, c)))

void rgb2hsv16(const uint16_t src_r, const uint16_t src_g, const uint16_t src_b, float *h, float *s, float *v)
{
    float r = (float)src_r / MAX16f;
    float g = (float)src_g / MAX16f;
    float b = (float)src_b / MAX16f;

    float max = max_f(r, g, b);
    float min = min_f(r, g, b);

    *v = max;

    if (max == 0.0f) {
        *s = 0;
        *h = 0;
    }
    else if (max - min == 0.0f) {
        *s = 0;
        *h = 0;
    }
    else {
        *s = (max - min) / max;

        if (max == r) {
            *h = 60 * ((g - b) / (max - min)) + 0;
        }
        else if (max == g) {
            *h = 60 * ((b - r) / (max - min)) + 120;
        }
        else {
            *h = 60 * ((r - g) / (max - min)) + 240;
        }
    }

    if (*h < 0) *h += 360.0f;
}

void hsv2rgb16(float h, float s, float v, uint16_t *dst_r, uint16_t *dst_g, uint16_t *dst_b)
{

    float r, g, b; // 0.0-1.0

    int   hi = (int)(h / 60.0f) % 6;
    float f  = (h / 60.0f) - hi;
    float p  = v * (1.0f - s);
    float q  = v * (1.0f - s * f);
    float t  = v * (1.0f - s * (1.0f - f));

    switch(hi) {
        case 0: r = v, g = t, b = p; break;
        case 1: r = q, g = v, b = p; break;
        case 2: r = p, g = v, b = t; break;
        case 3: r = p, g = q, b = v; break;
        case 4: r = t, g = p, b = v; break;
        case 5: r = v, g = p, b = q; break;
        default: r = 0.f, g = 0.f, b = 0.f; break;
    }

    *dst_r = (uint16_t)(r * MAX16); // dst_r : 0-255
    *dst_g = (uint16_t)(g * MAX16); // dst_r : 0-255
    *dst_b = (uint16_t)(b * MAX16); // dst_r : 0-255
}
