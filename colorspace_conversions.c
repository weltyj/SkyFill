#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <malloc.h>
#include "skyfill_tif.h"
#include "colorspace_conversions.h"

#define min_f(a, b, c)  (fminf(a, fminf(b, c)))
#define max_f(a, b, c)  (fmaxf(a, fmaxf(b, c)))

#define no_V_IS_REALLY_I
#ifdef V_IS_REALLY_I
// The call is still ...hsv..., but the color model is HSI
void rgb2hsv16(const uint16_t src_r, const uint16_t src_g, const uint16_t src_b, float *H, float *S, float *I)
{
    float r = (float)src_r / MAX16f;
    float g = (float)src_g / MAX16f;
    float b = (float)src_b / MAX16f;

    float max = max_f(r, g, b);
    float min = min_f(r, g, b);

    *I = (r+b+g)/3.;

    // handle completely unsaturated, avoiding floating point precision errors
    if(src_r == src_g && src_g == src_b) {
	*S = 0. ;
	*H = 0. ;
	return ;
    }

    *S = 1.-min / *I ;

    if (max == r) {
	*H = 60 * ((g - b) / (max - min)) + 0;
    }
    else if (max == g) {
	*H = 60 * ((b - r) / (max - min)) + 120;
    }
    else {
	*H = 60 * ((r - g) / (max - min)) + 240;
    }

    if (*H < 0) *H += 360.0f;
}

void hsv2rgb16(float H, float S, float I, uint16_t *dst_r, uint16_t *dst_g, uint16_t *dst_b)
{

    float r, g, b; // 0.0-1.0

    float Htag = H / 60. ;
    float Z = 1. - fabsf(fmod(Htag,2.) - 1.) ;
    float C = (3.*I*S)/(1.+Z) ;
    float X = C * Z ;

    if(Htag < 1.) {
	r = C  ; g = X  ; b = 0. ;
    } else if(Htag < 2.) {
	r = X  ; g = C  ; b = 0. ;
    } else if(Htag < 3.) {
	r = 0. ; g = C  ; b = X ;
    } else if(Htag < 4.) {
	r = 0. ; g = X  ; b = C ;
    } else if(Htag < 5.) {
	r = X  ; g = 0. ; b = C ;
    } else {
	r = C  ; g = 0. ; b = X ;
    }

    float m = I*(1.-S) ;
    r += m ;
    g += m ;
    b += m ;

/*      if(r > 1.0) r=1.0 ;  */
/*      if(g > 1.0) g=1.0 ;  */
/*      if(b > 1.0) b=1.0 ;  */

    *dst_r = (uint16_t)((r) * MAX16); // dst_r : 0-255
    *dst_g = (uint16_t)((g) * MAX16); // dst_g : 0-255
    *dst_b = (uint16_t)((b) * MAX16); // dst_b : 0-255
}
#else

#define R 0
#define G 1
#define B 2
#define H 0
#define S 1
#define V 2
void rgb2hsv16(const uint16_t *rgb, float *hsv)
{
    float r = (float)rgb[0] / MAX16f;
    float g = (float)rgb[1] / MAX16f;
    float b = (float)rgb[2] / MAX16f;

    float max = max_f(r, g, b);
    float min = min_f(r, g, b);

    hsv[V] = max;

    if (max == 0.0f) {
        hsv[H] = 0;
        hsv[S] = 0;
    }
    else if (max - min == 0.0f) {
        hsv[H] = 0;
        hsv[S] = 0;
    }
    else {
        hsv[S] = (max - min) / max;

        if (max == r) {
            hsv[H] = 60 * ((g - b) / (max - min)) + 0;
        }
        else if (max == g) {
            hsv[H] = 60 * ((b - r) / (max - min)) + 120;
        }
        else {
            hsv[H] = 60 * ((r - g) / (max - min)) + 240;
        }
    }

    if (hsv[H] < 0) hsv[H] += 360.0f;
}

void hsv2rgb16(float *hsv, uint16_t *rgb)
{

    float r, g, b; // 0.0-1.0

    int   hi = (int)(hsv[H] / 60.0f) % 6;
    float f  = (hsv[H] / 60.0f) - hi;
    float p  = hsv[V] * (1.0f - hsv[S]);
    float q  = hsv[V] * (1.0f - hsv[S] * f);
    float t  = hsv[V] * (1.0f - hsv[S] * (1.0f - f));

    switch(hi) {
        case 0: r = hsv[V], g = t, b = p; break;
        case 1: r = q, g = hsv[V], b = p; break;
        case 2: r = p, g = hsv[V], b = t; break;
        case 3: r = p, g = q, b = hsv[V]; break;
        case 4: r = t, g = p, b = hsv[V]; break;
        case 5: r = hsv[V], g = p, b = q; break;
        default: r = 0.f, g = 0.f, b = 0.f; break;
    }

    rgb[R] = (uint16_t)(r * MAX16); // dst_r : 0-255
    rgb[G] = (uint16_t)(g * MAX16); // dst_r : 0-255
    rgb[B] = (uint16_t)(b * MAX16); // dst_r : 0-255
}
#endif
