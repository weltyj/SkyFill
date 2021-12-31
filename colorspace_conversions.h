#ifndef COLORSPACE_CONVERSIONS_H
#define COLORSPACE_CONVERSIONS_H
#define min_f(a, b, c)  (fminf(a, fminf(b, c)))
#define max_f(a, b, c)  (fmaxf(a, fmaxf(b, c)))

void rgb2hsv16(const uint16_t src_r, const uint16_t src_g, const uint16_t src_b, float *h, float *s, float *v) ;
void hsv2rgb16(float h, float s, float v, uint16_t *dst_r, uint16_t *dst_g, uint16_t *dst_b) ;
#endif
