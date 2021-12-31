#ifndef PIXEL_TESTS_H
#define PIXEL_TESTS_H

int xy_is_opaque_pixel(tdata_t *image, int32_t x, int32_t y) ;
int xy_has_nonblack_pixel(tdata_t *image, int32_t x, int32_t y) ;
int wrong_hue_or_black(tdata_t *image, int32_t x, int32_t y, float hue_sky, float sat_sky, float val_sky) ;
int xy_is_transparent(tdata_t *image, int32_t x, int32_t y) ;
int xy_is_black_pixel(tdata_t *image, int32_t x, int32_t y) ;
#endif
