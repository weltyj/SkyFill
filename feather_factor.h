#ifndef FEATHER_FACTOR_H
#define FEATHER_FACTOR_H

// how much of the new estimate to use at pixel x, y
float raw_compute_feather_at_y(int x, int y, int feather_length, uint16_t b, float feather_factor, SKYFILL_DATA_t *pData) ;

int raw_compute_feather_length(int x, int *pFeather_end_y, float depth_of_fill, float feather_factor, int extra, SKYFILL_DATA_t *pData) ;

// how much of the new estimate to use at pixel x, y
float compute_feather_at_y(int x, int y, int feather_length, uint16_t b, float feather_factor, SKYFILL_DATA_t *pData) ;

int compute_feather_length(int x, int *pFeather_end_y, float depth_of_fill, float feather_factor, int extra, SKYFILL_DATA_t *pData) ;

// compute feather length using given end of sky
int compute_feather_length_with_eos(int x, int *pFeather_end_y, float depth_of_fill, float feather_factor, int extra, SKYFILL_DATA_t *pData, int end_of_sky) ;
#endif
