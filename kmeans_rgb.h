#ifndef KMEANS_RGB_H
#define KMEANS_RGB_H

#include "estimate_sky.h"

void kmeans_analysis(tdata_t *image, int n_pixels, int *x_list, int *y_list, struct sky_pixel_data *sp_list, int *classification, int first_cluster_pt_index) ;

#endif
