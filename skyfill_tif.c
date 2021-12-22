#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include "tiffio.h"
#include <math.h>
#include "amoeba_06.h"
#include "string.h"
#include "mstat.h"

/* coordinate system:
	TIF files are read.  data is stored in image[] array.  Top left is (0,0), bottom right is (width-1,height-1) ;
*/

// four macros to get or set rgb and alpha 
#define tif_get3c(image,x,y,r,g,b) {\
		r = ((uint16 *)(image[(y)]))[IMAGE_NSAMPLES*(x)+0] ; \
		g = ((uint16 *)(image[(y)]))[IMAGE_NSAMPLES*(x)+1] ; \
		b = ((uint16 *)(image[(y)]))[IMAGE_NSAMPLES*(x)+2] ; \
		}

#define tif_get4c(image,x,y,r,g,b,a) {\
		r = ((uint16 *)(image[(y)]))[IMAGE_NSAMPLES*(x)+0] ; \
		g = ((uint16 *)(image[(y)]))[IMAGE_NSAMPLES*(x)+1] ; \
		b = ((uint16 *)(image[(y)]))[IMAGE_NSAMPLES*(x)+2] ; \
		a = ((uint16 *)(image[(y)]))[IMAGE_NSAMPLES*(x)+3] ; \
		}

#define tif_set3c(image,x,y,r,g,b) {\
		((uint16 *)(image[(y)]))[IMAGE_NSAMPLES*(x)+0] = (uint16)r; \
		((uint16 *)(image[(y)]))[IMAGE_NSAMPLES*(x)+1] = (uint16)g; \
		((uint16 *)(image[(y)]))[IMAGE_NSAMPLES*(x)+2] = (uint16)b; \
		}

#define tif_set4c(image,x,y,r,g,b,a) {\
		((uint16 *)(image[(y)]))[IMAGE_NSAMPLES*(x)+0] = (uint16)r; \
		((uint16 *)(image[(y)]))[IMAGE_NSAMPLES*(x)+1] = (uint16)g; \
		((uint16 *)(image[(y)]))[IMAGE_NSAMPLES*(x)+2] = (uint16)b; \
		((uint16 *)(image[(y)]))[IMAGE_NSAMPLES*(x)+3] = (uint16)a; \
		}

#define MAX16 65535
#define HALF16 32767
#define MAX16f 65535.f
#define FMAX 1.e30f
// value which is black
#define BLK16 8621
float floatBLK=.04 ; // value at which we assume hugin has masked out a portion of the image and placed black there instead of a alpha of 0 ;
int verbose=0 ;
int fix_SOS_edges=0 ;  // by default, don't fill missing data on edges around the start of the sky
int fill_top_of_sky=0 ;  // by default, don't make start of the sky values equal in first stage of repairing the top of sky
float final_saturation_factor=1.0 ;

// global arrays, will indicate what was detected for each column (x) for these values
int16 *start_of_sky, *end_of_sky, *final_end_of_sky, *raw_start_of_sky ;  // raw start of sky is as detected before any repairs on image

#define MIN(a,b) ( (a) < (b) ? (a) : (b) )
#define MAX(a,b) ( (a) > (b) ? (a) : (b) )
void predict_sky_hsv(float px, float py, float *pH, float *pS, float *pV) ;
void predict_sky_huesat_from_val(float vhat, float *pH, float *pS, float px, float py) ;


// set to 1 to grid search for starting sun and horizon position
int estimate_only=0 ;
int show_raw_prediction=0 ;
int CIE_sky_index=-1 ; // if 3, will only use CIE coefs to predict given CIE index grid search

int fix_sky_hue =0 ;  // set all found sky hue pixels to the sky_hue_want ;
int min_sky_hue_mask = 0 ;
int max_sky_hue_mask = -1 ;
float sky_hue_want = 212. ;
float sky_hue_v_factor = 1./160. ;
int depth_of_fill_is_absolute = 0 ;
int depth_of_fill_absolute_y ;
float horizon_curvature=0. ; // a proportion, up or down, at the edge of the image relative to the center of the image


// defaults for detecting a change from sky pixels to end of sky pixels
float sky_hue_tolerance=360. ;
float sky_sat_tolerance=.3 ;
float sky_val_tolerance=.015 ;
float min_sky_end_p=1.0 ; // the minimum allowed for end of sky, 0 is bottom of image, 1 is top
float lowest_sky_py_found=1.0 ;  // after sky is detected, this is set
float repair_tos_thresh=0.02 ;  // how far off a r,g, or b pixel is off from mean before it is not used in local
                                // sky color calculation -- the get_mean_rgb() function.

float sky_gamma=-1. ; // applied to the luminosity estimate from the CIE sky model, reduces intensity near sun, negative value disables
float final_sky_gamma=1.0 ; // applied to luminosity estimate from the CIE sky model, but *AFTER* fitting has been, allows manual control of result
int model_is_being_fit=1 ;

float exposure_factor=1. ; // exposure factor

int have_pto_fov=0 ;
int n_strips ;
int32 IMAGE_HEIGHT ;
int32 IMAGE_WIDTH ;
uint16 IMAGE_NSAMPLES ; // must be 16 bit, needed by call to tif library
int IMAGE_HAS_ALPHA=0 ;

// horizon_py is at x=0 ;
float sun_x = 120., sun_py = .99, horizon_py = .5 ;
float hue_horizon = 200, sat_depth = .8, sun_lum=.9 ;

// as determined from the sample points
float hue_sky = 212 ;
float val_sky=1. ;
float sat_sky=.7 ;

float horizon_lift_angle_radians=0. ;

// constraints on sun x position
float sun_x_lo = -180. ;
float sun_x_hi = 180. ;

// constraints on horizon_py position
float horizon_py_lo = -0.5 ;
float horizon_py_hi = 1.5 ;

float horizon_slope=.05 ;
float FOV_horizontal=120. ;
float proportion_to_radian_factor_x ; // how many radians per proportion of full image in X dimension ?
float proportion_to_radian_factor_y ; // how many radians per proportion of full image in Y dimension ?
float maximum_vhat = 1. ; // maximum vhat in sky dome

#ifdef TRY_SUN_STRETCH_DIST
float sun_dist_stretch_factor=0.5 ; // adjusts px as a function of distance to sun_x for CIE calculation
#endif

int max_end_of_sky = -1 ;
int min_end_of_sky = -1 ;
int max_start_of_sky = -1 ;
int min_start_of_sky = -1 ;

// 2021, try nonlinear feathering
float nonlinear_feather=1.0 ;

//float perez_A=-1. ; // horizon-zenith gradient, -5 to 5
//float perez_B=-0.32 ; // gradient intensity, -10 to 0
//float perez_C=10. ; // circumsolar intensity, 0 to 25
//float perez_D=-3. ; // circumsolar radius, -10 to 0
//float perez_E=0.45 ; // backscattering effect, -1 to 5

//float perez_A=0.001 ; // horizon-zenith gradient, -5 to 5
//float perez_B=-1.00 ; // gradient intensity, -10 to 0
//float perez_C=2. ; // circumsolar intensity, 0 to 25
//float perez_D=-1.5 ; // circumsolar radius, -10 to 0
//float perez_E=0.15 ; // backscattering effect, -1 to 5

// for clear blue sky with some haze this will produce
// a result with high luminosity around the sun corona
// as well as some lightening at the horizon

float perez_A=-1.5 ; // horizon-zenith gradient, -5 to 5
float perez_B=-0.80 ; // gradient intensity, -10 to 0
float perez_C=2. ; // circumsolar intensity, 0 to 25
float perez_D=-2.4 ; // circumsolar radius, -10 to 0
float perez_E=-0.15 ; // backscattering effect, -1 to 5

#ifdef TRY_SUN_STRETCH_DIST
#define MAX_PARMS 17
#else
#define MAX_PARMS 16
#endif
struct OPT_PARM {
    float *variable_address ;
    char abreviation[8] ; // parameter name abreviation
    char name[80] ; // parameter name
    float default_value ; // parameter value
    float lo ; // parameter minimum allowed value
    float hi ; // parameter maximum allowed value
    int optimize_flag ;  // if 1, include in optimization, if 0 do not include ;
    int grid_optimize_flag ;  // if 1, include in grid search optimization, if 0 do not include ;
    char *units ;
} opt_parms[MAX_PARMS] = {
    {&sun_x, "sx", "sun_x", 0., -180., 180., 1, 1, "degrees: 0 is middle of image"}, // 0 is straight ahead, -180 or 180 is straight behind
    {&sun_py, "sy", "sun_py", .99, 0.0, 3., 1, 1, "proportion: 0 is bottom, 1 is top"}, // 0 is bottom of image, 1 is top of image
    {&horizon_py, "hy", "horizon_py", .5, -.2, 1.5, 1, 1, "proportion: 0 is bottom, 1 is top"}, // 0 is bottom of image, 1 is top of image
    {&horizon_curvature, "hc", "horizon curvature", 0.0, -.2,   .2, 1, 1, "proportion"},
    {&hue_sky, "hs", "hue sky", 212., 180., 250., 1, 0, "Hue from HSV model"},
    {&hue_horizon, "hh", "hue_horizon", 212., 180., 250., 1, 0, "Horizon hue, not currenly used"},
    {&sat_depth, "sd", "sat_depth", .8, 0.01, 1.1, 1, 1, "Amount of saturation to apply, 0 to 1"},
    {&sun_lum, "sl", "sky_lum", 1., 0.01, 20.1, 1, 1, "Factor applied to predicted sky intensity"},
    {&FOV_horizontal, "fov", "FOV_horizontal", 90., 1., 360., 1, 0, "The horizontal field of view of the image, degrees"},
    {&perez_A, "A", "perez_A", .001, -5., 5., 1, 0, "The A parameter of the perez CIE sky model"},
    {&perez_B, "B", "perez_B", -1., -10., -0.001, 1, 0, "The B parameter of the perez CIE sky model"},
    {&perez_C, "C", "perez_C", 2., 0., 25., 1, 0, "The C parameter of the perez CIE sky model"},
    {&perez_D, "D", "perez_D", -1.5, -10., -0.001, 1, 0, "The D parameter of the perez CIE sky model"},
    {&perez_E, "E", "perez_E", 0.15, -1., 5., 1, 0, "The E parameter of the perez CIE sky model"},
    {&sky_gamma, "sg", "sky_gamma", 1.0, 1.e-4, 1000, 1, 1, "Gamma value for sky intensity"},
    {&exposure_factor, "ef", "exposure factor",    1.0,  0.1, 10.0, 1, 1, "Final exposure appled to sky intensity"},
#ifdef TRY_SUN_STRETCH_DIST
    {&sun_dist_stretch_factor, "sun_dist_stretch_factor", 0.5, -2., 2., 1},
#endif
    } ;

// restrict CIE sky types used in optimization
#define ALL_CIE_SKIES 0x01
#define UNIFORM_CIE_SKIES 0x02
#define NONUNIFORM_CIE_SKIES 0x03
int allowed_sky_type = ALL_CIE_SKIES ;

// these are the "standard" CIE parameters for various types of skies.
#define N_CIE_STD_PARMS 16
struct CIE_STD_PARMS {
    int type ;
    int gradation ;
    float A ;
    float B ;
    float C ;
    float D ;
    float E ;
    char description[120] ;
    }
    cie_parms[N_CIE_STD_PARMS] = {
	{ 1,1,   4.0, -0.70,  0, -1.0, 0.00, "CIE standard overcast sky, alternative form.  Steep lum gradation towars zenith, azimuthal uniformity"},
	{ 2,1,   4.0, -0.70,  2, -1.5, 0.15, "Overcast, with steep lum gradation and slight brightening towards sun"},
	{ 3,2,   1.1, -0.80,  0, -1.0, 0.00, "Overcast, moderately graded with azimuthal uniformity"},
	{ 4,2,   1.1, -0.80,  2, -1.5, 0.15, "Overcast, moderately graded with slight brightening towards sun"},
	{ 5,3,   0.0, -1.00,  0, -1.0, 0.00, "Sky of uniform luminance"},
	{ 6,3,   0.0, -1.00,  2, -1.5, 0.15, "Partly cloudy sky, no gradation towards zenith, slight brigtening towards sun"},
	{ 7,3,   0.0, -1.00,  5, -2.5, 0.35, "Partly cloudy sky, no gradation towards zenith, brighter circumsolar region"},
	{ 8,3,   0.0, -1.00, 10, -3.0, 0.45, "Partly cloudy sky, no gradation towards zenith, distinct solar corona"},
	{ 9,4,  -1.0, -0.55,  2, -1.5, 0.15, "Partly cloudy, with the obscured sun"},
	{10,4,  -1.0, -0.55,  5, -2.5, 0.30, "Partly cloudy, with brighter circumsolar region"},
	{11,4,  -1.0, -0.55, 10, -3.0, 0.45, "White-blue sky with distinct solar corona"},
	{12,5,  -1.0, -0.32, 10, -3.0, 0.45, "CIE Standard Clear Sky, low illuminance"},
	{13,5,  -1.0, -0.32, 16, -3.0, 0.30, "CIE Standard Clear Sky, polluted atmosphere"},
	{14,6,  -1.0, -0.15, 16, -3.0, 0.30, "Cloudless turbid sky with broad solar corona"},
	{15,6,  -1.0, -0.15, 24, -2.8, 0.15, "White-blue turbid sky with broad solar corona"},
	{16,3,  -1.0, -0.30,  0, -1.0, 0.00, "Sky of uniform luminance, maximum horizon glow"},
	} ;

float p_half_image_width ; // scaled value of half the image width (use as X origin) ;

#define IMAGE_PIXEL_X_TO_RELATIVE(x) ((float)( ((float)(x)/(float)(IMAGE_WIDTH-1)) -p_half_image_width))
#define IMAGE_PIXEL_Y_TO_RELATIVE(y) ((float)(1.f-(float)(y)/(float)(IMAGE_HEIGHT-1)))

#define IMAGE_RELATIVE_TO_PIXEL_X(px) (int)(((px)+p_half_image_width)*(float)(IMAGE_WIDTH-1)+0.5) ;
#define IMAGE_RELATIVE_TO_PIXEL_Y(py) (int)((1.-(py))*(float)(IMAGE_HEIGHT-1)+0.5) ;


struct V3D {
    float A,B,C ;
    } ;
struct V3D V_zenith;
struct V3D V_sun;

unsigned char *column_mask ; // user requested masked columns
unsigned char *column_repair_mask ; // user requested masked columns

// set half space coefs for horizon
/*  void set_horizon_coefs()  */
/*  {  */
/*      float x0 = 0. ;  */
/*      float y0 = horizon_py ;  */
/*      float x1 = 1 ;  */
/*      float y1 = y0+horizon_slope ;  */
/*      float dx = x1-x0 ;  */
/*      float dy = y1-y0 ;  */
/*      horizon_A = dy ;  */
/*      horizon_B = -dx ;  */
/*      horizon_C = -(horizon_A*x1 + horizon_B*y1) ;  */
/*      float denom = sqrt(dy*dy+dx*dx) ;  */
/*      horizon_A /= denom ;  */
/*      horizon_B /= denom ;  */
/*      horizon_C /= denom ;  */
/*  }  */

void dump_parameters_to_stdout()
{
    int i ;

    for(i=0 ; i < MAX_PARMS ; i++) {
	printf("%d:%12s %f\n", opt_parms[i].optimize_flag, opt_parms[i].name, *opt_parms[i].variable_address) ;
    }

    printf("Maximum vhat=%f\n", maximum_vhat) ;

    printf("command line settings:\n  ") ;

    for(i=0 ; i < MAX_PARMS ; i++) {
	printf(" -%s %f", opt_parms[i].abreviation, *opt_parms[i].variable_address) ;
    }
    printf("\n") ;
}

#define MAX_SAT_FOR_SAMPLE 1.0

#define min_f(a, b, c)  (fminf(a, fminf(b, c)))
#define max_f(a, b, c)  (fmaxf(a, fmaxf(b, c)))

// given sun_x (as an angle from the center of the image), convert
// to px, where px=0 is left side of image, px=1. is right side of image
float sun_x_angle2px(float sun_x_angle)
{
    return sun_x_angle/FOV_horizontal ;
}

void rgb2hsv16(const uint16 src_r, const uint16 src_g, const uint16 src_b, float *h, float *s, float *v)
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

void hsv2rgb16(float h, float s, float v, uint16 *dst_r, uint16 *dst_g, uint16 *dst_b)
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

    *dst_r = (uint16)(r * MAX16); // dst_r : 0-255
    *dst_g = (uint16)(g * MAX16); // dst_r : 0-255
    *dst_b = (uint16)(b * MAX16); // dst_r : 0-255
}


void read_tif_image(tdata_t *image,TIFF *tif,int h,uint16 spp, uint16 tif_config)
{
    int32 y ;
    /* read a row */
    for(y = 0 ; y < h ; y++) {
	if(tif_config == PLANARCONFIG_CONTIG) {
	    TIFFReadScanline(tif,image[y],y,0) ;
	} else if(tif_config == PLANARCONFIG_SEPARATE) {
	    uint16 s ;

	    for(s = 0 ; s < spp ; s++)
		TIFFReadScanline(tif,image[y],y,s) ;

	} else {
	    fprintf(stderr, "Cannot read planar configuration of tiff image, config=%d!\n", (int)tif_config) ;
	    exit(1) ;
	}
    }
}

void copy_pixel(tdata_t *image, int32 xsrc, int32 ysrc, int32 xdest, int32 ydest)
{
    ((uint16 *)(image[ydest]))[IMAGE_NSAMPLES*xdest+0] = ((uint16 *)(image[ysrc]))[IMAGE_NSAMPLES*xsrc+0] ;
    ((uint16 *)(image[ydest]))[IMAGE_NSAMPLES*xdest+1] = ((uint16 *)(image[ysrc]))[IMAGE_NSAMPLES*xsrc+1] ;
    ((uint16 *)(image[ydest]))[IMAGE_NSAMPLES*xdest+2] = ((uint16 *)(image[ysrc]))[IMAGE_NSAMPLES*xsrc+2] ;

    if(IMAGE_HAS_ALPHA) {
	((uint16 *)(image[ydest]))[IMAGE_NSAMPLES*xdest+3] = ((uint16 *)(image[ysrc]))[IMAGE_NSAMPLES*xsrc+3] ;
    }
}

int xy_is_opaque_pixel(tdata_t *image, int32 x, int32 y)
{
    /* if alpha channel is nonzero then hugin has placed image data here */
	if(((uint16 *)(image[y]))[IMAGE_NSAMPLES*x+3] < HALF16)
	    return 0 ;

    return 1 ;
}

int xy_has_nonblack_pixel(tdata_t *image, int32 x, int32 y)
{
    /* if alpha channel is nonzero then hugin has placed image data here */
    if(IMAGE_HAS_ALPHA) {
	if(((uint16 *)(image[y]))[IMAGE_NSAMPLES*x+3] < HALF16)
	    return 0 ;
    }

    if(((uint16 *)(image[y]))[IMAGE_NSAMPLES*x+0] > BLK16)
	return 1 ;
    if(((uint16 *)(image[y]))[IMAGE_NSAMPLES*x+1] > BLK16)
	return 1 ;
    if(((uint16 *)(image[y]))[IMAGE_NSAMPLES*x+2] > BLK16)
	return 1 ;

    return 0 ;
}

void simple_find_start_of_sky(int16 *start_of_sky,uint16 w,uint16 h,tdata_t *image)
{
    int32 x, y ;

    for(x = 0 ; x < w ; x++) {
	start_of_sky[x] = 0 ;

	if(column_mask[x] == 1) {
	    end_of_sky[x] = 0 ;
	    continue ;
	}

	end_of_sky[x] = IMAGE_HEIGHT-1 ;

	for(y = 0 ; y < h ; y++) {
	    if(xy_has_nonblack_pixel(image, x, y) == 1) {
		start_of_sky[x] = y ;
		break ;
	    }
	}
    }

    fprintf(stderr, "Finished simple find start of sky\n") ;
}

void find_start_of_sky(int16 *start_of_sky,uint16 w,uint16 h,tdata_t *image,int fix_SOS_edges)
{
    int32 x, y ;

    simple_find_start_of_sky(start_of_sky,w,h,image) ;
    return ;

    /* look for vertical gaps in sky (alpha=0), which can happen after lens correction */

    for(x = 0 ; x < w ; x++) {
	int32 gap_start=-1, gap_end=-1 ;

	for(y = start_of_sky[x] ; y < h ; y++) {
	    if(xy_has_nonblack_pixel(image, x, y) == 0) {
		gap_start = y ;
		break ;
	    }
	}

	if(gap_start > -1) {
	    for(y = gap_start ; y < h ; y++) {
		/* if alpha channel is zero then hugin has placed image data here */
		if(xy_has_nonblack_pixel(image, x, y) == 1) {
		    gap_end = y-1 ;
		    break ;
		}
	    }
	}

	if(gap_end> -1 && gap_start > -1) {
	    /* interpolate gap values from start to end */
	    float start_r,start_g,start_b ;
	    float end_r,end_g,end_b ;
	    tif_get3c(image,x,gap_start-1,start_r,start_g,start_b) ;
	    tif_get3c(image,x,gap_end+1,start_r,start_g,start_b) ;

	    for(y = gap_start ; y < gap_end+1 ; y++) {
		float pStart = 1.f - ((float)y-((float)gap_start-1.f))/((float)gap_end-(float)gap_start+3.f) ;
		float pEnd = 1.f - pStart ;

		tif_set4c(image,x,y, (uint16)(pStart*start_r + pEnd*end_r), (uint16)(pStart*start_g + pEnd*end_g), (uint16)(pStart*start_b + pEnd*end_b), MAX16) ; 
	    }
	}

    }

    for(x = 0 ; x < w ; x++) {
	if(column_mask[x] == 1) start_of_sky[x]=0 ;
    }


    if(fix_SOS_edges == 1) {
	fprintf(stderr, "Fixing start of sky edges\n") ;

	// eliminate long vertical edges in the sky end, which can create
	// visible edges in the feathering result
	for(x = 0 ; x < w-1 ; x++) {
	    if(column_mask[x] == 1) continue ;
	    if(column_mask[x+1] == 1) continue ;

	    if( (start_of_sky[x] - start_of_sky[x+1]) > 2) {
		start_of_sky[x+1] = start_of_sky[x]-2 ;
	    }
	}

	for(x = w-1 ; x > 0 ; x--) {
	    if(column_mask[x] == 1) continue ;
	    if(column_mask[x-1] == 1) continue ;

	    if( (start_of_sky[x] - start_of_sky[x-1]) > 2) {
		start_of_sky[x-1] = start_of_sky[x]-2 ;
	    }
	}

    }
}


void get_sky_mean_var(tdata_t *image, int x, int y0, int n, float mean[], float var[])
{
    float sum[3] = {0., 0., 0.} ;  // for r,g,b
    float sumsq[3] = {0., 0., 0.} ;  // for r,g,b
    int y ;

    for(y = y0 ; y < y0+n ; y++) {
	uint16 r,g,b ;
	tif_get3c(image,x,y,r,g,b) ;

	sum[0] += r ;
	sum[1] += g ;
	sum[2] += b ;
	sumsq[0] += r*r ;
	sumsq[1] += g*g ;
	sumsq[2] += b*b ;
    }

    mean[0] = sum[0]/(float)n ;
    mean[1] = sum[1]/(float)n ;
    mean[2] = sum[2]/(float)n ;

    var[0] = sumsq[0] - sum[0]*mean[0] ; // equiv to sumsq - sum*sum/n ;
    var[1] = sumsq[1] - sum[1]*mean[1] ; // equiv to sumsq - sum*sum/n ;
    var[2] = sumsq[2] - sum[2]*mean[2] ; // equiv to sumsq - sum*sum/n ;
}


float *sky_hue, *sky_val, *sky_sat  ;
uint8 *is_clear_sky ;
int n_sky ;

int is_end_of_sky(uint16 y, int n_samples, float *hdiff, float *sdiff, float *vdiff)
{
    float m_hue0=0. ;
    float m_val0=0. ;
    float m_sat0=0. ;
    float m_hue1=0. ;
    float m_val1=0. ;
    float m_sat1=0. ;
    int iy ;

    if(1.f-(float)y/((float)IMAGE_HEIGHT-1.) > min_sky_end_p) return 0 ;  // is it greater than the min_sky_end_p requested?
    float sum_wgts=0. ;

    for(iy=y-n_samples+1 ; iy <= y ; iy++) {
	if(iy < 0 || iy > IMAGE_HEIGHT-1) {
	    fprintf(stderr, "iy=%d, y=%d n_samples=%d\n", (int)iy, (int)y, n_samples) ;
	    exit(1) ;
	}
	m_hue0 += sky_hue[iy]*is_clear_sky[iy] ;
	m_sat0 += sky_sat[iy]*is_clear_sky[iy] ;
	m_val0 += sky_val[iy]*is_clear_sky[iy] ;
	sum_wgts += is_clear_sky[iy] ;
    }

    m_hue0 /= sum_wgts ;
    m_sat0 /= sum_wgts ;
    m_val0 /= sum_wgts ;

    sum_wgts=0. ;
    for(iy=y+1 ; iy < y+1+n_samples/2. ; iy++) {
	m_hue1 += sky_hue[iy]*is_clear_sky[iy] ;
	m_sat1 += sky_sat[iy]*is_clear_sky[iy] ;
	m_val1 += sky_val[iy]*is_clear_sky[iy] ;
	sum_wgts += is_clear_sky[iy] ;
    }

    m_hue1 /= sum_wgts ;
    m_sat1 /= sum_wgts ;
    m_val1 /= sum_wgts ;

    float huediff = fabs(m_hue0 - m_hue1) ;
    float satdiff = fabs(m_sat0 - m_sat1) ;
    float valdiff = fabs(m_val0 - m_val1) ;

    huediff *= (m_sat0+m_sat1)/2. ;   // adjust hue difference by amount of saturation
    *hdiff = huediff ;
    *sdiff = satdiff ;
    *vdiff = valdiff ;


    if(huediff > sky_hue_tolerance) {
	//fprintf(stderr, "ieos %d %3.0f %0.3f %0.3f\n", y, huediff, satdiff, valdiff) ;
	return 1 ;
    }

    if(valdiff > sky_val_tolerance) {
	//fprintf(stderr, "ieos %d %3.0f %0.3f %0.3f\n", y, huediff, satdiff, valdiff) ;
	return 1 ;
    }

    if(satdiff > sky_sat_tolerance) {
	//fprintf(stderr, "ieos %d %3.0f %0.3f %0.3f\n", y, huediff, satdiff, valdiff) ;
	return 1 ;
    }

    return 0 ;
}

float sobel(tdata_t *image, int xc, int yc)
{
    uint16 r[3][3] ;
    uint16 g[3][3] ;
    uint16 b[3][3] ;

    for(int x = xc-1 ; x <= xc+1 ; x++) {
	int i = x - (xc-1) ;
	for(int y = yc-1 ; y <= yc+1 ; y++) {
	    int j = y - (yc-1) ;
	    tif_get3c(image,x,y,r[i][j],g[i][j],b[i][j]) ;
	}
    }

    float Gx_r = 
	-  r[0][0] +   r[0][2]
	-2*r[1][0] + 2*r[0][2]
	-  r[2][0] +   r[2][2] ;
    float Gx_g = 
	-  g[0][0] +   g[0][2]
	-2*g[1][0] + 2*g[0][2]
	-  g[2][0] +   g[2][2] ;
    float Gx_b = 
	-  b[0][0] +   b[0][2]
	-2*b[1][0] + 2*b[0][2]
	-  b[2][0] +   b[2][2] ;

    float Gy_r = 
	    r[0][0] + 2*r[0][1] + r[0][2]
	   -r[2][0] - 2*r[2][1] - r[2][2] ;
    float Gy_g = 
	    g[0][0] + 2*g[0][1] + g[0][2]
	   -g[2][0] - 2*g[2][1] - g[2][2] ;
    float Gy_b = 
	    b[0][0] + 2*b[0][1] + b[0][2]
	   -b[2][0] - 2*b[2][1] - b[2][2] ;

    return fabs(Gx_r) + fabs(Gx_g) + fabs(Gx_b) + (fabs(Gy_r) + fabs(Gy_g) + fabs(Gy_b))/2. ;
}


void print_sobel(tdata_t *image, int x)
{
    int16 y ;

    for(y = end_of_sky[x]-5 ; y < end_of_sky[x]+5 ; y++) {
	if( y == end_of_sky[x]) {
	    fprintf(stderr, "SOBEL col %d,%d * %f\n", x, y, sobel(image,x,y)) ;
	} else {
	    fprintf(stderr, "SOBEL col %d,%d   %f\n", x, y, sobel(image,x,y)) ;
	}
    }
    exit(1) ;
}

int get_sobel_eos(tdata_t *image, int x)
{
    int y ;

    if(0) {
	if(x < 1) x = 1 ;
	if(x > IMAGE_WIDTH-2) x = IMAGE_WIDTH-2 ;
	int best_y = 0 ;
	float best_sobel=0. ;

	for(y = end_of_sky[x]-5 ; y < end_of_sky[x]+5 ; y++) {
	    if(y < 1 || y > IMAGE_HEIGHT-2) {
		fprintf(stderr, "get_sobel: y:%d not between 1 and %d\n", y, IMAGE_HEIGHT-2) ;
		continue ;
	    }
	    float this = sobel(image,x,y) ;
	    if(this > best_sobel) {
		best_sobel = this ;
		best_y = y ;
	    }
	}

	if(best_y == 0) return end_of_sky[x] ;

	return best_y-1 ;
    }

    // trying out new end of sky detection
    float obuf[20] ; 
    //current eos-10 will be obuf index 0 ;
    int offset = end_of_sky[x]-10 ;
    int y0 = end_of_sky[x]-10 ;
    int y1 = y0+19 ;

    int sample_length=3 ;
    for(y = y0 ; y <= y1 ; y++) {
	float r0,g0,b0,a0 ;
	float r1,g1,b1,a1 ;
	r0=g0=b0=a0=0. ;
	r1=g1=b1=a1=0. ;
	float r,g,b ;

	if(y-sample_length > 0 && y+sample_length <= IMAGE_HEIGHT-1) {
	    for(int z = y-sample_length ; z < y ; z++) {
		tif_get3c(image,x,z,r,g,b) ;
		r0 += r ;
		g0 += g ;
		b0 += b ;
	    }

	    for(int z = y+1 ; z <= y+sample_length ; z++) {
		tif_get3c(image,x,z,r,g,b) ;
		r1 += r ;
		g1 += g ;
		b1 += b ;
	    }

	    r0 /= (float)sample_length ;
	    r1 /= (float)sample_length ;
	    g0 /= (float)sample_length ;
	    g1 /= (float)sample_length ;
	    b0 /= (float)sample_length ;
	    b1 /= (float)sample_length ;

	    float dr = fabs(r1-r0) ;
	    float dg = fabs(g1-g0) ;
	    float db = fabs(b1-b0) ;

	    obuf[y-offset] = dr+dg+db ;
	} else {
	    obuf[y-offset] = 0. ;
	}

    }

    int eos = 0 ;
    float val_tol = sky_val_tolerance/.015* 4000. ;
    for(y = y0 ; y <= y1  ; y++) {
	if(obuf[y-offset] > val_tol) {
	    eos = y+sample_length-1 ;
	    return eos ;
	}
    }

    return end_of_sky[x] ;

}



void find_end_of_sky(int16 *end_of_sky,int16 *start_of_sky,int w,int h,tdata_t *image,int fix_edges,int fix_slivers)
{
    int16 x, y ;

    is_clear_sky = (uint8 *)calloc(IMAGE_HEIGHT, sizeof(uint8)) ;
    sky_hue = (float *)calloc(IMAGE_HEIGHT, sizeof(float)) ;
    sky_val = (float *)calloc(IMAGE_HEIGHT, sizeof(float)) ;
    sky_sat = (float *)calloc(IMAGE_HEIGHT, sizeof(float)) ;

    // look for end of sky now

    for(x = 0 ; x < w ; x++) {

	if(column_mask[x] == 1) {
	    end_of_sky[x] = 0 ;
	    continue ;
	}

	float hdiff, sdiff, vdiff ;

	int first_y = raw_start_of_sky[x]+10 ;  // skip first 10 pixels

	//first_y = end_of_sky[x]-6 ;  // temporarily see how this looks 
	float top_mean[3], top_var[3] ;
	float bot_mean[3], bot_var[3] ;

	//fprintf(stderr, "\nNEW EOS x:%5d eos:%5d\n", x, end_of_sky[x]) ;

	for(y = first_y ; y < IMAGE_HEIGHT-10 ; y++) {
	    get_sky_mean_var(image, x, y-10, 10, top_mean, top_var) ;
	    get_sky_mean_var(image, x, y+5, 5, bot_mean, bot_var) ;
	    float rr = fabs(top_mean[0]/bot_mean[0]-1.) ;
	    float rg = fabs(top_mean[1]/bot_mean[1]-1.) ;
	    float rb = fabs(top_mean[2]/bot_mean[2]-1.) ;

	    // do we have the same color above and below the test area?
	    if(rr > .1 || rg > .1 || rb > .1) {
		for(; y < IMAGE_HEIGHT-10 ; y++) {
		    is_clear_sky[y]=1 ;
		}
		break ;
	    }

	    is_clear_sky[y]=1 ;

	    // if transparent or black, must be at bottom of image
	    uint16 r,g,b,a ;
	    tif_get4c(image,x,y+10,r,g,b,a) ;

	    if(a < HALF16) break ;
	    if( r < BLK16 && g < BLK16 && b < BLK16 ) break ;

	    int last_y = y ;


	    for(int y0 = y ; y0 < y+5 ; y0++) {
		// try to skip small highlights (jet trails) in sky
		tif_get3c(image,x,y,r,g,b) ;

		float rtop=(float)r/top_mean[0] ;
		float gtop=(float)g/top_mean[1] ;
		float btop=(float)b/top_mean[2] ;
		float rbot=(float)r/bot_mean[0] ;
		float gbot=(float)g/bot_mean[1] ;
		float bbot=(float)b/bot_mean[2] ;

		if( rtop > 1. && gtop > 1. && btop > 1. && rbot > 1. && gbot > 1. && bbot > 1.) {
		    is_clear_sky[y0] = 0 ;
		    last_y = y0 ;
		} else {
		    is_clear_sky[y0] = 1 ;
		}
	    }

	    y = last_y ;

	}


	// load up arrays

	for(y = first_y ; y < first_y+10 ; y++) {
	    uint16 r,g,b ;
	    tif_get3c(image,x,y,r,g,b) ;

	    float h,s,v ;
	    rgb2hsv16(r,g,b,&h,&s,&v) ;
	    float l = (.299*(float)r+.587*(float)g+.114*(float)b)/MAX16f ;
	    sky_hue[y] = h ;
	    sky_sat[y] = s ;
	    sky_val[y] = l ;
	}

	for(y = first_y ; y < h-10 ; y++) {
	    end_of_sky[x] = y+4 ;

	    if(is_end_of_sky(y+4,5,&hdiff,&sdiff,&vdiff)) {
		break ;
	    }

	    uint16 r,g,b ;
	    tif_get3c(image,x,y+10,r,g,b) ;

	    float h,s,v ;
	    rgb2hsv16(r,g,b,&h,&s,&v) ;
	    float l = (.299*(float)r+.587*(float)g+.114*(float)b)/MAX16f ;
	    sky_hue[y+10] = h ;
	    sky_sat[y+10] = s ;
	    sky_val[y+10] = l ;
	}

	// fine tune the detection -- increase y as long as vdiff becomes larger
	float vdiff_prev=vdiff ;

	for(y = end_of_sky[x]+1 ; y < h-10 ; y++) {
	    uint16 r,g,b ;
	    tif_get3c(image,x,y+5,r,g,b) ;

	    float h,s,v ;
	    rgb2hsv16(r,g,b,&h,&s,&v) ;

	    float l = (.299*(float)r+.587*(float)g+.114*(float)b)/MAX16f ;
	    sky_hue[y+5] = h ;
	    sky_sat[y+5] = s ;
	    sky_val[y+5] = l ;

	    tif_get3c(image,x,y+6,r,g,b) ;

	    rgb2hsv16(r,g,b,&h,&s,&v) ;

	    l = (.299*(float)r+.587*(float)g+.114*(float)b)/MAX16f ;
	    sky_hue[y+6] = h ;
	    sky_sat[y+6] = s ;
	    sky_val[y+6] = l ;

	    is_end_of_sky(y,5,&hdiff,&sdiff,&vdiff) ;

	    //if(x == 232)
		//printf("fty:sdiff:vdiff %d:%f:%f\n", y, sdiff, vdiff) ;

	    if(vdiff <= vdiff_prev)
		break ;

	    vdiff_prev = vdiff ;

	    end_of_sky[x]++ ;
	}

	end_of_sky[x]-=2 ;
	end_of_sky[x] = get_sobel_eos(image,x) ;
    }

    for(x = w ; x < w ; x++) {
	if(column_mask[x] == 1) continue ;

	// load up arrays
	int first_y = raw_start_of_sky[x]+10 ;  // skip first 10 pixels

	//first_y = end_of_sky[x]-6 ;  // temporarily see how this looks 
	float top_mean[3], top_var[3] ;
	float bot_mean[3], bot_var[3] ;
	float max_diff_mean=0 ;
	float n_diff_means=0. ;
	float sum_diff_means=0. ;
	int best_y=first_y ;

	//fprintf(stderr, "\nNEW EOS x:%5d eos:%5d\n", x, end_of_sky[x]) ;

	for(y = first_y ; y < IMAGE_HEIGHT-10 ; y++) {
	    get_sky_mean_var(image, x, y-10, 10, top_mean, top_var) ;
	    get_sky_mean_var(image, x, y+1, 10, bot_mean, bot_var) ;

	    // if transparent or black, must be at bottom of image
	    uint16 r,g,b,a ;
	    tif_get4c(image,x,y+10,r,g,b,a) ;

	    if(a < HALF16) break ;
	    if( r < BLK16 && g < BLK16 && b < BLK16 ) break ;

	    // try to skip small highlights (jet trails) in sky
	    tif_get3c(image,x,y,r,g,b) ;

	    float rtop=(float)r/top_mean[0] ;
	    float gtop=(float)g/top_mean[1] ;
	    float btop=(float)b/top_mean[2] ;
	    float rbot=(float)r/bot_mean[0] ;
	    float gbot=(float)g/bot_mean[1] ;
	    float bbot=(float)b/bot_mean[2] ;
	    if( rtop > 1. && rbot > 1. && gtop > 1. && gbot > 1. && btop > 1. && bbot > 1.)
		continue ;

	    //float diff_mean = .299*fabs(top_mean[0]-bot_mean[0]) + .587*fabs(top_mean[1]-bot_mean[1]) + .114*fabs(top_mean[2]-bot_mean[2]) ;
	    float diff_mean = fabs(top_mean[0]-bot_mean[0]) + fabs(top_mean[1]-bot_mean[1]) + fabs(top_mean[2]-bot_mean[2]) ;
	    float rdiff=1. ;

	    if(n_diff_means > 10.) rdiff = diff_mean/(sum_diff_means/n_diff_means) ;

	    float prob = rdiff - 1. ;
	    if(prob < 0.) prob=0 ;



	    // probability it's a true diff is quite low in top of image
	    //float prob = (float)y/(float)(IMAGE_HEIGHT/2.) ;
	    //if(prob > 1.) prob = 1. ;

	    //float tot_top_var = top_var[0] + top_var[1] + top_var[2] ;
	    //float tot_bot_var = bot_var[0] + bot_var[1] + bot_var[2] ;

	    if(max_diff_mean < diff_mean*prob) {
		max_diff_mean = diff_mean*prob ;
		best_y=y ;
	    }

	    sum_diff_means += diff_mean ;
	    n_diff_means += 1. ;

	    //fprintf(stderr, "        y:%5d %10.0f   %10.0f   %10.0f\n", y, diff_mean, tot_top_var, tot_bot_var) ;
	}

	end_of_sky[x] = best_y ;
	end_of_sky[x] = get_sobel_eos(image,x) ;
    }

    //look for and fix slivers
    fix_slivers=0 ;
    while(fix_slivers >0) {
	fprintf(stderr, "Fixing slivers\n") ;

	for(x = 1 ; x < w-3 ; x++) {
	    if(column_mask[x] == 1) continue ;
	    if(column_mask[x-1] == 1) continue ;

	    if( (end_of_sky[x-1] - end_of_sky[x]) > 10 && (end_of_sky[x+1] - end_of_sky[x]) > 10) {
		if(column_mask[x+1] == 1) continue ;
		//fprintf(stderr,"fix sliver at x=%d\n", x) ;
		end_of_sky[x] = (int)((end_of_sky[x-1]+end_of_sky[x+1])/2+0.5) ;

	    } else if( (end_of_sky[x-1] - end_of_sky[x]) > 10 && (end_of_sky[x+2] - end_of_sky[x]) > 10) {
		if(column_mask[x+2] == 1) continue ;
		//fprintf(stderr,"fix 2-wide sliver at x=%d\n", x) ;
		end_of_sky[x] = (int)((end_of_sky[x-1]+end_of_sky[x+2])/2+0.5) ;
		end_of_sky[x+1] = (int)((end_of_sky[x-1]+end_of_sky[x+2])/2+0.5) ;

	    } else if( (end_of_sky[x-1] - end_of_sky[x]) > 10 && (end_of_sky[x+3] - end_of_sky[x]) > 10) {
		if(column_mask[x+3] == 1) continue ;
		//fprintf(stderr,"fix 3-wide sliver at x=%d\n", x) ;
		end_of_sky[x] = (int)((end_of_sky[x-1]+end_of_sky[x+3])/2+0.5) ;
		end_of_sky[x+1] = (int)((end_of_sky[x-1]+end_of_sky[x+3])/2+0.5) ;
		end_of_sky[x+2] = (int)((end_of_sky[x-1]+end_of_sky[x+3])/2+0.5) ;
	    }
	}
	fix_slivers-- ;
    }

    if(0 && fix_edges == 1) {
	fprintf(stderr, "Fixing edges\n") ;

	// eliminate long vertical edges in the sky end, which can create
	// visible edges in the feathering result
	for(x = 0 ; x < w-1 ; x++) {
	    if(column_mask[x] == 1) continue ;
	    if(column_mask[x+1] == 1) continue ;

	    if( (end_of_sky[x+1] - end_of_sky[x]) > 2) {
		//fprintf(stderr, "set EOS %d to %d\n", x+1, end_of_sky[x]+2) ;
		end_of_sky[x+1] = end_of_sky[x]+2 ;
	    }
	}

	for(x = w-1 ; x > 0 ; x--) {
	    if(column_mask[x] == 1) continue ;
	    if(column_mask[x-1] == 1) continue ;

	    if( (end_of_sky[x-1] - end_of_sky[x]) > 2) {
		end_of_sky[x-1] = end_of_sky[x]+2 ;
	    }
	}
    }

    // set lowest_sky_py_found
    for(x = 0 ; x < w-1 ; x++) {
	if(column_mask[x] == 1) continue ;
	float py = 1.-(float)end_of_sky[x]/(float)(h-1) ;

	if(lowest_sky_py_found > py) {
	    lowest_sky_py_found = py ;
	}
    }

/*      print_sobel(image,w/2) ;  */

    free(is_clear_sky) ;
    free(sky_hue) ;
    free(sky_sat) ;
    free(sky_val) ;
}

int find_x1(int x0, int desired_width, int max_width, int w)
{

    int prev_x1=x0 ;
    int x1 ;

    for(x1=x0+1 ; x1 < w ; x1++) {
	if(column_mask[x1] == 0) {
	    if((x1-x0) >= desired_width || x1 == w-1) {
		if((x1-x0) > max_width && prev_x1 > x0) {
		    x1=prev_x1 ;
		}
		break ;
	    } else {
		prev_x1 = x1 ;
	    }
	}
    }

    return x1 ;
}

// return luminance for image data at pixel x,y
float get_xy_L(tdata_t *image, int x, int y)
{
    if(y < 0) y = 0 ;
    if(y > IMAGE_HEIGHT-1) y = IMAGE_HEIGHT-1 ;
    uint16 r,g,b ;
    tif_get3c(image,x,y,r,g,b) ;

    float fr = (float)r/MAX16f ;
    float fg = (float)g/MAX16f ;
    float fb = (float)b/MAX16f ;
    return 0.21126*fr + 0.7152*fg + 0.0722*fb ;
}

void repair_sky_hue(tdata_t *image,int16 *start_of_sky,int16 *end_of_sky)
{
    int x, y ;

    fprintf(stderr, "repair sky hue\n") ;

    // need some correction factors
    double s_coefs_left[2] ;
    double h_coefs_left[2] ;
    double s_coefs_right[2] ;
    double h_coefs_right[2] ;

    float sum_s=0. , sum_h=0. ;
    float sum_shat=0. , sum_hhat=0. ;
    int n=0, n_l=0, n_r=0 ;

    struct mstat m_fit_s_l = init_reg(1) ; // will be a simple linear in y
    struct mstat m_fit_h_l = init_reg(1) ; // will be a simple linear in y
    struct mstat m_fit_s_r = init_reg(1) ; // will be a simple linear in y
    struct mstat m_fit_h_r = init_reg(1) ; // will be a simple linear in y

    fprintf(stderr, "Repair sky hue from %d to %d\n", min_sky_hue_mask, max_sky_hue_mask) ;

    for(x=min_sky_hue_mask-10 ; x < min_sky_hue_mask ; x ++) {
	if(x < 0) continue ;
	if(x > IMAGE_WIDTH-1) continue ;
	if(column_mask[x] == 1) continue ;

	for(y = start_of_sky[x] ; y <= end_of_sky[x] ; y++) {
	    if(y > IMAGE_HEIGHT-1) break ;

	    uint16 r,g,b ;
	    tif_get3c(image,x,y,r,g,b) ;

	    float h,s,v ;
	    rgb2hsv16(r,g,b,&h,&s,&v) ;

	    float px = IMAGE_PIXEL_X_TO_RELATIVE(x) ;
	    float py = IMAGE_PIXEL_Y_TO_RELATIVE(y) ;
	    float hhat, shat ;
	    predict_sky_huesat_from_val(v, &hhat, &shat, px, py) ;

	    if(shat > .1) {
		sum_s += s ;
		sum_h += h ;
		double xv[1] = {y} ;
		sum_reg(&m_fit_s_l, xv, s/shat) ;
		sum_reg(&m_fit_h_l, xv, h/hhat) ;
		sum_shat += shat ;
		sum_hhat += hhat ;

		n_l++ ;

		n++ ;
	    }
	}
    }

    for(x=max_sky_hue_mask ; x < max_sky_hue_mask+10 ; x ++) {
	if(x < 0) continue ;
	if(x > IMAGE_WIDTH-1) break ;
	if(column_mask[x] == 1) continue ;

	for(y = start_of_sky[x] ; y <= end_of_sky[x] ; y++) {
	    if(y > IMAGE_HEIGHT-1) break ;

	    uint16 r,g,b ;
	    tif_get3c(image,x,y,r,g,b) ;

	    float h,s,v ;
	    rgb2hsv16(r,g,b,&h,&s,&v) ;

	    float px = IMAGE_PIXEL_X_TO_RELATIVE(x) ;
	    float py = IMAGE_PIXEL_Y_TO_RELATIVE(y) ;
	    float hhat, shat ;
	    predict_sky_huesat_from_val(v, &hhat, &shat, px, py) ;

	    if(shat > .1) {
		sum_s += s ;
		sum_h += h ;
		double xv[1] = {y} ;
		sum_reg(&m_fit_s_r, xv, s/shat) ;
		sum_reg(&m_fit_h_r, xv, h/hhat) ;
		sum_shat += shat ;
		sum_hhat += hhat ;

		n_r++ ;

		n++ ;
	    }
	}
    }

    print_reg_matrix(&m_fit_s_l, "m_fit_h_l") ;
    print_reg_matrix(&m_fit_s_l, "m_fit_s_l") ;
    print_reg_matrix(&m_fit_s_r, "m_fit_s_r") ;

    float s_correct=1. ;
    float h_correct=1. ;

    if(n > 0) {
	s_correct = sum_s/sum_shat ;
	h_correct = sum_h/sum_hhat ;
    }

    if(n_l > 2) {
	estimate_reg(&m_fit_s_l, s_coefs_left) ;
	estimate_reg(&m_fit_h_l, h_coefs_left) ;
    } else {
	s_coefs_left[0] = s_correct ;
	s_coefs_left[1] = 0.0 ;
    }

    if(n_r > 2) {
	estimate_reg(&m_fit_s_r, s_coefs_right) ;
	estimate_reg(&m_fit_h_r, h_coefs_right) ;
    } else {
	s_coefs_right[0] = s_correct ;
	s_coefs_right[1] = 0.0 ;
    }

    if(n_l > 2 && n_r <= 2) {
	s_coefs_right[0] = s_coefs_left[0] ;
	s_coefs_right[1] = s_coefs_left[1] ;
	h_coefs_right[0] = h_coefs_left[0] ;
	h_coefs_right[1] = h_coefs_left[1] ;
    }

    if(n_r > 2 && n_l <= 2) {
	s_coefs_left[0] = s_coefs_right[0] ;
	s_coefs_left[1] = s_coefs_right[1] ;
	h_coefs_left[0] = h_coefs_right[0] ;
	h_coefs_left[1] = h_coefs_right[1] ;
    }

    fprintf(stderr, "repair sky hue, h,s correct = %f, %f\n", h_correct, s_correct) ;
    fprintf(stderr, "repair sky hue, h_coefs(l,r) = (%f, %f) (%f,%f)\n", h_coefs_left[0], h_coefs_left[1], h_coefs_right[0], h_coefs_right[1]) ;
    fprintf(stderr, "repair sky hue, s_coefs(l,r) = (%f, %f) (%f,%f)\n", s_coefs_left[0], s_coefs_left[1], s_coefs_right[0], s_coefs_right[1]) ;



    for(x=min_sky_hue_mask ; x < max_sky_hue_mask ; x ++) {
	uint16 r,g,b ;
	if(column_mask[x] == 1) continue ;

	float pleft = 1.-(float)(x-min_sky_hue_mask)/(float)(max_sky_hue_mask-min_sky_hue_mask) ;


	for(y = start_of_sky[x] ; y <= end_of_sky[x] ; y++) {
	    if(y >= IMAGE_HEIGHT) break ;
	    tif_get3c(image,x,y,r,g,b) ;

	    float h,s,v ;
	    rgb2hsv16(r,g,b,&h,&s,&v) ;

	    float px = IMAGE_PIXEL_X_TO_RELATIVE(x) ;
	    float py = IMAGE_PIXEL_Y_TO_RELATIVE(y) ;

	    predict_sky_huesat_from_val(v, &h, &s, px, py) ;

	    //compute correction factors on left and right side
	    float h_correct_l = h_coefs_left[0] + h_coefs_left[1]*(float)y ;
	    float s_correct_l = s_coefs_left[0] + s_coefs_left[1]*(float)y ;
	    float h_correct_r = h_coefs_right[0] + h_coefs_right[1]*(float)y ;
	    float s_correct_r = s_coefs_right[0] + s_coefs_right[1]*(float)y ;

	    // interpolate correction factors
	    float pixel_h_correct = pleft*h_correct_l + (1.-pleft)*h_correct_r ;
	    float pixel_s_correct = pleft*s_correct_l + (1.-pleft)*s_correct_r ;


	    if(x == min_sky_hue_mask) {
		fprintf(stderr, "repair sky hue, y=%d, h,s correct = %f, %f\n", y, pixel_h_correct, pixel_s_correct) ;
	    }

	    h *= pixel_h_correct ;
	    s *= pixel_s_correct ;

	    hsv2rgb16(h,s,v,&r,&g,&b) ;
	    tif_set3c(image,x,y,r,g,b) ;

	}
    }
}

struct sample_point {
    float px, py ;
    float h,s,v ;
    float v_hat ;
    uint16 r,g,b ;
} ;

struct sample_point *samples = NULL ;
int n_samples_to_optimize = 0 ;

int sat_prediction_method = 1 ; // 1 means predict as function of (val,val*va), 0 means just a mean value is used
float S_from_V_coefs[3] = {0.5,0.,0.} ; // coefs to predict sat from value
float H_from_V_coefs[4] = {212,0.,0.,0.} ; // coefs to predict hue from value

int sample_sky_points(int n_per_column, int n_columns,tdata_t *image,int16 *start_of_sky,int16 *end_of_sky)
{
    int x, y ;
    float mean_column_length ;
    float sum_column_lengths=0. ;
    int n_columns_in_sum=0 ;

    for(x=0 ; x < IMAGE_WIDTH ; x++) {
	if(column_mask[x] == 1) continue ;

	 sum_column_lengths += end_of_sky[x] - start_of_sky[x] + 1 ;
	 n_columns_in_sum++ ;
    }

    mean_column_length = sum_column_lengths/(float)n_columns_in_sum ;

    // try to get 1000 sample points, each will be a mean of a 5x5 area around the actual point
    // do this with a rectangular grid over the image area, refining the grid to the sky area only

    n_columns = 50 ;
    int n_rows = (int)(1000./(float)n_columns * (float)IMAGE_HEIGHT/mean_column_length+0.5) ;
    int dx = IMAGE_WIDTH/n_columns ;
    int dy = IMAGE_HEIGHT/n_rows ;
    int n_samples=0 ;

    // how many valid samples in the grid?
    for(x=0 ; x < IMAGE_WIDTH ; x += dx) {
	if(column_mask[x] == 1) continue ;

	for(y=0 ; y < IMAGE_HEIGHT ; y += dy) {

	    if(y >= start_of_sky[x] && y <= end_of_sky[x])
		n_samples++ ;
	}
    }

    printf("Using a grid over the sky with %d samples\n", n_samples) ;

    if(samples == NULL) {
	samples = (struct sample_point *) calloc(n_samples+100, sizeof(struct sample_point)) ;
    }

    int max_samples = n_samples+100 ;

    n_samples=0 ;

    float sum_val_wgts = 0. ;
    float sum_val=0. ;
    float sum_hue_wgts = 0. ;
    float sum_hue=0. ;
    float sum_sat_wgts = 0. ;
    float sum_sat=0. ;
    int n_max_rows= IMAGE_HEIGHT/dy ;
    float sum_h_by_row[n_max_rows] ;
    int n_h_by_row[n_max_rows] ;

    for(int i = 0 ; i < n_max_rows ; i++) {
	sum_h_by_row[i] = 0. ;
	n_h_by_row[i] = 0 ;
    }


    for(x=0 ; x < IMAGE_WIDTH ; x += dx) {
	uint16 r,g,b ;

	if(column_mask[x] == 1) continue ; // don't attempt a mean sky sample inside a masked area
	if(fix_sky_hue && x>=min_sky_hue_mask && x <= min_sky_hue_mask) continue ;

	float sum_h=0. ;
	float sum_s=0. ;
	float sum_v=0. ;
	float sum_n=0. ;
	int row=0 ;

	for(y=0 ; y < IMAGE_HEIGHT ; y += dy) {

	    if(y >= start_of_sky[x] && y <= end_of_sky[x]) {
		int sx, sy ;
		int n=0 ;
		float sum_r=0. ;
		float sum_g=0. ;
		float sum_b=0. ;

		if(n_samples == max_samples) {
		    fprintf(stderr, "ACCCK!  sample_sky_points(), too many samples\n") ;
		    exit(0) ;
		}

		int x0=x-2 ;
		int x1=x+3 ;
		int y0=y-2 ;
		int y1=y+3 ;
		if(x0 < 0) x0 = 0 ;
		if(y0 < 0) y0 = 0 ;
		if(x1 > IMAGE_WIDTH) x1 = IMAGE_WIDTH ;
		if(y1 > IMAGE_HEIGHT) y1 = IMAGE_HEIGHT ;

		for(sx = x0 ; sx < x1 ; sx++) {
		    if(column_mask[x] == 1) continue ; // no sky sample from inside a masked area
		    for(sy = y0 ; sy < y1 ; sy++) {
			if(sy >= start_of_sky[sx] && sy <= end_of_sky[sx]) {
			    tif_get3c(image,x,y,r,g,b) ;

			    sum_r += r ;
			    sum_g += g ;
			    sum_b += b ;
			    n++ ;
			}
		    }
		}

		if(n == 0) continue ; // n == 0 could potentially happen near edges of detected sky area

		r = (uint16)(sum_r/(float)n) ;
		g = (uint16)(sum_g/(float)n) ;
		b = (uint16)(sum_b/(float)n) ;

		float py = IMAGE_PIXEL_Y_TO_RELATIVE(y) ;
		float px = IMAGE_PIXEL_X_TO_RELATIVE(x) ;

		float h,s,v ;
		rgb2hsv16(r,g,b,&h,&s,&v) ;
		samples[n_samples].px = px ;
		samples[n_samples].py = py ;
		samples[n_samples].r = r ;
		samples[n_samples].g = g ;
		samples[n_samples].b = b ;

		samples[n_samples].h = h ;
		samples[n_samples].s = s ;
		samples[n_samples].v = v ;
		sum_h_by_row[row] += h ;
		n_h_by_row[row]++ ;

		if(fix_sky_hue && (x < min_sky_hue_mask || x > max_sky_hue_mask)) {
		    float wgt = (1.-s) ;
		    sum_hue += h*wgt ;
		    sum_hue_wgts += wgt ;
		    sum_sat += s ;
		    sum_sat_wgts += 1.f ;
		    sum_val += v*wgt ;
		    sum_val_wgts += wgt ;
		} else {
		    float wgt = (1.-s) ;
		    sum_hue += h*wgt ;
		    sum_hue_wgts += wgt ;
		    sum_sat += s ;
		    sum_sat_wgts += 1.f ;
		    sum_val += v*wgt ;
		    sum_val_wgts += wgt ;
		}

		sum_h += h ;
		sum_s += s ;
		sum_v += v ;
		sum_n += 1.f ;

		n_samples++ ;

	    }

	    row++ ; // end of y for loop now
	}

	//fprintf(stderr, "column %d, h=%f, s=%f, v=%f\n", x, sum_h/sum_n, sum_s/sum_n, sum_v/sum_n) ;
    }

    fprintf(stderr, "\n") ;
    for(int i = 0 ; i < n_max_rows ; i++) {
	if(n_h_by_row[i] > 0) {
	    float mean = sum_h_by_row[i]/(float)n_h_by_row[i] ;
	    fprintf(stderr, "INFO:  Mean sky hue at y=%d is %f, n=%d\n", i*dy, mean, n_h_by_row[i]) ;
	}
    }
    fprintf(stderr, "\n") ;

    hue_sky = sum_hue/sum_hue_wgts ;
    sat_sky = sum_sat/sum_sat_wgts ;
    val_sky = sum_val/sum_val_wgts ;

    S_from_V_coefs[0] = sat_sky ;
    S_from_V_coefs[1] = 0. ;
    S_from_V_coefs[2] = 0. ;

    if(sat_prediction_method == 1) {
	// predict S and H as functions of V, given sample points

	struct mstat m_fit = init_reg(2) ;

	for(int i = 0 ; i < n_samples ; i++) {
	    double xv[2] = {samples[i].v,samples[i].v*samples[i].v} ;
	    double y = samples[i].s ;
	    sum_reg(&m_fit, xv, y) ;
	}

	double coefs[3] ;
	estimate_reg(&m_fit, coefs) ;
	S_from_V_coefs[0] = coefs[0] ;
	S_from_V_coefs[1] = coefs[1] ;
	S_from_V_coefs[2] = coefs[2] ;
    }

    fprintf(stderr, "Sat from v coefs %f %f %f\n",
    S_from_V_coefs[0],
    S_from_V_coefs[1],
    S_from_V_coefs[2]) ;

    // initially just set hue to a constant
    H_from_V_coefs[0] = hue_sky ;
    H_from_V_coefs[1] = 0.f ;
    H_from_V_coefs[2] = 0.f ;
    H_from_V_coefs[3] = 2.f ;

    if(0) {
	// fit H as a function of py
	// THIS IS NOT WORKING WELL AT ALL! 

	struct mstat m_fit = init_reg(1) ;
	double coefs[3] ;
	coefs[0] = coefs[1] = coefs[2] = 0. ;
	float max_py=0. ;

	for(int i = 0 ; i < n_samples ; i++) {
	    int ay=(1.-samples[i].py)*IMAGE_HEIGHT ;
	    float ty=1.-(float)(ay)/(float)max_end_of_sky ;
	    double t = log(ty*10.) ;
	    t = 1./(ty/10.+.0001) ;
	    double xv[2] = {t,t*t} ;
	    double y = samples[i].h ;
	    sum_reg(&m_fit, xv, y) ;
	}

	estimate_reg(&m_fit, coefs) ;
	H_from_V_coefs[0] = coefs[0] ;
	H_from_V_coefs[1] = coefs[1] ;
	H_from_V_coefs[2] = coefs[2] ;
	H_from_V_coefs[3] = max_py ;

	fprintf(stderr, "\n") ;
	for(int i = 0 ; i < n_max_rows ; i++) {
	    float mean=0. ;
	    if(n_h_by_row[i] > 0) {
		mean = sum_h_by_row[i]/(float)n_h_by_row[i] ;
	    }

	    float y = i*dy ;
	    if(y > max_end_of_sky) break ;

	    float ty=1.-(float)(y)/(float)max_end_of_sky ;
	    double t = log(ty*10.) ;
	    t = 1./(ty/10.+.0001) ;
	    float h_hat = H_from_V_coefs[0] + H_from_V_coefs[1]*t + H_from_V_coefs[2]*t*t ;
	    fprintf(stderr, "INFO:  Mean,pred sky hue at y,ty,t=%d,%4.2f,%8f is (%6.1f,%6.1f), n=%d\n", i*dy,ty,t, mean, h_hat, n_h_by_row[i]) ;
	}
	fprintf(stderr, "\n") ;
    }

    fprintf(stderr, "sample points, sky hue=%f\n", hue_sky) ;

    return n_samples ;
}

// convert point to normalized vector from origin (0.,0.,0.) to point
struct V3D P3toV3(float x, float y, float z)
{
    struct V3D v ;
    float sum = sqrt(x*x + y*y + z*z) ;
    v.A = x/sum ;
    v.B = y/sum ;
    v.C = z/sum ;

    return v ;
}

void set_FOV_factor()
{
    float FOV_horz_rad = M_PI/180.0*FOV_horizontal ;
    float FOV_vert_rad = FOV_horz_rad * (float)IMAGE_HEIGHT/(float)IMAGE_WIDTH ;

    proportion_to_radian_factor_x = FOV_horz_rad ;
    proportion_to_radian_factor_y = FOV_vert_rad ;
}

void set_xy_constraints()
{
    // sun_x, MAX,MIN in case a tighter constraint already
    opt_parms[0].lo = MAX(-180.f, opt_parms[0].lo) ;
    opt_parms[0].lo = MIN(180.f, opt_parms[0].lo) ;

    opt_parms[0].hi = MIN(180.f, opt_parms[0].hi) ;
    opt_parms[0].hi = MAX(-180.f, opt_parms[0].hi) ;

    // sun_py
/*      opt_parms[1].lo = MAX(min_y, opt_parms[1].lo) ;  */
/*      opt_parms[1].lo = MIN(max_y, opt_parms[1].lo) ;  */

/*      opt_parms[1].hi = MIN(max_y, opt_parms[1].hi) ;  */
/*      opt_parms[1].hi = MAX(min_y, opt_parms[1].hi) ;  */

    // horizon_py
/*      opt_parms[2].lo = MAX(min_y, opt_parms[2].lo) ;  */
/*      opt_parms[2].lo = MIN(max_y, opt_parms[2].lo) ;  */

/*      opt_parms[2].hi = MIN(max_y, opt_parms[2].hi) ;  */
/*      opt_parms[2].hi = MAX(min_y, opt_parms[2].hi) ;  */
}

float angle_between(struct V3D *v1, struct V3D *v2, float *dotproduct)
{
    *dotproduct = v1->A*v2->A + v1->B*v2->B + v1->C*v2->C ;

    // normally divide by magnitude of each vector
    // but these vectors have been already normalized so their magnitude is 1.0
    // dotproduct /= mag(v1)*mag(v2) ;
    if(*dotproduct > 1.f) *dotproduct=1.f ; // catch rounding error, acos(> 1.) is NAN
    if(*dotproduct < -1.f) *dotproduct=-1.f ; // catch rounding error, acos(> 1.) is NAN

    float angle = acos(*dotproduct) ;
    return angle ;
}

float py_adjusted_for_horizon_curvature(float px, float py)
{

#define OLD_HC
#ifdef OLD_HC
    float tmp = 1. - fabs(horizon_curvature) ;
    float f = 1. - tmp*tmp ;
    // *4 -> px is only -0.5 to 0.5, we want curvature for py as if px runs -1 to 1
    // but the user (or algorithm) considers the horizon curvature at px=0.5 
    tmp = 1. - f*px*px*4. ; 
    if(tmp < 0.) {
/*  	fprintf(stderr, "FATAL: pxpy_toV3D, negative adjusted px\n") ;  */
/*  	fprintf(stderr, "FATAL: px:%f\n", px) ;  */
/*  	fprintf(stderr, "FATAL: py:%f\n", py) ;  */
/*  	fprintf(stderr, "FATAL: horizon_curvature:%f\n", horizon_curvature) ;  */
/*  	fprintf(stderr, "FATAL: f:%f\n", f) ;  */
/*  	fprintf(stderr, "FATAL: tmp:%f\n", tmp) ;  */
	tmp = 0. ;
    }

    float hc_py = 1. - sqrt(tmp) ;
    if(horizon_curvature < 0.) hc_py = -hc_py ;

    // hc_py is where we expect the horizon to actually be at, in image coordinates
    // so if we find the distance between the image given py and hc_py, this will
    // correspond to a 'py' with curvature 0.  (with some obvious issues, but should
    // be good enough for very low curvatures

    return py-hc_py ;
#else

    // bug in this somewhere, causing NAN's in F_CIE2003()

    if(fabs(horizon_curvature) < .00001) {
	return py ;
    } else {
	float R = (horizon_curvature*horizon_curvature+.25)/(2.*horizon_curvature) ;
	R = fabs(R) ;
	float hc_py = R - sqrt(R*R-px*px) ;

	if(horizon_curvature > 0.) {
	    return py-hc_py ;
	} else {
	    return py+hc_py ;
	}
    }
#endif
}

struct V3D pxpy_toV3D(float px, float py) {
    float py_above_horizon=py_adjusted_for_horizon_curvature(px,py)  ;

    if(0) {
	// for debugging
	horizon_curvature = .1 ;
	py=.1 ;
	for(px = -0.5 ; px < 0.51 ; px += .05) {
	    py_above_horizon=py_adjusted_for_horizon_curvature(px,py)  ;
	    fprintf(stderr, "px:%5.2f py:%5.2f py_above:%5.2f\n", px, py, py_above_horizon) ;
	}
	exit(1) ;
    }

    float angle_about_Y = -px*proportion_to_radian_factor_x ;
    float angle_about_X = py_above_horizon*proportion_to_radian_factor_y ;


    float X = 0. ;
    float Y = 0. ;
    float Z = -1. ;

    float C = cos(angle_about_X) ;
    float S = sin(angle_about_X) ;
    float X1 = X ;
    float Y1 = Y*C-Z*S ;
    float Z1 = Y*S+Z*C ;

    C = cos(angle_about_Y) ;
    S = sin(angle_about_Y) ;
    Y = Y1 ;
    X = X1*C+Z1*S ;
    Z = -X1*S+Z1*C ;

    if(verbose) {
	struct V3D view = P3toV3(X,Y,Z) ;
	float dp_sun, dp_zenith ;
	fprintf(stderr, "pxpy: %5.2f %5.2f, AOX:%7.2f AOY:%7.2f ", px, py, angle_about_X*180./M_PI, angle_about_Y*180./M_PI) ;

	float Gamma = angle_between(&view,&V_sun,&dp_sun) ;
	float Theta = angle_between(&view,&V_zenith,&dp_zenith) ;
	fprintf(stderr, "Asun:%7.2f Azenith:%7.2f ", 180./M_PI*Gamma, 180./M_PI*Theta) ;

	return view ;
    }

    return P3toV3(X,Y,Z) ;
}

// give a image relative pixel (where the bottom of the image is py=0, top is py=1, left is -p_half_image_width, right is p_half_image_width)
// create a 3d vector
struct V3D image_relative_pixel_to_V3(float px, float py)
{
    py -= horizon_py ;
    if(py < 0.) py = 0. ;

    return pxpy_toV3D(px,py) ;
}

struct V3D image_relative_pixel_to_V3_noclip(float px, float py)
{
    py -= horizon_py ;

    return pxpy_toV3D(px,py) ;
}

float F_CIE2003(float px, float py, float *pGamma, float *pTheta)
{
    /* assume observer is at origin */
    struct V3D view= image_relative_pixel_to_V3(px, py) ;
    float dp_zenith, dp_sun ; // dot product of vectors to view

    *pGamma = angle_between(&view,&V_sun,&dp_sun) ;

    float cosGamma = cos(*pGamma) ;

    // check if we are within the sun diameter
    if(*pGamma < M_PI/180.f*0.533f/2.) {
	*pGamma = 0.f ;
	cosGamma = 1.f ;
    }

    *pTheta = angle_between(&view,&V_zenith,&dp_zenith) ;

    if(0) {
	// JJW reduce angle between zenith and horizon,
	// which will raise the level of horizon glow (via A & B coefs and Perez formula)
	float horizon_angle_factor = (M_PI/2.)/(M_PI/2.-horizon_lift_angle_radians) ;

	*pTheta *= horizon_angle_factor ;
    }

    if(*pTheta > M_PI/2.) *pTheta = M_PI/2. ;

    // clamp angle between view and zenith at PI/2 so prediction
    // values don't fall below horizon which is not good at all
    if(*pTheta > M_PI/2.-0.01) *pTheta=M_PI/2.-.01 ;

    if(isnan(*pGamma)) {
	fprintf(stderr, "NAN in F_CIE2003, sun_x=%f sun_py=%f px=%f py=%f horizon_py=%f theta=%f gamma=%f\n",
	    sun_x, sun_py, px, py, horizon_py, *pTheta, *pGamma) ;
	fprintf(stderr, "View: %f %f %f\n", view.A, view.B, view.C) ;
	fprintf(stderr, "Sun: %f %f %f\n", V_sun.A, V_sun.B, V_sun.C) ;
	float dotproduct = view.A*V_sun.A + view.B*V_sun.B + view.C*V_sun.C ;
	fprintf(stderr, "dotproduct %f\n", dotproduct) ;
	fprintf(stderr, "acos %f\n", acos(dotproduct)) ;
	exit(1) ;
    }

    float L = (1.+perez_A*exp(perez_B/cos(*pTheta)))*(1.+perez_C*(exp(perez_D*(*pGamma))-exp(perez_D*M_PI/2.)) + perez_E*cosGamma*cosGamma) ;

    if(isnan(L)) {
	fprintf(stderr, "NAN in F_CIE2003, sun_x=%f sun_py=%f px=%f py=%f horizon_py=%f theta=%f gamma=%f\n",
	    sun_x, sun_py, px, py, horizon_py, *pTheta, *pGamma) ;
	fprintf(stderr, "View: %f %f %f\n", view.A, view.B, view.C) ;
	fprintf(stderr, "Zenith: %f %f %f\n", V_zenith.A, V_zenith.B, V_zenith.C) ;
	float dotproduct = view.A*V_zenith.A + view.B*V_zenith.B + view.C*V_zenith.C ;
	fprintf(stderr, "dotproduct %f\n", dotproduct) ;
	fprintf(stderr, "acos %f\n", acos(dotproduct)) ;
	exit(1) ;
    }

    if(verbose) {
	fprintf(stderr, "sky L:%5.2f\n", L) ;
    }

    return L ;
}

float find_maximum_vhat(void)
{
    float gamma,theta ;
    float px, py ;

    float sun_px = sun_x_angle2px(sun_x) ;
    float max_vhat = 0. ;

    if(sun_px >= -0.5 && sun_px <= 0.5 && sun_py >= 0.0 && sun_py <= 1.0) {
	// if the sun is in the image, get vhat directly at the sun
	max_vhat = F_CIE2003(sun_x_angle2px(sun_x), sun_py, &gamma, &theta) ;
    }


    for(px = -0.5f ; px < 0.501f ; px += .05f) {
	for(py = 0.f ; py < 1.001f ; py += .05f) {
	    float vhat_here = F_CIE2003(px, py, &gamma, &theta) ;
	    if(vhat_here > max_vhat) max_vhat = vhat_here ;
	}
    }

    if(max_vhat < 1.0f) max_vhat = 1.0f ;

    return max_vhat ;
}

float sample_based_v_correction = 1.f ; // factor to make mean v from predicted "hsv" match mean v from actual

void predict_sky_huesat_from_val(float vhat, float *pH, float *pS, float px, float py)
{

//    Regression 'Sat=f(Value)' from a fit based on a sky image
//****************************************************************
//Name           	       Theta	      StdErr	 T-stat
//[0='B0']       	      5.7540	      0.0206	 279.70
//[1='B1']       	    -11.8315	      0.0589	-200.80
//[2='B2']       	      6.4355	      0.0419	 153.54

/*      shat = 5.7540-11.8315*vhat+6.4355*vhat*vhat ;  */

    float shat = S_from_V_coefs[0] + S_from_V_coefs[1]*vhat + S_from_V_coefs[2]*vhat*vhat ;


    if(shat > 1.0) shat = 1.0 ;
    if(shat < 0.0) shat = 0.0 ;

    shat *= sat_depth ;
    vhat *= sun_lum ;

    if(vhat > 1.f) vhat = 1.f ;
    if(vhat < 0.f) vhat = 0.f ;

    double t ;
    if(py > H_from_V_coefs[3]) {
	t = log(H_from_V_coefs[3]) ;
    } else {
	t = log(py) ;
    }

    int ay=(1.-py)*IMAGE_HEIGHT ;
    float ty=1.-(float)(ay)/(float)max_end_of_sky ;
    t = 1./(ty/10.+.0001) ;
    float hhat = H_from_V_coefs[0] + H_from_V_coefs[1]*t + H_from_V_coefs[2]*t*t ;

    *pH = hhat ;
    *pS = shat ;
}

void predict_sky_hsv(float px, float py, float *pH, float *pS, float *pV)
{
    float vhat ;

    float gamma ; // angle between view and sun vector
    float theta ; // angle between view and zenith vector

    vhat = F_CIE2003(px, py, &gamma, &theta)/maximum_vhat ;

    if(vhat > 1.f) vhat = 1.f ;
    if(vhat < 0.f) vhat = 0.f ;

    if(vhat < 0.) {
	fprintf(stderr, "WARNING: F_CIE2003 returns negative vhat, fixing!\n") ;
	vhat = 0. ;
    }

    float cos_gamma = fabs(cos(gamma)) ; // note, we can get this direct from F_CIE2003 if we want

    if(sky_gamma > 1.e-10)
	vhat *= (1.-pow(cos_gamma,1./sky_gamma)) ;

    //vhat = pow(vhat, 1./sky_gamma) ;

    vhat *= sample_based_v_correction ; // factor to make mean v from predicted "hsv" match mean v from actual, was computed during fits

    if(!model_is_being_fit)
	vhat = pow(vhat, 1./final_sky_gamma) ;

    predict_sky_huesat_from_val(vhat, pH, pS, px, py) ;

    *pV = vhat * exposure_factor ;
}

void predict_sky_color(float px, float py, uint16 *rhat, uint16 *ghat, uint16*bhat)
{
    float hhat, shat, vhat ;

    predict_sky_hsv(px,py,&hhat,&shat,&vhat) ;

    hsv2rgb16(hhat,shat,vhat,rhat,ghat,bhat) ;
}


#define V_from_HSV 0x01
#define FULL_RGB 0x02

int optimize_type = V_from_HSV ;
int calln=0 ;
int autoset_lumf_flag=1 ;

float samples_error(void)
{
    int i ;
    float sum_err_sq = 0. ;
    //optimize_type = FULL_RGB ;

    set_FOV_factor() ;


    V_sun= image_relative_pixel_to_V3_noclip(sun_x_angle2px(sun_x), sun_py) ;

    maximum_vhat = find_maximum_vhat() ;

    if(isinf(maximum_vhat)) {
	fprintf(stderr, "maximum_vhat is NAN in grid optimize: sun_x:%f sun_py%f\n", sun_x, sun_py) ;
	exit(1) ;
    }

    float lum_f=1. ;

    float sum_v=0 ;
    float sum_v_hat=0 ;

    // recalculate value correction using sky samples and
    // current prediction coefficients
    sample_based_v_correction = 1.f ;

    for(i = 0 ; i < n_samples_to_optimize ; i++) {
	float h,s,v ;
	predict_sky_hsv(samples[i].px, samples[i].py, &h, &s, &v) ;
	samples[i].v_hat = v ;
	sum_v += samples[i].v ;
	sum_v_hat += samples[i].v_hat ;
    }

    sample_based_v_correction = sum_v/sum_v_hat ;

    // XXXX
    // clamp this -- gets way to high if the sun is in the image and the CIE parms allow for a sun
/*      if(sample_based_v_correction > 4.0) sample_based_v_correction = 4.0 ;   */

    if(optimize_type == V_from_HSV) {

	float sum_vhat =0. ;

	for(i = 0 ; i < n_samples_to_optimize ; i++) {
	    float v_hat = samples[i].v_hat*lum_f ;
	    float err_v = (samples[i].v - v_hat) ;

	    if(isinf(err_v)) {
		fprintf(stderr, "err_v is inf in grid optimize: v:%f px:%f py%f\n", v_hat, samples[i].px, samples[i].py) ;
		exit(1) ;
	    }

	    sum_err_sq += err_v*err_v ;
	    sum_vhat += v_hat ;
	}

	if(!(calln % 1000)) fprintf(stderr, "mean vhat=%f\n", sum_vhat/(float)n_samples_to_optimize) ;
    } else {
	for(i = 0 ; i < n_samples_to_optimize ; i++) {
	    uint16 rhat,ghat,bhat ;
	    predict_sky_color(samples[i].px, samples[i].py, &rhat, &ghat, &bhat) ;
	    float err_r = (samples[i].r - rhat) ;
	    float err_g = (samples[i].g - ghat) ;
	    float err_b = (samples[i].b - bhat) ;
	    sum_err_sq += err_r*err_r + err_g*err_g + err_b*err_b ;
	}
    }

    calln++ ;


    return sum_err_sq ;
}

void set_cie_parms(int CIE_index)
{
    perez_A = cie_parms[CIE_index].A ;
    perez_B = cie_parms[CIE_index].B ;
    perez_C = cie_parms[CIE_index].C ;
    perez_D = cie_parms[CIE_index].D ;
    perez_E = cie_parms[CIE_index].E ;
}

float optimize_grid_function(int CIE_index)
{

    set_cie_parms(CIE_index) ;

    return samples_error() ;

}

float optimize_horizon_function(int CIE_index, float B)
{

    perez_A = cie_parms[CIE_index].A ;
    perez_B = B ;
    perez_C = cie_parms[CIE_index].C ;
    perez_D = cie_parms[CIE_index].D ;
    perez_E = cie_parms[CIE_index].E ;

    return samples_error() ;

}

float find_horizon_py(int n_samples)
{
    int i ;
    float best_feval = 1.e30 ;
    float best_hpy=lowest_sky_py_found ;
    float best_B=-.1 ;
    float best_hc=0. ;
    float hpy0 = lowest_sky_py_found ;
    if(hpy0 > 0.5) hpy0 = 0.5 ;
    float hpy1 = lowest_sky_py_found ;
    if(hpy1 > 0.99) hpy0 = 0.99 ;

    optimize_type = V_from_HSV ;

    autoset_lumf_flag=1 ;
    for(i=0 ; i < MAX_PARMS ; i++) {
	if(opt_parms[i].grid_optimize_flag == 0) {
	    if(!strcmp(opt_parms[i].name, "sun_lum")) {
		autoset_lumf_flag=0 ;
	    }
	}
    }

    n_samples_to_optimize=n_samples ;

    for(horizon_curvature=-.1 ; horizon_curvature < .101 ; horizon_curvature += .02) {
	for(perez_B=-.1 ; perez_B > -1.01 ; perez_B -= .1) {
	    for(horizon_py = hpy0 ; horizon_py < hpy1+.0001 ; horizon_py += .05) {
		float feval = optimize_horizon_function(15,perez_B) ;
		//float feval = optimize_grid_function(8) ;
		if(feval < best_feval) {
		    best_B = perez_B ;
		    best_hpy = horizon_py ;
		    best_hc = horizon_curvature ;
		    best_feval = feval ;
		}
	    }
	}
	printf("============ horizon curvature:%5.2lf, horizon:%f, best_feval %f ============\n", horizon_curvature, best_hpy, best_feval) ;
    }

    horizon_curvature = best_hc ;
    horizon_py = best_hpy ;
    perez_B = best_B ;
    //horizon_py = lowest_sky_py_found ;
    //best_hpy = lowest_sky_py_found ;


    printf("============ Grid, lowest_horizon is %5.2lf, Find horizon py returns ============\n", lowest_sky_py_found) ;
    dump_parameters_to_stdout() ;

    return best_hpy ;
}

float optimize_grid(int n_samples, int verbose)
{
    int i ;
    float best_feval = 1.e30 ;
    int best_i = -1 ;
    float feval ;
    optimize_type = V_from_HSV ;

    autoset_lumf_flag=1 ;
    for(i=0 ; i < MAX_PARMS ; i++) {
	if(opt_parms[i].grid_optimize_flag == 0) {
	    if(!strcmp(opt_parms[i].name, "sun_lum")) {
		autoset_lumf_flag=0 ;
	    }
	}
    }

    n_samples_to_optimize=n_samples ;

    if(CIE_sky_index > 0) {
	i=CIE_sky_index-1 ;

	feval = optimize_grid_function(i) ;
	if(feval < best_feval) {
	    best_i = i ;
	    best_feval = feval ;
	}
    } else {
	for(i = 0 ; i < N_CIE_STD_PARMS ; i++) {
	    int is_uniform_sky=0 ;
	    if(cie_parms[i].C < .01 && cie_parms[i].D < .01) is_uniform_sky=1;

	    if(is_uniform_sky && allowed_sky_type == NONUNIFORM_CIE_SKIES) continue ;
	    if(!is_uniform_sky && allowed_sky_type == UNIFORM_CIE_SKIES) continue ;

	    feval = optimize_grid_function(i) ;
	    if(feval < best_feval) {
		best_i = i ;
		best_feval = feval ;
	    }
	}
	feval = optimize_grid_function(best_i) ;
    }


    if(verbose) printf("GRID feval:%g CIE index %d, %s\n", feval, cie_parms[best_i].type, cie_parms[best_i].description) ; ;

    if(verbose) dump_parameters_to_stdout() ;

    return feval ;
}

float smart_optimize_function(float *pParams)
{
    int i ;

    int j=1 ;

    for(i=0 ; i < MAX_PARMS ; i++) {
	if(opt_parms[i].optimize_flag == 1) {
	    *opt_parms[i].variable_address = pParams[j] ;
	    j++ ;
	}
    }

    return samples_error() ;
}

void smart_optimize(int n_samples)
{
    float guess[MAX_PARMS+1] = {0,sun_x,sun_py,sun_lum,horizon_py,perez_A,perez_B,perez_C,perez_D,perez_E,FOV_horizontal} ;
    float guess_delta[MAX_PARMS+1], c_lo[MAX_PARMS+1], c_hi[MAX_PARMS+1] ;
    int i ;
    int nparms_to_optimize=0 ;

    n_samples_to_optimize=n_samples ;

    int j=1 ;

    for(i=0 ; i < MAX_PARMS ; i++) {
	if(opt_parms[i].optimize_flag == 1) {
	    guess[j] = *opt_parms[i].variable_address ;
	    c_lo[j] = opt_parms[i].lo ;
	    c_hi[j] = opt_parms[i].hi ;
	    if(fabs(guess[j]) > 1.e-5) {
		guess_delta[j] = fabs(guess[j])/10. ;
	    } else {
		guess_delta[j] = fabs(opt_parms[i].hi-opt_parms[i].lo)/10. ;
	    }

	    if(!strcmp(opt_parms[i].name, "sun_x")) {
		// special case for sun_x_angle
		guess_delta[j] = 5. ;
		c_lo[j] = guess[j]-20. ;
		c_hi[j] = guess[j]+20. ;
	    }

	    if(!strcmp(opt_parms[i].name, "horizon_py")) {
		// special case for horizon_py
		guess_delta[j] = .01 ;
		c_lo[j] = 0. ;
		c_hi[j] = ((float)(IMAGE_HEIGHT)-(float)max_end_of_sky)/(float)IMAGE_HEIGHT ;
	    }

	    j++ ;
	    printf("OPTIMIZE->%s\n", opt_parms[i].name) ;
	    nparms_to_optimize++ ;
	}
    }

    float feval = smart_optimize_function(guess) ;

    printf("Initial guess feval:%g\n", feval) ;

    dump_parameters_to_stdout() ;

    feval = AmoebaFit(smart_optimize_function,1.e-7,1000,nparms_to_optimize,guess,guess_delta,c_lo,c_hi,0) ;

    feval = smart_optimize_function(guess) ;

    printf("Optimized feval:%g\n", feval) ;

    dump_parameters_to_stdout() ;

}



// how much of the new estimate to use at pixel x, y
float raw_compute_feather_at_y(int x, int y, int feather_length, uint16 b, float feather_factor)
{
    //if(x == 0) {
	//fprintf(stderr, "INFO: CFAY 1, x:%d y:%d sos:%d eos:%d fl:%d ff:%f\n", x, y, raw_start_of_sky[x], end_of_sky[x], feather_length, feather_factor) ;
    //} 

    // compute feathering amount
    float fp0 = 1. ; // proportion of prediction to use

    if(y <= raw_start_of_sky[x]) return 1. ;
    if(y > raw_start_of_sky[x]+feather_length) return 0. ;

    if(y >= raw_start_of_sky[x])
	fp0 = 1.-((float)y-(float)raw_start_of_sky[x])/(float)feather_length ;

    if(b < floatBLK*MAX16) {
	// pixel is black, use 100% of predicted value
	fp0 = 1.0 ;
    }

    float fp1 = 1.-fp0 ;
    fp1 *= feather_factor ; // adjust amount of feathering by feather_factor (useful for debugging)
    fp0 = 1.-fp1 ;

    if(fp0 < 0.0) {
	fprintf(stderr, "FATAL: raw CFAY 1, fp0:%f x:%d y:%d sos:%d fl:%d ff:%f\n", fp0, x, y, raw_start_of_sky[x], feather_length, feather_factor) ;
	exit(1) ;
    }

    if(y == end_of_sky[x]) {
	// make the end of sky pixel half estimated
	fp0 = .5 ;
	fp1 = .5 ;
    }

    if(fp0 < 0.0) {
	fprintf(stderr, "FATAL: raw CFAY 2, fp0:%f x:%d y:%d sos:%d fl:%d ff:%f\n", fp0, x, y, raw_start_of_sky[x], feather_length, feather_factor) ;
	exit(1) ;
    }

    fp0 = powf(fp0,nonlinear_feather) ;
    return fp0 ;
}

int raw_compute_feather_length(int x, int *pFeather_end_y, float depth_of_fill, float feather_factor, int extra)
{
    float sky_height ;

    if(depth_of_fill_is_absolute) {
	sky_height = depth_of_fill_absolute_y - raw_start_of_sky[x] ;
    } else {
	sky_height = end_of_sky[x] - start_of_sky[x] ;
    }

    *pFeather_end_y = raw_start_of_sky[x] + (int)(sky_height*depth_of_fill+0.5) + extra ;

    if(*pFeather_end_y > end_of_sky[x]) *pFeather_end_y = end_of_sky[x] ;

    if(estimate_only) {
	*pFeather_end_y = IMAGE_HEIGHT-1 ;
	feather_factor=0. ;
    }

    return *pFeather_end_y - raw_start_of_sky[x] ;
}

// how much of the new estimate to use at pixel x, y
float compute_feather_at_y(int x, int y, int feather_length, uint16 b, float feather_factor)
{
    //if(x == 0) {
	//fprintf(stderr, "INFO: CFAY 1, x:%d y:%d sos:%d eos:%d fl:%d ff:%f\n", x, y, start_of_sky[x], end_of_sky[x], feather_length, feather_factor) ;
    //} 

    // compute feathering amount
    float fp0 = 1. ; // proportion of prediction to use

    if(y >= start_of_sky[x])
	fp0 = 1.-((float)y-(float)start_of_sky[x])/(float)feather_length ;

    if(b < floatBLK*MAX16) {
	// pixel is black, use 100% of predicted value
	fp0 = 1.0 ;
    }

    float fp1 = 1.-fp0 ;
    fp1 *= feather_factor ; // adjust amount of feathering by feather_factor (useful for debugging)
    fp0 = 1.-fp1 ;

    if(fp0 < 0.0) {
	fprintf(stderr, "FATAL: CFAY 1, fp0:%f x:%d y:%d sos:%d fl:%d ff:%f\n", fp0, x, y, start_of_sky[x], feather_length, feather_factor) ;
	exit(1) ;
    }

    if(y == end_of_sky[x]) {
	// make the end of sky pixel half estimated
	fp0 = .5 ;
	fp1 = .5 ;
    }

    if(fp0 < 0.0) {
	fprintf(stderr, "FATAL: CFAY 2, fp0:%f x:%d y:%d sos:%d fl:%d ff:%f\n", fp0, x, y, start_of_sky[x], feather_length, feather_factor) ;
	exit(1) ;
    }

    fp0 = powf(fp0,nonlinear_feather) ;
    return fp0 ;
}

int compute_feather_length(int x, int *pFeather_end_y, float depth_of_fill, float feather_factor, int extra)
{
    float sky_height ;

    if(depth_of_fill_is_absolute) {
	sky_height = depth_of_fill_absolute_y - start_of_sky[x] ;
    } else {
	sky_height = end_of_sky[x] - start_of_sky[x] ;
    }

    *pFeather_end_y = start_of_sky[x] + (int)(sky_height*depth_of_fill+0.5) + extra ;

    if(*pFeather_end_y > end_of_sky[x]) *pFeather_end_y = end_of_sky[x] ;

    if(estimate_only) {
	*pFeather_end_y = IMAGE_HEIGHT-1 ;
	feather_factor=0. ;
    }

    return *pFeather_end_y - start_of_sky[x] ;
}

int wrong_hue_or_black(tdata_t *image, int32 x, int32 y)
{
    uint16 r,g,b,a ;
    tif_get4c(image,x,y,r,g,b,a) ;

    /* if alpha channel is nonzero then hugin has placed image data here */
    if(a < HALF16) return 0 ;

    float h,s,v ;
    rgb2hsv16(r,g,b,&h,&s,&v) ;
    float e_h = fabs((h-hue_sky)/hue_sky) ;
    float e_s = fabs((s-sat_sky)/sat_sky) ;
    float e_v = fabs((v-val_sky)/val_sky) ;

    if(e_h > .05) return 1 ;
    if(e_s > .20) return 1 ;
    if(e_v > .20) return 1 ;

    return 0 ;
}


int xy_is_transparent(tdata_t *image, int32 x, int32 y)
{
    /* if alpha channel is nonzero then hugin has placed image data here */
    if(IMAGE_HAS_ALPHA) {
	if(((uint16 *)(image[y]))[IMAGE_NSAMPLES*x+3] < 256*64)
	    return 1 ;
    }

    return 0 ;
}


int xy_is_black_pixel(tdata_t *image, int32 x, int32 y)
{
    /* if alpha channel is nonzero then hugin has placed image data here */
    if(IMAGE_HAS_ALPHA) {
	if(((uint16 *)(image[y]))[IMAGE_NSAMPLES*x+3] < 256*64)
	    return 0 ;
    }

    if(((uint16 *)(image[y]))[IMAGE_NSAMPLES*x+0] > 256*15)
	return 0 ;
    if(((uint16 *)(image[y]))[IMAGE_NSAMPLES*x+1] > 256*15)
	return 0 ;
    if(((uint16 *)(image[y]))[IMAGE_NSAMPLES*x+2] > 256*15)
	return 0 ;

    return 1 ;
}

void repair_alpha(tdata_t *image)
{
    // hugin is creating pixels near the edges of the data with low r,g,b values
    // presumeably because they didn't represent a full pixel from an original
    // image but are part of a pixel.  These often also have really bad colors,
    // the max 16 bit alpha is correct at 65535
    // out of hugin, so we need to roll through the entire image and set the alpha
    // and r,g,b of those suspect pixels to 0,0,0,0
    int x,y ;

    for(x = 0 ; x < IMAGE_WIDTH ; x++) {
	for(y = 0 ; y < IMAGE_HEIGHT ; y++) {
	    uint16 r,g,b ;
	    tif_get3c(image,x,y,r,g,b) ;

	    // nearest neighbor search +/- 2 pixels to find max r,g,b in area
	    uint max_r=0 ;
	    uint max_g=0 ;
	    uint max_b=0 ;
	    int x0 = x-2 ; if(x0 < 0) x0 = 0 ;
	    int y0 = y-2 ; if(y0 < 0) y0 = 0 ;
	    int x1 = x+2 ; if(x1 > IMAGE_WIDTH-1) x1 = IMAGE_WIDTH-1 ;
	    int y1 = y+2 ; if(y1 > IMAGE_HEIGHT-1) y1 = IMAGE_HEIGHT-1 ;
	    int found_edge=0;  // are we near the edge of image data (i.e. there is a transparent pixel in the neighborhood?) ;

	    for(int dx=x0 ; dx<=x1 ; dx++) {
		for(int dy=y0 ; dy<=y1 ; dy++) {
		    uint16 r_neighbor,g_neighbor,b_neighbor,a_neighbor ;
		    tif_get4c(image,dx,dy,r_neighbor,g_neighbor,b_neighbor,a_neighbor) ;

		    if(a_neighbor > 0) {
			if(max_r < r_neighbor) max_r = r_neighbor ;
			if(max_g < g_neighbor) max_g = g_neighbor ;
			if(max_b < b_neighbor) max_b = b_neighbor ;
		    } else {
			found_edge=1 ;
		    }
		}
	    }

	    if(found_edge && (r < max_r/8 || g < max_g/8 || b < max_b/8)) {
		// this is a suspect pixel on the edge, with at least one very low color component
		tif_set4c(image,x,y,0,0,0,0) ;
	    }
	}
    }

}

void set_minmax_sky_values(void)
{
    max_end_of_sky=0 ;
    min_end_of_sky=IMAGE_HEIGHT-1 ;
    max_start_of_sky=0 ;
    min_start_of_sky=IMAGE_HEIGHT-1 ;

    for(int x = 0 ; x < IMAGE_WIDTH ; x++) {
	if(column_mask[x] == 1) continue ;
	if(max_end_of_sky < end_of_sky[x]) max_end_of_sky = end_of_sky[x] ;
	if(max_start_of_sky < start_of_sky[x]) max_start_of_sky = start_of_sky[x] ;
	if(min_end_of_sky > end_of_sky[x]) min_end_of_sky = end_of_sky[x] ;
	if(min_start_of_sky > start_of_sky[x]) min_start_of_sky = start_of_sky[x] ;
    }
}

#ifdef GMR_DEBUG
int gmr_debug = 0 ;
#endif

int get_mean_rgb(tdata_t *image, int xc,int yc,double rcoef[],double gcoef[], double bcoef[],int search_width, int clip_eos_flag, int repair_bad_points_flag, int recursion_level)
{
    int x0 = xc-search_width ;
    int x1 = xc+search_width ;

    if(x0 < 0) x0 = 0 ;
    if(x1 < 0) x1 = 0 ;
    if(x0 > IMAGE_WIDTH-1) x0 = IMAGE_WIDTH-1 ;
    if(x1 > IMAGE_WIDTH-1) x1 = IMAGE_WIDTH-1 ;

    rcoef[0] =  rcoef[1] =  rcoef[2] =  rcoef[3] = 0. ;
    gcoef[0] =  gcoef[1] =  gcoef[2] =  gcoef[3] = 0. ;
    bcoef[0] =  bcoef[1] =  bcoef[2] =  bcoef[3] = 0. ;

    int x,y ;
    int increment=2 ;

#ifdef GMR_DEBUG
    if(gmr_debug) {
	fprintf(stderr, "GMR: x,y:%d,%d\n", xc,yc) ;
    }
#endif

    // first find the coordinates of valid point and store them
    int max_points = (2*search_width)*(2*search_width)/2+16 ;

    if(repair_bad_points_flag) {
	max_points *= 4 ;
	increment=1 ;
    }

    uint16 *x_possible_coords, *y_possible_coords, n_possible_points ;
    uint16 *x_coords, *y_coords, n_points ;
    uint8 *point_status ;

    n_points=0 ;
    point_status = (uint8 *)calloc(max_points, sizeof(uint8)) ;
    x_coords = (uint16 *)calloc(max_points, sizeof(uint16)) ;
    y_coords = (uint16 *)calloc(max_points, sizeof(uint16)) ;

    n_possible_points=0 ;
    x_possible_coords = (uint16 *)calloc(max_points, sizeof(uint16)) ;
    y_possible_coords = (uint16 *)calloc(max_points, sizeof(uint16)) ;

    // look through all possible pixels
    for(x = x0 ; x <= x1 ; x += increment) {
	if(column_mask[x] == 1) continue ;

	int y1 = yc+search_width ;
	if(y1 < start_of_sky[x]) continue ;

	int y0 = yc-search_width ;
	if(y0 < start_of_sky[x]) y0 = start_of_sky[x] ;

	if(clip_eos_flag && y1 > end_of_sky[x]) y1 = end_of_sky[x] ;

#ifdef GMR_DEBUG
	if(gmr_debug) {
	    fprintf(stderr, "GMR: col search x:%d, y0,y1:%d,%d\n", x,y0,y1) ;
	}
#endif

	for(y = y0 ; y <= y1 ; y += increment) {
	    /* if alpha channel is small then hugin indicates no image data here */
#ifdef GMR_DEBUG
	    if(gmr_debug) {
		uint16 r,g,b,a ;
		tif_get4c(image,x,y,r,g,b,a) ;
		fprintf(stderr, "GMR: pixel x,y:%d,%d rgba=%5d,%5d,%5d,%5d\n", x,y,r,g,b,a) ;
	    }
#endif
	    if(((uint16 *)(image[y]))[IMAGE_NSAMPLES*x+3] < MAX16)
		continue ;

/*  	    if(((uint16 *)(image[y]))[IMAGE_NSAMPLES*x+0] <= BLK16) continue ;  */
/*  	    if(((uint16 *)(image[y]))[IMAGE_NSAMPLES*x+1] <= BLK16) continue ;  */
/*  	    if(((uint16 *)(image[y]))[IMAGE_NSAMPLES*x+2] <= BLK16) continue ;  */

	    if(n_possible_points == max_points) {
		fprintf(stderr, "FATAL: n_possible_points > max_points(%d) in get_mean_rgb(), search_width=%d\n", n_possible_points, search_width) ;
		exit(1) ;
	    }

#ifdef GMR_DEBUG
	    if(gmr_debug) {
		fprintf(stderr, "GMR: ADDED pixel x,y:%d,%d\n", x,y) ;
	    }
#endif

	    x_possible_coords[n_possible_points] = x ;
	    y_possible_coords[n_possible_points] = y ;
	    point_status[n_possible_points] = 1 ;
	    n_possible_points++ ;
	}
    }

    if(n_possible_points <  16) {
/*  	fprintf(stderr, "WARNING: n_possible_points is < 16 in get_mean_rgb(), x,y=%d,%d search_width=%d\n", xc, yc, search_width) ;  */

	free(point_status) ;
	free(x_possible_coords) ;
	free(y_possible_coords) ;
	free(x_coords) ;
	free(y_coords) ;

	recursion_level-- ;

	if(recursion_level <= 0)
	    return 0 ;

	search_width *= 2 ;
/*  	fprintf(stderr, "         increasing search_width to %d\n", search_width) ;  */

	return get_mean_rgb(image, xc,yc,rcoef,gcoef, bcoef,search_width, clip_eos_flag, repair_bad_points_flag, recursion_level) ;
    }

    int c ;

    // compute mean r,g,b to help look for outliers
    float sum[3] = {0.,0.,0.} ;
    float mean[3] = {0.,0.,0.} ;

    for(int i=0 ; i < n_possible_points ; i++) {
	x = x_possible_coords[i] ;
	y = y_possible_coords[i] ;
	sum[0] += (float)((uint16 *)(image[y]))[IMAGE_NSAMPLES*x+0] ; // R
	sum[1] += (float)((uint16 *)(image[y]))[IMAGE_NSAMPLES*x+1] ; // G
	sum[2] += (float)((uint16 *)(image[y]))[IMAGE_NSAMPLES*x+2] ; // B
    }

    mean[0] = sum[0] / (float)n_possible_points ;
    mean[1] = sum[1] / (float)n_possible_points ;
    mean[2] = sum[2] / (float)n_possible_points ;

#ifdef GMR_DEBUG
    if(gmr_debug) {
	fprintf(stderr, "\nGMR: mean rgb %5d,%5d,%5d\n", (int)mean[0],(int)mean[1],(int)mean[2]) ;
    }
#endif

    int min_x=IMAGE_WIDTH-1 ;
    int max_x=0 ;
    int min_y=IMAGE_HEIGHT-1 ;
    int max_y=0 ;

    // look at differences from the mean to throw out bad colored pixels (and mark them needing repair
    for(int i=0 ; i < n_possible_points ; i++) {
	x = x_possible_coords[i] ;
	y = y_possible_coords[i] ;

	float r,g,b ;
	tif_get3c(image,x,y,r,g,b) ;

	point_status[i] = 0 ; // this point might need repair

	float r_err = fabs( (r-mean[0])/mean[0]) ;
	float g_err = fabs( (g-mean[1])/mean[1]) ;
	float b_err = fabs( (b-mean[2])/mean[2]) ;
#ifdef GMR_DEBUG
	fprintf(stderr, "GMR:      err x,y:%d,%d %f %f %f\n", x,y,r_err,g_err,b_err) ;
#endif

	// throw out point of different by more than 3% from mean ;
	if(r_err > repair_tos_thresh) continue ;
	if(g_err > repair_tos_thresh) continue ;
	if(b_err > repair_tos_thresh) continue ;

	point_status[i] = 1 ; // mark this point as not needing repair, if requested

	x_coords[n_points] = x ;
	y_coords[n_points] = y ;
	n_points++ ;

	if(min_x > x) min_x=x ;
	if(max_x < x) max_x=x ;
	if(min_y > y) min_y=y ;
	if(max_y < y) max_y=y ;
    }

#ifdef GMR_DEBUG
    if(gmr_debug) {
	fprintf(stderr, "\nGMR: n_possible,n %d,%d\n", n_possible_points, n_points) ;
	if(xc == 1705 && yc == 91)
	    exit(1) ;
    }
#endif

    if(n_points < 16) {
	free(point_status) ;
	free(x_possible_coords) ;
	free(y_possible_coords) ;
	free(x_coords) ;
	free(y_coords) ;

	recursion_level-- ;

	if(recursion_level <= 0) {
	    fprintf(stderr, "WARNING: number of fitting points is < 16 in get_mean_rgb, coefs are all 0.0!\n") ;
	    return 0 ;
	}
	fprintf(stderr, "WARNING: number of fitting points is < 16 in get_mean_rgb, doubling search width\n") ;

	search_width *= 2 ;
/*  	fprintf(stderr, "         increasing search_width to %d\n", search_width) ;  */

	return get_mean_rgb(image, xc,yc,rcoef,gcoef, bcoef,search_width, clip_eos_flag, repair_bad_points_flag, recursion_level) ;
    }

// fit types, FIT_X_Y2 is color=B0 + B1*x + B2*y + B3*y*y ;
// fit types, FIT_X_Y1 is color=B0 + B1*x + B2*y + 0 *y*y ;
// fit types, FIT_0_Y2 is color=B0 + 0 *x + B2*y + B3*y*y ;
// fit types, FIT_0_Y1 is color=B0 + 0 *x + B2*y + 0 *y*y ;
// fit types, FIT_X_0  is color=B0 + B1*x + 0 *y + 0 *y*y ;
// fit types, FIT_0_0  is color=B0 + 0 *x + 0 *y + 0 *y*y ;

// note the first 2 places are the order of the "X" coord used, places 3&4 are the order of the "Y" coord used
#define FIT_X_Y2 0x0102
#define FIT_X_Y1 0x0101
#define FIT_X_Y0 0x0101
#define FIT_0_Y2 0x0002
#define FIT_0_Y1 0x0001
#define FIT_0_Y0 0x0000

    int fit_type=0x0000 ;

    if((max_x-min_x) > 0) fit_type |= 0x0100 ;

    if((max_y-min_y) > 1)
	fit_type |= 0x0002 ;
    else if((max_y-min_y) > 0)
	fit_type |= 0x0001 ;

    int order_x = (fit_type&0xff00) >> 8 ;
    int order_y = (fit_type&0x00ff) ;

    int total_vars = order_x + order_y ;

    struct mstat fit_struct[3] ;

    for(c=0 ; c < 3 ; c++) {
	fit_struct[c] = init_reg(total_vars) ;
    }

    for(int i=0 ; i < n_points ; i++) {
	x = x_coords[i] ;
	y = y_coords[i] ;

	double xv[3] ;
	int j=0 ;
	if(order_x) { xv[j] = x ; j++; }
	if(order_y) { xv[j] = y ; j++; }
	if(order_y > 1) { xv[j] = y*y ; j++; }

	uint16 pixel[3] ;
	tif_get3c(image,x,y,pixel[0],pixel[1],pixel[2]) ;

	sum_reg((&fit_struct[0]), xv, (double)pixel[0]) ;
	sum_reg((&fit_struct[1]), xv, (double)pixel[1]) ;
	sum_reg((&fit_struct[2]), xv, (double)pixel[2]) ;
    }

    double fit_coefs[4] ;
    double *coef_dest[3] = {rcoef,gcoef,bcoef} ;

    for(c=0 ; c < 3 ; c++) {
	double *coefs = coef_dest[c] ;
	estimate_reg((&fit_struct[c]), fit_coefs) ;

	int j=0 ;
	int k=0 ;
	coefs[0] = fit_coefs[0] ; //B0 is always used
	j++ ;
	k++ ;

	if(order_x) {  coefs[j] = fit_coefs[k] ; k++ ; j++; } else { coefs[j] = 0.0 ; j++;}
	if(order_y) {  coefs[j] = fit_coefs[k] ; k++ ; j++; } else { coefs[j] = 0.0 ; j++;}
	if(order_y > 1) {  coefs[j] = fit_coefs[k] ; k++ ; j++; } else { coefs[j] = 0.0 ; j++;}
    }

#ifdef GMR_DEBUG
    if(gmr_debug) {
	fprintf(stderr, "GMR: rcoef %f,%f,%f,%f\n", rcoef[0], rcoef[1], rcoef[2], rcoef[3]) ;
	fprintf(stderr, "GMR: gcoef %f,%f,%f,%f\n", gcoef[0], gcoef[1], gcoef[2], gcoef[3]) ;
	fprintf(stderr, "GMR: bcoef %f,%f,%f,%f\n", bcoef[0], bcoef[1], bcoef[2], bcoef[3]) ;
    }
#endif

    int n_repaired=0 ;
    if(repair_bad_points_flag == 1) {

	for(int i=0 ; i < n_possible_points ; i++) {

	    if(point_status[i] == 1)
		continue ;

	    x = x_possible_coords[i] ;

	    // pixels in column masks have already been eliminated in the very first loop
	    // now see if the repair mask is set for this column
	    if(column_repair_mask[x] == 1) continue ;
	   
	    y = y_possible_coords[i] ;

	    float rhat = rcoef[0] + rcoef[1]*(float)x + rcoef[2]*(float)y + rcoef[3]*(float)y*(float)y ;
	    float ghat = gcoef[0] + gcoef[1]*(float)x + gcoef[2]*(float)y + gcoef[3]*(float)y*(float)y ;
	    float bhat = bcoef[0] + bcoef[1]*(float)x + bcoef[2]*(float)y + bcoef[3]*(float)y*(float)y ;
#ifdef GMR_DEBUG
	    if(gmr_debug) {
		fprintf(stderr, "   Set pixel %d,%d to  %5d,%5d,%5d\n", x,y,(int)rhat,(int)bhat,(int)ghat) ;
	    }
#endif

	    tif_set4c(image,x,y,rhat,ghat,bhat,MAX16) ;

	    n_repaired++ ;
	}
    }

    free(point_status) ;
    free(x_possible_coords) ;
    free(y_possible_coords) ;
    free(x_coords) ;
    free(y_coords) ;
#ifdef GMR_DEBUG
    if(gmr_debug) {
	if(n_possible_points == 121 && n_points == 0) exit(1) ;
    }
#endif

    return n_repaired ;
}

void repair_top_of_sky(tdata_t *image, int end_of_sky_is_known_flag, int phase)
{
    int x, y ;
    fprintf(stderr, "RTOS::: sky hsv (%3.0f,%4.2f,%4.2f)\n", hue_sky, sat_sky, val_sky) ;

    if(phase == 1) {
	// check "edges of found sky for issues:
	// 1. are there black pixels (with alpha near MAX16 ?
	fprintf(stderr, "RTOS eliminate black pixels at top of sky\n") ;

	for(x = 0 ; x < IMAGE_WIDTH ; x++) {
	    if(column_mask[x] == 1) continue ;
	    int n_ok = 0 ;
	    for(y = raw_start_of_sky[x]-1 ; y > 0 ; y--) {
		if(xy_is_black_pixel(image,x,y) == 1) {
		    fprintf(stderr, "RTOS found black pixel at %d,%d\n",x,y) ;
		    ((uint16 *)(image[y]))[IMAGE_NSAMPLES*x+3] = 0 ;
		    n_ok=0 ;
		} else if(n_ok > 10) {
		    break ;
		} else {
		    n_ok++ ;
		}
	    }
	}

	fprintf(stderr, "RTOS eliminate black pixels on vertical edges of sky\n") ;
	set_minmax_sky_values() ;

	for(y = min_start_of_sky ; y < max_start_of_sky+1 ; y++) {

	    // left to right
	    int max_x=IMAGE_WIDTH-1 ;
	    for(x = 0 ; x < max_x ; x++) {
		if(xy_is_transparent(image,x,y) && xy_is_black_pixel(image,x+1,y) == 1) {
		    if(column_mask[x+1] == 0) continue ;
		    fprintf(stderr, "RTOS LR fix x,y %d,%d\n", x,y) ;
		    ((uint16 *)(image[y]))[IMAGE_NSAMPLES*(x+1)+3] = 0 ;
		    max_x = x+10 ;
		    if(max_x > IMAGE_WIDTH-1) max_x = IMAGE_WIDTH-1 ;
		}
	    }

	    // right to left
	    int min_x =0 ;
	    for(x = IMAGE_WIDTH-1 ; x > min_x ; x--) {
		if(xy_is_transparent(image,x,y) && xy_is_black_pixel(image,x-1,y) == 1) {
		    if(column_mask[x-1] == 0) continue ;
		    fprintf(stderr, "RTOS RL fix x,y %d,%d\n", x,y) ;
		    ((uint16 *)(image[y]))[IMAGE_NSAMPLES*(x-1)+3] = 0 ;
		    min_x = x-10 ;
		    if(min_x < 0) min_x = 0 ;
		}
	    }

	}

	fprintf(stderr, "RTOS first scan for bad pixels in top 5 pixels of sky\n") ;

	double rcoef[4], gcoef[4], bcoef[4] ;

	// are there wrong pixels in the top 50 % of the sky?

	int xc = 0 ;
	while(1) {
	    while(column_mask[xc] == 1 && xc < IMAGE_WIDTH)
		xc++ ;

	    if(xc >= IMAGE_WIDTH)
		break ;

	    int search_width = 5 ;
	    int yc = start_of_sky[x] ;

	    int n_repaired = get_mean_rgb(image,xc,yc,rcoef,gcoef,bcoef,search_width,end_of_sky_is_known_flag,1,1) ;

	    if(n_repaired > 0)
		fprintf(stderr, "Bad pixels near top of sky: Repaired %d pixel colors near %d,%d\n", n_repaired, xc, yc) ;
	    xc++ ;
	}

    }

    if(phase == 1) {
	// fill in vertical gaps where alpha is 0, happens
	// when there is a mask in the image that hugin used that
	// appears on the edge

	// first, if there are gaps at the edges of the image,
	// these can extend below the actual sky, so they need
	// to be detected.
	int left_edge_max_x=-1 ;

	for(x = 0 ; x < IMAGE_WIDTH ; x++) {
	    int found=0 ;
	    for(y = raw_start_of_sky[x] ; y <(IMAGE_HEIGHT*3)/4 ; y++) {
		if(xy_is_opaque_pixel(image,x,y) == 0) {
		    left_edge_max_x=x ;
		    found=1 ;
		    break ;
		}
	    }
	    if(found == 0)
		break ;
	}

	int right_edge_min_x=IMAGE_WIDTH ;

	for(x = IMAGE_WIDTH-1 ; x > -1 ; x--) {
	    int found=0 ;
	    for(y = raw_start_of_sky[x] ; y <(IMAGE_HEIGHT*3)/4 ; y++) {
		if(xy_is_opaque_pixel(image,x,y) == 0) {
		    right_edge_min_x=x ;
		    found=1 ;
		    break ;
		}
	    }
	    if(found == 0)
		break ;
	}

	fprintf(stderr, "RTOS fill vertical gap,  leftedgemask:%d,  rightedgemask,%d\n", left_edge_max_x, right_edge_min_x) ;

	for(x = left_edge_max_x+1 ; x < right_edge_min_x-1 ; x++) {
	    if(column_mask[x] == 1) continue ;
	    int maxy = end_of_sky_is_known_flag ? end_of_sky[x] : (IMAGE_HEIGHT*3)/4 ;
	    maxy = (IMAGE_HEIGHT*3)/4 ;

	    for(y = raw_start_of_sky[x] ; y < maxy ; y++) {

		if(xy_is_opaque_pixel(image,x,y) == 0) {
		    int ty = y ;
		    int by ;

		    for(by = y+1 ; by < maxy ; by++) {
			if(xy_is_opaque_pixel(image,x,by) == 1)
			    break ;
		    }

		    if(by == maxy)
			continue ; // can't fix transparent pixels at bottom of image

		    fprintf(stderr, "\nRTOS fill vertical gap at x:%d,  y %d to %d\n", x, ty, by) ;

		    int search_width = 50 ;

		    // the top and bottom pixels may be badly colored pixels from edges
		    // of images that hugin output, correct those while getting prediction
		    // coefs
		    double rcoef0[4], gcoef0[4], bcoef0[4] ;
		    double rcoef1[4], gcoef1[4], bcoef1[4] ;
		    int nt = get_mean_rgb(image,x,ty,rcoef0,gcoef0,bcoef0,search_width,end_of_sky_is_known_flag,1,2) ;
		    if(nt > 0)
			fprintf(stderr, "Repaired %d pixels near the top of the gap\n",nt) ;
		    int nb = get_mean_rgb(image,x,by,rcoef1,gcoef1,bcoef1,search_width,end_of_sky_is_known_flag,1,2) ;
		    if(nb > 0)
			fprintf(stderr, "Repaired %d pixels near the top of the gap\n",nb) ;


		    for(y = ty ; y <= by ; y++) {
			float p1 = (float)(y-ty)/(by-ty) ;
			float p0 = 1.f - p1 ;

			float rhat0 = rcoef0[0] + rcoef0[1]*(float)x + rcoef0[2]*(float)y + rcoef0[3]*(float)y*(float)y ;
			float ghat0 = gcoef0[0] + gcoef0[1]*(float)x + gcoef0[2]*(float)y + gcoef0[3]*(float)y*(float)y ;
			float bhat0 = bcoef0[0] + bcoef0[1]*(float)x + bcoef0[2]*(float)y + bcoef0[3]*(float)y*(float)y ;

			float rhat1 = rcoef1[0] + rcoef1[1]*(float)x + rcoef1[2]*(float)y + rcoef1[3]*(float)y*(float)y ;
			float ghat1 = gcoef1[0] + gcoef1[1]*(float)x + gcoef1[2]*(float)y + gcoef1[3]*(float)y*(float)y ;
			float bhat1 = bcoef1[0] + bcoef1[1]*(float)x + bcoef1[2]*(float)y + bcoef1[3]*(float)y*(float)y ;

			float rhat = p0*rhat0 + p1*rhat1 ;
			float ghat = p0*ghat0 + p1*ghat1 ;
			float bhat = p0*bhat0 + p1*bhat1 ;

			tif_set4c(image,x,y, (uint16)rhat, (uint16)ghat, (uint16)bhat, MAX16 );
		    }
		}

		// the search continues for the next greater values of y
	    }
	}

	set_minmax_sky_values() ;

	fprintf(stderr, "RTOS fill divots\n") ;

	// look for "divots" at the top of the sky to fill in
	// this is essentially filling horizontal gaps, where a gap is black pixels or alpha == 0 ;

	for(y = max_start_of_sky ; y >= min_start_of_sky ; y--) {
	    continue ;

	    int lx, rx ;

	    for(x = 0 ; x < IMAGE_WIDTH ;) {
		if(column_mask[x] == 1) {
		    x++ ;
		    continue ;
		}

		// find lx, which is next x where SOS[x] <= y and SOS[x+1] > y
		lx=-1 ;

		while(x < IMAGE_WIDTH-2) {
		    if(column_mask[x] == 1) {
			x++ ;
			continue ;
		    }
		    if(xy_has_nonblack_pixel(image, x, y) == 1)
			lx=x ;

		    if(start_of_sky[x+1] > y) {
			break ;
		    }

		    x++ ;
		}

		if(lx > -1) {
		    // find rx, which is next x where SOS[x] <= y ;
		    rx=-1 ;
		    x=lx+1 ;
		    while(x < IMAGE_WIDTH-1) {
			if(column_mask[x] == 1) {
			    x++ ;
			    continue ;
			}
			if(xy_has_nonblack_pixel(image, x, y) == 1) {
			    rx=x ;
			    break ;
			}
			x++ ;
		    }
		}

		if(lx > -1 && rx > -1) {
		    float width = rx-lx ;

		    if(column_mask[rx]) {
			fprintf(stderr, "RTOS: Fatal rx=%d is masked column!\n", rx) ;
			exit(1) ;
		    }

		    if(column_mask[lx]) {
			fprintf(stderr, "RTOS: Fatal lx=%d is masked column!\n", lx) ;
			exit(1) ;
		    }

		    if(width < IMAGE_WIDTH) {
			fprintf(stderr, "RTOS fill horizontal gap,  y:%d, l:%d, r:%d\n", y, lx, rx) ;

			double rcoef0[4], gcoef0[4], bcoef0[4] ;
			double rcoef1[4], gcoef1[4], bcoef1[4] ;

			int search_width=5 ;

			get_mean_rgb(image,lx,y,rcoef0,gcoef0,bcoef0,search_width,end_of_sky_is_known_flag,1,2) ;
			get_mean_rgb(image,rx,y,rcoef1,gcoef1,bcoef1,search_width,end_of_sky_is_known_flag,1,2) ;

			for(x = lx+1 ; x < rx ; x++) {
			    if(column_mask[x] == 1) {
				continue ;
			    }

			    float p1 = (float)(x-lx)/width ;
			    float p0 = 1. - p1 ;

			    float rhat0 = rcoef0[0] + rcoef0[1]*(float)x + rcoef0[2]*(float)y + rcoef0[3]*(float)y*(float)y ;
			    float ghat0 = gcoef0[0] + gcoef0[1]*(float)x + gcoef0[2]*(float)y + gcoef0[3]*(float)y*(float)y ;
			    float bhat0 = bcoef0[0] + bcoef0[1]*(float)x + bcoef0[2]*(float)y + bcoef0[3]*(float)y*(float)y ;

			    float rhat1 = rcoef1[0] + rcoef1[1]*(float)x + rcoef1[2]*(float)y + rcoef1[3]*(float)y*(float)y ;
			    float ghat1 = gcoef1[0] + gcoef1[1]*(float)x + gcoef1[2]*(float)y + gcoef1[3]*(float)y*(float)y ;
			    float bhat1 = bcoef1[0] + bcoef1[1]*(float)x + bcoef1[2]*(float)y + bcoef1[3]*(float)y*(float)y ;

			    float rhat = p0*rhat0 + p1*rhat1 ;
			    float ghat = p0*ghat0 + p1*ghat1 ;
			    float bhat = p0*bhat0 + p1*bhat1 ;

			    tif_set4c(image,x,y, (uint16)rhat, (uint16)ghat, (uint16)bhat, MAX16 );

			    start_of_sky[x] = y ;
			}
		    } else {
			x=rx ;
		    }

		} else {
		    x++ ;
		}
	    }

	}

    }

    // note, there is no phase 2, just so we can break at a debug level that matches the phase

    if(phase == 3) {
	fprintf(stderr, "RTOS phase 3: fill side edges\n") ;
	// do the side edges need to be fixed?
	int edge_width = IMAGE_WIDTH/40 ;
	int min_start_y=IMAGE_HEIGHT ;
	int max_start_y=0 ;

	// LEFT side
	for(x = 0 ; x < edge_width ; x++) {
	    if(min_start_y > start_of_sky[x]) min_start_y = start_of_sky[x] ;
	    if(max_start_y < start_of_sky[x]) max_start_y = start_of_sky[x] ;
	}

	for(int y = max_start_y ; y >= min_start_y ; y--) {
	    int xgood ;
	    for(xgood = 0 ; xy_has_nonblack_pixel(image, xgood, y) == 0 && xgood < edge_width ; xgood++) ;
	    if(xgood == 0) continue ; // left edge pixel is filled

	    if(xgood == edge_width) {
		fprintf(stderr, "WARNING: could not find good start pixel on left edge gap at y=%d\n", y) ;
		continue ;
	    }

	    double rcoef[4], gcoef[4], bcoef[4] ;
	    int search_width = (xgood-0)*10 ;
	    if(search_width > 20) search_width=20 ;
	    get_mean_rgb(image,xgood,y,rcoef,gcoef,bcoef,search_width,end_of_sky_is_known_flag,0,1) ;

	    fprintf(stderr, "RTOS fix sky left edge from x,y %d,%d to %d,%d\n", 0, y, xgood-1, y) ;
	    for(x = 0 ; x < xgood ; x++) {
		float rhat = rcoef[0] + rcoef[1]*(float)x + rcoef[2]*(float)y + rcoef[3]*(float)y*(float)y ;
		float ghat = gcoef[0] + gcoef[1]*(float)x + gcoef[2]*(float)y + gcoef[3]*(float)y*(float)y ;
		float bhat = bcoef[0] + bcoef[1]*(float)x + bcoef[2]*(float)y + bcoef[3]*(float)y*(float)y ;

		tif_set4c(image,x,y, (uint16)rhat, (uint16)ghat, (uint16)bhat, MAX16 );
		start_of_sky[x] = y ;
	    }

	}

	// RIGHT side
	min_start_y=IMAGE_HEIGHT ;
	max_start_y=0 ;
	int inside_x = IMAGE_WIDTH-1-edge_width ;

	for(x = inside_x ; x < IMAGE_WIDTH ; x++) {
	    if(min_start_y > start_of_sky[x]) min_start_y = start_of_sky[x] ;
	    if(max_start_y < start_of_sky[x]) max_start_y = start_of_sky[x] ;
	}

	for(int y = max_start_y ; y >= min_start_y ; y--) {
	    int xgood ;
	    for(xgood = IMAGE_WIDTH-1 ; xy_has_nonblack_pixel(image, xgood, y) == 0 && xgood > inside_x ; xgood--) ;
	    if(xgood == IMAGE_WIDTH-1) continue ; // right edge pixel is filled

	    if(xgood == inside_x) {
		fprintf(stderr, "WARNING: could not find good start pixel on right edge gap at y=%d\n", y) ;
		continue ;
	    }

	    double rcoef[4], gcoef[4], bcoef[4] ;
	    int search_width = (IMAGE_WIDTH-xgood)*10 ;
	    if(search_width > 20) search_width=20 ;
	    get_mean_rgb(image,xgood,y,rcoef,gcoef,bcoef,search_width,end_of_sky_is_known_flag,0,1) ;

	    fprintf(stderr, "RTOS fix sky right edge from x,y %d,%d to %d,%d\n", xgood+1, y, IMAGE_WIDTH-1,y) ;

	    for(x = xgood+1 ; x < IMAGE_WIDTH ; x++) {
		float rhat = rcoef[0] + rcoef[1]*(float)x + rcoef[2]*(float)y + rcoef[3]*(float)y*(float)y ;
		float ghat = gcoef[0] + gcoef[1]*(float)x + gcoef[2]*(float)y + gcoef[3]*(float)y*(float)y ;
		float bhat = bcoef[0] + bcoef[1]*(float)x + bcoef[2]*(float)y + bcoef[3]*(float)y*(float)y ;

		//fprintf(stderr, "    x:%d rgb:%d,%d,%d\n", x, (int)rhat,(int)ghat,(int)bhat) ;

		tif_set4c(image,x,y, (uint16)rhat, (uint16)ghat, (uint16)bhat, MAX16 );
		start_of_sky[x] = y ;
	    }

	}



    }

    if(phase == 4) {
#ifdef GMR_DEBUG
	gmr_debug = 1000 ;
#endif
	fprintf(stderr, "RTOS second scan for bad pixels in top 50%% of sky\n") ;

	double rcoef[4], gcoef[4], bcoef[4] ;
	// are there wrong pixels in the top 50 % of the sky?


	int search_width = 5 ;
	int n_repaired_gt=0 ; // grand total repaired ;

	for(int xc=search_width ; xc < IMAGE_WIDTH-5 ; xc += search_width) {
	    int y0 = 0 ;
	    int y1 = IMAGE_HEIGHT-1 ;
	    int x0 = xc-search_width ;
	    int x1 = xc+search_width ;
	    if(x1 > IMAGE_WIDTH-1) x1 = IMAGE_WIDTH-1 ;
	    int n_unmasked_columns=0 ;

	    for(int x = x0 ; x <= x1 ; x++) {
		if(column_mask[x] == 1) continue ;
		if(column_repair_mask[x] == 1) continue ;
		if(y0 < start_of_sky[x]) y0 = start_of_sky[x] ;
		if(y1 > end_of_sky[x]) y1 = end_of_sky[x] ;
		n_unmasked_columns++ ;
	    }

	    if(n_unmasked_columns == 0) continue ;

	    y1 = (y0 + y1)/2. ;
/*  	    for(int yc=y0-search_width ; yc < y1-search_width ; yc += search_width) {  */
	    for(int yc=y1-search_width ; yc > y0 ; yc -= search_width) {
		int n_repaired = get_mean_rgb(image,xc,yc,rcoef,gcoef,bcoef,search_width,end_of_sky_is_known_flag,1,2) ;
		n_repaired_gt += n_repaired ;
/*  		if(n_repaired > 0)  */
/*  		    fprintf(stderr, "Bad pixels in top 50%% of sky: clip EOS=%d,repaired %d pixel colors near %d,%d\n", end_of_sky_is_known_flag, n_repaired, xc, yc) ;  */

#ifdef GMR_DEBUG
		gmr_debug-- ;
		if(gmr_debug == 0)
		    exit(1) ;
#endif
	    }
	}

	fprintf(stderr, "Bad pixels in top 50%% of sky: repaired %d pixel colors\n", n_repaired_gt) ;
    }

    // are there sky slopes so steep they will affect the result?
    if(phase == 5) {

        // run back and forth until there is no slope too steep
	int needs_more=1 ;

        while(needs_more) {
	    needs_more=0 ;
	    fprintf(stderr, "RTOS slope correct\n") ;
	    int search_width = 15 ;

	    double rcoef[4], gcoef[4], bcoef[4] ;

	    for(x = 0 ; x < IMAGE_WIDTH-4 ; x++) {
		int x2 = x+3 ;
		if(start_of_sky[x] > start_of_sky[x2]+1) {
		    needs_more=1 ;
		    y=start_of_sky[x] ;
		    get_mean_rgb(image,x,start_of_sky[x],rcoef,gcoef,bcoef,search_width,end_of_sky_is_known_flag,0,2) ;

		    float rhat = rcoef[0] + rcoef[1]*(float)x + rcoef[2]*(float)y + rcoef[3]*(float)y*(float)y ;
		    float ghat = gcoef[0] + gcoef[1]*(float)x + gcoef[2]*(float)y + gcoef[3]*(float)y*(float)y ;
		    float bhat = bcoef[0] + bcoef[1]*(float)x + bcoef[2]*(float)y + bcoef[3]*(float)y*(float)y ;

		    tif_set4c(image,x,y, (uint16)(rhat+0.5), (uint16)(ghat+0.5), (uint16)(bhat+0.5), MAX16 );

		    start_of_sky[x]-- ;
		}
		if(start_of_sky[x]+1 < start_of_sky[x2]) {
		    needs_more=1 ;
		    y=start_of_sky[x2] ;
		    get_mean_rgb(image,x2,start_of_sky[x2],rcoef,gcoef,bcoef,search_width,end_of_sky_is_known_flag,0,2) ;

		    float rhat = rcoef[0] + rcoef[1]*(float)x2 + rcoef[2]*(float)y + rcoef[3]*(float)y*(float)y ;
		    float ghat = gcoef[0] + gcoef[1]*(float)x2 + gcoef[2]*(float)y + gcoef[3]*(float)y*(float)y ;
		    float bhat = bcoef[0] + bcoef[1]*(float)x2 + bcoef[2]*(float)y + bcoef[3]*(float)y*(float)y ;

		    tif_set4c(image,x2,y, (uint16)(rhat+0.5), (uint16)(ghat+0.5), (uint16)(bhat+0.5), MAX16 );
		    start_of_sky[x2]-- ;
		}
	    }
	}
    }

    if(fill_top_of_sky == 1 && phase == 6) {
	fprintf(stderr, "RTOS fill gaps to min start of sky\n") ;

	for(x = 0 ; x < IMAGE_WIDTH ; x++) {
	    if(column_mask[x] == 1) continue ;

	    if(start_of_sky[x] > min_start_of_sky) {
		int search_width = 30 ; // need a robust estimate, so wide search width
		double rcoef[4], gcoef[4], bcoef[4] ;
		get_mean_rgb(image,x,start_of_sky[x],rcoef,gcoef,bcoef,search_width,end_of_sky_is_known_flag,1,2) ;

		for(y = start_of_sky[x] ; y > min_start_of_sky ; y--) {
		    float rhat = rcoef[0] + rcoef[1]*(float)x + rcoef[2]*(float)y + rcoef[3]*(float)y*(float)y ;
		    float ghat = gcoef[0] + gcoef[1]*(float)x + gcoef[2]*(float)y + gcoef[3]*(float)y*(float)y ;
		    float bhat = bcoef[0] + bcoef[1]*(float)x + bcoef[2]*(float)y + bcoef[3]*(float)y*(float)y ;

		    tif_set4c(image,x,y, (uint16)(rhat+0.5), (uint16)(ghat+0.5), (uint16)(bhat+0.5), MAX16 );
		}

		start_of_sky[x] = min_start_of_sky ;


	    }
	}
    }

    fprintf(stderr, "Finished RTOS::: sky hsv (%3.0f,%4.2f,%4.2f)\n", hue_sky, sat_sky, val_sky) ;

}

void estimate_sky(int x0,int x1,tdata_t *image,int16 *start_of_sky,int16 *end_of_sky,int16 *final_end_of_sky,
		    float depth_of_sample, int h, int w, float feather_factor, int debug, float depth_of_fill, int extra_feather_length)
{
    int x, y ;
    fprintf(stderr, "estimate sky x:%d to %d, feather_factor=%f, nonlinear_feather=%f\n", x0,x1,feather_factor, nonlinear_feather) ;

    float sun_px = sun_x_angle2px(sun_x) ;

    if(x0 < 0) x0=0 ;
    if(x1 > w-1) x1=w-1 ;

    V_sun= image_relative_pixel_to_V3(sun_x_angle2px(sun_x), sun_py) ;

    maximum_vhat = find_maximum_vhat() ;

    float *vcorrect_at_x, *hcorrect_at_x, *scorrect_at_x ;
    float *avg_vcorrect_at_x, *avg_hcorrect_at_x, *avg_scorrect_at_x ;

    // actual h,s,v divided by predicted at column x
    hcorrect_at_x = (float *)calloc(IMAGE_WIDTH, sizeof(float)) ;
    scorrect_at_x = (float *)calloc(IMAGE_WIDTH, sizeof(float)) ;
    vcorrect_at_x = (float *)calloc(IMAGE_WIDTH, sizeof(float)) ;

    // moving average actual h,s,v divided by predicted at column x
    avg_hcorrect_at_x = (float *)calloc(IMAGE_WIDTH, sizeof(float)) ;
    avg_scorrect_at_x = (float *)calloc(IMAGE_WIDTH, sizeof(float)) ;
    avg_vcorrect_at_x = (float *)calloc(IMAGE_WIDTH, sizeof(float)) ;

    // is the sky saturation prediction not right, saturation should increases from max start of sky to y=0 ;
    {
	struct mstat m_fit = init_reg(1) ; // will be a simple linear in y
	for(x = x0 ; x <= x1 ; x += (x1-x0)/20) {
	    if(column_mask[x] == 1) continue ;
	    for(y = 0 ; y < max_start_of_sky ; y++) {
		float px = IMAGE_PIXEL_X_TO_RELATIVE(x) ;
		float py = IMAGE_PIXEL_Y_TO_RELATIVE(y) ;

		float hhat,shat,vhat ;
		predict_sky_hsv(px, py, &hhat, &shat, &vhat) ;
		double xv[1] = {py} ;
		sum_reg(&m_fit, xv, shat) ;
	    }
	}

	double coefs[2] ;
	estimate_reg(&m_fit, coefs) ;

	if(coefs[1] < 0.0) {
	    fprintf(stderr, "*********************************************************************************\n") ;
	    fprintf(stderr, "*********************************************************************************\n") ;
	    fprintf(stderr, "*****************  WARNING WARNING WARNING **************************************\n") ;
	    fprintf(stderr, "It appears that the saturation prediction is not correct\n") ;

	    if(sat_prediction_method == 1) {
		fprintf(stderr, " --reverting to a simpler\n") ;
		fprintf(stderr, "prediction method -- but a better fit could be had by adding \"-spm 0\" to the command line\n") ;
	    } else {
		fprintf(stderr, " results may look unnatural.\n") ;
	    }

	    fprintf(stderr, "*********************************************************************************\n") ;
	    fprintf(stderr, "*********************************************************************************\n") ;

	    S_from_V_coefs[0] = sat_sky ;
	    S_from_V_coefs[1] = 0. ;
	    S_from_V_coefs[2] = 0. ;
	} else {
	    fprintf(stderr, "*********************************************************************************\n") ;
	    fprintf(stderr, "Modeled sky saturation coefs are b0:%f b1:%f\n", coefs[0], coefs[1]) ;
	    fprintf(stderr, "*********************************************************************************\n") ;
	}

    }


    if(show_raw_prediction || estimate_only) {
	for(x = x0 ; x <= x1 ; x++) {
	    if(column_mask[x] == 1) continue ;
	    int ymax = show_raw_prediction ? start_of_sky[x] : IMAGE_HEIGHT ;

	    for(y = 0 ; y < ymax ; y++) {
		// completely model sky color
		float px = IMAGE_PIXEL_X_TO_RELATIVE(x) ;
		float py = IMAGE_PIXEL_Y_TO_RELATIVE(y) ;

		float hhat,shat,vhat ;
		uint16 rhat, ghat, bhat ;

		predict_sky_hsv(px, py, &hhat, &shat, &vhat) ;
		shat *= final_saturation_factor ;

		hsv2rgb16(hhat,shat,vhat,&rhat,&ghat,&bhat) ;

		tif_set4c(image,x,y, (uint16)(rhat+0.5), (uint16)(ghat+0.5), (uint16)(bhat+0.5), MAX16 );
	    }

	}

	goto estimate_sky_cleanup ;
    }

    // need to generate correction factors by column

    for(x = x0 ; x <= x1 ; x++) {
	if(column_mask[x] == 1) continue ;


	// to collect correction data for this column only look 85% into the column
	float sky_height = end_of_sky[x] - start_of_sky[x] ;
	int feather_end_y ;
/*  	feather_end_y = start_of_sky[x] + (int)(sky_height*.10+0.5) ;  */
/*  	feather_end_y = start_of_sky[x] + 30 ;  */
	feather_end_y = start_of_sky[x] + (int)(sky_height*.25+0.5) ;

	if(feather_end_y > start_of_sky[x] + 10) feather_end_y = start_of_sky[x] + 10 ;

	if(feather_end_y > IMAGE_HEIGHT-1) feather_end_y = IMAGE_HEIGHT-1 ;
	if(feather_end_y > end_of_sky[x]) feather_end_y = end_of_sky[x] ;

	float sum_h=0., sum_s=0., sum_v=0. ;
	float sum_hhat=0., sum_shat=0., sum_vhat=0. ;

	float sum_wgts=0. ;

	for(y = start_of_sky[x]+5 ; y < feather_end_y ; y++) {
	    float hhat,shat,vhat ;

	    if(xy_has_nonblack_pixel(image, x, y) == 0) continue ;
	    sum_wgts=0. ;

	    // completely model sky color
	    float px = IMAGE_PIXEL_X_TO_RELATIVE(x) ;
	    float py = IMAGE_PIXEL_Y_TO_RELATIVE(y) ;

	    predict_sky_hsv(px, py, &hhat, &shat, &vhat) ;

	    uint16 r,g,b ;
	    tif_get3c(image,x,y,r,g,b) ;

	    float h,s,v ;
	    rgb2hsv16(r,g,b,&h,&s,&v) ;

	    sum_h += h*s ;
	    sum_s += s ;
	    sum_v += v*s ;

	    sum_hhat += hhat*s ;
	    sum_shat += shat ;
	    sum_vhat += vhat*s ;
	    sum_wgts += s ;
	}

	if(sum_wgts == 0) {
	    // can happen at edges
	    hcorrect_at_x[x] = 1.f ;
	    scorrect_at_x[x] = 1.f ;
	    vcorrect_at_x[x] = 1.f ;
	} else {
	    hcorrect_at_x[x] = sum_h/sum_hhat ;
	    scorrect_at_x[x] = sum_s/sum_shat ;
	    vcorrect_at_x[x] = sum_v/sum_vhat ;
	}
    }

    if(0 && allowed_sky_type == UNIFORM_CIE_SKIES) {
	float sum_h=0., sum_s=0., sum_v=0. ;
	int n_correct = 0 ;

	for(x = x0 ; x <= x1 ; x++) {
	    if(column_mask[x0] == 1) continue ;

	    // vcorrect is 0. if there were no samples to calculate it
	    if(vcorrect_at_x[x0] > 1.e-10) {
		sum_h += hcorrect_at_x[x] ;
		sum_s += scorrect_at_x[x] ;
		sum_v += vcorrect_at_x[x] ;
		n_correct++ ;
	    }
	}

	if(n_correct > 0) {
	    for(x = x0 ; x <= x1 ; x++) {
		hcorrect_at_x[x] = (hcorrect_at_x[x] + sum_h / (float)n_correct)/2. ;
		scorrect_at_x[x] = (scorrect_at_x[x] + sum_s / (float)n_correct)/2. ;
		vcorrect_at_x[x] = (vcorrect_at_x[x] + sum_v / (float)n_correct)/2. ;
	    }
	}
    }

    for(x = x0 ; x <= x1 ; x++) {
	if(column_mask[x] == 1) continue ;
	if(show_raw_prediction) continue ;

	float sum_h=0., sum_s=0., sum_v=0. ;
	int n_correct = 0 ;

	avg_hcorrect_at_x[x] = 1. ;
	avg_scorrect_at_x[x] = 1. ;
	avg_vcorrect_at_x[x] = 1. ;

	for(int dx=-126 ; dx < 127 ; dx++) {
	    if(x+dx < 0 || x+dx > IMAGE_WIDTH-1) continue ;
	    if(column_mask[x+dx] == 1) continue ;

	    // vcorrect is 0. if there were no samples to calculate it
	    if(vcorrect_at_x[x+dx] > 1.e-10) {
		sum_h += hcorrect_at_x[x+dx] ;
		sum_s += scorrect_at_x[x+dx] ;
		sum_v += vcorrect_at_x[x+dx] ;
		n_correct++ ;
	    }
	}

	if(feather_factor > .01 && n_correct > 0) {
	    avg_hcorrect_at_x[x] = sum_h / (float)n_correct ;
	    avg_scorrect_at_x[x] = sum_s / (float)n_correct ;
	    avg_vcorrect_at_x[x] = sum_v / (float)n_correct ;
	}
    }



    for(x = x0 ; x <= x1 ; x++) {
	if(column_mask[x] == 1) continue ;

	int feather_end_y ;
	int feather_length = compute_feather_length(x, &feather_end_y, depth_of_fill, feather_factor, extra_feather_length) ;


	for(y = 0 ; y < feather_end_y ; y++) {
	    uint16 a = ((uint16 *)(image[y]))[IMAGE_NSAMPLES*x+3] ;

	    if(a < HALF16 && y > raw_start_of_sky[x]) continue ; // do not adjust transparent pixels below original start of sky

	    uint16 rhat, ghat, bhat ;

	    // current pixel value
	    uint16 r,g,b ;
	    tif_get3c(image,x,y,r,g,b) ;

	    // completely model sky color
	    float px = IMAGE_PIXEL_X_TO_RELATIVE(x) ;
	    float py = IMAGE_PIXEL_Y_TO_RELATIVE(y) ;

	    float fp0 = compute_feather_at_y(x, y, feather_length, b, feather_factor) ;
	    float fp1 = 1.-fp0 ;

	    float hhat,shat,vhat ;
	    predict_sky_hsv(px, py, &hhat, &shat, &vhat) ;

	    shat *= final_saturation_factor ;

	    // are we at or above the sun?
	    if(sun_py < 1. && sun_px > -0.5 && sun_px < 0.5 && py >= sun_py && fabs(px-sun_px) < py-sun_py) {
		float px_dist = (py-sun_py) ;
		float fdist_from_center = fabs(px-sun_px)/(py-sun_py) ;

		// x0, x1 are left and right extents where we want to
		// modify values more directly 
		int x0 = IMAGE_RELATIVE_TO_PIXEL_X(sun_px-px_dist) ;
		int x1 = IMAGE_RELATIVE_TO_PIXEL_X(sun_px+px_dist) ;
		float px1=(float)(x-x0)/(float)(x1-x0) ;
		float px0=1.-px1 ;

/*  		if(fdist_from_center > .99) {  */
/*  		    fprintf(stderr, "sun est x,y:%d,%d, x0,x1:%d,%d px0, px1 %f,%f\n", x,y, x0,x1, px0, px1) ;  */
/*  		}  */

		hhat *= px0*avg_hcorrect_at_x[x0] + px1*avg_hcorrect_at_x[x1] ;
		shat *= px0*avg_scorrect_at_x[x0] + px1*avg_scorrect_at_x[x1] ;
		vhat *= px0*avg_vcorrect_at_x[x0] + px1*avg_vcorrect_at_x[x1] ;

		fp0 = 1. ;
		fp1 = 1.-fp0 ;
	    } else {
		hhat *= avg_hcorrect_at_x[x] ;
		shat *= avg_scorrect_at_x[x] ;
		vhat *= avg_vcorrect_at_x[x] ;
	    }

	    if(vhat < 0.f) vhat=0.f ;
	    if(vhat > 1.f) vhat=1.f ;
	    if(shat < 0.f) shat=0.f ;
	    if(shat > 1.f) shat=1.f ;
	    hsv2rgb16(hhat,shat,vhat,&rhat,&ghat,&bhat) ;



	    r = (uint16)(fp0*(float)rhat + fp1*(float)r+0.5) ;
	    g = (uint16)(fp0*(float)ghat + fp1*(float)g+0.5) ;
	    b = (uint16)(fp0*(float)bhat + fp1*(float)b+0.5) ;
	    tif_set4c(image,x,y, r,g,b,MAX16) ;
	}

	// did end of sky not go far enough?   See if there is an increase in L
	// at end of sky +1 followed by a decrease in L at end of sky +2

	int yeos ;
	final_end_of_sky[x] = end_of_sky[x] ;

	for(yeos = end_of_sky[x]-2 ; yeos < end_of_sky[x]+2 ; yeos++) {

	    float L_eos0 = get_xy_L(image, x, yeos) ;
	    float L_eos1 = get_xy_L(image, x, yeos+1) ;
	    float L_eos2 = get_xy_L(image, x, yeos+2) ;

	    if(L_eos1 > L_eos0 && L_eos1 > L_eos2) {
		float fac = L_eos0/L_eos1 ; // scale rgb values back to L found at end of sky
		y = yeos+1 ;
		final_end_of_sky[x] = y ;

		uint16 r,g,b,a ;
		tif_get4c(image,x,y,r,g,b,a) ;

		if(a < HALF16 && y > raw_start_of_sky[x]) continue ; // do not adjust transparent pixels below original start of sky

		tif_set3c(image,x,y, (uint16)((float)r*fac+0.5), (uint16)((float)g*fac+0.5), (uint16)((float)b*fac+0.5));
	    }
	}
    }

    //for(x = x0 ; x <= x1 ; x++) {
//	float sky_height = end_of_sky[x] - start_of_sky[x] ;
//	int feather_end_y = start_of_sky[x] + (int)(sky_height*depth_of_fill+0.5) ;
 //   }

estimate_sky_cleanup:

    free(hcorrect_at_x) ;
    free(scorrect_at_x) ;
    free(vcorrect_at_x) ;

    free(avg_hcorrect_at_x) ;
    free(avg_scorrect_at_x) ;
    free(avg_vcorrect_at_x) ;

    fprintf(stderr, "done with estimate sky, column=%d\n", x) ;

}

// a function to reduce horizontal variation of sky r,g,b values
void make_sky_uniform(tdata_t *image, int flag, int x0, int x1)
{
    // flag=0 => run smoother from top of image(y0) to lowest end of sky(y1)
    // flag=1 => run smoother from top of image(y0) to lowest start of sky(y1)
    // flag=2 => run smoother from highest start of sky(y0) to lowest end of sky(y1)

    // if flag is 0 or 1, this smoother is run after the entire sky is predicted
    // if flag is 2, this smoother is run before the optimization and sky prediction

/*      int min_esky_y=IMAGE_HEIGHT-1 ;  */
    int min_ssky_y=IMAGE_HEIGHT-1 ;
    int max_esky_y=0 ;
    int max_ssky_y=0 ;

    for(int x = x0 ; x <= x1 ; x++) {
/*  	if(min_esky_y > end_of_sky[x]) min_esky_y = end_of_sky[x] ;  */
	if(min_ssky_y > start_of_sky[x]) min_ssky_y = start_of_sky[x] ;
	if(max_esky_y < end_of_sky[x]) max_esky_y = end_of_sky[x] ;
	if(max_ssky_y < start_of_sky[x]) max_ssky_y = start_of_sky[x] ;
    }

    int y0 = 0 ;
    if(flag == 2) y0 = min_ssky_y ;

    int y1=max_esky_y ;
    if(flag == 1) y1 = max_ssky_y ;

    fprintf(stderr, "make sky uniform, y %d to %d\n", y0, y1) ;

    for(int y = y0 ; y <= y1 ; y++) {
	float sum[3] = {0.,0.,0.} ;
	int n=0 ;
	for(int x = x0 ; x <= x1 ; x++) {
	    if(end_of_sky[x] >= y) {
		if(flag == 2 && y < start_of_sky[x]) continue ;

		uint16 r,g,b ;
		tif_get3c(image,x,y,r,g,b) ;

		sum[0] += r ;
		sum[1] += g ;
		sum[2] += b ;
		n++ ;
	    }
	}

	if(n > 0) {
	    sum[0] /= (float)n ;
	    sum[1] /= (float)n ;
	    sum[2] /= (float)n ;

    // flag == 0, blend is 0 at eos[x], to 1 at sos[x]
    // flag == 1, blend is 0 y1, to 1 at sos[x]
    // flag == 2, blend is 0 y1, to 1 at (eos[x]+sos[x])/2

	    for(int x = x0 ; x <= x1 ; x++) {
		if(y <= start_of_sky[x]) {
		    tif_set4c(image,x,y, sum[0], sum[1], sum[2], MAX16);
		    continue ;
		}

		if(flag == 0) {
		    int blend_y0 = (int)((start_of_sky[x]+end_of_sky[x])/2.+0.5) ;
		    int blend_y1 = end_of_sky[x] ;
		    if(y > end_of_sky[x]) continue ;

		    if(y <= blend_y0) {
			tif_set4c(image,x,y, sum[0], sum[1], sum[2], MAX16);
			continue ;
		    }

		    float pmean = 1.-(float)(y-blend_y0)/(float)(blend_y1-blend_y0) ;
		    float r = ((uint16 *)(image[y]))[IMAGE_NSAMPLES*x+0] ;
		    float g = ((uint16 *)(image[y]))[IMAGE_NSAMPLES*x+1] ;
		    float b = ((uint16 *)(image[y]))[IMAGE_NSAMPLES*x+2] ;
		    r = pmean*sum[0] + (1.-pmean)*r ;
		    g = pmean*sum[1] + (1.-pmean)*g ;
		    b = pmean*sum[2] + (1.-pmean)*b ;
		    tif_set4c(image,x,y,r,g,b, MAX16) ;
		    continue ;
		}

		float blend_y0 = start_of_sky[x] ;
		float blend_y1 = y1 ;

		if(flag == 2) blend_y0 = (float)(end_of_sky[x]+start_of_sky[x])/2. ;
		if(flag == 1 && blend_y1 > (float)end_of_sky[x]) blend_y1 = (float)end_of_sky[x] ;

		if(blend_y1 > end_of_sky[x]) blend_y1 = (float)end_of_sky[x] ;
		if(flag == 1) blend_y1 = (float)(end_of_sky[x]+start_of_sky[x])/2. ;
		if(flag == 1) blend_y1 = end_of_sky[x] ;

		float pmean = 1.-((float)y-blend_y0)/(float)(blend_y1-blend_y0) ; //at y == y0, use 100% of mean value, linear change to 0% of mean at y == y1

		if(pmean > 0.) {
		    float r = ((uint16 *)(image[y]))[IMAGE_NSAMPLES*x+0] ;
		    float g = ((uint16 *)(image[y]))[IMAGE_NSAMPLES*x+1] ;
		    float b = ((uint16 *)(image[y]))[IMAGE_NSAMPLES*x+2] ;
		    r = pmean*sum[0] + (1.-pmean)*r ;
		    g = pmean*sum[1] + (1.-pmean)*g ;
		    b = pmean*sum[2] + (1.-pmean)*b ;
		    tif_set4c(image,x,y,r,g,b, MAX16) ;
		}
	    }

	}
    }

    //should we move start of sky up?
    if(flag == 2) {
	// since this happened *before* sky sampling/optimization
	for(int x = x0 ; x <= x1 ; x++) {
	    start_of_sky[x] = y0 ;
	}
    }
}

void fit_sky_model(int force_grid_optimization, int n_custom_optimizations, int optimization_var[10][MAX_PARMS+1], int n_samples)
{
    model_is_being_fit=1 ;

    // limit high value of horizon_py
    opt_parms[2].hi = ((float)(IMAGE_HEIGHT)-(float)max_end_of_sky)/(float)IMAGE_HEIGHT ;

    if(force_grid_optimization || n_custom_optimizations == 0) {


	// default case, do a grid search for sun position
	float best_horizon_py ;
	float best_s_x ;
	float best_s_y ;
	float best_err=FMAX ;
	int n_tot = 0 ;
	float sun_dx = 1.e30 ;
	float sun_dy = 1.e30 ;
	float horizon_dy = 1.e30 ;
	opt_parms[0].lo = -180. ;
	opt_parms[0].hi = 180. ;
	opt_parms[1].lo = 0. ;
	opt_parms[1].hi = 1.5 ;
	opt_parms[2].hi = 1. ;

	fprintf(stderr, "1  Sunx,suny,horz_y=%f,%f,%f\n", sun_x, sun_py, horizon_py) ;

	if(opt_parms[0].grid_optimize_flag == 1) {
	    sun_x=opt_parms[0].lo ;
	    sun_dx = 10. ;
	    fprintf(stderr, "Grid optimize sun_x\n") ;
	}
	if(opt_parms[1].grid_optimize_flag == 1) {
	    sun_py=opt_parms[1].lo ;
	    sun_dy = .5 ;
	    fprintf(stderr, "Grid optimize sun_py\n") ;
	}

	if(opt_parms[2].grid_optimize_flag == 1) {
	    horizon_py = find_horizon_py(n_samples) ;
	    // now, horizon_curvature is also set

	    // if(n_custom_optimizations == 0) goto estimate_sky ;
	    //horizon_py=opt_parms[2].lo ;
	    //horizon_dy = .1 ;
	    fprintf(stderr, "Grid optimize horizon_py found is %f\n",horizon_py) ;
	}

	float start_sun_x = sun_x ;
	float start_sun_py = sun_py ;
	float start_horizon_py = horizon_py ;

	fprintf(stderr, "Grid optimize sun_x from %f to %f by %f\n", start_sun_x, opt_parms[0].hi, sun_dx) ;
	fprintf(stderr, "Grid optimize sun_py from %f to %f by %f\n", start_sun_py, opt_parms[1].hi, sun_dy) ;
	fprintf(stderr, "Grid optimize horizon_py from %f to %f by %f\n", start_horizon_py, opt_parms[2].hi, horizon_dy) ;

#define FTOL .0001
	for(sun_x=start_sun_x; sun_x < opt_parms[0].hi+FTOL ; sun_x += sun_dx) {
	    for(sun_py=start_sun_py; sun_py < opt_parms[1].hi+FTOL ; sun_py += sun_dy) {
		for(horizon_py=start_horizon_py; horizon_py < opt_parms[2].hi+FTOL ; horizon_py += horizon_dy) {
		    //fprintf(stderr, "grid evaluations %d, sx:%f sy:%f hy:%f\n", n_tot, sun_x, sun_py, horizon_py) ;
		    n_tot++ ;
		}
	    }
	}

	int i_grid_search_no=0 ;
	fprintf(stderr, "%d grid evaluations will be done\n", n_tot) ;

	for(sun_x=start_sun_x; sun_x < opt_parms[0].hi+FTOL ; sun_x += sun_dx) {
	    for(sun_py=start_sun_py; sun_py < opt_parms[1].hi+FTOL ; sun_py += sun_dy) {
		for(horizon_py=start_horizon_py; horizon_py < opt_parms[2].hi+FTOL ; horizon_py += horizon_dy) {

		    float err = optimize_grid(n_samples,0) ;
		    i_grid_search_no++ ;
		    printf("Completed %d of %d grid evaluations for sun position\n", i_grid_search_no, n_tot) ;

		    if(err < best_err) {
			best_err = err ;
			best_horizon_py = horizon_py ;
			best_s_x = sun_x ;
			best_s_y = sun_py ;
		    }

		}
	    }
	}

	/* reset for best values */
	sun_x = best_s_x ;
	sun_py = best_s_y ;
	horizon_py = best_horizon_py ;
	optimize_grid(n_samples,1) ;

	if(n_custom_optimizations == 0) {

    //-O sx hy B C sl sg -O sx sy sl sg A D E

	    char *opt0[6] ;
	    char *opt1[7] ;

	    char *opt0_a[6] = {"sx","hy","B","C","sl", "sd"} ;
	    char *opt1_a[7] = {"sx","sy","sl","A","D", "E", "sd"} ;

	    char *opt0_u[5] = {"sx","sy","hy","A","B"} ;
	    char *opt1_u[4] = {"sx","sy","sl","B"} ;
	    int i,j ;

	    int np0, np1 ;

	    if(allowed_sky_type == UNIFORM_CIE_SKIES) {
		np0 = 5 ;
		np1 = 4 ;
		for(j = 0 ; j < np0 ; j++) {
		    opt0[j] = opt0_u[j] ;
		}
		for(j = 0 ; j < np1 ; j++) {
		    opt1[j] = opt1_u[j] ;
		}
	    } else {
		np0 = 6 ;
		np1 = 7 ;
		for(j = 0 ; j < np0 ; j++) {
		    opt0[j] = opt0_a[j] ;
		}
		for(j = 0 ; j < np1 ; j++) {
		    opt1[j] = opt1_a[j] ;
		}
	    }


	    int n_vars=0 ;
	    optimization_var[n_custom_optimizations][0] = -1 ;

	    for(j = 0 ; j < np0 ; j++) {
		for(i = 0 ; i < MAX_PARMS ; i++) {
		    if(!strcmp(opt0[j], opt_parms[i].abreviation)) {
			if(opt_parms[i].optimize_flag == 1) {
			    optimization_var[n_custom_optimizations][n_vars] = i ;
			    n_vars++ ;
			    optimization_var[n_custom_optimizations][n_vars] = -1 ;
			}
		    }
		}
	    }
	    n_custom_optimizations++ ;

	    n_vars=0 ;
	    optimization_var[n_custom_optimizations][0] = -1 ;

	    for(j = 0 ; j < np1 ; j++) {
		for(i = 0 ; i < MAX_PARMS ; i++) {
		    if(!strcmp(opt1[j], opt_parms[i].abreviation)) {
			if(opt_parms[i].optimize_flag == 1) {
			    optimization_var[n_custom_optimizations][n_vars] = i ;
			    n_vars++ ;
			    optimization_var[n_custom_optimizations][n_vars] = -1 ;
			}
		    }
		}
	    }
	    n_custom_optimizations++ ;
	    optimization_var[n_custom_optimizations][0] = -1 ;
	}

    }

    {
	// perform custom optimizations either from defaults or specified on the command line
	int i,o ;

	for(o = 0 ; o < n_custom_optimizations ; o++) {
	    optimize_type = V_from_HSV ;

	    for(i = 0 ; i < MAX_PARMS ; i++) {
		opt_parms[i].optimize_flag=0 ;
	    }

	    for(i = 0 ; i < MAX_PARMS ; i++) {
		int j = optimization_var[o][i] ;
		if( j < 0) break ;
		opt_parms[j].optimize_flag=1 ;
		if(!strcmp(opt_parms[j].name, "sat_depth")) optimize_type = FULL_RGB ;
		if(!strcmp(opt_parms[j].name, "hue_sky")) optimize_type = FULL_RGB ;
		if(!strcmp(opt_parms[j].name, "hue_horizon")) optimize_type = FULL_RGB ;
	    }

	    smart_optimize(n_samples) ;

	}

	free(samples) ;

    }
}

void usage(char *msg, char *msg2)
{
    fprintf(stderr, "skyfill_tif: %s %s\n", msg, msg2) ;

    fprintf(stderr, "usage: skyfill_tif <inputfile> <outputfile> [options]\n") ;
    fprintf(stderr, "options are:\n") ;
    fprintf(stderr, "    -ds <DOS> -- sample sky values from sky start to (sky end-sky start)*DOS+sky_start in a column, DOS is a proportion, default 0.5\n") ;
    fprintf(stderr, "    -cdf <DOF> -- fill sky values from O to (min sky end- max sky start)*DOF+sky_start in all columns, DOF is a proportion, default 0.5\n") ;
    fprintf(stderr, "              useful when there is a lot of small peaks at the skyline which are visible in the sky merge.  Try if -df doesn't look good\n") ;
    fprintf(stderr, "    -df <DOF> -- fill sky values from O to (sky end-sky_start)*DOF+sky_start in a column, DOF is a proportion, default 0.5\n") ;
    fprintf(stderr, "    -ds <depth_of_sample> -- proportion, control how far in a column to find samples to compare to predicted sky\n") ;
    fprintf(stderr, "             -- affects final blending of predicted sky to original sky at that column\n") ;
    fprintf(stderr, "    -ff <FF> -- feather factor.  Default 1.0 Normally replaced sky is completely estimated at sky start, and\n") ;
    fprintf(stderr, "              uses 100%% of the orginal value at the final row (determined by DOF above\n") ;
    fprintf(stderr, "              setting FF to 0.0 would cause 100%% of the estimated value to be use on every\n") ;
    fprintf(stderr, "              pixel replaced\n") ;
    fprintf(stderr, "    -fse <0|1>  -- turn off|on algorithm to repair vertical edges with missing pixels near the start of the sky\n") ;
    fprintf(stderr, "    -fsg <gamma>  -- a typical image gamma function applied to predicted value of sky (v = v**(1./gamma))\n") ;
    fprintf(stderr, "    -fsh <left> <right> -- attempt to repair sky hue in a portion of the image, left and right are proportions\n") ;
    fprintf(stderr, "                       fsh 0.5 1.0 would attempt the repair on the right half of the image\n") ;
    fprintf(stderr, "    -ftos <0|1> -- brings top of sky up to a constant value before modeling sky values if set to 1\n") ;
    fprintf(stderr, "    -efl <number_pixels>  -- extend computed feather length of all columns by this many pixels \n") ;
    fprintf(stderr, "    -tol <th> <ts> <tv> -- tolerance for determining edge of skyline at horizon, default 360 .06 0.015\n") ;
    fprintf(stderr, "    -im -- interpolate sky in masked columns\n") ;
    fprintf(stderr, "    -hf -- run horizontal filter on sky area to smooth artifacts\n") ;


    fprintf(stderr, "    -fm <l:q><l:q><l:q> -- fit model type for hsv, l is linear, q is quadratic, default is qqq (quadratic fit for h,s,v components\n") ;
    fprintf(stderr, "    -m <l> <r> -- mask out columns from l-r from sky sampling or filling (up to 20 masks may be specified)\n") ;
    fprintf(stderr, "    -d -- debug, output image will have start and end of sky marked\n") ;
    fprintf(stderr, "    -G -- run the grid optimizer, using fixed parameters supplied on command line\n") ;
    fprintf(stderr, "    -rm <l> <r> -- mask out columns from l-r from having pixels replaced i.e. repaird because they appear to be outliers\n") ;
    fprintf(stderr, "    -r_tos_thresh <thresh> -- sets the tolerance for finding \"outlier\" pixels in a small section of the sky, if the pixel\n") ;
    fprintf(stderr, "                       is more than this fraction away from the mean, it is classified as an outlier.  Default thresh is 0.02\n") ;
    fprintf(stderr, "    -spm <0|1> -- saturation prediction method 0 will just use a mean value, 1 will predict saturation as a function of value\n") ;
    fprintf(stderr, "    -verbose -- adds a little information at the end of the run, really only useful for debuggin\n") ;
    fprintf(stderr, "    -t -- quick test most, scales the output image by 50%% to speed up processing\n") ;
    fprintf(stderr, "    -EO -- show estimated sky only, will replace the entire image with the optimized estimate of the sky\n") ;
    fprintf(stderr, "    -SR -- show estimated sky only for the sky area to be filled -- it will not have final blending applied\n") ;



    fprintf(stderr, "\n") ;
    fprintf(stderr, "    -O <one or more parms> optimize given parameters of sky dome, possible parameters are:\n") ;

    int i ;

    for(i = 0 ; i < MAX_PARMS ; i++) {
	fprintf(stderr, "               %s (%s) range is %lf to %lf --  %s\n",
	    opt_parms[i].abreviation,
	    opt_parms[i].name,
	    opt_parms[i].lo,
	    opt_parms[i].hi,
	    opt_parms[i].units) ;
    }
    fprintf(stderr, "\n               More than one -O may be specified, causing multiple optimizations to be performed\n") ;
    fprintf(stderr, "               sequentially with individual control over which parameters are optimized\n\n") ;

    fprintf(stderr, "    -any parameter may be specified by prepending with a '-', i.e. \"-sx 90\" sets the sun x position 90 degrees to the right\n") ;

    exit(1) ;
}

int main(int argc, char* argv[])
{
    float depth_of_sample = 0.5 ; // how far into the sky to sample for hue,val and sat
    float depth_of_fill = 0.5 ; // how far into existing sky to fill
    float feather_factor = 1.0 ; // how much feathering to do
    int extra_feather_length = 0 ; // extend feathering length by this many pixels
    int run_hf=0 ; // run horizontal filter
    int debug=0 ;

    int n_masks=0 ;
    int mask_l[20] ;
    int mask_r[20] ;

    int n_repair_masks=0 ;
    int repair_mask_l[20] ;
    int repair_mask_r[20] ;

    char infilename[512] ;
    char outfilename[512] ;
    char ptofilename[512] ;
    int interpolate_masked_columns=0 ;
    int verbose_flag=0 ;
    int use_exiftool=1 ;
    int make_sky_uniform_flag=-1 ;
    int quick_test=0 ;
    float fmin_sky_hue_mask = 2. ;
    float fmax_sky_hue_mask = 2. ;


    if(argc < 3)
	usage("Input or output file not specified", "") ;

    // strip final \n off the last input argument
    {
	int l = strlen(argv[argc-1]) ;

	if(argv[argc-1][l-1] == '\n')
	    argv[argc-1][l-1] = '\0' ;
    }

    // make sure input tif is not the same as output tif
    if(! strcmp(argv[1], argv[2]) ) {
	usage("Input file name must not be same as output filename", "") ;
    }


    strcpy(infilename, argv[1]) ;
    strcpy(outfilename, argv[2]) ;
    strcpy(ptofilename, argv[1]) ;
    int pto_fov=-1 ;

    // can we open a corresponding hugin pto file, and get the horizontal FOV from it?
    {
	int i ;

	for(i=0 ; i < strlen(ptofilename) ; i++) {
	    if(ptofilename[i] == '.')
		break ;
	}

	ptofilename[i+1] = 'p' ;
	ptofilename[i+2] = 't' ;
	ptofilename[i+3] = 'o' ;
	ptofilename[i+4] = '\0' ;
	printf("ptofilename is %s\n", ptofilename) ;

	FILE *fp = fopen(ptofilename,"r") ;

	if(fp != NULL) {
	    char s[1000] ;
	    while(fgets(s, 999, fp)) {
		int j ;

		if(s[0] != 'p')
		    continue ;

		for(j=0 ; j < strlen(s) ; j++) {
		    if(s[j] == 'v')
			break ;
		}

		if(s[j] == 'v' && s[j-1] == ' ') {
		    sscanf(&s[j+1], "%d", &pto_fov) ;
		    printf("PTO fov is %d\n", pto_fov) ;
		}
	    }

	    fclose(fp) ;
	}
    }


    TIFF* tif = TIFFOpen(infilename, "r");
    argc -= 3 ;
    argv += 3 ;
    int n_custom_optimizations=0 ;
    int force_grid_optimization = 0 ;

    int optimization_var[10][MAX_PARMS+1] ;

    TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &IMAGE_NSAMPLES);

    if(IMAGE_NSAMPLES < 4) {
	fprintf(stderr, "FATAL:  skyfill only works with 4 channel (rgba) images produced by hugin\n") ;
	exit(1) ;
    }

    uint32 input_W, input_H ;
    uint16 BPS, SPP, TIF_CONFIG ;

    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &input_W);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &input_H);
    TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &BPS);
    TIFFGetField(tif, TIFFTAG_PLANARCONFIG, &TIF_CONFIG);

    if(BPS != 16) {
	fprintf(stderr, "FATAL:  skyfill only works with 16 bit tiff images\n") ;
	exit(1) ;
    }


    while(argc > 0) {
	int found_set_parm=0 ;
	int i ;

	for(i = 0 ; i < MAX_PARMS ; i++) {
	    char set_string[20]  ;
	    snprintf(set_string, 19, "-%s", opt_parms[i].abreviation) ;

	    if(!strcmp(argv[0], set_string)) {
		*opt_parms[i].variable_address = atof(argv[1]) ;
		opt_parms[i].grid_optimize_flag = 0 ;
		opt_parms[i].optimize_flag = 0 ;
		printf("setting %s to %f\n", opt_parms[i].name, *opt_parms[i].variable_address) ;
		found_set_parm=1 ;
		argv++ ;
		argc-- ;
		argv++ ;
		argc-- ;
		break ;
	    }
	}

	if(found_set_parm) continue ;


	if(!strcmp(argv[0], "-O")) {
	    // optimization requested, which variables?
	    argc-- ;
	    argv++ ;
	    int n_vars=0 ;

	    while(argc > 0) {
		if(argv[0][0] != '-') {
		    for(i = 0 ; i < MAX_PARMS ; i++) {
			if(!strcmp(argv[0], opt_parms[i].abreviation)) {
			    optimization_var[n_custom_optimizations][n_vars] = i ;
			    n_vars++ ;
			    optimization_var[n_custom_optimizations][n_vars] = -1 ;
			    printf("optimization %d, optimize %s\n", n_custom_optimizations, opt_parms[i].name) ;
			    argc-- ;
			    argv++ ;
			    break ;
			}
		    }
		} else {
		    break ;
		}
	    }
	    n_custom_optimizations++ ;
	    continue ;
	}

	if(!strcmp(argv[0], "-CIE")) {
	    CIE_sky_index = atoi(argv[1]) ;
	    set_cie_parms(CIE_sky_index-1) ;
	    argc -= 2 ;
	    argv += 2 ;
	    continue ;
	}

	if(!strcmp(argv[0], "-CIE_types")) {

	    if(argv[1][0] == 'a')
		allowed_sky_type = ALL_CIE_SKIES ;
	    else if(argv[1][0] == 'n')
		allowed_sky_type = NONUNIFORM_CIE_SKIES ;
	    else if(argv[1][0] == 'u')
		allowed_sky_type = UNIFORM_CIE_SKIES ;
	    else {
		usage("unknown sky type for -CIE_types, allowed are 'a', 'n' and 'u'", "") ;
	    }

	    argc -= 2 ;
	    argv += 2 ;
	    continue ;
	}

	if(!strcmp(argv[0], "-d")) {
	    if(debug < 1) debug=1 ;
	    argc -= 1 ;
	    argv += 1 ;
	    continue ;
	}

	if(!strcmp(argv[0], "-t")) {
	    quick_test=1 ;
	    argc -= 1 ;
	    argv += 1 ;
	    continue ;
	}

	if(strlen(argv[0]) == 3) {
	    if(argv[0][0] == '-' && argv[0][1] == 'd') {
		if(argv[0][2] >= '0' && argv[0][2] <= '9') {
		    debug = atoi(&(argv[0][2])) ;
		    fprintf(stderr, "debug set to %d\n", debug) ;
		    argc -= 1 ;
		    argv += 1 ;
		    continue ;
		}
	    }
	}

	if(!strcmp(argv[0], "-cdf")) {
	    depth_of_fill = atof(argv[1]) ;
	    depth_of_fill_is_absolute = 1 ;
	    argc -= 2 ;
	    argv += 2 ;
	    continue ;
	}

	if(!strcmp(argv[0], "-df")) {
	    depth_of_fill = atof(argv[1]) ;
	    argc -= 2 ;
	    argv += 2 ;
	    continue ;
	}

	if(!strcmp(argv[0], "-ds")) {
	    depth_of_sample = atof(argv[1]) ;
	    argc -= 2 ;
	    argv += 2 ;
	    continue ;
	}

	if(!strcmp(argv[0], "-efl")) {
	    extra_feather_length = atoi(argv[1]) ;
	    argc -= 2 ;
	    argv += 2 ;
	    continue ;
	}

	if(!strcmp(argv[0], "-ff")) {
	    feather_factor = atof(argv[1]) ;
	    argc -= 2 ;
	    argv += 2 ;
	    continue ;
	}

	if(!strcmp(argv[0], "-fs")) {
	    final_saturation_factor = atof(argv[1]) ;
	    argc -= 2 ;
	    argv += 2 ;
	    continue ;
	}

	if(!strcmp(argv[0], "-fse")) {
	    fix_SOS_edges = atoi(argv[1]) ;
	    argc -= 2 ;
	    argv += 2 ;
	    continue ;
	}

	if(!strcmp(argv[0], "-fsg")) {
	    final_sky_gamma = atof(argv[1]) ;
	    argc -= 2 ;
	    argv += 2 ;
	    continue ;
	}

	if(!strcmp(argv[0], "-fsh")) {
	    fix_sky_hue=1 ;
	    fmin_sky_hue_mask = atof(argv[1]) ;
	    fmax_sky_hue_mask = atof(argv[2]) ;
	    argc -= 3 ;
	    argv += 3 ;
	    continue ;
	}

	if(!strcmp(argv[0], "-ftos")) {
	    fill_top_of_sky = atoi(argv[1]) ;
	    argc -= 2 ;
	    argv += 2 ;
	    continue ;
	}

	if(!strcmp(argv[0], "-hf")) {
	    run_hf=1 ;
	    argc -= 1 ;
	    argv += 1 ;
	    continue ;
	}
	if(!strcmp(argv[0], "-im")) {
	    interpolate_masked_columns = 1 ;
	    argc -= 1 ;
	    argv += 1 ;
	    continue ;
	}

	if(!strcmp(argv[0], "-m")) {
	    mask_l[n_masks] = atoi(argv[1]) ;
	    mask_r[n_masks] = atoi(argv[2]) ;
	    n_masks++ ;
	    argc -= 3 ;
	    argv += 3 ;
	    continue ;
	}

	if(!strcmp(argv[0], "-mse")) {
	    min_sky_end_p = atof(argv[1]) ;
	    argc -= 2 ;
	    argv += 2 ;
	    continue ;
	}
	if(!strcmp(argv[0], "-msu")) {
	    make_sky_uniform_flag = atoi(argv[1]) ;
	    argc -= 2 ;
	    argv += 2 ;
	    continue ;
	}

	if(!strcmp(argv[0], "-nf")) {
	    nonlinear_feather = atof(argv[1]) ;
	    argc -= 2 ;
	    argv += 2 ;
	    continue ;
	}

	if(!strcmp(argv[0], "-r_tos_thresh")) {
	    repair_tos_thresh = atof(argv[1]) ;
	    argc -= 2 ;
	    argv += 2 ;
	    continue ;
	}

	if(!strcmp(argv[0], "-rm")) {
	    repair_mask_l[n_repair_masks] = atoi(argv[1]) ;
	    repair_mask_r[n_repair_masks] = atoi(argv[2]) ;
	    n_repair_masks++ ;
	    argc -= 3 ;
	    argv += 3 ;
	    continue ;
	}

	if(!strcmp(argv[0], "-spm")) {
	    sat_prediction_method = atoi(argv[1]) ;
	    argc -= 2 ;
	    argv += 2 ;
	    continue ;
	}

	if(!strcmp(argv[0], "-tol")) {
	    sky_hue_tolerance = atof(argv[1]) ;
	    sky_sat_tolerance = atof(argv[2]) ;
	    sky_val_tolerance = atof(argv[3]) ;
	    argc -= 4 ;
	    argv += 4 ;
	    continue ;
	}
	if(!strcmp(argv[0], "-verbose")) {
	    verbose_flag = atoi(argv[1]) ;
	    argc -= 2 ;
	    argv += 2 ;
	    continue ;
	}

	if(!strcmp(argv[0], "-EO")) {
	    estimate_only=1 ;
	    argc -= 1 ;
	    argv += 1 ;
	    continue ;
	}
	if(!strcmp(argv[0], "-G")) {
	    force_grid_optimization=1 ;
	    argc -= 1 ;
	    argv += 1 ;
	    continue ;
	}
	if(!strcmp(argv[0], "-SR")) {
	    show_raw_prediction=1 ;
	    argc -= 1 ;
	    argv += 1 ;
	    continue ;
	}
	usage("Unknown option", argv[0]) ;
    }

    if(pto_fov > 1) {
	FOV_horizontal = (float)pto_fov ;
	have_pto_fov=1 ;
    }
    if (!tif) {
	fprintf(stderr, "failed to open %s\n", infilename) ;
	exit(1) ;
    }


    TIFF* tifout = TIFFOpen(outfilename, "w");

    tsize_t ROWSIZE ;
    int32 x, y, w, h ;

    // prepare output image to have same tags as input image
    uint16 alltags[] = {
			TIFFTAG_IMAGEWIDTH,
			TIFFTAG_IMAGELENGTH,
			TIFFTAG_BITSPERSAMPLE,
			TIFFTAG_SAMPLESPERPIXEL,
			TIFFTAG_ROWSPERSTRIP,
			TIFFTAG_ORIENTATION,
			TIFFTAG_PLANARCONFIG,
			TIFFTAG_PHOTOMETRIC,
			TIFFTAG_SAMPLEFORMAT,
			TIFFTAG_COMPRESSION,
			} ;
    int i ;

    for(i = 0 ; i < 10 ; i++) {
	uint16 tag = alltags[i] ;
	uint32 tmp ;

	if(tag == TIFFTAG_IMAGEWIDTH && quick_test) {
	    if(TIFFGetField(tif, tag, &tmp)) {
		tmp /= 2 ;
		TIFFSetField(tifout, tag, tmp);
	    }
	} else if(tag == TIFFTAG_IMAGELENGTH && quick_test) {
	    if(TIFFGetField(tif, tag, &tmp)) {
		tmp /= 2 ;
		TIFFSetField(tifout, tag, tmp);
	    }
	} else {
	    if(TIFFGetField(tif, tag, &tmp)) {
		fprintf(stderr, "setting tag number %d in output file\n", tag) ;
		TIFFSetField(tifout, tag, tmp);
	    }
	}
    }

    if(!use_exiftool) {
	{
	    uint16 tag = TIFFTAG_ICCPROFILE ;

	    uint32 tmp[4] ;
	    unsigned char buf[4][2000] ;
	    if(TIFFGetField(tif, tag, tmp, buf)) {
		fprintf(stderr, "setting ICCPROFILE tag in output file, tmp is %d,%d,%d\n", (int)tmp[0],(int)tmp[1],(int)tmp[2]) ;
		TIFFSetField(tifout, tag, tmp[0], buf);
	    }
	    //exit(1) ;
	}
    }


    tdata_t *image ;


    ROWSIZE = TIFFScanlineSize(tif) ;
    SPP = IMAGE_NSAMPLES ;

    IMAGE_HAS_ALPHA = IMAGE_NSAMPLES > 3 ? 1 : 0 ;


    // make signed int values for width and height
    h = (int32)input_H ;
    w = (int32)input_W ;

    if(quick_test) {
	h /= 2 ;
	w /= 2 ;
    }

    min_sky_hue_mask = fmin_sky_hue_mask*w ;
    max_sky_hue_mask = fmax_sky_hue_mask*w ;

    IMAGE_HEIGHT = (int32)h ; // set global variable
    IMAGE_WIDTH = (int32)w ; // set global variable

    p_half_image_width = 0.5 ;
    V_zenith= P3toV3(0., 1., 0.) ;

    // if sun x position was supplied via the command line, it was given as a proportion of the image width
    // now it must be converted to the uniform coordinates

#ifdef OLD_SUN_X
#define XP2U(x) ((float)(x)*((float)IMAGE_WIDTH)/((float)IMAGE_HEIGHT))

    if(sun_x_pwidth > -9999.) {
	sun_x = XP2U(sun_x_pwidth) ;
    }

    if(sun_x_pwidth_hi > -9999.) {
	sun_x_hi = XP2U(sun_x_pwidth_hi) ;
    }

    if(sun_x_pwidth_lo > -9999.) {
	sun_x_lo = XP2U(sun_x_pwidth_lo) ;
    }

    opt_parms[0].lo = sun_x_lo ;
    opt_parms[0].hi = sun_x_hi ;
#endif

    fprintf(stderr, "input image is %d x %d, %d bits, %d channels, %d = %d\n", (int)input_W, (int)input_H, (int)BPS, (int)SPP, (int)(input_W*BPS/8*SPP), (int)ROWSIZE) ;

    // apply column masks to column_mask array

    column_mask = (unsigned char *)calloc(w, sizeof(unsigned char)) ;
    for(x = 0 ; x < w ; x++) {
	column_mask[x] = 0 ;
    }
    for(int m = 0 ; m < n_masks ; m++) {
	for(x = mask_l[m] ; x <= mask_r[m] ; x++) {
	    if(x < 0) continue ;
	    if(x > input_W-1) continue ;
	    column_mask[x/(quick_test+1)] = 1 ;
	}
    }

    int x0 ;

    for(x0 = 0 ; x0 < w ; x0++) {
	if(column_mask[x0] == 0) break ;
    }

    if(x0 >= w) {
	fprintf(stderr, "All columns are masked!, no output generated\n") ;
	exit(1) ;
    }

    // apply column repair masks to column_repair_mask array

    column_repair_mask = (unsigned char *)calloc(w, sizeof(unsigned char)) ;
    for(x = 0 ; x < w ; x++) {
	column_repair_mask[x/(quick_test+1)] = 0 ;
    }
    for(int m = 0 ; m < n_repair_masks ; m++) {
	for(x = repair_mask_l[m] ; x <= repair_mask_r[m] ; x++) {
	    if(x < 0) continue ;
	    if(x > input_W-1) continue ;
	    column_repair_mask[x/(quick_test+1)] = 1 ;
	}
    }

    // Read the input tif image

    if((int)ROWSIZE != (int)(input_W*BPS/8*SPP)) {
	fprintf(stderr, "ROWSIZE is WRONG!, fixing...\n") ;
	ROWSIZE = input_W*BPS/8*SPP ;
	fprintf(stderr, "ROWSIZE is now %d\n", (int)ROWSIZE) ;
    }

    /* allocate memory for image */
    image = (tdata_t *)calloc(input_H, sizeof(tdata_t *)) ;

    for(y = 0 ; y < input_H ; y++) {
	image[y] = (tdata_t) _TIFFmalloc(ROWSIZE);
	if(image[y] == NULL) {
	    fprintf(stderr, "Failed to alloc memory for image at row %d\n", y) ;
	    exit(1) ;
	}
    }

    read_tif_image(image,tif,input_H,SPP, TIF_CONFIG) ;


    if(quick_test) {
	//scale image by 1/2
	fprintf(stderr, "Scale image (%d,%d) by half using nearest pixel method to (%d,%d)\n", input_W,input_H,IMAGE_WIDTH,IMAGE_HEIGHT) ;

	for(x = 0 ; x < IMAGE_WIDTH  ; x++) {
	    int x2 = x*2 ;
	    if(x2 > input_W) {
		fprintf(stderr, "uh-oh, x:%d * 2 (%d) > input width %d\n", x, x*2, input_W) ;
		exit(1) ;
	    }


	    for(y = 0 ; y < IMAGE_HEIGHT ; y++) {
		int y2 = y*2 ;
		if(y2 > input_H) {
		    fprintf(stderr, "uh-oh, y:%d * 2 (%d) > input height %d\n", y, y*2, input_H) ;
		    exit(1) ;
		}
/*  		fprintf(stderr, "y:%d * 2 (%d) > input height %d\n", y, y*2, input_H) ;  */
		uint r,g,b,a ;
		tif_get4c(image,x2,y2,r,g,b,a) ;
		tif_set4c(image,x,y,r,g,b,a) ;
	    }
	}

    }

    // input image has read, and possibly scaled by 1/2, ready for image processing


    fprintf(stderr, "repair alpha\n") ;
    repair_alpha(image) ;

    /* find start of real image data */
    raw_start_of_sky = (int16 *)calloc(w, sizeof(int16 *)) ;
    start_of_sky = (int16 *)calloc(w, sizeof(int16 *)) ;
    end_of_sky = (int16 *)calloc(w, sizeof(int16 *)) ;
    final_end_of_sky = (int16 *)calloc(w, sizeof(int16 *)) ;


    simple_find_start_of_sky(start_of_sky,w,h,image) ;

    for(x = 0 ; x < w ; x++) {
	raw_start_of_sky[x] = start_of_sky[x] ;
    }

    find_start_of_sky(start_of_sky,w,h,image,fix_SOS_edges) ;

    int fix_slivers=1 ;
    int fix_edges=fix_SOS_edges ;
    find_end_of_sky(end_of_sky,start_of_sky,w,h,image,fix_edges,fix_slivers) ;

    repair_top_of_sky(image, 1, 1) ;  // 1-> end of sky is known, 1 means phase 1
    if(debug == 1) goto writeout ;

    find_start_of_sky(start_of_sky,w,h,image,fix_SOS_edges) ;
    find_end_of_sky(end_of_sky,start_of_sky,w,h,image,fix_edges,fix_slivers) ;
    if(debug == 2) goto writeout ;

    set_minmax_sky_values() ;

    if(depth_of_fill_is_absolute) {
	float length = min_end_of_sky - max_start_of_sky ;
	depth_of_fill_absolute_y = (int)((float)max_start_of_sky + depth_of_fill*length+0.5) ;
    }

    set_FOV_factor() ;
    set_xy_constraints() ;

    int n_samples = sample_sky_points(10, 200,image,start_of_sky,end_of_sky) ;
    fprintf(stderr, "sampled %d points in image\n", n_samples) ;

    if(fix_sky_hue) {
	// mainly corrects for blown out regions where hue is messed up,
	// ** Note, the starting hue_sky is used as estimate of actual hue
	repair_sky_hue(image,start_of_sky,end_of_sky) ;
/*  	fix_sky_hue=0 ;  */
/*  	n_samples = sample_sky_points(10, 200,image,start_of_sky,end_of_sky) ;  */
    }

    repair_top_of_sky(image, 1, 3) ;  // 1 means end of sky is known
    if(debug == 3) goto writeout ;
    repair_top_of_sky(image, 1, 4) ;  // 1 means end of sky is known
    if(debug == 4) goto writeout ;
    repair_top_of_sky(image, 1, 5) ;  // 1 means end of sky is known
    if(debug == 5) goto writeout ;
    repair_top_of_sky(image, 1, 6) ;  // 1 means end of sky is known
    if(debug == 6) goto writeout ;
    if(debug == 9) goto writeout ;

    if(make_sky_uniform_flag == 2) make_sky_uniform(image,make_sky_uniform_flag,0,IMAGE_WIDTH-1) ;

    if(debug == 7) goto writeout ;

    fit_sky_model(force_grid_optimization, n_custom_optimizations, optimization_var, n_samples) ;

estimate_sky:
    maximum_vhat = find_maximum_vhat() ;
    model_is_being_fit=0 ;

    estimate_sky(0,IMAGE_WIDTH-1,image,start_of_sky,end_of_sky,final_end_of_sky,depth_of_sample,h,w,feather_factor,debug,depth_of_fill,extra_feather_length) ;

    fprintf(stderr, "back in main\n") ;

    int min_sky=h ;

    if(interpolate_masked_columns || run_hf) {
	// find minimum end of sky for all x
	for(x = 0 ; x < w ; x++) {
	    if(column_mask[x] == 1) continue ;
	    if(min_sky > end_of_sky[x]) min_sky = end_of_sky[x] ;
	}
	fprintf(stderr, "min sky is %d\n", min_sky) ;
    }

    if(interpolate_masked_columns) {
	fprintf(stderr, "main:interpolate masked columns\n") ;
	for(int m = 0 ; m < n_masks ; m++) {
	    if(mask_l[m] > 0 && mask_r[m] < w-1) {
		int32 x0 = mask_l[m]-1 ;
		int32 x1 = mask_r[m]+1 ;

		int max_y = start_of_sky[x0] ;

		for(x = x0+1 ; x < x1 ; x++) {
		    if(max_y < start_of_sky[x]) max_y = start_of_sky[x] ;
		}

		max_y += 10 ;

		for(y = 0 ; y < max_y ; y++) {
		    /* interpolate value from start to end */
		    float start_r, start_g, start_b ;
		    float   end_r, end_g, end_b ;

		    tif_get3c(image,x0,y,start_r,start_g,start_b) ;
		    tif_get3c(image,x1,y,end_r,end_g,end_b) ;

		    for(x = x0+1 ; x < x1 ; x++) {
			float pStart = 1.f - ((float)x-((float)x0))/((float)x1-(float)x0+1.f) ;
			float pEnd = 1.f - pStart ;
			tif_set4c(image,x,y, (uint16)(pStart*start_r + pEnd*end_r), (uint16)(pStart*start_g + pEnd*end_g), (uint16)(pStart*start_b + pEnd*end_b), MAX16) ;
		    }
		}

	    }
	}
    }

    if(make_sky_uniform_flag == 0 || make_sky_uniform_flag == 1) make_sky_uniform(image,make_sky_uniform_flag, 0,IMAGE_WIDTH-1) ;

    if(run_hf) {
	// run a Horizontal filter on the sky to smooth edges
	fprintf(stderr, "main:horizontal filter\n") ;


	// alloc a buffer to hold output data
#define HBW 50
	float *buf = (float *) calloc(w, sizeof(float)) ;

	for(int ch=0 ; ch < 3 ; ch++) {
	    for(y = 0 ; y < min_sky ; y++) {
		// read values into buf
		for(x = 0 ; x < w ; x++) {
		    buf[x] = (float)((uint16 *)(image[y]))[IMAGE_NSAMPLES*x+ch] ;
		}

		for(x=HBW ; x < w-HBW ; x++) {
		    float sum=0. ;
		    float n =0.;
		    int dx ;

		    for(dx=x-HBW ; dx < x+HBW ; dx++) {
			if(column_mask[dx] == 1) continue ;
			sum += buf[dx] ;
			n += 1. ;
		    }
		    ((uint16 *)(image[y]))[IMAGE_NSAMPLES*x+ch] = (uint16)(sum/n) ;

		    if(x == HBW) {
			for(dx=0 ; dx < HBW ; dx++) {
			    if(column_mask[dx] == 1) continue ;
			    ((uint16 *)(image[y]))[IMAGE_NSAMPLES*dx+ch] = (uint16)(sum/(float)(HBW*2+1)) ;
			}
		    } else if(x == w-HBW-1) {
			for(dx=x ; dx < x+HBW ; dx++) {
			    if(column_mask[dx] == 1) continue ;
			    ((uint16 *)(image[y]))[IMAGE_NSAMPLES*dx+ch] = (uint16)(sum/(float)(HBW*2+1)) ;
			}
		    } else {
			if(column_mask[x] == 1) continue ;
			((uint16 *)(image[y]))[IMAGE_NSAMPLES*x+ch] = (uint16)(sum/(float)(HBW*2+1)) ;
		    }

		}

	    }
	}

	free(buf) ;
    }



#define DBW 4
#define DBH 6
    int df_passes=2 ;
    {
	// this needs to be in {}'s because the "goto writeout" isn't happy when
	// other variables are declared in the scope between the goto and the label
	uint16 *output_buf[IMAGE_HEIGHT] ;
	float save_depth_of_fill=depth_of_fill ;
	depth_of_fill = 1. ;
	//extra_feather_length=5 ;

	// run a Disc filter on the sky to smooth edges
	if(df_passes > 0) {

	    fprintf(stderr, "main:disc filter size %d X %d\n", DBW, DBH) ;
	    int i ;

	    // alloc a buffer to hold output data
	    for(i = 0 ; i < IMAGE_HEIGHT ; i++) {
		output_buf[i] = (uint16 *) calloc(IMAGE_WIDTH, sizeof(uint16)) ;
	    }
	}

	while(df_passes > 0) {
	    int i ;
	    df_passes-- ;

	    fprintf(stderr, "main:disc filter\n") ;

	    for(int ch=0 ; ch < 3 ; ch++) {
		int ox, oy ;

		for(ox = 0 ; ox < IMAGE_WIDTH ; ox++) {
		    if(column_mask[ox] == 1) continue ;
		    int feather_end_y ;
		    int feather_length = raw_compute_feather_length(ox, &feather_end_y, depth_of_fill, feather_factor, extra_feather_length) ;
		    int maxy = start_of_sky[ox]+feather_length ;

		    if(maxy < raw_start_of_sky[ox]) maxy = raw_start_of_sky[ox] ;
		    if(maxy > IMAGE_HEIGHT-1) maxy = IMAGE_HEIGHT-1 ;

		    if(ox == 0 && ch == 0) fprintf(stderr, "main:disc filter, start y loop\n") ;

		    for(oy = 0 ; oy < maxy ; oy++) {

			float sum=0. ;
			float n =0.;

			for(x=ox-DBW ; x <= ox+DBW ; x++) {
			    if(column_mask[x] == 1) continue ;
			    if(x < 0) continue ;
			    if(x > IMAGE_WIDTH-1) continue ;

			    for(y=oy-DBH ; y <= oy+DBH ; y++) {
				if(y < 0) continue ;
				if(y > maxy) break ;
				if(y > IMAGE_HEIGHT-1) break ;
				sum += ((uint16 *)(image[y]))[IMAGE_NSAMPLES*x+ch] ;
				n += 1. ;
			    }
			}

			output_buf[oy][ox] = (uint16)(sum/n+0.5) ;
		    }
		}

		for(ox = 0 ; ox < IMAGE_WIDTH ; ox++) {
		    if(column_mask[ox] == 1) continue ;
		    int feather_end_y ;
		    int feather_length = raw_compute_feather_length(ox, &feather_end_y, depth_of_fill, feather_factor, extra_feather_length) ;
		    int maxy = start_of_sky[ox]+feather_length ;

		    if(maxy < raw_start_of_sky[ox]) maxy = raw_start_of_sky[ox] ;
		    if(maxy > IMAGE_HEIGHT-1) maxy = IMAGE_HEIGHT-1 ;
		    if(ox == 0 && ch == 0) fprintf(stderr, "main:disc filter, start second y loop\n") ;

		    for(oy = 0 ; oy < maxy ; oy++) {

			uint16 b = ((uint16 *)(image[oy]))[IMAGE_NSAMPLES*ox+2] ;
			float fp0 = raw_compute_feather_at_y(ox, oy, feather_length, b, feather_factor) ;
			uint16 new = output_buf[oy][ox] ;
			uint16 old = ((uint16 *)(image[oy]))[IMAGE_NSAMPLES*ox+ch] ;
			((uint16 *)(image[oy]))[IMAGE_NSAMPLES*ox+ch] = (uint16)(fp0*(float)new + (1.-fp0)*(float)old + 0.5) ;

/*  			uint16 old = ((uint16 *)(image[oy]))[IMAGE_NSAMPLES*ox+ch] = output_buf[oy][ox] ;  */
		    }
		}
	    }


	    if(df_passes == 0) {
		for(i = 0 ; i <= max_end_of_sky ; i++) {
		    free(output_buf[i]) ;
		}
	    }

	}

	depth_of_fill = save_depth_of_fill ;
    }


    int dither_sky=1 ;
    if(dither_sky) {
	srand(0) ;
	/* dither sky */
	for(x = 0 ; x < IMAGE_WIDTH  ; x++) {
	    if(column_mask[x] == 1) continue ;
	    for(y = 0 ; y < raw_start_of_sky[x] ; y++) {
		int16 dither = ((rand()%128) - 64)*8 ;
		((uint16 *)(image[y]))[IMAGE_NSAMPLES*x+0] += dither ;
		dither = ((rand()%128) - 64)*8 ;
		((uint16 *)(image[y]))[IMAGE_NSAMPLES*x+1] += dither ;
		dither = ((rand()%128) - 64)*8 ;
		((uint16 *)(image[y]))[IMAGE_NSAMPLES*x+2] += dither ;
	    }
	}
    }

writeout:
    if(debug == 9) {
	// trying out new end of sky detection
	float *obuf ;
	obuf = (float *) calloc(IMAGE_HEIGHT, sizeof(float)) ;
	for(x = 0 ; x < IMAGE_WIDTH  ; x++) {
	    if(column_mask[x] == 1) continue ;

	    int width=3 ;
	    for(y = IMAGE_HEIGHT-width-1 ; y > width ; y--) {
		float r0,g0,b0,a0 ;
		float r1,g1,b1,a1 ;
		r0=g0=b0=a0=0. ;
		r1=g1=b1=a1=0. ;
		float r,g,b ;

		for(int z = y-width ; z < y ; z++) {
		    tif_get3c(image,x,z,r,g,b) ;
		    r0 += r ;
		    g0 += g ;
		    b0 += b ;
		}

		for(int z = y+1 ; z <= y+width ; z++) {
		    tif_get3c(image,x,z,r,g,b) ;
		    r1 += r ;
		    g1 += g ;
		    b1 += b ;
		}

		r0 /= (float)width ;
		r1 /= (float)width ;
		g0 /= (float)width ;
		g1 /= (float)width ;
		b0 /= (float)width ;
		b1 /= (float)width ;

		float gamma = 1./2.2 ;
		gamma = 1. ;
		float dr = powf(fabs(r1-r0),1./gamma) ;
		float dg = powf(fabs(g1-g0),1./gamma) ;
		float db = powf(fabs(b1-b0),1./gamma) ;
/*  		float da = powf(fabs(a1-a0),1./gamma) ;  */

		dr = dr + dg + db ;
		dg = dr ;
		db = dr ;
/*  		dr = (g1+1.)/(r1+1.) * (float)MAX16 ;  */
/*  		dg = dr ;  */
/*  		db = dr ;  */
/*  		tif_set4c(image,x,y,dr,dg,db,MAX16) ;  */
		obuf[y] = dr ;
	    }

	    int sos = 0 ;
	    for(y = width+1 ; y < IMAGE_HEIGHT-width-1 ; y++) {
		if(obuf[y] > 8000) {
		    sos = y+width ;
		    break ;
		}
	    }

	    int eos = 0 ;
	    float val_tol = sky_val_tolerance/.015* 4000. ;
	    for(y = sos+width*2 ; y < IMAGE_HEIGHT-width-1 ; y++) {
		if(obuf[y] > val_tol) {
		    eos = y+width-1 ;
		    break ;
		}
	    }

	    tif_set4c(image,x,sos-1,MAX16,0,0,MAX16) ;
	    tif_set4c(image,x,eos+2,0,MAX16,0,MAX16) ;

/*  	    for(y = IMAGE_HEIGHT-width-1 ; y > width ; y--) {  */
/*  		tif_set4c(image,x,y,obuf[y],obuf[y],obuf[y],MAX16) ;  */
/*  	    }  */
	}
	free(obuf) ;
    } else if(debug) {
	for(x = 0 ; x < IMAGE_WIDTH  ; x++) {
	    if(column_mask[x] == 1) continue ;

	    float sky_height = end_of_sky[x] - start_of_sky[x] ;
	    int feather_end_y = start_of_sky[x] + (int)(sky_height*depth_of_fill+0.5) ;

	    // indicate start of sky as red line
	    for(y = start_of_sky[x] ; y < start_of_sky[x]+1 ; y++) {
		if(y < 0 || y >=h) {
		    fprintf(stderr, "debug: sos y=%d\n", y) ;
		    exit(1) ;
		}
		tif_set4c(image,x,y,MAX16,0,0,MAX16) ;
	    }

	    // indicate raw_start of sky as yellow line
	    for(y = raw_start_of_sky[x] ; y < raw_start_of_sky[x]+1 ; y++) {
		if(y < 0 || y >=h) {
		    fprintf(stderr, "debug: sos y=%d\n", y) ;
		    exit(1) ;
		}

		tif_set4c(image,x,y,HALF16,HALF16,0,MAX16) ;
	    }

	    // indicate end of fill as black
	    for(y = feather_end_y ; y < feather_end_y+1 ; y++) {
		if(y < 0 || y >=h) {
		    fprintf(stderr, "debug: feather y=%d\n", y) ;
		    exit(1) ;
		}

		tif_set4c(image,x,y,0,0,0,MAX16) ;
	    }

	    // indicate end of detected sky as green
	    //for(y = end_of_sky[x]-1 ; y < end_of_sky[x]+1 ; y++) {
	    for(y = end_of_sky[x] ; y < end_of_sky[x]+1 ; y++) {
		if(y < 0 || y >=h) {
		    fprintf(stderr, "debug: eos x=%d, y=%d\n", x, y) ;
		    exit(1) ;
		}

		int yset = y+2 ; // set the color of the pixel two below EOS to green, so you can see the transition pixels above

		if(yset > h-1) yset = h-1 ;

		tif_set4c(image,x,yset,0,MAX16,0,MAX16) ;
	    }
	}
    }



    fprintf(stderr, "main:save\n") ;
    fprintf(stderr, "sky detect tolerances (hsv) are %lf, %lf, %lf\n", sky_hue_tolerance, sky_sat_tolerance, sky_val_tolerance) ;
    fprintf(stderr, "Value correction factor is %lf\n", sample_based_v_correction) ;

    verbose=verbose_flag ;

    if(verbose) {
	//some useful debugging information
	{
	    float sun_px = sun_x_angle2px(sun_x) ;
	    float gamma,theta ;
	    printf("horizon_py:%f, sun_x_angle:%f sun(px,py):(%6.2f,%6.2f)\n\n", horizon_py, sun_x, sun_px, sun_py) ;
	    printf("======== image_relative pixel_to_V3 for sun ============\n") ;
	    F_CIE2003(sun_px, sun_py, &gamma, &theta) ;
	    printf("======== image_relative pixel_to_V3 for sun ============\n\n") ;
	}

	float px, py ;
	py = sun_py ;
	for(px = -0.5 ; px < 0.51 ; px += .05) {
	    float gamma,theta ;
	    float vhat = F_CIE2003(px, py, &gamma, &theta) ;
	}
    }

    srand(0) ;
    /* save image data */
    for(y = 0 ; y < h ; y++) {
/*  	if(! (y%200)) {  */
/*  	    fprintf(stderr, "main:save row:%d\n", (int)y) ;  */
/*  	}  */
	if(TIF_CONFIG == PLANARCONFIG_CONTIG) {
	    TIFFWriteScanline(tifout,image[y],y,0) ;
	} else if(TIF_CONFIG == PLANARCONFIG_SEPARATE) {
	    uint16 s ;
	    for(s = 0 ; s < SPP ; s++)
		TIFFWriteScanline(tifout,image[y],y,s) ;
	} else {
	    fprintf(stderr, "Cannot read planar configuration of tiff image, config=%d!\n", (int)TIF_CONFIG) ;
	    exit(1) ;
	}
    }

    TIFFClose(tif);
    TIFFClose(tifout);

    fprintf(stderr, "free-ing image memory\n") ;

    /* free memory for image */

    for(y = 0 ; y < h ; y++) {
	_TIFFfree(image[y]) ;
    }

    fprintf(stderr, "free-ing image\n") ;
    free(image) ;
    free(raw_start_of_sky) ;
    free(start_of_sky) ;
    free(end_of_sky) ;
    free(final_end_of_sky) ;
    free(column_mask) ;

    // finally, copy exif tags from original to output image
    // The ICC profile seems to only be recognized in exif tags, not tiff tags, by at least
    // Geegie, and the linux image viewer "xviewer"
    if(use_exiftool) {
	char cmd[1000] ;
	if(0) {
	    //this doesn't work, maybe a bug in exiftool, it seems to create a malformed
	    //tag
	    snprintf(cmd,999,"exiftool -tagsFromFile %s %s\n", infilename, outfilename) ;
	    printf("%s", cmd) ;
	    system(cmd) ;
	} else {
	    // extract the ICC profile to external file
	    snprintf(cmd,999,"exiftool -icc_profile -b %s > skyfill_tmp.icc\n", infilename) ;
	    printf("%s", cmd) ;
	    system(cmd) ;

	    // add the ICC profile from the external file to the new tiff image
	    snprintf(cmd,999,"exiftool \"-icc_profile<=skyfill_tmp.icc\" %s\n", outfilename) ;
	    printf("%s", cmd) ;
	    system(cmd) ;

	    // system command to delete the file
	    remove("skyfill_tmp.icc") ;
	}
    }

    exit(0) ;
}
