#ifndef SAMPLE_AND_FIT_SKY_MODEL_H
#define SAMPLE_AND_FIT_SKY_MODEL_H
struct sample_point {
    float px, py ;
    float normalized_v ;
    float normalized_s ;
    float h,s,v ;
    float v_hat ;
    float abs_error[2] ; // absolute error from model of r+g+b, by location index
    float abs_h_error[2] ; // absolute error from model of hue, by location index
    float abs_s_error[2] ; // absolute error from model of sat, by location index
    float abs_v_error[2] ; // absolute error from model of val, by location index
    float cv_r, cv_g, cv_b ; // coefficiant of variation (proportion not, %) for each color component
    uint16_t n_close_neighbors ;
    uint16_t r,g,b,x,y ;
    uint8_t is_outlier ;
} ;

#define HUE_OUTLIER 0x01
#define SAT_OUTLIER 0x02
#define VAL_OUTLIER 0x04
// these have to be global
// because the nonlinear fitting function doesn't pass data, only
// parameters to the function to be minimized
extern struct sample_point *samples ;
extern int n_samples_to_optimize ;

// the value model is a single model, not 2 (or more) fit with different weights based on px
#define VAL_MODEL_IS_SINGLE

#define FULL_RAW_HSV_MODEL

// coefs to model sky HSV from input file sky area
// model is H|S|V = B0 + B1*px + B2*py + B3*px*py ;
#define MAX_COEFS 10
struct HSV_model_coefs {
    float h_coefs[MAX_COEFS];
    float s_coefs[MAX_COEFS];
    float v_coefs[MAX_COEFS];
    } ;

// two arrays will hold the underlying curve vs horizontal angle shape
#define N_BASE_KNOTS 101
extern float *val_base_knots ;
extern float *sat_base_knots ;

// for getting an interpolated point given px
#define N_BASE_KNOTSM1 100
#define KNOT_PWIDTH (1./(float)N_BASE_KNOTSM1)
#define PX2KNOT(px) ( (int)((px+0.5)/KNOT_PWIDTH) )

#define PX2_base_v(px,base_v) {\
    int k0 = PX2KNOT(px) ;\
    int k1 = k0+1 ;\
    float px_k0=(float)k0*KNOT_PWIDTH-.5;\
    float prop_k1 = (px-px_k0)/KNOT_PWIDTH ;\
    base_v = (1.-prop_k1)*val_base_knots[k0] + prop_k1*val_base_knots[k1] ;\
    }

#define PX2_base_s(px,base_s) {\
    int k0 = PX2KNOT(px) ;\
    int k1 = k0+1 ;\
    float px_k0=(float)k0*KNOT_PWIDTH-.5;\
    float prop_k1 = (px-px_k0)/KNOT_PWIDTH ;\
    base_s = (1.-prop_k1)*sat_base_knots[k0] + prop_k1*sat_base_knots[k1] ;\
    }


extern struct HSV_model_coefs raw_hsv_coefs[2] ;

// weighting for all points when fitting model for single function for entire sky
static inline float px_wgt_single(float px, float abs_error)
{
    return 1./(1.+abs_error)  ;
}

static inline float px_wgt_left(float px, float abs_error)
{
    float wgt = (.5-px) ;
    return wgt / (1.+abs_error) ;
}

static inline float px_wgt_right(float px, float abs_error)
{
    float wgt = (px+.5) ;
    return wgt / (1.+abs_error) ;
}

#ifdef SDLFKJSD
#define LOWER_PEAK

// weighting for all points when fitting model for left half of sky
static inline float px_wgt_left(float px, float abs_error)
{
#ifdef LOWER_PEAK
    px *= 1.175 ;
    px -= .25 ;
    if(px < 0.) px = 0. ;
    return (1.-3.*px*px+2.*px*px*px) / (1.+abs_error) ;
#else
    px = px*1.5 -.25;

    if(px < 0.) px = 0. ;
    if(px > 1.) px = 1. ;
    float lr_wgt = (1.-3.*px*px+2.*px*px*px) ;
    return lr_wgt*lr_wgt*lr_wgt*lr_wgt / (1.+abs_error) ;
#endif
}
    float right_wgt = 1.-left_wgt ;\

// weighting for all points when fitting model for right half of sky
static inline float px_wgt_right(float px, float abs_error)
{
#ifdef LOWER_PEAK
    px = (1.-px)*1.175 -.25;
    if(px < 0.) px = 0. ;
    return (1.-3.*px*px+2.*px*px*px) / (1.+abs_error) ;
#else
    px = (1.-px)*1.5 -.25;

    if(px < 0.) px = 0. ;
    if(px > 1.) px = 1. ;
    float lr_wgt = (1.-3.*px*px+2.*px*px*px) ;
    return lr_wgt*lr_wgt*lr_wgt*lr_wgt / (1.+abs_error) ;
#endif
}

#endif


static inline float local_h_from_pxpy(float px, float py_hc, float *coefs)
{
    float half=coefs[0] ;
    float spread=coefs[1] ;
    float speed=coefs[2] ;
    float center=coefs[3]+coefs[4]*px ;
    return spread*(py_hc-center)*speed/(1.+fabs(py_hc-center)*speed) + half ;

}

static inline float local_s_from_pxpy(float px, float py_hc, float *coefs)
{
#ifdef QUADRATIC_SAT_VERSION
    float x = py_hc-(coefs[0]+coefs[3]*px+coefs[4]*px*px) ;
    if(x < 0.0) x = 0.0 ;
    return coefs[5] + x/(coefs[1]+coefs[2]*x) ;
#else
    float x = py_hc-(coefs[0]+coefs[3]*px) ;
    if(x < 0.0) x = 0.0 ;
    return coefs[4] + x/(coefs[1]+coefs[2]*x) ;
#endif
}

static inline float modelled_local_v_from_pxpy(float px, float py_hc, float *coefs)
{
    float v_hat ;

    if(pData_fit->val_model_is_linear) {
	float mapped_x = pData_fit->FOV_horizontal*px*M_PI/180.*pData_fit->angle_factor ;
	return coefs[0] + coefs[1]*cos(mapped_x) + coefs[2]*sin(mapped_x)
		+ coefs[3]*py_hc + coefs[4]*py_hc*py_hc ;

    } else {
	float x_shifted = (pData_fit->FOV_horizontal*px+coefs[6]) ;
	float mapped_x = (x_shifted)*M_PI/180. ;
	float amplitude_by_x = (1.-coefs[4]*fabs(x_shifted)/180.) ;

	v_hat = coefs[0] + coefs[1]/(1.+exp((py_hc-coefs[2])*coefs[3])) // this is just a logistic
		+ coefs[5]*cos(mapped_x)*amplitude_by_x // this will move the logistic result up and down along the x-axis
		;
    }

    if(v_hat > 1.) v_hat=1. ;
    return v_hat ;
}

static inline float local_v_from_pxpy(float px, float py_hc, float *coefs)
{
    float v_hat ;

    float base_v ;
    PX2_base_v(px,base_v) ;

    if(0 && pData_fit->val_model_is_linear) {
/*  	float mapped_x = pData_fit->FOV_horizontal*px*M_PI/180.*pData_fit->angle_factor ;  */
/*  	return coefs[0] + coefs[1]*cos(mapped_x) + coefs[2]*sin(mapped_x)  */
/*  		+ coefs[3]*py_hc + coefs[4]*py_hc*py_hc ;  */

	return coefs[0] + coefs[1]*base_v + coefs[2]*py_hc ;

    } else {
/*  	float x_shifted = (pData_fit->FOV_horizontal*px+coefs[6]) ;  */
/*  	float mapped_x = (x_shifted)*M_PI/180. ;  */
/*  	float amplitude_by_x = (1.-coefs[4]*fabs(x_shifted)/180.) ;  */
	float vert_angle = (py_hc-pData_fit->horizon_py-coefs[3]*base_v) * pData_fit->proportion_to_radian_factor_y ;
	float cosa = cos(coefs[2]*vert_angle) ;
	v_hat = base_v*coefs[1] + coefs[0]*(exp(cosa)/exp(1.)-exp(-2.))*1.156517642 ;

#ifdef OLD
	v_hat = coefs[0] + coefs[1]/(1.+exp((py_hc-coefs[2])*coefs[3])) // this is just a logistic
		+ base_v // this will move the logistic result up and down along the x-axis
/*  		+ coefs[5]*cos(mapped_x)*amplitude_by_x // this will move the logistic result up and down along the x-axis  */
		;
#endif
    }

/*      if(v_hat > 1.) v_hat=1. ;  */
    return v_hat ;
}

#define hsv_blended(f,px,py,v) {\
    float left_wgt = px_wgt_left(px,1.) ;\
    float right_wgt = px_wgt_right(px,1.) ;\
    float lvalue = f(px,py,&raw_hsv_coefs[0]) ;\
    float rvalue = f(px,py,&raw_hsv_coefs[1]) ;\
    v = (left_wgt*lvalue + right_wgt*rvalue)/(left_wgt+right_wgt) ;\
    }

#define v_blended(f,px,py,v) {\
    float left_wgt = px_wgt_left(px,1.) ;\
    float right_wgt = px_wgt_right(px,1.) ;\
    float lvalue = f(px,py,raw_hsv_coefs[0].v_coefs) ;\
    float rvalue = f(px,py,raw_hsv_coefs[1].v_coefs) ;\
    v = (left_wgt*lvalue + right_wgt*rvalue)/(left_wgt+right_wgt) ;\
    }

#define h_blended(f,px,py,h) {\
    float left_wgt = px_wgt_left(px,1.) ;\
    float right_wgt = px_wgt_right(px,1.) ;\
    float lvalue = f(px,py,raw_hsv_coefs[0].h_coefs) ;\
    float rvalue = f(px,py,raw_hsv_coefs[1].h_coefs) ;\
    h = (left_wgt*lvalue + right_wgt*rvalue)/(left_wgt+right_wgt) ;\
    }

/*      float left_wgt = px_wgt_left(px,1.) ;\  */
/*      float right_wgt = px_wgt_right(px,1.) ;\  */

#define s_blended(f,px,py,s) {\
    float left_wgt = (.5-px) ;\
    float right_wgt = 1.-left_wgt ;\
    float lvalue = f(px,py,raw_hsv_coefs[0].s_coefs) ;\
    float rvalue = f(px,py,raw_hsv_coefs[1].s_coefs) ;\
    s = (left_wgt*lvalue + right_wgt*rvalue)/(left_wgt+right_wgt) ;\
    }

// saturation needs to pass horizon_py as an additional argument
#define hsv_blended_hpy(f,px,py,v) {\
    float left_wgt = px_wgt_left(px,1.) ;\
    float right_wgt = px_wgt_right(px,1.) ;\
    float lvalue = f(px,py,&raw_hsv_coefs[0]) ;\
    float rvalue = f(px,py,&raw_hsv_coefs[1]) ;\
    v = (left_wgt*lvalue + right_wgt*rvalue)/(left_wgt+right_wgt) ;\
    }


static inline float h_from_pxpy(float px, float py_hc, int full_hsv_model_level)
{
    if(full_hsv_model_level==1)
	return local_h_from_pxpy(px,py_hc,raw_hsv_coefs[0].h_coefs) ;

    float h ;
    h_blended(local_h_from_pxpy, px, py_hc, h) ;
    return h ;
}

static inline float s_from_pxpy(float px, float py_hc, int full_hsv_model_level, float horizon_py)
{
    if(full_hsv_model_level==1)
	return local_s_from_pxpy(px,py_hc,raw_hsv_coefs[0].s_coefs) ;

    float s ;
    s_blended(local_s_from_pxpy, px, py_hc, s) ;
    return s ;
}

static inline float modelled_v_from_pxpy(float px, float py_hc, int full_hsv_model_level)
{
    if(full_hsv_model_level==1)
	return modelled_local_v_from_pxpy(px,py_hc,raw_hsv_coefs[0].v_coefs) ;

    float v ;
    v_blended(modelled_local_v_from_pxpy, px, py_hc, v) ;
    return v ;
}

static inline float v_from_pxpy(float px, float py_hc, int full_hsv_model_level)
{
#ifdef VAL_MODEL_IS_SINGLE
	return local_v_from_pxpy(px,py_hc,raw_hsv_coefs[0].v_coefs) ;
#else
    if(full_hsv_model_level==1)
	return local_v_from_pxpy(px,py_hc,raw_hsv_coefs[0].v_coefs) ;

    float v ;
    v_blended(local_v_from_pxpy, px, py_hc, v) ;
    return v ;
#endif
}

float fit_full_HSV_model(int n_samples, int location_index, int calc_error, SKYFILL_DATA_t *pData) ;

void reset_sample_errors(int n_samples) ;

int sample_sky_points(int n_per_column, int n_columns,tdata_t *image,int16_t *start_of_sky,int16_t *end_of_sky, SKYFILL_DATA_t *pData, int need_model_fit) ;
int read_samples_and_fit_model(char *samplefile, SKYFILL_DATA_t *pData) ;
void output_sample_data(int n_samples, char filename[]) ;

#endif
