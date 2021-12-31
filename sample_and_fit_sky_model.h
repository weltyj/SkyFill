#ifndef SAMPLE_AND_FIT_SKY_MODEL_H
#define SAMPLE_AND_FIT_SKY_MODEL_H
struct sample_point {
    float px, py ;
    float h,s,v ;
    float v_hat ;
    float abs_error[2] ; // absolute error from model of r+g+b, by location index
    float abs_h_error[2] ; // absolute error from model of hue, by location index
    float abs_s_error[2] ; // absolute error from model of sat, by location index
    float abs_v_error[2] ; // absolute error from model of val, by location index
    uint16_t r,g,b,x,y ;
} ;

// these have to be global
// because the nonlinear fitting function doesn't pass data, only
// parameters to the function to be minimized
extern struct sample_point *samples ;
extern int n_samples_to_optimize ;

#define FULL_RAW_HSV_MODEL

// coefs to model sky HSV from input file sky area
// model is H|S|V = B0 + B1*px + B2*py + B3*px*py ;
struct HSV_model_coefs {
    float raw_HSV_coefs[3][5];
    float S_inv_coef ; // coef for saturation on the inverse function
    } ;

extern struct HSV_model_coefs raw_hsv_coefs[2] ;

// weighting for all points when fitting model for single function for entire sky
static inline float px_wgt_single(float px, float abs_error)
{
    return 1./(1.+abs_error)  ;
}

// weighting for all points when fitting model for left half of sky
static inline float px_wgt_left(float px, float abs_error)
{
    px *= 1.175 ;
    px -= .25 ;
    if(px < 0.) px = 0. ;
    return (1.-3.*px*px+2.*px*px*px) / (1.+abs_error) ;
}

// weighting for all points when fitting model for right half of sky
static inline float px_wgt_right(float px, float abs_error)
{
    px = (1.-px)*1.175 -.25;
    if(px < 0.) px = 0. ;
    return (1.-3.*px*px+2.*px*px*px) / (1.+abs_error) ;
}


static inline float local_h_from_pxpy(float px, float py_hc, struct HSV_model_coefs *p)
{
    return p->raw_HSV_coefs[0][0] + p->raw_HSV_coefs[0][1]*px + p->raw_HSV_coefs[0][2]*log(1.+py_hc) ;
}

static inline float local_ty_from_pxpy(float px, float py_hc, struct HSV_model_coefs *p, float horizon_py)
{
    py_hc -= horizon_py ;

    return 1./(1.+exp(-py_hc*p->S_inv_coef)) ; // really a simple logistic centered on horizon_py
}

static inline float local_s_from_pxpy(float px, float py_hc, struct HSV_model_coefs *p)
{
    // nonlinear model now, fit with amoeba (Nelder Meade)
    return p->raw_HSV_coefs[1][0]/(1.+exp(-(py_hc-p->raw_HSV_coefs[1][1])*p->raw_HSV_coefs[1][2]))
	    +p->raw_HSV_coefs[1][3] +p->raw_HSV_coefs[1][4]*py_hc ;
}

static inline float local_v_from_pxpy(float px, float py_hc, struct HSV_model_coefs *p)
{
    float horz_angle = pData_fit->FOV_horizontal*px*M_PI/180. ;
/*      fprintf(stderr, "LVFP: HA:%f\n", horz_angle) ;  */
    return p->raw_HSV_coefs[2][0] + p->raw_HSV_coefs[2][1]*cos(2.*horz_angle) + p->raw_HSV_coefs[2][2]*sin(2.*horz_angle) + p->raw_HSV_coefs[2][3]*py_hc ;
}

#define hsv_blended(f,px,py,v) {\
    float left_wgt = px_wgt_left(px,1.) ;\
    float right_wgt = px_wgt_right(px,1.) ;\
    float lvalue = f(px,py,&raw_hsv_coefs[0]) ;\
    float rvalue = f(px,py,&raw_hsv_coefs[1]) ;\
    v = (left_wgt*lvalue + right_wgt*rvalue)/(left_wgt+right_wgt) ;\
    }

// saturation needs to pass horizon_py as an additional argument
#define hsv_blended_hpy(f,px,py,v) {\
    float left_wgt = px_wgt_left(px,1.) ;\
    float right_wgt = px_wgt_right(px,1.) ;\
    float lvalue = f(px,py,&raw_hsv_coefs[0]) ;\
    float rvalue = f(px,py,&raw_hsv_coefs[1]) ;\
    v = (left_wgt*lvalue + right_wgt*rvalue)/(left_wgt+right_wgt) ;\
    }


static inline float h_from_pxpy(float px, float py, int full_hsv_model_level)
{
    if(full_hsv_model_level==1)
	return local_h_from_pxpy(px,py,&raw_hsv_coefs[0]) ;

    float h ;
    hsv_blended(local_h_from_pxpy, px, py, h) ;
    return h ;
}

static inline float s_from_pxpy(float px, float py, int full_hsv_model_level, float horizon_py)
{
    if(full_hsv_model_level==1)
	return local_s_from_pxpy(px,py,&raw_hsv_coefs[0]) ;

    float s ;
    hsv_blended_hpy(local_s_from_pxpy, px, py, s) ;
    return s ;
}

static inline float v_from_pxpy(float px, float py, int full_hsv_model_level)
{
    if(full_hsv_model_level==1)
	return local_v_from_pxpy(px,py,&raw_hsv_coefs[0]) ;

    float v ;
    hsv_blended(local_v_from_pxpy, px, py, v) ;
    return v ;
}

float fit_full_HSV_model(int n_samples, int location_index, int calc_error, SKYFILL_DATA_t *pData) ;

void reset_sample_errors(int n_samples) ;

int sample_sky_points(int n_per_column, int n_columns,tdata_t *image,int16_t *start_of_sky,int16_t *end_of_sky, SKYFILL_DATA_t *pData) ;

#endif
