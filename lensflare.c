#ifndef LENSFLARE_MAIN
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <malloc.h>
#include <strings.h>
#include <string.h>
#include <ctype.h>

#include "skyfill_tif.h"
#include "colorspace_conversions.h"
#else

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <malloc.h>
#include <strings.h>
#include <string.h>
#include <ctype.h>
#include "tiffio.h"
#include "colorspace_conversions.h"

// half space coefficients for 2D
struct HS {
    float A, B, C ;
} ;

// in case Windows MSVC, or 
#define _USE_MATH_DEFINES
#include <math.h>

// A few truly global variables
int32_t IMAGE_HEIGHT ;
int32_t IMAGE_WIDTH ;
uint16_t IMAGE_NSAMPLES ; // must be 16 bit, needed by call to tif library
int IMAGE_HAS_ALPHA ;
float p_half_image_width ; // scaled value of half the image width (use as X origin) ;

/* coordinate system:
	TIF files are read.  data is stored in image[] array.  Top left is (0,0), bottom right is (width-1,height-1) ;
*/

#define tif_b(image,x,y) (((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+2])

// macros to get or set rgb and alpha 
#define tif_get3c(image,x,y,r,g,b) {\
		r = ((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+0] ; \
		g = ((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+1] ; \
		b = ((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+2] ; \
		}
#define tif_get3cv(image,x,y,array) {\
		array[0] = ((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+0] ; \
		array[1] = ((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+1] ; \
		array[2] = ((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+2] ; \
		}
#define tif_get4cv(image,x,y,array) {\
		array[0] = ((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+0] ; \
		array[1] = ((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+1] ; \
		array[2] = ((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+2] ; \
		array[3] = ((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+3] ; \
		}
#define dtif_get3c(image,x,y,r,g,b) {\
		fprintf(stderr, "get image data at %d, %d\n", x, y) ; \
		r = ((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+0] ; \
		g = ((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+1] ; \
		b = ((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+2] ; \
		}
#define dtif_get4c(image,x,y,r,g,b,a) {\
		fprintf(stderr, "get image data at %d, %d\n", x, y) ; \
		r = ((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+0] ; \
		g = ((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+1] ; \
		b = ((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+2] ; \
		a = ((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+3] ; \
		}

#define tif_get4c(image,x,y,r,g,b,a) {\
		r = ((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+0] ; \
		g = ((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+1] ; \
		b = ((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+2] ; \
		a = ((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+3] ; \
		}

#define tif_set3c(image,x,y,r,g,b) {\
		((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+0] = (uint16_t)r; \
		((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+1] = (uint16_t)g; \
		((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+2] = (uint16_t)b; \
		}

#define tif_set4c(image,x,y,r,g,b,a) {\
		((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+0] = (uint16_t)r; \
		((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+1] = (uint16_t)g; \
		((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+2] = (uint16_t)b; \
		((uint16_t *)(image[(y)]))[IMAGE_NSAMPLES*(x)+3] = (uint16_t)a; \
		}

#define MAX16 65535
#define HALF16 32767
#define MAX16f 65535.f
#define FMAX 1.e30f
// value which is black
#define BLK16 8621

#define MIN(a,b) ( (a) < (b) ? (a) : (b) )
#define MAX(a,b) ( (a) > (b) ? (a) : (b) )

#define IMAGE_PIXEL_X_TO_RELATIVE(x) ((float)( ((float)(x)/(float)(IMAGE_WIDTH-1)) -p_half_image_width))
#define IMAGE_PIXEL_Y_TO_RELATIVE(y) ((float)(1.f-(float)(y)/(float)(IMAGE_HEIGHT-1)))

#define IMAGE_RELATIVE_TO_PIXEL_X(px) (int)(((px)+p_half_image_width)*(float)(IMAGE_WIDTH-1)+0.5)
#define IMAGE_RELATIVE_TO_PIXEL_Y(py) (int)((1.-(py))*(float)(IMAGE_HEIGHT-1)+0.5)

#endif

// construct half space equation orthoganal to two points
// x2,y2 is on the line, x1,y1 is a postive distance from the line
struct HS half_space_from_pts(float x1, float y1, float x2, float y2)
{
    struct HS s ;
    float dx = x2-x1 ;
    float dy = y2-y1 ;
    s.A = -dx ;
    s.B = -dy ;
    s.C = -(s.A*x2 + s.B*y2) ;

    float denom = sqrt(dx*dx + dy*dy) ;
    s.A /= denom ;
    s.B /= denom ;
    s.C /= denom ;

    return s ;
}

#define fdither(f16) {\
    float dither = (float)(((rand()%128) - 64)*8) ; \
    f16 += dither ; \
    if(f16 < 0.0) f16 = 0.0 ; \
    if(f16 > MAX16f-1.) f16 = MAX16f-1. ; \
}

void new_map_rfunc(float r_func, float *e, float *b, float *c)
{

    //if(r_func < 0.5) r_func = 0.5 - r_func ;
    // b is radius cutoff thresh
    // c is base
    float ep5 = exp((1.-r_func)*2.5)-1. ;

    // b is thresh
    float bp5 = -0.86*r_func*r_func-0.11*r_func+0.97 ;

    float cp5 = r_func ;

    if(r_func < .5) {
	*e = exp((1.-r_func)*2.5)-1. ;

	// b is thresh
	*b = -0.86*r_func*r_func-0.11*r_func+0.97 ;

	// c is base
	*c = r_func ;
    } else if(r_func < 1.) {
	float p1 = (r_func-.5)/.5 ;
	float pp5 = 1.-p1 ;
	*e = ep5*pp5 +1.00*p1 ; // e => ep5 to 1
	*b = bp5*pp5 +0.92*p1 ;
	*c = cp5*pp5 +0.90*p1 ;
    } else if(r_func < 2.) {
	*e = 1. ;
	*b = (2.-r_func)*.92 ;
	*c = .9 ;
    } else {
	*b=0 ;
	*e=r_func-1. ;
    }

/*      fprintf(stderr, "Rfunc:%f e:%f b:%f c:%f\n", r_func, *e, *b, *c) ;  */
}

void best_map_rfunc(float r_func, float *e, float *b, float *c)
{
    if(r_func < 1.) {
	*e = exp((1.-r_func)*2.5)-1. ;

	// b is thresh
	*b = -0.86*r_func*r_func-0.11*r_func+0.97 ;

	// c is base
	*c = r_func ;
    } else {
	*b=0 ;
	*e=r_func ;
    }

/*      fprintf(stderr, "Rfunc:%f e:%f b:%f c:%f\n", r_func, *e, *b, *c) ;  */
}

void map_rfunc(float r_func, float *e, float *b, float *c)
{
    new_map_rfunc(r_func, e, b, c) ;
/*      fprintf(stderr, "Rfunc:%f e:%f b:%f c:%f\n", r_func, *e, *b, *c) ;  */
    return ;
    if(r_func < 1.) {
	best_map_rfunc(1.-r_func, e, b, c) ;
    } else if(r_func < 2.) {
	*e = exp(2.5)-1. ;
	*b = 1. ;
	*c = r_func-1. ;
    } else {
	*b=0 ;
	*e=r_func-2. ;
    }

    fprintf(stderr, "Rfunc:%f e:%f b:%f c:%f\n", r_func, *e, *b, *c) ;
}
// special function to map hsv to rgb, inplace, output is scaled 0 to 1

void hsv2rgb_1(float *in)
{
    uint16_t rgb[3] ;
    hsv2rgb16(in,rgb) ;
    in[0] = (float)rgb[0] / MAX16f ;
    in[1] = (float)rgb[1] / MAX16f ;
    in[2] = (float)rgb[2] / MAX16f ;
}

#define MAX_GHOSTS_IN_LINE 25

struct ghost_ring_struct {
    int n_ghosts ;
    float dist ;
    float radius ;
    float MAX_CA ;
    float angle[MAX_GHOSTS_IN_LINE] ;
    float I[MAX_GHOSTS_IN_LINE][3] ;
    int n_leaves ;
    float pcircle[MAX_GHOSTS_IN_LINE] ;
    float r_func[MAX_GHOSTS_IN_LINE] ;
    } ;

struct ghost_line_struct {
    int n_ghosts ;
    float rotate_angle ;
    float MAX_CA ;
    float pcircle[MAX_GHOSTS_IN_LINE] ;
    float r_func[MAX_GHOSTS_IN_LINE] ;
    float dist[MAX_GHOSTS_IN_LINE] ;
    float radius[MAX_GHOSTS_IN_LINE] ;
    float I[MAX_GHOSTS_IN_LINE][3] ;
    int n_leaves ;
    } ;

struct ghost_image_struct {
    float x,y ;
    float radius ;
    float I[3] ;
    float e ;
    float b ;
    float c ;
    float MAX_CA ;
    int n_leaves ;
    float pcircle ;
    float r_func ;
    } ;


#define MAX_GHOST_LINES 10
#define MAX_GHOST_RINGS 10
struct lens_flare_metadata
{

    int n_ghost_lines ;
    int n_ghost_rings ;
    int n_leaves ;

    int flare_type ;
    int flare_radius ;
    float flare_reverse_radius ;
    float sun_x ;
    float sun_y ;

    float TRUE_CENTER_X ;
    float TRUE_CENTER_Y ;

    struct ghost_ring_struct ghost_rings[MAX_GHOST_LINES] ;

    struct ghost_line_struct ghost_lines[MAX_GHOST_RINGS] ;

    int n_ghosts ;

    float scale_hsv[3] ;
    float scale_radius ;
    float scale_dist ;
    float rotate ; // about sun center
    float D_per_angle ;
    float E_per_angle ;
} ;

void render_simple_lens_flare(char *filename, struct lens_flare_metadata *md, tdata_t *image) ;
void render_lens_flare(char *filename, struct lens_flare_metadata *md, tdata_t *image) ;

void render_ghosts(struct ghost_image_struct *gi, int n_ghosts, tdata_t *image, struct lens_flare_metadata *md)
{
    for(int i=0 ; i < n_ghosts ; i++) {
	// render_ghost_flare(&gi[i]) ;
	// if n_leaves was 4, would have to increase radius by sqrt(2) to make sure
	// the entire area was covered
	int x0=gi[i].x-gi[i].radius*sqrt(2.)-1-gi[i].MAX_CA ;
	int x1=gi[i].x+gi[i].radius*sqrt(2.)+1+gi[i].MAX_CA ;
	int y0=gi[i].y-gi[i].radius*sqrt(2.)-1-gi[i].MAX_CA ;
	int y1=gi[i].y+gi[i].radius*sqrt(2.)+1+gi[i].MAX_CA ;

	fprintf(stderr, "Render flare ghost %d of %d at x,y %f,%f\n", i+1, n_ghosts, gi[i].x, gi[i].y) ;

	if(x0 < 0) x0=0 ;
	if(x0 > IMAGE_WIDTH-1) x0=IMAGE_WIDTH-1 ;
	if(x1 < 0) x1=0 ;
	if(x1 > IMAGE_WIDTH-1) x1=IMAGE_WIDTH-1 ;

	if(y0 < 0) y0=0 ;
	if(y0 > IMAGE_HEIGHT-1) y0=IMAGE_HEIGHT-1 ;
	if(y1 < 0) y1=0 ;
	if(y1 > IMAGE_HEIGHT-1) y1=IMAGE_HEIGHT-1 ;

	// create the normal equations that bound the space equivalent
	// to the leaves for the apeture.

	struct HS ne[12][3] ;

	float ca_offset[3] = {gi[i].MAX_CA, 0., -gi[i].MAX_CA} ;
	float ca_center_x[3], ca_center_y[3] ;

	for(int c=0 ; c < 3 ; c++) {
	    float dx_sun=(gi[i].x-md->sun_x) ;
	    float dy_sun=(gi[i].y-md->sun_y) ;
	    float r=sqrt(dx_sun*dx_sun+dy_sun*dy_sun) ;
	    dx_sun /= 1000. ;
	    dy_sun /= 1000. ;
	    ca_center_x[c] = ca_offset[c]*dx_sun ;
	    ca_center_y[c] = ca_offset[c]*dy_sun ;
	}

	for(int l=0 ; l < gi[i].n_leaves ; l++) {
	    for(int c=0 ; c < 3 ; c++) {
		float angle=0+360./gi[i].n_leaves*(float)l ;
		float radians = angle*M_PI/180. ;
		float dx = cosf(radians)*gi[i].radius ;
		float dy = sinf(radians)*gi[i].radius ;
		float ca_x = ca_center_x[c] ;
		float ca_y = ca_center_y[c] ;
		ne[l][c] = half_space_from_pts(gi[i].x+ca_x, gi[i].y+ca_y, gi[i].x+ca_x+dx, gi[i].y+ca_y+dy) ;
	    }
	}

	if(gi[i].n_leaves > 2) {
	    int midx = (x0+x1)/2 ;
	    int midy = (y0+y1)/2 ;

	    for(int x=x0 ; x <= x1 ; x++) {
		for(int y=y0 ; y <= y1 ; y++) {
		    float flare_rgb[3]={0.,0.,0.} ;
		    int found=0 ;

		    for(int c=0 ; c < 3 ; c++) {
			// find the minimum d, where d > 0.
			float d=1.e30 ;
			int inside=1 ;

			for(int l=0 ; l < gi[i].n_leaves ; l++) {
			    float this_d = (ne[l][c].A*x+ne[l][c].B*y+ne[l][c].C) ;
			    if(this_d <= 0.) {
				inside=0 ;
				break ;
			    }
			    if(this_d > 0. && this_d < d) {
				d = this_d ;
			    }
			}

			if(inside && d < gi[i].radius*2) {
			    // what is the distance from the edge of an equivalent circle?
			    float dx=x-(ca_center_x[c]+gi[i].x) ;
			    float dy=y-(ca_center_y[c]+gi[i].y) ;
			    float r_center = sqrt(dx*dx+dy*dy) ;
			    float d_circle = gi[i].radius-r_center ;
			    d = gi[i].pcircle*d_circle + (1.-gi[i].pcircle)*d ;

			    if(d <= 0.) inside=0 ;
			}

			if(inside && d < gi[i].radius*2) {

			    float f = d ;

			    if(f < 0.) continue ;
			    if(f > 1.) f=1. ;

			    // not technically correct except
			    // where a leaf plane is tangent to the
			    // enclosed circle, but let's try it
			    float pr= 1.-d/gi[i].radius ;

			    if(pr < 0.) continue ;
			    if(pr > 1.) continue ;

			    // intensity of as function of radius
			    float I ;
			    if(gi[i].b > 0.) {
				if(pr < gi[i].b) {
				    I=gi[i].c+exp(-fabsf(gi[i].b-pr)*gi[i].e)*(1.-gi[i].c) ;
				} else {
				    I=1. - (pr-gi[i].b)/(1.-gi[i].b) ;
				}
			    } else {
				I = powf(1.-pr, gi[i].e) ;
			    }


#define gamma(x) powf((x),3.0)
#define ungamma(x) powf((x),1./3.0)

			    I *= gi[i].I[c]*f ;
			    flare_rgb[c] += gamma(I) ;
			    found=1 ;
			}
		    }

		    if(found) {

			float r,g,b ;
			tif_get3c(image,x,y,r,g,b) ;

			r += ungamma(flare_rgb[0])*MAX16 ;
			g += ungamma(flare_rgb[1])*MAX16 ;
			b += ungamma(flare_rgb[2])*MAX16 ;

			fdither(r) ;
			fdither(g) ;
			fdither(b) ;

			tif_set3c(image,x,y,r,g,b) ;
		    }
		}
	    }

	} else {
#ifdef CIRCLES
	    for(int x=x0 ; x <= x1 ; x++) {
		float dx = x-ai[i].x ;
		for(int y=y0 ; y <= y1 ; y++) {
		    float dy = y-ai[i].y ;
		    float r_xy = sqrt(dx*dx+dy*dy) ;
		    float f = ai[i].radius-r_xy ;

		    if(f < 0.) continue ;
		    if(f > 1.) f=1. ;

		    float pr= r_xy/ai[i].radius ;
		    float I ;
		    if(pr < .8) {
			I=ai[i].b+exp(-fabsf(0.8-pr)*ai[i].e)*(1.-ai[i].b) ;
		    } else {
			I=1. - (pr-.8)/.2 ;
		    }

		    I *= ai[i].I*f ;

		    float r,g,b ;
		    tif_get3c(image,x,y,r,g,b) ;

		    r += I*MAX16 ;
		    g += I*MAX16 ;
		    b += I*MAX16 ;

		    fdither(r) ;
		    fdither(g) ;
		    fdither(b) ;

		    tif_set3c(image,x,y,r,g,b) ;
		}
	    }
#endif
	}
    } // end of for() render ghosts loop
}

int string_is_numeric(char *buf)
{
    int i=0 ;
    int found_digit=0 ;
    while(buf[i] != '\0') {

	if(isspace(buf[i])) {
	    i++ ;
	} else if(isdigit(buf[i])) {
	    found_digit=1 ;
	    i++ ;
	} else if(buf[i] == '+') {
	    i++ ;
	} else if(buf[i] == '-') {
	    i++ ;
	} else if(buf[i] == '.') {
	    i++ ;
	} else {
	    return 0 ;
	}
    }

    if(found_digit)
	return 1 ;

    return 0 ;
}

#define BUFSIZE 512
char * get_next_line(FILE *fp, char *buf, char *filename)
{
    while(1) {
	if(fgets(buf, BUFSIZE-1, fp) == NULL) {
	    fprintf(stderr, "EOF while reading %s for data\n", filename) ;
	    return NULL ;
	}

	int i=0 ;
	while(isspace(buf[i]) && i < BUFSIZE-1) {
	    i++ ;
	}

	if(i >= BUFSIZE-1 || buf[i] == '#' || buf[i] == '\0')
	    continue ;
	break ;
    }

    return buf ;
}


int render_ghost_lines_and_rings(tdata_t *image, struct lens_flare_metadata *md)
{

	// convert data read into individual ghosts
	int n_ghosts=0 ;

	float sun_dx = (md->sun_x-md->TRUE_CENTER_X) ;
	float sun_dy = (md->sun_y-md->TRUE_CENTER_Y) ;
	float sun_dist = sqrt(sun_dx*sun_dx+sun_dy*sun_dy) ;
	float sun_pdx = sun_dx/sun_dist ;
	float sun_pdy = sun_dy/sun_dist ;
	float line_length = sun_dist*2. ;

#define MAX_GHOSTS 200
	struct ghost_image_struct *gi = (struct ghost_image_struct *)calloc(MAX_GHOSTS, sizeof(struct ghost_image_struct)) ;

	for(int gl=0 ; gl < md->n_ghost_lines ; gl++) {
	    struct ghost_line_struct *pGL = &md->ghost_lines[gl] ;
	    float c=cosf((pGL->rotate_angle)*M_PI/180.) ;
	    float s=sinf((pGL->rotate_angle)*M_PI/180.) ;
	    float px = -sun_pdx*c + sun_pdy*s ;
	    float py = -sun_pdx*s - sun_pdy*c ;

	    fprintf(stderr, "Convert ghost line %d with %d ghosts\n", gl, pGL->n_ghosts) ;


	    for(int i=0 ; i < pGL->n_ghosts ; i++) {
		if(n_ghosts == MAX_GHOSTS) {
		    fprintf(stderr, "LENSFLARE: exceeded maximum allowed total number of ghosts (%d)\n", MAX_GHOSTS) ;
		    return 1 ;
		}
		gi[n_ghosts].x = md->sun_x + pGL->dist[i]*px*line_length ;
		gi[n_ghosts].y = md->sun_y + pGL->dist[i]*py*line_length ;
		gi[n_ghosts].pcircle = pGL->pcircle[i] ;
		gi[n_ghosts].r_func = pGL->r_func[i] ;
		gi[n_ghosts].radius = pGL->radius[i] ;
		gi[n_ghosts].I[0] = pGL->I[i][0] ;
		gi[n_ghosts].I[1] = pGL->I[i][1] ;
		gi[n_ghosts].I[2] = pGL->I[i][2] ;
		gi[n_ghosts].MAX_CA = pGL->MAX_CA ;
		map_rfunc(gi[n_ghosts].r_func, &gi[n_ghosts].e, &gi[n_ghosts].b, &gi[n_ghosts].c) ;
		//if(gi[n_ghosts].pcircle > 1.) gi[n_ghosts].pcircle = 1.;
		//if(gi[n_ghosts].pcircle < 0.) gi[n_ghosts].pcircle = 0.;
		gi[n_ghosts].n_leaves = md->n_leaves ;


/*  		fprintf(stderr, "ghosts %d, pcircle:%f, radius:%f, rfunc:%f ebc:%f,%f,%f\n", n_ghosts,  */
/*  				gi[n_ghosts].pcircle, gi[n_ghosts].radius, gi[n_ghosts].r_func, gi[n_ghosts].e, gi[n_ghosts].b, gi[n_ghosts].c) ;  */

		n_ghosts++ ;
	    }
	    fprintf(stderr, "final n_ghosts %d\n", n_ghosts) ;
	}

	for(int gr=0 ; gr < md->n_ghost_rings ; gr++) {
	    struct ghost_ring_struct *pGR = &md->ghost_rings[gr] ;
	    for(int i=0 ; i < pGR->n_ghosts ; i++) {
		if(n_ghosts == MAX_GHOSTS) {
		    fprintf(stderr, "LENSFLARE: exceeded maximum allowed total number of ghosts (%d)\n", MAX_GHOSTS) ;
		    return 1 ;
		}
		float c=cosf((pGR->angle[i])*M_PI/180.) ;
		float s=sinf((pGR->angle[i])*M_PI/180.) ;
		float px = -sun_pdx*c + sun_pdy*s ;
		float py = -sun_pdx*s - sun_pdy*c ;
		gi[n_ghosts].x = md->sun_x + pGR->dist*px*line_length ;
		gi[n_ghosts].y = md->sun_y + pGR->dist*py*line_length ;
		gi[n_ghosts].radius = pGR->radius ;
		gi[n_ghosts].I[0] = pGR->I[i][0] ;
		gi[n_ghosts].I[1] = pGR->I[i][1] ;
		gi[n_ghosts].I[2] = pGR->I[i][2] ;
		gi[n_ghosts].MAX_CA = pGR->MAX_CA ;
		gi[n_ghosts].r_func = pGR->r_func[i] ;
		gi[n_ghosts].pcircle = pGR->pcircle[i] ;
		map_rfunc(gi[n_ghosts].r_func, &gi[n_ghosts].e, &gi[n_ghosts].b, &gi[n_ghosts].c) ;
		gi[n_ghosts].n_leaves = md->n_leaves ;
		//if(gi[n_ghosts].pcircle > 1.) gi[n_ghosts].pcircle = 1.;
		//if(gi[n_ghosts].pcircle < 0.) gi[n_ghosts].pcircle = 0.;

		n_ghosts++ ;
	    }
	}

	render_ghosts(gi, n_ghosts, image, md) ;
	free(gi) ;

	return 0 ;
}

int parse_lens_flare_file_and_render_ghosts(char *filename, struct lens_flare_metadata *md, tdata_t *image)
{
    FILE *fp = fopen(filename, "r") ;

    if(fp != NULL) {
	char cmd[100] ;
	int rval ;
	int str_is_filled=0 ;
	while(1) {
	    char str[BUFSIZE] ;
	    if(str_is_filled == 0) {
		if(get_next_line(fp, str, filename) == NULL) {
		    break ;
		}
	    }

	    str_is_filled=0 ;

	    rval = sscanf(str, "%s", cmd) ;

	    fprintf(stderr, "lensflare rval %d, string %s\n", rval, cmd) ;
	    if(rval < 1) break ;

	    if(cmd[0] == '#') {
		// get rest of line and continue
		if(fgets(str, BUFSIZE-1, fp) == NULL) {
		    fprintf(stderr, "FATAL: died while reading %s for data\n", filename) ;
		}
		continue ;
	    }

	    if(!strcasecmp("BLACKSKY", cmd)) {
		float r,g,b ;
		//int rval = sscanf(str, "%s%f%f%f", cmd, &r, &g, &b) ;
		//printf("BLACKSKY rval=%d\n", rval) ;
		if(sscanf(str, "%s%f%f%f", cmd, &r, &g, &b) != 4) {
		    r=g=b=0. ;
		}
		r *= MAX16f ;
		g *= MAX16f ;
		b *= MAX16f ;
		for(int x=0 ; x < IMAGE_WIDTH ; x++) {
		    for(int y=0 ; y < IMAGE_HEIGHT ; y++) {
			tif_set3c(image,x,y,r,g,b) ;
		    }
		}
		continue ;
	    }

	    if(!strcasecmp("FILE", cmd)) {
		//get_next_line(fp, str, filename) ;
		char newfilename[512] ;
		sscanf(str, "%*s%s", newfilename) ;
		struct lens_flare_metadata new_md ;
		new_md = *md ;
		parse_lens_flare_file_and_render_ghosts(newfilename, &new_md, image) ;
		continue ;
	    }

	    if(!strcasecmp("LEAVES", cmd)) {
		//get_next_line(fp, str, filename) ;
		sscanf(str, "%*s%d", &md->n_leaves) ;
		continue ;
	    }

	    if(!strcasecmp("FLARE", cmd)) {
		//get_next_line(fp, str, filename) ;
		sscanf(str, "%*s%d%d%f%f%f", &md->flare_type, &md->flare_radius,
					&md->flare_reverse_radius, &md->D_per_angle, &md->E_per_angle) ;

		render_lens_flare(filename, md, image) ;

		continue ;
	    }

	    if(!strcasecmp("SCALE", cmd)) {
		//get_next_line(fp, str, filename) ;
		sscanf(str, "%*s%f%f%f%f%f",
				    &md->scale_hsv[0], &md->scale_hsv[1], &md->scale_hsv[2],
				    &md->scale_radius, &md->scale_dist) ;
		continue ;
	    }

	    if(!strcasecmp("ROTATE", cmd)) {
		//get_next_line(fp, str, filename) ;
		sscanf(str, "%*s%f", &md->rotate) ;
		continue ;
	    }

	    if(!strcasecmp("IMAGE_CENTER", cmd)) {
		//get_next_line(fp, str, filename) ;
		sscanf(str, "%*s%f%f", &md->TRUE_CENTER_X, &md->TRUE_CENTER_Y) ;
		md->TRUE_CENTER_X *= IMAGE_WIDTH ;
		md->TRUE_CENTER_Y *= IMAGE_HEIGHT ;
		continue ;
	    }

	    if(!strcasecmp("SUN_CENTER", cmd)) {
		float px,py ;
		sscanf(str, "%*s%f%f", &px, &py) ;
		fprintf(stderr, "image center is at %0.0lf,%0.0lf\n", md->TRUE_CENTER_X, md->TRUE_CENTER_Y) ;
		fprintf(stderr, "Sun center was at %0.0lf,%0.0lf", md->sun_x, md->sun_y) ;
		md->sun_x = (px*(float)IMAGE_WIDTH+0.5) ;
		md->sun_y = (py*(float)IMAGE_HEIGHT+0.5) ;
		fprintf(stderr, " now at %0.0lf,%0.0lf\n", md->sun_x, md->sun_y) ;
		continue ;
	    }

	    if(!strcasecmp("MOVE_CENTER", cmd)) {
		float delta_px,delta_py ;
		sscanf(str, "%*s%f%f", &delta_px, &delta_py) ;
		fprintf(stderr, "image center was at %0.0lf,%0.0lf\n", md->TRUE_CENTER_X, md->TRUE_CENTER_Y) ;
		fprintf(stderr, "Sun center was at %0.0lf,%0.0lf\n", md->sun_x, md->sun_y) ;
		md->sun_x += (delta_px*(float)IMAGE_WIDTH+0.5) ;
		md->sun_y += (delta_py*(float)IMAGE_HEIGHT+0.5) ;
		md->TRUE_CENTER_X += (delta_px*(float)IMAGE_WIDTH+0.5) ;
		md->TRUE_CENTER_Y += (delta_py*(float)IMAGE_HEIGHT+0.5) ;
		fprintf(stderr, "image center now at %0.0lf,%0.0lf\n", md->TRUE_CENTER_X, md->TRUE_CENTER_Y) ;
		fprintf(stderr, "Sun center now at %0.0lf,%0.0lf\n", md->sun_x, md->sun_y) ;
		continue ;
	    }

	    if(!strcasecmp("GHOST_LINE", cmd)) {
		struct ghost_line_struct *pGL = &md->ghost_lines[md->n_ghost_lines] ;
		float p_inv_circle ;
		if(md->n_ghost_lines == MAX_GHOST_LINES) {
		    fprintf(stderr, "LENSFLARE: exceeded maximum allowed number of ghost lines (%d)\n", MAX_GHOST_LINES) ;
		    return 1 ;
		}
		md->n_ghost_lines++ ;
		//get_next_line(fp, str, filename) ;
		sscanf(str, "%*s%f%f", &pGL->rotate_angle, &pGL->MAX_CA) ;
		pGL->rotate_angle += md->rotate ;
		pGL->n_ghosts=0 ;
		int i=0 ;

		while(1) {
		    if(get_next_line(fp, str, filename) == NULL)
			break ;
		    if(!string_is_numeric(str)) {
			// set flag that the next line is already loaded into str
			str_is_filled=1 ;
			break ;
		    }

		    if(i == MAX_GHOSTS_IN_LINE) {
			fprintf(stderr, "LENSFLARE: exceeded maximum allowed number of ghosts in a line (%d)\n", MAX_GHOSTS_IN_LINE) ;
			return 1 ;
		    }

		    sscanf(str, "%f%f%f%f%f%f%f", &pGL->dist[i], &pGL->radius[i], &pGL->I[i][0], &pGL->I[i][1], &pGL->I[i][2],
			&p_inv_circle, &pGL->r_func[i]) ;
		    pGL->pcircle[i] = 1.-p_inv_circle ;
		    pGL->dist[i] *= md->scale_dist ;
		    pGL->radius[i] *= md->scale_radius ;
		    pGL->I[i][0] += md->scale_hsv[0] ;
		    pGL->I[i][1] *= md->scale_hsv[1] ;
		    pGL->I[i][2] *= md->scale_hsv[2] ;
		    hsv2rgb_1(pGL->I[i]) ;
		    pGL->n_ghosts++ ;
		    i++ ;

		    //fprintf(stderr, "GL: %f %f %f %f %f\n", pGL->dist[i], pGL->radius[i], pGL->I[i][0], pGL->I[i][1], pGL->I[i][2]) ;
		}

		if(render_ghost_lines_and_rings(image, md) == 1)
		    return 1 ;

		md->n_ghost_lines-- ;
		continue ;
	    }
	    if(!strcasecmp("GHOST_RING", cmd)) {
		struct ghost_ring_struct *pGR = &md->ghost_rings[md->n_ghost_rings] ;
		float p_inv_circle ;
		if(md->n_ghost_rings == MAX_GHOST_RINGS) {
		    fprintf(stderr, "LENSFLARE: exceeded maximum allowed number of ghost rings (%d)\n", MAX_GHOST_RINGS) ;
		    return 1 ;
		}
		md->n_ghost_rings++ ;
		//get_next_line(fp, str, filename) ;
		sscanf(str, "%*s%d%f%f%f", &pGR->n_ghosts, &pGR->dist, &pGR->radius, &pGR->MAX_CA) ;
		pGR->dist *= md->scale_dist ;
		pGR->radius *= md->scale_radius ;
		//fprintf(stderr, "%f %f %f\n", pGR->dist, pGR->radius, pGR->MAX_CA) ;
		pGR->n_ghosts=0 ;

		int i=0 ;

		while(1) {
		    if(get_next_line(fp, str, filename) == NULL)
			break ;
		    if(!string_is_numeric(str)) {
			// set flag that the next line is already loaded into str
			str_is_filled=1 ;
			break ;
		    }

		    if(i == MAX_GHOSTS_IN_LINE) {
			fprintf(stderr, "LENSFLARE: exceeded maximum allowed number of ghosts in a line (%d)\n", MAX_GHOSTS_IN_LINE) ;
			return 1 ;
		    }

		    sscanf(str, "%f%f%f%f%f%f", &pGR->angle[i], &pGR->I[i][0], &pGR->I[i][1], &pGR->I[i][2],
			&p_inv_circle, &pGR->r_func[i]) ;
		    pGR->pcircle[i] = 1.-p_inv_circle ;
		    pGR->angle[i] += md->rotate ;
		    pGR->I[i][0] += md->scale_hsv[0] ;
		    pGR->I[i][1] *= md->scale_hsv[1] ;
		    pGR->I[i][2] *= md->scale_hsv[2] ;
		    hsv2rgb_1(pGR->I[i]) ;
		    pGR->n_ghosts++ ;
		    i++ ;
		    //fprintf(stderr, "GR: %f %f %f %f\n", pGR->angle[i], pGR->I[i][0], pGR->I[i][1], pGR->I[i][2]) ;
		}

		if(render_ghost_lines_and_rings(image, md) == 1)
		    return 1 ;

		md->n_ghost_rings-- ;
		continue ;
	    }

	    fprintf(stderr, "FATAL: unrecognnized line \'%s\' in commandfile %s\n", str, filename) ;
	    return 1 ;
	}

	fclose(fp) ;

	//render_ghost_lines_and_rings(image, md) ;

    } else {
	fprintf(stderr, "... Unable to open requested lens flare data file \'%s\', returning ...\n", filename) ;
	return 1 ;
    }

    return 0 ;
}

void render_lens_flare(char *filename, struct lens_flare_metadata *md, tdata_t *image)
{
    if(md->flare_type <= 2) {
	render_simple_lens_flare(filename, md, image) ;
	return ;
    }

    struct flare_struct {
	float A,B,C,D,E ;
	float radius ;
	float r,g,b ;
	float intensity ;
	float width ;
	int hue_based ;
	float hue_start ;
	float hue_speed ;
	float falloff_A ;
	float falloff_B ;
	float sina, cosa ;
	float angle ;
	int group ;
    } fs[180] ;

    int n_flares=0 ;

    struct group_data {
	float radius ;
	float r,g,b ;
	float intensity ;
	float width ;
	float hue_start ;
	float hue_speed ;
	float angle_delta ;
	float falloff_A ;
	float falloff_B ;
    } gd[20] ;

    int flares_per_group=1 ;

    float group_angle_delta=2.*M_PI/(float)md->n_leaves ;

    if(md->flare_type >= 3) {
	// RX100ii like
	group_angle_delta=2.*M_PI/(float)md->n_leaves ;

	flares_per_group=12 ;
	float angle_delta=group_angle_delta/(float)flares_per_group/16. ;
	float mid_i=(flares_per_group-1.)/2. ;
	for(int i=0 ; i < flares_per_group ; i++) {
	    gd[i].radius = md->flare_radius ;
	    float pedge = abs((float)i-mid_i)/mid_i ;
	    gd[i].r = .9 + pedge*.1 ;
	    gd[i].g = .9 + pedge*.08 ;
	    gd[i].b = .9 ;
	    gd[i].intensity = .125 ;
	    gd[i].width = 6.25 ;
	    gd[i].width = 12.5 ;
	    gd[i].hue_start = 90 ;
	    gd[i].hue_speed = 360 ;
	    gd[i].angle_delta = i > 0 ? angle_delta : 0. ;
	    gd[i].falloff_A = .65 ;
	    gd[i].falloff_B = 10 ;
	}
    }

    // add lens flare
    int group=0 ;
    for(float group_flare_angle=group_angle_delta*.4 ; group_flare_angle < 2.*M_PI-group_angle_delta/2. ; group_flare_angle += group_angle_delta) {
	group++ ;
	float angle=group_flare_angle ;
	fprintf(stderr, "Group flare angle=%f\n", group_flare_angle*180./M_PI) ;

	float start_angle=angle*180./M_PI ;
	float final_angle=angle*180./M_PI ;
	int first_flare=n_flares ;

	for(int g=0 ; g < flares_per_group ; g++) {
	    angle += gd[g].angle_delta + (float)rand()/(float)RAND_MAX*M_PI/180 ;
	    fs[n_flares].angle=angle*180./M_PI ;
	    fprintf(stderr, "   flare %d in group, angle=%f\n", g, angle*180./M_PI) ;
	    final_angle=angle*180./M_PI ;
	    float sina = sinf(angle) ;
	    float cosa = cosf(angle) ;

	    if(1) {
		float x1=0. ;
		float y1=0. ;
		float x2=1. ;
		float y2=0. ;
		float dx=x2-x1 ;
		float dy=y2-y1 ;
		fs[n_flares].sina = sina ;
		fs[n_flares].cosa = cosa ;
		fs[n_flares].A = -dy ;
		fs[n_flares].B = dx ;
		fs[n_flares].C = -(fs[n_flares].A*x1 + fs[n_flares].B*y1) ;
		float denom = sqrt(dx*dx+dy*dy) ;
		fs[n_flares].A /= denom ;
		fs[n_flares].B /= denom ;
		fs[n_flares].C /= denom ;
	    } else {
		// since we are working with a unit circle, this is shorter, faster
		float x2=1. ;
		float y2=0. ;
		fs[n_flares].A = -y2 ;
		fs[n_flares].B = x2 ;
		fs[n_flares].C = 0. ;
	    }

	    fs[n_flares].group = group ;
	    fs[n_flares].radius = gd[g].radius ;

	    uint16_t rgb[3] = { gd[g].r*MAX16f, gd[g].g*MAX16f, gd[g].b*MAX16f} ;
	    float hsv[3] ;
	    rgb2hsv16(rgb,hsv) ;
	    fprintf(stderr, "ray %d, hue=%f + %f, s=%f * %f, v=%f * %f\n", n_flares,
									    hsv[0], md->scale_hsv[0],
									    hsv[1], md->scale_hsv[1],
									    hsv[2], md->scale_hsv[2]) ;
	    hsv[0] += md->scale_hsv[0] ;
	    hsv[1] *= md->scale_hsv[1] ;
	    hsv[2] *= md->scale_hsv[2] ;
	    if(hsv[1] > 1.) hsv[1] = 1. ;
	    if(hsv[2] > 1.) hsv[2] = 1. ;
	    if(hsv[1] < 0.) hsv[1] = 0. ;
	    if(hsv[2] < 0.) hsv[2] = 0. ;
	    hsv2rgb16(hsv,rgb) ;

	    fs[n_flares].r = (float)rgb[0]/MAX16f ;
	    fs[n_flares].g = (float)rgb[1]/MAX16f ;
	    fs[n_flares].b = (float)rgb[2]/MAX16f ;

	    fs[n_flares].hue_based = 0 ;
	    fs[n_flares].intensity=gd[g].intensity ;
	    fs[n_flares].width=gd[g].width ;
	    fs[n_flares].hue_start=gd[g].hue_start ;
	    fs[n_flares].hue_speed=gd[g].hue_speed ;
	    fs[n_flares].falloff_A=gd[g].falloff_A ;
	    fs[n_flares].falloff_B=gd[g].falloff_B ;
	    n_flares++ ;
	    if(md->flare_reverse_radius > 0.) {
		fs[n_flares] = fs[n_flares-1] ;
		fs[n_flares].A = -fs[n_flares-1].A ;
		fs[n_flares].B = -fs[n_flares-1].B ;
		fs[n_flares].radius *= md->flare_reverse_radius ;
		n_flares++ ;
	    }
	}

	float mid_angle=(start_angle+final_angle)/2. ;
	for(int i=first_flare ; i < n_flares ; i++) {
	    fs[i].D = md->D_per_angle/fs[i].radius/100. * (fs[i].angle-mid_angle) ;
	    fs[i].E = md->E_per_angle * (fs[i].angle-mid_angle) ;
	}
    }

    fprintf(stderr, "Modeling %d flares\n", n_flares) ;

    int x0=md->sun_x-md->flare_radius ;
    if(x0 < 0) x0=0 ;
    if(x0 > IMAGE_WIDTH-1) x0=IMAGE_WIDTH-1 ;
    int x1=md->sun_x+md->flare_radius ;
    if(x1 < 0) x1=0 ;
    if(x1 > IMAGE_WIDTH-1) x1=IMAGE_WIDTH-1 ;

    int y0=md->sun_y-md->flare_radius ;
    if(y0 < 0) y0=0 ;
    if(y0 > IMAGE_HEIGHT-1) y0=IMAGE_HEIGHT-1 ;
    int y1=md->sun_y+md->flare_radius ;
    if(y1 < 0) y1=0 ;
    if(y1 > IMAGE_HEIGHT-1) y1=IMAGE_HEIGHT-1 ;

    for(int x=x0 ; x <= x1 ; x++) {
	int rx = x-md->sun_x ;
	for(int y=y0 ; y <= y1 ; y++) {
	    int ry = y-md->sun_y ;

	    float sum_r=0., sum_g=0., sum_b=0. ;
	    float sum_wgts=0. ;

	    // for adjusting the ray length, current set to constant of down
	    float length_x=0. ;
	    float length_y=1. ;

	    for(int i=0 ; i < n_flares ; i++)
	    {
		// need to rotate the point, the flares all are aligned on the positive X axis
		float px = -(float)rx*fs[i].cosa + (float)ry*fs[i].sina ;
		float py = -(float)rx*fs[i].sina - (float)ry*fs[i].cosa ;

		float point_A=px ;
		float point_B=py ;
		float mag=px*px+py*py ;
		point_A /= sqrt(mag) ;
		point_B /= sqrt(mag) ;

		// get the dot product of this flare vs the point
		float dp = point_A*fs[i].B - point_B*fs[i].A ; // remember, for flare, A=-dy, B=dx ;

		// dp is the cos of the angle between the vectors, so if the dp < 0, the point x,y must be on the
		// opposite side of the flare
		if(dp < 0.0) continue ;

		float dp_length = length_x*fs[i].B - length_y*fs[i].A ;

		float angle_effect = (dp_length +1)/2. * 0.5 + 0.5 ;
		float radius_at_angle = fs[i].radius * angle_effect ;

		// no angle effect on rays
		float width_at_angle = fs[i].width ;

		// FIXME -- here is where a quadratic formulation can make the rays curve.
		float d=fs[i].A*px + fs[i].B*py ;

		d = fs[i].D*px*px + fs[i].B*py + fs[i].E ;

		float hypot_sq = px*px + py*py ;
		//h**2 = r**2 + d**2 ;
		float radius = sqrt(hypot_sq-d*d) ;
		float pr = radius/radius_at_angle ;


		if(pr > 1.) continue ;

		float hw=(0.25+0.75*pr)*width_at_angle ;

		float I = exp(6.0-6.0*pr)/exp(6.0) ;

		if(radius_at_angle-radius <= hw) {
		    float p = 1.-(radius_at_angle-radius)/hw ;
		    hw = hw*sqrt(1.-p*p) ;
		}

		float f = hw-fabsf(d) ;


		if(f > 0.) {
		    float pf = 1.-f/hw ;
		    float pI=1./(1.+expf((pf-fs[i].falloff_A)*fs[i].falloff_B)) ;
		    if(f > 1.) f=1. ;
		    uint16_t flare_rgb[3] ;

		    if(fs[i].hue_based) {
			float hs_for_group=((fs[i].group % 2) > 0 ? 360 : 720) ;
			float hue = fs[i].hue_start+hs_for_group*pr ;
			if(hue >= 360.) hue -= 360. ;
			if(hue >= 360.) hue -= 360. ;

			float s=1./(1.+exp(-(pr-.4)*10.)) ;
			s=1. ;

			float hsv[3] = {hue, s, fs[i].intensity*I*f*pI} ;

			hsv[1] = s ;
			hsv2rgb16(hsv, flare_rgb) ;
		    } else {
			flare_rgb[0] = fs[i].r*MAX16f*fs[i].intensity*I*f*pI ;
			flare_rgb[1] = fs[i].g*MAX16f*fs[i].intensity*I*f*pI ;
			flare_rgb[2] = fs[i].b*MAX16f*fs[i].intensity*I*f*pI ;
		    }

		    sum_r += flare_rgb[0] ;
		    sum_g += flare_rgb[1] ;
		    sum_b += flare_rgb[2] ;
		    sum_wgts += f ;
		}
	    }

	    if(sum_wgts > 0.) {
		float r,g,b ;
		tif_get3c(image,x,y,r,g,b) ;

		r += sum_r ;
		g += sum_g ;
		b += sum_b ;

		fdither(r) ;
		fdither(g) ;
		fdither(b) ;

		tif_set3c(image,x,y,r,g,b) ;
	    }


	}
    }
}


int lens_flare(char *filename, tdata_t *image, int sun_x_pixel, int sun_y_pixel)
{
    struct lens_flare_metadata md ;

    // initialize metadata
    md.n_ghosts=0 ;
    md.n_ghost_lines=0 ;
    md.n_ghost_rings=0 ;
    md.n_leaves=7 ;

    md.flare_type = 1 ;
    md.flare_radius=650 ;
    md.flare_reverse_radius=0.5 ;
    md.sun_x = sun_x_pixel ;
    md.sun_y = sun_y_pixel ;

    md.TRUE_CENTER_X = IMAGE_WIDTH/2 ;
    md.TRUE_CENTER_Y = IMAGE_HEIGHT/2 ;


    md.scale_hsv[0] = 0. ;
    md.scale_hsv[1] = 1. ;
    md.scale_hsv[2] = 1. ;
    md.scale_radius=1. ;
    md.scale_dist=1. ;
    md.rotate=0. ;

    if(parse_lens_flare_file_and_render_ghosts(filename, &md, image) == 1)
	return 1 ;

    return 0 ;
}

void render_simple_lens_flare(char *filename, struct lens_flare_metadata *md, tdata_t *image)
{

    struct flare_struct {
	float A,B,C ;
	float radius ;
	float r,g,b ;
	float intensity ;
	float width ;
	float hue_start ;
	float hue_speed ;
	float falloff_A ;
	float falloff_B ;
	int group ;
    } fs[180] ;

    int n_flares=0 ;

    struct group_data {
	float radius ;
	float r,g,b ;
	float intensity ;
	float width ;
	float hue_start ;
	float hue_speed ;
	float angle_delta ;
	float falloff_A ;
	float falloff_B ;
    } gd[10] ;

    int flares_per_group=1 ;

    float group_angle_delta=2.*M_PI/(float)md->n_leaves ;

    if(md->flare_type == 2) {
	// RX100ii like
	group_angle_delta=2.*M_PI/(float)md->n_leaves ;

	flares_per_group=6 ;
	float angle_delta=group_angle_delta/(float)flares_per_group/8. ;

#define FLARE_W 12.5
	struct group_data gd0 = {md->flare_radius*.8, 1.0, 0.8, 0.8, .125, FLARE_W,  90, 360, 0, .65, 10} ;
	struct group_data gd1 = {md->flare_radius*.9, 0.5, 0.5, 1.0, .125, FLARE_W,  90, 360, angle_delta, .65, 10} ;
	struct group_data gd2 = {md->flare_radius, 0.5, 0.5, 1.0, .125, FLARE_W,  90, 360, angle_delta, .65, 10} ;
	struct group_data gd3 = {md->flare_radius, 0.5, 0.5, 1.0, .125, FLARE_W,  90, 360, angle_delta, .65, 10} ;
	struct group_data gd4 = {md->flare_radius*.9, 0.5, 0.5, 1.0, .125, FLARE_W,  90, 360, angle_delta, .65, 10} ;
	struct group_data gd5 = {md->flare_radius*.8, 0.5, 0.5, 1.0, .125, FLARE_W,  90, 360, angle_delta, .65, 10} ;
	gd[0] = gd0 ;
	gd[1] = gd1 ;
	gd[2] = gd2 ;
	gd[3] = gd3 ;
	gd[4] = gd4 ;
	gd[5] = gd5 ;

    } else {

	flares_per_group=1 ;
	float angle_delta=group_angle_delta/(float)flares_per_group/4. ;

	struct group_data gd0 = {md->flare_radius, 1.0, 0.8, 0.8, .25, 15,  90, 360, 0, .50, 5} ;
	struct group_data gd1 = {md->flare_radius, 1.0, 0.8, 0.8, .5, 45,  90, 360, 0, .65, 10} ;
	gd[0] = gd0 ;
	gd[1] = gd1 ;
    }

    for(int i=0 ; i < flares_per_group ; i++) {
	gd[i].intensity *= md->scale_hsv[2] ;
	if(gd[i].intensity > 1.) gd[i].intensity = 1. ;
	if(gd[i].intensity < 0.) gd[i].intensity = 0. ;
    }

    // add lens flare
    int group=0 ;
    for(float group_flare_angle=group_angle_delta*.4 ; group_flare_angle < 2.*M_PI-group_angle_delta/2. ; group_flare_angle += group_angle_delta) {
	group++ ;
	float angle=group_flare_angle ;
	fprintf(stderr, "Group flare angle=%f\n", group_flare_angle*180./M_PI) ;

	for(int g=0 ; g < flares_per_group ; g++) {
	    angle += gd[g].angle_delta + (float)rand()/(float)RAND_MAX*M_PI/90 ;
	    float sina = sinf(angle) ;
	    float cosa = cosf(angle) ;

	    if(1) {
		float x1=0. ;
		float y1=0. ;
		float x2=cosa ;
		float y2=sina ;
		float dx=x2-x1 ;
		float dy=y2-y1 ;
		fs[n_flares].A = -dy ;
		fs[n_flares].B = dx ;
		fs[n_flares].C = -(fs[n_flares].A*x1 + fs[n_flares].B*y1) ;
		float denom = sqrt(dx*dx+dy*dy) ;
		fs[n_flares].A /= denom ;
		fs[n_flares].B /= denom ;
		fs[n_flares].C /= denom ;
	    } else {
		// since we are working with a unit circle, this is shorter, faster
		float x2=cosa ;
		float y2=sina ;
		fs[n_flares].A = -y2 ;
		fs[n_flares].B = x2 ;
		fs[n_flares].C = 0. ;
	    }

	    // get the dot product of this flare vs the point
	    float dp = 0.*fs[n_flares].B - 1.*fs[n_flares].A ; // remember, for flare, A=-dy, B=dx ;

	    fs[n_flares].group = group ;
	    fs[n_flares].radius = gd[g].radius ;
	    fs[n_flares].r = gd[g].r ;
	    fs[n_flares].g = gd[g].g ;
	    fs[n_flares].b = gd[g].b ;
	    fs[n_flares].r = (float)rand()/(float)RAND_MAX ;
	    fs[n_flares].g = (float)rand()/(float)RAND_MAX ;
	    fs[n_flares].b = (float)rand()/(float)RAND_MAX ;
	    fs[n_flares].intensity=gd[g].intensity ;
	    fs[n_flares].width=gd[g].width ;
	    fs[n_flares].hue_start=gd[g].hue_start ;
	    fs[n_flares].hue_speed=gd[g].hue_speed ;
	    fs[n_flares].falloff_A=gd[g].falloff_A ;
	    fs[n_flares].falloff_B=gd[g].falloff_B ;
	    n_flares++ ;
	    if(md->flare_reverse_radius > 0.) {
		fs[n_flares] = fs[n_flares-1] ;
		fs[n_flares].A = -fs[n_flares-1].A ;
		fs[n_flares].B = -fs[n_flares-1].B ;
		fs[n_flares].radius *= md->flare_reverse_radius ;
		n_flares++ ;
	    }
	}
    }

    fprintf(stderr, "Modeling %d flares\n", n_flares) ;

    int x0=md->sun_x-md->flare_radius ;
    if(x0 < 0) x0=0 ;
    if(x0 > IMAGE_WIDTH-1) x0=IMAGE_WIDTH-1 ;
    int x1=md->sun_x+md->flare_radius ;
    if(x1 < 0) x1=0 ;
    if(x1 > IMAGE_WIDTH-1) x1=IMAGE_WIDTH-1 ;

    int y0=md->sun_y-md->flare_radius ;
    if(y0 < 0) y0=0 ;
    if(y0 > IMAGE_HEIGHT-1) y0=IMAGE_HEIGHT-1 ;
    int y1=md->sun_y+md->flare_radius ;
    if(y1 < 0) y1=0 ;
    if(y1 > IMAGE_HEIGHT-1) y1=IMAGE_HEIGHT-1 ;

    for(int x=x0 ; x <= x1 ; x++) {
	int rx = x-md->sun_x ;
	for(int y=y0 ; y <= y1 ; y++) {
	    int ry = y-md->sun_y ;
	    float point_A=rx ;
	    float point_B=ry ;
	    float mag=rx*rx+ry*ry ;
	    point_A /= sqrt(mag) ;
	    point_B /= sqrt(mag) ;

	    float sum_r=0., sum_g=0., sum_b=0. ;
	    float sum_wgts=0. ;

	    // for adjusting the ray length, current set to constant pointing down
	    float length_x=0. ;
	    float length_y=1. ;

#ifdef USE_CLOSEST_FLARE_ONLY
	    // find the closest flare of interest
	    float closest_d=1.e30 ;
	    int closest_i=0 ;

	    for(int i=0 ; i < n_flares ; i++) {
		// get the dot product of this flare vs the point
		float dp = point_A*fs[i].B - point_B*fs[i].A ; // remember, for flare, A=-dy, B=dx ;

		// dp is the cos of the angle between the vectors, so if the dp < 0, the point x,y must be on the
		// opposite side of the flare
		if(dp < 0.0) continue ;

		float dp_length = length_x*fs[i].B - length_y*fs[i].A ;

		float angle_effect = (dp_length +1)/2. * 0.5 + 0.5 ;
		float radius_at_angle = fs[i].radius * angle_effect ;

		// no angle effect on rays
		float width_at_angle = fs[i].width ;

		float d=fs[i].A*(float)rx + fs[i].B*(float)ry ;
		float hypot_sq = rx*rx + ry*ry ;
		//h**2 = r**2 + d**2 ;
		float radius = sqrt(hypot_sq-d*d) ;
		float pr = radius/radius_at_angle ;

		/*  		if(i == 0 && ry == 0) {  */
		/*  		    fprintf(stderr, "rx,ry:%d,%d d:%f hypotsq:%f radius:%f pr:%f\n", rx, ry, d, hypot_sq, radius, pr) ;  */
		/*  		}  */

		if(pr > 1.) continue ;
		if(d < 0.) d=-d ;

		if(d < closest_d) {
		    closest_d = d ;
		    closest_i = i ;
		}
	    }

	    //for(int i=closest_i ; i < closest_i+1 ; i++)
#endif
	    for(int i=0 ; i < n_flares ; i++)
	    {
		// get the dot product of this flare vs the point
		float dp = point_A*fs[i].B - point_B*fs[i].A ; // remember, for flare, A=-dy, B=dx ;

		// dp is the cos of the angle between the vectors, so if the dp < 0, the point x,y must be on the
		// opposite side of the flare
		if(dp < 0.0) continue ;

		float dp_length = length_x*fs[i].B - length_y*fs[i].A ;

		float angle_effect = (dp_length +1)/2. * 0.5 + 0.5 ;
		float radius_at_angle = fs[i].radius * angle_effect ;

		// no angle effect on rays
		float width_at_angle = fs[i].width ;

		float d=fs[i].A*(float)rx + fs[i].B*(float)ry ;
		float hypot_sq = rx*rx + ry*ry ;
		//h**2 = r**2 + d**2 ;
		float radius = sqrt(hypot_sq-d*d) ;
		float pr = radius/radius_at_angle ;

		/*  		if(i == 0 && ry == 0) {  */
		/*  		    fprintf(stderr, "rx,ry:%d,%d d:%f hypotsq:%f radius:%f pr:%f\n", rx, ry, d, hypot_sq, radius, pr) ;  */
		/*  		}  */

		if(pr > 1.) continue ;

/*  		    float I=exp(-pr*pr*60.)*0.25+exp(-pr*3.)*0.5 ;  */
/*  		    I = powf(1.-pr,0.75) ;  */
/*  		    I=exp(-pr*3.)-exp(-3.) ;  */

		//float I=exp(-pr*pr*90.)*0.75+(exp(-pr*2.) - exp(-2.))*.25;
		//I=1 + pr*(-4.126984 + pr*(6.407407 + pr*(-3.280423))) ;
		float hw=(0.25+0.75*pr)*width_at_angle ;

		// max I = 1./0.25 * (1.-0) = 4.
		//I = (1./(0.25+0.75*pr)) / 4. * exp(-4.*pr);
		//I = exp(-4.*pr);
		float I = exp(6.0-6.0*pr)/exp(6.0) ;

		if(radius_at_angle-radius <= hw) {
		    float p = 1.-(radius_at_angle-radius)/hw ;
		    hw = hw*sqrt(1.-p*p) ;
		}

		float f = hw-fabsf(d) ;


		if(f > 0.) {
		    float pf = 1.-f/hw ;
		    float pI=1./(1.+expf((pf-fs[i].falloff_A)*fs[i].falloff_B)) ;
		    if(f > 1.) f=1. ;
		    float hs_for_group=((fs[i].group % 2) > 0 ? 360 : 720) ;
		    float hue = fs[i].hue_start+hs_for_group*pr ;
		    if(hue >= 360.) hue -= 360. ;
		    if(hue >= 360.) hue -= 360. ;

		    float s=1./(1.+exp(-(pr-.4)*10.)) ;

		    s *= md->scale_hsv[1] ;
		    if(s < 0.) s=0. ;
		    if(s > 1.) s=1. ;
		    float hsv[3] = {hue, s, fs[i].intensity*I*f*pI} ;

		    uint16_t flare_rgb[3] ;
		    hsv2rgb16(hsv, flare_rgb) ;
		    sum_r += flare_rgb[0] ;
		    sum_g += flare_rgb[1] ;
		    sum_b += flare_rgb[2] ;
		    sum_wgts += f ;
		}
	    }

	    if(sum_wgts > 0.) {
		float r,g,b ;
		tif_get3c(image,x,y,r,g,b) ;

		r += sum_r ;
		g += sum_g ;
		b += sum_b ;
		if(r > MAX16f) r=MAX16f ;
		if(g > MAX16f) g=MAX16f ;
		if(b > MAX16f) b=MAX16f ;

		fdither(r) ;
		fdither(g) ;
		fdither(b) ;

		tif_set3c(image,x,y,r,g,b) ;
	    }


	}
    }
}

#ifdef LENSFLARE_MAIN
int main(int argc, char* argv[])
{

    char outfilename[512] ;
    char lensflarefilename[512] ;
    lensflarefilename[0] = 0 ; // this one needs initialization here


    if(argc < 4) {
	fprintf(stderr, "%s <output_TIF_file> -LF <lensflare_command_file>", argv[0]) ;
	exit(1) ;
    }

    strcpy(outfilename, argv[1]) ;

    if(!strcmp(argv[2], "-LF")) {
	strcpy(lensflarefilename, argv[3]) ;
	argc -= 2 ;
	argv += 2 ;
    } else {
	fprintf(stderr, "%s <output_TIF_file> -LF <lensflare_command_file>", argv[0]) ;
	exit(1) ;
    }

    // Read the input tif image
    int32_t x, y ;
    int input_W=800 ;
    int input_H=800 ;
    int is_16_bit_image=1 ;
    int BPS=8+8*is_16_bit_image ;
    IMAGE_NSAMPLES=3 ;
    int SPP=IMAGE_NSAMPLES ;
    tdata_t *image ;
    tsize_t ROWSIZE ;
    ROWSIZE = input_W*BPS/8*SPP ;
    IMAGE_HAS_ALPHA = IMAGE_NSAMPLES > 3 ? 1 : 0 ;

    IMAGE_HEIGHT = (int32_t)input_H ; // set global variable
    IMAGE_WIDTH = (int32_t)input_W ; // set global variable

    if((int)ROWSIZE != (int)(input_W*BPS/8*SPP)) {
	fprintf(stderr, "ROWSIZE is WRONG!, fixing...\n") ;
	ROWSIZE = input_W*BPS/8*SPP ;
    }

    if(is_16_bit_image) {
	fprintf(stderr, "16 bit image, ROWSIZE is now %d\n", (int)ROWSIZE) ;
    } else {
	ROWSIZE *= 2 ;
	fprintf(stderr, "8 bit image, ROWSIZE is now %d\n", (int)ROWSIZE) ;
    }

    /* allocate memory for image */
    image = (tdata_t *)calloc(input_H, sizeof(tdata_t *)) ;

    for(int y = 0 ; y < input_H ; y++) {
	image[y] = (tdata_t) _TIFFmalloc(ROWSIZE);
	if(image[y] == NULL) {
	    fprintf(stderr, "Failed to alloc memory for image at row %d\n", y) ;
	    exit(1) ;
	}
    }

    for(x = 0 ; x < IMAGE_WIDTH  ; x++) {
	for(y = 0 ; y < IMAGE_HEIGHT ; y++) {

	    if(IMAGE_HAS_ALPHA) {
		tif_set4c(image,x,y,0,0,0,1) ;
	    } else {
		tif_set3c(image,x,y,0,0,0) ;
	    }

	}
    }

    int sun_x = IMAGE_WIDTH/2 ;
    int sun_y = IMAGE_HEIGHT/2 ;

    // render lens flare if requested
    if(lensflarefilename[0] != 0) {
	lens_flare(lensflarefilename ,image, sun_x, sun_y) ;
    }

    TIFF* tifout = NULL ; // may not be used if just refitting models for analysis
    tifout = TIFFOpen(outfilename, "w");
    uint16_t TIF_CONFIG = PLANARCONFIG_CONTIG ;

    // prepare output image to have same tags as input image
    uint16_t alltags[] = {
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



/*  setting i:0, tag number 256 to 2810, in output file  */
/*  setting i:1, tag number 257 to 867, in output file  */
/*  setting i:2, tag number 258 to 16, in output file  */
/*  setting i:3, tag number 277 to 4, in output file  */
/*  setting i:4, tag number 278 to 93, in output file  */
/*  setting i:5, tag number 274 to 1, in output file  */
/*  setting i:6, tag number 284 to 1, in output file  */
/*  setting i:7, tag number 262 to 2, in output file  */
/*  setting i:8, tag number 339 to 1, in output file  */
/*  setting i:9, tag number 259 to 5, in output file  */


    TIFFSetField(tifout, TIFFTAG_IMAGEWIDTH, (uint32_t)IMAGE_WIDTH);
    TIFFSetField(tifout, TIFFTAG_IMAGELENGTH, (uint32_t)IMAGE_HEIGHT);
    TIFFSetField(tifout, TIFFTAG_BITSPERSAMPLE, (uint32_t)BPS);
    TIFFSetField(tifout, TIFFTAG_SAMPLESPERPIXEL, (uint32_t)SPP);
/*      TIFFSetField(tifout, TIFFTAG_ROWSPERSTRIP, (uint32_t)93);  */
    TIFFSetField(tifout, TIFFTAG_ROWSPERSTRIP, TIFFDefaultStripSize(tifout,IMAGE_WIDTH*SPP));
    TIFFSetField(tifout, TIFFTAG_ORIENTATION, (uint32_t)1);
    TIFFSetField(tifout, TIFFTAG_PLANARCONFIG, (uint32_t)1);

    TIFFSetField(tifout, TIFFTAG_PHOTOMETRIC, (uint32_t)2);
    TIFFSetField(tifout, TIFFTAG_SAMPLEFORMAT, (uint32_t)1);
    TIFFSetField(tifout, TIFFTAG_COMPRESSION, (uint32_t)5);

    /* save image data */

    if(TIF_CONFIG == PLANARCONFIG_CONTIG) {
	printf("TIF samples stored contiguously\n") ;
    } else {
	printf("TIF samples stored separately\n") ;
    }

    uint16_t r,g,b ;
    tif_get3c(image,100,100,r,g,b) ;


    printf("value of 100,100  is %d,%d,%d\n", r, g, b) ;


    for(y = 0 ; y < IMAGE_HEIGHT ; y++) {
	if(TIF_CONFIG == PLANARCONFIG_CONTIG) {
	    if(TIFFWriteScanline(tifout,image[y],y,0) < 0) {
		printf("Error writing contigous row %d\n", y) ;
		exit(1) ;
	    }
	} else if(TIF_CONFIG == PLANARCONFIG_SEPARATE) {
	    uint16_t s ;
	    for(s = 0 ; s < SPP ; s++) {
		if(TIFFWriteScanline(tifout,image[y],y,s) < 0) {
		    printf("Error writing separate row %d, component %d\n", y, s) ;
		    exit(1) ;
		}
	    }
	} else {
	    fprintf(stderr, "Cannot write planar configuration of tiff image, config=%d!\n", (int)TIF_CONFIG) ;
	    exit(1) ;
	}
    }


    TIFFClose(tifout);

    fprintf(stderr, "free-ing image memory\n") ;

    /* free memory for image */

    for(y = 0 ; y < IMAGE_HEIGHT ; y++) {
	_TIFFfree(image[y]) ;
    }

}
#endif
