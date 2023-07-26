CFLAGS=-Wall -O3 -Wno-unused-variable
#CFLAGS=-Wall -g -Wno-unused-variable
SRC=skyfill_tif.c repair_sky.c amoeba_06.c mstat.c colorspace_conversions.c pixel_tests.c feather_factor.c find_sky.c estimate_sky.c sample_and_fit_sky_model.c optimize.c kmeans_rgb.c mpfit.c lensflare.c
OBJS=$(SRC:.c=.o)
LIBS=-ltiff -lm

all: skyfill_tif lensflare

skyfill_tif : $(OBJS) Makefile
	gcc $(CFLAGS) $(OBJS) $(LIBS) -o skyfill_tif

lensflare : lensflare.c colorspace_conversions.o Makefile
	gcc $(CFLAGS) -DLENSFLARE_MAIN lensflare.c colorspace_conversions.o $(LIBS) -o lensflare

$(OBJS) : skyfill_tif.h sample_and_fit_sky_model.h feather_factor.h Makefile colorspace_conversions.h estimate_sky.h kmeans_rgb.h mpfit.h
 
.c.o:
	gcc -c $(CFLAGS) $< -o $@

backup:
	tar cvzf ../skyfill.tgz *.c *.h Makefile CMakeLists.txt
