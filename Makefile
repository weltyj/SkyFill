CFLAGS=-Wall -g
SRC=skyfill_tif.c repair_sky.c amoeba_06.c mstat.c colorspace_conversions.c pixel_tests.c feather_factor.c find_sky.c estimate_sky.c sample_and_fit_sky_model.c optimize.c
OBJS=$(SRC:.c=.o)
LIBS=-ltiff -lm

skyfill_tif : $(OBJS) Makefile
	gcc $(CFLAGS) $(OBJS) $(LIBS) -o skyfill_tif

$(OBJS) : skyfill_tif.h sample_and_fit_sky_model.h feather_factor.h Makefile colorspace_conversions.h
 
.c.o:
	gcc -c $(CFLAGS) $< -o $@

backup:
	tar cvzf ../skyfill.tgz *.c *.h Makefile CMakeLists.txt
