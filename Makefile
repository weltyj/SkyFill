CFLAGS=-Wall -g
SRC=skyfill_tif.c repair_sky.c amoeba_06.c mstat.c colorspace_conversions.c pixel_tests.c feather_factor.c find_sky.c estimate_sky.c sample_and_fit_sky_model.c optimize.c
OBJS=$(SRC:.c=.o)
LIBS=-ltiff -lm

skyfill_tif : $(OBJS) Makefile sample_and_fit_sky_model.h
	gcc $(CFLAGS) $(OBJS) $(LIBS) -o skyfill_tif
 
.c.o: Makefile sample_and_fit_sky_model.h
	gcc -c $(CFLAGS) $< -o $@

backup:
	tar cvzf ../skyfill.tgz *.c *.h Makefile CMakeLists.txt
