skyfill_tif : skyfill_tif.c amoeba_06.c mstat.c Makefile
	#gcc -Wall -O5 skyfill_tif.c amoeba_06.c mstat.c -ltiff -lm -o skyfill_tif
	gcc -Wall -g skyfill_tif.c amoeba_06.c mstat.c -ltiff -lm -o skyfill_tif
