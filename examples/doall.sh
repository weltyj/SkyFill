FHML="-fhml 2"
#FLAGS=-SR
FLAGS=
SF=../git/SkyFill/skyfill_tif

#$SF 2021_GC_lunch.tif 2021_GC_lunchsf2.tif -fsr $FHML $FLAGS $1

#$SF 2021_pan02.tif 2021_pan02sf.tif -fsr $FHML -CIE 6 -hy .4 -df 1.0 -m 0 493 -rm 0 1626 -tol 360 .2 .15 -C 2 -E 0.0 -sx -5 -sy .87 -O sx sy F G $FLAGS -sl 1 -D -18 -F .64 -G .89 $1 $2

#$SF 2021_pan19.tif 2021_pan19sf.tif -fsr $FHML -fsh 0.8 1.0 -df .9 -tol 360 .3 .05 $FLAGS $1

#$SF 2021_pan24.tif 2021_pan24sf.tif -fsr -df .5 -CIE 6 $FHML -sx 60 -sy 1.05 -sl .2 -D -10 -E 0.0 -O sd C D $FLAGS $1 $2

#$SF 2021_pan37.tif 2021_pan37sf.tif -fsr -hy .48 -CIE 6 -df 1.0 -nf 0.5 -sx -55 -sy .825 -C .29 -D -9.3 -F 3.5 -G 1.0 -O sl sd -sm 300 1000 -rm 457 794 -tol 5 .3 .08 $FHML $FLAGS $1 $2

#  This looks better without the -fsr
#$SF 2021_ShoshonePt.tif 2021_ShoshonePtsf2.tif $FHML -hy .65 -tol 360 .3 .03 -O sd -fs 1.0 -df 0.9 $FLAGS $1 

#$SF DSC03507-DSC03512.tif DSC03507-DSC03512sf.tif $FHML -hy .42 -fsh 0.65 1.0 -tol 360 .3 .03 $FLAGS $1

#$SF meditation_rock.tif meditation_rocksf.tif -fsr -hy .4 $FHML $FLAGS $1
