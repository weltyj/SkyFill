FHML="-fhml 2"
#FLAGS=-SR
FLAGS=
SF=../git/SkyFill/skyfill_tif

#$SF 2021_GC_lunch.tif 2021_GC_lunchsf2.tif $FHML $FLAGS $1

#$SF 2021_pan02.tif 2021_pan02sf.tif -fsr $FHML -CIE 6 -df 1.0 -m 0 493 -rm 0 1626 -tol 360 .2 .15 -C 2 -E 0.0 -sx -5 -sy .87 -O sx sy F G $FLAGS -sl 1 -D -18 -F .64 -G .89 $1 $2

#$SF 2021_pan19.tif 2021_pan19sf.tif $FHML -fsh 0.8 1.0 -df .9 -tol 360 .3 .05 $FLAGS $1

#$SF 2021_pan24.tif 2021_pan24sf.tif -df .5 -CIE 6 $FHML -sx 60 -sy 1.05 -sl .2 -D -10 -E 0.0 -O sd C D $FLAGS $1 $2
#$SF 2021_pan24.tif 2021_pan24sf_fsr.tif -fsr -df .5 -CIE 6 $FHML -sx 60 -sy 1.05 -sl .2 -D -10 -E 0.0 -O sd C D $FLAGS $1 $2

#$SF 2021_pan37.tif 2021_pan37sf.tif -fsr -hy .48 -CIE 6 -df 1.0 -nf 0.5 -sx -55 -sy .825 -C .29 -D -9.3 -F 3.5 -G 1.0 -sm 300 1000 -rm 457 794 -tol 5 .3 .08 $FHML $FLAGS $1 $2

#$SF 2021_ShoshonePt.tif 2021_ShoshonePtsf2.tif $FHML -tol 360 .3 .03 -O sd -fs 0.7 -df 0.8 $FLAGS $1 

#$SF DSC03507-DSC03512.tif DSC03507.tif $FHML -fsh 0.65 1.0 -tol 360 .3 .03 $FLAGS $1

$SF meditation_rock.tif meditation_rocksf.tif $FHML $FLAGS $1
