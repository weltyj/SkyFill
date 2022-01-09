FHML="-fhml 2"
#FLAGS=-SR
FLAGS=
#SF=../git/SkyFill/skyfill_tif
SF=../git/merge/skyfill_tif

#$SF 2021_GC_lunch.tif 2021_GC_lunchsf2.tif $FHML $FLAGS $1

$SF 2021_pan02_v2.tif 2021_pan02_v2sf.tif -fsr -df 1.0 -nceos $FHML -hy .35 -rm 0 1626 -tm 490 799 0 284 -tm 899 948 398 429 -tm 1330 1364 406 430 -tm 1024 1600 0 339 -m 0 493 -nsrm 0 493 -tol 10 .1 .05 -D 2.8 -F 1.0 -G 1000 -sx -2.45 -sy .87 -O sd $FLAGS -sl 1  $1 -fsrt .95 10

#$SF 2021_pan19.tif 2021_pan19sf.tif -fsr $FHML -fsh 0.8 1.0 -df .9 -tol 360 .3 .05 $FLAGS $1 -fsrt .8 10

#$SF 2021_pan24.tif 2021_pan24sf.tif -tol 360 .3 .03 -fsr -nceos -df .5 $FHML -sx 60 -sy 1.05 -sl 1. -D -10 -O sd C D $FLAGS $1 $2

#$SF 2021_pan37.tif 2021_pan37sf.tif -fsr -hy .48 -fsrt .9 20 -df 0.5 -nf 1.0 -sx -55 -sy .825 -D 4 -F 1.0 -G 1000 -O sl sd -sm 300 1000 -rm 457 794 -tol 5 .3 .08 $FHML $FLAGS $1 $2 $3

#  This looks better without the -fsr
#$SF 2021_ShoshonePt.tif 2021_ShoshonePtsf2.tif $FHML -hy .65 -tol 360 .3 .03 -O sd -fs 1.0 -df 0.9 $FLAGS $1 

#$SF bruno_test16.tif bruno_test16sf.tif -tm 429 443 24 62 -fsr -nsrm 425 727 -df 0.1 -ff 0.0 -fsrt .95 100 -hy .4 -tol 1 .02 .02 -m 608 727 $FHML $FLAGS $1

#$SF DSC03507-DSC03512.tif DSC03507-DSC03512sf.tif -fsr $FHML -hy .42 -fsh 0.65 1.0 -tol 360 .3 .03 $FLAGS $1

#$SF meditation_rock.tif meditation_rocksf.tif -hy .4 $FHML $FLAGS $1

gnuplot < samples.gnuplot
