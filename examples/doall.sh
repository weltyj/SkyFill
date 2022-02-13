
GDB=""

# uncomment next line to run debugger
#GDB="gdb --args"

rm eos_test.prms
#echo "2 10 .01" > eos_test.prms

EXE=../git/SkyFill/skyfill_tif
#EXE=../git/new_CIE_test/skyfill_tif
#EXE=../git/Jan15_before_spline_curves/skyfill_tif

SF="$GDB $EXE"

#echo "6 -1 3.5 200 4" > eos_test.prms
#$SF 2021_GC_lunch.tif 2021_GC_lunchsf2.tif $SUN $1

MASK="-rm 0 1626 -tm 490 799 0 284 -tm 899 948 398 429 -tm 1330 1364 406 430 -tm 1024 1600 0 339 -sm 0 493 -nsrm 0 493"
SUN="-C 6.5 -D 1.0 -F 1.0 -G 0.12 -sx -2.45 -sy .87"
#$SF 2021_pan02_v2.tif 2021_pan02_v2sf.tif -fsrt .95 10 -fsr -df 0.5 $MASK $SUN $1

# a version to try full CIE2003 modeling, not as good...
#$SF 2021_pan02_v2.tif 2021_pan02_v2sf.tif -fsrt .95 10 -fsr -df 0.5 $MASK $SUN -CIE 12 $1

#$SF 2021_pan19.tif 2021_pan19sf.tif -fsr -fsh 0.8 1.0 $1

SUN="-sx 50 -sy 1.01 -D 3.5 -C 2 -F 1.0 -G .1 -sl .2"
#$SF 2021_pan24.tif 2021_pan24sf.tif -fsr -nceos -df .5 $SUN $1

MASK=" -tm 480 800 0 250 -sm 300 1000 -rm 457 794"
SUN="-sx -46.9 -sy .83 -D 3.5 -F 1.0 -G .1"
$SF 2021_pan37.tif 2021_pan37sf.tif -fsr $SUN $MASK $1

#  This looks better without the -fsr
#  -msu 0 => make sky smooth, after all processing run a smoother from the highest end sky to top of image
#         effect is 100% at the top of the image, and 0% a the highest end of sky
#$SF 2021_ShoshonePt.tif 2021_ShoshonePtsf2.tif -msu 0 -r_tos_thresh .15 -rsv 0.0 -sd 1.0      -fs 1.0 -df 0.9 $1 
#$SF test8.tif test8sf.tif -msu 0 -r_tos_thresh .15 -rsv 0.0 -sd 1.0      -fs 1.0 -df 0.9 $1 

MASK="-tm 427 443 23 62 -tm 166 174 254 266 -tm 337 366 179 214 -nsrm 0 727 -rm 580 727 -m 580 727"
#$SF bruno_test16.tif bruno_test16sf.tif -fov 80 -fsr $MASK $1

#$SF DSC03507-DSC03512.tif DSC03507-DSC03512sf.tif -ef .9 -fsr -hy .42 -fsh 0.65 1.0 $1

#echo "6 -1 5. 200 4" > eos_test.prms
#$SF LipanPoint2.tif LipanPoint2sf.tif -fsr -df 0.75 -rm 0 10000 $1

#echo "4 1457 0.999 3.9" > eos_test.prms
#echo "6 -1 3.5 200 4" > eos_test.prms
#$SF meditation_rock.tif meditation_rocksf.tif -hy .4 $1


#rm *.png
#~/bin/skyfill_plot.pl
