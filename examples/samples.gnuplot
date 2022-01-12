set title 'test'
unset hidden3d
set ticslevel 0.5
#set view 0,0
set view map
set autoscale
set parametric
set style data points
# 1  2 3 4 5 6    7    8    9    10   11   12    13        14
#px py h s v hhat shat vhat herr serr verr angle sumdistsv sumdistvy
set key box
set terminal pngcairo size 800,600 enhanced font 'Verdana,10'

#set title 'pred sat vs hue, colored by py'
#set output 'samples_s3val.png'
#splot "samples.dat" u 3:7:2 with points palette pointsize 1 pointtype 7

set title 'sumdist2 SV vs sumdist2 VY'
set output 'samples_outlier.png'
plot "samples.dat" u 14:13 with points pointsize 1 pointtype 7

set title 'Actual Saturation vs Value, colored by py'
set output 'samples_s_vs_I.png'
splot "samples.dat" u 5:4:2 with points palette pointsize 1 pointtype 7

set title 'hue err vs predicted hue, colored by py'
set output 'samples_h3err.png'
splot "samples.dat" u 6:9:2 with points palette pointsize 1 pointtype 7

set title 'sat err vs predicted sat, colored by py'
set output 'samples_s3err.png'
splot "samples.dat" u 7:10:2 with points palette pointsize 1 pointtype 7

set title 'val err vs py, colored by predicted val'
set output 'samples_v3err.png'
splot "samples.dat" u 2:11:8 with points palette pointsize 1 pointtype 7

set yrange [120:230]
set title 'Hue vs py, colored by px'
set output 'samples_h3pya.png'
splot "samples.dat" u 2:3:1 with points palette pointsize 1 pointtype 7

set title 'Hue vs horizontal angle, colored by py'
set output 'samples_h3pxa.png'
splot "samples.dat" u 12:3:2 with points palette pointsize 1 pointtype 7

set title 'Predicted Hue vs py, colored by px'
set output 'samples_h3pyp.png'
splot "samples.dat" u 2:6:1 with points palette pointsize 1 pointtype 7

set title 'Predicted Hue vs horizontal angle, colored by py'
set output 'samples_h3pxp.png'
splot "samples.dat" u 12:6:2 with points palette pointsize 1 pointtype 7

#set title 'Predicted Hue0 vs py, colored by px'
#set output 'samples_h3pyp0.png'
#splot "samples.dat" u 2:13:1 with points palette pointsize 1 pointtype 7

#set title 'Predicted Hue1 vs py, colored by px'
#set output 'samples_h3pyp1.png'
#splot "samples.dat" u 2:14:1 with points palette pointsize 1 pointtype 7


set yrange [0.0:0.7]
set title 'Saturation vs py, colored by px'
set output 'samples_s3pya.png'
splot "samples.dat" u 2:4:1 with points palette pointsize 1 pointtype 7

set title 'Saturation vs horizontal angle, colored by py'
set output 'samples_s3pxa.png'
splot "samples.dat" u 12:4:2 with points palette pointsize 1 pointtype 7

set title 'Predicted Saturation vs py, colored by px'
set output 'samples_s3pyp.png'
splot "samples.dat" u 2:7:1 with points palette pointsize 1 pointtype 7

set title 'Predicted Saturation vs horizontal angle, colored by py'
set output 'samples_s3pxp.png'
splot "samples.dat" u 12:7:2 with points palette pointsize 1 pointtype 7



set yrange [0.0:1.0]
set title 'py vs px, colored by value'
set output 'samples_v3pypx.png'
splot "samples.dat" u 1:2:5 with points palette pointsize 1 pointtype 7

set title 'Value vs py, colored by px'
set output 'samples_v3pya.png'
splot "samples.dat" u 2:5:1 with points palette pointsize 1 pointtype 7

set title 'Value vs horizontal angle, colored by py'
set output 'samples_v3pxa.png'
splot "samples.dat" u 12:5:2 with points palette pointsize 1 pointtype 7

set title 'Predicted Value vs py, colored by px'
set output 'samples_v3pyp.png'
splot "samples.dat" u 2:8:1 with points palette pointsize 1 pointtype 7

set title 'Predicted Value vs horizontal angle, colored by py'
set output 'samples_v3pxp.png'
splot "samples.dat" u 12:8:2 with points palette pointsize 1 pointtype 7
