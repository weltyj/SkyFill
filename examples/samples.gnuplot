set title 'test'
unset hidden3d
set ticslevel 0.5
#set view 0,0
set view map
set autoscale
set parametric
set style data points
#1 2 3 4 5 6    7    8    9    10   11   12
#x y h s v hhat shat vhat herr serr verr angle
set key box
set terminal pngcairo size 800,600 enhanced font 'Verdana,10'

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


set title 'Saturation vs py, colored by px'
set output 'samples_s3pya.png'
splot "samples.dat" u 2:4:1 with points palette pointsize 1 pointtype 7

set title 'Saturation vs horizontal angle, colored by py'
set output 'samples_s3pya.png'
splot "samples.dat" u 12:4:2 with points palette pointsize 1 pointtype 7

set title 'Predicted Saturation vs py, colored by px'
set output 'samples_s3pyp.png'
splot "samples.dat" u 2:7:1 with points palette pointsize 1 pointtype 7

set title 'Predicted Saturation vs horizontal angle, colored by py'
set output 'samples_s3pyp.png'
splot "samples.dat" u 12:7:2 with points palette pointsize 1 pointtype 7



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

#set output 'samples_h0.png'
#plot "samples.dat" u 1:3 with circles lc rgb "red"
#set output 'samples_s0.png'
#plot "samples.dat" u 1:4 with circles lc rgb "red"
#set output 'samples_v0.png'
#plot "samples.dat" u 1:5 with circles lc rgb "red"

#set output 'samples_h1hat.png'
#plot "samples.dat" u 1:6 with circles lc rgb "red"
#set output 'samples_s1hat.png'
#plot "samples.dat" u 1:7 with circles lc rgb "red"
#set output 'samples_v1hat.png'
#plot "samples.dat" u 1:8 with circles lc rgb "red"

#set output 'samples_h2err.png'
#plot "samples.dat" u 6:9 with circles lc rgb "red"
#set output 'samples_s2err.png'
#plot "samples.dat" u 7:10 with circles lc rgb "red"
#set output 'samples_v2err.png'
#plot "samples.dat" u 8:11 with circles lc rgb "red"

