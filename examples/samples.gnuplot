set title 'test'
unset hidden3d
set ticslevel 0.5
#set view 0,0
set view map
set autoscale
set parametric
set style data points

#column names
#x y ha va px py hue sat val nv hue_hat sat_hat val_hat hue_err sat_err val_err cust1 cust2

# old column numbers
# 1  2 3 4 5 6    7    8    9    10   11   12     13        14        15            16
#px py h s v hhat shat vhat herr serr verr hangle sumdistsv sumdistvy normalized_v  vangle

set key box
set terminal pngcairo size 800,600 enhanced font 'Verdana,10'

#set title 'pred sat vs hue, colored by py'
#set output 'samples_s3val.png'
#splot "samples.dat" u 3:7:2 with points palette pointsize 1 pointtype 7

set title 'sumdist2 SV vs sumdist2 VY'
set output 'samples_outlier.png'
plot "samples.dat" u (column("cust2")):(column("cust1")) with points pointsize 1 pointtype 7

set title 'Actual Saturation vs Value, colored by vertical angle'
set output 'samples_s_vs_I.png'
splot "samples.dat" u (column("val")):(column("sat")):(column("py")):(column("gpsymbol")) with points palette pointsize 1 pointtype variable

set title 'hue err vs predicted hue, colored by vertical angle'
set output 'samples_h3err.png'
splot "samples.dat" u (column("hue_hat")):(column("hue_err")):(column("py")):(column("gpsymbol")) with points palette pointsize 1 pointtype variable

set title 'sat err vs predicted sat, colored by vertical angle'
set output 'samples_s3err.png'
splot "samples.dat" u (column("sat_hat")):(column("sat_err")):(column("py")):(column("gpsymbol")) with points palette pointsize 1 pointtype variable

set title 'val err vs va, colored by predicted val'
set output 'samples_v3err.png'
splot "samples.dat" u (column("py")):(column("val_err")):(column("val_hat")):(column("gpsymbol")) with points palette pointsize 1 pointtype variable

set yrange [120:230]
set title 'Hue vs va, colored by px'
set output 'samples_h3pya.png'
splot "samples.dat" u (column("py")):(column("hue")):(column("px")):(column("gpsymbol")) with points palette pointsize 1 pointtype variable

set title 'Hue vs horizontal angle, colored by va'
set output 'samples_h3pxa.png'
splot "samples.dat" u (column("ha")):(column("hue")):(column("py")):(column("gpsymbol")) with points palette pointsize 1 pointtype variable

set title 'Predicted Hue vs va, colored by px'
set output 'samples_h3pyp.png'
splot "samples.dat" u (column("py")):(column("hue_hat")):(column("px")):(column("gpsymbol")) with points palette pointsize 1 pointtype variable

set title 'Predicted Hue vs horizontal angle, colored by va'
set output 'samples_h3pxp.png'
splot "samples.dat" u (column("ha")):(column("hue_hat")):(column("py")):(column("gpsymbol")) with points palette pointsize 1 pointtype variable

#set title 'Predicted Hue0 vs py, colored by px'
#set output 'samples_h3pyp0.png'
#splot "samples.dat" u 2:13:1 with points palette pointsize 1 pointtype 7

#set title 'Predicted Hue1 vs py, colored by px'
#set output 'samples_h3pyp1.png'
#splot "samples.dat" u 2:14:1 with points palette pointsize 1 pointtype 7


set yrange [0.0:0.7]
set title 'Saturation vs va, colored by px'
set output 'samples_s3pya.png'
splot "samples.dat" u (column("py")):(column("sat")):(column("px")):(column("gpsymbol")) with points palette pointsize 1 pointtype variable

set title 'Normalized Saturation vs horizontal angle'
set output 'samples_s3pxn.png'
splot "normalized_samples.dat" u (column("ha")):(column("ns")):(column("py")):(column("gpsymbol")) with points palette pointsize 1 pointtype variable

set title 'Saturation vs horizontal angle, colored by va'
set output 'samples_s3pxa.png'
splot "samples.dat" u (column("px")):(column("sat")):(column("py")):(column("gpsymbol")) with points palette pointsize 1 pointtype variable

set title 'Predicted Saturation vs va, colored by px'
set output 'samples_s3pyp.png'
splot "samples.dat" u (column("py")):(column("sat_hat")):(column("px")):(column("gpsymbol")) with points palette pointsize 1 pointtype variable

set title 'Predicted Saturation vs horizontal angle, colored by va'
set output 'samples_s3pxp.png'
splot "samples.dat" u (column("px")):(column("sat_hat")):(column("py")):(column("gpsymbol")) with points palette pointsize 1 pointtype variable

set yrange [-0.3:0.3]
set title 'Error Saturation vs va, colored by px'
set output 'samples_s3pye.png'
splot "samples.dat" u (column("py")):(column("sat_err")):(column("px")):(column("gpsymbol")) with points palette pointsize 1 pointtype variable

set title 'Error Saturation vs horizontal angle, colored by va'
set output 'samples_s3pxe.png'
splot "samples.dat" u (column("px")):(column("sat_err")):(column("py")):(column("gpsymbol")) with points palette pointsize 1 pointtype variable



set yrange [0.0:1.0]
set title 'py vs px, colored by value'
set output 'samples_v3pypx.png'
splot "samples.dat" u (column("px")):(column("py")):(column("val")) with points palette pointsize 1 pointtype 7

set title 'Value vs va, colored by px'
set output 'samples_v3pya.png'
splot "samples.dat" u (column("py")):(column("val")):(column("px")):(column("gpsymbol")) with points palette pointsize 1 pointtype variable

set title 'Normalized Value vs horizontal angle, colored by va'
set output 'samples_v3pxn.png'
splot "normalized_samples.dat" u (column("px")):(column("nv")):(column("py")):(column("gpsymbol")) with points palette pointsize 1 pointtype variable

set title 'Value vs horizontal angle, colored by va'
set output 'samples_v3pxa.png'
splot "samples.dat" u (column("px")):(column("val")):(column("py")):(column("gpsymbol")) with points palette pointsize 1 pointtype variable

set title 'Predicted Value vs va, colored by px'
set output 'samples_v3pyp.png'
splot "samples.dat" u (column("py")):(column("val_hat")):(column("px")):(column("gpsymbol")) with points palette pointsize 1 pointtype variable

set title 'Predicted Value vs horizontal angle, colored by va'
set output 'samples_v3pxp.png'
splot "samples.dat" u (column("px")):(column("val_hat")):(column("py")):(column("gpsymbol")) with points palette pointsize 1 pointtype variable
