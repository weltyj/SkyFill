SUNFLARE="-C 6.5 -D 0.5 -F 1.0 -G 0.12 -sx -2.45 -sy .87"
#../skyfill_tif ../../../tests/Tutorial/pan02_v2.tif lensflare_example_plain.tif $SUNFLARE -LF lf3_plain.dat
#convert lensflare_example_plain.tif lensflare_example_plain.jpg
rm lensflare_example_plain.tif
cp ../../../tests/Tutorial/lf3.jpg lensflare_example.jpg
exit

../lensflare ghost_samples1.tif -LF ghostsamples.dat
../lensflare flare_samples1.tif -LF flaresamples.dat
convert ghost_samples1.tif ghost_samples1.jpg
convert flare_samples1.tif flare_samples1.jpg
rm ghost_samples1.tif
rm flare_samples1.tif
