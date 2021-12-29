# SkyFill
Fill the top of clear skies in stitched panoramas

Currently, it is designed to work with 4 channel, 16 bit TIFF files from Hugin.

Having the corresponding .pto for the project is very helpful as skyfill
is able to parse the .pto file for relavent information (FOV) to help predict
the sky image needed to blend with the stitched image coming from Hugin

*************************************************************************
THIS PROJECT IS STILL UNDER HEAVY DEVELOPMENT, EXPECT LOTS OF WARNING MESSAGES FROM THE COMPILER
AND CHANGES TO THE FUNCTIONALITY AND RESULTS.  HOWEVER, IT MAY STILL BE USEFUL
IN MANY CASES
*************************************************************************


Dependencies:

* libtiff
* exiftool

TODO (Dec 29, 2021):
* Better implementation of blending CIE sun model into the new hsv sky model
* fix sky hue (-fsh) leaves noticeable edge in repaired area
