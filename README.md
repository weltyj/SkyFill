# SkyFill
Fill the top of clear skies in stitched panoramas

Currently, it is designed to work with 4 channel, 16 bit TIFF files from Hugin.

Having the corresponding .pto for the project is *extremely* helpful as skyfill
is able to parse the .pto file for relavent information (FOV) to help predict
the sky image needed to blend with the stitched image coming from Hugin

*************************************************************************
THIS PROJECT IS STILL UNDER HEAVY DEVELOPMENT, EXPECT LOTS OF WARNING MESSAGES FROM THE COMPILER
AND CHANGES TO THE FUNCTIONALITY AND RESULTS.  HOWEVER, IT MAY STILL BE USEFUL
IN MANY CASES
*************************************************************************

Dec 31, 2021
NEWS:
* The new full sky replacement mode (-fsr switch on the commandline) is functional.
* The sky HSV model has been improved, made more robust, and will likely change again.

CODE:
WARNING -- it will compile cleanly using the enclose Makefile, BUT:
* there are no dependencies for the required header files for each source module
* The CMakeLists.txt file is out of date and needs to be updated.
* Still a lot of cleanup to do, but the breaking apart into separate source modules was a huge first step

USING:
* Create a 16 bit stitched panorama with hugin, it must be 4 channel.  Hugin sets alpha to 0 for areas with
  no image data.  You can aquire your image via another method but this is specifically designed for Hugin.  If you
  don't have Hugin and know the Horizontal FOV, specify it on the commandline with (-FOV <degrees>)

* Run it with skyfill_tif *inputimage* *outputimage*

Dependencies:

* libtiff
* exiftool
