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

USING:
* Create a 16 bit stitched panorama with hugin, it must be 4 channel.  Hugin sets alpha to 0 for areas with
  no image data.  You can aquire your image via another method but this is specifically designed for Hugin.  If you
  don't have Hugin and know the Horizontal FOV, specify it on the commandline with (-FOV \<degrees\>)

* Run it with skyfill_tif \<inputimage\> \<outputimage\>

* Skyfill detects the start of sky, and end of sky for each column (X value) in the image.  *If it does not get the
  end of sky correct (usually by marking the end of sky below the actual end of sky), the subsequent estimate of HSV will be 
  corrupted by using non-sky pixels in the estimate of actual sky HSV.*   Clouds are also a potential problem.

  - Run skyfill with the -d2 flag to see where it thinks the start and end of sky is at.  It will put a green line two pixels
    *below* where it has identified the end of sky, so you can visually see if the actual pixel color at the detected
    end of sky looks correct.
  - If the end of sky looks incorrect, the only option to fix it is by changing the tolerances for end of sky detection
     "-tol \<th\> \<ts\> \<tv\>" -- tolerance for determining edge of skyline at horizon, default 360 .06 0.015

     th is the allowed change in Hue in the actual sky pixels, any change in hue from one pixel to another greater than this value
     will trigger end of sky to be detected.
     
     ts and tv work the same but for Saturation and Value
 

DEPENDENCIES:

* libtiff
* exiftool

WARNING:
* it will compile cleanly using the enclose Makefile, BUT:
  - there are no dependencies for the required header files for each source module
  - The CMakeLists.txt file is out of date and needs to be updated.
* Still a lot of cleanup to do, but the breaking apart into separate source modules was a huge first step

UPDATES:
Jan 1, 2022
* Fixed bug in sky HSV model process
* Automatic detection and attempted repair when top of sky is too dark
* Made the output less verbose by default
Dec 31, 2021
* The new full sky replacement mode (-fsr switch on the commandline) is functional.
* The sky HSV model has been improved, made more robust, and will likely change again.

TODO:
* Probably isn't going to work for skies at dusk/dawn, need to test
* Better blending at edges of masked columns
* Run a smoother along start of sky in from original image to smooth the seam lines
