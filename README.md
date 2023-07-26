# SkyFill
Fill the top of clear skies in stitched panoramas

Currently, it is designed to work with 4 channel, 8 or 16 bit TIFF files from Hugin.
Because of artifacts created in jpeg compression, there is no attempt to provide I/O for jpeg files

Having the corresponding .pto for the project is *extremely* helpful as skyfill
is able to parse the .pto file for relavent information (FOV) to help predict
the sky image needed to blend with the stitched image coming from Hugin

*************************************************************************
As of Feb 12, 2022 This is alpha level code.  Major development has ended, and documentation has begun
Only bugs found will be addressed.  Bugs do not include very hard cases of images where the sky is
too complex for SkyFill to get a good estimate of where the sky is, and how to fill various areas of
the sky
*************************************************************************

**CREDITS**
***Levenberg Marquardt Fitting code***

*   Markwardt, C. B. 2009, "Non-Linear Least Squares Fitting in IDL with MPFIT," in proc. Astronomical Data Analysis Software and Systems XVIII, Quebec, Canada, ASP Conference Series, Vol. 411, eds. D. Bohlender, P. Dowler & D. Durand (Astronomical Society of the Pacific: San Francisco), p. 251-254 (ISBN: 978-1-58381-702-5; Link to ASP title listing)
    ADS Bibcode: 2009ASPC..411..251M (click for Bibtex and other citation formats)
    Arxiv preprint: arXiv:0902.2850v1
    Bibtex entry:

*   \bibitem[Markwardt(2009)]{2009ASPC..411..251M} Markwardt, C.~B.\ 2009,
    Astronomical Data Analysis Software and Systems XVIII, 411, 251

*   http://purl.com/net/mpfit

*   Moré, J. 1978, "The Levenberg-Marquardt Algorithm: Implementation and Theory," in Numerical Analysis, vol. 630, ed. G. A. Watson (Springer-Verlag: Berlin), p. 105

*The code I found was at this location in February of 2022*
    http://cow.physics.wisc.edu/~craigm/idl/cmpfit.html

***Thanks to***

* Thomas Modes for providing CMake functionality, and finding compilation issues on Windows.

**Tutorials**

There are three [tutorials](Tutorial/tutorial_index.md) describing the usage of SkyFill, and some of the underlying concepts

**Getting Started**

* Create a 8 or 16 bit stitched panorama with hugin, it must be 4 channel.  Hugin sets alpha to 0 for areas with
  no image data.  You can aquire your image via another method but this is specifically designed for Hugin.  If you
  don't have Hugin and know the Horizontal FOV, specify it on the commandline with (-FOV \<degrees\>)

* Run it with skyfill\_tif \<inputimage\> \<outputimage\>

* Skyfill detects the start of sky, and end of sky for each column (X value) in the image.  *If it does not get the end of sky correct (usually by marking the end of sky below the actual end of sky), the subsequent estimate of HSV will be 
  corrupted by using non-sky pixels in the estimate of actual sky HSV.*   Vegetation, clouds are examples of potential problems.
  Skyfill analyzes the sky in columns (i.e. the "y") dimension.  It starts at the first opaque pixels in an image column and moves downward until 2 or more adjacent pixels in the column do not match the expected values for the sky in that localize
  area.

  - Run skyfill with the -d3 flag to see where it identified the start and end of sky is at.  It will put a green line two pixels *below* where it has identified the end of sky, so you can visually see if the actual pixel color at the detected end of sky looks correct.
  - If the end of sky looks incorrect, the only option to fix it is by changing the tolerances for end of sky detection
     "-tol \<t\>" -- tolerance for determining edge of skyline at horizon, default is 3.5  Increasing the tolerance will cause
     the algorithm to move farther down the sky in a column.  Decreasing the tolerance will have the opposit effect.

**Masks**

* Masks are an important option for handling images that are not simple situations with unobstructed clear sky at the top of the image
* Column Mask
  "-m \<l\> \<r\>"  -- this flag prevents any analysis, or changes to image columns from *l* to *r*
* Sample Mask
  "-sm \<l\> \<r\>"  -- this flag prevents using image columns from *l* to *r* as estimates of clear sky for purposes of modelling the sky HSV values
* Repair Mask
  "-rm \<l\> \<r\>"  -- this flag prevents image columns from *l* to *r* from having changes make to pixels in the sky, though the area may still be used for sampling actual sky values
* Non Sky Repair Mask
  "-nsrm \<l\> \<r\>"  -- this flag prevents image columns from *l* to *r* from having changes make to pixels in the sky which are classsified as likely to be not a sky pixel (i.e. a cloud, or a piece of vegetation)
* Test Mask
  "-tm \<left\> \<right\> \<top\> \<bottom\>"  -- this flag identifies a rectangle that is skipped over for end of sky detection, and no sky samples are taken in the rectangle.  In the examples, this has been used to exclude regions with direct sun and lens flare that are in the actual sky area, but would be bad to use for detecting end of sky and for modelling the sky color

**Normal sky replacement**

* After a sky HSV model has been created by SkyFill, it then processes the image one column at a time. It blends the modelled sky (now starting at y=0), and blending to a y value that is 50% of the way from the first opaque pixel (top of sky), to the identified end of sky.  

* The "-df \<D\>" flag determines how far down to stop the blending.  The default is 50%.    "-df 0", would leave the original sky pixels completely unchanged and only fill in the pixels above -- (leaving a noticeable seam due to errors from the sky model).

* "-df 1." will fill sky pixels all the way down to the detected end of sky.  The blending uses 100% of the modeled sky at the original top of the sky in the image, and linearly changes to using 100% of the original pixel value at the lowest (largest y) pixel changed.

* If the end of sky is highly variable, there may be some visible artifacts in the resulting sky from this method and it would be worth trying Full Sky Replacement Mode

**Full Sky Replacement Mode**

* When you run it in full sky replacement mode (-fsr).  The placement of the horizon determines the lowest part of the image where any pixel replacement will occur. To see how far down in the image replacement is going to occur, run it with an undocumented   flag (-SRE).  The output will be grayscale down to the horizon (which may be curved).  White pixels will be fully replaced, black pixels will not be replaced.  You can move the horizon line with the -hy \<p\> flag, where \<p\>=1.0 is the top of the image, \<p\>=0.5 is halfway between the top and bottom...

* In this mode, you can supply parameters for putting a modeled sun in the sky.  The actual sun diameter is much smaller in an image than you might first guess, because so much of the area around the sun has been blown out and exceeded the sensor maximum  value.   The parameters for defining the sun are:
  * -sx \<Xangle\> -- in degrees -- Xangle of 0 is the middle of the image, Xangle=-FOV/2 is left side of image, Xangle=FOV/2 is right side of image
  * -sy \<Py\> -- as a proportion, 1 is top of image, 0 is bottom of image
  * --- note, you can put the sun outside of the image, and get just some glow from the modelled sun to appear in the image
  * -C \<Cvalue\> -- overall brighness of the sun.  Default is 2
  * -D \<Dvalue\> -- overall diameter of the sun.  Default is .01  You *must* set this to a large value, try the range 1 to 10.  The "size" of the generated sun will be somewhat linear with respect to the Dvalue chosen
  * -F \<Fvalue\> -- Blending amount of the modeled sun into the normal modelled sky, default is 0.5
  * -G \<Gvalue\> -- This controls how much of the sun+modelled sky to use in replacing pixel values.  Think of this as a "blending factor" that overrides the normal top of sky to end of sky blending factor.  Gvalue is given as a proportional distance (where the proportions are 0 to 1 for both image width and image height)  Default is 0, so you *must* set this to a nonzero value to have an effect.  Try small values, .1 is a good first try

**Miscellaneous Flags**

  * -nceos -- "no clip end of sky".  In full sky replacement mode, will disable feathering for any pixels up to end of sky.  If the sky model is different from the actual sky near the horizon, this may produce a more consistent result
  * -fsrt \<thresh\> \<ramp\> -- Sets the threshold and ramp for probability a pixel is identified as a sky pixel (and will be replaced).  Default thresh is 0.9, default ramp is 200.  Reasonable thresh values are 0.5 to 0.99.  Higher values result in less pixels identified as sky.  ramp changes how fast probability changes around thresh.  ramp of 10 gives a very slow change, while ramp of 200 gives a reasonably fast change.
  *  -fsh \<pl\> \<pr\> -- Fix Sky Hue.  If a small area of the sky has been blown out, the hue will be incorrect.  If there is also an area of sky adjacent to the blown out area, a repair can be attempted.  pl, and pr are the left and right markers for the area of the sky to attempt a repair.  They are expressed as a proportion, i.e "-fsh 0.75 1.0" will attempt to repair the right 25% of the image where sky is detected.
  * -msu \<0|1|2\>  -- run a smoothing filter over the sky after all processing is complete.  Recommended mode is "-msu 0", the other modes are experimental.
  * -sd \<S\> -- a factor applied to the esimated sky saturation model, default is 1.0
  * -r_tos_thresh \<thresh\> -- threshold to ignore pixels in very localized areas of the sky for estimating local sky rgb in preliminary steps repairing errant sky pixels.  Default is 0.02.  Reasonable values are 0.01 to 0.05.   The default will be sufficient in most cases.
  * -LF <filename> -- add lens flare and aperture ghosts created by specular light (i.e. the sun) defined in *filename*

**Adding lensflare**

Using a simple model, lens flare (rays, aperture ghosts) can be rendered onto the image
There is [documentation](Lensflare.md) for a lens flare data file


**DEPENDENCIES**

* libtiff
* exiftool

**INSTALLATION**

* it will compile cleanly using the enclosed Makefile on Linux
* Thomas Modes has supplied a CMakeLists.txt file which allows use of the CMake system to compile on Linux or Windows or ???

**UPDATES**

Jan 1, 2022
* Fixed bug in sky HSV model process
* Automatic detection and attempted repair when top of sky is too dark
* Made the output less verbose by default

Dec 31, 2021
* The new full sky replacement mode (-fsr switch on the commandline) is functional.
* The sky HSV model has been improved, made more robust, and will likely change again.

**TODO**

* Probably isn't going to work for skies at dusk/dawn, need to test
