# some examples
The file doall.sh shows the command line parameters used

Generally the filename prefixes that end in "sf" are the skyfilled images

These are jpegs created from the original tiffs -- the skyfill_tif program will
not work with these images.  The original tiffs are too large to put on github, so
these examples were created to show the output from skyfill_tif.

---

For the extremely curious:
* When you run skyfill, a file named "samples.dat" will also be created.  It is the samples
that skyfill found in the sky, and used to model the HSV across the sky in the X and Y dimensions

* "samples.gnuplot" is used to create some visual scatterplots to help analyze
how well the sky HSV model worked.  Simply run "gnuplot < samples.gnuplot" to create the scatterplots
