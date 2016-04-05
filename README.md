Reduce Gemini/GMOS long-slit spectra for The Best & Brightest Metal-Poor Stars 
==============================================================================

These scripts reduce Gemini/GMOS longslit spectra. The scripts are Python
but the guts are actually IRAF (ew, I know).

1) Clone this repository:

    git clone git@github.com:andycasey/gemini-gmos-reduction-script.git

2) Install [Ureka](http://ssb.stsci.edu/ureka/1.5.1/) from STScI.

3) Download all your Gemini science and calibration images into the `raw_data`
   directory.

4) Run the preparation script:

    python prepare.py
   
5) The preparation script will produce a `reduce_all.sh` script. Run it:

    source reduce_all.sh


The one-dimensional, wavelength-calibrated, extracted, stacked spectra will 
live in `reduced_data/<OBJECT_NAME>/<OBJECT_NAME>.fits`
