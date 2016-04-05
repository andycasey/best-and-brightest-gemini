#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
An `in-folder` reduction script for The Best and Brightest Metal-Poor Stars
observed through Gemini.

Note: You should really be running prepare.py, which will generate a `reduce.sh`
      script that calls this file.
"""

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"
__version__ = "June, 2015"


import logging
import os
import time
from glob import glob

import numpy as np
from astropy.io import fits
from astropy.table import Table

# I am ashamed:
try:
    from pyraf import iraf

except ImportError:
    raise ImportError("No pyraf module -- did you forget to do run `ur_setup`?")

# Set up logging.
logging.basicConfig(level=logging.DEBUG, filename="reduction.log",
    format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s",
    datefmt="%m-%d %H:%M", filemode="w")

console = logging.StreamHandler()
console.setLevel(logging.INFO)
formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
console.setFormatter(formatter)
logging.getLogger("reduction").addHandler(console)

logger = logging.getLogger("reduction")

def log_iraf_result(result, raise_errors=False):
    message = "\n".join(result)
    logger.info("Result:\n{0}".format(message))
    if raise_errors and "error" in message:
        raise IrafError(message)

if not os.path.exists("login.cl"):
    raise IOError("no login.cl file -- have you run mkiraf? [select xgterm]")

logger.info("Learning..")

iraf.fitsutil()
iraf.gemini()
iraf.gmos()

logger.info("Unlearning..")

iraf.unlearn("gemini")
iraf.unlearn("gmos")
iraf.unlearn("gemtools")
iraf.set(stdimage="imtgmos2", homedir=os.path.expanduser("~/"))

logger.info("Cleaning up from any previous reductions..")
os.system("rm -Rf master_flat.fits mosaic_master_flat.fits stgscg*.fits "
    "cg*.fits tgs*.fits g*.fits estgsc*.fits database J??????.?????????.?.fits")
                                                      
files = np.array(glob("*.fits"))

logger.info("Getting object types..")
obstypes = np.array([iraf.hselect(filename + "[0]", fields="OBSTYPE",
    expr="1=1", Stdout=1)[0] for filename in files])

obsclass = np.array([iraf.hselect(filename + "[0]", fields="OBSCLASS",
    expr="1=1", Stdout=1)[0] for filename in files])

folder = os.path.dirname(files[0])

# Check that we have one of each
if "OBJECT" not in obstypes:
    raise IOError("no OBJECT image found in {0}".format(folder))

if "ARC" not in obstypes:
    raise IOError("no ARC image found in {0}".format(folder))

if "FLAT" not in obstypes:
    raise IOError("no FLAT image found in {0}".format(folder))

arc_filenames = files[obstypes == "ARC"]
flat_filename = files[obstypes == "FLAT"][0]
object_filenames = files[(obstypes == "OBJECT") * (obsclass == "science")]

# Prepare the data
logger.info("Preparing..")
prepare_result = iraf.gprepare("*.fits", fl_addmdf=True, Stdout=1)
log_iraf_result(prepare_result)


# Do cosmic ray rejection on the science frames
for filename in object_filenames:
    logger.info("Doing cosmic ray rejection on {}..".format(filename))
    cosmic_ray_result = iraf.gscrrej("g{0}".format(filename),
        "cg{0}".format(filename), Stdout=1)
    log_iraf_result(cosmic_ray_result)


# Create a master flat
logger.info("Creating master flat..")
flat_result = iraf.gsflat(inflats="g{0}".format(flat_filename),
    specflat="master_flat.fits", fl_over=False, fl_trim=True, fl_dark=False,
    fl_fixpix=False, fl_inter=False, function="chebyshev", order=15,
    fl_detec=True, ovs_flinter=False, fl_vardq=False, fl_bias=False, Stdout=1)
log_iraf_result(flat_result)


# Create a mosaic of the master flat
logger.info("Creating master flat mosaic..")
mosaic_result = iraf.gmosaic("master_flat.fits", outpref="mosaic_", Stdout=1)
log_iraf_result(mosaic_result)

# Reduce the science frame
logger.info("Reducing science frame(s)..")
for object_filename in object_filenames:
    logger.info("Reducing {}".format(object_filename))
    reduce_science_result = iraf.gsreduce("cg{0}".format(object_filename),
        fl_inter=False, fl_over=False, fl_trim=True, fl_dark=False,
        fl_flat=True, flatim="mosaic_master_flat.fits", fl_gmosaic=True,
        fl_fixpix=True, fl_bias=False, fl_cut=True, fl_gsappwave=True,
        ovs_flinter=False, fl_vardq=False, yoffset=5.0, Stdout=1)
log_iraf_result(reduce_science_result)


# Reduce the arc frames
logger.info("Reducing arc frames..")
for filename in arc_filenames:
    reduce_arc_result = iraf.gsreduce("g{0}".format(filename), 
        fl_over=False, fl_trim=True, fl_bias=False, fl_dark=False, 
        fl_flat=False, fl_cut=True, fl_gsappwave=True, yoffset=5.0, Stdout=1)
    log_iraf_result(reduce_arc_result)

    logger.info("Running gswavelength..")
    #gswavelength_result = iraf.gswavelength("gsg{0}".format(filename),
    #    fl_inter="NO", nsum=5, step=5, function='chebyshev', order=6, 
    #    fitcxord=5, fitcyord=4, Stdout=1)
    gswavelength_result = iraf.gswavelength("gsg{0}".format(filename),
         fl_inter="NO", nsum=5, step=5, function="chebyshev", order="6",
         fitcxord=5, fitcyord=4, Stdout=1, niterate=10, low_reject=1.5,
         high_reject=1.5, fl_dbwrite="YES", fl_overwrite="yes")
    log_iraf_result(gswavelength_result)

    # Apply transformations
    logger.info("Applying wavelength transformations to arc..")
    gstransform_arc_result = iraf.gstransform("gsg{0}".format(filename),
        wavtraname="gsg{0}".format(filename), Stdout=1)
    log_iraf_result(gstransform_arc_result)


logger.info("Doing wavelength transformations, sky and extraction:")
for object_filename in object_filenames:
    logger.info("Working on filename {}".format(object_filename))
    logger.info("Applying wavelength transformations to object..")

    # Get nearest arc by time.
    arc_times = []
    for arc_filename in arc_filenames:
        image = fits.open(arc_filename)

        time_string = "{0}T{1}".format(
            image[0].header["DATE"], image[0].header["UT"]) \
            if "T" not in image[0].header["DATE"] else image[0].header["DATE"]
        arc_times.append(time.mktime(time.strptime(time_string.split(".")[0],
            "%Y-%m-%dT%H:%M:%S")))

    image = fits.open(object_filename)
    time_string = "{0}T{1}".format(
        image[0].header["DATE"], image[0].header["UT"]) \
        if "T" not in image[0].header["DATE"] else image[0].header["DATE"]
    obs_time = time.mktime(time.strptime(time_string.split(".")[0],
        "%Y-%m-%dT%H:%M:%S"))

    # Get closest arc.
    index = np.argmin(np.abs(np.array(arc_times) - obs_time))
    arc_filename = arc_filenames[index]

    log_iraf_result(["ASSOCIATING ARC {0} WITH FILENAME {1}".format(
        arc_filename, object_filename)])

    log_iraf_result(["ASSOCIATING FLAT {0} WITH FILENAME {1}".format(
        flat_filename, object_filename)])

    gstransform_object_result = iraf.gstransform(
        "gscg{0}".format(object_filename),
        wavtraname="gsg{0}".format(arc_filename), Stdout=1)
    log_iraf_result(gstransform_object_result)

    logger.info("Subtracting sky..")
    sky_subtraction_result = iraf.gsskysub("tgscg{0}".format(object_filename),
        fl_inter=False, Stdout=1)
    log_iraf_result(sky_subtraction_result)

    logger.info("Extracting..")
    extract_result = iraf.gsextract("stgscg{0}".format(object_filename),
        fl_inter=False, find=True, back="fit", bfunc="chebyshev", border=1,
        tfunct="spline3", torder=5, tnsum=20, tstep=50, refimage="",
        apwidth=1.3, recent=True, trace=True, fl_vardq=False, Stdout=1)
    log_iraf_result(extract_result)

logger.info("Producing 1D spectrum..")
reduced_filenames = glob("estgscg*.fits")

image = fits.open(reduced_filenames[0])
dispersion = image[2].header["CRVAL1"] \
    + np.arange(image[2].header["NAXIS1"]) * image[2].header["CD1_1"]

# Stack the flux (this assumes sequential exposures)
flux = np.zeros_like(dispersion)
for filename in reduced_filenames:
    with fits.open(filename) as science_image:
        flux += science_image[2].data.flatten()

    # And create an easy save file:
    primary_hdu = fits.PrimaryHDU(header=image[0].header)
    disp = image[2].header["CRVAL1"] \
        + np.arange(image[2].header["NAXIS1"]) * image[2].header["CD1_1"]
    data_hdu = fits.new_table([
        fits.Column(name="disp", format="1D", array=disp),
        fits.Column(name="flux", format="1D", array=flux),
        fits.Column(name="inv_var", format="1D", array=1.0/flux)])

    hdu_list = fits.HDUList([primary_hdu, data_hdu])
    hdu_list.writeto("{0}-{1}".format(image[0].header["OBJECT"],
        filename))

flux[0 >= flux] = np.nan

primary_hdu = fits.PrimaryHDU(header=image[0].header)
data_hdu = fits.new_table([
    fits.Column(name="disp", format="1D", array=dispersion),
    fits.Column(name="flux", format="1D", array=flux),
    # I am further ashamed:
    fits.Column(name="variance", format="1D", array=flux)
    ])
hdu_list = fits.HDUList([primary_hdu, data_hdu])
hdu_list.writeto("{}.fits".format(image[0].header["OBJECT"]))
logger.info("Created extracted spectrum {}.fits".format(image[0].header["OBJECT"]))

