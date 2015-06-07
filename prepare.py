#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Create a reduction script for The Best and Brightest Metal-Poor Stars observed
through Gemini in 2015A.

Note: The reduction script uses Ureka from STScI, so you will need to install
      that before running the `reduce.sh` file produced by this script.
"""

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"
__version__ = "June, 2015"


import os
import logging
import time
from glob import glob

import numpy as np
from astropy.io import fits
from astropy.table import Table

# Set up logging.
logging.basicConfig(level=logging.DEBUG, filename="preparation.log",
    format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s",
    datefmt="%m-%d %H:%M", filemode="w")

console = logging.StreamHandler()
console.setLevel(logging.INFO)
formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
console.setFormatter(formatter)
logging.getLogger("preparation").addHandler(console)

logger = logging.getLogger("preparation")


# Specify the raw and reduced data folders.
REDUCTION_SCRIPT = "reduce_all.sh"
RAW_DATA_FOLDER = "raw_data/"
REDUCED_DATA_FOLDER = "reduced_data/"


def observation_summary(filenames, additional_headers=None):
    """
    Create a table to serve as an observing summary.
    """

    headers = ["INSTRUME", "OBJECT", "OBSTYPE", "OBSCLASS", "GEMPRGID", "OBSID",
    "DATALAB", "RA", "DEC", "UT", "DATE", "TIME-OBS", "GRATING", "EXPOSURE"]
    if additional_headers is not None:
        headers += additional_headers

    rows = []
    for filename in filenames:
        with fits.open(filename) as image:
            rows.append(dict(zip(headers,
                [image[0].header.get(h, None) for h in headers])))

        # Add the filename.
        rows[-1]["FILENAME"] = filename

        # Create a time stamp from the header information.
        time_string = "{0}T{1}".format(rows[-1]["DATE"], rows[-1]["UT"]) \
            if "T" not in rows[-1]["DATE"] else rows[-1]["DATE"]
        rows[-1]["TIME"] = time.mktime(time.strptime(time_string.split(".")[0],
            "%Y-%m-%dT%H:%M:%S"))
    
    headers.insert(0, "FILENAME")
    headers.append("TIME")

    return Table(rows=rows, names=headers)


# IT BEGINS...
logger.info("Raw data folder: {}".format(RAW_DATA_FOLDER))
logger.info("Reduced data folder: {}".format(REDUCED_DATA_FOLDER))
if not os.path.exists(REDUCED_DATA_FOLDER):
    os.mkdir(REDUCED_DATA_FOLDER)

# Gunzip any files in the raw image directory.
gzipped_raw_filenames = glob(os.path.join(RAW_DATA_FOLDER, "*.fits.gz"))
if len(gzipped_raw_filenames) > 0:
    logger.info("Gunzipping {} raw data images..".format(
        len(gzipped_raw_filenames)))

    for filename in gzipped_raw_filenames:
        os.system("gunzip -fv {}".format(filename))


logger.info("Creating summary table.")
summary = observation_summary(glob("{}/*.fits".format(RAW_DATA_FOLDER)))

science_objects = set(summary["OBJECT"]).difference(["Bias", "CuAr", "GCALflat"])
logger.info("{0} science objects".format(len(science_objects)))

commands = ["ur_setup"] # These are the commands for the final shell script.
for science_object in science_objects:
    logger.info("Creating structure for {}".format(science_object))

    # Create the folder for this object.
    folder = os.path.join(REDUCED_DATA_FOLDER, science_object)
    if not os.path.exists(folder):
        os.mkdir(folder)
        logger.info("Created folder {}".format(folder))

    indices = np.where(summary["OBJECT"] == science_object)[0]

    # Create symbolic links for science data.
    for filename in summary["FILENAME"][indices]:
        basename = os.path.basename(filename)
        source_filename = os.path.abspath(filename)
        destination_filename = os.path.join(folder, basename)

        if not os.path.exists(destination_filename):
            os.symlink(source_filename, destination_filename)

            logger.info("Created symbolic link {0} --> {1}".format(source_filename,
                destination_filename))

    # Use the first observation to associate an ARC and FLAT.
    obs = summary[indices[0]]
    instrument_match = (summary["INSTRUME"] == obs["INSTRUME"])
    flats = np.where((summary["OBSTYPE"] == "FLAT") * instrument_match)[0]
    arcs = np.where((summary["OBSTYPE"] == "ARC") * instrument_match)[0]

    # Get the closest arc and flat.
    flat_index = flats[np.argmin(np.abs(obs["TIME"] - summary["TIME"][flats]))]
    arc_index = arcs[np.argmin(np.abs(obs["TIME"] - summary["TIME"][arcs]))]

    logger.info("Found arc {0} for {1}".format(summary["FILENAME"][arc_index],
        science_object))
    logger.info("Found flat {0} for {1}".format(summary["FILENAME"][flat_index],
        science_object))

    # Create symbol links for arc and flat data.
    for filename \
    in (summary["FILENAME"][flat_index], summary["FILENAME"][arc_index]):

        basename = os.path.basename(filename)
        source_filename = os.path.abspath(filename)
        destination_filename = os.path.join(folder, basename)

        if not os.path.exists(destination_filename):
            os.symlink(source_filename, destination_filename)

            logger.info("Created symbolic link {0} --> {1}".format(
                source_filename, destination_filename))

    # Create symbolic link for the reduction script.
    source_filename = os.path.abspath("reduce.py")
    destination_filename = os.path.join(folder, os.path.basename("reduce.py"))
    if not os.path.exists(destination_filename):
        logger.info("Created symbolic link {0} --> {1}".format(source_filename,
            destination_filename))
        os.symlink(source_filename, destination_filename)

    # Add the commands that we will need to perform the reductions.
    commands.append("cd {}".format(os.path.abspath(folder)))
    commands.append("mkiraf -f")
    commands.append("python {}".format(os.path.basename("reduce.py")))

# Change directories back to where we started, and forget Ureka.
commands.append("cd {}".format(os.getcwd()))
commands.append("ur_forget")

# Create the reduction script, and make it executable.
with open(REDUCTION_SCRIPT, "w+") as fp:
    fp.write("\n".join(commands))

logger.info("Created reduction script {}".format(REDUCTION_SCRIPT))
logger.info("Now you should run `source {}`".format(REDUCTION_SCRIPT))

