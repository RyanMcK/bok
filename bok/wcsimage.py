#!/usr/bin/env python

import os.path
from subprocess import call

import numpy as np

import astropy.io.fits as fits

class WCSImage(object):
    def __init__(self, imgpath):
        self.imgpath = imgpath
        basedir, filename = os.path.split(imgpath)
        
        weightname = filename.replace(".pw.", ".p.")
        weightname = weightname.replace(".fits", ".weight.fits")
        self.weight = os.path.join(basedir, weightname)

        catname = filename.replace(".fits", ".cat")
        self.cat = os.path.join(basedir, catname)

    def run_sextractor(self):
        params = {'CATALOG_TYPE': 'FITS_1.0',
                  'DETECT_TYPE': 'CCD',
                  'DETECT_MINAREA': '5',
                  'DETECT_THRESH': '1.0',
                  'ANALYSIS_THRESH': '5.0',
                  'FILTER': 'Y',
                  'SATUR_LEVEL': '50000',
                  'MAG_ZEROPOINT': '27.5',
                  'BACK_SIZE': '128',
                  'BACK_FILTERSIZE': '3',
                  'BACKPHOTO_TYPE': 'GLOBAL',
                  'WEIGHT_TYPE': 'MAP_VAR',
                  'FILTER_NAME': '/project/projectdirs/deepsky/ptf/soft/cvs/ptf/wcs/default.conv',
                  'STARNNW_NAME': '/project/projectdirs/deepsky/ptf/soft/cvs/ptf/wcs/default.nnw',
                  'WEIGHT_IMAGE': self.weight,
                  'CATALOG_NAME': self.cat,
                  'PARAMETERS_NAME': '/global/project/projectdirs/cosmo/work/bok/proc/photometry/photometry.param'}

        params = " ".join("-%s %s" % (k, v) for k, v in params.iteritems())
        cmd = ["sex", self.imgpath] + params.split()
        call(cmd)
