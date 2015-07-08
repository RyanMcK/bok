#!/usr/bin/env python

import os

import numpy as np

from astropy.io import fits
import astropy.stats as stats

from rawimage import RawImage
import util

class StackedBias(object):
    def __init__(self, filenames, hdu_num):
        self.filenames = filenames
        self.nfiles = len(filenames)
        self.hdu_num = hdu_num
        self.images = []
        for filename in filenames:
            img = RawImage(filename, hdu_num)
            img.subtract_overscan()
            self.images.append(img)

        nx, ny = self.images[0].data.shape
        cube = np.zeros((self.nfiles, nx, ny))
        cube_weights = np.zeros_like(cube)
        for j in range(self.nfiles):
            cube[j, :, :] = self.images[j].data
            cube_weights[j, :, :] = self.images[j].data_err**-2.0

        cubeclipped = stats.sigma_clip(cube, sig=3.0, cenfunc=np.mean, axis=0)

        # Compute the weighted mean for each pixel in the stack, using
        # only those weights for counts that weren't clipped.
        cubemean = np.ma.average(cubeclipped, axis=0, weights=cube_weights)
        # In the very unlikely (impossible?) event that all pixels in the 
        # stack were masked, assume zero counts.
        cubemean = np.ma.filled(cubemean, 0.0)
        self.data = cubemean

        # Use pairwise differences of overscan-subtracted biases to estimate
        # the error.
        rmses = []
        for i in range(0, self.nfiles-1, 2):
            diff = self.images[i+1].data[self.images[i+1].py_datasec] - \
                   self.images[i].data[self.images[i].py_datasec]
            mn = np.mean(diff)
            rms = np.sqrt(np.mean((diff - mn)**2.0))
            rmses.append(rms)
        mean_rms = np.mean(rmses)
        self.data_err = np.full_like(self.data, mean_rms)
        print "hdu %d, bias diff rms %f, scatter %f" % (hdu_num, mean_rms, np.sqrt(np.mean((np.array(rmses)-mean_rms)**2.0)))
        
        # Mask those pixels whose rms in the stack of unclipped pixels is
        # greater than five times the rms measured from pairwise differences.
        rms_arr = np.ma.average((cubeclipped-cubemean)**2.0, axis=0, weights=cube_weights)
        rms_arr = np.sqrt(np.ma.filled(rms_arr,  0.0))

        pixmask = np.zeros_like(cubemean, dtype=int)
        bad_pix = ((rms_arr == 0.0) | \
                   (rms_arr > 5.0*mean_rms)) & \
                  (pixmask == 0)
        pixmask[bad_pix] = 1
        self.data_mask = pixmask

    def save(self):
        header = fits.Header()
        header['OBJECT'] = "stackedbias"
        header['COMMENT'] = "Stackedbias FITS file."
        header['DETSEC'] = self.images[0].detsec
        header['CCDSEC'] = self.images[0].ccdsec
        header['DATASEC'] = self.images[0].datasec
        header['BIASSEC'] = self.images[0].biassec
        header['CCDNAME'] = self.images[0].ccdname
        header['AMP-CFG'] = self.images[0].ampcfg
        
        # Transpose because of the way FITS handles indexing.
        imgHDU = fits.PrimaryHDU(data=self.data.T, header=header)
        imgHDU_err = fits.PrimaryHDU(data=self.data_err.T, header=header)
        imgHDU_mask = fits.PrimaryHDU(data=self.data_mask.T, header=header)

        hdulist = fits.HDUList([imgHDU])
        basedir, filename = os.path.split(self.images[0].imgpath)
        basename, ext = os.path.splitext(filename)
        newpath = os.path.join(basedir, "stackedbias.%02d.fits" % self.hdu_num)
        hdulist.writeto(newpath, clobber=True)

        hdulist_err = fits.HDUList([imgHDU_err])
        newpath_err = os.path.join(basedir, "stackedbias.err.%02d.fits" % self.hdu_num)
        hdulist_err.writeto(newpath_err, clobber=True)

        hdulist_mask = fits.HDUList([imgHDU_mask])
        newpath_mask = os.path.join(basedir, "stackedbias.mask.%02d.fits" % self.hdu_num)
        hdulist_mask.writeto(newpath_mask, clobber=True)
