#!/usr/bin/env python

import os

import numpy as np

from astropy.io import fits

from rawimage import RawImage
import util

class StackedFlat(object):
    def __init__(self, filenames, stackedbias, stackedbias_err,
                 stackedbias_mask, hdu_num):
        self.filenames = filenames
        self.nfiles = len(filenames)
        self.stackedbias = fits.open(stackedbias)
        self.stackedbias_err = fits.open(stackedbias_err)
        self.stackedbias_mask = fits.open(stackedbias_mask)
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
            cube[j, :, :] = self.images[j].data - self.stackedbias[1].data.T
            cube_err = np.sqrt(self.images[j].data_err**2 +
                               self.stackedbias_err[1].data.T**2)
            
            med = np.median(cube[j, :, :])
            cube[j, :, :] /= med
            # If we're normalizing the counts, we should do the same for
            # the errors?
            cube_err /= med
            
            cube_weights[j, :, :] = cube_err**-2.0

        cubemedian = np.median(cube, axis=0)
        # We took the median, not the weighted mean, for the counts, but
        # we'll approximate the error using the weights.
        # See http://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Dealing_with_variance.
        sum_weights = np.sum(cube_weights, axis=0)
        cubestd = np.sqrt(sum_weights**-1.0)

        med = np.median(cubemedian)
        cubemedian /= med
        # Again, normalize the errors whenever we normalize the counts.
        cubestd /= med
        self.data = cubemedian
        self.data_err = cubestd

        # Use the superbias and superflat errors to create a mask, where
        # 0 = good pixel, 1 = bad bias pixel, and 2 = bad flat pixel.
        pixmask = self.stackedbias_mask[1].data.T
        flat_err = cubestd
        cflat_err = flat_err - np.median(flat_err)
        ns = np.percentile(cflat_err, 15.87)
        ps = np.percentile(cflat_err, 84.13)
        sigma_err = max(np.abs(ns), np.abs(ps)) 
        bad_pix = ((flat_err == 0.0) | \
                   (np.abs(flat_err - np.median(flat_err)) > 5.0*sigma_err)) & \
                  (pixmask == 0)
        pixmask[bad_pix] = 2

        self.data_mask = pixmask

    def save(self):
        header = fits.Header()
        header['OBJECT'] = "stackedflat"
        header['COMMENT'] = "Stackedflat FITS file."
        priHDU = fits.PrimaryHDU(header=header)

        header = fits.Header()
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
        newpath = os.path.join(basedir, "stackedflat.%02d.fits" % self.hdu_num)
        hdulist.writeto(newpath, clobber=True)
 
        hdulist_err = fits.HDUList([imgHDU_err])
        newpath_err = os.path.join(basedir, "stackedflat.err.%02d.fits" % self.hdu_num)
        hdulist_err.writeto(newpath_err, clobber=True)

        hdulist_mask = fits.HDUList([imgHDU_mask])
        newpath_mask = os.path.join(basedir, "stackedflat.mask.%02d.fits" % self.hdu_num)
        hdulist_mask.writeto(newpath_mask, clobber=True)
