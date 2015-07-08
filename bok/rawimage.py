#!/usr/bin/env python

import os

import numpy as np

from astropy.io import fits
import astropy.stats as stats

import util

def parse_fits_range(fits_range, to_python=False):
    """Convert, e.g., ' [1:2016,1:2048]' to ((1, 2016), (1, 2048)) if to_python
    is False, and a tuple of slices if true.  From Yu Feng's frontloader."""
    x, y = fits_range.strip()[1:-1].split(',')
    x = tuple([int(i) for i in x.split(':')])
    y = tuple([int(i) for i in y.split(':')])
    if to_python:
        stepy = 1 if y[1] > y[0] else -1
        stepx = 1 if x[1] > x[0] else -1
        return (slice(x[0] - 1, x[1] + stepx - 1, stepx),
                slice(y[0] - 1, y[1] + stepy - 1, stepy))
    return x, y

class RawImage(object):
    def __init__(self, imgpath, hdu_num):
        self.imgpath = imgpath
        self.hdu_num = hdu_num

    def subtract_overscan(self):
        print "Subtracting overscan for %s" % self.imgpath
        hdulist = fits.open(self.imgpath)

        # To ensure we always traverse the amps in the same order, let's use
        # header information to figure out the proper indexing.  We'll define
        # index = 4*CCD + AMP - 3, where CCD in [1, 2, 3, 4] and AMP
        # in [0, 1, 2, 3] to keep consistent with header numbering.
        self.i = None
        for j in range(1, len(hdulist)):
            ccd = int(hdulist[j].header['CCDNAME'][-1])
            amp = int(hdulist[j].header['AMP-CFG'])
            if util.pair_to_hdu(ccd, amp) == self.hdu_num:
                self.i = j
        assert(self.i is not None)
        
        header = hdulist[self.i].header
        self.pHDU = hdulist[0].header
        self.iHDU = header
        self.detsec = header['DETSEC']
        self.ccdsec = header['CCDSEC']
        self.datasec = header['DATASEC']
        self.biassec = header['BIASSEC']
        self.ccdname = header['CCDNAME']
        self.ampcfg = header['AMP-CFG']

        newdata = hdulist[self.i].data.copy().T
        print "HDU %d's data has shape %s" % (self.i, newdata.shape)
        
        self.py_datasec = parse_fits_range(header['DATASEC'], to_python=True)
        self.py_biassec = parse_fits_range(header['BIASSEC'], to_python=True)
       
        print "Performing overscan subtraction with"
        print "CCDSEC = %s" % header['CCDSEC'] 
        print "DATASEC = %s" % header['DATASEC']
        print "BIASSEC = %s" % header['BIASSEC']

        biasreg = newdata[self.py_biassec][:, :]
        biasclip = stats.sigma_clip(biasreg, sig=3.0, cenfunc=np.mean, axis=0)
        biasmean = np.ma.filled(np.mean(biasclip, axis=0))
        newdata[self.py_datasec] -= biasmean[None, :]

        # The way we're estimating things, each row that shares the overscan
        # region will have the same error.
        row = np.sqrt(biasmean) / np.sqrt(biasreg.shape[0])
        newerr = np.repeat(row[np.newaxis, :], newdata.shape[0], 0)
        assert(newdata.shape == newerr.shape)
        
        self.data = newdata
        self.data_err = newerr

    def subtract_stackedbias(self, stackedbias, stackedbias_err):
        sb = fits.open(stackedbias)
        sb_err = fits.open(stackedbias_err)
        self.data -= sb[1].data.T
        #self.data_err = np.sqrt(self.data_err**2.0 + sb_err[1].data.T**2.0)
        # TODO: use real RMS error
        self.data_err = np.sqrt(self.data_err**2.0 + 7.0**2.0)

    def divide_stackedflat(self, stackedflat, stackedflat_err):
        sf = fits.open(stackedflat)
        sf_err = fits.open(stackedflat_err)
        self.data /= sf[1].data.T
        proc_err = sf_err[1].data.T * self.data
        self.data_err = np.sqrt(self.data_err**2.0 + proc_err**2.0)

    def update_mask(self, stackedflat_mask):
        sf = fits.open(stackedflat_mask)
        self.data_mask = sf[1].data.T

        #cdata_err = self.data_err - np.median(self.data_err)
        #ns = np.percentile(cdata_err, 15.87)
        #ps = np.percentile(cdata_err, 84.13)
        #sigma_err = max(np.abs(ns), np.abs(ps))
        #bad_pix = ((self.data_err == 0.0) | \
        #           (np.abs(self.data_err - np.median(self.data_err)) > 5.0*sigma_err)) & \
        #          (self.data_mask == 0)
        #self.data_mask[bad_pix] = 3

    def save(self):
        # For final, raw images, we'll only write out the data section.
        header = fits.Header()
        header['OBJECT'] = "processed_image"
        header['COMMENT'] = "Processed FITS file."
        header['DETSEC'] = self.detsec
        header['CCDSEC'] = self.ccdsec
        header['DATASEC'] = self.datasec
        header['CCDNAME'] = self.ccdname
        header['AMP-CFG'] = self.ampcfg
        header['RA'] = self.pHDU['RA']
        header['DEC'] = self.pHDU['DEC']
        header['EPOCH'] = self.pHDU['EPOCH']
        header['EQUINOX'] = self.pHDU['EQUINOX']
        header['FILTER'] = self.pHDU['FILTER']
        header['PIXSCAL1'] = self.pHDU['PIXSCAL1']
        header['PIXSCAL2'] = self.pHDU['PIXSCAL2']
        header['EXPTIME'] = self.pHDU['EXPTIME']
        header['DATE'] = self.pHDU['DATE']
        header['DATE-OBS'] = self.pHDU['DATE-OBS']
        header['CRVAL1'] = self.iHDU['CRVAL1']
        header['CRVAL2'] = self.iHDU['CRVAL2']
        header['CRPIX1'] = self.iHDU['CRPIX1']
        header['CRPIX2'] = self.iHDU['CRPIX2']
        header['CTYPE1'] = self.iHDU['CTYPE1']
        header['CTYPE2'] = self.iHDU['CTYPE2']

        # Transpose because of the way FITS handles indexing.
        imgHDU = fits.PrimaryHDU(data=self.data[self.py_datasec].T, header=header)
        weights = self.data_err**-2.0
        weights[self.data_mask != 0] = 0.0
        imgHDU_weight = fits.PrimaryHDU(data=weights[self.py_datasec].T, header=header)
        imgHDU_mask = fits.PrimaryHDU(data=self.data_mask[self.py_datasec].T, header=header)

        hdulist = fits.HDUList([imgHDU])
        basedir, filename = os.path.split(self.imgpath)
        basename, ext = os.path.splitext(filename)
        newpath = os.path.join(basedir, "%s.p.%02d.fits" % (basename, self.hdu_num))
        hdulist.writeto(newpath, clobber=True)

        hdulist_weight = fits.HDUList([imgHDU_weight])
        newpath_weight = os.path.join(basedir, "%s.p.%02d.weight.fits" % (basename, self.hdu_num))
        hdulist_weight.writeto(newpath_weight, clobber=True)

        hdulist_mask = fits.HDUList([imgHDU_mask])
        newpath_mask = os.path.join(basedir, "%s.p.%02d.mask.fits" % (basename, self.hdu_num))
        hdulist_mask.writeto(newpath_mask, clobber=True)
