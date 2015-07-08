#!/usr/bin/env python

import argparse
import glob
import os

import numpy as np

from database import Database
from rawimage import RawImage
from stackedbias import StackedBias
from stackedflat import StackedFlat

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--import_images", help="import images from directory")
    parser.add_argument("--stacked_bias", help="make stacked_bias from directory")
    parser.add_argument("--stacked_image", help="make stacked image from directory")
    parser.add_argument("--processed_image", help="make processed images from directory")

    parser.add_argument("--hdu", help="hdu number", required=True)
    args = parser.parse_args()
    hdu_num = int(args.hdu)

    if args.import_images:
        bok_dir = "/project/projectdirs/cosmo/staging/bok"
        db = Database()
        db.import_images(bok_dir)

    elif args.stacked_bias:
        filenames = glob.glob(os.path.join(args.stacked_bias, "d*.fits"))
        if len(filenames) > 100:
            np.random.shuffle(filenames)
            filenames = filenames[0:100]
        sb = StackedBias(filenames, hdu_num)
        sb.save()

    elif args.stacked_image:
        filenames = glob.glob(os.path.join(args.stacked_image, "d*.fits"))
        np.random.shuffle(filenames)
        si = StackedFlat(filenames[0:100], "/scratch2/scratchdirs/nugent/ryan_counts/mar_biases/stackedbias.%02d.fits" % hdu_num, "/scratch2/scratchdirs/nugent/ryan_counts/mar_biases/stackedbias.err.%02d.fits" % hdu_num, hdu_num)
        si.save()

    # At some point will want to rename the below data products to use %02d.
    elif args.processed_image:
        filenames = sorted(glob.glob("/scratch2/scratchdirs/nugent/ryan_counts/20150313/d*.fits"))
        sb = "/scratch2/scratchdirs/nugent/ryan_counts/mar_biases/stackedbias.%02d.fits" % hdu_num
        sb_err = "/scratch2/scratchdirs/nugent/ryan_counts/mar_biases/stackedbias.err.%02d.fits" % hdu_num
        sf = "/scratch2/scratchdirs/nugent/ryan_counts/mar_1000/stackedflat.%02d.fits" % hdu_num
        sf_err = "/scratch2/scratchdirs/nugent/ryan_counts/mar_1000/stackedflat.err.%02d.fits" % hdu_num
        sf_mask = "/scratch2/scratchdirs/nugent/ryan_counts/mar_1000/stackedflat.mask.%02d.fits" % hdu_num
        for filename in filenames:
            if ".p." in filename: continue
            ri = RawImage(filename, hdu_num)
            ri.subtract_overscan()
            ri.subtract_stackedbias(sb, sb_err)
            ri.divide_stackedflat(sf, sf_err)
            ri.update_mask(sf_mask)
            ri.save()

