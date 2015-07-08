#!/usr/bin/env python

import numpy as np

def mad(data, axis=None):
    return np.median(np.abs(data - np.median(data, axis)), axis)

# To ensure we always traverse the amps in the same order, let's use
# header information to figure out the proper indexing.  We'll define
# index = 4*CCD + AMP - 3, where CCD in [1, 2, 3, 4] and AMP
# in [0, 1, 2, 3] to keep consistent with header numbering.
def pair_to_hdu(ccd, amp):
    """Returns an integer between 1 and 16."""
    return 4*ccd + amp - 3

def hdu_to_pair(hdu):
    amp = (hdu + 3) % 4
    ccd = (hdu + 3 - amp) / 4
    return (ccd, amp)
