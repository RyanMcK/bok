#!/usr/bin/env python

import argparse
import json
import os

from astropy.coordinates import SkyCoord
from astropy.io import fits
import astropy.units as u
import psycopg2
from psycopg2.extras import DictCursor

class Database(object):
    def __init__(self):
        # TODO: Maybe there's an alternative to hardcoded settings path?
        info = json.load(open("/global/u1/m/mckinnon/bok/database.json"))
        self.conn = psycopg2.connect(database=info['database'],
                                     user=info['user'],
                                     host=info['host'],
                                     password=info['password'])
        self.cursor = self.conn.cursor(cursor_factory=DictCursor)

    def import_images(self, rootdir):
        exts = ('.fits', '.fits.gz')
        for dirname, subdirs, files in os.walk(rootdir):
            imgs = [os.path.join(dirname, f) for f in files if f.endswith(exts)]
            map(self._import_single_image, imgs)
            self.conn.commit()   
 
    def _import_single_image(self, imgpath, overwrite=False):
        basedir, filename = os.path.split(imgpath)
        self.cursor.execute("""SELECT count(*) FROM bok_images WHERE basedir=%s
                            AND filename=%s""", (basedir, filename))
        img_exists = self.cursor.fetchone()[0]

        if img_exists and not overwrite:
            print "Skipping %s" % imgpath
            return
        elif img_exists and overwrite:
            print "Overwriting %s" % imgpath
            self.cursor.execute("""DELETE FROM bok_images WHERE basedir=%s
                                AND filename=%s""", (basedir, filename))

        print "Reading %s" % imgpath

        hdulist = fits.open(imgpath)
        header = hdulist[0].header
        header1 = hdulist[1].header

        origra = header.get('RA')
        origdec = header.get('DEC')
        if origra and origdec:
            coords = SkyCoord(origra, origdec, unit=(u.hour, u.degree))
            ra = coords.ra.degree
            dec = coords.dec.degree
        else:
            ra = dec = None
        filt = header.get('FILTER')
        obs_id = header.get('TELESCOP')
        exptime = header.get('EXPTIME')
        airmass = header.get('AIRMASS')
        try:
            date_obs = header.get('DATE-OBS') + 'T' + header.get('UTC-OBS')
        except TypeError as e:
            date_obs = None
        jd_obs = header.get('JULIAN')
        cd1_1 = header1.get('CD1_1')
        cd1_2 = header1.get('CD1_2')
        cd2_1 = header1.get('CD2_1')
        cd2_2 = header1.get('CD2_2')
        detsize = header1.get('DETSIZE')
        ccdsize = header1.get('CCDSIZE')
        biassec = header1.get('BIASSEC')
        datasec = header1.get('DATASEC')
        trimsec = header1.get('TRIMSEC')
        ampsec = header1.get('AMPSEC') 
        detsec = header1.get('DETSEC') 
        ccdsec = header1.get('CCDSEC') 
        ccdsec1 = header1.get('CCDSEC1') 

        self.cursor.execute("""INSERT INTO bok_images (ra, dec, filter,
                            observatory_id, exptime, airmass, date_obs, jd_obs, 
                            cd1_1, cd1_2, cd2_1, cd2_2, detsize, ccdsize, 
                            biassec, datasec, trimsec, ampsec, detsec, ccdsec, 
                            ccdsec1, filename, basedir) VALUES (%s, %s, %s, %s,
                            %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, 
                            %s, %s, %s, %s, %s, %s)""",
                           (ra, dec, filt, obs_id, exptime, airmass, date_obs,
                            jd_obs, cd1_1, cd1_2, cd2_1, cd2_2, detsize,
                            ccdsize, biassec, datasec, trimsec, ampsec, detsec,
                            ccdsec, ccdsec1, filename, basedir))
