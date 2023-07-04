#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  3 18:10:15 2023

@author: Damien Turpin
"""

from astropy.io import fits
import astropy.units as u
from astropy.coordinates import SkyCoord

def get_cat_data(cat_path):
    """
    This code reads the x-ray source catalog and provide the cat info (data and
    header) in the cat and header arrays. It removes the duplicates based on
    the XMM name comparison

    Parameters
    ----------
    cat_path : string
        Path to the customed x-ray source catalog to be crossmatched with 
        QPO_MXT sources.

    Returns
    -------
    cat : astropy.io.fits.fitsrec.FITS_rec
        Catalog data.
    header : astropy.io.fits.header.Header
        Header of the catalog.

    """
    
    cat = fits.getdata(cat_path)
    indexes = [cat.XMM_NAME.tolist().index(x) for x in set(cat.XMM_NAME.tolist())]
    cat_noduplicate = cat[indexes]
    header = fits.getheader(cat_path)
    
    return cat_noduplicate, header


def make_qpomxt_conesearch(mxt_skycoords,mxt_r90s,cat):
    """
    This code performs the crossmatching of the mxt transient candidate
    position with the catalogued x-ray sources positions. It gives back the
    list cat sources that are compatible with errors with the MXT candidate 

    Parameters
    ----------
    mxt_skycoords : astropy.coordinates.sky_coordinate.SkyCoord
        SkyCoordinates of the MXT transient candidates.
    mxt_r90s : list
        List of the MXT transient candidate error boxes.
    cat : astropy.io.fits.fitsrec.FITS_rec
        Catalog data.

    Returns
    -------
    cat_match
        List of the catalogued sources with a positive match with MXT 
        candidates (json format).

    """
    
    cat_skycoords = SkyCoord(cat.RA, cat.DEC, unit='degree', frame='icrs')
    cat_match = {}
    for i in range(len(mxt_r90s)):
        ang_dist = mxt_skycoords[i].separation(cat_skycoords)
        mask_match = ang_dist<mxt_r90s[i]*u.arcmin+cat.ERROR*u.arcsec
        if mask_match.any():
            print('MXT candidate #',str(i),
                  '. There is at least one positive match with the x-ray catalog')
        else:
            print('MXT candidate #',str(i),
                  '. There is no match with the x-ray catalog')
        k = 0
        cat_src_matchs = []
        for cat_src in cat[mask_match]:
            cat_src_info = {"cat_src_name_xmm":cat_src[0],
                            "cat_src_name":cat_src[8],
                            "cat_name": cat_src[9],
                            "cat_src_ra":cat_src[1],
                            "cat_src_dec":cat_src[2],
                            "cat_src_r90":cat_src[3],
                            "cat_src_r90_unit":"arcsec",
                            "cat_src_l":cat_src[4],
                            "cat_src_b":cat_src[5],
                            "cat_src_disttomxt": ang_dist[mask_match][k].arcsec,
                            "cat_src_disttomxt_unit": "arcsec"
                            }
            cat_src_matchs.append(cat_src_info)
            k = k+1
            
        cat_match['MXT_cand'+str(i)] = cat_src_matchs
    return cat_match
    