#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  3 18:10:15 2023

@author: Damien Turpin
"""
import numpy as np
from astropy.io import fits
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.simbad import Simbad
Simbad.add_votable_fields('otype(opt)','plx', 'distance','typed_id')

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


def make_cat_crossmatch(mxt_skycoords,mxt_r90s,cat):
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
    cat_match : json
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
                            "ra":cat_src[1],
                            "dec":cat_src[2],
                            "r90":cat_src[3],
                            "r90_unit":"arcsec",
                            "l":cat_src[4],
                            "b":cat_src[5],
                            "ang_dist_tomxt": ang_dist[mask_match][k].arcsec,
                            "ang_dist_tomxt_unit": "arcsec"
                            }
            cat_src_matchs.append(cat_src_info)
            k = k+1
            
        cat_match['MXT_cand'+str(i)] = cat_src_matchs
    return cat_match

def get_simbad_galaxy_types():
    """
    Return the list of the galaxy taxonomy according to the Simbad DB

    Returns
    -------
    list
        DESCRIPTION.

    """
    
    return ["LSB","bCG","SBG","H2G","EmG","AGN", "SyG", "Sy1", "Sy2", "rG", 
            "LIN", "QSO", "Bla", "BLL", "IG", "PaG", "GrG", "CGG", "ClG", 
            "PCG", "SCG"]

def make_simbad_gal_conesearch(mxt_skycoords,mxt_r90s):
    """
    Perform the conesearch query to the Simbad database and select the nearby
    galaxies inside the MXT R90 errorbox

    Parameters
    ----------
    mxt_skycoords : astropy.coordinates.sky_coordinate.SkyCoord
        SkyCoordinates of the MXT transient candidates.
    mxt_r90s : list
        List of the MXT transient candidate error boxes.

    Returns
    -------
    simbad_galaxy_match : json
        List of the nearby galaxies having a positive match with MXT 
        candidates.

    """
    
    simbad_galaxy_match = {}
    for i in range(len(mxt_r90s)):
        sim_match = Simbad.query_region(mxt_skycoords[i], radius=mxt_r90s[i] * u.arcmin)
        #---- First step we mask all the objects for which we don't know
        #     the distance
        mask_dist = ~sim_match['PLX_VALUE'].mask | ~sim_match['Distance_distance'].mask
        sim_match_sel = sim_match[mask_dist]
        if not sim_match_sel:
            print("There is no known nearby galaxy within the MXT R90 error box")
            cat_gal_matchs = []
        else:
           #---- Second step we mask all the Simbad objects for which the OTYPE is not 
           # compatible with a galaxy type from the SIMBAD galaxy taxonomy
           mask_galaxy = []
           for obj in sim_match_sel:
               mask_galaxy.append(obj['OTYPE_opt_1'] in get_simbad_galaxy_types())
           gal_match = sim_match_sel[mask_galaxy]
           cat_gal_matchs = []
           for gal in gal_match:
               # compute the galaxy skycoords and the angular distance to the 
               # MXT position
               gal_coords = SkyCoord(gal["RA"],gal["DEC"],
                                   unit=(u.hourangle,u.deg))
               ang_dist = mxt_skycoords[i].separation(gal_coords)
               if type(gal['PLX_VALUE']) is np.float64:
                    dist = (1000./gal['PLX_VALUE'])/1e6 #in Mpc
                    dist_origin = 'parallax'
               else:
                   dist_origin = gal['Distance_bibcode']
                   if gal['Distance_unit'] == 'Mpc':
                       dist = gal['Distance_distance']
                   elif gal['Distance_unit'] == 'kpc':
                       dist = gal['Distance_distance']/1e3
                   elif gal['Distance_unit'] == 'pc':
                       dist = gal['Distance_distance']/1e6
                   elif gal['Distance_unit'] == 'Gpc':
                       dist = gal['Distance_distance']*1e3
               cat_gal_info = {"galaxy_name":gal["MAIN_ID"],
                               "cat_name":"Simbad",
                               "ra": gal_coords.ra.degree,
                               "dec":gal_coords.dec.degree,
                               "l":gal_coords.galactic.l,
                               "b":gal_coords.galactic.b,
                               "distance":dist,
                               "distance_unit":"Mpc",
                               "distance_source":dist_origin,
                               "ang_dist_tomxt": ang_dist.arcsec,
                               "ang_dist_tomxt_unit": "arcsec"
                               }
               cat_gal_matchs.append(cat_gal_info)
               print("The ",gal["MAIN_ID"]," galaxy ("+"%.2f" % dist+
                     " Mpc) is located inside the MXT R90 error box!")

        simbad_galaxy_match['MXT_cand'+str(i)] = cat_gal_matchs
    
    return simbad_galaxy_match
    