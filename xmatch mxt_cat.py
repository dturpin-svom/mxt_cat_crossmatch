#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  3 18:13:05 2023

@author: Damien Turpin (CEA-Saclay)
"""

import crossmatch_utils
import astropy.units as u
from astropy.coordinates import SkyCoord

if __name__ == "__main__":
    # loading cat path information
    cat_path = "/media/dt270490/9632505A32504203/workspace/SVOM/MXT/"
    cat_name = "MIC-X-ray-input-catalogue-v0.1.fits"
    # loading cat data with no duplicate
    cat, header = crossmatch_utils.get_cat_data(cat_path+cat_name) 
    # getting MXT transient candidate skyCoords and r90
    mxt_skycoords=SkyCoord([233.75265,148.8882194],[53.3437465,69.06529514],unit='degree')
    mxt_r90s = [1.0,1.0] #in arcmin
    #Perform the crossmatch with the x-ray source catalog 
    cat_src_match = crossmatch_utils.make_cat_crossmatch(mxt_skycoords, mxt_r90s, cat)
    
    # Perform the Simbad cone search for nearby galaxies
    simbad_gal_match = crossmatch_utils.make_simbad_gal_conesearch(mxt_skycoords,mxt_r90s)
    
