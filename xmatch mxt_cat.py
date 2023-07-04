#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  3 18:13:05 2023

@author: dt270490
"""

import crossmatch_utils
import astropy.units as u
from astropy.coordinates import SkyCoord

if __name__ == "__main__":

    cat_path = "/media/dt270490/9632505A32504203/workspace/SVOM/MXT/"
    cat_name = "MIC-X-ray-input-catalogue-v0.1.fits"
    
    cat, header = crossmatch_utils.get_cat_data(cat_path+cat_name) 
    
    mxt_skycoords=SkyCoord([233.75265,180.],[53.3437465,45.0],unit='degree')
    mxt_r90 = [1.0,5.0] #in arcmin
    data = crossmatch_utils.make_crossmatch(mxt_skycoords, mxt_r90, cat)
