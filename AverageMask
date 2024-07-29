#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 09:29:46 2024

@author: Ethan Bass

Used to create a composite image of each user mask and find a representative mask along with a measurement of disagreement among users. 
Saves the output data to a csv file. Requires all mosaics downloaded to run. Ethan has a copy of the output csv. 


"""

import traceback 
import ast 

import numpy as np
import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

plt.ion()
#get_ipython().run_line_magic('matplotlib', 'inline')
# some plots
#get_ipython().run_line_magic('matplotlib', 'qt')
# the interactive plot
from matplotlib.patches import Circle
from matplotlib.colors import LogNorm
from matplotlib.colors import SymLogNorm

#import astropy.units as u
from astropy.io import fits
from astropy import wcs
#from astropy.wcs import WCS
from astropy.io import ascii
#from astropy.coordinates import SkyCoord
#from astropy.coordinates import ICRS, Galactic, FK4, FK5
from astropy.visualization import make_lupton_rgb
#from astropy.modeling import models, fitting
from scipy import interpolate
#import itertools
import sys
import math
import csv
#import pylab as py
import copy
import cv2
import os
import pandas as pd
from astropy.nddata import Cutout2D

# These lines supress warnings
import warnings
warnings.filterwarnings('ignore')

#For beta testing "Pick up where you left off?" option

#You will need to make sure the following packages have been installed first:

#from tkinter import *
#import tkinter
#conda install -c anaconda tk

#from photutils import centroid_com
#https://photutils.readthedocs.io/en/stable/install.html
#conda install photutils -c astropy


path = '.'
catalog_name = os.path.join(path, 'USE_THIS_CATALOG_ybcat_MWP_with_ID.csv')

#Open the catalog file to get YB names, l, b, radius
data = ascii.read(catalog_name, delimiter=',')



#class choose_image: use l range to choose and open right mosaics
# returns .um8, .um8data, .um8w - full info, data array, wcs header at 8 um
# returns .um12, .um12data, .um12w - full info, data array, wcs header at 12 um
# returns .um24, .um24data, .um124w - full info, data array, wcs header at 24 um
# GWC adding mosaics beyond the pilot region starting 02/22/22    
class choose_image():
    def __init__(self, l, b):
        #currently don't need the WCS files for 12, 24 um because they
        #are reprojections onto 8um coordinate grid
        #GWC added 'b' on 4/5/22. 
        #Adding mosaics 021, 024, 027 on 10/17/23.
        if l > 19.5 and l <= 22.5:
            path8 = os.path.join(path,
                                 'mosaics/GLM_02100+0000_mosaic_I4.fits')
            path12 = os.path.join(path,
                                  'mosaics/WISE_12um_02100_mosaic.fits')
            path24 = os.path.join(
                path,
                'mosaics/MIPSGAL_24um_02100_mosaic.fits')
            path70 = os.path.join(
                path,
                'mosaics/PACS_70um_02100_mosaic.fits')    
        elif l > 22.5 and l <= 25.5:
            path8 = os.path.join(path,
                                 'mosaics/GLM_02400+0000_mosaic_I4.fits')
            path12 = os.path.join(path,
                                  'mosaics/WISE_12um_02400_mosaic.fits')
            path24 = os.path.join(
                path,
                'mosaics/MIPSGAL_24um_02400_mosaic.fits')
            path70 = os.path.join(
                path,
                'mosaics/PACS_70um_02400_mosaic.fits')    
        elif l > 25.5 and l <= 28.5:
            path8 = os.path.join(path,
                                 'mosaics/GLM_02700+0000_mosaic_I4.fits')
            path12 = os.path.join(path,
                                  'mosaics/WISE_12um_02700_mosaic.fits')
            path24 = os.path.join(
                path,
                'mosaics/MIPSGAL_24um_02700_mosaic.fits')
            path70 = os.path.join(
                path,
                'mosaics/PACS_70um_02700_mosaic.fits')    
        elif l > 28.5 and l <= 31.5:
            path8 = os.path.join(path,
                                 'mosaics/GLM_03000+0000_mosaic_I4.fits')
            path12 = os.path.join(path,
                                  'mosaics/WISE_12um_03000_mosaic.fits')
            path24 = os.path.join(
                path,
                'mosaics/MIPSGAL_03000_mosaic_reprojected.fits')
            path70 = os.path.join(
                path,
                'mosaics/PACS_70um_03000_mosaic.fits')
        elif l > 31.5 and l <= 34.5:
            path8 = os.path.join(path,
                                 'mosaics/GLM_03300+0000_mosaic_I4.fits')
            path12 = os.path.join(path,
                                  'mosaics/WISE_12um_03300_mosaic.fits')
            path24 = os.path.join(
                path,
                'mosaics/MIPSGAL_03300_mosaic_reprojected.fits')
            path70 = os.path.join(
                path,
                'mosaics/PACS_70um_03300_mosaic.fits')
        elif l > 34.5 and l <= 37.5:
            path8 = os.path.join(path,
                                 'mosaics/GLM_03600+0000_mosaic_I4.fits')
            path12 = os.path.join(path,
                                  'mosaics/WISE_12um_03600_mosaic.fits')
            path24 = os.path.join(
                path,
                'mosaics/MIPSGAL_03600_mosaic_reprojected.fits')
            path70 = os.path.join(
                path,
                'mosaics/PACS_70um_03600_mosaic.fits')
        elif l > 37.5 and l <= 40.5:
            path8 = os.path.join(path,
                                 'mosaics/GLM_03900+0000_mosaic_I4.fits')
            path12 = os.path.join(path,
                                  'mosaics/WISE_12um_03900_mosaic.fits')
            path24 = os.path.join(
                path,
                'mosaics/MIPSGAL_03900_mosaic_reprojected.fits')
            path70 = os.path.join(
                path,
                'mosaics/PACS_70um_03900_mosaic.fits')
        elif l > 40.5 and l <= 43.5:
            path8 = os.path.join(path,
                                 'mosaics/GLM_04200+0000_mosaic_I4.fits')
            path12 = os.path.join(path,
                                  'mosaics/WISE_12um_04200_mosaic.fits')
            path24 = os.path.join(
                path,
                'mosaics/MIPSGAL_24um_04200_mosaic.fits')
            path70 = os.path.join(
                path,
                'mosaics/PACS_70um_04200_mosaic.fits')
        elif l > 43.5 and l <= 46.5:
            path8 = os.path.join(path,
                                 'mosaics/GLM_04500+0000_mosaic_I4.fits')
            path12 = os.path.join(path,
                                  'mosaics/WISE_12um_04500_mosaic.fits')
            path24 = os.path.join(
                path,
                'mosaics/MIPSGAL_24um_04500_mosaic.fits')
            path70 = os.path.join(
                path,
                'mosaics/PACS_70um_04500_mosaic.fits')       
        elif l > 46.5 and l <= 49.5:
            path8 = os.path.join(path,
                                 'mosaics/GLM_04800+0000_mosaic_I4.fits')
            path12 = os.path.join(path,
                                  'mosaics/WISE_12um_04800_mosaic.fits')
            path24 = os.path.join(
                path,
                'mosaics/MIPSGAL_24um_04800_mosaic.fits')
            path70 = os.path.join(
                path,
                'mosaics/PACS_70um_04800_mosaic.fits') 
        elif l > 49.5 and l <= 52.5:
            path8 = os.path.join(path,
                                 'mosaics/GLM_05100+0000_mosaic_I4.fits')
            path12 = os.path.join(path,
                                  'mosaics/WISE_12um_05100_mosaic.fits')
            path24 = os.path.join(
                path,
                'mosaics/MIPSGAL_24um_05100_mosaic.fits')
            path70 = os.path.join(
                path,
                'mosaics/PACS_70um_05100_mosaic.fits')
        elif l > 52.5 and l <= 55.5:  
            path8 = os.path.join(path,
                                 'mosaics/GLM_05400+0000_mosaic_I4.fits')
            path12 = os.path.join(path,
                                  'mosaics/WISE_12um_05400_mosaic.fits')
            path24 = os.path.join(
                path,
                'mosaics/MIPSGAL_24um_05400_mosaic.fits')
            path70 = os.path.join(
                path,
                'mosaics/PACS_70um_05400_mosaic.fits')
        elif l > 55.5 and l <= 58.5:
            path8 = os.path.join(path,
                                 'mosaics/GLM_05700+0000_mosaic_I4.fits')
            path12 = os.path.join(path,
                                  'mosaics/WISE_12um_05700_mosaic.fits')
            path24 = os.path.join(
                path,
                'mosaics/MIPSGAL_24um_05700_mosaic.fits')
            path70 = os.path.join(
                    path,
                    'mosaics/PACS_70um_05700_mosaic.fits')   
        elif l > 58.5 and l <= 61.5:
            path8 = os.path.join(path,
                                     'mosaics/GLM_06000+0000_mosaic_I4.fits')
            path12 = os.path.join(path,
                                      'mosaics/WISE_12um_06000_mosaic.fits')
            path24 = os.path.join(
                    path,
                    'mosaics/MIPSGAL_24um_06000_mosaic.fits')
            path70 = os.path.join(
                    path,
                    'mosaics/PACS_70um_06000_mosaic.fits')  
        elif l > 61.5 and l <= 64.5:
            path8 = os.path.join(path,
                                     'mosaics/GLM_06300+0000_mosaic_I4.fits')
            path12 = os.path.join(path,
                                      'mosaics/WISE_12um_06300_mosaic.fits')
            path24 = os.path.join(
                    path,
                    'mosaics/MIPSGAL_24um_06300_mosaic.fits')
            path70 = os.path.join(
                    path,
                    'mosaics/PACS_70um_06300_mosaic.fits')                   
        elif l > 64.5 and l <= 65.5:
            path8 = os.path.join(path,
                                     'mosaics/GLM_06600+0000_mosaic_I4.fits')
            path12 = os.path.join(path,
                                      'mosaics/WISE_12um_06600_mosaic.fits')
            path24 = os.path.join(
                    path,
                    'mosaics/MIPSGAL_24um_06600_mosaic.fits')
            path70 = os.path.join(
                    path,
                    'mosaics/PACS_70um_06600_mosaic.fits')     
        elif l > 294.8 and l <= 295.5:
            path8 = os.path.join(path,
                                     'mosaics/GLM_29400+0000_mosaic_I4.fits')
            path12 = os.path.join(path,
                                      'mosaics/WISE_12um_29400_mosaic.fits')
            path24 = os.path.join(
                    path,
                    'mosaics/MIPSGAL_24um_29400_mosaic.fits')
            path70 = os.path.join(
                    path,
                    'mosaics/PACS_70um_29400_mosaic.fits')              
        elif l > 295.5 and l <= 298.5:
            path8 = os.path.join(path,
                                     'mosaics/GLM_29700+0000_mosaic_I4.fits')
            path12 = os.path.join(path,
                                      'mosaics/WISE_12um_29700_mosaic.fits')
            path24 = os.path.join(
                    path,
                    'mosaics/MIPSGAL_24um_29700_mosaic.fits')
            path70 = os.path.join(
                    path,
                    'mosaics/PACS_70um_29700_mosaic.fits')
        elif l > 298.5 and l <= 301.5:
            path8 = os.path.join(path,
                                     'mosaics/GLM_30000+0000_mosaic_I4.fits')
            path12 = os.path.join(path,
                                      'mosaics/WISE_12um_30000_mosaic.fits')
            path24 = os.path.join(
                    path,
                    'mosaics/MIPSGAL_24um_30000_mosaic.fits')
            path70 = os.path.join(
                    path,
                    'mosaics/PACS_70um_30000_mosaic.fits')
        elif l > 301.5 and l <= 304.5:
            path8 = os.path.join(path,
                                     'mosaics/GLM_30300+0000_mosaic_I4.fits')
            path12 = os.path.join(path,
                                      'mosaics/WISE_12um_30300_mosaic.fits')
            path24 = os.path.join(
                    path,
                    'mosaics/MIPSGAL_24um_30300_mosaic.fits')
            path70 = os.path.join(
                    path,
                    'mosaics/PACS_70um_30300_mosaic.fits')
        elif l > 304.5 and l <= 307.5:
            path8 = os.path.join(path,
                                     'mosaics/GLM_30600+0000_mosaic_I4.fits')
            path12 = os.path.join(path,
                                      'mosaics/WISE_12um_30600_mosaic.fits')
            path24 = os.path.join(
                    path,
                    'mosaics/MIPSGAL_24um_30600_mosaic.fits')
            path70 = os.path.join(
                    path,
                    'mosaics/PACS_70um_30600_mosaic.fits')       
        #Added these mosaics on 7/11/23. For some reason, many of the elif statements
        #for regions I've done are not here. I added these above on 7/26/23. I had to
        #copy them over from ExpertPhotom.py
        #Adding more as I complete photometry for eacj sector.
        elif l > 307.5 and l <= 310.5:
            path8 = os.path.join(path,
                                     'mosaics/GLM_30900+0000_mosaic_I4.fits')
            path12 = os.path.join(path,
                                      'mosaics/WISE_12um_30900_mosaic.fits')
            path24 = os.path.join(
                    path,
                    'mosaics/MIPSGAL_24um_30900_mosaic.fits')
            path70 = os.path.join(
                    path,
                    'mosaics/PACS_70um_30900_mosaic.fits')
        elif l > 310.5 and l <= 313.5:
            path8 = os.path.join(path,
                                 'mosaics/GLM_31200+0000_mosaic_I4.fits')
            path12 = os.path.join(path,
                                  'mosaics/WISE_12um_31200_mosaic.fits')
            path24 = os.path.join(
                        path,
                        'mosaics/MIPSGAL_24um_31200_mosaic.fits')
            path70 = os.path.join(
                        path,
                        'mosaics/PACS_70um_31200_mosaic.fits')
        elif l > 313.5 and l <= 316.5:
            path8 = os.path.join(path,
                                  'mosaics/GLM_31500+0000_mosaic_I4.fits')
            path12 = os.path.join(path,
                                   'mosaics/WISE_12um_31500_mosaic.fits')
            path24 = os.path.join(
                         path,
                         'mosaics/MIPSGAL_24um_31500_mosaic.fits')
            path70 = os.path.join(
                         path,
                         'mosaics/PACS_70um_31500_mosaic.fits')    
        elif l > 316.5 and l <= 319.5:
            path8 = os.path.join(path,
                                   'mosaics/GLM_31800+0000_mosaic_I4.fits')
            path12 = os.path.join(path,
                                    'mosaics/WISE_12um_31800_mosaic.fits')
            path24 = os.path.join(
                          path,
                          'mosaics/MIPSGAL_24um_31800_mosaic.fits')
            path70 = os.path.join(
                          path,
                          'mosaics/PACS_70um_31800_mosaic.fits')      
        elif l > 319.5 and l <= 322.5:
            path8 = os.path.join(path,
                                   'mosaics/GLM_32100+0000_mosaic_I4.fits')
            path12 = os.path.join(path,
                                    'mosaics/WISE_12um_32100_mosaic.fits')
            path24 = os.path.join(
                          path,
                          'mosaics/MIPSGAL_24um_32100_mosaic.fits')
            path70 = os.path.join(
                          path,
                          'mosaics/PACS_70um_32100_mosaic.fits')   
        elif l > 322.5 and l <= 325.5:
            path8 = os.path.join(path,
                                   'mosaics/GLM_32400+0000_mosaic_I4.fits')
            path12 = os.path.join(path,
                                    'mosaics/WISE_12um_32400_mosaic.fits')
            path24 = os.path.join(
                          path,
                          'mosaics/MIPSGAL_24um_32400_mosaic.fits')
            path70 = os.path.join(
                          path,
                          'mosaics/PACS_70um_32400_mosaic.fits')          
        elif l > 325.5 and l <= 328.5:
            path8 = os.path.join(path,
                                   'mosaics/GLM_32700+0000_mosaic_I4.fits')
            path12 = os.path.join(path,
                                    'mosaics/WISE_12um_32700_mosaic.fits')
            path24 = os.path.join(
                          path,
                          'mosaics/MIPSGAL_24um_32700_mosaic.fits')
            path70 = os.path.join(
                          path,
                          'mosaics/PACS_70um_32700_mosaic.fits')         
        elif l > 328.5 and l <= 331.5:
            path8 = os.path.join(path,
                                       'mosaics/GLM_33000+0000_mosaic_I4.fits')
            path12 = os.path.join(path,
                                        'mosaics/WISE_12um_33000_mosaic.fits')
            path24 = os.path.join(
                              path,
                              'mosaics/MIPSGAL_24um_33000_mosaic.fits')
            path70 = os.path.join(
                              path,
                              'mosaics/PACS_70um_33000_mosaic.fits')   
        elif l > 331.5 and l <= 334.5:
            path8 = os.path.join(path,
                                       'mosaics/GLM_33300+0000_mosaic_I4.fits')
            path12 = os.path.join(path,
                                        'mosaics/WISE_12um_33300_mosaic.fits')
            path24 = os.path.join(
                              path,
                              'mosaics/MIPSGAL_24um_33300_mosaic.fits')
            path70 = os.path.join(
                              path,
                              'mosaics/PACS_70um_33300_mosaic.fits')           
        elif l > 334.5 and l <= 337.5:
            path8 = os.path.join(path,
                                    'mosaics/GLM_33600+0000_mosaic_I4.fits')
            path12 = os.path.join(path,
                                  'mosaics/WISE_12um_33600_mosaic.fits')
            path24 = os.path.join(
                                  path,
                                  'mosaics/MIPSGAL_24um_33600_mosaic.fits')
            path70 = os.path.join(
                                  path,
                                  'mosaics/PACS_70um_33600_mosaic.fits')  
        elif l > 337.5 and l <= 340.5:
            path8 = os.path.join(path,
                                     'mosaics/GLM_33900+0000_mosaic_I4.fits')
            path12 = os.path.join(path,
                                   'mosaics/WISE_12um_33900_mosaic.fits')
            path24 = os.path.join(
                                   path,
                                   'mosaics/MIPSGAL_24um_33900_mosaic.fits')
            path70 = os.path.join(
                                   path,
                                   'mosaics/PACS_70um_33900_mosaic.fits')        
        elif l > 340.5 and l <= 343.5:
            path8 = os.path.join(path,
                                     'mosaics/GLM_34200+0000_mosaic_I4.fits')
            path12 = os.path.join(path,
                                   'mosaics/WISE_12um_34200_mosaic.fits')
            path24 = os.path.join(
                                   path,
                                   'mosaics/MIPSGAL_24um_34200_mosaic.fits')
            path70 = os.path.join(
                                   path,
                                   'mosaics/PACS_70um_34200_mosaic.fits') 
        elif l > 343.5 and l <= 346.5:
            path8 = os.path.join(path,
                                 'mosaics/GLM_34500+0000_mosaic_I4.fits')
            path12 = os.path.join(path,
                                  'mosaics/WISE_12um_34500_mosaic.fits')
            path24 = os.path.join(
                                   path,
                                      'mosaics/MIPSGAL_24um_34500_mosaic.fits')
            path70 = os.path.join(
                                   path,
                                       'mosaics/PACS_70um_34500_mosaic.fits') 
        elif l > 346.5 and l <= 349.5:
            path8 = os.path.join(path,
                                'mosaics/GLM_34800+0000_mosaic_I4.fits')
            path12 = os.path.join(path,
                                 'mosaics/WISE_12um_34800_mosaic.fits')
            path24 = os.path.join(
                                  path,
                                     'mosaics/MIPSGAL_24um_34800_mosaic.fits')
            path70 = os.path.join(
                                  path,
                                      'mosaics/PACS_70um_34800_mosaic.fits') 
        elif l > 349.5 and l <= 352.5:
            path8 = os.path.join(path,
                                'mosaics/GLM_35100+0000_mosaic_I4.fits')
            path12 = os.path.join(path,
                                 'mosaics/WISE_12um_35100_mosaic.fits')
            path24 = os.path.join(
                                  path,
                                     'mosaics/MIPSGAL_24um_35100_mosaic.fits')
            path70 = os.path.join(
                                  path,
                                      'mosaics/PACS_70um_35100_mosaic.fits') 
        elif l > 352.5 and l <= 355.5:
            path8 = os.path.join(path,
                                'mosaics/GLM_35400+0000_mosaic_I4.fits')
            path12 = os.path.join(path,
                                 'mosaics/WISE_12um_35400_mosaic.fits')
            path24 = os.path.join(
                                  path,
                                     'mosaics/MIPSGAL_24um_35400_mosaic.fits')
            path70 = os.path.join(
                                  path,
                                      'mosaics/PACS_70um_35400_mosaic.fits') 
        elif l > 355.5 and l <= 358.5:
            path8 = os.path.join(path,
                                 'mosaics/GLM_35700+0000_mosaic_I4.fits')
            path12 = os.path.join(path,
                                  'mosaics/WISE_12um_35700_mosaic.fits')
            path24 = os.path.join(
                                   path,
                                      'mosaics/MIPSGAL_24um_35700_mosaic.fits')
            path70 = os.path.join(
                                   path,
                                       'mosaics/PACS_70um_35700_mosaic.fits')    
        elif (l > 358.5 and l <= 360.1) or (l > -0.1 and l <= 1.5):
            path8 = os.path.join(path,
                                 'mosaics/GLM_00000+0000_mosaic_I4.fits')
            path12 = os.path.join(path,
                                  'mosaics/WISE_12um_00000_mosaic.fits')
            path24 = os.path.join(
                                   path,
                                      'mosaics/MIPSGAL_24um_00000_mosaic.fits')
            path70 = os.path.join(
                                   path,
                                       'mosaics/PACS_70um_00000_mosaic.fits')
        elif l > 1.5 and l <= 4.5:
            path8 = os.path.join(path,
                                 'mosaics/GLM_00300+0000_mosaic_I4.fits')
            path12 = os.path.join(path,
                                  'mosaics/WISE_12um_00300_mosaic.fits')
            path24 = os.path.join(
                                   path,
                                      'mosaics/MIPSGAL_24um_00300_mosaic.fits')
            path70 = os.path.join(
                                   path,
                                       'mosaics/PACS_70um_00300_mosaic.fits')
        elif l > 4.5 and l <= 7.5:
            path8 = os.path.join(path,
                                  'mosaics/GLM_00600+0000_mosaic_I4.fits')
            path12 = os.path.join(path,
                                   'mosaics/WISE_12um_00600_mosaic.fits')
            path24 = os.path.join(
                                    path,
                                       'mosaics/MIPSGAL_24um_00600_mosaic.fits')
            path70 = os.path.join(
                                    path,
                                        'mosaics/PACS_70um_00600_mosaic.fits')
        elif l > 7.5 and l <= 10.5:
            path8 = os.path.join(path,
                                  'mosaics/GLM_00900+0000_mosaic_I4.fits')
            path12 = os.path.join(path,
                                   'mosaics/WISE_12um_00900_mosaic.fits')
            path24 = os.path.join(
                                    path,
                                       'mosaics/MIPSGAL_24um_00900_mosaic.fits')
            path70 = os.path.join(
                                    path,
                                        'mosaics/PACS_70um_00900_mosaic.fits')
        elif l > 10.5 and l <= 13.5:
            path8 = os.path.join(path,
                                  'mosaics/GLM_01200+0000_mosaic_I4.fits')
            path12 = os.path.join(path,
                                   'mosaics/WISE_12um_01200_mosaic.fits')
            path24 = os.path.join(
                                    path,
                                       'mosaics/MIPSGAL_24um_01200_mosaic.fits')
            path70 = os.path.join(
                                    path,
                                        'mosaics/PACS_70um_01200_mosaic.fits')
        elif l > 13.5 and l <= 16.5:
            path8 = os.path.join(path,
                                  'mosaics/GLM_01500+0000_mosaic_I4.fits')
            path12 = os.path.join(path,
                                   'mosaics/WISE_12um_01500_mosaic.fits')
            path24 = os.path.join(
                                    path,
                                       'mosaics/MIPSGAL_24um_01500_mosaic.fits')
            path70 = os.path.join(
                                    path,
                                        'mosaics/PACS_70um_01500_mosaic.fits')
        elif l > 16.5 and l <= 19.5:
            path8 = os.path.join(path,
                                  'mosaics/GLM_01800+0000_mosaic_I4.fits')
            path12 = os.path.join(path,
                                   'mosaics/WISE_12um_01800_mosaic.fits')
            path24 = os.path.join(
                                    path,
                                       'mosaics/MIPSGAL_24um_01800_mosaic.fits')
            path70 = os.path.join(
                                    path,
                                        'mosaics/PACS_70um_01800_mosaic.fits')
        #The following are for the SMOG region.  
        #GWC: Something went wonky on 2/7/24 -- need to revisit how to cover SMOG.
        elif l > 101.0 and l <= 105.59 and b < 3.06:
            path8 = os.path.join(path,
                                 'mosaics/SMOG_08um_10300_mosaic.fits')
            path12 = os.path.join(path,
                                  'mosaics/SMOG_12um_10300_mosaic.fits')
            path24 = os.path.join(
                path,
                'mosaics/SMOG_24um_10300_mosaic.fits')
            path70 = os.path.join(
                path,
                'mosaics/SMOG_PACS_70um_10300_mosaic.fits')
            # Replaced 'mosaics/SMOG_70um_10300_mosaic.fits') with PACS on 7/7/23
        elif l > 101.0 and l <= 105.59 and b >= 3.06:
            path8 = os.path.join(path,
                                 'mosaics/SMOG_08um_10300_mosaic.fits')
            path12 = os.path.join(path,
                                  'mosaics/SMOG_12um_10300_mosaic.fits')
            path24 = os.path.join(
                path,
                'mosaics/SMOG_24um_10300_mosaic_high_b.fits')
            path70 = os.path.join(
                path,
                'mosaics/SMOG_PACS_70um_10300_mosaic.fits')
        elif l > 105.59 and l <= 110.2:
            path8 = os.path.join(path,
                                 'mosaics/SMOG_08um_10700_mosaic.fits')
            path12 = os.path.join(path,
                                  'mosaics/SMOG_12um_10700_mosaic.fits')
            path24 = os.path.join(
                path,
                'mosaics/SMOG_24um_10700_mosaic.fits')
            path70 = os.path.join(
                path,
                'mosaics/SMOG_PACS_70um_10700_mosaic.fits')
            # Replaced 'mosaics/SMOG_70um_10700_mosaic.fits') with PACS on 7/7/23
        # The following were added for Cyg-X by GW-C on 2/7/24.
        elif l > 75.5 and l <= 76.5:
            path8 = os.path.join(path,
                                 'mosaics/CYGX_08um_07500+0050_mosaic.fits')
            path12 = os.path.join(path,
                                  'mosaics/CYGX_12um_07500+0050_mosaic.fits')
            path24 = os.path.join(
                path,
                'mosaics/CYGX_24um_07500+0050_mosaic.fits')
            path70 = os.path.join(
                path,
                'mosaics/CYGX_70um_07500+0050_mosaic.fits')
        elif l > 76.5 and l <= 79.5 and b < 0.82:
            path8 = os.path.join(path,
                                 'mosaics/CYGX_08um_07800-0085_mosaic.fits')
            path12 = os.path.join(path,
                                  'mosaics/CYGX_12um_07800-0085_mosaic.fits')
            path24 = os.path.join(
                path,
                'mosaics/CYGX_24um_07800-0085_mosaic.fits')
            path70 = os.path.join(
                path,
                'mosaics/CYGX_70um_07800-0085_mosaic.fits')
        elif l > 76.5 and l <= 79.5 and b >= 0.82:
            path8 = os.path.join(path,
                                 'mosaics/CYGX_08um_07800+0250_mosaic.fits')
            path12 = os.path.join(path,
                                  'mosaics/CYGX_12um_07800+0250_mosaic.fits')
            path24 = os.path.join(
                path,
                'mosaics/CYGX_24um_07800+0250_mosaic.fits')
            path70 = os.path.join(
                path,
                'mosaics/CYGX_70um_07800+0250_mosaic.fits')
        elif l > 79.5 and l <= 82.5 and b < 0.82:
            path8 = os.path.join(path,
                                 'mosaics/CYGX_08um_08100-0070_mosaic.fits')
            path12 = os.path.join(path,
                                  'mosaics/CYGX_12um_08100-0070_mosaic.fits')
            path24 = os.path.join(
                path,
                'mosaics/CYGX_24um_08100-0070_mosaic.fits')
            path70 = os.path.join(
                path,
                'mosaics/CYGX_70um_08100-0070_mosaic.fits')
        elif l > 79.5 and l <= 82.5 and b >= 0.82:
            path8 = os.path.join(path,
                                 'mosaics/CYGX_08um_08100+0235_mosaic.fits')
            path12 = os.path.join(path,
                                  'mosaics/CYGX_12um_08100+0235_mosaic.fits')
            path24 = os.path.join(
                path,
                'mosaics/CYGX_24um_08100+0235_mosaic.fits')
            path70 = os.path.join(
                path,
                'mosaics/CYGX_70um_08100+0235_mosaic.fits')
        elif l > 82.5 and l <= 83.0:
            path8 = os.path.join(path,
                                 'mosaics/CYGX_08um_08400+0005_mosaic.fits')
            path12 = os.path.join(path,
                                  'mosaics/CYGX_12um_08400+0005_mosaic.fits')
            path24 = os.path.join(
                path,
                'mosaics/CYGX_24um_08400+0005_mosaic.fits')
            path70 = os.path.join(
                path,
                'mosaics/CYGX_70um_08400+0005_mosaic.fits')
        else:
            # GWC revised print statement from "outside the pilot..."
            print('Your YB is outside the region.')
            print('Please try again.')
            sys.exit()

        temp = fits.open(path8)[0]
        self.um8 = temp
        self.um8data = temp.data
        self.um8w = wcs.WCS(temp.header)
        temp = fits.open(path12)[0]
        self.um12 = temp
        self.um12data = temp.data
        self.um12w = wcs.WCS(temp.header)
        temp = fits.open(path24)[0]
        self.um24 = temp
        self.um24data = temp.data
        self.um24w = wcs.WCS(temp.header)
        temp = fits.open(path70)[0]
        self.um70 = temp
        self.um70data = temp.data
        self.um70w = wcs.WCS(temp.header)

#class that does the masking and interpolation, returns masked, blanked, interpolated, and residual
class do_interp():
    def __init__(self, img, verts):

        #use the clicked values from the user to create a NaN mask
        vertices = np.array([verts], dtype=np.int32)
        xyvals = np.array(verts, dtype=np.int32)
        xmin = min(xyvals[:, 0]) - 5
        xmax = max(xyvals[:, 0]) + 5
        ymin = min(xyvals[:, 1]) - 5
        ymax = max(xyvals[:, 1]) + 5
        #print(xmin, xmax, ymin, ymax)
        mask = np.zeros_like(img)
        inverse_mask = np.zeros_like(img)
        region_mask = np.zeros_like(img)
        cutout = np.zeros_like(img)

        # filling pixels inside the polygon defined by "vertices" with the fill color
        cv2.fillPoly(mask, vertices, 255)
        #TURN ALL Non-ZERO to NaN
        inverse_mask[np.nonzero(mask)] = int(
            1)  # ones inside poly, zero outside

        mask[np.nonzero(mask)] = float('nan')
        #TURN ALL ZERO to 1
        mask[np.where(mask == 0)] = int(1)  # nan inside poly, 1 outside
        region_mask = mask
        region_mask = np.nan_to_num(region_mask)  # zero in poly, 1 outside
        cutout[ymin:ymax, xmin:xmax] = mask[ymin:ymax, xmin:xmax]
        #TURN EVERYTHING OUTSIDE THAT RANGE to NaN
        cutout[np.where(cutout == 0)] = float('nan')

        #TAKE image=workask*mask will make a image with original values but NaN in polygon
        #blank = img*mask
        blank = img * cutout
        self.masked = mask
        self.blanked = blank
        goodvals = np.where(np.isfinite(blank))

        #perform the interpolation over the masked coordinates
        x = goodvals[1]  # x values of finite coordinates
        y = goodvals[0]  # y values of finite coordinates

        def get_fvals(x, y):
            range_array = np.arange(x.size)
            vals = np.zeros(x.size)
            for (i, xi, yi) in zip(range_array, x, y):
                vals[i] = img[yi][xi]
            return vals

        fvals = get_fvals(x, y)

        newfunc = interpolate.Rbf(
            x, y, fvals,
            function='multiquadric')  # the function that does interpolation
        allvals = np.where(img)  # whole region to interpolate over
        xnew = allvals[1]
        ynew = allvals[0]
        fnew = newfunc(xnew, ynew)

        #put the interpolated values back into a 2D array for display and other uses
        def make_2D(fnew, xnew, ynew, img):
            new_array = np.zeros(
                (int(xnew.size /
                     ((img.shape)[0])), int(ynew.size / ((img.shape)[1]))),
                dtype=float)
            #print(new_array)
            range_array = np.arange(fnew.size)
            #print("rangearay:",range_array)

            for (i, x, y) in zip(range_array, xnew, ynew):
                new_array[y][x] = fnew[i]

            return new_array

        fnew_2D = make_2D(fnew, xnew, ynew, img)

        self.interp = img * region_mask + fnew_2D * inverse_mask

        #generate the residual image (original - interpolated background)
        self.resid = img - (img * region_mask + fnew_2D * inverse_mask)


class get_flux():
    def __init__(self, d8, d12, d24, d70):
        #reset the total flux value
        #flux_tot8 = 0
        #flux_tot12 = 0
        #flux_tot24 = 0
        #flux_tot70 = 0
        #go through 100 x 100 pixel residual image and add up the pixel values
        flux_tot8 = sum(map(sum,d8))
        flux_tot12 = sum(map(sum,d12))
        flux_tot24 = sum(map(sum,d24))
        flux_tot70 = sum(map(sum,d70))
        
        #for ROW in range(0, 100):
        #    for column in range(0, 100):

        #        flux_tot8 = flux_tot8 + d8[ROW][
        #            column]  #gives the value of photometry
        #        flux_tot12 = flux_tot12 + d12[ROW][
        #            column]  #gives the value of photometry
        #        flux_tot24 = flux_tot24 + d24[ROW][
        #            column]  #gives the value of photometry
        #        flux_tot70 = flux_tot70 + d70[ROW][
        #            column]  #gives the value of photometry

            #convert units of flux total.  MIPS/IRAC in MJy/Sr*Pixel, want Jy
            #conversion: (MJy/Sr*Pixel)*(Sr/Pix)*(Jy/MJy)
            #WISE is in units of data number per pixel
            #Conversion: DN/Pix*Pixel*Jy/DN

        pixscale8 = 0.000333333
        #pixscale8 = abs(image.um8w.wcs.cdelt1)
        #print(pixscale8)
        #8 um square degree to square pixel conversion-- x*y
        #GWC edit 4/1/22: Set sqdeg_tosqpix8 directly to pixscale8 * pixscale8
        #sqdeg_tosqpix8 = abs(image.um8w.wcs.cd[0][0]) * abs(
        #    image.um8w.wcs.cd[1][1])
        sqdeg_tosqpix8 = pixscale8 * pixscale8
        #8 um steradian to pixel conversion (used to convert MJy/Pixel)
        #     will be used to convert to MJy
        str_to_pix8 = sqdeg_tosqpix8 * 0.0003046118
        dn_to_jy12 = 1.83 * 10**-6
        
        
        flux_tot8 = flux_tot8 * str_to_pix8 * 10**6
        flux_tot12 = flux_tot12 * dn_to_jy12
        flux_tot24 = flux_tot24 * str_to_pix8 * 10**6
        flux_tot70 = flux_tot70 #No conversion necesary

        self.um8 = flux_tot8
        self.um12 = flux_tot12
        self.um24 = flux_tot24
        self.um70 = flux_tot70


class AverageMask:
    def __init__(self):
        self.results = pd.DataFrame(columns=["YB", 
                                             "AvMaskFlux70", "AvMaskFlux24", "AvMaskFlux12", "AvMaskFlux8",
                                             "AvMaskEr70", "AvMaskEr24", "AvMaskEr12", "AvMaskEr8"])
    
    def show_mask(self, YBnum):
        
        
        #get the YB's location and radius
        print(YBnum, data['YB'][YBnum-1])
        YB_long = data[YBnum-1]['l']
        YB_lat = data[YBnum-1]['b']
        YB_rad = data[YBnum-1]['r']
    
        #Use the location to determine the correct image files to use
        #GWC added YB_lat on 4/5/22.
        image = choose_image(YB_long, YB_lat)
    
        ##          Unit Conversions Info         ##
        # 8 um:
        #GLIMPSE 8 um flux units are MJy/steradian
        #obtain the 8 um pixel scale from the image header
        #this gives the degree/pixel conversion factor used for overdrawing circle
        #and flux conversions (this is x direction only but they are equal in this case)
        #GWC edit 4/1/22: Set pixscale directly to 0.000333333
        #pixscale8 = abs(image.um8w.wcs.cd[0][0])
        pixscale8 = 0.000333333
        #pixscale8 = abs(image.um8w.wcs.cdelt1)
        #print(pixscale8)
        #8 um square degree to square pixel conversion-- x*y
        #GWC edit 4/1/22: Set sqdeg_tosqpix8 directly to pixscale8 * pixscale8
        #sqdeg_tosqpix8 = abs(image.um8w.wcs.cd[0][0]) * abs(
        #    image.um8w.wcs.cd[1][1])
        sqdeg_tosqpix8 = pixscale8 * pixscale8
        #8 um steradian to pixel conversion (used to convert MJy/Pixel)
        #     will be used to convert to MJy
        str_to_pix8 = sqdeg_tosqpix8 * 0.0003046118
        #WISE Units are Data Number per Pixel, Jy/DN is 1.8326*10^-6
        #See http://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec2_3f.html
        dn_to_jy12 = 1.83 * 10**-6
    
        #convert YB l and b and radius to pixel coordinates
        ybwcs = np.array([[YB_long, YB_lat]], np.float_)
        pixcoords = image.um8w.wcs_world2pix(ybwcs, 1)
        YB_long_pix = pixcoords[0][0]
        YB_lat_pix = pixcoords[0][1]
        YB_rad_pix = YB_rad / pixscale8
    
        # This is for the added points to show the user what other YBs are in the images
        # Read in the l, b, and r values for all the YBs and convert them to pixels
        YBloc = pd.read_csv(catalog_name, usecols = ['l', 'b', 'r'])
        # Convert l, b, and r from YBloc into pixels
        for i in range(0, len(YBloc)):
            yblocwcs = np.array([[YBloc['l'][i], YBloc['b'][i]]], np.float_)
            pixcoordsloc = image.um8w.wcs_world2pix(yblocwcs, 1)
            YB_l = pixcoordsloc[0][0]
            YBloc['l'][i] = YB_l
            YB_b = pixcoordsloc[0][1]
            YBloc['b'][i] = YB_b
            YB_radloc = YBloc['r'][i]
            YB_r = YB_radloc / pixscale8
            YBloc['r'][i] = YB_r
    
        #Create cropped 100 x 100 pixel image arrays centered on YB
        #orig = image.um8data[y1:y2,x1:x2]
        #orig12 = image.um12data[y1:y2,x1:x2]
        #orig24 = image.um24data[y1:y2,x1:x2]
    
        #use Cutout2D to make the zoomed windows
        position = (YB_long_pix + 0.5, YB_lat_pix + 0.5)
        size = (100, 100)
    
        cut8 = Cutout2D(data=image.um8data,
                        position=position,
                        size=size,
                        wcs=image.um8w)
        cut12 = Cutout2D(data=image.um12data,
                         position=position,
                         size=size,
                         wcs=image.um12w)
        cut24 = Cutout2D(data=image.um24data,
                         position=position,
                         size=size,
                         wcs=image.um24w)
        cut70 = Cutout2D(data=image.um70data,
                         position=position,
                         size=size,
                         wcs=image.um70w)
    
        fitcopy8 = image.um8
        fitcopy8.data = cut8.data
        fitcopy8.header.update(cut8.wcs.to_header())
    
        fitcopy12 = image.um12
        fitcopy12.data = cut12.data
        fitcopy12.header.update(cut12.wcs.to_header())
    
        fitcopy24 = image.um24
        fitcopy24.data = cut24.data
        fitcopy24.header.update(cut24.wcs.to_header())
    
        fitcopy70 = image.um70
        fitcopy70.data = cut70.data
        fitcopy70.header.update(cut70.wcs.to_header())
    
        orig = cut8.data
        orig12 = cut12.data
        orig24 = cut24.data
        orig70 = cut70.data
    
        #create copies of cropped images called workmasks
        workmask = copy.deepcopy(orig)
        workmask12 = copy.deepcopy(orig12)
        workmask24 = copy.deepcopy(orig24)
        workmask70 = copy.deepcopy(orig70)
        
        LinearThreshold = 0.001
        
        ###################################
        #
        ## Here is where we start to actually differ from the main code
        #
        ###################################
        
        #fig = plt.figure(figsize=(8, 8))

        ####################################################################################################################################################
        #
        #This is the list of the names of all the csv's of points that we are using
        #
        #Update this list with the csv files you want to use for the average mask. Note that this section is to generate the average mask for a specific YB.  
        #
        ####################################################################################################################################################
        
        csvlist = ['YBphotometry_results_EthanBass.csv', 'YBphotometry_results_Ignasi.csv', 'YBphotometry_results_colemaya.csv', 'YBphotometry_results_IgnasiBonmatiGonzalvez.csv', 'YBphotometry_results_KORRA.csv', 'YBphotometry_results_WolfChase1-ALL_YBs.csv', 'YBphotometry_results_Asta.csv']
        
    
    
        #csvlist = ['YBphotometry_results_EthanBass.csv']
        
        headerlist = ['vertices 70', 'vertices 24', 'vertices 12', 'vertices 8']
        imagelist = [workmask70, workmask24, workmask12, workmask]
        averagelist = []
        residlist = []

        for i in range(0, 4):
            average_mask = np.zeros_like(workmask70)
            img = imagelist[i]
            header = headerlist[i]
            validlength = 0
            for csvname in csvlist:
                userdata = ascii.read(csvname, delimiter=',')
                if userdata[YBnum][header] != '':
                    vertices = ast.literal_eval(userdata[YBnum][header])
                    interp = do_interp(img, vertices)

                    im = interp.masked
                    im[np.where(im != 1)] = int(0)
                    average_mask += im
                    validlength += 1
                #else:
                    #print(csvname)
                    
                
            print(f"Num of sources: {validlength}, for YB {YBnum}, at {headerlist[i]}um")
            if (validlength == 0):
                return 0
            average_mask /= validlength
            #print(average_mask[50])
            averagelist.append(np.copy(average_mask)*-1 +1)
    
                #here average_mask is 1 outside, 0 inside of all shapes, and some fractions between
            sensitivity = 0.5 #determines the breakoff point
                # a low sensitivity means that we take the smaller shape
                # a high sensitivity means that we take bigger shapes'', 'outer', 'inner', 'flatinner'
                # thus a low sensitivity will be more 'zoomed in' on the YB
            
            

            average_mask[np.where(average_mask >= sensitivity)] = 1
            average_mask[np.where(average_mask <  sensitivity)] = float('nan')
            #print(average_mask[50])
            region_mask = np.copy(average_mask)
            region_mask[np.where(region_mask != 1)] = 0
                
            inverse_mask = region_mask*-1 +1
                
                #In order to recreate goodvals, I need to find the maximum 0 values in order to create the rectangle
            shapevals = np.where(average_mask != 1)
                #am using 1 for x values and 0 for y values since that is what is done in the main code
            xmin = min(shapevals[1]) - 5
            xmax = max(shapevals[1]) + 5
            ymin = min(shapevals[0]) - 5
            ymax = max(shapevals[0]) + 5
            
            if xmin < 0:
                xmin = 5
            if ymin < 0:
                ymin = 5
            if xmax > 99:
                xmax = 95
            if ymax > 99:
                ymax = 95
            
            #print(average_mask[50])
            #print(xmin, xmax, ymin, ymax)
            cutout = np.zeros_like(average_mask)
            cutout[ymin:ymax, xmin:xmax] = average_mask[ymin:ymax, xmin:xmax]
            cutout[np.where(cutout == 0)] = float('nan')
                #this yields a cutout with nan inside our shape and outside the bounding rectangle
                #else it is 1 
            
                #this gives us the values inside rect but outside of the shape
            goodvals = np.where(np.isfinite(cutout))
            #print('GOOOOD\n', goodvals)
            x = goodvals[1]  # x values of finite coordinates
            y = goodvals[0]  # y values of finite coordinates
                
                #define fvals without making a function
                #fvals are then
            range_array = np.arange(x.size)
            fvals = np.zeros(x.size)
            for (i, xi, yi) in zip(range_array, x, y):
                fvals[i] = img[yi][xi] #get proper image here
            #print('fvals', fvals)
            #print(fvals)
            
            

                
            newfunc = interpolate.Rbf(
                x, y, fvals,
                function='multiquadric')  # the function that does interpolation
            allvals = np.where(img)  # whole region to interpolate over
            xnew = allvals[1]
            ynew = allvals[0]
            fnew = newfunc(xnew, ynew)
                    
            
            #Reduced Make_2D function
            #put the interpolated values back into a 2D array for display and other uses
            fnew_2D = np.zeros(
                (int(xnew.size /
                     ((img.shape)[0])), int(ynew.size / ((img.shape)[1]))),
                dtype=float)
            
            range_array = np.arange(fnew.size)
                    
            for (i, x, y) in zip(range_array, xnew, ynew):
                fnew_2D[y][x] = fnew[i]
                        
                
                
            resid = img - (img * region_mask + fnew_2D * inverse_mask) 
            residlist.append(resid)
            #interp = img * region_mask + fnew_2D * inverse_mask
            #averagelist.append(average_mask)
                        
                        
          
       
            
        flux = get_flux(residlist[3], residlist[2], residlist[1], residlist[0])
        fluxvalues = [flux.um70, flux.um24, flux.um12, flux.um8]
       
        
        fig, axs = plt.subplots(3, 4, figsize=(8, 16)) 
        plt.title(YBnum)
        
        for i, ax in enumerate(axs.flatten()):
            if i < 4:
                ax.imshow(averagelist[i], cmap='hot')
            elif i < 8:
                ax.imshow(residlist[i-4], cmap='hot')
            elif i < 12:
                ax.imshow(imagelist[i-12], cmap='hot', norm=SymLogNorm(linthresh= LinearThreshold))
    
    
            ax.axis('off')
    
        plt.tight_layout()
        plt.show()             
        #plt.pause(-1)
        
        
        return 0
    
    
    def photom_overlay(self, YBnum):
        
        
        #get the YB's location and radius
        print(YBnum, data['YB'][YBnum-1])
        YB_long = data[YBnum-1]['l']
        YB_lat = data[YBnum-1]['b']
        YB_rad = data[YBnum-1]['r']
    
        #Use the location to determine the correct image files to use
        #GWC added YB_lat on 4/5/22.
        image = choose_image(YB_long, YB_lat)
    
        ##          Unit Conversions Info         ##
        # 8 um:
        #GLIMPSE 8 um flux units are MJy/steradian
        #obtain the 8 um pixel scale from the image header
        #this gives the degree/pixel conversion factor used for overdrawing circle
        #and flux conversions (this is x direction only but they are equal in this case)
        #GWC edit 4/1/22: Set pixscale directly to 0.000333333
        #pixscale8 = abs(image.um8w.wcs.cd[0][0])
        pixscale8 = 0.000333333
        #pixscale8 = abs(image.um8w.wcs.cdelt1)
        #print(pixscale8)
        #8 um square degree to square pixel conversion-- x*y
        #GWC edit 4/1/22: Set sqdeg_tosqpix8 directly to pixscale8 * pixscale8
        #sqdeg_tosqpix8 = abs(image.um8w.wcs.cd[0][0]) * abs(
        #    image.um8w.wcs.cd[1][1])
        sqdeg_tosqpix8 = pixscale8 * pixscale8
        #8 um steradian to pixel conversion (used to convert MJy/Pixel)
        #     will be used to convert to MJy
        str_to_pix8 = sqdeg_tosqpix8 * 0.0003046118
        #WISE Units are Data Number per Pixel, Jy/DN is 1.8326*10^-6
        #See http://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec2_3f.html
        dn_to_jy12 = 1.83 * 10**-6
    
        #convert YB l and b and radius to pixel coordinates
        ybwcs = np.array([[YB_long, YB_lat]], np.float_)
        pixcoords = image.um8w.wcs_world2pix(ybwcs, 1)
        YB_long_pix = pixcoords[0][0]
        YB_lat_pix = pixcoords[0][1]
        YB_rad_pix = YB_rad / pixscale8
    
        # This is for the added points to show the user what other YBs are in the images
        # Read in the l, b, and r values for all the YBs and convert them to pixels
        YBloc = pd.read_csv(catalog_name, usecols = ['l', 'b', 'r'])
        # Convert l, b, and r from YBloc into pixels
        for i in range(0, len(YBloc)):
            yblocwcs = np.array([[YBloc['l'][i], YBloc['b'][i]]], np.float_)
            pixcoordsloc = image.um8w.wcs_world2pix(yblocwcs, 1)
            YB_l = pixcoordsloc[0][0]
            YBloc['l'][i] = YB_l
            YB_b = pixcoordsloc[0][1]
            YBloc['b'][i] = YB_b
            YB_radloc = YBloc['r'][i]
            YB_r = YB_radloc / pixscale8
            YBloc['r'][i] = YB_r
    
        #Create cropped 100 x 100 pixel image arrays centered on YB
        #orig = image.um8data[y1:y2,x1:x2]
        #orig12 = image.um12data[y1:y2,x1:x2]
        #orig24 = image.um24data[y1:y2,x1:x2]
    
        #use Cutout2D to make the zoomed windows
        position = (YB_long_pix + 0.5, YB_lat_pix + 0.5)
        size = (100, 100)
    
        cut8 = Cutout2D(data=image.um8data,
                        position=position,
                        size=size,
                        wcs=image.um8w)
        cut12 = Cutout2D(data=image.um12data,
                         position=position,
                         size=size,
                         wcs=image.um12w)
        cut24 = Cutout2D(data=image.um24data,
                         position=position,
                         size=size,
                         wcs=image.um24w)
        cut70 = Cutout2D(data=image.um70data,
                         position=position,
                         size=size,
                         wcs=image.um70w)
    
        fitcopy8 = image.um8
        fitcopy8.data = cut8.data
        fitcopy8.header.update(cut8.wcs.to_header())
    
        fitcopy12 = image.um12
        fitcopy12.data = cut12.data
        fitcopy12.header.update(cut12.wcs.to_header())
    
        fitcopy24 = image.um24
        fitcopy24.data = cut24.data
        fitcopy24.header.update(cut24.wcs.to_header())
    
        fitcopy70 = image.um70
        fitcopy70.data = cut70.data
        fitcopy70.header.update(cut70.wcs.to_header())
    
        orig = cut8.data
        orig12 = cut12.data
        orig24 = cut24.data
        orig70 = cut70.data
    
        #create copies of cropped images called workmasks
        workmask = copy.deepcopy(orig)
        workmask12 = copy.deepcopy(orig12)
        workmask24 = copy.deepcopy(orig24)
        workmask70 = copy.deepcopy(orig70)
        
        LinearThreshold = 0.001
        
        ###################################
        #
        ## Here is where we start to actually differ from the main code
        #
        ###################################
        
        #fig = plt.figure(figsize=(8, 8))
        
        ####################################################################################################################################################
        #
        #This is the list of the names of all the csv's of points that we are using
        #
        #Update this list with the csv files you want to use for the average mask. 
        #
        ####################################################################################################################################################

        csvlist = ['YBphotometry_results_EthanBass.csv', 'YBphotometry_results_Ignasi.csv', 'YBphotometry_results_colemaya.csv', 'YBphotometry_results_IgnasiBonmatiGonzalvez.csv', 'YBphotometry_results_KORRA.csv', 'YBphotometry_results_WolfChase1-ALL_YBs.csv', 'YBphotometry_results_Asta.csv']
        #csvlist = ['YBphotometry_results_EthanBass.csv']
    
    
        #csvlist = ['YBphotometry_results_EthanBass.csv']
        
        headerlist = ['vertices 70', 'vertices 24', 'vertices 12', 'vertices 8']
        imagelist = [workmask70, workmask24, workmask12, workmask]
        averagelist = []
        residlist = []
        baselist = []
        errorlist = []

        for i in range(0, 4):
            average_mask = np.zeros_like(workmask70)
            img = imagelist[i]
            header = headerlist[i]
            validlength = 0
            for csvname in csvlist:
                userdata = ascii.read(csvname, delimiter=',')
                if userdata[YBnum][header] != '':
                    vertices = ast.literal_eval(userdata[YBnum][header])
                    interp = do_interp(img, vertices)
                    if csvname == 'YBphotometry_results_EthanBass.csv':
                        baselist.append(interp.resid)
                    im = interp.masked
                    im[np.where(im != 1)] = int(0)
                    average_mask += im
                    validlength += 1
                #else:
                    #print(csvname)
                    
                
            print(f"Num of sources: {validlength}, for YB {YBnum}, at {headerlist[i]}um")
            if (validlength == 0):
                return 0
            average_mask /= validlength
            #print(average_mask[50])
            averagelist.append(np.copy(average_mask)*-1 +1)
    
                #here average_mask is 1 outside, 0 inside of all shapes, and some fractions between
            sensitivity = 0.5 #determines the breakoff point
                # a low sensitivity means that we take the smaller shape
                # a high sensitivity means that we take bigger shapes
                # thus a low sensitivity will be more 'zoomed in' on the YB

            errorlist.append((average_mask[np.where(average_mask <  sensitivity)]).sum().sum() + (1 - average_mask[np.where(average_mask >= sensitivity)]).sum().sum())

            average_mask[np.where(average_mask >= sensitivity)] = 1
            average_mask[np.where(average_mask <  sensitivity)] = float('nan')
            #print(average_mask[50])
            region_mask = np.copy(average_mask)
            region_mask[np.where(region_mask != 1)] = 0
                
            inverse_mask = region_mask*-1 +1
                
                #In order to recreate goodvals, I need to find the maximum 0 values in order to create the rectangle
            shapevals = np.where(average_mask != 1)
                #am using 1 for x values and 0 for y values since that is what is done in the main code
            xmin = min(shapevals[1]) - 5
            xmax = max(shapevals[1]) + 5
            ymin = min(shapevals[0]) - 5
            ymax = max(shapevals[0]) + 5
            
            if xmin < 0:
                xmin = 5
            if ymin < 0:
                ymin = 5
            if xmax > 99:
                xmax = 95
            if ymax > 99:
                ymax = 95
            '''if average_mask.sum() == 0:
                print('this line works')'''
            #print(average_mask[50])
            #print(xmin, xmax, ymin, ymax)
            cutout = np.zeros_like(average_mask)
            cutout[ymin:ymax, xmin:xmax] = average_mask[ymin:ymax, xmin:xmax]
            cutout[np.where(cutout == 0)] = float('nan')
                #this yields a cutout with nan inside our shape and outside the bounding rectangle
                #else it is 1 
            
                #this gives us the values inside rect but outside of the shape
            goodvals = np.where(np.isfinite(cutout))
            #print('GOOOOD\n', goodvals)
            x = goodvals[1]  # x values of finite coordinates
            y = goodvals[0]  # y values of finite coordinates
                
                #define fvals without making a function
                #fvals are then
            range_array = np.arange(x.size)
            fvals = np.zeros(x.size)
            for (i, xi, yi) in zip(range_array, x, y):
                fvals[i] = img[yi][xi] #get proper image here
            #print('fvals', fvals)
            #print(fvals)
            
            

                
            newfunc = interpolate.Rbf(
                x, y, fvals,
                function='multiquadric')  # the function that does interpolation
            allvals = np.where(img)  # whole region to interpolate over
            xnew = allvals[1]
            ynew = allvals[0]
            fnew = newfunc(xnew, ynew)
                    
            
            #Reduced Make_2D function
            #put the interpolated values back into a 2D array for display and other uses
            fnew_2D = np.zeros(
                (int(xnew.size /
                     ((img.shape)[0])), int(ynew.size / ((img.shape)[1]))),
                dtype=float)
            
            range_array = np.arange(fnew.size)
                    
            for (i, x, y) in zip(range_array, xnew, ynew):
                fnew_2D[y][x] = fnew[i]
                        
                
                
            resid = img - (img * region_mask + fnew_2D * inverse_mask) 
            residlist.append(resid)
            #interp = img * region_mask + fnew_2D * inverse_mask
            #averagelist.append(average_mask)
                        
                        
          
       
            
        flux = get_flux(residlist[3], residlist[2], residlist[1], residlist[0])
        fluxvalues = [flux.um70, flux.um24, flux.um12, flux.um8]
       
        
        '''fig, axs = plt.subplots(4, 4, figsize=(8, 16)) 
        plt.title(YBnum)
        
        for i, ax in enumerate(axs.flatten()):
            if i < 4:
                ax.imshow(averagelist[i], cmap='hot')
            elif i < 8:
                ax.imshow(residlist[i-4], cmap='hot')
            elif i < 12:
                ax.imshow(imagelist[i-12], cmap='hot', norm=SymLogNorm(linthresh= LinearThreshold))
    
    
            ax.axis('off')
    
        plt.tight_layout()
        plt.show()   '''            
        #plt.pause(-1)
        self.results.loc[len(self.results)] = [int(YBnum)]+fluxvalues+errorlist
        
        
        return 0
    
    
#TODO change photom_overlay to be an object and have results be an item of the object. 
Mask = AverageMask()

for i in range(1672, 2009):
    Mask.photom_overlay(i)

fluxvals = userdata = ascii.read("YB_frequency.csv", delimiter=',')

graphvals = pd.merge(fluxvals.to_pandas(), Mask.results, on='YB')

'''for i in [70, 24, 12, 8]: 
    plt.figure(figsize=(8, 6))  # Adjust the figure size if needed
    plt.plot(Mask.results['YB'], Mask.results[f'AvMaskFlux{i}'], marker='o', linestyle='None')
    plt.plot(fluxvals["YB"], fluxvals[f"flux{i}_median"], marker='o', linestyle='None')
    plt.plot(fluxvals["YB"], fluxvals[f"flux{i}_mean"], marker='o', linestyle='None')
    plt.xlim(Mask.results["YB"].min()-0.5, Mask.results["YB"].max()+0.5)
    plt.title(f"Flux values at {i}um")
    plt.yscale('log')
    plt.xticks(Mask.results['YB'])
    plt.show()'''

from scipy.stats import pearsonr
from scipy import stats
Mask = AverageMask() 
graphvals = pd.read_csv('AverageMaskAllData.csv')

outlierslist = []
threshold = 5
for x in ['']:
    for i in [70, 24, 12, 8]:
        
        posvals = graphvals[graphvals[f"flux{i}_mean"]>0]
        
        # Calculate Z-scores
        z_scores = stats.zscore(posvals[f'AvMaskEr{i}{x}'])

        # Find outliers
        outliers = (np.abs(z_scores) > threshold)

        # Print outliers
        print("Outliers:", posvals['YB'][outliers])
        for i in posvals['YB'][outliers]:
            if i not in outlierslist:
                outlierslist.append(i)
    
for i in [70, 24, 12, 8]:
    
    posvals = graphvals[graphvals[f"flux{i}_mean"]>0]
    
    # Calculate Z-scores
    z_scores = stats.zscore(posvals[f'flux{i}_std']/posvals[f'flux{i}_mean'])

    # Find outliers
    outliers = (np.abs(z_scores) > threshold)

    # Print outliers
    print("Outliers:", posvals['YB'][outliers])
    for i in posvals['YB'][outliers]:
        if i not in outlierslist:
            outlierslist.append(i)

for i in outlierslist:
    Mask.show_mask(i)


for i in [70, 24, 12, 8]:
    posvals = graphvals[graphvals[f"flux{i}_mean"]>0]
    
    plt.figure(figsize=(8, 6)) 
    #plt.plot(graphvals[f'AvMaskFlux{i}'], graphvals[f"flux{i}_median"], marker='o', linestyle='None')
    plt.plot(posvals[f'AvMaskFlux{i}'], posvals[f"flux{i}_mean"], marker='o', linestyle='None')
    plt.plot(posvals[f'AvMaskFlux{i}'][posvals['YB'].isin(outlierslist)], posvals[f"flux{i}_mean"][posvals['YB'].isin(outlierslist)], marker='o', linestyle='None')
    plt.title(f"Average Mask Flux vs Mean Flux at {i}μm")
    plt.axline((0, 0), slope=1, color='black', linestyle='--')
    plt.yscale('log')
    plt.xscale('log')
    plt.show()

for x in ['']: 
    for i in [70, 24, 12, 8]:
        
        posvals = graphvals[graphvals[f"flux{i}_mean"]>0]
        
        for j in outlierslist:
            posvals = posvals[posvals["YB"] != j]
        
        #posvals = posvals[posvals['YB'] != 1693]
        corr_coeff, p_value = pearsonr(posvals[f"AvMaskEr{i}{x}"], posvals[f"flux{i}_std"]/posvals[f"flux{i}_mean"])
        
        # Calculate standard error of the estimate
        
        std_err_estimate = np.sqrt((1 - corr_coeff**2) * np.var(graphvals[f"flux{i}_std"]/graphvals[f"flux{i}_mean"]))
        
        # Output the results
        
        print(f"\nPearson correlation coefficient: {corr_coeff:.4f} {x}")
        print(f"P-value: {p_value:.4f} {x}")
        print(f"Standard error of the estimate: {std_err_estimate:.4f} {x}\n")
        
        plt.figure(figsize=(8, 6)) 
        plt.plot(posvals[f"AvMaskEr{i}{x}"], posvals[f"flux{i}_std"]/posvals[f"flux{i}_mean"], marker='o', linestyle='None')
        plt.title(f"Average Mask Error indicator vs standard deviation at {i}μm {x}")
        plt.show()
 
    
 
graphvals.to_csv('AverageMaskAllData.csv', index=True)  
 
    
 
    



