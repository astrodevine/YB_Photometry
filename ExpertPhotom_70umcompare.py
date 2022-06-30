'''
Version Last updated 7/9/2021 by Makenzie Stapley
-Includes 70 micron images
-Skips images that are saturated and writes out ' ' for those coordinate values
-Axis labels removed to make the 6 and 4-panel images pretty
-Compactness class and calculations have been removed
-If there are multiple YBs in the cropped image, the rg image has a blue circle with
the MWP location and radius
'''

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
#pip install opencv-python
import os
import pandas as pd
from astropy.nddata import Cutout2D

# These lines supress warnings
import warnings
warnings.filterwarnings('ignore')

#For beta testing "Pick up where you left off?" option
import pickle

#You will need to make sure the following packages have been installed first:

#from tkinter import *
#import tkinter
#conda install -c anaconda tk

#from photutils import centroid_com
#https://photutils.readthedocs.io/en/stable/install.html
#conda install photutils -c astropy

from itertools import chain, repeat
#import threading
#import mynormalize
#import mycolorbar
""" Katie Devine/Anupa Poudyal

Opens images of YBs and prompts users to click around YB images at 8, 12, and 24 um.
Polygon mask is then applied, and background interpolation done. 
Residual (original-background) is used to calculate total flux density in Jy 
Output files are a table of the photometric results and images showing the mask, interpolation, and residual.

To use, make sure you have your file system set up according to the file i/o section
including the output folder for the images, change the settings and flags to whatever you
would like them to be, enter the range of sources you want to analyze, and run from the 
command line.

This code requires that the 8, 12, and 24 um mosaics all have the same dimensions
(i.e. that the 12 and 24 um mosaics have been reprojected to match the 8 um pixel size)
Conventions:
Mips in red, Glimpse in green
Longer dashes, longer wavelengths (8um gets dots, 24um gets dashes)
fpa is flux per area

Files and directory structure needed: 
    -all of the mosaic files listed below in path+/mosaics
    -directory called path+/photom_images/ to store output images
    -"YBphotometry_results.csv" will be created when program is run for the first time and appended each time program is run
#"""

######################################################
#Paths and File Locations                            #
######################################################

#EDIT THIS PATH FOR THE FILE LOCATION ON YOUR MACHINE
# '.' means current directory
path = '.'
path1 = '.'
image_name = os.path.join(path, 'GLM_03000+0000_mosaic_I4.fits')
catalog_name = os.path.join(path, 'USE_THIS_CATALOG_ybcat_MWP_with_ID.csv')
#out_name = os.path.join(path1, 'YBphotometry_results.csv')
instID = 'Devine_70umcompare' #Change to be your ID
out_name = os.path.join(path, 'YB_70um_photometry_results_' + instID + '.csv')

######################################################
# Define my functions and classes                    #
######################################################


#function that shows image and collects the coordinates for interpolation
def get_coords(img, imgw, wave, ybid):
    #Plot the cropped image, this will be for user to draw polygon on
    global clkfig
    #clkfig = plt.figure()
    #clkfig.add_subplot(111, projection = imgw)
    #clkfig = plt.subplot(1,2,2, title = 'select the coordinates for polygon in this image.', projection = imgw)
    #clkfig.suptitle('select the coordinates for polygon in this image.', fontsize=8)
    #plt.imshow(img, cmap = 'hot', norm = LogNorm())
    #plt.xlabel('Longitude')
    #plt.ylabel('Latitude')

    circle = Circle((dxw / 2, dyw / 2), YB_rad_pix, fill=False)
    circle1 = Circle((dxw / 2, dyw / 2), YB_rad_pix, fill=False)
    circle2 = Circle((dxw / 2, dyw / 2), YB_rad_pix, fill=False)
    circle3 = Circle((dxw / 2, dyw / 2), YB_rad_pix, fill=False)
    circle4 = Circle((dxw / 2, dyw / 2), YB_rad_pix, fill=False)
    circle5 = Circle((dxw / 2, dyw / 2), YB_rad_pix, fill=False)

    clkfig = plt.figure(figsize=(18, 27))

    # Create a df of plottable points from YBloc
    plotloc = []
    for i in range(0, len(YBloc)):
        if 0 <  abs(YBloc['l'][i] - YB_long_pix) < (dxw/2 - 5) and 0 < abs(YBloc['b'][i] - YB_lat_pix) < (dyw/2 -5):
            plotloc.append((YBloc['l'][i] - YB_long_pix + (dxw/2 - 5), YBloc['b'][i] - YB_lat_pix + (dyw/2 -5)))

    #Plot rg image,
    axrg = plt.subplot(3, 3, 1, title='RG (24 + 8 um)', projection=imgw)
    r = orig24
    g = orig
    b = np.zeros(orig.shape)
    axrg.axis('off')
    axrg.imshow(make_lupton_rgb(r, g, b, stretch=200, Q=0), norm=LogNorm())
    axrg.add_artist(circle)
    for point in range(0, len(plotloc)):
        #axrg.plot(plotloc[point][0], plotloc[point][1], marker='+', markersize=40)
        circleloc = Circle((plotloc[point][0], plotloc[point][1]), YBloc['r'][point], fill = False, color = 'Blue')
        axrg.add_artist(circleloc)

    #Add subtitles to the plot
    clkfig.add_subplot(3,
                       3,
                       2,
                       title='Select coordinates in this ' + wave + ' image',
                       projection=imgw)
    clkfig.suptitle('You are examining the ' + wave + ' image for YB' +
                    str(ybid))
    plt.axis('off')
    plt.imshow(img, cmap='gray', norm=LogNorm())

    #Plot the 70 um
    if math.isnan(orig70m.min()) == False and math.isnan(orig70m.max()) == False:
        fee = plt.subplot(3, 3, 3, title='70 um MIPS', projection=imgw)
        plt.axis('off')
        plt.imshow(orig70m,
                   cmap='hot',
                   norm=LogNorm(vmin=orig70m.min(), vmax=orig70m.max()))
        fee.add_artist(circle1)
    else:
        fee = plt.subplot(3, 3, 3, title='70 um MIPS', projection=imgw)
        plt.axis('off')
        plt.imshow(orig70,
                   cmap='hot',
                   norm=LogNorm())
        fee.add_artist(circle1)
        

    #Plot the 24 um
    if math.isnan(orig24.min()) == False and math.isnan(orig24.max()) == False:
        foo = plt.subplot(3, 3, 4, title='24 um', projection=imgw)
        plt.axis('off')
        plt.imshow(orig24,
                   cmap='hot',
                   norm=LogNorm(vmin=orig24.min(), vmax=orig24.max()))
        foo.add_artist(circle2)
    else:
        foo = plt.subplot(3, 3, 4, title='24 um', projection=imgw)
        plt.axis('off')
        plt.imshow(orig24,
                   cmap='hot',
                   norm=LogNorm())
        foo.add_artist(circle2)
        
    #Plot the 12 um
    if math.isnan(orig12.min()) == False and math.isnan(orig12.max()) == False:
        faa = plt.subplot(3, 3, 5, title='12 um', projection=imgw)
        plt.axis('off')
        plt.imshow(orig12,
                   cmap='hot',
                   norm=LogNorm(vmin=orig12.min(), vmax=orig12.max()))
        faa.add_artist(circle3)
    else:
        faa = plt.subplot(3, 3, 5, title='12 um', projection=imgw)
        plt.axis('off')
        plt.imshow(orig12,
                   cmap='hot',
                   norm=LogNorm())
        faa.add_artist(circle3)

    #Plot the 8um
    if math.isnan(orig.min()) == False and math.isnan(orig.max()) == False:
        fum = plt.subplot(3, 3, 6, title='8 um', projection=imgw)
        plt.axis('off')
        plt.imshow(orig,
                   cmap='hot',
                   norm=LogNorm(vmin=orig.min(), vmax=orig.max()))
        fum.add_artist(circle4)
    else:
        fum = plt.subplot(3, 3, 6, title='8 um', projection=imgw)
        plt.axis('off')
        plt.imshow(orig,
                   cmap='hot',
                   norm=LogNorm())
        fum.add_artist(circle4)
        
        
    #Plot the PACS 70 um
    if math.isnan(orig70.min()) == False and math.isnan(orig70.max()) == False:
        fum = plt.subplot(3, 3, 7, title='70 um PACS', projection=imgw)
        plt.axis('off')
        plt.imshow(orig70,
                   cmap='hot',
                   norm=LogNorm(vmin=orig70.min(), vmax=orig70.max()))
        fum.add_artist(circle5)
    else:
        fum = plt.subplot(3, 3, 7, title='70 um PACS', projection=imgw)
        plt.axis('off')
        plt.imshow(orig70,
                   cmap='hot',
                   norm=LogNorm())
        fum.add_artist(circle5)
    #cbar = plt.colorbar(format='%05.2f')
    #cbar.set_norm(mynormalize.MyNormalize(vmin=orig.min(),vmax=orig.max(),stretch='linear'))
    #cbar = mycolorbar.DraggableColorbar(cbar,orig)
    #cbar.connect()

    #Plot instructions
    plt.subplot(3, 3, 9)
    plt.axis([0, 10, 0, 10])
    plt.axis('off')
    text = ("*Click in the grey-scale image. ")
    text1 = ("*Left click to add points.")
    text2 = ("*Middle click when finished to exit.")
    text3 = ("*All other images are for reference.")
    text4 = ("*Circles indicate MWP user radius.")
    text5 = ("*You will be prompted to inspect results;")
    text6 = ("type 'n' to continue, anything else to redo.")
    plt.text(1, 8, text, ha='left', wrap=True)
    plt.text(1, 7, text1, ha='left', wrap=True)
    plt.text(1, 6, text2, ha='left', wrap=True)
    plt.text(1, 5, text3, ha='left', wrap=True)
    plt.text(1, 4, text4, ha='left', wrap=True)
    plt.text(1, 3, text5, ha='left', wrap=True)
    plt.text(1, 2, text6, ha='left', wrap=True)

    #interactive clicking to fill up coords
    coords = clkfig.ginput(n=-1, timeout=300, show_clicks=True, mouse_stop=2)

    plt.close('all')
    #print(coords)
    return coords


#generates and saves the four panel images to the folder photom_images
#call with (image, masked image, interpolated image, resid)
def make_figs(im1, im2, im3, im4, fitfile, imw, um):

    ############Generate the figures for each source##################
    #note I'm being lazy here and calling from the code things that aren't defined in function
    #this is pretty ugly and should maybe get cleaned up
    fig = plt.figure(figsize=(8, 8))

    #Plot the original image and MWP User YB circle
    circle = Circle((dxw / 2, dyw / 2), YB_rad_pix, fill=False)
    fig1 = plt.subplot(2, 2, 1, title='Cropped image', projection=imw)
    plt.imshow(im1, cmap='hot', norm=LogNorm(vmin=im1.min(), vmax=im1.max()))
    plt.axis('off')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    fig1.add_artist(circle)
    fig1.text(dxw / 2,
              dyw / 2 - 5,
              'MWP size',
              color='white',
              ha='center',
              va='top',
              weight='bold')

    #Plot the mask
    plt.subplot(2, 2, 2, title='Masked Image', projection=imw)
    plt.imshow(im2, cmap='hot', norm=LogNorm(vmin=im1.min(), vmax=im1.max()))
    plt.axis('off')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')

    #Plot the interpolated background
    plt.subplot(2, 2, 3, title='Interpolated image', projection=imw)
    plt.imshow(im3, cmap='hot', norm=LogNorm(vmin=im1.min(), vmax=im1.max()))
    plt.axis('off')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')

    #Plot the residual
    plt.subplot(2, 2, 4, title='Residual(Image-Interp)', projection=imw)
    plt.imshow(im4, cmap='hot')
    plt.axis('off')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')

    #Make the plots pretty
    plt.tight_layout(pad=0.2, w_pad=5, h_pad=6)
    plt.subplots_adjust(top=0.9)
    fig.suptitle('Interpolation of YB %s at ' % (YB) + um, fontsize=10)

    plt.pause(1)
    
    # Save this result as a new png
    figurename = os.path.join(path1,
                              f"photom_images/{um}interpolation_YB_{YB}.png")
    fig.savefig(figurename)

    # Save the fits cut-outs for future use if needed
    im1name = os.path.join(path1, f"fits_cutouts/{um}cropped_YB_{YB}.fits")
    im2name = os.path.join(path1, f"fits_cutouts/{um}masked_YB_{YB}.fits")
    im3name = os.path.join(path1, f"fits_cutouts/{um}interp_YB_{YB}.fits")
    im4name = os.path.join(path1, f"fits_cutouts/{um}resid_YB_{YB}.fits")

    fitfile.data=im1
    fitfile.writeto(im1name, overwrite=True)

    fitfile.data=im2
    fitfile.writeto(im2name, overwrite=True)

    fitfile.data=im3
    fitfile.writeto(im3name, overwrite=True)

    fitfile.data=im4
    fitfile.writeto(im4name, overwrite=True)

######################################################
# Define my classes                                  #
######################################################
#class choose_image: use l range to choose and open right mosaics
# returns .um8, .um8data, .um8w - full info, data array, wcs header at 8 um
# returns .um12, .um12data, .um12w - full info, data array, wcs header at 12 um
# returns .um24, .um24data, .um124w - full info, data array, wcs header at 24 um
# returns .um70, .um70data, .um70w  (from Herschel PACS)
# returns .um70M, .um70dataM, .um70wM (same as above but for MIPS data)
class choose_image():
    def __init__(self, l, b):
        #currently don't need the WCS files for 12, 24 um because they
        #are reprojections onto 8um coordinate grid
        #GWC added 'b' on 4/5/22. 

        if l >= 28.5 and l <= 31.5:
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
            path70m = os.path.join(
                path,
                'mosaics/SMOG_70um_10300_mosaic.fits')
        elif l > 101.0 and l <= 105.59 and b >= 3.06:
            path8 = os.path.join(path,
                                 'mosaics/SMOG_08um_10300_mosaic.fits')
            path12 = os.path.join(path,
                                  'mosaics/SMOG_12um_10300_mosaic.fits')
            path24 = os.path.join(
                path,
                'mosaics/SMOG_24um_10300_mosaic_high_b.fits')
            path70m = os.path.join(
                path,
                'mosaics/SMOG_70um_10300_mosaic.fits')    
            path70 = os.path.join(
                path, 'mosaics/SMOG_PACS_70um_10300_mosaic.fits')
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
            path70m = os.path.join(
                path,
                'mosaics/SMOG_70um_10700_mosaic.fits')
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
        temp = fits.open(path70m)[0]
        self.um70m = temp
        self.um70datam = temp.data
        self.um70wm = wcs.WCS(temp.header)


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

        for i in x:
            for j in y:

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


#class that gets the flux from residual images. Unit conversions are applied here.
class get_flux():
    def __init__(self, d70, d70m):
        #reset the total flux value
        
        flux_tot70 = 0
        flux_tot70m = 0
        
        #go through 100 x 100 pixel residual image and add up the pixel values
        for ROW in range(0, 100):
            for column in range(0, 100):

                flux_tot70 = flux_tot70 + d70[ROW][
                    column]  #gives the value of photometry
                flux_tot70m = flux_tot70m + d70m[ROW][
                    column]  #gives the value of photometry
                
            #convert units of flux total.  MIPS/IRAC in MJy/Sr*Pixel, want Jy
            #conversion: (MJy/Sr*Pixel)*(Sr/Pix)*(Jy/MJy)
            #WISE is in units of data number per pixel
            #Conversion: DN/Pix*Pixel*Jy/DN


        flux_tot70 = flux_tot70 #No conversion necesary for Herschel data
        flux_tot70m = flux_tot70m * str_to_pix8 * 10**6


        self.um70 = flux_tot70
        self.um70m = flux_tot70m


######################################################
# Actual Program Begins Here                         #
######################################################

#Open the catalog file to get YB names, l, b, radius
data = ascii.read(catalog_name, delimiter=',')

#Open the output file and write the column headers

if os.path.exists(out_name):
    append_write = 'a'  # append if already exists
    output_file = open(out_name,
                       append_write)  #####opens up files for creating csv file
    headers = [
        'YB', 'YB_long', 'YB_lat', 'vertices 70', '70umphotom MIPS', '70umphotom PACS']
    writer = csv.DictWriter(output_file, fieldnames=headers)
    output_file.close()

else:
    output_file = open(out_name,
                       'w')  #####opens up files for creating csv file
    writer = csv.DictWriter(
        output_file,
        fieldnames=[
        'YB', 'YB_long', 'YB_lat', 'vertices 70', '70umphotom MIPS', '70umphotom PACS'],
        lineterminator='\n')
    writer.writeheader()
    writer.writerow({
        'YB': 'ID Number',
        'YB_long': 'degree',
        'YB_lat': 'degree',
        'vertices 70': 'pixel coords',
        '70umphotom MIPS': 'Jy',
        '70umphotom PACS': 'Jy',
    })

    for k in range(0,6176): ## THIS IS THE FULL CATALOG RANGE
        YB = data[k]['YB']
        YB_long = data[k]['l']
        YB_lat = data[k]['b']
        #output_file = open(out_name,append_write) #####opens up file for writing data
        #writer=csv.DictWriter(output_file,fieldnames=['YB','YB_long','YB_lat','8umphotom','12umphotom','24umphotom'])

        writer.writerow({'YB': YB, 'YB_long': YB_long, 'YB_lat': YB_lat})

    output_file.close()

######################################################
# Begin the loop through YBs in the catalog (k=YB)   #
######################################################

#Set Pre-chosen range
BegYB = 3074
YBlast = 3100

#The below allows the user to start at the beginning of the range or from where they left off last time the program was ran
print('Welcome! Where would you like to start?')
startpos = input(
    "Enter 'a' to start from the beginning or anything else to start where you left off! "
)
if startpos == 'a':
    YBfirst = BegYB
else:
    pickle_in = open("leftYB.pickle", "rb")
    leftYB = pickle.load(
        pickle_in)  #loading in the last YB they completed last time
    YBfirst = leftYB + 1

YB1 = int(YBfirst) - 1
YB2 = int(YBlast)
k = YB1
currentYB = YB1
while (k < YB2):
    #get the YB's location and radius
    YB = data[k]['YB']
    YB_long = data[k]['l']
    YB_lat = data[k]['b']
    YB_rad = data[k]['r']

    #Use the location to determine the correct image files to use
    image = choose_image(YB_long, YB_lat)

    ##          Unit Conversions Info         ##
    # 8 um:
    #GLIMPSE 8 um flux units are MJy/steradian
    #obtain the 8 um pixel scale from the image header
    #this gives the degree/pixel conversion factor used for overdrawing circle
    #and flux conversions (this is x direction only but they are equal in this case)
    pixscale8 = abs(image.um8w.wcs.cdelt[0])
    #8 um square degree to square pixel conversion-- x*y
    #MODIFIED FOR SMOG DATA HEADERS!
    sqdeg_tosqpix8 = abs(image.um8w.wcs.cdelt[0]) * abs(
        image.um8w.wcs.cdelt[1])
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
    
    #define a window to zoom in on the YB
    xw = xs = int(YB_long_pix +
                  0.5)  # x coordinates of source and window center
    yw = ys = int(YB_lat_pix +
                  0.5)  # y coordinates of source and window center
    dxw = 200  # x width of window
    dyw = 200  # y width of window

    #find the pixel coordinates LLH and URH corner of zoomed window
    x1 = int(xw - 0.5 * dxw)
    y1 = int(yw - 0.5 * dyw)
    x2 = int(xw + 0.5 * dxw)
    y2 = int(yw + 0.5 * dyw)
        
    #use Cutout2D to make the zoomed windows
    position = (YB_long_pix + 0.5, YB_lat_pix + 0.5)
    size = (dxw, dyw)

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
    cut70m = Cutout2D(data=image.um70datam,
                     position=position,
                     size=size,
                     wcs=image.um70wm)
    
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

    fitcopy70m = image.um70m
    fitcopy70m.data = cut70m.data
    fitcopy70m.header.update(cut70m.wcs.to_header())


    orig = cut8.data
    orig12 = cut12.data
    orig24 = cut24.data
    orig70 = cut70.data
    orig70m = cut70m.data
    
    wcs8 = cut8.wcs
    wcs12 = cut12.wcs
    wcs24 = cut24.wcs
    wcs70 = cut70.wcs
    wcs70m = cut70m.wcs
    
    #create empty residuals to fill up later
    diff8 = orig * 0
    diff12 = orig12 * 0
    diff24 = orig24 * 0
    diff70 = orig70 * 0
    diff70m = orig70m * 0

    #create copies of cropped images called workmasks
    workmask = copy.deepcopy(orig)
    workmask12 = copy.deepcopy(orig12)
    workmask24 = copy.deepcopy(orig24)
    workmask70 = copy.deepcopy(orig70)
    workmask70m = copy.deepcopy(orig70m)

    ###################################################################
    # Call the classes to draw polygons and perform interpolation, use MIPS for polygon and apply to PACS    #
    ###################################################################
    if np.isnan(orig70m.min()) == False and np.isnan(orig70m.max()) == False:
        check = 'y'
        while check != 'n':
            # 70 um image analysis
            #reset global list coords that gets creaeted in get_coords
            print('######################################################')
            print('Beginning the 70 um analysis -- MIPS data')
            coords = []
            #get the coordinates on 70um image
            coordinates = get_coords(workmask70m, wcs70m, '70 um MIPS', YB)

            #if no clicks, exit
            if coordinates == []:
                print(
                    "Timeout limit of 5 minutes reached. Exiting program. Results for most recent YB will not be saved. Please pick 'start where you left off' option next time to re-start this YB."
                )
                pickle_out = open("leftYB.pickle", "wb")
                pickle.dump(currentYB, pickle_out)
                pickle_out.close()
                sys.exit()

            print('got coords')
            #do the masking and interpolation on 70um image
            print('starting interp')
            interp70m = do_interp(workmask70m, coordinates)
            diff70m = interp70m.resid
            #display and save the images for the 70um image
            make_figs(workmask70m, interp70m.blanked, interp70m.interp,
                      interp70m.resid, fitcopy70m, wcs70m, '70_um_MIPS')

            #Prompt user to accept image or redo
            plt.pause(1)
            check = input(
                'Please consult the residual image. Would you like to redo? Type n to continue, anything else to redo:  '
            )
            if check != 'n':
                plt.close('all')
        plt.close('all')
        coord70 = str(coordinates)
        plt.close('all')
    else:
        print('70 MIPS micron image is saturated.')
        coord70 = ' '
        
    #comparing requires no saturation in either image
    if np.isnan(orig70m.min()) == False and np.isnan(orig70m.max()) == False and np.isnan(orig70.min()) == False and np.isnan(orig70.max()) == False:
        check = 'y'
        while check != 'n':
            # 70 um image analysis
            #reset global list coords that gets creaeted in get_coords
            print('######################################################')
            print('Beginning the 70 um analysis -- PACS image')
            
            #commenting out all the getting coordinates-- only using the MIPS coordinates for this
            #coords = []
            #coordinates = get_coords(workmask70, wcs70, '70 um', YB)

            #if no clicks, exit
            #if coordinates == []:
            #    print(
            #        "Timeout limit of 5 minutes reached. Exiting program. Results for most recent YB will not be saved. Please pick 'start where you left off' option next time to re-start this YB."
            #    )
            #    pickle_out = open("leftYB.pickle", "wb")
            #    pickle.dump(currentYB, pickle_out)
            #    pickle_out.close()
            #    sys.exit()

            #print('got coords')
            print('using coordinates from MIPS image')
            #do the masking and interpolation on 70um image
            print('starting interp')
            interp70 = do_interp(workmask70, coordinates)
            diff70 = interp70.resid
            #display and save the images for the 70um image
            make_figs(workmask70, interp70.blanked, interp70.interp,
                      interp70.resid, fitcopy70, wcs70, '70_um')

            #Prompt user to accept image or redo
            plt.pause(1)
            check = input(
                'Please consult the residual image. Would you like to redo? Type n to continue, anything else to redo:  '
            )
            if check != 'n':
                plt.close('all')
        plt.close('all')
    else:
        print('PACS 70 micron image is saturated.')

    ##############################################################################
    # Use residual images to perform photometry and write out to table with flags#
    ##############################################################################

    #call the get_flux class
    flux_tot = get_flux(diff70,diff70m)

    #'8umphotom':flux_tot.um8,'12umphotom':flux_tot.um12,'24umphotom':flux_tot.um24}

    df = pd.read_csv(out_name)

    kk = str(YB)

    df.loc[df["YB"] == kk, '70umphotom PACS'] = round(flux_tot.um70, 5)
    df.loc[df["YB"] == kk, '70umphotom MIPS'] = round(flux_tot.um70m, 5)

    df.loc[df["YB"] == kk, 'vertices 70'] = coord70

    df.to_csv(out_name, index=False)
    k = k + 1
    currentYB = currentYB + 1
    pickle_out = open("leftYB.pickle", "wb")
    pickle.dump(currentYB, pickle_out)
    pickle_out.close()
    if k < YB2:
        #Allow the user to safely exit the program between YBs
        cont = input(
            "Would you like to continue? Type 'y' for yes or anything else for no: "
        )
        if cont == 'y':
            print('Okay! Continuing photometry...')
        else:
            print('Goodbye! See you next time!')
            #Save current YB, so the user can pick up from here next time they run the program
            pickle_out = open("leftYB.pickle", "wb")
            pickle.dump(currentYB, pickle_out)
            pickle_out.close()
            sys.exit()

plt.close('all')
print(
    'Congratulations! You have completed the photometry for all the YBs in your range!'
)
