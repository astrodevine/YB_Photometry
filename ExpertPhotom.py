'''
Version Last updated 7/15/2021 by Makenzie Stapley
-Includes 70 micron images
-Skips images that are saturated and writes out ' ' for those coordinate values
-Axis labels removed to make the 6 and 4-panel images pretty
-Compactness class and calculations have been removed
-If there are multiple YBs in the cropped image, the rg image has a blue circle with
the MWP location and radius
'''

import traceback #used for debugging to see the actual errors traceback.print_exc() prints error log

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
instID = 'EthanBass' #'WolfChase1' #Change to be your ID
out_name = os.path.join(path, 'YBphotometry_results_' + instID + '.csv')

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

    clkfig = plt.figure(figsize=(18, 27))

    # Create a df of plottable points from YBloc
    plotloc = []
    for i in range(0, len(YBloc)):
        if 0 <  abs(YBloc['l'][i] - YB_long_pix) < (dxw/2 - 5) and 0 < abs(YBloc['b'][i] - YB_lat_pix) < (dyw/2 -5):
            plotloc.append((YBloc['l'][i] - YB_long_pix + (dxw/2 - 5), YBloc['b'][i] - YB_lat_pix + (dyw/2 -5)))

    #Plot rg image,
    axrg = plt.subplot(2, 4, 1, title='RG (24 + 8 um)', projection=imgw)
    r = orig24
    g = orig
    b = np.zeros(orig.shape)
    axrg.axis('off')
    axrg.imshow(make_lupton_rgb(r, g, b, stretch=200, Q=0))
                #, norm=LogNorm())
    axrg.add_artist(circle)
    for point in range(0, len(plotloc)):
        #axrg.plot(plotloc[point][0], plotloc[point][1], marker='+', markersize=40)
        circleloc = Circle((plotloc[point][0], plotloc[point][1]), YBloc['r'][point], fill = False, color = 'Blue')
        axrg.add_artist(circleloc)

    #Add subtitles to the plot
    clkfig.add_subplot(2,
                       4,
                       2,
                       title='Select coordinates in this ' + wave + ' image',
                       projection=imgw)
    clkfig.suptitle('You are examining the ' + wave + ' image for YB' +
                    str(ybid))
    plt.axis('off')
    plt.imshow(img, cmap='gray', norm=SymLogNorm(linthresh= LinearThreshold))

    #Plot the 70 um
    if math.isnan(orig70.min()) == False and math.isnan(orig70.max()) == False:
        fee = plt.subplot(2, 4, 5, title='70 um', projection=imgw)
        plt.axis('off')
        plt.imshow(orig70,
                   cmap='hot',
                   norm=SymLogNorm(linthresh= LinearThreshold, vmin=orig70.min(), vmax=orig70.max()))
        fee.add_artist(circle1)
    else:
        fee = plt.subplot(2, 4, 5, title='70 um', projection=imgw)
        plt.axis('off')
        plt.imshow(orig70,
                   cmap='hot',
                   norm=SymLogNorm(linthresh= LinearThreshold))
        fee.add_artist(circle1)

    #Plot the 24 um
    if math.isnan(orig24.min()) == False and math.isnan(orig24.max()) == False:
        foo = plt.subplot(2, 4, 6, title='24 um', projection=imgw)
        plt.axis('off')
        plt.imshow(orig24,
                   cmap='hot',
                   norm=SymLogNorm(linthresh= LinearThreshold, vmin=orig24.min(), vmax=orig24.max()))
        foo.add_artist(circle2)
    else:
        foo = plt.subplot(2, 4, 6, title='24 um', projection=imgw)
        plt.axis('off')
        plt.imshow(orig24,
                   cmap='hot',
                   norm=SymLogNorm(linthresh= LinearThreshold))
        foo.add_artist(circle2)
        
    #Plot the 12 um
    if math.isnan(orig12.min()) == False and math.isnan(orig12.max()) == False:
        faa = plt.subplot(2, 4, 7, title='12 um', projection=imgw)
        plt.axis('off')
        plt.imshow(orig12,
                   cmap='hot',
                   norm=SymLogNorm(linthresh= LinearThreshold, vmin=orig12.min(), vmax=orig12.max()))
        faa.add_artist(circle3)
    else:
        faa = plt.subplot(2, 4, 7, title='12 um', projection=imgw)
        plt.axis('off')
        plt.imshow(orig12,
                   cmap='hot',
                   norm=SymLogNorm(linthresh= LinearThreshold))
        faa.add_artist(circle3)

    #Plot the 8um
    if math.isnan(orig.min()) == False and math.isnan(orig.max()) == False:
        fum = plt.subplot(2, 4, 8, title='8 um', projection=imgw)
        plt.axis('off')
        plt.imshow(orig,
                   cmap='hot',
                   norm=SymLogNorm(linthresh= LinearThreshold, vmin=orig.min(), vmax=orig.max()))
        fum.add_artist(circle4)
    else:
        fum = plt.subplot(2, 4, 8, title='8 um', projection=imgw)
        plt.axis('off')
        plt.imshow(orig,
                   cmap='hot',
                   norm=SymLogNorm(linthresh= LinearThreshold))
        fum.add_artist(circle4)
    #cbar = plt.colorbar(format='%05.2f')
    #cbar.set_norm(mynormalize.MyNormalize(vmin=orig.min(),vmax=orig.max(),stretch='linear'))
    #cbar = mycolorbar.DraggableColorbar(cbar,orig)
    #cbar.connect()

    #Plot instructions
    plt.subplot(2, 4, 3)
    plt.axis([0, 10, 0, 10])
    plt.axis('off')
    text = ("*Click in the grey-scale image. ")
    text1 = ("*Left click to add points.")
    text2 = ("*Right click to undo most recent point.")
    text3 = ("*Middle click when finished to exit.")
    text4 = ("*All other images are for reference.")
    text5 = ("*Circles indicate MWP user radius.")
    text6 = ("*You will be prompted to inspect results;")
    text7 = ("type 'n' to continue, anything else to redo.")
    plt.text(1, 8, text, ha='left', wrap=True)
    plt.text(1, 7, text1, ha='left', wrap=True)
    plt.text(1, 6, text2, ha='left', wrap=True)
    plt.text(1, 5, text3, ha='left', wrap=True)
    plt.text(1, 4, text4, ha='left', wrap=True)
    plt.text(1, 3, text5, ha='left', wrap=True)
    plt.text(1, 2, text6, ha='left', wrap=True)
    plt.text(1, 1, text7, ha='left', wrap=True)

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
    plt.imshow(im1, cmap='hot', norm=SymLogNorm(linthresh = LinearThreshold, vmin=im1.min(), vmax=im1.max()))
    plt.axis('off')
    #plt.xlabel('Longitude')
    #plt.ylabel('Latitude')
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
    plt.imshow(im2, cmap='hot', norm=SymLogNorm(linthresh=LinearThreshold, vmin=im1.min(), vmax=im1.max()))
    plt.axis('off')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')

    #Plot the interpolated background
    plt.subplot(2, 2, 3, title='Interpolated image', projection=imw)
    plt.imshow(im3, cmap='hot', norm=SymLogNorm(linthresh=LinearThreshold, vmin=im1.min(), vmax=im1.max()))
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


#Function to examine images and input flags for output file
def make_flags(fim1, fim2, um):
    plt.figure(figsize=(6, 3))

    plt.subplot(1, 2, 1, title='Original Data')
    plt.imshow(fim1,
               cmap='hot',
               norm=SymLogNorm(linthresh=LinearThreshold, vmin=fim1.min(), vmax=fim1.max()))

    plt.subplot(1, 2, 2, title='Bkgrnd Removed')
    plt.imshow(fim2, cmap='hot')

    plt.pause(1)

    flag = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

    print('####################################')
    print('Do you want to note any flags in the output file?')
    print('Select the following flag(s) to apply')
    print('####################################')

    if um == '70' or um == '24' or um == '12':
        foo = 0

        while foo != 9:

            flag1 = "Saturated Image"
            flag2 = "Multiple sources within masked area"
            flag4 = "No obvious source at this wavelength"
            flag6 = "Diffraction Pattern/Star"
            flag7 = "Poor Confidence in Photometry"
            flag8 = "Other/Revisit this source"
            print('flag options:')
            print('[1] ' + flag1)
            print('[2] ' + flag2)
            print('[4] ' + flag4)
            print('[6] ' + flag6)
            print('[7] ' + flag7)
            print('[8] ' + flag8)
            print('[9] Done Flagging')
            print('[10] Clear Flags and Start Over')

            prompts = chain(["Enter a number from the flagging options:"],
                            repeat("Not a flagging option! Try again:"))
            replies = map(input, prompts)
            numeric_strings = filter(str.isnumeric, replies)
            numbers = map(float, numeric_strings)
            is_positive = (0).__lt__
            valid_response = next(filter(is_positive, numbers))
            foo = valid_response

            if foo == 1:
                flag[0] = 1
            if foo == 2:
                flag[1] = 1
            if foo == 4:
                flag[3] = 1
            if foo == 6:
                flag[5] = 1
            if foo == 7:
                flag[6] = 1
            if foo == 8:
                flag[7] = 1
            if foo == 10:
                flag == [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            if foo == 9:
                print("done flagging")
            else:
                foo == 0

    if um == '8':
        foo = 0

        while foo != 9:
            flag1 = "Saturated Image"
            flag2 = "Multiple sources within masked area"
            flag3 = "Filamentary or bubble rim structure"
            flag4 = "No obvious source at this wavelength"
            flag5 = "IRDC Association"
            flag6 = "Diffraction Pattern/Star"
            flag7 = "Poor Confidence in Photometry"
            flag8 = "Other/Revisit this source"

            print('flag options:')
            print('[1] ' + flag1)
            print('[2] ' + flag2)
            print('[3] ' + flag3)
            print('[4] ' + flag4)
            print('[5] ' + flag5)
            print('[6] ' + flag6)
            print('[7] ' + flag7)
            print('[8] ' + flag8)
            print('[9] Done Flagging')
            print('[10] Clear Flags and Start Over')

            prompts = chain(["Enter a number from the flagging options:"],
                            repeat("Not a flagging option! Try again:"))
            replies = map(input, prompts)
            numeric_strings = filter(str.isnumeric, replies)
            numbers = map(float, numeric_strings)
            is_option = (
                0
            ).__lt__  #This provides the condition thst the input should be greater than 0.
            valid_response = next(filter(is_option, numbers))
            foo = valid_response

            if foo == 1:
                flag[0] = 1
            if foo == 2:
                flag[1] = 1
            if foo == 3:
                flag[2] = 1
            if foo == 4:
                flag[3] = 1
            if foo == 5:
                flag[4] = 1
            if foo == 6:
                flag[5] = 1
            if foo == 7:
                flag[6] = 1
            if foo == 8:
                flag[7] = 1
            if foo == 10:
                flag = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            if foo == 9:
                print("done flagging")
            else:
                foo == 0
    return flag


######################################################
# Define my classes                                  #
######################################################
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
        
        ########################################
        #
        ## OOB bug fix
        #
        ########################################
        
        if xmin < 0:
            print('Too small of x values, please click closer to the center')
            print('interpolation will be done with the closest safest value')
            xmin = 0
        
        if ymin < 0:
            print('Too small of y values, please click closer to the center')
            print('interpolation will be done with the closest safest value')
            ymin = 0
        
          
        
        
        #define the 4 areas
        mask = np.zeros_like(img)
        inverse_mask = np.zeros_like(img)
        region_mask = np.zeros_like(img)
        cutout = np.zeros_like(img)

        # filling pixels inside the polygon defined by "vertices" with the fill color
        cv2.fillPoly(mask, vertices, 255)
        #This makes the mask have 0's and 255 in the interior region
        
        #TURN ALL Non-ZERO to NaN
        #This makes the interior region be positive while exterior is 0
        inverse_mask[np.nonzero(mask)] = int(1)  # ones inside poly, zero outside

        
        mask[np.nonzero(mask)] = float('nan') #nan inside, 0 outside
        
        #TURN ALL ZERO to 1
        mask[np.where(mask == 0)] = int(1)  # nan inside poly, 1 outside
        
        #region mask is 0 inside poly, 1 outside
        region_mask = mask
        region_mask = np.nan_to_num(region_mask)  # zero in poly, 1 outside
        
        #cutout is the smaller area that is 5 pixels outside of the polygon
        #result is nan in poly, 1 outside poly, nan outside of rectangle
        cutout[ymin:ymax, xmin:xmax] = mask[ymin:ymax, xmin:xmax]
        #TURN EVERYTHING OUTSIDE THAT RANGE to NaN
        cutout[np.where(cutout == 0)] = float('nan')

        #blank is then nan outside of rectangle, 0 inside polygon, 
        #and normal image values elsewhere
        blank = img * cutout
        self.masked = mask #nan inside poly, 1 outside
        self.blanked = blank #nan inside poly, normal image outside, nan outside rect
        
        #goodvals is list of values that are inside of the rectangle not in polygon
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

        #inverse_mask is ones inside poly, 0 out. 
        #region_mask is 0 inside poly, 1 outside. 
        
        
        self.interp = img * region_mask + fnew_2D * inverse_mask

        #generate the residual image (original - interpolated background)
        self.resid = img - (img * region_mask + fnew_2D * inverse_mask)
        
        #Weird coloration on images is from there exiting negative values inside of the inverse mask
        #ie when fnew_2d > img for values, or that the interpolated values are greater than the source extracted
        #thus if we dont extract the extra light from a star then we should get a bright background


#class that gets the flux from residual images. Unit conversions are applied here.
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

        flux_tot8 = flux_tot8 * str_to_pix8 * 10**6
        flux_tot12 = flux_tot12 * dn_to_jy12
        flux_tot24 = flux_tot24 * str_to_pix8 * 10**6
        flux_tot70 = flux_tot70 #No conversion necesary

        self.um8 = flux_tot8
        self.um12 = flux_tot12
        self.um24 = flux_tot24
        self.um70 = flux_tot70


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
        'YB', 'YB_long', 'YB_lat', 'vertices 8', 'vertices 12', 'vertices 24',
        'vertices 70', '8umphotom', '8flag1', '8flag2', '8flag3', '8flag4',
        '8flag5', '8flag6', '8flag7', '8flag8', '12umphotom', '12flag1', '12flag2', '12flag4',
        '12flag6', '12flag7', '12flag8', '24umphotom', '24flag1', '24flag2', '24flag4', '24flag6',
        '24flag7', '24flag8', '70umphotom', '70flag1', '70flag2', '70flag4', '70flag6', '70flag7',
        '70flag8'
    ]
    writer = csv.DictWriter(output_file, fieldnames=headers)
    output_file.close()

else:
    output_file = open(out_name,
                       'w')  #####opens up files for creating csv file
    writer = csv.DictWriter(
        output_file,
        fieldnames=[
        'YB', 'YB_long', 'YB_lat', 'vertices 8', 'vertices 12', 'vertices 24',
        'vertices 70', '8umphotom', '8flag1', '8flag2', '8flag3', '8flag4',
        '8flag5', '8flag6', '8flag7', '8flag8', '12umphotom', '12flag1', '12flag2', '12flag4',
        '12flag6', '12flag7', '12flag8', '24umphotom', '24flag1', '24flag2', '24flag4', '24flag6',
        '24flag7', '24flag8', '70umphotom', '70flag1', '70flag2', '70flag4', '70flag6', '70flag7',
        '70flag8'
        ],
        lineterminator='\n')
    writer.writeheader()
    writer.writerow({
        'YB': 'ID Number',
        'YB_long': 'degree',
        'YB_lat': 'degree',
        'vertices 8': 'pixel coords',
        'vertices 12': 'pixel coords',
        'vertices 24': 'pixel coords',
        'vertices 70': 'pixel coords',
        '8umphotom': 'Jy',
        '12umphotom': 'Jy',
        '24umphotom': 'Jy',
        '70umphotom': 'Jy',
        '70flag1': 'Saturated',
        '70flag2': 'Multiple sources within YB',
        '70flag4': 'No obvious source at this wavelength',
        '70flag6': 'Star/Diffraction Pattern',
        '70flag7': 'Poor Confidence',
        '70flag8': 'Other/Follow Up',
        '24flag1': 'Saturated',
        '24flag2': 'Multiple sources within YB',
        '24flag4': 'No obvious source at this wavelength',
        '24flag6': 'Star/Diffraction Pattern',
        '24flag7': 'Poor Confidence',
        '24flag8': 'Other/Follow Up',
        '12flag1': 'Saturated',
        '12flag2': 'Multiple sources within YB',
        '12flag4': 'No obvious source at this wavelength',
        '12flag6': 'Star/Diffraction Pattern',
        '12flag7': 'Poor Confidence',
        '12flag8': 'Other/Follow Up',
        '8flag1': 'Saturated',
        '8flag2': 'Multiple sources within YB',
        '8flag3': 'Filament or Bubble Rim',
        '8flag4': 'No obvious source at this wavelength',
        '8flag5': 'IRDC Association',
        '8flag6': 'Star/Diffraction Pattern',
        '8flag7': 'Poor Confidence',
        '8flag8': 'Other/Follow Up'
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
BegYB = 1667

YBlast = 3034

#Determine how much of the scalar negative numbers get
LinearThreshold = 0.001

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
    
    #define a window to zoom in on the YB
    xw = xs = int(YB_long_pix +
                  0.5)  # x coordinates of source and window center
    yw = ys = int(YB_lat_pix +
                  0.5)  # y coordinates of source and window center
    dxw = 100  # x width of window
    dyw = 100  # y width of window

    #find the pixel coordinates LLH and URH corner of zoomed window
    x1 = int(xw - 0.5 * dxw)
    y1 = int(yw - 0.5 * dyw)
    x2 = int(xw + 0.5 * dxw)
    y2 = int(yw + 0.5 * dyw)

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

    wcs8 = cut8.wcs
    wcs12 = cut12.wcs
    wcs24 = cut24.wcs
    wcs70 = cut70.wcs

    #create empty residuals to fill up later
    diff8 = orig * 0
    diff12 = orig12 * 0
    diff24 = orig24 * 0
    diff70 = orig70 * 0

    #create copies of cropped images called workmasks
    workmask = copy.deepcopy(orig)
    workmask12 = copy.deepcopy(orig12)
    workmask24 = copy.deepcopy(orig24)
    workmask70 = copy.deepcopy(orig70)
    ###################################################################
    # Call the classes to draw polygons and perform interpolation     #
    ###################################################################
    try:
        if np.isnan(orig70.min()) == False and np.isnan(orig70.max()) == False:
            check = 'y'
            try:
                while check != 'n':
                    # 70 um image analysis
                    #reset global list coords that gets creaeted in get_coords
                    print('######################################################')
                    print('Beginning the 70 um analysis')
                    coords = []
                    #get the coordinates on 70um image
                    coordinates = get_coords(workmask70, wcs70, '70 um', YB)
        
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
                #flag70 = make_flags(workmask70, interp70.resid, '70') 
                flag70 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
                coord70 = str(coordinates)
                plt.close('all')
            except(ValueError):
                print("There was a problem with the 70 micron image.")
                traceback.print_exc()
                
                coord70 = ' '
        else:
            print('70 micron image is saturated.')
            coord70 = ' '
            flag70 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    
        if math.isnan(orig24.min()) == False and math.isnan(orig24.max()) == False:
            check = 'y'
            
            ##########################################################
            #
            ###  Temporary code to see if we can reuse prior points
            #
            ##########################################################
            if coordinates != []:
                interp24 = do_interp(workmask24, coordinates)
                diff24 = interp24.resid
                #display and save the images for the 24um image
                make_figs(workmask24, interp24.blanked, interp24.interp,
                      interp24.resid, fitcopy24, wcs24, '24_um')

                #Prompt user to accept image or redo
                plt.pause(0.1)
                check = input(
                    'Please consult the residual image. Would you like to reuse the points? Type n to continue, anything else to redo:  '
                    )
                if check != 'n':
                    plt.close('all')
            
            
            ##########################################################
            
            while check != 'n':
                # 24 um image analysis
                #reset global list coords that gets creaeted in get_coords
                print('######################################################')
                print('Beginning the 24 um analysis')
                coords = []
                #get the coordinates on 24um image
                coordinates = get_coords(workmask24, wcs24, '24 um', YB)
    
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
                #do the masking and interpolation on 24um image
                print('starting interp')
                interp24 = do_interp(workmask24, coordinates)
                diff24 = interp24.resid
                #display and save the images for the 24um image
                make_figs(workmask24, interp24.blanked, interp24.interp,
                          interp24.resid, fitcopy24, wcs24, '24_um')
    
                #Prompt user to accept image or redo
                plt.pause(1)
                check = input(
                    'Please consult the residual image. Would you like to redo? Type n to continue, anything else to redo:  '
                )
                if check != 'n':
                    plt.close('all')
            plt.close('all')
            #flag24 = make_flags(workmask24, interp24.resid, '24')
            flag24 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            coord24 = str(coordinates)
            plt.close('all')
        else:
            print('24 micron image is saturated.')
            coord24 = ' '
            flag24 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    
        if math.isnan(orig12.min()) == False and math.isnan(orig12.max()) == False:
            check = 'y'
            
            ##########################################################
            #
            ###  Temporary code to see if we can reuse prior points
            #
            ##########################################################
            if coordinates != []:
                #do the masking and interpolation on 12um image
                interp12 = do_interp(workmask12, coordinates)
                diff12 = interp12.resid
                #display and save the images for the 12um image
                make_figs(workmask12, interp12.blanked, interp12.interp,
                          interp12.resid, fitcopy12, wcs12, '12_um')
                
                #Prompt user to accept image or redo
                plt.pause(0.1)
                check = input(
                    'Please consult the residual image. Would you like to reuse? Type n to continue, anything else to redo:  '
                    )
                if check != 'n':
                    plt.close('all')
            
            
            ##########################################################
            while check != 'n':
                # 12 um image analysis
                #reset global list coords that gets created in get_coords
                print('######################################################')
                print('Beginning the 12 um analysis')
                coords = []
                #get the coordinates on 12um image
                coordinates = get_coords(workmask12, wcs12, '12 um', YB)
                #if no clicks, exit
                if coordinates == []:
                    print(
                        "Timeout limit of 5 minutes reached. Exiting program. Results for most recent YB will not be saved. Please pick 'start where you left off' option next time to re-start this YB."
                    )
                    pickle_out = open("leftYB.pickle", "wb")
                    pickle.dump(currentYB, pickle_out)
                    pickle_out.close()
                    sys.exit()
    
                #do the masking and interpolation on 12um image
                interp12 = do_interp(workmask12, coordinates)
                diff12 = interp12.resid
                #display and save the images for the 12um image
                make_figs(workmask12, interp12.blanked, interp12.interp,
                          interp12.resid, fitcopy12, wcs12, '12_um')
    
                #Prompt user to accept image or redo
                plt.pause(1)
                check = input(
                    'Please consult the residual image. Would you like to redo? Type n to continue, anything else to redo:  '
                )
                if check != 'n':
                    plt.close('all')
            plt.close('all')
            #flag12 = make_flags(workmask12, interp12.resid, '12')
            flag12 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            coord12 = str(coordinates)
            plt.close('all')
        else:
            print('12 micron image is saturated.')
            coord12 = ' '
            flag12 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    
        if math.isnan(orig.min()) == False and math.isnan(orig.max()) == False:
            check = 'y'
            
            ##########################################################
            #
            ###  Temporary code to see if we can reuse prior points
            #
            ##########################################################
            if coordinates != []:
                interp8 = do_interp(workmask, coordinates)
                diff8 = interp8.resid
                #display and save the images for the 8um image
                make_figs(workmask, interp8.blanked, interp8.interp, interp8.resid,
                          fitcopy8, wcs8, '8_um')
                
                #Prompt user to accept image or redo
                plt.pause(1)
                check = input(
                    'Please consult the residual image. Would you like to reuse? Type n to continue, anything else to redo:  '
                    )
                if check != 'n':
                    plt.close('all')
            
            
            ##########################################################
            
            while check != 'n':
                # 8 um image analysis
                #reset global list coords that gets creaeted in get_coords
                print('######################################################')
                print('Beginning the 8 um analysis')
                coords = []
                #get the coordinates on 8um image
                coordinates = get_coords(workmask, wcs8, '8 um', YB)
                #if no clicks, exit
                if coordinates == []:
                    print(
                        "Timeout limit of 5 minutes reached. Exiting program. Results for most recent YB will not be saved. Please pick 'start where you left off' option next time to re-start this YB."
                    )
                    pickle_out = open("leftYB.pickle", "wb")
                    pickle.dump(currentYB, pickle_out)
                    pickle_out.close()
                    sys.exit()
    
                #do the masking and interpolation on 8um image
                interp8 = do_interp(workmask, coordinates)
                diff8 = interp8.resid
                #display and save the images for the 8um image
                make_figs(workmask, interp8.blanked, interp8.interp, interp8.resid,
                          fitcopy8, wcs8, '8_um')
    
                #Prompt user to accept image or redo
                plt.pause(1)
                check = input(
                    'Please consult the residual image. Would you like to redo? Type n to continue, anything else to redo:  '
                )
                if check != 'n':
                    plt.close('all')
            plt.close('all')
            #flag8 = make_flags(workmask, interp8.resid, '8')
            flag8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            coord8 = str(coordinates)
            plt.close('all')
        else:
            print('8 micron image is saturated.')
            coord8 = ' '
            flag8 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    
        ##############################################################################
        # Use residual images to perform photometry and write out to table with flags#
        ##############################################################################
    
        #call the get_flux class
        flux_tot = get_flux(diff8, diff12, diff24, diff70)
    
        #'8umphotom':flux_tot.um8,'12umphotom':flux_tot.um12,'24umphotom':flux_tot.um24}
    
        df = pd.read_csv(out_name)
    
        kk = str(YB)
    
        df.loc[df["YB"] == kk, "8umphotom"] = round(flux_tot.um8, 5)
        df.loc[df["YB"] == kk, '12umphotom'] = round(flux_tot.um12, 5)
        df.loc[df["YB"] == kk, '24umphotom'] = round(flux_tot.um24, 5)
        df.loc[df["YB"] == kk, '70umphotom'] = round(flux_tot.um70, 5)
    
        df.loc[df["YB"] == kk, 'vertices 8'] = coord8
        df.loc[df["YB"] == kk, 'vertices 12'] = coord12
        df.loc[df["YB"] == kk, 'vertices 24'] = coord24
        df.loc[df["YB"] == kk, 'vertices 70'] = coord70
    
        df.loc[df["YB"] == kk, '24flag1'] = flag24[0]
        df.loc[df["YB"] == kk, '24flag2'] = flag24[1]
        df.loc[df["YB"] == kk, '24flag4'] = flag24[3]
        df.loc[df["YB"] == kk, '24flag6'] = flag24[5]
        df.loc[df["YB"] == kk, '24flag7'] = flag24[6]
        df.loc[df["YB"] == kk, '24flag8'] = flag24[7]
    
        df.loc[df["YB"] == kk, '70flag1'] = flag70[0]
        df.loc[df["YB"] == kk, '70flag2'] = flag70[1]
        df.loc[df["YB"] == kk, '70flag4'] = flag70[3]
        df.loc[df["YB"] == kk, '70flag6'] = flag70[5]
        df.loc[df["YB"] == kk, '70flag7'] = flag70[6]
        df.loc[df["YB"] == kk, '70flag8'] = flag70[7]
    
        df.loc[df["YB"] == kk, '12flag1'] = flag12[0]
        df.loc[df["YB"] == kk, '12flag2'] = flag12[1]
        df.loc[df["YB"] == kk, '12flag4'] = flag12[3]
        df.loc[df["YB"] == kk, '12flag6'] = flag12[5]
        df.loc[df["YB"] == kk, '12flag7'] = flag12[6]
        df.loc[df["YB"] == kk, '12flag8'] = flag12[7]
    
        df.loc[df["YB"] == kk, '8flag1'] = flag8[0]
        df.loc[df["YB"] == kk, '8flag2'] = flag8[1]
        df.loc[df["YB"] == kk, '8flag3'] = flag8[2]
        df.loc[df["YB"] == kk, '8flag4'] = flag8[3]
        df.loc[df["YB"] == kk, '8flag5'] = flag8[4]
        df.loc[df["YB"] == kk, '8flag6'] = flag8[5]
        df.loc[df["YB"] == kk, '8flag7'] = flag8[6]
        df.loc[df["YB"] == kk, '8flag8'] = flag8[7]
    
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
    except(ValueError):
        k = k+1
        currentYB = currentYB + 1
        pickle_out = open("leftYB.pickle", "wb")
        pickle.dump(currentYB, pickle_out)
        pickle_out.close()
plt.close('all')
print(
    'Congratulations! You have completed the photometry for all the YBs in your range!'
)
