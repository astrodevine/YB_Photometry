# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 11:37:12 2024

@author: colemaya
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, CheckButtons, TextBox
from matplotlib.patches import Circle
import matplotlib as mpl

plt.ion()
from matplotlib.colors import SymLogNorm
from astropy.io import fits
from astropy import wcs
from astropy.io import ascii
from scipy import interpolate
import sys
import copy
import cv2
import os
import pandas as pd
from astropy.nddata import Cutout2D
import ast

# These lines supress warnings
import warnings
warnings.filterwarnings('ignore')

#For "Pick up where you left off" option
import pickle

#You will need to make sure the following packages have been installed first:

#from tkinter import *
#import tkinter
#conda install -c anaconda tk

#from photutils import centroid_com
#https://photutils.readthedocs.io/en/stable/install.html
#conda install photutils -c astropy

######################################################
#Paths and File Locations                            #
######################################################

#EDIT THIS PATH FOR THE FILE LOCATION ON YOUR MACHINE
# '.' means current directory
path= '.'
path1= '.'
image_name= os.path.join(path, 'GLM_03000+0000_mosaic_I4.fits')
catalog_name= os.path.join(path, 'USE_THIS_CATALOG_ybcat_MWP_with_ID.csv')
#out_name= os.path.join(path1, 'YBphotometry_results.csv')
instID= 'flagcsvfinal' #Change to be your ID
out_name= os.path.join(path, 'YBphotometry_results_' + instID + '.csv')

######################################################
# Define my functions                                #
######################################################

#continue to next image
def cont(event):
    if event.key == 'enter' or event.key == None:
        global forward
        forward= True
        plt.close()
    
#quit program
def savequit(event):
    global done
    done= True
    plt.close()

def entrynumber(entry):
    global startpos
    startpos = int(entry)
    plt.close()
    return startpos

def leftoff(event):
    global startpos
    pickle_in = open("latestYB.pickle", "rb")
    latestYB = pickle.load(
        pickle_in)  #loading in the last YB they completed last time
    startpos = latestYB + 1
    plt.close()
    return startpos

def make_flags(imlist1, imlist2, imlist3, imlist4, imlist5, imgwlist, YB):
    
    umlist= ['70 um', '24 um', '12 um', '8 um']
    
    #flag options
    flagnames= ["Multiple sources within masked area", "No Obvious Source at 70 um",
                "No Obvious Source at 24 um", "No Obvious Source at 12 um", "No Obvious Source at 8 um",
                 "Poor Confidence in 70 um Photometry", "Poor Confidence in 24 um Photometry",
                 "Poor Confidence in 12 um Photometry", "Poor Confidence in 8 um Photometry", "Very Circular and Extended", "Review Source"]
    finish= "Done Flagging"
            
    fig= plt.figure(figsize=(32, 32))
    
    #gather points for other YBs in image
    plotloc = []
    #for i in range(len(YBloc)):
    #    if 0 <  abs(YBloc['l'][i] - YB_long_pix) < dxw and 0 < abs(YBloc['b'][i] - YB_lat_pix) < dyw:
    #        plotloc.append((YBloc['l'][i] - YB_long_pix + dxw, YBloc['b'][i] - YB_lat_pix + dyw, YBloc['r'][i]))
            
    cmap = mpl.colormaps.get_cmap('hot')
    cmap.set_bad(color='purple')
    
    cmap2 = mpl.colormaps.get_cmap('hot')
    cmap2.set_bad(color='blue')
    
    #plot original images
    for i in range(len(imlist1)):
                
        img1= imlist1[i]
        circle = Circle((dxw, dyw), YB_rad_pix, fill=False)
        
        origplot = plt.subplot(4, 6, 1 + 6*i, projection = imgwlist[i])
        if i == 0:
            plt.title('Original Image', fontsize= 15)
        #plot image as normal if not saturated
        if ~np.isnan(img1.min()) and ~np.isnan(img1.max()):
            plt.imshow(img1, cmap='hot',
               norm=SymLogNorm(linthresh= LinearThreshold, vmin=img1.min(), vmax=img1.max()))
        #remove vmin and vmax values if image is saturated
        else:
            plt.imshow(img1, cmap=cmap,
               norm=SymLogNorm(linthresh= LinearThreshold))
        origplot.add_artist(circle)
        
        #mark focused YB with black x, other YBs with blue x
        plt.plot(dxw, dyw, color= 'black', marker= 'x', markersize= 12)
        for point in range(len(plotloc)):
            plt.plot(plotloc[point][0], plotloc[point][1], color= 'blue', marker= 'x', markersize= 10),
            circleloc = Circle((plotloc[point][0], plotloc[point][1]), plotloc[point][2], fill = False, color = 'Blue')
            origplot.add_artist(circleloc)
        
        plt.axis('off')
        plt.text(-30, dyw, umlist[i], fontsize= 15)
    
    #plot average mask images
    for i in range(len(imlist2)):
        img2= imlist2[i]
        plt.subplot(4, 6, 2 + 6*i, projection= imgwlist[i])
        if i == 0:
            plt.title('Average Mask', fontsize= 15)
        #plot background subtracted image if not saturated
        if np.max(img2) != 0:
            plt.imshow(imlist2[i], cmap= cmap2)
        #"Saturated Image" text over black image if saturated
        else:
            plt.imshow(imlist2[i], cmap= cmap2)
            plt.text(8, dyw, 'Saturated', fontsize= 20, color= 'white')
            plt.text(23, dyw - 20, 'Image', fontsize= 20, color= 'white')
        plt.axis('off')
        
    #plot residual images
    for i in range(len(imlist3)):
        img3= imlist3[i]
        plt.subplot(4, 6, 3 + 6*i, projection= imgwlist[i])
        if i == 0:
            plt.title('Background Removed', fontsize= 15)
        #plot background subtracted image if not saturated
        if np.max(img3) != 0:
            plt.imshow(img3, cmap= cmap)
               #norm=SymLogNorm(linthresh= LinearThreshold, vmin=img3.min(), vmax=img3.max()))
        #"Saturated Image" text over black image if saturated
        else:
            plt.imshow(imlist3[i], cmap= cmap)
            plt.text(8, dyw, 'Saturated', fontsize= 20, color= 'white')
            plt.text(23, dyw - 20, 'Image', fontsize= 20, color= 'white')
        plt.axis('off')
        
    #plot interpolated images
    for i in range(len(imlist4)):
        img4= imlist4[i]
        plt.subplot(4, 6, 4 + 6*i, projection= imgwlist[i])
        if i == 0:
            plt.title('Background Only', fontsize= 15)
        #plot background subtracted image if not saturated
        if np.max(img4) != 0:
            plt.imshow(img4, cmap=cmap)
               #norm=SymLogNorm(linthresh= LinearThreshold, vmin=img4.min(), vmax=img4.max()))
        #"Saturated Image" text over black image if saturated
        else:
            plt.imshow(imlist4[i], cmap= cmap)
            plt.text(8, dyw, 'Saturated', fontsize= 20, color= 'white')
            plt.text(23, dyw - 20, 'Image', fontsize= 20, color= 'white')
        plt.axis('off')
        
    #plot heatmaps
    for i in range(len(imlist5)):
         img5= imlist5[i]
         plt.subplot(4, 6, 5 + 6*i, projection= imgwlist[i])
         if i == 0:
             plt.title('Selection Heatmap', fontsize= 15)
         #plot background subtracted image if not saturated
         if np.max(img5) != 0:
             plt.imshow(img5, cmap=cmap)
                #norm=SymLogNorm(linthresh= LinearThreshold, vmin=img5.min(), vmax=img5.max()))
         #"Saturated Image" text over black image if saturated
         else:
             plt.imshow(imlist5[i], cmap= cmap)
             plt.text(8, dyw, 'Saturated', fontsize= 20, color= 'white')
             plt.text(23, dyw - 20, 'Image', fontsize= 20, color= 'white')
         plt.axis('off')
    
    fig.suptitle('Flagging for YB ' + str(YB), fontsize= 15)
    
    #space for flagging options
    plt.subplot(4, 6, 6)
    plt.title('Flagging Options', fontsize= 15)
    plt.axis('off')
    
    #adds flags to list when boxes are checked
    def marked(label):
        index= flagnames.index(label)
        checked[index]= not checked[index]
        if checked[index] == True:
            flag[index]= 1
        elif checked[index] == False:
            flag[index]= 0
        
    flag= [0] * len(flagnames)
    
    #checkboxes for flags
    checked= [False] * len(flagnames)
    checkaxes= plt.axes([.775, .575, .21, .305])
    checks= CheckButtons(checkaxes, flagnames, checked)
    checks.on_clicked(marked)
    checks.rectangles
    for i in range(len(flagnames)):
        checks.labels[i].set_fontsize(10)
    
    #continue button
    buttonaxes= plt.axes([.775, .475, .1, .075])
    contbutton= Button(buttonaxes, finish)
    contbutton.on_clicked(cont)
    contbutton.label.set_fontsize(15)
    
    #quit button
    button2axes= plt.axes([.885, .475, .1, .075])
    exitbutton= Button(button2axes, 'Quit')
    exitbutton.on_clicked(savequit)
    exitbutton.label.set_fontsize(15)
    
    #allows you to press enter to continue instead of clicking the button
    plt.connect('key_press_event', cont)
    
    global forward
    forward= None
    plt.show()
    # Wait for user to press a button
    while forward is None and not done:
        plt.pause(0.1)
            
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
        if l > 1.5 and l <= 4.5:
            path8= os.path.join(path,
                                 'mosaics/GLM_00300+0000_mosaic_I4.fits')
            path12= os.path.join(path,
                                  'mosaics/WISE_12um_00300_mosaic.fits')
            path24= os.path.join(
                                   path,
                                      'mosaics/MIPSGAL_24um_00300_mosaic.fits')
            path70= os.path.join(
                                   path,
                                       'mosaics/PACS_70um_00300_mosaic.fits')
        elif l > 4.5 and l <= 7.5:
            path8= os.path.join(path,
                                  'mosaics/GLM_00600+0000_mosaic_I4.fits')
            path12= os.path.join(path,
                                   'mosaics/WISE_12um_00600_mosaic.fits')
            path24= os.path.join(
                                    path,
                                       'mosaics/MIPSGAL_24um_00600_mosaic.fits')
            path70= os.path.join(
                                    path,
                                        'mosaics/PACS_70um_00600_mosaic.fits')
        elif l > 7.5 and l <= 10.5:
            path8= os.path.join(path,
                                  'mosaics/GLM_00900+0000_mosaic_I4.fits')
            path12= os.path.join(path,
                                   'mosaics/WISE_12um_00900_mosaic.fits')
            path24= os.path.join(
                                    path,
                                       'mosaics/MIPSGAL_24um_00900_mosaic.fits')
            path70= os.path.join(
                                    path,
                                        'mosaics/PACS_70um_00900_mosaic.fits')
        elif l > 10.5 and l <= 13.5:
            path8= os.path.join(path,
                                  'mosaics/GLM_01200+0000_mosaic_I4.fits')
            path12= os.path.join(path,
                                   'mosaics/WISE_12um_01200_mosaic.fits')
            path24= os.path.join(
                                    path,
                                       'mosaics/MIPSGAL_24um_01200_mosaic.fits')
            path70= os.path.join(
                                    path,
                                        'mosaics/PACS_70um_01200_mosaic.fits')
        elif l > 13.5 and l <= 16.5:
            path8= os.path.join(path,
                                  'mosaics/GLM_01500+0000_mosaic_I4.fits')
            path12= os.path.join(path,
                                   'mosaics/WISE_12um_01500_mosaic.fits')
            path24= os.path.join(
                                    path,
                                       'mosaics/MIPSGAL_24um_01500_mosaic.fits')
            path70= os.path.join(
                                    path,
                                        'mosaics/PACS_70um_01500_mosaic.fits')
        elif l > 16.5 and l <= 19.5:
            path8= os.path.join(path,
                                  'mosaics/GLM_01800+0000_mosaic_I4.fits')
            path12= os.path.join(path,
                                   'mosaics/WISE_12um_01800_mosaic.fits')
            path24= os.path.join(
                                    path,
                                       'mosaics/MIPSGAL_24um_01800_mosaic.fits')
            path70= os.path.join(
                                    path,
                                        'mosaics/PACS_70um_01800_mosaic.fits')
        #Adding mosaics 021, 024, 027 on 10/17/23.
        elif l > 19.5 and l <= 22.5:
            path8= os.path.join(path,
                                 'mosaics/GLM_02100+0000_mosaic_I4.fits')
            path12= os.path.join(path,
                                  'mosaics/WISE_12um_02100_mosaic.fits')
            path24= os.path.join(
                path,
                'mosaics/MIPSGAL_24um_02100_mosaic.fits')
            path70= os.path.join(
                path,
                'mosaics/PACS_70um_02100_mosaic.fits')    
        elif l > 22.5 and l <= 25.5:
            path8= os.path.join(path,
                                 'mosaics/GLM_02400+0000_mosaic_I4.fits')
            path12= os.path.join(path,
                                  'mosaics/WISE_12um_02400_mosaic.fits')
            path24= os.path.join(
                path,
                'mosaics/MIPSGAL_24um_02400_mosaic.fits')
            path70= os.path.join(
                path,
                'mosaics/PACS_70um_02400_mosaic.fits')    
        elif l > 25.5 and l <= 28.5:
            path8= os.path.join(path,
                                 'mosaics/GLM_02700+0000_mosaic_I4.fits')
            path12= os.path.join(path,
                                  'mosaics/WISE_12um_02700_mosaic.fits')
            path24= os.path.join(
                path,
                'mosaics/MIPSGAL_24um_02700_mosaic.fits')
            path70= os.path.join(
                path,
                'mosaics/PACS_70um_02700_mosaic.fits')    
        elif l > 28.5 and l <= 31.5:
            path8= os.path.join(path,
                                 'mosaics/GLM_03000+0000_mosaic_I4.fits')
            path12= os.path.join(path,
                                  'mosaics/WISE_12um_03000_mosaic.fits')
            path24= os.path.join(
                path,
                'mosaics/MIPSGAL_03000_mosaic_reprojected.fits')
            path70= os.path.join(
                path,
                'mosaics/PACS_70um_03000_mosaic.fits')
        elif l > 31.5 and l <= 34.5:
            path8= os.path.join(path,
                                 'mosaics/GLM_03300+0000_mosaic_I4.fits')
            path12= os.path.join(path,
                                  'mosaics/WISE_12um_03300_mosaic.fits')
            path24= os.path.join(
                path,
                'mosaics/MIPSGAL_03300_mosaic_reprojected.fits')
            path70= os.path.join(
                path,
                'mosaics/PACS_70um_03300_mosaic.fits')
        elif l > 34.5 and l <= 37.5:
            path8= os.path.join(path,
                                 'mosaics/GLM_03600+0000_mosaic_I4.fits')
            path12= os.path.join(path,
                                  'mosaics/WISE_12um_03600_mosaic.fits')
            path24= os.path.join(
                path,
                'mosaics/MIPSGAL_03600_mosaic_reprojected.fits')
            path70= os.path.join(
                path,
                'mosaics/PACS_70um_03600_mosaic.fits')
        elif l > 37.5 and l <= 40.5:
            path8= os.path.join(path,
                                 'mosaics/GLM_03900+0000_mosaic_I4.fits')
            path12= os.path.join(path,
                                  'mosaics/WISE_12um_03900_mosaic.fits')
            path24= os.path.join(
                path,
                'mosaics/MIPSGAL_03900_mosaic_reprojected.fits')
            path70= os.path.join(
                path,
                'mosaics/PACS_70um_03900_mosaic.fits')
        elif l > 40.5 and l <= 43.5:
            path8= os.path.join(path,
                                 'mosaics/GLM_04200+0000_mosaic_I4.fits')
            path12= os.path.join(path,
                                  'mosaics/WISE_12um_04200_mosaic.fits')
            path24= os.path.join(
                path,
                'mosaics/MIPSGAL_24um_04200_mosaic.fits')
            path70= os.path.join(
                path,
                'mosaics/PACS_70um_04200_mosaic.fits')
        elif l > 43.5 and l <= 46.5:
            path8= os.path.join(path,
                                 'mosaics/GLM_04500+0000_mosaic_I4.fits')
            path12= os.path.join(path,
                                  'mosaics/WISE_12um_04500_mosaic.fits')
            path24= os.path.join(
                path,
                'mosaics/MIPSGAL_24um_04500_mosaic.fits')
            path70= os.path.join(
                path,
                'mosaics/PACS_70um_04500_mosaic.fits')       
        elif l > 46.5 and l <= 49.5:
            path8= os.path.join(path,
                                 'mosaics/GLM_04800+0000_mosaic_I4.fits')
            path12= os.path.join(path,
                                  'mosaics/WISE_12um_04800_mosaic.fits')
            path24= os.path.join(
                path,
                'mosaics/MIPSGAL_24um_04800_mosaic.fits')
            path70= os.path.join(
                path,
                'mosaics/PACS_70um_04800_mosaic.fits') 
        elif l > 49.5 and l <= 52.5:
            path8= os.path.join(path,
                                 'mosaics/GLM_05100+0000_mosaic_I4.fits')
            path12= os.path.join(path,
                                  'mosaics/WISE_12um_05100_mosaic.fits')
            path24= os.path.join(
                path,
                'mosaics/MIPSGAL_24um_05100_mosaic.fits')
            path70= os.path.join(
                path,
                'mosaics/PACS_70um_05100_mosaic.fits')
        elif l > 52.5 and l <= 55.5:  
            path8= os.path.join(path,
                                 'mosaics/GLM_05400+0000_mosaic_I4.fits')
            path12= os.path.join(path,
                                  'mosaics/WISE_12um_05400_mosaic.fits')
            path24= os.path.join(
                path,
                'mosaics/MIPSGAL_24um_05400_mosaic.fits')
            path70= os.path.join(
                path,
                'mosaics/PACS_70um_05400_mosaic.fits')
        elif l > 55.5 and l <= 58.5:
            path8= os.path.join(path,
                                 'mosaics/GLM_05700+0000_mosaic_I4.fits')
            path12= os.path.join(path,
                                  'mosaics/WISE_12um_05700_mosaic.fits')
            path24= os.path.join(
                path,
                'mosaics/MIPSGAL_24um_05700_mosaic.fits')
            path70= os.path.join(
                    path,
                    'mosaics/PACS_70um_05700_mosaic.fits')   
        elif l > 58.5 and l <= 61.5:
            path8= os.path.join(path,
                                     'mosaics/GLM_06000+0000_mosaic_I4.fits')
            path12= os.path.join(path,
                                      'mosaics/WISE_12um_06000_mosaic.fits')
            path24= os.path.join(
                    path,
                    'mosaics/MIPSGAL_24um_06000_mosaic.fits')
            path70= os.path.join(
                    path,
                    'mosaics/PACS_70um_06000_mosaic.fits')  
        elif l > 61.5 and l <= 64.5:
            path8= os.path.join(path,
                                     'mosaics/GLM_06300+0000_mosaic_I4.fits')
            path12= os.path.join(path,
                                      'mosaics/WISE_12um_06300_mosaic.fits')
            path24= os.path.join(
                    path,
                    'mosaics/MIPSGAL_24um_06300_mosaic.fits')
            path70= os.path.join(
                    path,
                    'mosaics/PACS_70um_06300_mosaic.fits')                   
        elif l > 64.5 and l <= 65.5:
            path8= os.path.join(path,
                                     'mosaics/GLM_06600+0000_mosaic_I4.fits')
            path12= os.path.join(path,
                                      'mosaics/WISE_12um_06600_mosaic.fits')
            path24= os.path.join(
                    path,
                    'mosaics/MIPSGAL_24um_06600_mosaic.fits')
            path70= os.path.join(
                    path,
                    'mosaics/PACS_70um_06600_mosaic.fits') 

        # The following were added for Cyg-X by GW-C on 2/7/24.
        elif l > 75.5 and l <= 76.5:
            path8= os.path.join(path,
                                 'mosaics/CYGX_08um_07500+0050_mosaic.fits')
            path12= os.path.join(path,
                                  'mosaics/CYGX_12um_07500+0050_mosaic.fits')
            path24= os.path.join(
                path,
                'mosaics/CYGX_24um_07500+0050_mosaic.fits')
            path70= os.path.join(
                path,
                'mosaics/CYGX_70um_07500+0050_mosaic.fits')
        elif l > 76.5 and l <= 79.5 and b < 0.82:
            path8= os.path.join(path,
                                 'mosaics/CYGX_08um_07800-0085_mosaic.fits')
            path12= os.path.join(path,
                                  'mosaics/CYGX_12um_07800-0085_mosaic.fits')
            path24= os.path.join(
                path,
                'mosaics/CYGX_24um_07800-0085_mosaic.fits')
            path70= os.path.join(
                path,
                'mosaics/CYGX_70um_07800-0085_mosaic.fits')
        elif l > 76.5 and l <= 79.5 and b >= 0.82:
            path8= os.path.join(path,
                                 'mosaics/CYGX_08um_07800+0250_mosaic.fits')
            path12= os.path.join(path,
                                  'mosaics/CYGX_12um_07800+0250_mosaic.fits')
            path24= os.path.join(
                path,
                'mosaics/CYGX_24um_07800+0250_mosaic.fits')
            path70= os.path.join(
                path,
                'mosaics/CYGX_70um_07800+0250_mosaic.fits')
        elif l > 79.5 and l <= 82.5 and b < 0.82:
            path8= os.path.join(path,
                                 'mosaics/CYGX_08um_08100-0070_mosaic.fits')
            path12= os.path.join(path,
                                  'mosaics/CYGX_12um_08100-0070_mosaic.fits')
            path24= os.path.join(
                path,
                'mosaics/CYGX_24um_08100-0070_mosaic.fits')
            path70= os.path.join(
                path,
                'mosaics/CYGX_70um_08100-0070_mosaic.fits')
        elif l > 79.5 and l <= 82.5 and b >= 0.82:
            path8= os.path.join(path,
                                 'mosaics/CYGX_08um_08100+0235_mosaic.fits')
            path12= os.path.join(path,
                                  'mosaics/CYGX_12um_08100+0235_mosaic.fits')
            path24= os.path.join(
                path,
                'mosaics/CYGX_24um_08100+0235_mosaic.fits')
            path70= os.path.join(
                path,
                'mosaics/CYGX_70um_08100+0235_mosaic.fits')
        elif l > 82.5 and l <= 83.0:
            path8= os.path.join(path,
                                 'mosaics/CYGX_08um_08400+0005_mosaic.fits')
            path12= os.path.join(path,
                                  'mosaics/CYGX_12um_08400+0005_mosaic.fits')
            path24= os.path.join(
                path,
                'mosaics/CYGX_24um_08400+0005_mosaic.fits')
            path70= os.path.join(
                path,
                'mosaics/CYGX_70um_08400+0005_mosaic.fits')

        #The following are for the SMOG region.  
        #GWC: Something went wonky on 2/7/24 -- need to revisit how to cover SMOG.
        elif l > 101.0 and l <= 105.59 and b < 3.06:
            path8= os.path.join(path,
                                 'mosaics/SMOG_08um_10300_mosaic.fits')
            path12= os.path.join(path,
                                  'mosaics/SMOG_12um_10300_mosaic.fits')
            path24= os.path.join(
                path,
                'mosaics/SMOG_24um_10300_mosaic.fits')
            path70= os.path.join(
                path,
                'mosaics/SMOG_PACS_70um_10300_mosaic.fits')
            # Replaced 'mosaics/SMOG_70um_10300_mosaic.fits') with PACS on 7/7/23
        elif l > 101.0 and l <= 105.59 and b >= 3.06:
            path8= os.path.join(path,
                                 'mosaics/SMOG_08um_10300_mosaic.fits')
            path12= os.path.join(path,
                                  'mosaics/SMOG_12um_10300_mosaic.fits')
            path24= os.path.join(
                path,
                'mosaics/SMOG_24um_10300_mosaic_high_b.fits')
            path70= os.path.join(
                path,
                'mosaics/SMOG_PACS_70um_10300_mosaic.fits')
        elif l > 105.59 and l <= 110.2:
            path8= os.path.join(path,
                                 'mosaics/SMOG_08um_10700_mosaic.fits')
            path12= os.path.join(path,
                                  'mosaics/SMOG_12um_10700_mosaic.fits')
            path24= os.path.join(
                path,
                'mosaics/SMOG_24um_10700_mosaic.fits')
            path70= os.path.join(
                path,
                'mosaics/SMOG_PACS_70um_10700_mosaic.fits')
            # Replaced 'mosaics/SMOG_70um_10700_mosaic.fits') with PACS on 7/7/23
        
        elif l > 294.8 and l <= 295.5:
            path8= os.path.join(path,
                                     'mosaics/GLM_29400+0000_mosaic_I4.fits')
            path12= os.path.join(path,
                                      'mosaics/WISE_12um_29400_mosaic.fits')
            path24= os.path.join(
                    path,
                    'mosaics/MIPSGAL_24um_29400_mosaic.fits')
            path70= os.path.join(
                    path,
                    'mosaics/PACS_70um_29400_mosaic.fits')              
        elif l > 295.5 and l <= 298.5:
            path8= os.path.join(path,
                                     'mosaics/GLM_29700+0000_mosaic_I4.fits')
            path12= os.path.join(path,
                                      'mosaics/WISE_12um_29700_mosaic.fits')
            path24= os.path.join(
                    path,
                    'mosaics/MIPSGAL_24um_29700_mosaic.fits')
            path70= os.path.join(
                    path,
                    'mosaics/PACS_70um_29700_mosaic.fits')
        elif l > 298.5 and l <= 301.5:
            path8= os.path.join(path,
                                     'mosaics/GLM_30000+0000_mosaic_I4.fits')
            path12= os.path.join(path,
                                      'mosaics/WISE_12um_30000_mosaic.fits')
            path24= os.path.join(
                    path,
                    'mosaics/MIPSGAL_24um_30000_mosaic.fits')
            path70= os.path.join(
                    path,
                    'mosaics/PACS_70um_30000_mosaic.fits')
        elif l > 301.5 and l <= 304.5:
            path8= os.path.join(path,
                                     'mosaics/GLM_30300+0000_mosaic_I4.fits')
            path12= os.path.join(path,
                                      'mosaics/WISE_12um_30300_mosaic.fits')
            path24= os.path.join(
                    path,
                    'mosaics/MIPSGAL_24um_30300_mosaic.fits')
            path70= os.path.join(
                    path,
                    'mosaics/PACS_70um_30300_mosaic.fits')
        elif l > 304.5 and l <= 307.5:
            path8= os.path.join(path,
                                     'mosaics/GLM_30600+0000_mosaic_I4.fits')
            path12= os.path.join(path,
                                      'mosaics/WISE_12um_30600_mosaic.fits')
            path24= os.path.join(
                    path,
                    'mosaics/MIPSGAL_24um_30600_mosaic.fits')
            path70= os.path.join(
                    path,
                    'mosaics/PACS_70um_30600_mosaic.fits')       
        #Added these mosaics on 7/11/23. For some reason, many of the elif statements
        #for regions I've done are not here. I added these above on 7/26/23. I had to
        #copy them over from ExpertPhotom.py
        #Adding more as I complete photometry for eacj sector.
        elif l > 307.5 and l <= 310.5:
            path8= os.path.join(path,
                                     'mosaics/GLM_30900+0000_mosaic_I4.fits')
            path12= os.path.join(path,
                                      'mosaics/WISE_12um_30900_mosaic.fits')
            path24= os.path.join(
                    path,
                    'mosaics/MIPSGAL_24um_30900_mosaic.fits')
            path70= os.path.join(
                    path,
                    'mosaics/PACS_70um_30900_mosaic.fits')
        elif l > 310.5 and l <= 313.5:
            path8= os.path.join(path,
                                 'mosaics/GLM_31200+0000_mosaic_I4.fits')
            path12= os.path.join(path,
                                  'mosaics/WISE_12um_31200_mosaic.fits')
            path24= os.path.join(
                        path,
                        'mosaics/MIPSGAL_24um_31200_mosaic.fits')
            path70= os.path.join(
                        path,
                        'mosaics/PACS_70um_31200_mosaic.fits')
        elif l > 313.5 and l <= 316.5:
            path8= os.path.join(path,
                                  'mosaics/GLM_31500+0000_mosaic_I4.fits')
            path12= os.path.join(path,
                                   'mosaics/WISE_12um_31500_mosaic.fits')
            path24= os.path.join(
                         path,
                         'mosaics/MIPSGAL_24um_31500_mosaic.fits')
            path70= os.path.join(
                         path,
                         'mosaics/PACS_70um_31500_mosaic.fits')    
        elif l > 316.5 and l <= 319.5:
            path8= os.path.join(path,
                                   'mosaics/GLM_31800+0000_mosaic_I4.fits')
            path12= os.path.join(path,
                                    'mosaics/WISE_12um_31800_mosaic.fits')
            path24= os.path.join(
                          path,
                          'mosaics/MIPSGAL_24um_31800_mosaic.fits')
            path70= os.path.join(
                          path,
                          'mosaics/PACS_70um_31800_mosaic.fits')      
        elif l > 319.5 and l <= 322.5:
            path8= os.path.join(path,
                                   'mosaics/GLM_32100+0000_mosaic_I4.fits')
            path12= os.path.join(path,
                                    'mosaics/WISE_12um_32100_mosaic.fits')
            path24= os.path.join(
                          path,
                          'mosaics/MIPSGAL_24um_32100_mosaic.fits')
            path70= os.path.join(
                          path,
                          'mosaics/PACS_70um_32100_mosaic.fits')   
        elif l > 322.5 and l <= 325.5:
            path8= os.path.join(path,
                                   'mosaics/GLM_32400+0000_mosaic_I4.fits')
            path12= os.path.join(path,
                                    'mosaics/WISE_12um_32400_mosaic.fits')
            path24= os.path.join(
                          path,
                          'mosaics/MIPSGAL_24um_32400_mosaic.fits')
            path70= os.path.join(
                          path,
                          'mosaics/PACS_70um_32400_mosaic.fits')          
        elif l > 325.5 and l <= 328.5:
            path8= os.path.join(path,
                                   'mosaics/GLM_32700+0000_mosaic_I4.fits')
            path12= os.path.join(path,
                                    'mosaics/WISE_12um_32700_mosaic.fits')
            path24= os.path.join(
                          path,
                          'mosaics/MIPSGAL_24um_32700_mosaic.fits')
            path70= os.path.join(
                          path,
                          'mosaics/PACS_70um_32700_mosaic.fits')         
        elif l > 328.5 and l <= 331.5:
            path8= os.path.join(path,
                                       'mosaics/GLM_33000+0000_mosaic_I4.fits')
            path12= os.path.join(path,
                                        'mosaics/WISE_12um_33000_mosaic.fits')
            path24= os.path.join(
                              path,
                              'mosaics/MIPSGAL_24um_33000_mosaic.fits')
            path70= os.path.join(
                              path,
                              'mosaics/PACS_70um_33000_mosaic.fits')   
        elif l > 331.5 and l <= 334.5:
            path8= os.path.join(path,
                                       'mosaics/GLM_33300+0000_mosaic_I4.fits')
            path12= os.path.join(path,
                                        'mosaics/WISE_12um_33300_mosaic.fits')
            path24= os.path.join(
                              path,
                              'mosaics/MIPSGAL_24um_33300_mosaic.fits')
            path70= os.path.join(
                              path,
                              'mosaics/PACS_70um_33300_mosaic.fits')           
        elif l > 334.5 and l <= 337.5:
            path8= os.path.join(path,
                                    'mosaics/GLM_33600+0000_mosaic_I4.fits')
            path12= os.path.join(path,
                                  'mosaics/WISE_12um_33600_mosaic.fits')
            path24= os.path.join(
                                  path,
                                  'mosaics/MIPSGAL_24um_33600_mosaic.fits')
            path70= os.path.join(
                                  path,
                                  'mosaics/PACS_70um_33600_mosaic.fits')  
        elif l > 337.5 and l <= 340.5:
            path8= os.path.join(path,
                                     'mosaics/GLM_33900+0000_mosaic_I4.fits')
            path12= os.path.join(path,
                                   'mosaics/WISE_12um_33900_mosaic.fits')
            path24= os.path.join(
                                   path,
                                   'mosaics/MIPSGAL_24um_33900_mosaic.fits')
            path70= os.path.join(
                                   path,
                                   'mosaics/PACS_70um_33900_mosaic.fits')        
        elif l > 340.5 and l <= 343.5:
            path8= os.path.join(path,
                                     'mosaics/GLM_34200+0000_mosaic_I4.fits')
            path12= os.path.join(path,
                                   'mosaics/WISE_12um_34200_mosaic.fits')
            path24= os.path.join(
                                   path,
                                   'mosaics/MIPSGAL_24um_34200_mosaic.fits')
            path70= os.path.join(
                                   path,
                                   'mosaics/PACS_70um_34200_mosaic.fits') 
        elif l > 343.5 and l <= 346.5:
            path8= os.path.join(path,
                                 'mosaics/GLM_34500+0000_mosaic_I4.fits')
            path12= os.path.join(path,
                                  'mosaics/WISE_12um_34500_mosaic.fits')
            path24= os.path.join(
                                   path,
                                      'mosaics/MIPSGAL_24um_34500_mosaic.fits')
            path70= os.path.join(
                                   path,
                                       'mosaics/PACS_70um_34500_mosaic.fits') 
        elif l > 346.5 and l <= 349.5:
            path8= os.path.join(path,
                                'mosaics/GLM_34800+0000_mosaic_I4.fits')
            path12= os.path.join(path,
                                 'mosaics/WISE_12um_34800_mosaic.fits')
            path24= os.path.join(
                                  path,
                                     'mosaics/MIPSGAL_24um_34800_mosaic.fits')
            path70= os.path.join(
                                  path,
                                      'mosaics/PACS_70um_34800_mosaic.fits') 
        elif l > 349.5 and l <= 352.5:
            path8= os.path.join(path,
                                'mosaics/GLM_35100+0000_mosaic_I4.fits')
            path12= os.path.join(path,
                                 'mosaics/WISE_12um_35100_mosaic.fits')
            path24= os.path.join(
                                  path,
                                     'mosaics/MIPSGAL_24um_35100_mosaic.fits')
            path70= os.path.join(
                                  path,
                                      'mosaics/PACS_70um_35100_mosaic.fits') 
        elif l > 352.5 and l <= 355.5:
            path8= os.path.join(path,
                                'mosaics/GLM_35400+0000_mosaic_I4.fits')
            path12= os.path.join(path,
                                 'mosaics/WISE_12um_35400_mosaic.fits')
            path24= os.path.join(
                                  path,
                                     'mosaics/MIPSGAL_24um_35400_mosaic.fits')
            path70= os.path.join(
                                  path,
                                      'mosaics/PACS_70um_35400_mosaic.fits') 
        elif l > 355.5 and l <= 358.5:
            path8= os.path.join(path,
                                 'mosaics/GLM_35700+0000_mosaic_I4.fits')
            path12= os.path.join(path,
                                  'mosaics/WISE_12um_35700_mosaic.fits')
            path24= os.path.join(
                                   path,
                                      'mosaics/MIPSGAL_24um_35700_mosaic.fits')
            path70= os.path.join(
                                   path,
                                       'mosaics/PACS_70um_35700_mosaic.fits')    
        elif (l > 358.5 and l <= 360.1) or (l > -0.1 and l <= 1.5):
            path8= os.path.join(path,
                                 'mosaics/GLM_00000+0000_mosaic_I4.fits')
            path12= os.path.join(path,
                                  'mosaics/WISE_12um_00000_mosaic.fits')
            path24= os.path.join(
                                   path,
                                      'mosaics/MIPSGAL_24um_00000_mosaic.fits')
            path70= os.path.join(
                                   path,
                                       'mosaics/PACS_70um_00000_mosaic.fits')
        else:
            # GWC revised print statement from "outside the pilot..."
            print('Your YB is outside the region.')
            print('Please try again.')
            sys.exit()

        temp= fits.open(path8)[0]
        self.um8= temp
        self.um8data= temp.data
        self.um8w= wcs.WCS(temp.header)
        temp= fits.open(path12)[0]
        self.um12= temp
        self.um12data= temp.data
        self.um12w= wcs.WCS(temp.header)
        temp= fits.open(path24)[0]
        self.um24= temp
        self.um24data= temp.data
        self.um24w= wcs.WCS(temp.header)
        temp= fits.open(path70)[0]
        self.um70= temp
        self.um70data= temp.data
        self.um70w= wcs.WCS(temp.header)


#class that does the masking and interpolation, returns masked, blanked, interpolated, and residual
class do_interp():
    def __init__(self, img, verts):

        #use the clicked values from the user to create a NaN mask
        vertices= np.array([verts], dtype=np.int32)
        xyvals= np.array(verts, dtype=np.int32)
        xmin= min(xyvals[:, 0]) - 5
        xmax= max(xyvals[:, 0]) + 5
        ymin= min(xyvals[:, 1]) - 5
        ymax= max(xyvals[:, 1]) + 5
        #print(xmin, xmax, ymin, ymax)
        mask= np.zeros_like(img)
        inverse_mask= np.zeros_like(img)
        #region_mask= np.zeros_like(img)
        cutout= np.zeros_like(img)

        # filling pixels inside the polygon defined by "vertices" with the fill color
        cv2.fillPoly(mask, vertices, 255)
        #TURN ALL Non-ZERO to NaN
        inverse_mask[np.nonzero(mask)]= int(
            1)  # ones inside poly, zero outside

        mask[np.nonzero(mask)]= float('nan')
        #TURN ALL ZERO to 1
        mask[np.where(mask == 0)]= int(1)  # nan inside poly, 1 outside
        region_mask= inverse_mask*(-1)+1
        #region_mask= np.nan_to_num(region_mask)  # zero in poly, 1 outside
        cutout[ymin:ymax, xmin:xmax]= mask[ymin:ymax, xmin:xmax]
        #TURN EVERYTHING OUTSIDE THAT RANGE to NaN
        cutout[np.where(cutout == 0)]= float('nan')

        #TAKE image=workask*mask will make a image with original values but NaN in polygon
        #blank= img*mask
        blank= img * cutout
        self.masked= mask
        self.blanked= blank
        goodvals= np.where(np.isfinite(blank))

        #perform the interpolation over the masked coordinates
        x= goodvals[1]  # x values of finite coordinates
        y= goodvals[0]  # y values of finite coordinates

        range_array= np.arange(x.size)
        fvals= np.zeros(x.size)
        for (i, xi, yi) in zip(range_array, x, y):
            fvals[i]= img[yi][xi]

        newfunc= interpolate.Rbf(
            x, y, fvals,
            function='multiquadric')  # the function that does interpolation
        allvals= np.where(img)  # whole region to interpolate over
        xnew= allvals[1]
        ynew= allvals[0]
        fnew= newfunc(xnew, ynew)

        #put the interpolated values back into a 2D array for display and other uses
        #def make_2D(fnew, xnew, ynew, img):
        fnew_2D= np.zeros(
            (int(xnew.size /
                 ((img.shape)[0])), int(ynew.size / ((img.shape)[1]))),
            dtype=float)
        #print(new_array)
        range_array2= np.arange(fnew.size)
        #print("rangearay:",range_array)
        for (i, x, y) in zip(range_array2, xnew, ynew):
            fnew_2D[y][x]= fnew[i]

        #fnew_2D= make_2D(fnew, xnew, ynew, img)
        interp= img * region_mask + fnew_2D * inverse_mask
        self.interp= interp

        #generate the residual image (original - interpolated background)
        self.resid= img - interp
        
######################################################
# Code Begins Here                                   #
######################################################

#Open the catalog file to get YB names, l, b, radius
data= ascii.read(catalog_name, delimiter=',')

#Set Pre-chosen range
BegYB= 1

YBlast= 6176

#Set linear threshold of the SymLogNorm
LinearThreshold= 0.001

#The below allows the user to start at the beginning of the range or from where they left off last time the program was ran

#The below allows the user to start at a specified YB number or from where they left off last time the program was ran
entry = plt.figure(figsize= (5,2.5))
entry.suptitle('Welcome! Where would you like to start?')
startbox_axes = plt.axes([.2,.6,.7,.25])
entrybox = TextBox(startbox_axes, '')
entrybox.on_submit(entrynumber)
startbutton_axes = plt.axes([.2,.2,.7,.25])
startbutton = Button(startbutton_axes, 'Start where you left off')
startbutton.on_clicked(leftoff)

startpos = 0
plt.show()
while startpos == 0:
    plt.pause(0.1)

YB1= int(startpos-1)
YB2= int(YBlast)
done= False
wavelengths= ['70', '24', '12', '8']

while (YB1 < YB2) and not done:
    #get the YB's location and radius
    YB= data[YB1]['YB']
    YB_long= data[YB1]['l']
    YB_lat= data[YB1]['b']
    YB_rad= data[YB1]['r']
    YBcoords= pd.read_csv(out_name, usecols= ['vertices 8', 'vertices 12', 'vertices 24', 'vertices 70'])
    
    #Use the location to determine the correct image files to use
    #image= choose_image(YB_long, YB_lat)
    
    ##          Unit Conversions Info         ##
    # 8 um:
    #GLIMPSE 8 um flux units are MJy/steradian
    #obtain the 8 um pixel scale from the image header
    #this gives the degree/pixel conversion factor used for overdrawing circle
    #and flux conversions (this is x direction only but they are equal in this case)
    #GWC edit 4/1/22: Set pixscale directly to 0.000333333
    #pixscale8= abs(image.um8w.wcs.cd[0][0])
    pixscale8= 0.000333333
    #pixscale8= abs(image.um8w.wcs.cdelt1)
    #print(pixscale8)
    #8 um square degree to square pixel conversion-- x*y
    #GWC edit 4/1/22: Set sqdeg_tosqpix8 directly to pixscale8 * pixscale8
    #sqdeg_tosqpix8= abs(image.um8w.wcs.cd[0][0]) * abs(
    #    image.um8w.wcs.cd[1][1])
    sqdeg_tosqpix8= pixscale8 * pixscale8
    #8 um steradian to pixel conversion (used to convert MJy/Pixel)
    #     will be used to convert to MJy
    str_to_pix8= sqdeg_tosqpix8 * 0.0003046118
    #WISE Units are Data Number per Pixel, Jy/DN is 1.8326*10^-6
    #See http://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec2_3f.html
    dn_to_jy12= 1.83 * 10**-6

    #convert YB l and b and radius to pixel coordinates
    ybwcs= np.array([[YB_long, YB_lat]], np.float_)
    #pixcoords= image.um8w.wcs_world2pix(ybwcs, 1)
    #YB_long_pix= pixcoords[0][0]
    #YB_lat_pix= pixcoords[0][1]
    YB_rad_pix= YB_rad / pixscale8
    dxw = 50
    dyw = 50
    YBloc= pd.read_csv(catalog_name, usecols= ['l', 'b', 'r'])
    '''
    # This is for the added points to show the user what other YBs are in the images
    # Read in the l, b, and r values for all the YBs and convert them to pixels
    YBloc= pd.read_csv(catalog_name, usecols= ['l', 'b', 'r'])
    YBcoords= pd.read_csv(out_name, usecols= ['vertices 8', 'vertices 12', 'vertices 24', 'vertices 70'])
    # Convert l, b, and r from YBloc into pixels
    for i in range(len(YBloc)):
        yblocwcs= np.array([[YBloc['l'][i], YBloc['b'][i]]], np.float_)
        pixcoordsloc= image.um8w.wcs_world2pix(yblocwcs, 1)
        YB_l= pixcoordsloc[0][0]
        YBloc['l'][i]= YB_l
        YB_b= pixcoordsloc[0][1]
        YBloc['b'][i]= YB_b
        YB_radloc= YBloc['r'][i]
        YB_r= YB_radloc / pixscale8
        YBloc['r'][i]= YB_r
    
    #define a window to zoom in on the YB
    xw= YB_long_pix + 0.5 # x coordinates of source and window center
    yw= YB_lat_pix + 0.5 # y coordinates of source and window center
    dxw= 50  # x width of window
    dyw= 50  # y width of window

    #find the pixel coordinates LLH and URH corner of zoomed window
    x1= int(xw - dxw)
    y1= int(yw - dyw)
    x2= int(xw + dxw)
    y2= int(yw + dyw)
    
    #use Cutout2D to make the zoomed windows
    position= (xw, yw)
    size= (2 * dxw, 2 * dyw)

    cut8= Cutout2D(data=image.um8data,
                    position=position,
                    size=size,
                    wcs=image.um8w)
    cut12= Cutout2D(data=image.um12data,
                     position=position,
                     size=size,
                     wcs=image.um12w)
    cut24= Cutout2D(data=image.um24data,
                     position=position,
                     size=size,
                     wcs=image.um24w)
    cut70= Cutout2D(data=image.um70data,
                     position=position,
                     size=size,
                     wcs=image.um70w)

    fitcopy8= image.um8
    fitcopy8.data= cut8.data
    fitcopy8.header.update(cut8.wcs.to_header())

    fitcopy12= image.um12
    fitcopy12.data= cut12.data
    fitcopy12.header.update(cut12.wcs.to_header())

    fitcopy24= image.um24
    fitcopy24.data= cut24.data
    fitcopy24.header.update(cut24.wcs.to_header())

    fitcopy70= image.um70
    fitcopy70.data= cut70.data
    fitcopy70.header.update(cut70.wcs.to_header())

    orig8= cut8.data
    orig12= cut12.data
    orig24= cut24.data
    orig70= cut70.data
    origlist= [orig70, orig24, orig12, orig8]

    wcslist= [cut70.wcs, cut24.wcs, cut12.wcs, cut8.wcs]

    #create empty residuals to fill up later
    diffslist= [orig70 * 0, orig24 * 0, orig12 * 0, orig8 * 0]

    #create copies of cropped images called workmasks
    workmasks= [copy.deepcopy(orig70), copy.deepcopy(orig24), copy.deepcopy(orig12), copy.deepcopy(orig8)]
    '''
    
    origlist = []
    avmasks = []
    interps= []
    resids= []
    wcslist= []
    heatmaps= []
    
    for i in range(len(wavelengths)):
        origpath = os.path.join(path, f"Flagging_cutouts/{wavelengths[i]}cropped_YB_{YB1+1}.fits")
        avmaskpath = os.path.join(path, f"Flagging_cutouts/{wavelengths[i]}masked_YB_{YB1+1}.fits")
        interppath = os.path.join(path, f"Flagging_cutouts/{wavelengths[i]}interp_YB_{YB1+1}.fits")
        residpath = os.path.join(path, f"Flagging_cutouts/{wavelengths[i]}resid_YB_{YB1+1}.fits")
        heatpath = os.path.join(path, f"Flagging_cutouts/{wavelengths[i]}heatmap_YB_{YB1+1}.fits")
        address= 'vertices ' + wavelengths[i]
        address2= wavelengths[i] + 'flag1'
        coordinates= YBcoords[address][YB1 + 1]
        #if coordinates available, do flags
        if coordinates != ' ':
            orig = fits.open(origpath)
            origlist.append(orig[0].data)
            avmask= fits.open(avmaskpath)
            avmasks.append(avmask[0].data)
            interp= fits.open(interppath)
            tempinterp = interp[0].data
            np.nan_to_num(tempinterp, copy=False, nan=0.0, posinf=None, neginf=None)
            interps.append(tempinterp)
            resid= fits.open(residpath)
            resids.append(resid[0].data)
            heatmap= fits.open(heatpath)
            heatmaps.append(heatmap[0].data)
            w = wcs.WCS(orig[0].header)
            wcslist.append(w)
            
        #check if saturated
        else:
            orig = fits.open(origpath)
            origlist.append(orig[0].data)
            avmasks.append('Saturated Image')
            interps.append('Saturated Image')
            resids.append('Saturated Image')
            heatmaps.append('Saturated Image')
            w = wcs.WCS(orig[0].header)
            wcslist.append(w)
    
    flag= make_flags(origlist, avmasks, resids, interps, heatmaps, wcslist, YB)

    #####################################################
    # Writing to CSV                                    #        
    #####################################################
    
    flagheaders= ['8flag1', '12flag1', '24flag1', '70flag1', 'flag2', '8flag3', '12flag3', '24flag3', '70flag3', '8flag4', '12flag4',
                  '24flag4', '70flag4', 'VCEflag', 'Review']
    
    #write flags to csv
    df= pd.read_csv(out_name)
    
    k= str(YB)
    
    for j in range(len(flag)):
        df.loc[df['YB'] == k, flagheaders[j]]= flag[j]
                
    df.to_csv(out_name, index= False)
    
    YB1 += 1
    pickle_out= open("latestYB.pickle", "wb")
    pickle.dump(YB1, pickle_out)
    pickle_out.close()
    orig.close()
    avmask.close()
    interp.close()
    resid.close()
    heatmap.close()
    
    
    if done:
        print('Goodbye! See you next time!')
        #Save current YB, so the user can pick up from here next time they run the program
        pickle_out= open("latestYB.pickle", "wb")
        pickle.dump(YB1, pickle_out)
        pickle_out.close()
        sys.exit()
        
plt.close('all')
print(
    'Congratulations! You have finished flagging all the YBs in your range!'
)
    
