#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 13:24:30 2020

@author: anupapoudyal
"""
import numpy as np
#import matplotlib
import matplotlib.pyplot as plt
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
#from PIL import Image
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
#You will need to make sure the following packages have been installed first:
import ast
from photutils import centroid_com
#https://photutils.readthedocs.io/en/stable/install.html
#conda install photutils -c astropy


#import threading
#import mynormalize
#import mycolorbar

"""
    Katie Devine/Anupa Poudyal
    
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
    
#######################################################
# Paths and File Locations                           #
######################################################

#EDIT THIS PATH FOR THE FILE LOCATION ON YOUR MACHINE
path='/Users/anupapoudyal/Desktop/images/'
path1='/Users/anupapoudyal/Desktop/images/2020Summer/'
#image_name = path+'GLM_03000+0000_mosaic_I4.fits'
catalog_name = path+'USE_THIS_CATALOG_ybcat_MWP_with_ID.csv' 
out_name = path1+ 'YBphotometry_results.csv'
#######################################################
# Define my functions and classes                    # 
######################################################

#function that shows image and collects the coords for interp
def get_coords(wave, ybid):
        
        df= pd.read_csv(out_name)
        
        print (ybid)
        
        yb = ybid - 1152
        ybid = str(ybid)
        
        row = df[df['YB']== ybid]
        colname = 'vertices '+ wave
        
    
        
        coord = row[colname]
        
        if (coord[yb]) != (coord[yb]):
            coords = []
        
        else:
            coordss= ast.literal_eval(coord[yb]) 
        
        #converting string to list
        #coordss= ast.literal_eval(d[0].append(1))
        
            coords = list(coordss)
        
        return coords

#generates and saves the images
#call with (image, masked image, interpolated image, resid)
def make_figs(im1, im2, im3, im4, fitfile, imw, um, cin):
    
############Generate the figures for each source##################
#note I'm being lazy here and calling from the code things that aren't defined in function
#this is pretty ugly and should maybe get cleaned up 
    fig = plt.figure(figsize=(8, 8))
    
    #Plot the original image and MWP User YB circle

    circle = Circle((dxw/2,dyw/2),YB_rad_pix, fill=False) 
    fig1=plt.subplot(2,2,1,title='Cropped image',projection = imw)
    plt.imshow(im1, cmap = 'hot', norm = LogNorm(vmin=im1.min(), vmax=im1.max()))
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    fig1.add_artist(circle)    
    fig1.text(dxw/2,dyw/2-5, 'MWP size', color='white', ha='center', va='top', weight='bold')
   
    #Plot the mask 
    plt.subplot(2,2,2,title='Masked Image', projection = imw)
    plt.imshow(im2, cmap = 'hot', norm = LogNorm(vmin=im1.min(), vmax=im1.max()))
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    
    #Plot the interpolated background
    plt.subplot(2,2,3,title='Interpolated image', projection = imw)
    plt.imshow(im3, cmap = 'hot', norm = LogNorm(vmin=im1.min(), vmax=im1.max()))
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')

    #Plot the residual
    text='Compactness Index: ' +str(round(cin, 4))
    plt.subplot(2,2,4,title = 'Residual(Image-Interp)', projection = imw)
    plt.imshow(im4, cmap = 'hot')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.text(10,10, text, color='white', ha='left')
#Make the plots pretty 
    plt.tight_layout(pad=0.2, w_pad=5, h_pad=6)    
    plt.subplots_adjust(top=0.9)
    fig.suptitle('Interpolation of YB %s at '%(YB) + um , fontsize=10)
    
    # Save this result as a new png
    #figurename=path1+'compactimages/'+um+'interpolation_YB_%s.png' %(YB)
    #fig.savefig(figurename)


    # Save the fits cut-outs for fututre use if needed
    #im1name=path1+'compact_fits/'+um+'cropped_YB_%s.fits' %(YB)
    #im2name=path1+'compact_fits/'+um+'masked_YB_%s.fits' %(YB)
    #im3name=path1+'compact_fits/'+um+'interp_YB_%s.fits' %(YB)
    #im4name=path1+'compact_fits/'+um+'resid_YB_%s.fits' %(YB)

    #fitfile.data=im1 
    #fitfile.writeto(im1name, overwrite=True)

    #fitfile.data=im2
    #fitfile.writeto(im2name, overwrite=True)

    #fitfile.data=im3 
    #fitfile.writeto(im3name, overwrite=True)

    #fitfile.data=im4
    #fitfile.writeto(im4name, overwrite=True)

    
#Function to examine images and input flags for output file
def make_flags(fim1, fim2, um):
    plt.figure(figsize=(6,3))

    plt.subplot(1,2,1,title='Original Data')
    plt.imshow(fim1, cmap = 'hot', norm = LogNorm(vmin=fim1.min(), vmax=fim1.max()))

    plt.subplot(1,2,2,title='Bkgrnd Removed')
    plt.imshow(fim2, cmap = 'hot')
    
    plt.pause(1)
    
    flag=[0,0,0,0,0,0,0,0,0,0,0]
            
                               
    return flag

#######################################################
# Define my classes              # 
######################################################
#class choose_image: use l range to choose and open right mosaics
# returns .um8, .um8data, .um8w - full info, data array, wcs header at 8 um
# returns .um12, .um12data, .um12w - full info, data array, wcs header at 12 um
# returns .um24, .um24data, .um124w - full info, data array, wcs header at 24 um
class choose_image():        
    def __init__(self, l):
        #currently don't need the WCS files for 12, 24 um because they
        #are reprojections onto 8um coordinate grid
        
        if l >= 28.5 and l<=31.5:
            path8 = path+'mosaics/GLM_03000+0000_mosaic_I4.fits'
            path12 = path+'mosaics/WISE_12um_03000_mosaic.fits'
            path24 = path+'mosaics/MIPSGAL_03000_mosaic_reprojected.fits'
        elif l > 31.5 and l <=34.5:
            path8 = path+'mosaics/GLM_03300+0000_mosaic_I4.fits'
            path12 = path+'mosaics/WISE_12um_03300_mosaic.fits'
            path24 = path+'mosaics/MIPSGAL_03300_mosaic_reprojected.fits'
        elif l > 34.5 and l <=37.5:
            path8 = path+'mosaics/GLM_03600+0000_mosaic_I4.fits'
            path12 = path+'mosaics/WISE_12um_03600_mosaic.fits'
            path24 = path+'mosaics/MIPSGAL_03600_mosaic_reprojected.fits'
        elif l > 37.5 and l <=40.5:
            path8 = path+'mosaics/GLM_03900+0000_mosaic_I4.fits'
            path12 = path+'mosaics/WISE_12um_03900_mosaic.fits'
            path24 = path+'mosaics/MIPSGAL_03900_mosaic_reprojected.fits'
        else:
            print('Your YB is outside the pilot region l=30-40 degrees.')
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


#class that does the masking and interpolation, returns masked, blanked, interp, resid
class do_interp():
        def __init__(self, img, verts):
                       
            #use the clicked values from the user to create a NaN mask 
            
            vertices = np.array([verts], dtype=np.int32)  
            xyvals = np.array(verts, dtype=np.int32)  
            xmin = min(xyvals[:,0])-5
            xmax = max(xyvals[:,0])+5
            ymin = min(xyvals[:,1])-5
            ymax = max(xyvals[:,1])+5
            #print(xmin, xmax, ymin, ymax)
            mask = np.zeros_like(img)
            inverse_mask = np.zeros_like(img)
            region_mask = np.zeros_like(img)
            cutout = np.zeros_like(img)
            
            # filling pixels inside the polygon defined by "vertices" with the fill color
            cv2.fillPoly(mask, vertices, 255)          
            #TURN ALL Non-ZERO to NaN
            inverse_mask[np.nonzero(mask)] = int(1)  # ones inside poly, zero outside  

            mask[np.nonzero(mask)]=float('nan') 
            #TURN ALL ZERO to 1
            mask[np.where(mask == 0)] = int(1) # nan inside poly, 1 outside
            region_mask=mask
            region_mask=np.nan_to_num(region_mask) # zero in poly, 1 outside
            cutout[ymin:ymax,xmin:xmax]=mask[ymin:ymax,xmin:xmax]
            #TURN EVERYTHING OUTSIDE THAT RANGE to NaN
            cutout[np.where(cutout == 0)] = float('nan')
             
            #TAKE image=workask*mask will make a image with original values but NaN in polygon
            #blank = img*mask  
            blank = img*cutout    
            self.masked = mask
            self.blanked = blank
            goodvals = np.where(np.isfinite(blank))

            #perform the interpolation over the masked coordinates
            x = goodvals[1] # x values of finite coordinates
            y = goodvals[0]# y values of finite coordinates


            for i in x:
                for j in y:
        
                    def get_fvals(x,y):
                        range_array = np.arange(x.size)
                        vals = np.zeros(x.size)
                        for (i, xi, yi) in zip(range_array, x, y):
                             vals[i] = img[yi][xi]
                        return vals
                    
            fvals = get_fvals(x,y)
  
            newfunc = interpolate.Rbf(x, y, fvals, function='multiquadric') # the function that does interpolation
            allvals = np.where(img)  # whole region to interpolate over
            xnew = allvals[1]
            ynew = allvals[0]
            fnew = newfunc(xnew, ynew)

        #put the interpolated values back into a 2D array for display and other uses
            def make_2D(fnew, xnew, ynew, img):
                new_array = np.zeros((int(xnew.size/((img.shape)[0])), int(ynew.size/((img.shape)[1]))), dtype = float)
                #print(new_array)
                range_array = np.arange(fnew.size)
                #print("rangearay:",range_array)
        
                for (i, x, y) in zip(range_array, xnew, ynew):
                    new_array[y][x] = fnew[i]
        
                return new_array

            fnew_2D = make_2D(fnew, xnew, ynew, img)
            
            self.interp = img*region_mask+fnew_2D*inverse_mask
            
            #generate the residual image (original - interpolated background)
            self.resid = img - (img*region_mask+fnew_2D*inverse_mask)
            
            #turn the interpolated image into a FITS file   
            #hdu_new = fits.PrimaryHDU(fnew_2D)  

#class that gets the flux from residual images. Unit conversions are applied here.
class get_flux():
    def __init__(self, d8, d12, d24):
        #reset the total flux value
        flux_tot8=0
        flux_tot12=0            
        flux_tot24=0        
        #go through 100 x 100 pixel residual image and add up the pixel values
        for ROW in range (0,100):       
            for column in range(0,100):
                    
                flux_tot8= flux_tot8+ d8[ROW][column] #gives the value of photometry
                flux_tot12= flux_tot12+ d12[ROW][column] #gives the value of photometry
                flux_tot24= flux_tot24+ d24[ROW][column] #gives the value of photometry 
                
            #convert units of flux total.  MIPS/IRAC in MJy/Sr*Pixel, want Jy
            #conversion: (MJy/Sr*Pixel)*(Sr/Pix)*(Jy/MJy)
            #WISE is in units of data number per pixel
            #Conversion: DN/Pix*Pixel*Jy/DN
            
        flux_tot8=flux_tot8*str_to_pix8*10**6
        flux_tot12=flux_tot12*dn_to_jy12
        flux_tot24=flux_tot24*str_to_pix8*10**6               
 
        self.um8 = flux_tot8
        self.um12 = flux_tot12            
        self.um24 = flux_tot24    

#def crop_center(pil_img, crop_width, crop_height):
        #img_width, img_height = pil_img.shape
        #return pil_img.crop(((int(img_width - crop_width) // 2),
                             #(int(img_height - crop_height) // 2),
                             #(int(img_width + crop_width) // 2),
                             #int((img_height + crop_height) // 2)))
#class to calculate the "compactness index"
class compactness():
    def __init__(self, inresid, inrad,fitfile):
        dif_flux = 0
        flux = 0
        
        #new = inresid[(100-50)/2): ((100-50)/2), ((100+50)/2):((100+50)/2)]
        #cropped_img = inresid.crop(((100-50)//2, (100-50)//2, (100+50)//2, (100+50)//2))
        x, y = centroid_com(inresid)
        xdif = abs(x-50)
        ydif = abs(y-50)
        offset = math.sqrt(xdif**2+ydif**2)
        
        #height, width = inresid.shape
        #print(img.shape)
        #wi=(width/2)
        #he=(height/2)
        #cropped_img = inresid.crop((wi-50)//2, (he-50)//2, (wi+50)//2, (he+50)//2)

        fitfile.data=inresid
        
        position1=(x,y)
        size1=(51,51)  #Using odd sixe to use the middle pixel as the center
        
        cut = Cutout2D(data=inresid.data, position=position1, size=size1)
        cuto = cut.data
        
        #generate a rotated by 90 deg image and get flux difference in Jy
        rotated_image = [[] for x in range(len(cuto))]
        for i in range(len(cuto)):
            for j in range(len(cuto[i])):
                rotated_image[len(cuto) - j - 1].append(cuto[i][j])
        
        dif_im = rotated_image-cuto
        
        for ROW in range (0,51):       #sum the values of the pixels
            for column in range(0,51):
                dif_flux = dif_flux + abs(dif_im[ROW][column]) 
                flux = flux + abs(inresid[ROW][column]) 
                
        #normalize the flux difference by original input
        asym = abs(dif_flux)/abs(flux)
                    
        #find the offset between the user identified center (at 50,50) and maximum
        #x, y = centroid_com(inresid)
        #xdif = abs(x-50)
        #ydif = abs(y-50)
        #offset = math.sqrt(xdif**2+ydif**2)
        #normalize by radius and then increase by some scaling factor, currently 10
        
        offset = offset/inrad*10
        
        #use the information above to calucate the "compactness index" 
        cindex=math.sqrt((0.25*offset)**2+(0.75*asym)**2)
        
        self.index = cindex
 #############
        #plt.figure(figsize=(6,3))

        #plt.subplot(1,3,1,title='Resid')
       # plt.plot(25.5, 25.5, "og", markersize=5,marker = "x")
        #plt.imshow(cuto, cmap = 'hot')#, norm = LogNorm(vmin=rotated_image.min(), vmax=rotated_image.max()))
        
        #plt.subplot(1,3,2,title='Rotated')
        #plt.plot(25.5, 25.5, "og", markersize=5, marker = "x")
        #plt.imshow(rotated_image, cmap = 'hot')
        
        #plt.subplot(1,3,3,title='diff')
        #plt.imshow(dif_im, cmap = 'hot')
        
        #plt.pause(1)
 #############
        #figurename=path1+'compactness_rotated/'+'interpolation_YB_%s.png' %(YB)
        #plt.savefig(figurename)
#######################################################
# Actual Program Begins Here   #
######################################################

#Open the catalog file to get YB names, l, b, radius
data = ascii.read(catalog_name, delimiter = ',')

#Open the output file and write the column headers

if os.path.exists(out_name):
    append_write = 'a'# append if already exists
    output_file = open(out_name,append_write) #####opens up files for creating csv file 
    headers = ['YB', 'YB_long', 'YB_lat','vertices 8','vertices 12','vertices 24','8umphotom','8cin', '8flag1', '8flag2', 
               '8flag3', '8flag4','8flag5', '8flag6','8flag7','8flag8','12umphotom','12cin', '12flag1',
               '12flag2','12flag3','12flag4','24umphotom','24cin','24flag1', '24flag2', '24flag3', '24flag4']
    writer = csv.DictWriter(output_file,fieldnames=headers)
    output_file.close()
    
else:
    output_file = open(out_name,'w') #####opens up files for creating csv file 
    writer=csv.DictWriter(output_file,fieldnames=['YB','YB_long','YB_lat','vertices','vertices 8','vertices 12','vertices 24','8umphotom','8cin', '8flag1', 
                                                  '8flag2','8flag3', '8flag4','8flag5', '8flag6','8flag7', '8flag8',
                                                  '12umphotom','12cin', '12flag1', '12flag2','12flag3', '12flag4',
                                                  '24umphotom','24cin','24flag1','24flag2','24flag3',
                                                  '24flag4'], lineterminator = '\n')
    writer.writeheader()
    writer.writerow({'YB':'ID Number','YB_long':'degree','YB_lat':'degree','vertices 8':'pixel coords','vertices 12':'pixel coords','vertices 24':'pixel coords','8umphotom':'Jy', 
                     '12umphotom':'Jy', '24umphotom':'Jy', '24flag1':'Saturated', 
                     '24flag2':'Star/Diffraction Pattern', '24flag3':'Poor Confidence', '24flag4':'Other/Follow Up', 
                     '12flag1':'Saturated', '12flag2':'Star/Diffraction Pattern', '12flag3':'Poor Confidence',
                     '12flag4':'Other/Follow Up','8flag1':'Saturated', 
                     '8flag2':'Multiple sources within YB','8flag3':'Filament or Bubble Rim',
                     '8flag4':'Not a YB or Star','8flag5':'IRDC Association',
                     '8flag6':'Star/Diffraction Pattern', '8flag7':'Poor Confidence','8flag8':'Other/Follow Up'})
       
    for k in range (1152,1671):
        YB=data[k]['YB']
        YB_long=data[k]['l']
        YB_lat=data[k]['b']
        writer.writerow({'YB':YB,'YB_long':YB_long,'YB_lat':YB_lat})
    
    output_file.close()

#######################################################
# Begin the loop through YBs in the catalog (k=YB)   #
######################################################
    
YBfirst=input("Which YB do you want to start with: ")
YBlast=input('Which YB do you want to end loop with: ')

YB1=int(YBfirst)-1
YB2=int(YBlast)
for k in range (YB1,YB2):    
    #get the YB's location and radius
    
    YB = data[k]['YB']
    YB_long = data[k]['l']
    YB_lat = data[k]['b']
    YB_rad = data[k]['r']
    
    #Use the location to determine the correct image files to use
    image=choose_image(YB_long)
    
    ##          Unit Conversions Info         ##
    # 8 um:
    #GLIMPSE 8 um flux units are MJy/steradian
    #obtain the 8 um pixel scale from the image header 
    #this gives the degree/pixel conversion factor used for overdrawing circle
    #and flux conversions (this is x direction only but they are equal in this case)
    pixscale8=abs(image.um8w.wcs.cd[0][0])
    #8 um square degree to square pixel conversion-- x*y 
    sqdeg_tosqpix8=abs(image.um8w.wcs.cd[0][0])*abs(image.um8w.wcs.cd[1][1])
    #8 um steradian to pixel conversion (used to convert MJy/Pixel)
    #     will be used to convert to MJy
    str_to_pix8=sqdeg_tosqpix8*0.0003046118
    #WISE Units are Data Number per Pixel, Jy/DN is 1.8326*10^-6
    #See http://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec2_3f.html
    dn_to_jy12=1.83*10**-6
    
    #convert YB l and b and radius to pixel coordinates
    ybwcs = np.array([[YB_long,YB_lat]], np.float_)
    pixcoords = image.um8w.wcs_world2pix(ybwcs,1)
    YB_long_pix = pixcoords[0][0]
    YB_lat_pix = pixcoords[0][1]
    YB_rad_pix = YB_rad/pixscale8

#    #define a window to zoom in on the YB
    xw = xs = int(YB_long_pix + 0.5) # x coordinates of source and window center
    yw = ys = int(YB_lat_pix + 0.5)  # y coordinates of source and window center
    dxw = 100                        # x width of window
    dyw = 100                        # y width of window

    #find the pixel coordinates LLH and URH corner of zoomed window
    x1 = int(xw-0.5*dxw)
    y1 = int(yw-0.5*dyw)
    x2 = int(xw+0.5*dxw)
    y2 = int(yw+0.5*dyw)

#    #Create cropped 100 x 100 pixel image arrays centered on YB
#    orig = image.um8data[y1:y2,x1:x2]
#    orig12 = image.um12data[y1:y2,x1:x2]
#    orig24 = image.um24data[y1:y2,x1:x2]
#    
    
    #use Cutout2D to make the zoomed windows
    position=(YB_long_pix+0.5, YB_lat_pix+0.5)
    size=(100,100)
    
    cut8=Cutout2D(data=image.um8data, position=position, size=size, wcs=image.um8w)
    cut12=Cutout2D(data=image.um12data, position=position, size=size, wcs=image.um12w)
    cut24=Cutout2D(data=image.um24data, position=position, size=size, wcs=image.um24w)
  
    fitcopy8=image.um8
    fitcopy8.data=cut8.data
    fitcopy8.header.update(cut8.wcs.to_header())  
    
    fitcopy12=image.um12
    fitcopy12.data=cut12.data
    fitcopy12.header.update(cut12.wcs.to_header())  

    fitcopy24=image.um24
    fitcopy24.data=cut24.data
    fitcopy24.header.update(cut24.wcs.to_header())  
    
    orig=cut8.data
    orig12=cut12.data
    orig24=cut24.data
    
    wcs8=cut8.wcs
    wcs12=cut12.wcs
    wcs24=cut24.wcs
    
    #create empty residuals to fill up later
    diff8 = orig*0
    diff12 = orig12*0
    diff24 = orig24*0
    
    #create copies of cropped images called workmasks
    workmask = copy.deepcopy(orig)
    workmask12 = copy.deepcopy(orig12) 
    workmask24 = copy.deepcopy(orig24) 
    
    
    

    ###################################################################
    # Call the classes to draw polygons and perform interp #
    ###################################################################
    #check='y'  
    #while check != 'n':
        
        # 24 um image analysis
        #reset global list coords that gets creaeted in get_coords   
    print('######################################################')
    print('Beginning the 24 um analysis')
        
    coords=[]  
        #get the coordinates on 24um image
        #print(coords)
        
    coordinates= get_coords('24', YB)
    

        #if no clicks, exit
        
        #if coordinates == []:
            #print("Timeout limit of 5 minutes reached. Exiting program. Results for most recent YB will not be saved.")
            #sys.exit()
        
        #do the masking and interpolation on 8um image
    if coordinates != []:
        interp24 = do_interp(workmask24, coordinates)
        
        diff24=interp24.resid
        #calculate the compactness index
        compact = compactness(diff24, YB_rad_pix,fitcopy24)
        compact24 = compact.index
        #display and save the images for the 8um image
        make_figs(workmask24, interp24.blanked, interp24.interp, interp24.resid, fitcopy24, wcs24, '24_um', compact24)
    else:
        compact24 = 0
        #ACCEPT IMAGE OR REDO
    plt.pause(1)
    plt.close('all')
    #flag24 = make_flags(workmask24, interp24.resid, '24')
    #flag24=[0,0,0,0,0,0,0,0,0,0,0]
    coord24=str(coordinates)
    plt.close('all')   
    
    
        # 12 um image analysis
        #reset global list coords that gets creaeted in get_coords 
    print('######################################################')
    print('Beginning the 12 um analysis')    
    coords=[]  
        #get the coordinates on 24um image
    coordinates=get_coords('12', YB)
    if coordinates != []:           
        #do the masking and interpolation on 8um image
        interp12 = do_interp(workmask12, coordinates)
        diff12=interp12.resid
        #calculate the compactness index
        compact = compactness(diff12, YB_rad_pix,fitcopy12)
        compact12 = compact.index
        #display and save the images for the 8um image
        make_figs(workmask12, interp12.blanked, interp12.interp, interp12.resid, fitcopy12, wcs12, '12_um', compact12)
    else:
        compact12 = 0
        #PROMPT TO ACCEPT IMAGE OR REDO 
    plt.pause(1) 
    #flag12 = make_flags(workmask12, interp12.resid, '12')
    #flag12=[0,0,0,0,0,0,0,0,0,0,0]
    coord12=str(coordinates)
    plt.close('all')  
          
        # 8 um image analysis
        #reset global list coords that gets creaeted in get_coords    
    print('######################################################')
    print('Beginning the 8 um analysis')
    coords=[]  
        #get the coordinates on 8um image
    coordinates=get_coords('8', YB)
               
        #do the masking and interpolation on 8um image
    if coordinates != []:
        interp8 = do_interp(workmask, coordinates)
        diff8=interp8.resid
        #calculate the compactness index
        compact = compactness(diff8, YB_rad_pix,fitcopy8)
        compact8 = compact.index
        #display and save the images for the 8um image
        make_figs(workmask, interp8.blanked, interp8.interp, interp8.resid, fitcopy8, wcs8, '8_um', compact8)
        #PROMPT TO ACCEPT IMAGE OR REDO 
    else:
        compact8 = 0
    plt.pause(1)
    plt.close('all')
    #flag8 = make_flags(workmask, interp8.resid, '8')
    #flag8=[0,0,0,0,0,0,0,0,0,0,0]
    coord8=str(coordinates)
    plt.close('all')  
    
    
    
    ###################################################################
    # Use residual images to perform photometry and write out to table with flags#
    ###################################################################

    #call the get_flux class
    flux_tot=get_flux(diff8, diff12, diff24)
    
    #'8umphotom':flux_tot.um8,'12umphotom':flux_tot.um12,'24umphotom':flux_tot.um24}               
    
    df = pd.read_csv(out_name)
    
    kk=str(YB)
    
    #df.loc[df["YB"]== kk , "8umphotom"] = round(flux_tot.um8,5)
    #df.loc[df["YB"]== kk ,'12umphotom'] = round(flux_tot.um12,5)
    #df.loc[df["YB"]== kk ,'24umphotom'] = round(flux_tot.um24,5)
    
    df.loc[df["YB"]== kk ,'8cin'] = round(compact8,4)
    df.loc[df["YB"]== kk ,'12cin'] = round(compact12,4)    
    df.loc[df["YB"]== kk ,'24cin'] = round(compact24,4)

    #df.loc[df["YB"]== kk ,'vertices 8'] = coord8
    #df.loc[df["YB"]== kk ,'vertices 12'] = coord12
    #df.loc[df["YB"]== kk ,'vertices 24'] = coord24
    
    #df.loc[df["YB"]== kk ,'24flag1'] = flag24[0]
    #df.loc[df["YB"]== kk ,'24flag2'] = flag24[1] 
    #df.loc[df["YB"]== kk ,'24flag3'] = flag24[2]
    #df.loc[df["YB"]== kk ,'24flag4'] = flag24[3]
    
    #df.loc[df["YB"]== kk ,'12flag1'] = flag12[0]
    #df.loc[df["YB"]== kk ,'12flag2'] = flag12[1]
    #df.loc[df["YB"]== kk ,'12flag3'] = flag12[2]
    #df.loc[df["YB"]== kk ,'12flag4'] = flag12[3] 
    
    #df.loc[df["YB"]== kk ,'8flag1'] = flag8[0]
    #df.loc[df["YB"]== kk ,'8flag2'] = flag8[1] 
    #df.loc[df["YB"]== kk ,'8flag3'] = flag8[2]
    #df.loc[df["YB"]== kk ,'8flag4'] = flag8[3] 
    #df.loc[df["YB"]== kk ,'8flag5'] = flag8[4]
    #df.loc[df["YB"]== kk ,'8flag6'] = flag8[5] 
    #df.loc[df["YB"]== kk ,'8flag7'] = flag8[6] 
    #df.loc[df["YB"]== kk ,'8flag8'] = flag8[7] 
        
    df.to_csv(out_name, index=False)
    
plt.close('all')



