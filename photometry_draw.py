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
#import itertools  
import sys
import math
import csv
#import pylab as py
import copy
import cv2
import os
import pandas as pd

#You will need to make sure the following packages have been installed first:

#from tkinter import *
#import tkinter
#conda install -c anaconda tk


from photutils import centroid_com
#https://photutils.readthedocs.io/en/stable/install.html
#conda install photutils -c astropy


#import threading
#import mynormalize
#import mycolorbar

"""
    Katie Devine/Anupa Poudyal
    Last Update: 3 October 2019
    
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
path='/Users/katiedevine/YB/code/'
#image_name = path+'GLM_03000+0000_mosaic_I4.fits'
catalog_name = path+'USE_THIS_CATALOG_ybcat_MWP_with_ID.csv' 
out_name = path+ 'YBphotometry_results.csv'
#######################################################
# Define my functions and classes                    # 
######################################################

#function that shows image and collects the coords for interp
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

        circle = Circle((dxw/2,dyw/2),YB_rad_pix, fill=False) 
        circle2 = Circle((dxw/2,dyw/2),YB_rad_pix, fill=False) 
        circle3 = Circle((dxw/2,dyw/2),YB_rad_pix, fill=False) 
        circle4 = Circle((dxw/2,dyw/2),YB_rad_pix, fill=False) 
        
        clkfig = plt.figure(figsize=(18,27))
        
        #Plot rg image,
        axrg = plt.subplot(2,3,1, title = 'RG (24 + 8 um)', projection = imgw)
        r = orig24
        g = orig
        b = np.zeros(orig.shape)
        axrg.imshow(make_lupton_rgb(r, g, b, stretch=200, Q=0), norm=LogNorm())
        axrg.add_artist(circle)   
    

        clkfig.add_subplot(2,3,2, title='Select coordinates in this '+wave+' image', projection = imgw)
        clkfig.suptitle('You are examining the '+wave+' image for YB' + str(ybid))
        plt.imshow(img, cmap = 'gray', norm = LogNorm())
        

        #Plot the 24 um
        foo = plt.subplot(2,3,4,title='24 um', projection = imgw)
        plt.imshow(orig24, cmap = 'hot', norm = LogNorm(vmin=orig24.min(), vmax=orig24.max()))
        foo.add_artist(circle2)  
        
        #Plot the 12 um
        faa = plt.subplot(2,3,5,title='12 um', projection = imgw)
        plt.imshow(orig12, cmap = 'hot', norm = LogNorm(vmin=orig12.min(), vmax=orig12.max()))
        faa.add_artist(circle3)  
    
        #Plot the 8um
        fum = plt.subplot(2,3,6,title = '8 um', projection = imgw)
        plt.imshow(orig, cmap = 'hot', norm = LogNorm(vmin=orig.min(), vmax=orig.max()))   
        fum.add_artist(circle4) 
#        cbar = plt.colorbar(format='%05.2f')
#        cbar.set_norm(mynormalize.MyNormalize(vmin=orig.min(),vmax=orig.max(),stretch='linear'))
#        cbar = mycolorbar.DraggableColorbar(cbar,orig)
#        cbar.connect()
        
        #Plot instructions
        plt.subplot(2,3,3)
        plt.axis([0,10,0,10])
        text = ("*Click in the grey-scale image. ")
        text1 = ("*Left click to add points.")
        text2 = ("*Middle click when finished to exit.")
        text3 = ("*All other images are for reference.")
        text4 = ("*Circles indicate MWP user radius.")
        text5 = ("*You will be prompted to inspect results;")
        text6 = ("type 'n' to continue, anything else to redo.")
        plt.text(1,8, text, ha='left', wrap=True)
        plt.text(1,7, text1, ha='left', wrap=True)
        plt.text(1,6, text2, ha='left', wrap=True)
        plt.text(1,5, text3, ha='left', wrap=True)        
        plt.text(1,4, text4, ha='left', wrap=True) 
        plt.text(1,3, text5, ha='left', wrap=True) 
        plt.text(1,2, text6, ha='left', wrap=True)         
        
        #interactive clicking to fill up coords
        coords=clkfig.ginput(n=-1, timeout=60, show_clicks=True, mouse_stop=2)


        #if no clicks, prompt to continue or exit
        if coords == []:
            answer = input("Do you wish to continue? Press y for yes. Else any other key, which will save nothing in your program. ")

            if answer == 'y':
               print ('You have chosen to continue on')

            else:
               print ("You have chosen to quit this program")
               sys.exit()
                   
        plt.close('clkfig')          
        return coords

#generates and saves the images
#call with (image, masked image, interpolated image, resid)
def make_figs(im1, im2, im3, im4, imw, um, cin):
    
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
    figurename=path+'photom_images/'+um+'interpolation_YB_%s.png' %(YB)
    fig.savefig(figurename)

#Function to examine images and input flags for output file
def make_flags(fim1, fim2, um):
    plt.figure(figsize=(6,3))

    plt.subplot(1,2,1,title='Original Data')
    plt.imshow(fim1, cmap = 'hot', norm = LogNorm(vmin=fim1.min(), vmax=fim1.max()))

    plt.subplot(1,2,2,title='Bkgrnd Removed')
    plt.imshow(fim2, cmap = 'hot')
    
    plt.pause(1)
    
    flag=[0,0,0,0,0,0,0,0,0,0,0]
            
    print('####################################')
    print('Do you want to note any flags in the output file?')
    print('Select the following flag(s) to apply')
    print('####################################')
          
    if um == '24' or um == '12':
        foo=0
        
        while foo != 5:
            flag1="Saturated Image"
            flag2="Diffraction Pattern/Star"
            flag3="Poor Confidence in Photometry"
            flag4="Other/Revisit this source"
            print('flag options:')
            print('[1] '+flag1)
            print('[2] '+flag2)
            print('[3] '+flag3)
            print('[4] '+flag4)
            print('[5] Done Flagging')
            print('[6] Clear Flags and Start Over')
            foo = int(input("select option:  "))
            if foo == 1:
                flag[0]=1
            if foo == 2:
                flag[1]=1                
            if foo == 3:
                flag[2]=1 
            if foo == 4:
                flag[3]=1 
            if foo == 6:
                flag==[0,0,0,0,0,0,0,0,0,0,0]  
                
    if um == '8':
        foo=0
        while foo != 9:
            flag1="Saturated Image"
            flag2="Multiple sources within YB"
            flag3="Filamentary or bubble rim structure"
            flag4="Not a YB or a star"
            flag5="IRDC Association"
            flag6="Diffraction Pattern/Star"
            flag7="Poor Confidence in Photometry"
            flag8="Other/Revisit this source"   
            
            print('flag options:')
            print('[1] '+flag1)
            print('[2] '+flag2)
            print('[3] '+flag3)
            print('[4] '+flag4)
            print('[5] '+flag5)
            print('[6] '+flag6)
            print('[7] '+flag7)
            print('[8] '+flag8)            
            print('[9] Done Flagging')
            print('[10] Clear Flags and Start Over')
            
            foo = int(input("select option:  "))
            if foo == 1:
                flag[0]=1
            if foo == 2:
                flag[1]=1   
            if foo == 3:
                flag[2]=1
            if foo == 4:
                flag[3]=1    
            if foo == 5:
                flag[4]=1
            if foo == 6:
                flag[5]=1     
            if foo == 7:
                flag[6]=1   
            if foo == 8:
                flag[7]=1                   
            if foo == 10:
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

        temp = fits.open(path8)
        self.um8 = temp
        self.um8data = temp[0].data
        self.um8w = wcs.WCS(temp[0].header) 
        temp = fits.open(path12)    
        self.um12 = temp
        self.um12data = temp[0].data
        self.um12w = wcs.WCS(temp[0].header) 
        temp = fits.open(path24)    
        self.um24 = temp
        self.um24data = temp[0].data
#        self.um24w = wcs.WCS(temp[0].header)


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
        
#class to calculate the "compactness index"
class compactness():
    def __init__(self, inresid, inrad):
        dif_flux = 0
        flux = 0
        
        #generate a rotated by 90 deg image and get flux difference in Jy
        rotated_image = [[] for x in range(len(inresid))]
        for i in range(len(inresid)):
            for j in range(len(inresid[i])):
                rotated_image[len(inresid) - j - 1].append(inresid[i][j])
        
        dif_im = rotated_image-inresid
        
        for ROW in range (0,100):       #sum the values of the pixels
            for column in range(0,100):
                dif_flux = dif_flux + abs(dif_im[ROW][column]) 
                flux = flux + abs(inresid[ROW][column]) 
                
        #normalize the flux difference by original input
        asym = abs(dif_flux)/abs(flux)
                    
        #find the offset between the user identified center (at 50,50) and maximum
        x, y = centroid_com(inresid)
        xdif = abs(x-50)
        ydif = abs(y-50)
        offset = math.sqrt(xdif**2+ydif**2)
        #normalize by radius and then increase by some scaling factor, currently 10
        offset = offset/inrad*10
        
        #use the information above to calucate the "compactness index" 
        cindex=math.sqrt((0.25*offset)**2+(0.75*asym)**2)
        
        self.index = cindex
        
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
        #output_file = open(out_name,append_write) #####opens up file for writing data   
        #writer=csv.DictWriter(output_file,fieldnames=['YB','YB_long','YB_lat','8umphotom','12umphotom','24umphotom'])
    
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

    #define a window to zoom in on the YB
    xw = xs = int(YB_long_pix + 0.5) # x coordinates of source and window center
    yw = ys = int(YB_lat_pix + 0.5)  # y coordinates of source and window center
    dxw = 100                        # x width of window
    dyw = 100                        # y width of window

    #find the pixel coordinates LLH and URH corner of zoomed window
    x1 = int(xw-0.5*dxw)
    y1 = int(yw-0.5*dyw)
    x2 = int(xw+0.5*dxw)
    y2 = int(yw+0.5*dyw)

    #Create cropped 100 x 100 pixel image arrays centered on YB
    orig = image.um8data[y1:y2,x1:x2]
    orig12 = image.um12data[y1:y2,x1:x2]
    orig24 = image.um24data[y1:y2,x1:x2]
    
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
    check='y'  
    while check != 'n':
        # 24 um image analysis
        #reset global list coords that gets creaeted in get_coords   
        print('######################################################')
        print('Beginning the 24 um analysis')
        coords=[]  
        #get the coordinates on 24um image
        coordinates=get_coords(workmask24, image.um8w, '24 um', YB)
        print('got coords')
        #do the masking and interpolation on 8um image
        print('starting interp')
        interp24 = do_interp(workmask24, coordinates)
        diff24=interp24.resid
        #calculate the compactness index
        compact = compactness(diff24, YB_rad_pix)
        compact24 = compact.index
        #display and save the images for the 8um image
        make_figs(workmask24, interp24.blanked, interp24.interp, interp24.resid, image.um8w, '24_um', compact24)

        #ACCEPT IMAGE OR REDO
        plt.pause(1)
        check = input('Please consult the residual image. Would you like to redo? Type n to continue, anything else to redo:  ')
    plt.close('all')
    #flag24 = make_flags(workmask24, interp24.resid, '24')
    flag24=[0,0,0,0,0,0,0,0,0,0,0]
    coord24=str(coordinates)
    plt.close('all')   
    
    check='y'
    while check != 'n':
        # 12 um image analysis
        #reset global list coords that gets creaeted in get_coords 
        print('######################################################')
        print('Beginning the 12 um analysis')    
        coords=[]  
        #get the coordinates on 24um image
        coordinates=get_coords(workmask12, image.um12w, '12 um', YB)
        #do the masking and interpolation on 8um image
        interp12 = do_interp(workmask12, coordinates)
        diff12=interp12.resid
        #calculate the compactness index
        compact = compactness(diff12, YB_rad_pix)
        compact12 = compact.index
        #display and save the images for the 8um image
        make_figs(workmask12, interp12.blanked, interp12.interp, interp12.resid, image.um12w, '12_um', compact12)

        #PROMPT TO ACCEPT IMAGE OR REDO 
        plt.pause(1) 
        check = input('Please consult the residual image. Would you like to redo? Type n to continue, anything else to redo:  ')
    plt.close('all')
#    flag12 = make_flags(workmask12, interp12.resid, '12')
    flag12=[0,0,0,0,0,0,0,0,0,0,0]
    coord12=str(coordinates)
    plt.close('all')  
    
    check='y'
    while check != 'n':       
        # 8 um image analysis
        #reset global list coords that gets creaeted in get_coords    
        print('######################################################')
        print('Beginning the 8 um analysis')
        coords=[]  
        #get the coordinates on 8um image
        coordinates=get_coords(workmask, image.um8w, '8 um', YB)
        #do the masking and interpolation on 8um image
        interp8 = do_interp(workmask, coordinates)
        diff8=interp8.resid
        #calculate the compactness index
        compact = compactness(diff8, YB_rad_pix)
        compact8 = compact.index
        #display and save the images for the 8um image
        make_figs(workmask, interp8.blanked, interp8.interp, interp8.resid, image.um8w, '8_um', compact8)

        #PROMPT TO ACCEPT IMAGE OR REDO 
        plt.pause(1)
        check = input('Please consult the residual image. Would you like to redo? Type n to continue, anything else to redo:  ')
    plt.close('all')
#    flag8 = make_flags(workmask, interp8.resid, '8')
    flag8=[0,0,0,0,0,0,0,0,0,0,0]
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
    
    df.loc[df["YB"]== kk , "8umphotom"] = round(flux_tot.um8,5)
    df.loc[df["YB"]== kk ,'12umphotom'] = round(flux_tot.um12,5)
    df.loc[df["YB"]== kk ,'24umphotom'] = round(flux_tot.um24,5)
    df.loc[df["YB"]== kk ,'8cin'] = round(compact8,4)
    df.loc[df["YB"]== kk ,'12cin'] = round(compact12,4)    
    df.loc[df["YB"]== kk ,'24cin'] = round(compact24,4)

    df.loc[df["YB"]== kk ,'vertices 8'] = coord8
    df.loc[df["YB"]== kk ,'vertices 12'] = coord12
    df.loc[df["YB"]== kk ,'vertices 24'] = coord24
    
    df.loc[df["YB"]== kk ,'24flag1'] = flag24[0]
    df.loc[df["YB"]== kk ,'24flag2'] = flag24[1] 
    df.loc[df["YB"]== kk ,'24flag3'] = flag24[2]
    df.loc[df["YB"]== kk ,'24flag4'] = flag24[3]
    
    df.loc[df["YB"]== kk ,'12flag1'] = flag12[0]
    df.loc[df["YB"]== kk ,'12flag2'] = flag12[1]
    df.loc[df["YB"]== kk ,'12flag3'] = flag12[2]
    df.loc[df["YB"]== kk ,'12flag4'] = flag12[3] 
    
    df.loc[df["YB"]== kk ,'8flag1'] = flag8[0]
    df.loc[df["YB"]== kk ,'8flag2'] = flag8[1] 
    df.loc[df["YB"]== kk ,'8flag3'] = flag8[2]
    df.loc[df["YB"]== kk ,'8flag4'] = flag8[3] 
    df.loc[df["YB"]== kk ,'8flag5'] = flag8[4]
    df.loc[df["YB"]== kk ,'8flag6'] = flag8[5] 
    df.loc[df["YB"]== kk ,'8flag7'] = flag8[6] 
    df.loc[df["YB"]== kk ,'8flag8'] = flag8[7] 
        
    df.to_csv(out_name, index=False)
    
plt.close('all')


