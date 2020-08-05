#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 16:07:43 2020

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
import ast
from photutils import centroid_com
from itertools import chain, repeat
import shutil
from PIL import Image

#EDIT THIS PATH FOR THE FILE LOCATION ON YOUR MACHINE
path='/Users/anupapoudyal/Desktop/images/'
path1='/Users/anupapoudyal/Desktop/images/2020Summer/'
#image_name = path+'GLM_03000+0000_mosaic_I4.fits'
catalog_name = path1 +'YBcompactness.csv' 
out_name = path1+ 'YBcompactness.csv'

#######################################################
# Define my functions and classes                    # 
######################################################

#function that shows image and collects the coords for interp

#Function to examine images and input flags for output file
def make_flags(um):
    
    flag=[0,0,0,0,0,0,0,0,0,0,0]
            
    print('####################################')
    print('Do you want to note any flags in the output file?')
    print('Select the following flag(s) to apply')
    print('####################################')
          
    if um == '24' or um == '12':
       foo = 0
       
       while foo != 9:
           
            flag1="Saturated Image"
            flag2="Diffraction Pattern/Star"
            flag3="Poor Confidence in Photometry"
            flag4="Other/Revisit this source"
            print('flag options:')
            print('[1] '+flag1)
            print('[2] '+flag2)
            print('[3] '+flag3)
            print('[4] '+flag4)
            print('[9] Done Flagging')
            print('[10] Clear Flags and Start Over')
            
                
            prompts = chain(["Enter a number from the flagging options:"], repeat("Not a flagging option! Try again:"))
            replies = map(input, prompts)
            numeric_strings = filter(str.isnumeric, replies)
            numbers = map(float, numeric_strings)
            is_positive = (0).__lt__
            valid_response = next(filter(is_positive, numbers))
            foo = valid_response
            
            if foo == 1:
                flag[0]=1
            if foo == 2:
                flag[1]=1                
            if foo == 3:
                flag[2]=1 
            if foo == 4:
                flag[3]=1 
            if foo == 10:
                flag==[0,0,0,0,0,0,0,0,0,0,0] 
            if foo == 9:
                print ("done flagging")
            else:
                foo == 0
        
                
    if um == '8':
        foo=0
        
        while foo != 9:
            flag1="Saturated Image"
            flag2="Multiple sources within masked area"
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
            
            prompts = chain(["Enter a number from the flagging options:"], repeat("Not a flagging option! Try again:"))
            replies = map(input, prompts)
            numeric_strings = filter(str.isnumeric, replies)
            numbers = map(float, numeric_strings)
            is_option = (0).__lt__  #This provides the condition thst the input should be greater than 0.
            valid_response = next(filter(is_option, numbers)) 
            foo = valid_response
            
            
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
            if foo == 9:
                print ("done flagging")
            else:
                foo == 0
    return flag

data = ascii.read(out_name, delimiter = ',')

if os.path.exists(out_name):
    append_write = 'a'# append if already exists
    output_file = open(out_name,append_write) #####opens up files for creating csv file 
    headers = ['YB', 'YB_long', 'YB_lat','vertices 8','vertices 12','vertices 24','8umphotom','8cin', '8flag1', '8flag2', 
               '8flag3', '8flag4','8flag5', '8flag6','8flag7','8flag8','12umphotom','12cin', '12flag1',
               '12flag2','12flag3','12flag4','24umphotom','24cin','24flag1', '24flag2', '24flag3', '24flag4']
    writer = csv.DictWriter(output_file,fieldnames=headers)
    output_file.close()
    
        
YBfirst=input("Which YB do you want to start with: ")
YBlast=input('Which YB do you want to end loop with: ')

YB1=int(YBfirst)-1152
YB2=int(YBlast)-1151

for k in range (YB1,YB2):    
    #get the YB's location and radius
    
    YB = data[k]['YB']
    
    um = '24'
    im24 = Image.open(path1+'photom_images_2/'+um+'_uminterpolation_YB_%s.png' %(YB))
    im24.show()
    flag24 = make_flags(um)
    
    
    um = '12'
    im12 = Image.open(path1+'photom_images_2/'+um+'_uminterpolation_YB_%s.png' %(YB))
    im12.show()
    flag12 = make_flags(um)
    
    um = '8'
    im8 = Image.open(path1+'photom_images_2/'+um+'_uminterpolation_YB_%s.png' %(YB))
    im8.show()
    flag8 = make_flags(um)
    im8.close()

    
    df = pd.read_csv(out_name)
    
    kk=str(YB)
    
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
