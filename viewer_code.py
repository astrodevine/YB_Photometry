#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 14:52:22 2020

@author: anupapoudyal
"""
import numpy as np
#import matplotlib
import matplotlib.pyplot as plt
import csv
import cv2
import os
from PIL import Image
from astropy.io import ascii
import shutil

path='/Users/anupapoudyal/Desktop/images/'
path1='/Users/anupapoudyal/Desktop/images/2020Summer/'
out_name = path1+ 'Second_YBphotometry.csv'
images = path1 + 'photom_images_2/'
compact = path1 + 'compact<2'

data = ascii.read(out_name, delimiter = ',')

#Open the output file and write the column headers

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
    Cin_8 = data[k]['8cin']
    Cin_12 = data[k]['12cin']
    Cin_24 = data[k]['24cin']
    
    um = '24'
    
    if Cin_24 < 2.0:
        #im = Image.open(path1+'photom_images_2/'+um+'_uminterpolation_YB_%s.png' %(YB))
        #im.show()
        imageNames = [path1+'photom_images_2/'+um+'_uminterpolation_YB_%s.png' %(YB)]
        for imageName in imageNames:
            shutil.copy(os.path.join(images, imageName), compact)
            
    um = '12'
    if Cin_12 < 2.0:
        #im = Image.open(path1+'photom_images_2/'+um+'_uminterpolation_YB_%s.png' %(YB))
        #im.show()
        imageNames = [path1+'photom_images_2/'+um+'_uminterpolation_YB_%s.png' %(YB)]
        for imageName in imageNames:
            shutil.copy(os.path.join(images, imageName), compact)
    
    um = '8'
    if Cin_8 < 2.0:
        imageNames = [path1+'photom_images_2/'+um+'_uminterpolation_YB_%s.png' %(YB)]
        #print(YB)
        for imageName in imageNames:
            shutil.copy(os.path.join(images, imageName), compact)
