# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 09:52:58 2024

@author: colemaya
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

plt.ion()
import os
import pandas as pd

# These lines supress warnings
import warnings
warnings.filterwarnings('ignore')


#Paths and File Locations#

#EDIT THIS PATH FOR THE FILE LOCATION ON YOUR MACHINE
path= '.'
catalog_name= os.path.join(path, 'YBXmatch.csv')
instID= 'colemaya' #Change to be your ID
out_name= os.path.join(path, 'YBphotometry_results_' + instID + '.csv')

#Retrieve data from csv files

ybfluxdata= pd.read_csv(out_name, usecols= ['YB', '8umphotom','12umphotom','24umphotom', '70umphotom'])
ybxdata= pd.read_csv(catalog_name, usecols= ['YB', 'RMS ID', 'WISE Cl', 'CORNISH N ID', 'CORNISH S ID'])

yb8um= ybfluxdata['8umphotom']
yb12um= ybfluxdata['12umphotom']
yb24um= ybfluxdata['24umphotom']
yb70um= ybfluxdata['70umphotom']

numrejects12_8= 0
numrejects24_8= 0
numrejects70_24= 0

colorvals12_8= []
colorvals24_8= []
colorvals70_24= []

#calculate color values and store them into lists #

for i in range(1, len(ybfluxdata)):
    color12_8= np.log10(float(yb12um[i])/float(yb8um[i]))
    color24_8= np.log10(float(yb24um[i])/float(yb8um[i]))
    color70_24= np.log10(float(yb70um[i])/float(yb24um[i]))
    
    if str(color12_8) == 'nan':
        numrejects12_8 += 1
    colorvals12_8.append(color12_8)
        
    if str(color24_8) == 'nan':
        numrejects24_8 += 1
    colorvals24_8.append(color24_8)
        
    if str(color70_24) == 'nan':
        numrejects70_24 += 1
    colorvals70_24.append(color70_24)

#############################################
#Write data into lists based on association #
#############################################

withassoc= ybxdata['YB'].tolist()
RMSassoc= ybxdata['RMS ID']
WISEassoc= ybxdata['WISE Cl']
CORNSassoc= ybxdata['CORNISH S ID']
CORNNassoc= ybxdata['CORNISH N ID']

noassocdata= []
noassocdata12_8 = []

RMSdata= []
RMSdata12_8 = []

WISEQdata= []
WISEQdata12_8= []

WISEelsedata= []
WISEelsedata12_8 = []

CORNSdata= []
CORNSdata12_8= []

CORNNdata= []
CORNNdata12_8= []

for i in range(1, len(ybfluxdata)):
    if i not in withassoc:
        noassocdata.append([colorvals70_24[i-1], colorvals24_8[i-1]])
        if str(colorvals12_8[i-1]) != 'nan':
            noassocdata12_8.append(colorvals12_8[i-1])
    else:
        j= withassoc.index(i)
        if str(RMSassoc[j]) != 'nan':
            RMSdata.append([colorvals70_24[i-1], colorvals24_8[i-1]])
            if str(colorvals12_8[i-1]) != 'nan':
                RMSdata12_8.append(colorvals12_8[i-1])
        if WISEassoc[j] == 'Q':
            WISEQdata.append([colorvals70_24[i-1], colorvals24_8[i-1]])
            if str(colorvals12_8[i-1]) != 'nan':
                WISEQdata12_8.append(colorvals12_8[i-1])
        if WISEassoc[j] == 'K' or WISEassoc[j] == 'G' or WISEassoc[j] == 'C':
            WISEelsedata.append([colorvals70_24[i-1], colorvals24_8[i-1]])
            if str(colorvals12_8[i-1]) != 'nan':
                WISEelsedata12_8.append(colorvals12_8[i-1])
        if str(CORNSassoc[j]) != 'nan':
            CORNSdata.append([colorvals70_24[i-1], colorvals24_8[i-1]])
            if str(colorvals12_8[i-1]) != 'nan':
                CORNSdata12_8.append(colorvals12_8[i-1])
        if str(CORNNassoc[j]) != 'nan':
            CORNNdata.append([colorvals70_24[i-1], colorvals24_8[i-1]])
            if str(colorvals12_8[i-1]) != 'nan':
                CORNNdata12_8.append(colorvals12_8[i-1])

######################
#make the histograms #  
######################

colorvals12_8 = [x for x in colorvals12_8 if str(x) != 'nan']
colorvals12_8.sort()
ave12_8= sum(colorvals12_8)/len(colorvals12_8)
 
fig= plt.figure(figsize=(32,32))
fig.canvas.header_visible= False
plt.subplot(2,4,1, title=r'log(F12/F8) Color Histogram, all data')
plt.hist(colorvals12_8[2:-2], facecolor='dodgerblue', bins= 20)
plt.xlabel('log(F12/F8)')
plt.ylabel('Number')
HIIline= plt.axvline(x=-0.09, color='black', linestyle='-.', label= 'HII Region Average')
PAVEline= plt.axvline(x=-0.43, color='black', label= 'Pilot Region Average')
FAVEline= plt.axvline(x=ave12_8, color='hotpink', label= 'Full Catalog Average')
plt.legend(handles=[HIIline, PAVEline, FAVEline])

plt.subplot(2,4,2,title=r'log(F12/F8) Color Histogram, RMS Association')
plt.hist(RMSdata12_8, facecolor='red', bins= 20)
plt.xlabel('log(F12/F8)')
plt.ylabel('Number')
HIIline= plt.axvline(x=-0.09, color='black', linestyle='-.', label= 'HII Region Average')
PAVEline= plt.axvline(x=-0.43, color='black', label= 'Pilot Region Average')
FAVEline= plt.axvline(x=sum(RMSdata12_8)/len(RMSdata12_8), color='hotpink', label= 'RMS Average')
plt.legend(handles=[HIIline, PAVEline, FAVEline])

plt.subplot(2,4,3,title=r'log(F12/F8) Color Histogram, WISE Q Association')
plt.hist(WISEQdata12_8, facecolor='orange', bins= 20)
plt.xlabel('log(F12/F8)')
plt.ylabel('Number')
HIIline= plt.axvline(x=-0.09, color='black', linestyle='-.', label= 'HII Region Average')
PAVEline= plt.axvline(x=-0.43, color='black', label= 'Pilot Region Average')
FAVEline= plt.axvline(x=sum(WISEQdata12_8)/len(WISEQdata12_8), color='hotpink', label= 'WISE Q Average')
plt.legend(handles=[HIIline, PAVEline, FAVEline])

plt.subplot(2,4,4,title=r'log(F12/F8) Color Histogram, WISE C,G,K Association')
plt.hist(WISEelsedata12_8, facecolor='yellow', bins= 20)
plt.xlabel('log(F12/F8)')
plt.ylabel('Number')
HIIline= plt.axvline(x=-0.09, color='black', linestyle='-.', label= 'HII Region Average')
PAVEline= plt.axvline(x=-0.43, color='black', label= 'Pilot Region Average')
FAVEline= plt.axvline(x=sum(WISEelsedata12_8)/len(WISEelsedata12_8), color='hotpink', label= 'WISE C,G,K Average')
plt.legend(handles=[HIIline, PAVEline, FAVEline])

plt.subplot(2,4,5,title=r'log(F12/F8) Color Histogram, CORNISH S Association')
plt.hist(CORNSdata12_8, facecolor='green', bins= 20)
plt.xlabel('log(F12/F8)')
plt.ylabel('Number')
HIIline= plt.axvline(x=-0.09, color='black', linestyle='-.', label= 'HII Region Average')
PAVEline= plt.axvline(x=-0.43, color='black', label= 'Pilot Region Average')
FAVEline= plt.axvline(x=sum(CORNSdata12_8)/len(CORNSdata12_8), color='hotpink', label= 'CORNISH S Average')
plt.legend(handles=[HIIline, PAVEline, FAVEline])

plt.subplot(2,4,6,title=r'log(F12/F8) Color Histogram, CORNISH N Association')
plt.hist(CORNNdata12_8, facecolor='blue')
plt.xlabel('log(F12/F8)')
plt.ylabel('Number')
HIIline= plt.axvline(x=-0.09, color='black', linestyle='-.', label= 'HII Region Average')
PAVEline= plt.axvline(x=-0.43, color='black', label= 'Pilot Region Average')
FAVEline= plt.axvline(x=sum(CORNNdata12_8)/len(CORNNdata12_8), color='hotpink', label= 'CORNISH N Average')
plt.legend(handles=[HIIline, PAVEline, FAVEline])

plt.subplot(2,4,7,title=r'log(F12/F8) Color Histogram, No Association Data')
plt.hist(noassocdata12_8, facecolor='purple', bins=20)
plt.xlabel('log(F12/F8)')
plt.ylabel('Number')
HIIline= plt.axvline(x=-0.09, color='black', linestyle='-.', label= 'HII Region Average')
PAVEline= plt.axvline(x=-0.43, color='black', label= 'Pilot Region Average')
FAVEline= plt.axvline(x=sum(noassocdata12_8)/len(noassocdata12_8), color='hotpink', label= 'No Association Average')
plt.legend(handles=[HIIline, PAVEline, FAVEline])

############################
#make the color-color plot #
############################

fig= plt.figure(figsize=(10,10))
fig.canvas.header_visible= False
plt.subplot(1,1,1,title=r'YB log(F24/F8) v. log(F70/F24)  Color-Color Plot')
#plt.xlim(0,4.5)
#plt.ylim(-1.2,1.2)

for i in range(len(noassocdata)):
    ccdata= noassocdata[i]
    plt.plot(ccdata[0],ccdata[1], 'x', color= 'black')

for i in range(len(RMSdata)):
    ccdata= RMSdata[i]
    plt.plot(ccdata[0],ccdata[1], 's', fillstyle= 'none', color= 'red')
    
for i in range(len(WISEQdata)):
    ccdata= WISEQdata[i]
    plt.plot(ccdata[0],ccdata[1], '+', color= 'blue')
    
for i in range(len(WISEelsedata)):
    ccdata= WISEelsedata[i]
    plt.plot(ccdata[0],ccdata[1], '*', color= 'red')

for i in range(len(CORNSdata)):
    ccdata= CORNSdata[i]
    plt.plot(ccdata[0],ccdata[1], '^', fillstyle= 'none', color= 'purple')

for i in range(len(CORNNdata)):
    ccdata= CORNNdata[i]
    plt.plot(ccdata[0],ccdata[1], 'o', fillstyle= 'none', color= 'green')

plt.xlabel('log(F70/F24)')
plt.ylabel('log(F24/F8)')
plt.axhline(y=1.0,color='k', linestyle='dashed', label='')
plt.axvline(x=0.8, color='k', linestyle='dashed', label='')

plt.vlines(1.05, 0.26, 0.84, colors='k', linestyles='solid', label='')
plt.vlines(1.47, 0.26, 0.84, colors='k', linestyles='solid', label='')

plt.hlines(0.26, 1.05,1.47, colors='k', linestyles='solid', label='')
plt.hlines(0.84, 1.05, 1.47, colors='k', linestyles='solid', label='')

WQ = mlines.Line2D([], [], color='blue', marker='+', markersize=5, linestyle='none', label='WISE Q')
WALL = mlines.Line2D([], [], color='red', marker='*', markersize=5, linestyle='none', label='WISE K,G,C')
RRMS = mlines.Line2D([], [], color='red', marker='s', markersize=5, fillstyle='none',linestyle='none', label='RMS')
NO = mlines.Line2D([], [], color='black', marker='x', markersize=5,linestyle='none', label='No Assoc.')
CORNS = mlines.Line2D([], [], color='purple', marker='^', markersize=5, fillstyle = 'none', linestyle='none', label='CORNISH S')
CORNN = mlines.Line2D([], [], color='green', marker='o', markersize=5, fillstyle = 'none', linestyle='none', label='CORNISH N')

plt.legend(handles=[WQ, WALL, RRMS, NO, CORNS, CORNN])
