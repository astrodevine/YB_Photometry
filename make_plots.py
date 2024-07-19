# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 09:52:58 2024

@author: colemaya
"""
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

plt.ion()
import os
import pandas as pd

# These lines supress warnings
import warnings
warnings.filterwarnings('ignore')


#Paths and File Locations#

#EDIT THIS PATH FOR THE FILE LOCATION ON YOUR MACHINE
path= '.'
catalog_name= os.path.join(path, 'YBXmatch.csv') #Crossmatch catalog csv name
instID= 'WolfChase1-ALL_YBS' #Change to be ID of csv you want data from
out_name= os.path.join(path, 'YBphotometry_results_' + instID + '.csv')

#Retrieve data from csv files

ybfluxdata= pd.read_csv(out_name, usecols= ['YB', '8umphotom','12umphotom','24umphotom', '70umphotom'])
ybxdata= pd.read_csv(catalog_name, usecols= ['YB', 'RMS ID', 'WISE Cl', 'CORNISH N ID', 'CORNISH S ID'])

yb8um= ybfluxdata['8umphotom']
yb12um= ybfluxdata['12umphotom']
yb24um= ybfluxdata['24umphotom']
yb70um= ybfluxdata['70umphotom']

colorvals12_8= []
colorvals24_8= []
colorvals70_24= []

#calculate color values and store them into lists #

for i in range(1, len(ybfluxdata)):
    color12_8= np.log10(float(yb12um[i])/float(yb8um[i]))
    color24_8= np.log10(float(yb24um[i])/float(yb8um[i]))
    color70_24= np.log10(float(yb70um[i])/float(yb24um[i]))
    
    colorvals12_8.append(color12_8)
    colorvals24_8.append(color24_8)
    colorvals70_24.append(color70_24)

#############################################
#Write data into lists based on association #
#############################################

data1 = []
data2 = []
data3 = []
data4 = []

numplottedlist = [0] * 6
withassoc= ybxdata['YB'].tolist()
RMSassoc= ybxdata['RMS ID']
WISEassoc= ybxdata['WISE Cl']
CORNSassoc= ybxdata['CORNISH S ID']
CORNNassoc= ybxdata['CORNISH N ID']

#these lists are to make the color-color point plot work.
noassocdata= []
RMSdata= []
WISEQdata= []
WISEelsedata= []
CORNSdata= []
CORNNdata= []
badsmogdata= []

for i in range(1, len(ybfluxdata)):
    #Delete this first if statement to include the SMOG data
    if 3035 <= i <= 3296:
        badsmogdata.append([colorvals70_24[i-1], colorvals24_8[i-1]])
        data1.append(colorvals70_24[i-1])
        data2.append(colorvals24_8[i-1])
        data3.append(colorvals12_8[i-1])
        data4.append('SMOG')
    else:
        data1.append(colorvals70_24[i-1])
        data2.append(colorvals24_8[i-1])
        data3.append(colorvals12_8[i-1])
        data4.append('Full Catalog')
        if str(colorvals70_24[i-1]) != 'nan' and str(colorvals24_8[i-1]) != 'nan':
            numplottedlist[0] += 1
        if i not in withassoc:
            noassocdata.append([colorvals70_24[i-1], colorvals24_8[i-1]])
            data1.append(colorvals70_24[i-1])
            data2.append(colorvals24_8[i-1])
            data3.append(colorvals12_8[i-1])
            data4.append('No Association')
            if str(colorvals70_24[i-1]) != 'nan' and str(colorvals24_8[i-1]) != 'nan':
                numplottedlist[5] += 1
        else:
            j= withassoc.index(i)
            if str(RMSassoc[j]) != 'nan':
                RMSdata.append([colorvals70_24[i-1], colorvals24_8[i-1]])
                data1.append(colorvals70_24[i-1])
                data2.append(colorvals24_8[i-1])
                data3.append(colorvals12_8[i-1])
                data4.append('RMS')
                if str(colorvals70_24[i-1]) != 'nan' and str(colorvals24_8[i-1]) != 'nan':
                    numplottedlist[1] += 1
            if WISEassoc[j] == 'Q':
                WISEQdata.append([colorvals70_24[i-1], colorvals24_8[i-1]])
                data1.append(colorvals70_24[i-1])
                data2.append(colorvals24_8[i-1])
                data3.append(colorvals12_8[i-1])
                data4.append('WISE Q')
                if str(colorvals70_24[i-1]) != 'nan' and str(colorvals24_8[i-1]) != 'nan':
                    numplottedlist[2] += 1
            if WISEassoc[j] == 'K' or WISEassoc[j] == 'G' or WISEassoc[j] == 'C':
                WISEelsedata.append([colorvals70_24[i-1], colorvals24_8[i-1]])
                data1.append(colorvals70_24[i-1])
                data2.append(colorvals24_8[i-1])
                data3.append(colorvals12_8[i-1])
                data4.append('WISE C,G,K')
                if str(colorvals70_24[i-1]) != 'nan' and str(colorvals24_8[i-1]) != 'nan':
                    numplottedlist[3] += 1
            if str(CORNSassoc[j]) != 'nan' or str(CORNNassoc[j]) != 'nan':
                CORNSdata.append([colorvals70_24[i-1], colorvals24_8[i-1]])
                data1.append(colorvals70_24[i-1])
                data2.append(colorvals24_8[i-1])
                data3.append(colorvals12_8[i-1])
                data4.append('CORNISH')
                if str(colorvals70_24[i-1]) != 'nan' and str(colorvals24_8[i-1]) != 'nan':
                    numplottedlist[4] += 1

######################
#make the histograms #  
######################
df = pd.DataFrame({'log(F70/F24)': data1, 'log(F24/F8)': data2, 'log(F12/F8)': data3, 'Association': data4})
datatitles= ['Full Catalog', 'RMS', 'WISE Q', 'WISE C,G,K', 'CORNISH', 'No Association']
colorsets= ['green','red', 'blue', 'orange', 'purple', 'gray']

fig= plt.figure(figsize=(30,30))
fig.canvas.header_visible= False
fig.suptitle('log(F12/F8) Color Histograms', fontsize = 25)
plt.axis('off')
plt.text(.5, -.075, 'log(F12/F8)', fontsize=20, ha='center')
plt.text(-.075, 0.5, 'Number', fontsize=20, rotation='vertical', va='center')

for i in range(len(datatitles)):
    data= df[(df == datatitles[i]).any(axis=1)]['log(F12/F8)'].tolist()
    vals12_8 = [x for x in data if str(x) != 'nan']
    minifig= plt.subplot(2,3,i+1)
    plt.title(datatitles[i], fontsize= 15)
    plt.hist(vals12_8, facecolor=colorsets[i], bins= 20)
    minifig.tick_params(labelsize=15)
    ave12_8= sum(vals12_8)/len(vals12_8)
    HIIline= plt.axvline(x=-0.09, color='black', linestyle='-.', label= 'HII Region Average')
    PAVEline= plt.axvline(x=-0.43, color='black', label= 'Pilot Region Average')
    FAVEline= plt.axvline(x=ave12_8, color='hotpink', label= datatitles[i] + ' Average')
    plt.legend(handles=[HIIline, PAVEline, FAVEline])

##################
# Make KDE plots #
##################
#Uncertainty variables
xvar= 0.1
yvar= 0.15

colors =  ['Greens', 'Reds', 'Blues', 'Oranges', 'Purples', 'Greys']

fig= plt.figure(figsize=(30,30))
fig.canvas.header_visible= False
fig.suptitle('YB log(F24/F8) v. log(F70/F24)  Color-Color Density', fontsize=25)
plt.axis('off')
plt.text(.5, -.075, 'log(F70/F24)', fontsize=20, ha='center')
plt.text(-.075, 0.5, 'log(F24/F8)', fontsize=20, rotation='vertical', va='center')

for i in range(len(datatitles)):
    minifig= plt.subplot(2,3,i + 1)
    plt.title(datatitles[i], fontsize=15)
    plt.xlim(-.75,2.5)
    plt.ylim(-1.25,2)
    sns.kdeplot(data= df[(df == datatitles[i]).any(axis=1)], x='log(F70/F24)', y='log(F24/F8)', 
                shade=True, cmap=colors[i], shade_lowest=False)
    minifig.set(xlabel=None, ylabel=None)
    minifig.tick_params(labelsize=15)
    plt.axhline(y=1.0,color='k', linestyle='dashed', label='')
    plt.axvline(x=0.8, color='k', linestyle='dashed', label='')
    #Average HII box
    plt.vlines(1.05, 0.26, 0.84, colors='k', linestyles='solid', label='')
    plt.vlines(1.47, 0.26, 0.84, colors='k', linestyles='solid', label='')
    plt.hlines(0.26, 1.05,1.47, colors='k', linestyles='solid', label='')
    plt.hlines(0.84, 1.05, 1.47, colors='k', linestyles='solid', label='')
    #Uncertainty crosshatch
    plt.hlines(-1, 2.25-xvar, 2.25+xvar, colors='k', linestyles='dotted')
    plt.vlines(2.25, -1-yvar, -1+yvar, colors='k', linestyles='dotted')
    plt.text(2,1.8,'N=' + str(numplottedlist[i]), fontsize=15)
    
################################################
#make the color-color plot (replaced with KDE) #
################################################
'''
We are using KDE plotting instead, but you can uncomment this section
and the code will still work if you want to see pt plotting of your color-color data
'''
'''
fig= plt.figure(figsize=(10,10))
fig.canvas.header_visible= False
plt.subplot(1,1,1,title=r'YB log(F24/F8) v. log(F70/F24)  Color-Color Plot')
plt.xlim(-.75,2.5)
plt.ylim(-1,2)

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

#plt.subplot(1,2,2,title=r'YB log(F24/F8) v. log(F70/F24) SMOG Color-Color Plot')
#for i in range(len(badsmogdata)):
#    ccdata= badsmogdata[i]
#    plt.plot(ccdata[0],ccdata[1], 'x', color= 'black')

#plt.xlabel('log(F70/F24)')
#plt.ylabel('log(F24/F8)')
#plt.axhline(y=1.0,color='k', linestyle='dashed', label='')
#plt.axvline(x=0.8, color='k', linestyle='dashed', label='')
#
#plt.vlines(1.05, 0.26, 0.84, colors='k', linestyles='solid', label='')
#plt.vlines(1.47, 0.26, 0.84, colors='k', linestyles='solid', label='')
#
#plt.hlines(0.26, 1.05,1.47, colors='k', linestyles='solid', label='')
#plt.hlines(0.84, 1.05, 1.47, colors='k', linestyles='solid', label='')'''
