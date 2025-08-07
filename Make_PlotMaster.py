# -*- coding: utf-8 -*-
"""
Created on Sat Jun 28 16:26:39 2025

@author: colemaya
"""

import os
import pandas as pd
import csv

path= '.'
xmat_name= os.path.join(path, 'MIRION-XMAT.csv')
catalog_name= os.path.join(path, 'YBphotometry_results_compiled_6-12-25.csv')
out_name= os.path.join(path, 'Plot Master.csv')
flag_name= os.path.join(path, 'YBphotometry_results_flagcsvfinal.csv')

headers = ['YB', 'F8', 'u_F8', 'e_F8', 'F12','u_F12', 'e_F12', 'F24','u_F24', 'e_F24', 
           'F70','u_F70','e_F70', 'Multiple Sources', 'No Obvious Source 70', 'No Obvious Source 24','No Obvious Source 12',
           'No Obvious Source 8','Poor Confidence 70', 'Poor Confidence 24','Poor Confidence 12',
           'Poor Confidence 8', 'Very Circular',
           'RMS', 'WISE Q', 'WISE C,G,K', 'CORNISH', 'No Assoc.']
data = pd.read_csv(flag_name).fillna('')
data2 = pd.read_csv(catalog_name).fillna('')
data3 = pd.read_csv(xmat_name).fillna('none')

flags = data.values.tolist()
obs = data2.values.tolist()

output_file = open(out_name, 'w')

writer = csv.DictWriter(
    output_file,
    fieldnames=headers,
    lineterminator='\n')
writer.writeheader()

for k in range(1,6177): ## THIS IS THE FULL CATALOG RANGE
    YB = data2['YB'][k]
    YB_8 = data2['8umaveflux'][k]
    YB_8u = data2['8umstdev'][k]
    YB_8e = data2['8umFE'][k]
    YB_12 = data2['12umaveflux'][k]
    YB_12u = data2['12umstdev'][k]
    YB_12e = data2['12umFE'][k]
    YB_24 = data2['24umaveflux'][k]
    YB_24u = data2['24umstdev'][k]
    YB_24e = data2['24umFE'][k]
    YB_70 = data2['70umaveflux'][k]
    YB_70u = data2['70umstdev'][k]
    YB_70e = data2['70umFE'][k]
    YB_Flag = []
    for i in range(10):
        YB_Flag.append(int(flags[k][i+15]))
    if data3['RMS ID'][k-1]=='none':
        rms = 0
    else:
        rms = 1
    if data3['WISE Cl'][k-1]=='none':
        wiseq = 0
        wisecgk = 0
    elif data3['WISE Cl'][k-1]=='Q':
        wiseq = 1
        wisecgk = 0
    else:
        wiseq = 0
        wisecgk = 1
    if data3['CORNISH ID'][k-1]=='none':
        cornish = 0
    else:
         cornish = 1
    if [rms, wiseq, wisecgk, cornish] == [0]*4:
        noassoc=1
    else:
        noassoc = 0
    writer.writerow({'YB': YB, 'F8': YB_8, 'u_F8': YB_8u, 'e_F8': YB_8e,
                     'F12': YB_12, 'u_F12': YB_12u,'e_F12': YB_12e, 'F24': YB_24,'u_F24': YB_24u, 'e_F24': YB_24e,
                     'F70': YB_70, 'u_F70': YB_70u,'e_F70': YB_70e, 'Multiple Sources': YB_Flag[0],
                     'No Obvious Source 8': YB_Flag[1], 'No Obvious Source 12': YB_Flag[2],
                     'No Obvious Source 24': YB_Flag[3], 'No Obvious Source 70': YB_Flag[4],
                     'Poor Confidence 8': YB_Flag[5], 'Poor Confidence 12': YB_Flag[6],
                     'Poor Confidence 24': YB_Flag[7], 'Poor Confidence 70': YB_Flag[8],
                     'Very Circular': YB_Flag[9], 'RMS': rms,
                     'WISE Q': wiseq, 'WISE C,G,K': wisecgk, 'CORNISH': cornish, 
                     'No Assoc.': noassoc})
output_file.close()

