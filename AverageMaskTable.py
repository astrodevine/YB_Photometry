from google_auth_oauthlib.flow import InstalledAppFlow
from googleapiclient.discovery import build
from googleapiclient.http import MediaIoBaseDownload
import io
import traceback 
import ast 

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
from matplotlib.path import Path
from skimage.measure import find_contours
from ast import literal_eval

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
import os
import math
import pandas as pd
from astropy.nddata import Cutout2D

# These lines supress warnings
import warnings
warnings.filterwarnings('ignore')

import io
import pandas as pd
from google_auth_oauthlib.flow import InstalledAppFlow
from googleapiclient.discovery import build
from googleapiclient.http import MediaIoBaseDownload

# Authenticate
SCOPES = ['https://www.googleapis.com/auth/drive.readonly']

flow = InstalledAppFlow.from_client_secrets_file(
    'client_secret_527671004811-mov0s2s96a9gqsdim7h0qtl3od2dnuos.apps.googleusercontent.com.json',
    SCOPES
)

creds = flow.run_local_server(port=0)
service = build('drive', 'v3', credentials=creds)

# Folder ID
folder_id = "1E9w5Irzj0FP9VUbCVYHM3PQtuMmLZtJU2FQTkaqCcNnv_ZANHuTlFqbayGdb0v9jZYQQwz3C"

# Only return files whose MIME type is text/csv
query = (
    f"'{folder_id}' in parents "
    f"and mimeType='text/csv'"
)

results = service.files().list(
    q=query,
    fields="files(id,name)"
).execute()

files = results.get("files", [])

print(f"Found {len(files)} CSV files")

all_data = []

for file in files:

    file_id = file["id"]
    filename = file["name"]

    print(f"Loading: {filename}")

    request = service.files().get_media(fileId=file_id)

    fh = io.BytesIO()
    downloader = MediaIoBaseDownload(fh, request)

    done = False
    while not done:
        status, done = downloader.next_chunk()
        if status:
            print(f"{status.progress() * 100:.1f}%")

    fh.seek(0)

    df = pd.read_csv(fh, header=0, skiprows=[1])

    # Remove units row if present
    if len(df) > 0 and str(df.iloc[0]["YB"]) == "ID Number":
        df = df.iloc[1:].reset_index(drop=True)

    all_data.append(df)

print("\nDONE LOADING ALL CSV FILES")


SCOPES = ['https://www.googleapis.com/auth/drive.readonly']

flow = InstalledAppFlow.from_client_secrets_file(
    'client_secret_527671004811-mov0s2s96a9gqsdim7h0qtl3od2dnuos.apps.googleusercontent.com.json',
    SCOPES
)

creds = flow.run_local_server(port=0)

service = build('drive', 'v3', credentials=creds)

# File ID of the CSV
file_id = "1aNrUvi-GaKBhbK3k0QWrJzryypa2s34FQDtr3cuu9OY"

# Download the file into memory
request = service.files().export_media(
    fileId=file_id,
    mimeType="text/csv"
)
fh = io.BytesIO()

downloader = MediaIoBaseDownload(fh, request)

done = False
while not done:
    status, done = downloader.next_chunk()
    if status:
        print(f"{status.progress() * 100:.1f}%")

# Read with pandas
fh.seek(0)
control = pd.read_csv(fh, skiprows=[1])

def polygon_to_mask(vertices, shape=(100, 100)):
    y, x = np.mgrid[:shape[0], :shape[1]]
    points = np.column_stack((x.ravel(), y.ravel()))
    path = Path(vertices)
    mask = path.contains_points(points)
    return mask.reshape(shape)

yb_ids = control["YB"]
all_rows = [[
    "ID Number", "degree", "degree", "pixel coords", "pixel coords",
    "pixel coords", "pixel coords", "Jy",
    "Saturated", "Multiple sources within YB", "Filament or Bubble Rim",
    "No obvious source at this wavelength", "IRDC Association",
    "Star/Diffraction Pattern", "Poor Confidence", "Other/Follow Up",
    "Jy", "Saturated", "Multiple sources within YB",
    "No obvious source at this wavelength", "Star/Diffraction Pattern",
    "Poor Confidence", "Other/Follow Up",
    "Jy", "Saturated", "Multiple sources within YB",
    "No obvious source at this wavelength", "Star/Diffraction Pattern",
    "Poor Confidence", "Other/Follow Up",
    "Jy", "Saturated", "Multiple sources within YB",
    "No obvious source at this wavelength", "Star/Diffraction Pattern",
    "Poor Confidence", "Other/Follow Up"
]]
yb_head = control.head().columns.tolist()
lengths = ["vertices 8","vertices 12","vertices 24","vertices 70"]
umlst = ["8umphotom","12umphotom","24umphotom","70umphotom"]
flag_cols = ["8flag1","8flag2","8flag3","8flag4","8flag5","8flag6","8flag7","8flag8", '12flag1', '12flag2', '12flag4', '12flag6', '12flag7', '12flag8', '24flag1', '24flag2', '24flag4', '24flag6', '24flag7', '24flag8', '70flag1', '70flag2', '70flag4', '70flag6', '70flag7', '70flag8']
old_flag_cols = ["8flag1","8flag2","8flag3","8flag4","8flag5","8flag6","8flag7","8flag8", '12flag1', '12flag2', '12flag4', '12flag6', '12flag7', '12flag8', '24flag1', '24flag2', '24flag4', '24flag6', '24flag7', '24flag8']

for yb in yb_ids:
    um8_tot = um12_tot = um24_tot = um70_tot = 0
    um8ave = um12ave = um24ave = um70ave = 0
    count8 = count12 = count24 = count70 = 0
    count1 = 0
    count_70 = 0
    avgs_vertices = []
    for csv in all_data:
        row = csv[csv["YB"] == yb].iloc[0]
        if row["8umphotom"] != "Saturated" and not pd.isna(row["8umphotom"]):
            um8_tot += float(row["8umphotom"])
            count8 += 1
        if row["12umphotom"] != "Saturated" and not pd.isna(row["12umphotom"]):
            um12_tot += float(row["12umphotom"])
            count12 += 1
        if row["24umphotom"] != "Saturated" and not pd.isna(row["24umphotom"]):
            um24_tot += float(row["24umphotom"])
            count24 += 1
        if "70umphotom" in row.index:
            if row["70umphotom"] != "Saturated" and not pd.isna(row["70umphotom"]):
                um70_tot += float(row["70umphotom"])
                count70 += 1 
    for length in lengths:
        shape_lst = []
        masks= []
        for csv in all_data:
            if length not in csv.columns:
                continue
           # row = csv.iloc[yb-1]
            row = csv[csv["YB"] == yb].iloc[0]
           # if row[length]
            verts_str = row[length]
            
            if pd.isna(verts_str):
                continue
            try:
                verts = np.array(literal_eval(verts_str))
                shape_lst.append(verts)
            except:
                continue
        if len(shape_lst) == 0:
            avgs_vertices.append("")
            continue

        shape=(100,100)

        for verti in shape_lst:
            mask = polygon_to_mask(verti,shape)
            masks.append(mask)

        masks = np.array(masks)
            
        test = masks.mean(axis=0)
        average_mask = test >= .5

        contours = find_contours(average_mask.astype(float), level=.5)
        if len(contours) > 0:
            avg_vertices = contours[0]
            avg_vertices_str = str([(float(x), float(y))for x, y in avg_vertices])
            avgs_vertices.append(avg_vertices_str)
            count1 = 1
    count = count8+count12+count24+count70
    if count > 0:
        if count8 != 0:
            um8ave = round((um8_tot / count8), 4)
        if count12 != 0:
            um12ave = round((um12_tot / count12), 4)
        if count24 != 0:
            um24ave = round((um24_tot / count24), 4)
        if count70 != 0:
            um70ave = round((um70_tot / count70), 4)
    else:
        um8ave = um12ave = um24ave = um70ave = 0

    for head in yb_head[3:]:
        if head in flag_cols:
            avgs_vertices.append("")
        elif count > 0:
            if head == "8umphotom":
                if count8 != 0:
                    avgs_vertices.append(str(um8ave))
                else:
                    avgs_vertices.append("")
            elif head == "12umphotom":
                if count12 != 0:
                    avgs_vertices.append(str(um12ave))
                else:
                    avgs_vertices.append("")
            elif head == "24umphotom":
                if count24 != 0:
                    avgs_vertices.append(str(um24ave))
                else:
                    avgs_vertices.append("")
            elif head == "70umphotom":
                if count70 != 0:
                    avgs_vertices.append(str(um70ave))
                else:
                    avgs_vertices.append("")

    big_row = [int(yb),float(control['YB_long'].iloc[int(yb)-1]),float(control['YB_lat'].iloc[int(yb)-1])]+avgs_vertices           
    all_rows.append(big_row)
    print(yb)

connect = pd.DataFrame(all_rows, columns = yb_head)
connect.to_csv(f"/Users/wadevining/YellowBall/YB_CSV_Sorting/MasterTable.csv", index=False)    
    