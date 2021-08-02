import numpy as np
import glob
import geo_QY
import time
import pdb
import os, sys
import pickle
import shutil
import netCDF4 as nc4
from datetime import datetime

import json

start_t = time.time()
print ('current dir: ', os.getcwd())
if True:
  dataDir2='./'
  dataDir4='./'
    
  if True:
    nc_files = sorted(glob.glob(dataDir2+'*.nc'))
    ### print ('nc_files: ', nc_files)
    print ('len(nc_files): ', len(nc_files))

    t_size = 0
    for file1 in nc_files:
      print(os.path.getsize(file1))
      t_size += os.path.getsize(file1)

    print('t_size (GB): ', t_size/1000000000.0)
