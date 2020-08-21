import numpy as np
import glob
import geo_QY
import time
import pdb
import os
import pickle
import netCDF4 as nc4

start_time = time.time()
"""
dataDir1='/peate_archive/.data5/Ops/npp/noaa/op/2012/05/15/scris/'
dataDir2='/peate_archive/.data5/Ops/npp/noaa/op/2012/05/15/gcrso/'
dataDir3='/peate_archive/.data5/Ops/npp/noaa/op/2012/05/15/svm15/'
dataDir4='/peate_archive/.data5/Ops/npp/noaa/op/2012/05/15/gmodo/'
"""
### for iday in range(15,23,1):
for iday in range(15,16,1):
#dataDir2='/peate_archive/.data6/Ops/snpp/gdisc/2/2015/06/01/crisl1b/'
    dataDir2='/peate_archive/.data1/Ops/snpp/gdisc/2/2015/01/'+str(iday).zfill(2)+'/crisl1b/'
    ### dataDir4='/raid15/qyue/VIIRS/VIIRS/201501/'
    dataDir4='/raid15/qyue/VIIRS/VIIRS/201501/VNP03MOD/'
    
    ### for iloop in range(0,239,10):
    for iloop in range(0,9,10):
        print(iloop)   
# get CrIS files 
#cris_sdr_files = sorted(glob.glob(dataDir1+'SCRIS*d2012*'))[21:40]
        cris_geo_files = sorted(glob.glob(dataDir2+'SNDR.SNPP.CRIS*'))[iloop:iloop+10]
        print ('cris_geo_files: ', cris_geo_files)
# get VIIRS files 
#viirs_sdr_files = sorted(glob.glob(dataDir3+'SVM15*d2012*'))[31:59]
        viirs_geo_files = sorted(glob.glob(dataDir4+'VNP03MOD*A2015'+str(iday).zfill(3)+'*'))[iloop:iloop+10]
        print ('viirs_geo_files: ', viirs_geo_files)
# read VIIRS data 
        viirs_lon, viirs_lat, viirs_satAzimuth, viirs_satRange, viirs_satZenith, viirs_height,viirs_time = geo_QY.read_nasa_viirs_geo(viirs_geo_files)
#viirs_bt, viirs_rad, viirs_sdrQa = geo.read_viirs_sdr(viirs_sdr_files)i


# read CrIS data 
        cris_lon, cris_lat, cris_satAzimuth, cris_satRange, cris_satZenith, cris_time, cris_realLW = geo_QY.read_nasa_cris_geo(cris_geo_files)

#cris_realLW = geo.read_nasa_cris_sdr(cris_sdr_files , sdrFlag=True)

# compute CrIS Pos Vector in EFEC on the Earth Surface 
        cris_pos= np.zeros(np.append(cris_lat.shape, 3))
        cris_pos[:, :, :, 0], cris_pos[:, :, :, 1], cris_pos[:, :, :, 2] \
	    = geo_QY.LLA2ECEF(cris_lon, cris_lat, np.zeros_like(cris_lat))

# compute CrIS LOS Vector in ECEF 
        cris_east, cris_north, cris_up = geo_QY.RAE2ENU(cris_satAzimuth, cris_satZenith, cris_satRange)

        cris_los= np.zeros(np.append(cris_lat.shape, 3))
        cris_los[:, :, :, 0], cris_los[:, :, :, 1], cris_los[:, :, :, 2] = \
	    geo_QY.ENU2ECEF(cris_east, cris_north, cris_up, cris_lon, cris_lat)

# compute viirs POS vector in ECEF
        viirs_pos= np.zeros(np.append(viirs_lat.shape, 3))
        viirs_pos[:, :, 0], viirs_pos[:, :, 1], viirs_pos[:, :, 2] = \
	    geo_QY.LLA2ECEF(viirs_lon, viirs_lat, np.zeros_like(viirs_lat))

# cris_los is pointing from pixel to satellite, we need to
#   change from satellite to pixel
        cris_los = -1.0*cris_los

# using Kd-tree to find the closted pixel of VIIRS for each CrIS FOV
# Set fake viirs_sdrQa to be zero: good quality everywhere since not for calibration
#viirs_sdrQa=np.zeros(viirs_lon.shape)

#remove the sdrqa, but adding time requirement (less than 600S difference)
        dy, dx = geo_QY.match_cris_viirs_QY(cris_los, cris_pos, viirs_pos, cris_time, viirs_time)
#	print("collocation are done in --- %s seconds --- for %d files " % (time.time() - start_time, len(cris_geo_files)))

        dy_flatten = np.array([item for lst in dy.reshape(-1) for item in lst])
        dy_size = np.array([len(lst) for lst in dy.reshape(-1)]).reshape(dy.shape)
        dx_flatten = np.array([item for lst in dx.reshape(-1) for item in lst])

        ### f = nc4.Dataset('/raid15/qyue/VIIRS/VIIRS/201501/Index/IND_CrIS_VIIRSMOD_201501'+str(iday)+'_'+str(iloop)+'.nc','w', format='NETCDF4') #'w' stands for write
        f = nc4.Dataset('/raid15/leipan/VIIRS/VIIRS/201501/Index/IND_CrIS_VIIRSMOD_201501'+str(iday)+'_'+str(iloop)+'.nc','w', format='NETCDF4') #'w' stands for write

        f.createDimension('m',dy_flatten.size)
        f.createDimension('x', dy.shape[0])
        f.createDimension('y', dy.shape[1])
        f.createDimension('z', dy.shape[2])

        y_flatten = f.createVariable('dy', 'i4', ('m',))
        y_size=f.createVariable('dy_size','i4',('x', 'y', 'z',))
        x_flatten = f.createVariable('dx', 'i4', ('m',))

        ### attr1 = f.createAttribute('granules', 'the info of the 3 viirs granules involved')

        y_size[:]=dy_size
        y_flatten[:]=dy_flatten
        x_flatten[:]=dx_flatten

        f.description="Demo Data for 2015 Jan"

        f.close()
        
    print("done in --- %s seconds --- for 2015 Jan %d files" % (time.time() - start_time, iday))
    start_time=time.time()
# collocation is done
"""
##############################################################################
# showing the collocated images 
#############################################################################
start_time = time.time()

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.colors as colors
import matplotlib.cm as cmx

print(cris_lon.min(),cris_lat.min(),cris_lon.max(),cris_lat.max())

m = Basemap(resolution='l', projection='cyl',  \
		llcrnrlon=cris_lon.min(), llcrnrlat=cris_lat.min(),  
        urcrnrlon=cris_lon.max(), urcrnrlat=cris_lat.max())
m.drawcoastlines()
m.drawcountries()
m.drawstates()

# meridians on bottom and left
parallels = np.arange(0.,81,10.)
m.drawparallels(parallels,labels=[False,True,True,False])
meridians = np.arange(10.,351.,20.)
m.drawmeridians(meridians,labels=[True,False,False,True])

# create color map 
jet = cm = plt.get_cmap('jet') 
#cNorm  = colors.Normalize(vmin=220, vmax=310)
cNorm  = colors.Normalize(vmin=0, vmax=1000)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

# show collocated pixels 
for k, j, i in np.ndindex(cris_lat.shape):
	
	ix=dx[k,j,i]
	iy=dy[k,j,i]
	vcolorVal = np.squeeze(scalarMap.to_rgba(viirs_height[iy, ix]))
	vx, vy = m(viirs_lon[iy, ix], viirs_lat[iy, ix])
	cs1 = m.scatter(vx, vy, s=0.5, c=vcolorVal, edgecolor='none', cmap='jet', marker='.')

plt.savefig('myfig_20150601', dpi=300)    

print("making plots is using --- %s seconds " % (time.time() - start_time))
"""


 
