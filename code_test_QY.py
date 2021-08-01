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
"""
dataDir1='/peate_archive/.data5/Ops/npp/noaa/op/2012/05/15/scris/'
dataDir2='/peate_archive/.data5/Ops/npp/noaa/op/2012/05/15/gcrso/'
dataDir3='/peate_archive/.data5/Ops/npp/noaa/op/2012/05/15/svm15/'
dataDir4='/peate_archive/.data5/Ops/npp/noaa/op/2012/05/15/gmodo/'
"""
### for iday in range(15,23,1):
if True:
#dataDir2='/peate_archive/.data6/Ops/snpp/gdisc/2/2015/06/01/crisl1b/'
    ### dataDir2='/peate_archive/.data1/Ops/snpp/gdisc/2/2015/01/'+str(iday).zfill(2)+'/crisl1b/'
    dataDir2='./'
    ### dataDir4='/raid15/qyue/VIIRS/VIIRS/201501/'
    ### dataDir4='/raid15/qyue/VIIRS/VIIRS/201501/VNP03MOD/'
    dataDir4='./'
    
    ### for iloop in range(0,239,10):
    ### for iloop in range(0,9,10):
    if True:
        ### print(iloop)   
# get CrIS files 
#cris_sdr_files = sorted(glob.glob(dataDir1+'SCRIS*d2012*'))[21:40]
        ### cris_geo_files = sorted(glob.glob(dataDir2+'SNDR.SNPP.CRIS*'))[iloop:iloop+10]
        cris_geo_files = sorted(glob.glob(dataDir2+'SNDR.SNPP.CRIS*'))
        print ('cris_geo_files: ', cris_geo_files)
# get VIIRS files 
#viirs_sdr_files = sorted(glob.glob(dataDir3+'SVM15*d2012*'))[31:59]
        ### viirs_geo_files = sorted(glob.glob(dataDir4+'VNP03MOD*A2015'+str(iday).zfill(3)+'*'))[iloop:iloop+10]
        ### viirs_geo_files = sorted(glob.glob(dataDir4+'VNP03MOD*A2015'+'*'))
        viirs_geo_files = sorted(glob.glob(dataDir4+'VNP03MOD*A*'+'*'))
        print ('viirs_geo_files: ', viirs_geo_files)


# read VIIRS data 
        viirs_lon, viirs_lat, viirs_satAzimuth, viirs_satRange, viirs_satZenith, viirs_height, viirs_time = geo_QY.read_nasa_viirs_geo(viirs_geo_files)
        ### print ('viirs_time: ', viirs_time)
        ### print ('type(viirs_time): ', type(viirs_time))
        ### print ('viirs_time: ', viirs_time)
        print ('viirs_time.shape: ', viirs_time.shape)
        ### print ('viirs_time.min(): ', viirs_time.min())
        ### print ('viirs_time.max(): ', viirs_time.max())

        ### print ('viirs_lon: ', viirs_lon)
        ### print ('type(viirs_lon): ', type(viirs_lon))
        print ('viirs_lon.shape: ', viirs_lon.shape)

        start_time = viirs_time.min()
        end_time = viirs_time.max()


#viirs_bt, viirs_rad, viirs_sdrQa = geo.read_viirs_sdr(viirs_sdr_files)i


# read CrIS data 
        cris_lon, cris_lat, cris_satAzimuth, cris_satRange, cris_satZenith, cris_time, cris_realLW = geo_QY.read_nasa_cris_geo(cris_geo_files)
        ### print ('cris_time: ', cris_time)
        print ('cris_time.min(): ', cris_time.min())
        print ('cris_time.max(): ', cris_time.max())

        if start_time < cris_time.min():
          start_time = cris_time.min()

        if end_time > cris_time.max():
          end_time = cris_time.max()

        print ('start_time: ', start_time)
        print ('end_time: ', end_time)

        # CrIS and VIIRS use epoch time since 1/1/1993 (1993TAI),
        # and unix epoch time is since 1/1/1970
        # there is a 23 year difference
        diff = (datetime(1993,1,1,0,0) - datetime(1970,1,1)).total_seconds()
        print('diff: ', diff)
        start_time += diff
        end_time += diff

        os.environ['TZ'] = 'GMT'
        time.tzset()
        start_date = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.localtime(start_time))
        end_date = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.localtime(end_time))

        print ('start_date: ', start_date)
        print ('end_date: ', end_date)

        start_date2 = time.strftime('%Y%m%dT%H%M%S', time.localtime(start_time))
        end_date2 = time.strftime('%Y%m%dT%H%M%S', time.localtime(end_time))

        print ('start_date2: ', start_date2)
        print ('end_date2: ', end_date2)

        # changed per Qing's suggestion
        ### output_filename = 'IND_CrIS_VIIRSMOD_' + start_date2 + '_' + end_date2

        # Qing's suggestion:
        # use this name: IND_CrIS_VIIRSMOD_SNDR.SNPP.20150601T1548.g159.nc
        # for this CRIS granule: SNDR.SNPP.CRIS.20150601T1548.m06.g159.L1B_NSR.std.v02_05.G.180904193415.nc
        cris_geo_file = os.path.basename(cris_geo_files[0])
        print ('cris_geo_file: ', cris_geo_file)
        split1 = cris_geo_file.split('.')
        print ('split1: ', split1)
        # split1:  ['SNDR', 'SNPP', 'CRIS', '20150603T1836', 'm06', 'g187', 'L1B_NSR', 'std', 'v02_05', 'G', '180905023033', 'nc']
        output_filename = 'IND_CrIS_VIIRSMOD_' + split1[0] + '.' + split1[1] + '.' + split1[3] + '.' + split1[5]
        print ('output_filename: ', output_filename)

        if os.path.exists(output_filename):
          shutil.rmtree(output_filename)

        ### sys.exit(0)
        os.mkdir(output_filename)

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

        print ('cris_los.shape: ', cris_los.shape)

# compute viirs POS vector in ECEF
        viirs_pos= np.zeros(np.append(viirs_lat.shape, 3))
        viirs_pos[:, :, 0], viirs_pos[:, :, 1], viirs_pos[:, :, 2] = \
	    geo_QY.LLA2ECEF(viirs_lon, viirs_lat, np.zeros_like(viirs_lat))

        print ('viirs_pos.shape: ', viirs_pos.shape)

# cris_los is pointing from pixel to satellite, we need to
#   change from satellite to pixel
        cris_los = -1.0*cris_los

# using Kd-tree to find the closted pixel of VIIRS for each CrIS FOV
# Set fake viirs_sdrQa to be zero: good quality everywhere since not for calibration
#viirs_sdrQa=np.zeros(viirs_lon.shape)

#remove the sdrqa, but adding time requirement (less than 600S difference)
        dy, dx = geo_QY.match_cris_viirs_QY(cris_los, cris_pos, viirs_pos, cris_time, viirs_time)

        ### print ('dy: ', dy)
        print ('dy.shape: ', dy.shape)
        ### print ('dx: ', dx)
        print ('dx.shape: ', dx.shape)

#	print("collocation are done in --- %s seconds --- for %d files " % (time.time() - start_time, len(cris_geo_files)))

        dy_flatten = np.array([item for lst in dy.reshape(-1) for item in lst])
        dy_size = np.array([len(lst) for lst in dy.reshape(-1)]).reshape(dy.shape)
        dx_flatten = np.array([item for lst in dx.reshape(-1) for item in lst])

        ### f = nc4.Dataset('/raid15/qyue/VIIRS/VIIRS/201501/Index/IND_CrIS_VIIRSMOD_201501'+str(iday)+'_'+str(iloop)+'.nc','w', format='NETCDF4') #'w' stands for write
        ### f = nc4.Dataset('/raid15/leipan/VIIRS/VIIRS/201501/Index/IND_CrIS_VIIRSMOD_201501'+str(iday)+'_'+str(iloop)+'.nc','w', format='NETCDF4') #'w' stands for write
        ### f = nc4.Dataset('./IND_CrIS_VIIRSMOD_201501'+'.nc','w', format='NETCDF4') #'w' stands for write

        # make it a real standard netcdf product
        # e.g., include long name of all the variables

        f = nc4.Dataset(output_filename+'/'+output_filename+'.nc','w', format='NETCDF4') #'w' stands for write

        f.createDimension('m',dy_flatten.size)
        f.createDimension('x', dy.shape[0])
        f.createDimension('y', dy.shape[1])
        f.createDimension('z', dy.shape[2])

        y_flatten = f.createVariable('dy', 'i4', ('m',))
        y_flatten.setncatts({'long_name':u'y coords', 'units':u'meter', 'var_desc':u'Y coordnates'})
        y_size=f.createVariable('dy_size','i4',('x', 'y', 'z',))
        y_size.setncatts({'long_name':u'y size', 'units':u'none', 'var_desc':u'Y dimension size'})
        x_flatten = f.createVariable('dx', 'i4', ('m',))
        x_flatten.setncatts({'long_name':u'x coords', 'units':u'meter', 'var_desc':u'X coordnates'})

        print ('dx_flatten.shape: ', dx_flatten.shape)

        y_size[:]=dy_size
        y_flatten[:]=dy_flatten
        x_flatten[:]=dx_flatten

        # add global attributes

        viirs_str = ''
        for item1 in viirs_geo_files:
          viirs_str += ',' + os.path.basename(item1)

        viirs_str = viirs_str[1:]
        f.viirs_file_names = viirs_str

        ### f.cris_file_name = os.path.basename(cris_geo_file)
        cris_str = ''
        for item1 in cris_geo_files:
          cris_str += ',' + os.path.basename(item1)

        cris_str = cris_str[1:]
        f.cris_file_names = cris_str

        f.cris_start_time = start_date
        f.cris_end_time = end_date

        f.cris_min_lat = cris_lat.min()
        f.cris_min_lon = cris_lon.min()
        f.cris_max_lat = cris_lat.max()
        f.cris_max_lon = cris_lon.max()


        f.description="Co-location Data for 2015 Jan"

        f.close()


# datetime object containing current date and time
now = datetime.now()
print("now: ", now)

# dd/mm/YY H:M:S
### dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
dt_string = now.strftime("%Y-%m-%dT%H:%M:%SZ")
print("date and time =", dt_string)	

d1 = \
{
    "creation_timestamp": dt_string,
    "version": "v1.0",
    "starttime": start_date,
    "endtime": end_date,
    "label": "matchup_cris_viirs_"+ start_date2 + '_' + end_date2
}
with open(output_filename+'/'+output_filename+'.dataset.json', 'w') as datasetf:
    json.dump(d1, datasetf, indent=2)

d2 = {}
with open(output_filename+'/'+output_filename+'.met.json', 'w') as metf:
    json.dump(d2, metf, indent=2)


print("done in --- %.2f seconds --- " % (float(time.time() - start_t)))

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


 
