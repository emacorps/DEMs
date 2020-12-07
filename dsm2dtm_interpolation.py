# -*- coding: utf-8 -*-
"""
Created on Tues, July 23 2019

@author: E. Macorps
Interpolation of elevation values using forest mask created from supervised classification - To be used in QGIS-2

"""
import sys
import time
import datetime
import os, os.path
from getpass import getuser
import csv
import math
import numpy as np
import ogr,osr
from osgeo import gdal
from osgeo import gdal_array
from scipy import interpolate
import matplotlib.pyplot as plt
# QGIS modules
from qgis.core import *
# Warning: If modules not detected, run: export PYTHONPATH=/Applications/QGIS.app/Contents/Resources/python
import qgis.utils
from PyQt4.QtCore import *
from qgis.PyQt.QtCore import QVariant
from PyQt4.QtGui import QApplication
import processing

# Filenames
proc_path = 'Users/emacorps/workspace/dsm2dtm_interp/'
dem = proc_path + 'pleiades_30m_resamp.tif'
mask = proc_path + 'reclass_landsat_ground_trees.tif'
#coordinates for clip
xmin = 641605.000
xmax = 647899.000
ymin = 2148680.000
ymax = 2161036.000
#Other variables
interp_method = 'linear' #choose from 'linear' 'cubic' 'nearest'

# Open with Qgis to read resolution
mask_qgis = QgsRasterLayer(mask,'forest_mask_from_supervised_classification')
mask_res = mask_qgis.rasterUnitsPerPixelX()

# Resample DEM to MASK resolution
processing.runalg("grass7:r.resample",dem,"%f,%f,%f,%f"%(xmin,xmax,ymin,ymax),mask_res,proc_path+'dem_30m_resamp.tif')
dem_resamp = proc_path+'dem_30m_resamp.tif'

# Clip DEM and MASK to same extent
processing.runalg("gdalogr:cliprasterbyextent",mask,-9999,"%f,%f,%f,%f"%(xmin,xmax,ymin,ymax),5,0,75,6,1,False,0,True,None,proc_path+'mask_clip.tif')
processing.runalg("gdalogr:cliprasterbyextent",dem_resamp,-9999,"%f,%f,%f,%f"%(xmin,xmax,ymin,ymax),5,0,75,6,1,False,0,True,None,proc_path+'dem_resamp_clip.tif')

DEM = proc_path+'dem_resamp_clip.tif'
MASK = proc_path+'mask_resamp_clip.tif'

def dsm2dtm_interpolation (DEM, MASK, interp_method):
    # open dem and mask with GDAL
    src_dem = gdal.Open(DEM)
    src_mask = gdal.Open(MASK)
    # get dimensions of dem
    pixels = src_dem.RasterXSize
    lines = src_dem.RasterYSize
    # get band gdal
    band_dem = src_dem.GetRasterBand(1)
    band_mask = src_mask.GetRasterBand(1)
    # get band datatype and numpy type
    gdal_band_type = band_dem.DataType
    np_band_type = gdal_array.GDALTypeCodeToNumericTypeCode(gdal_band_type)
    # read gdal raster into array
    array_dem = band_dem.ReadAsArray(0,0,pixels,lines).astype(np_band_type)
    array_mask = band_mask.ReadAsArray(0,0,pixels,lines).astype(np_band_type)
    # get geotransform and projection information
    geotrans = src_dem.GetGeoTransform()
    projInfo = src_dem.GetProjection()
    # close band dataset
    src_dem = None
    src_mask = None

    #Create empty arrays
    temp_dem = np.zeros ((lines,pixels), dtype = np_band_type)
    orig_dem = np.zeros ((lines,pixels), dtype = np_band_type)
    interp_dtm = np.zeros ((lines,pixels), dtype = np_band_type)
    forest_mask = np.zeros ((lines,pixels), dtype = np_band_type)

    # get data to fill empty arrays
    temp_dem = array_dem
    orig_dem = array_dem
    forest_mask = array_mask

    #Use mask to remove values where mask == 0 (vegetation). When mask == 1: bare ground, keep elevation values
    trees_coord = np.where(forest_mask==0) #matrix with row and col coordinates where forest_mask = 0
    row_trees = trees_coord[0] #row coordinates
    col_trees = trees_coord[1] #column coordinates
    
    for i in trees_coord: #Masking DEM for forested area, assigning -1 values where "trees" are; leaving elevation value from DEM in others
        temp_dem[row_trees[i],col_trees[i]] = -1
   
    # Create empty meshgrid with dimensions lines x pixels:
    x = range (0,pixels,1)
    y = range (0,lines,1)

    X, Y = np.meshgrid(x, y)

    # Interpolation
    interp_dtm = interpolate.griddata((X[temp_dem!=-1],Y[temp_dem!=-1]),orig_dem[temp_dem!=-1], (X, Y), interp_method)
    
    # Create output and apply projection and transform info:
    driver = gdal.GetDriverByName('GTiff')
    dst_ds = driver.Create(proc_path+'/interp_dtm_'+interp_method+'.tif',pixels,lines,1,gdal.GDT_Float64)
    dst_ds.SetGeoTransform(geotrans)
    sp_ref=osr.SpatialReference()
    sp_ref.ImportFromWkt(projInfo)
    dst_ds.SetProjection(sp_ref.ExportToWkt())

    #Write interpolated array to band:
    dst_ds.GetRasterBand(1).WriteArray(temp_dem)    
    
    # Display rasters
    plt.subplot (2,2,1)
    plt.imshow(forest_mask)
    plt.title('forest mask')
    plt.subplot (2,2,2)
    plt.imshow(array_dem)
    plt.title('DSM')
    plt.subplot (2,2,3)
    plt.imshow(interp_dtm)
    plt.title('interpolated DTM')
    plt.show
    
    #Flush to disk and close file:
    dst_ds.FlushCache()
    dst_ds=None

#Call dsm2dtm_interpolation function to run interpolation on DEM with Forest mask created using supervised classification in QGIS:    
dsm2dtm_interpolation (DEM, MASK, interp_method)


    
    
