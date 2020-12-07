#! /usr/bin/env python
# -*- coding: utf-8 -*-
# ----------------------------------------------------------------------------------------------------------
#   Created on Mar 19 12:54:18 2019

#   Author: Elodie Macorps

#   This script computes the elevation difference between 2 DEMs by:
#   (1) resampling to the same resolution
#   (2) clipping both rasters to the same extent
#   (3) converting to array and substracting cell by cell


# --------------------------------------------- Module Import ---------------------------------------------

import os, os.path
import numpy as np
import gdal
import osr
from osgeo import gdal_array
import argparse
import matplotlib.pyplot as plt
from imview.lib import pltlib
from pygeotools.lib import malib
from mpl_toolkits.axes_grid1 import make_axes_locatable
import rasterio
from rasterio.plot import show

# ------------------------------------------------ Function ------------------------------------------------

def getparser():
    parser = argparse.ArgumentParser(description="Compute Raster Difference for a given extent", \
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('DEM_A', type=str, help='DEM A Filename')
    parser.add_argument('DEM_B', type=str, help='DEM B Filename to be substracted to DEM A')
    parser.add_argument('-res', type=float, default=5.0, help='output resolution for raster difference')
    parser.add_argument('-xmin', type=str, help='Xmin coordinate for clipping extent')
    parser.add_argument('-ymin', type=str, help='Ymin coordinate for clipping extent')
    parser.add_argument('-xmax', type=str, help='Xmax coordinate for clipping extent')
    parser.add_argument('-ymax', type=str, help='Ymax coordinate for clipping extent')
    parser.add_argument('-outname', type=str, help='Output filename')
    parser.add_argument('-outdir', type=str, default=None, help='Output directory')
    return parser

def main(argv=None):
    parser = getparser()
    args = parser.parse_args()
    #Should check that files exist
    DEM_A = args.DEM_A
    DEM_B = args.DEM_B

    dem_res = args.res
    xmin = args.xmin
    ymin = args.ymin
    xmax = args.xmax
    ymax = args.ymax
    outname = args.outname
    outdir = args.outdir
    if outdir is None:
        outdir = os.path.splitext(DEM_A)[0] + '_raster_diff/'
    
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    print("\nRASTER A: %s" % DEM_A)
    print("RASTER B: %s" % DEM_B)
    print("Output: %s\n" % os.path.join(outdir,outname))

    # Rasters MUST HAVE the same extent AND the same resolution
    # 1. Warp and clip to same resolution
    os.system('gdal_translate -of GTiff -r cubicspline -a_srs EPSG:32613 -a_nodata 0 -ot Float32 -tr '+str(dem_res)+' '+str(dem_res)+' -projwin ' + xmin +' '+ ymax +' '+ xmax +' '+ ymin +' '+'-co COMPRESS=NONE'+' '+ DEM_A +' '+ outdir+'/dem_A_clip.tif')
    os.system('gdal_translate -of GTiff -r cubicspline -a_srs EPSG:32613 -a_nodata 0 -ot Float32 -tr '+str(dem_res)+' '+str(dem_res)+' -projwin ' + xmin +' '+ ymax +' '+ xmax +' '+ ymin +' '+'-co COMPRESS=NONE'+' '+ DEM_B +' '+ outdir+'/dem_B_clip.tif')
    DEM_A_clip = outdir+'/dem_A_clip.tif'
    DEM_B_clip = outdir+'/dem_B_clip.tif'

    # Compute DEM difference between clipped and warpped rasters
    # Open DEM with Gdal
    gdal_DEM_A = gdal.Open(DEM_A_clip)
    gdal_DEM_B = gdal.Open(DEM_B_clip)
    # Get array dimensions
    pixels = gdal_DEM_A.RasterXSize
    lines = gdal_DEM_A.RasterYSize
    # Get band GDAL
    gdal_DEM_A_band = gdal_DEM_A.GetRasterBand(1)
    gdal_DEM_B_band = gdal_DEM_B.GetRasterBand(1)
    # Get band GDAL datatype and map to Numpy Array
    gdal_band_type = gdal_DEM_A_band.DataType
    np_band_type = gdal_array.GDALTypeCodeToNumericTypeCode(gdal_band_type)
    # Read the DEM rasters into arrays
    DEM_A_array = gdal_DEM_A_band.ReadAsArray(0,0,pixels,lines).astype(np_band_type)
    DEM_B_array = gdal_DEM_B_band.ReadAsArray(0,0,pixels,lines).astype(np_band_type)
    # Get Geotransform infos
    geotrans = gdal_DEM_A.GetGeoTransform()
    projInfo = gdal_DEM_A.GetProjection()
        
    # Create output and apply projection and transform info
    driver = gdal.GetDriverByName("GTiff")
    dst_ds = driver.Create(outdir+outname,pixels,lines,1,gdal.GDT_Float32)
    dst_ds.SetGeoTransform(geotrans)
    sp_ref = osr.SpatialReference()
    sp_ref.ImportFromWkt(projInfo)
    dst_ds.SetProjection(sp_ref.ExportToWkt())
    dst_ds.GetRasterBand(1).SetNoDataValue(-9999)

    # Create empty array
    diff_dem = np.zeros((lines,pixels), dtype = np_band_type)
        
    # Set i and j as Col and Row coordinates
    i = range(0,lines,1) #Row coord
    j = range(0,pixels,1) #Col coord
    # Compute DEM difference between top surface and channel
    for i in range(0,lines,1):
        for j in range(0,pixels,1):
            diff_dem[i,j] = DEM_A_array[i,j] - DEM_B_array[i,j]

    # Write diff_dem array to band
    dst_ds.GetRasterBand(1).WriteArray(diff_dem)
    #dst_ds = driver.Create(proc_path+output_name,pixels,lines,1,gdal.GDT_Float32)
    #dst_ds.SetGeoTransform(geotrans)
    #sp_ref = osr.SpatialReference()
    #sp_ref.ImportFromWkt(projInfo)
    #dst_ds.SetProjection(sp_ref.ExportToWkt())
    #dst_ds.GetRasterBand(1).SetNoDataValue(-9999)
    
    #Close files
    gdal_DEM_A = None
    gdal_DEM_B = None
    gdal_DEM_A_band = None
    gdal_DEM_B_band = None
    dst_ds.FlushCache()
    dst_ds = None

    new_tiff_file = outdir+outname
    with rasterio.open(new_tiff_file) as dem_diff_tif:
        dem_diff = dem_diff_tif.read(1)
        
    #Create Figure
    print ("Creating final plot")
    dz_clim = malib.calcperc_sym(dem_diff, (5, 95))
    fig, ax = plt.subplots(figsize=(10,5))
    im = ax.imshow(dem_diff, cmap='coolwarm', clim=dz_clim)
    show(dem_diff, ax=ax, cmap='coolwarm', clim=dz_clim, transform=dem_diff_tif.transform)
    #plt.axis('off')
    divider = make_axes_locatable(ax)
    cax1 = divider.append_axes("right", size="2%", pad=0.05)
    fig.colorbar(im, cax=cax1, extend='both',ticks=range(-10,15,5),label='Elevation difference (m)',aspect=5).ax.set_yticklabels(['< -10','-5','0','5','> 10'])
    ax.set_title('Raster Difference') 
    xstart, xend = ax.get_xlim()
    ax.xaxis.set_ticks(np.arange(xstart,xend,1000))
    ystart, yend = ax.get_ylim()
    ax.yaxis.set_ticks(np.arange(ystart,yend,1000))

    #Save Figure to .png file
    fig_fn = outdir + outname + '.png'
    print("Writing out figure")
    fig.savefig(fig_fn, dpi=300)

if __name__ == "__main__":
    main()

# ------------------------------------------ Set up variables ------------------------------------------
"""
proc_path = 'Users/emacorps/workspace/SPOT6_DEM/merge_dems/spot_5m_clip_ext_dem_align/raster_diff/'
DEM_A = proc_path+'spot_5m_align_tdx_full_ext.tif'
DEM_B = proc_path+'spot_5m_align_tdx_clip_ext.tif'
dem_res = 5.0000
output_name = 'spot_full_clip_diff.tif'
xmin = 641592.000
xmax = 647577.000
ymin = 2151931.000
ymax = 2153902.000
"""