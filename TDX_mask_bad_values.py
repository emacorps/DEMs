#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 15:14:56 2020

@author: emacorps
"""

import os
import numpy as np
import glob
from pygeotools.lib import iolib
import scipy.ndimage as ndimage
import argparse

"""
Script to mask TanDEM-X DEMs to isolate good elevation values
Initially written by D. Shean DOI:10.5281/zenodo.3451468
Documentation of AUX products is here: https://tandemx-science.dlr.de/pdfs/TD-GS-PS-0021_DEM-Product-Specification_v3.1.pdf
"""
def getparser():
    parser = argparse.ArgumentParser(description="Script to mask TanDEM-X DEMs to isolate good elevation values. Initially written by D. Shean DOI:10.5281/zenodo.3451468 Outputs in initial coordinate system and UTM", \
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('indir', type=str, help='TanDEM-X tile directory with DEM and AUX products. Directory must end with "_" and all files must have a common name. Example is "/TDM1_DEM__30_N19W104_"')
    parser.add_argument('-epsg', type=int, default=None ,help='EPSG code for UTM coordinate system transformation')
    return parser
    
def main(argv=None):
    parser = getparser()
    args = parser.parse_args()

    tiledir = args.indir
    epsg = args.epsg
    
    print(tiledir)
    dem_fn = glob.glob(os.path.join(tiledir, 'DEM/*DEM.tif'))[0]
    print(dem_fn)
    dem_ds = iolib.fn_getds(dem_fn)
    dem = iolib.ds_getma(dem_ds)
    
    #Get original mask, True where masked
    mask = np.ma.getmaskarray(dem)
    
    #Theoretical height error
    err_fn = glob.glob(os.path.join(tiledir, 'AUXFILES/*HEM.tif'))[0]
    err = iolib.fn_getma(err_fn)
    max_err_multi = 1.5
    mask = np.logical_or(mask, (err.data > max_err_multi))
    
    #Water mask
    wam_fn = glob.glob(os.path.join(tiledir, 'AUXFILES/*WAM.tif'))[0]
    wam = iolib.fn_getma(wam_fn)
    wam_clim = (33,127)
    #wam_clim = (3,127)
    mask = np.logical_or(mask, (wam >= wam_clim[0]) & (wam <= wam_clim[1]))
    
    #Consistency mask
    com_fn = glob.glob(os.path.join(tiledir, 'AUXFILES/*COM.tif'))[0]
    com = iolib.fn_getma(com_fn)
    com_valid = (8,9,10)
    #4 is only one obs
    #com_invalid = (0,1,2,4)
    com_invalid = (0,1,2)
    mask = np.logical_or(mask, np.isin(com.data, com_invalid))
    
    #More stringent error threshold for single strip pixels
    #max_err_single = 1.0
    #mask = np.logical_or(mask, (com.data == 4) & (err.data > max_err_single))
    
    #Apply
    dem_masked = np.ma.array(dem, mask=mask)
    out_fn = os.path.splitext(dem_fn)[0]+'_masked.tif'
    iolib.writeGTiff(dem_masked, out_fn, dem_ds)
    
    #Dilate mask by n_iter px to remove isolated pixels and values around nodata 
    n_iter = 1
    mask = ndimage.morphology.binary_dilation(mask, iterations=n_iter)
    #To keep valid edges, do subsequent erosion 
    mask = ndimage.morphology.binary_erosion(mask, iterations=n_iter)
    #(dilation of inverted mask, to avoid maasking outer edge)
    #mask = ~(ndimage.morphology.binary_dilation(~mask, iterations=n_iter))
    
    #Apply
    dem_masked = np.ma.array(dem, mask=mask)
    out_fn = os.path.splitext(dem_fn)[0]+'_masked_erode.tif'
    iolib.writeGTiff(dem_masked, out_fn, dem_ds)
    
    ref_dem = out_fn
    print (ref_dem)
    
    # Warp to UTM coord system to match High-res DEMs
    if epsg is not None:
        ref_dem_utm = os.path.splitext(ref_dem)[0]+'_utm.tif'
        os.system('gdalwarp -ot Float32 -t_srs EPSG:'+str(epsg)+' '+'-r cubicspline -of GTiff -co COMPRESS=NONE -wo OPTIMIZE_SIZE=TRUE' + ' ' + ref_dem + ' ' + ref_dem_utm)

if __name__ == "__main__":
    main()
