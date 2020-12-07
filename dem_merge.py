# -*- coding: utf-8 -*-
"""
Created by E. Macorps

Merge two DEM datasets with overlapping surface
Requirements:
    1. pre-process reference TDX DEM 90-m with TDX_mask_bad_values.py (modified from D. Shean)
    2. align original DEMs to reference TDX DEM using dem_align_EM_edits.py (modified from D. Shean)
    3. make sure that dem_align_EM_edits.py is in the same directory as working directory
"""
#------------------------------------------------------------------------------

import csv
import os
import time
import numpy as np
import argparse
import pathlib

#--------------------------------- FUNCTIONS ----------------------------------
def SaveNewCSV(easting, northing, z_elev, output_csv):
    with open(output_csv, 'w') as csvout: # output csv file to write in
        wr = csv.writer(csvout, lineterminator='\n')
        for i in range(len(easting)):
            wr.writerow([easting[i],northing[i],z_elev[i]])
    #Exit
    csvout.close()

# Get Arguments
def getparser():
    parser = argparse.ArgumentParser(description="Merge two DEMs after performing dem_align on a processed TanDEM-X DEM", \
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('ref_dem_fn', type=str, help='reference DEM Filename')
    parser.add_argument('src_dem_fn', type=str, help='source DEM Filename')
    parser.add_argument('-xmin', type=str, help='xmin coordinate for clipping mask')
    parser.add_argument('-xmax', type=str, help='xmax coordinate for clipping mask')
    parser.add_argument('-ymin_ref', type=str, help='ymin coordinate for clipping mask of reference DEM')
    parser.add_argument('-ymin_src', type=str, help='ymin coordinate for clipping mask of source DEM')    
    parser.add_argument('-ymax_ref', type=str, help='ymax coordinate for clipping mask of reference DEM')
    parser.add_argument('-ymax_src', type=str, help='ymax coordinate for clipping mask of source DEM')
    parser.add_argument('-outres', type=str, help='output spatial resolution in meters')
    parser.add_argument('-epsg', type=str, default='32613', help='EPSG value for geographic coordinate system')
    parser.add_argument('-outname', type=str, help='output filename for merged DEM without extension')
    parser.add_argument('-outdir', type=str, default=None, help='output directory for merged DEM')
    parser.add_argument('-workdir', type=str, help='path to working directory with scripts')
    parser.add_argument('-mask_shp_fn', type=str, default=None, help='non-static surfaces shapefile for DEM align')
    return parser
    
def main(argv=None):
    start_time = time.time()

    parser = getparser()
    args = parser.parse_args()
    
    ref_dem_fn = args.ref_dem_fn
    src_dem_fn = args.src_dem_fn
    
    xmin = args.xmin
    xmax = args.xmax
    ymin_ref = args.ymin_ref
    ymin_src = args.ymin_src
    ymax_ref = args.ymax_ref
    ymax_src = args.ymax_src
    outres = args.outres
    epsg = args.epsg
    
    workdir = args.workdir
    if workdir is None:
        workdir = pathlib.Path(__file__).parent.absolute()
    mask_shp_fn = args.mask_shp_fn
    outname = args.outname
    outdir = args.outdir
    if outdir is not None:
        if not os.path.exists(outdir):
            os.makedirs(outdir)
    
    ref_dem = os.path.splitext(os.path.basename(ref_dem_fn))[0]
    src_dem = os.path.splitext(os.path.basename(src_dem_fn))[0]


    """
    Use Gdal_translate to clip reference and source DEM to respective mask extents
    """
    os.system('gdal_translate -of GTiff -a_srs EPSG:'+epsg+' '+'-tr '+outres+' '+outres+' -a_nodata 0 -r cubicspline -projwin '+xmin+' '+ymax_ref+' '+xmax+' '+ymin_ref+' -co COMPRESS=NONE '+ref_dem_fn+' '+os.path.join(outdir,ref_dem)+'_clip2ext.tif')
    os.system('gdal_translate -of GTiff -a_srs EPSG:'+epsg+' '+'-tr '+outres+' '+outres+' -a_nodata 0 -r cubicspline -projwin '+xmin+' '+ymax_src+' '+xmax+' '+ymin_src+' -co COMPRESS=NONE '+src_dem_fn+' '+os.path.join(outdir,src_dem)+'_clip2ext.tif')   
    
    # Update path and filename for clipped reference DEM
    ref_dem_clip = ref_dem+'_clip2ext'
    
    """
    Use dem_align to re-align clipped source DEM to reference DEM
    """
    os.system('python '+workdir+'/dem_align_EM_edits.py -mode nuth -mask_shp_fn '+mask_shp_fn+' -outdir '+outdir+' '+os.path.join(outdir,ref_dem_clip)+'.tif '+os.path.join(outdir,src_dem)+'_clip2ext.tif')
    
    # Get path and filename for aligned clipped source DEM
    for file in os.listdir(outdir):
        if file.endswith('_align.tif'):
            filepath = os.path.join(outdir, file)
            src_dem_clip = os.path.splitext(os.path.basename(filepath))[0]
    
    
    """
    Use Gdal_translate to convert High-res DEMs to CSV file
    Process to merge DEM products
    """
    
    os.system('gdal_translate -of AAIGrid -a_srs EPSG:32613 -a_nodata 0 -r cubicspline' + ' ' + os.path.join(outdir,ref_dem_clip) + '.tif ' + os.path.join(outdir,ref_dem_clip) + '.asc')
    os.system('gdal_translate -of xyz -a_srs EPSG:32613 -a_nodata 0 -r cubicspline'+ ' ' + os.path.join(outdir,ref_dem_clip) + '.asc ' + os.path.join(outdir,ref_dem_clip) + '.csv')
    os.system('gdal_translate -of AAIGrid -a_srs EPSG:32613 -a_nodata 0 -r cubicspline'+ ' ' + os.path.join(outdir,src_dem_clip) + '.tif ' + os.path.join(outdir,src_dem_clip) +'.asc')
    os.system('gdal_translate -of xyz -a_srs EPSG:32613 -a_nodata 0 -r cubicspline'+ ' ' + os.path.join(outdir,src_dem_clip) + '.asc ' + os.path.join(outdir,src_dem_clip) + '.csv')
    
    ref_dst = os.path.join(outdir,ref_dem_clip) + '.csv'
    print (ref_dst)
    src_dst = os.path.join(outdir,src_dem_clip) + '.csv'
    print (src_dst)
    
    with open (ref_dst) as Input_Table:
        readCSV = csv.reader (Input_Table, delimiter=' ')
        x_A, y_A, z_A = [], [], []
        for row in readCSV:
            x_A.append(float(row[0]))
            y_A.append(float(row[1]))
            z_A.append(float(row[2]))
    
    with open (src_dst) as Input_Table:
        readCSV = csv.reader (Input_Table, delimiter=' ')
        x_B, y_B, z_B = [], [], []
        for row in readCSV:
            x_B.append(float(row[0]))
            y_B.append(float(row[1]))
            z_B.append(float(row[2]))     
    
    # Need to compute displacement between two dataset:
            
    # get spatial resolution
    res = x_A[1] - x_A[0]
    
    # compute displacement in x
    x_displ = x_A[0] - x_B[0]
    # correct src_dst for displacement in x
    for i in range (len(x_B)):
        x_B[i] = x_B[i] + x_displ
    
    # Compute displacement in y
    y_A_last = y_A[len(y_A)-1]
    
    search_arr = np.where(np.round(y_B)==np.round(y_A_last))[0]
    if len(search_arr) != 0:
        index = search_arr[0]
    else:
        for n in range (1,int(res)+1):
            search_arr = np.where(np.round(y_B) == np.round(y_A_last) + n)[0]
            if len(search_arr) == 0:
                n += 1
            else:
                index = search_arr[0]
    
    y_displ = y_A_last - y_B[index] # re-compute displacement with correct index to know if positive or negative value
    
    # Correct src_dst for displacements in y
    for i in range (len(y_B)):
        y_B[i] = y_B[i] + y_displ
    
    x_A_arr = np.array(x_A)
    y_A_arr = np.array(y_A)
    x_B_arr = np.array(x_B)
    y_B_arr = np.array(y_B)
    
    # overlap area extent
    xmin = x_A[0]
    xmax = x_A[len(x_A)-1]
    ymin = y_A[len(y_A)-1]
    ymax = y_B[0]
    
    # find indexes for each point file for overlap area
    list_index_overlap4arr_A = np.where((ymin <= y_A_arr) & (y_A_arr <= ymax) & (xmin <= x_A_arr) & (x_A_arr <= xmax))[0]
    list_index_overlap4arr_B = np.where((ymin <= y_B_arr) & (y_B_arr <= ymax) & (xmin <= x_B_arr) & (x_B_arr <= xmax))[0]
    # find indexes for DEM A only
    list_index_arr_A_only = []
    for i in range(list_index_overlap4arr_A[0]):
        list_index_arr_A_only.append(i)
    
    # find indexes for DEM B only    
    list_index_arr_B_only = []
    for i in range(list_index_overlap4arr_B[len(list_index_overlap4arr_B)-1]+1,len(x_B)):
        list_index_arr_B_only.append(i)
    
    # Process to merge
    
    easting, northing, z_elev = [], [], []
    
    for i in list_index_arr_A_only: # Use only Pleiades
        easting.append(x_A[i])
        northing.append(y_A[i])
        z_elev.append(z_A[i])
     
    for j in range(len(list_index_overlap4arr_A)): # Use weighted average of SPOT and Pleiades elevations
        z_average = (0.2 * z_B[list_index_overlap4arr_B[j]])+(0.8 * z_A[list_index_overlap4arr_A[j]])
        z_elev.append(z_average)
        easting.append(x_A[list_index_overlap4arr_A[j]])
        northing.append(y_A[list_index_overlap4arr_A[j]])
    
    for k in list_index_arr_B_only: # Use SPOT only
        easting.append(x_B[k])
        northing.append(y_B[k])
        z_elev.append(z_B[k])
               
    # Save output into CSV file
    output_csv = os.path.join(outdir,outname) + '.csv'
    SaveNewCSV(easting,northing,z_elev,output_csv)
    print ('merged point file created in '+ output_csv)
    print ('now converting into GTiff')
    
    
#------------------------------------------------------------------------------
# CONVERT CSV FILE TO RASTER using gdal_translate
    
    os.system('gdal_translate -of xyz -a_srs EPSG:32613 -a_nodata 0 -r cubicspline' + ' ' + os.path.join(outdir,outname)+'.csv' + ' ' + os.path.join(outdir,outname)+'.xyz')
    os.system('gdal_translate -of AAIGrid -a_srs EPSG:32613 -a_nodata 0 -r cubicspline' + ' ' + os.path.join(outdir,outname)+'.xyz'+ ' ' + os.path.join(outdir,outname) + '.asc')
    os.system('gdal_translate -of AAIGrid -a_srs EPSG:32613 -a_nodata 0 -r cubicspline' + ' ' + os.path.join(outdir,outname) + '.asc'+ ' ' + os.path.join(outdir,outname) + '.tif')
    
    
    elapsed_time = time.time() - start_time
    print ('Elapsed time after processing is ' + str(elapsed_time) + ' seconds')
    
if __name__ == "__main__":
    main()


#------------------------------------------------------------------------------
# xmin = 641597.00  
# xmax = 647712.00
# ymin_ref = 2151930.00
# ymax_ref = 2161128.00
# ymin_src = 2147270.00
# ymax_src = 2153900.00
# epsg = 32613
# outname = merge_dem_outfile_newtest
# outdir = /Users/emacorps/workspace/SPOT6_DEM/merge_dems/merge_dem_newtest
# workdir = /Users/emacorps/workspace/SPOT6_DEM/merge_dems
# mask_shp_fn = /Users/emacorps/workspace/SPOT6_DEM/merge_dems/non_static_surface.shp
#


    