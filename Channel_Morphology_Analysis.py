# -*- coding: utf-8 -*-
# ----------------------------------------------------------------------------------------------------------
#   Created on Mar 19 12:54:18 2019
    
#   @author: Elodie Macorps
    
#   This script creates:
#   1- transects across the valley to calculate cross-sections
#   2- segments of centerline of the valleys to calculate slope
#   3- polygons within each valleys to calculate channel capacity
#   4- changes polygon color depending on % change in channel capacity

    
#   Only inputs required are the polylines of the valley edges


# --------------------------------------------- Module Import ---------------------------------------------
# Python modules
import sys
import time
import datetime
import os, os.path
import glob
from getpass import getuser
import csv
import math
import numpy as np
import gdal
import ogr, osr
from osgeo import gdal_array
# QGIS modules
from qgis.core import *
# Warning: If modules not detected, run: export PYTHONPATH=/Applications/QGIS.app/Contents/Resources/python
import qgis.utils
from PyQt4.QtCore import *
from qgis.PyQt.QtCore import QVariant
from PyQt4.QtGui import QApplication, QColor
import processing

# --------------------- [parameter_name] = [parameter_type] [default_value_optional] ---------------------

##proc_path = string C:/workspace/Channel_Morpho_tests
##dem_name = string /PCI_Pleiades1m_2013.tif
##dem_source = string Pleiades
##year_dem = integer 2013
##equidistance = integer 50
##valley_edge_left = vector
##valley_edge_right = vector
##display = boolean True

# ---------------------------------- Set your workspace processing path ----------------------------------
proc_path = 'C:/workspace/Channel_Morpho_tests'
if not os.path.exists(proc_path+'/results'):
    os.mkdir(proc_path+'/results')
# -------------------------------------------- Set variables ---------------------------------------------
dem_name = '/PCI_Pleiades1m_2013.tif'
year_dem = '2014'
display = 'True'
equidistance = 50
dem_source = 'Pleiades'
dem_pre_erupt = QgsRasterLayer(proc_path+dem_name,dem_name)
dem_res = dem_pre_erupt.rasterUnitsPerPixelX()

if display == 'True':
    # Load layers to map
    iface.addRasterLayer(proc_path+dem_name,dem_source+'_dem_'+year_dem)
# ---------------------------------------- Cross-sectional areas -----------------------------------------
def ComputeCrossSectionalAreas(csvfile,profiles):
    with open(csvfile,"a") as output:
        writer = csv.writer(output,lineterminator='\n')
        writer.writerow(['Transect_ID','Distance along channel','cross-sectional area(m2)'])
    # Calculate cross-sectional areas for each profile and save to excel file
    dist = 0
    features = [feat for feat in profiles.getFeatures()]
    list = []
    for i in range (1,len(features)):
        line_id_i = features[i].attributes()[0]
        line_id_i_1 = features[i-1].attributes()[0]
        z_values = features[i-1].attributes()[6]
        if line_id_i == line_id_i_1: #save z_values to list
            list.append(z_values)
        else: #save to file and reset list
            print i
            list.append(z_values)
            print z_values
            max_z = max(list)
            print max_z
            sum = 0
            for value in list:
                cross_sec = (max_z-value)*dem_res
                sum = sum + cross_sec
            print sum
            with open(csvfile,"a") as file:
                writer=csv.writer(file,lineterminator='\n')
                dist += equidistance
                writer.writerow([line_id_i_1,dist,sum])
            #start new list
            list=[]
        if i == len(features)-1:
            list.append(z_values)
            max_z = max(list)
            sum = 0
            for value in list:
                cross_sec = (max_z-value)*dem_res
                sum = sum + cross_sec
            with open(csvfile,"a") as file:
                writer=csv.writer(file,lineterminator='\n')
                dist += equidistance
                writer.writerow([line_id_i,dist,sum])
        else: 
            continue

# -------------------------------------- Channel slope and bearing ---------------------------------------
def ComputeChannelSlopeAndBearing(csvfile,centerline):
    with open(csvfile,"a") as output:
        writer = csv.writer(output,lineterminator='\n')
        writer.writerow(['Segment_ID','Distance along channel','Slope','Azimuth'])
    dist = 0
    for feat in centerline.getFeatures():
        x_start = feat.geometry().asPolyline()[0][0]
        x_end = feat.geometry().asPolyline()[1][0]
        y_start = feat.geometry().asPolyline()[0][1]
        y_end = feat.geometry().asPolyline()[1][1]
        z_start = feat.attributes()[2]
        z_end = feat.attributes()[3]
        length = feat.attributes()[1]
        slope = math.degrees (math.atan ((math.fabs(z_start-z_end))/length))
        bearing = math.degrees (math.atan ((y_end-y_start)/(x_end-x_start)))
        if bearing > 0 :
            degree_change = (90 - bearing)
        elif bearing < 0 :
            degree_change = -(90 + bearing)
        with open(csvfile,"a") as file:
            writer = csv.writer(file,lineterminator='\n')
            dist += equidistance
            writer.writerow([feat.attributes()[0],dist,slope,degree_change])

# ------------------------------------------- Channel Capacity -------------------------------------------
def ComputeChannelCapacity(csvfile, zonal_stats):
    list = []
    capacity = 0.0
    for zones in zonal_stats.getFeatures():
        capacity = zones.attributes()[4]
        list.append (capacity)
    # Write all output diff_dem data
    with open(csvfile,"a") as output:
        writer = csv.writer(output,lineterminator='\n')
        writer.writerow(['Polygon_ID','Distance along channel','channel capacity(m^3)'])
    with open(csvfile,"a") as output:
        writer = csv.writer(output,lineterminator='\n')
        i=1
        dist=0
        for value in list:
            writer.writerow([i,dist,value])
            i+=1
            dist+=equidistance

# ---------------------------------------- Change Polygon Colors -----------------------------------------
def ChangeColorPolygons(csvfile,polygons):
    list_capacity = []
    list_values = [0]
    with open(csvfile,"r") as f:
        csv_reader = csv.reader(f, delimiter=',')
        for row in csv_reader:
            capacity = row[2]
            try:
                list_capacity.append(float(capacity))
            except ValueError:
                list_capacity.append(-9999)
    
    for i in range (0,len(list_capacity)-1):
        if list_capacity[i] == -9999 or list_capacity[i+1] == -9999:
            value = 0
        else:
            perc_change = ((list_capacity[i+1]-list_capacity[i])/list_capacity[i])*100
            if perc_change < -50.0:
                value = 2
            elif 0 > perc_change > -50.0:
                value = 1
            elif perc_change >= 0:
                value = 0
        list_values.append(value)
    print (len(list_values))
    
    lyr = QgsVectorLayer(polygons,'polygons',"ogr")
    provider = lyr.dataProvider()
    provider.addAttributes([QgsField("risk_level", QVariant.Int)])
    lyr.updateFields()
    lyr.startEditing()
    i = 0
    for features in lyr.getFeatures():
        level_index = list_values[i]
        id = features.id()
        attr_value = {2:level_index}
        provider.changeAttributeValues({id:attr_value})
        i += 1
    lyr.commitChanges()
    risk_level = {"0":(184,246,173,"low"),"1":(246,213,120,"medium"),"2":(246,113,103,"high")}
    categories = []
    for zones, (R,G,B,label) in risk_level.items():
        sym = QgsSymbolV2.defaultSymbol(lyr.geometryType())
        sym.setColor(QColor(R,G,B))
        category = QgsRendererCategoryV2(zones,sym,label)
        categories.append(category)
    field = "risk_level"
    renderer = QgsCategorizedSymbolRendererV2(field, categories)
    lyr.setRendererV2(renderer)
    QgsMapLayerRegistry.instance().addMapLayer(lyr)

# ---------------------------------------- Change Centerline Colors -----------------------------------------
def ChangeColorCenterline(csvfile,centerline):
    list_azim = []
    list_azim_values = [0]
    with open(csvfile,"r") as f:
        csv_reader = csv.reader(f, delimiter=',')
        for row in csv_reader:
            azimuth = row[3]
            try:
                list_azim.append(float(azimuth))
            except ValueError:
                continue
        for i in range (0,len(list_azim)-1):
            dir_change = (math.fabs(list_azim[i+1]-list_azim[i])/math.fabs(list_azim[i]))*100
            if dir_change > 100.0:
                value = 2
            elif 100.0 > dir_change > 80.0:
                value = 1
            elif dir_change < 80.0:
                value = 0
            list_azim_values.append(value)
    print (len(list_azim_values))
    lyr = QgsVectorLayer(centerline,'centerline',"ogr")
    provider = lyr.dataProvider()
    provider.addAttributes([QgsField("risk_level", QVariant.Int)])
    lyr.updateFields()
    lyr.startEditing()
    i = 0
    for features in lyr.getFeatures():
        level_index = list_azim_values[i]
        id = features.id()
        attr_value = {4:level_index}
        provider.changeAttributeValues({id:attr_value})
        i += 1
    lyr.commitChanges()
    risk_level = {"0":(96,156,82,"low"),"1":(197,108,6,"medium"),"2":(190,5,5,"high")}
    categories = []
    for zones, (R,G,B,label) in risk_level.items():
        sym = QgsSymbolV2.defaultSymbol(lyr.geometryType())
        sym.setColor(QColor(R,G,B))
        sym.setWidth(0.66)
        category = QgsRendererCategoryV2(zones,sym,label)
        categories.append(category)
    field = "risk_level"
    renderer = QgsCategorizedSymbolRendererV2(field, categories)
    lyr.setRendererV2(renderer)
    QgsMapLayerRegistry.instance().addMapLayer(lyr)


# ------------------------------------- Channel Morphology Analysis --------------------------------------
def ChannelMorphologyAnalysis(proc_path,dem_name,dem_source,year_dem,valley_edge_left,valley_edge_right,equidistance,display,n):
    if not os.path.exists(proc_path+'/'+str(n)):
        os.mkdir(proc_path+'/'+str(n))

    dem_pre_erupt = QgsRasterLayer(proc_path+dem_name,dem_name)
    dem_res = dem_pre_erupt.rasterUnitsPerPixelX()
    valley_W = QgsVectorLayer(valley_edge_left,'valley_W',"ogr")
    valley_E = QgsVectorLayer(valley_edge_right,'valley_E',"ogr")
    epsg = dem_pre_erupt.crs().postgisSrid()
    processing.runalg("qgis:definecurrentprojection",valley_E,"epsg:"+str(epsg))
    processing.runalg("qgis:definecurrentprojection",valley_W,"epsg:"+str(epsg))

    if display == 'True':
        # Load layers to map
        iface.addVectorLayer(valley_edge_left,'valley_W_'+str(n),"ogr")
        iface.addVectorLayer(valley_edge_right,'valley_E_'+str(n),"ogr")

    ################################### Create transects across valley ####################################

    # Create equidistant points along one valley line
    for feature in valley_W.getFeatures():
        length = feature.geometry().length()

    processing.runalg("qgis:createpointsalonglines",valley_W,equidistance,0,length,proc_path+'/'+str(n)+'/pts_along_valleyW')
    pts_along_W = QgsVectorLayer(proc_path+'/'+str(n)+'/pts_along_valleyW.shp','pts_along_W',"ogr")

    # create perpendicular lines (i.e. transects) from points on one valley side to the opposite valley edge:
    uri = "LineString?crs=epsg:" +str(epsg)+"&field=id:integer""&field=distance:double(20,2)&index=yes"
    transects = QgsVectorLayer(uri,
                               "transects",
                               "memory")
    prov = transects.dataProvider()
    line_feature = [ line_feat for line_feat in valley_E.getFeatures() ]
    points_feature = [ pt_feature for pt_feature in pts_along_W.getFeatures() ]
    feats = []
    
    for p in points_feature:
        feat = QgsFeature(transects.fields())
        minDistPoint = min([l.geometry().closestSegmentWithContext(p.geometry().asPoint() ) for l in line_feature])[1]
        feat.setGeometry(QgsGeometry.fromPolyline([p.geometry().asPoint(), minDistPoint]))
        feat.setAttributes([points_feature.index(p),feat.geometry().length()])
        feats.append(feat)
    
    prov.addFeatures(feats)
    
    crs = QgsCoordinateReferenceSystem("epsg:"+str(epsg))
    transects_shpfile = QgsVectorFileWriter.writeAsVectorFormat(transects, proc_path+'/'+str(n)+'/transects.shp', "UTF-8", crs , "ESRI Shapefile")
    transects_shp = QgsVectorLayer(proc_path+'/'+str(n)+'/transects.shp','transects',"ogr")
    
    
    ################################## Create centerline of the valley #####################################
    
    # Use the transects to extract mid-points:
    def mid (pt1, pt2):
        x = (pt1.x() + pt2.x())/2
        y = (pt1.y() + pt2.y())/2
        return QgsPoint(x,y)
    
    def pair(list):
        ''' iterate over pairs in a list '''
        for i in range (1, len(list)):
            yield list[i-1], list[i]
    
    def create_geometry(point,pr):
        seg = QgsFeature()
        seg.setGeometry(QgsGeometry.fromPoint(point))
        pr.addFeatures( [seg] )
    
    # create memory layer for mid-points
    uri = "MultiPoint?crs=epsg:" +str(epsg)+"&field=id:integer""&field=X:double(20,2)&field=Y:double(20,2)&index=yes"
    pt_layer = QgsVectorLayer(uri,
                              "pt_layer",
                              "memory")
    pr = pt_layer.dataProvider()
    
    for elem in transects.getFeatures():
        line = elem.geometry()
        i = int(elem.id())
        for seg_start, seg_end in pair(line.asPolyline()):
            line_start = QgsPoint(seg_start)
            line_end = QgsPoint(seg_end)
            #midpoint
            midpt = mid(line_start, line_end)
            #add midpoint to layer
            create_geometry(midpt,pr)
            pt_layer.updateExtents()
    
    # create a second point layer to add xy coordinates to mid-points
    uri = "MultiPoint?crs=epsg:" +str(epsg)+"&field=id:integer""&field=X:double(20,2)&field=Y:double(20,2)&index=yes"
    midpoints = QgsVectorLayer(uri,
                               "midpoints",
                               "memory")
    prov = midpoints.dataProvider()
    
    pnts = [ pts for pts in pt_layer.getFeatures() ]
    feats = []
    len_list = int(pt_layer.featureCount())
    
    for i in range(0,len_list):
        fet = QgsFeature(midpoints.fields())
        gPnt = QgsGeometry.fromPoint(QgsPoint(pnts[i].geometry().asPoint()))
        fet.setGeometry(gPnt)
        fet.setAttributes([pnts[i].id(),fet.geometry().asPoint()[0],fet.geometry().asPoint()[1]])
        feats.append(fet)
    
    prov.addFeatures(feats)
    
    crs = QgsCoordinateReferenceSystem("epsg:"+str(epsg))
    midpoints_shpfile = QgsVectorFileWriter.writeAsVectorFormat(midpoints, proc_path+'/'+str(n)+'/midpoints.shp', "UTF-8", crs , "ESRI Shapefile")
    midpoints_shp = QgsVectorLayer(proc_path+'/'+str(n)+'/midpoints.shp','midpoints',"ogr")
    
    # Add z values to mid-points
    processing.runalg("saga:addgridvaluestoshapes",midpoints_shp,proc_path+dem_name,0,0,proc_path+'/'+str(n)+'/midpoints_z')
    midpoints_z = QgsVectorLayer(proc_path+'/'+str(n)+'/midpoints_z.shp','midpoints_z',"ogr")
    
    # Create centerline from mid-points and save to memory
    uri = "LineString?crs=epsg:" +str(epsg)+"&field=id:integer""&field=distance:double(20,2)&field=pt_start:double(20,2)&field=pt_end:double(20,2)&index=yes"
    centerline = QgsVectorLayer(uri,
                                "centerline",
                                "memory")
    provider = centerline.dataProvider()
    pts = [ pt_feature for pt_feature in midpoints.getFeatures() ]
    pts_z = [ pt_z for pt_z in midpoints_z.getFeatures() ]
    len_list = int(midpoints.featureCount())
    feats = []
    
    for i in range(0,(len_list-1)):
        feat = QgsFeature(centerline.fields())
        line_start = QgsPoint(pts[i].geometry().asPoint())
        line_end = QgsPoint(pts[i+1].geometry().asPoint())
        feat.setGeometry(QgsGeometry.fromPolyline([line_start,line_end]))
        feat.setAttributes([pts[i].id(),feat.geometry().length(),pts_z[i].attributes()[3],pts_z[i+1].attributes()[3]])
        feats.append(feat)
    
    provider.addFeatures(feats)
    
    # Write and save centerline to permanent shapefile
    centerline_shpfile = QgsVectorFileWriter.writeAsVectorFormat(centerline, proc_path+'/'+str(n)+'/centerline.shp', "UTF-8", crs , "ESRI Shapefile")
    centerline_shp = QgsVectorLayer(proc_path+'/'+str(n)+'/centerline.shp','centerline',"ogr")
 
    ############################ Create Cross-Profiles along valley centerline #############################    

    # Use centerline to create cross profiles across the valley
    length_line = 0.0
    for line in centerline_shp.getFeatures():
        length_line = line.geometry().length() + length_line
    
    trans_len_lst = []
    for transect in transects_shp.getFeatures():
        trans_len = transect.geometry().length()
        trans_len_lst.append(trans_len)
    
    profile_length = np.round(max(trans_len_lst),0)+10.0
    
    processing.runalg("saga:crossprofiles",dem_pre_erupt,centerline_shp,length_line,profile_length,10,0,proc_path+'/'+str(n)+'/cross_profiles')
    cross_profiles = QgsVectorLayer(proc_path+'/'+str(n)+'/cross_profiles.shp','x-profiles',"ogr")
    processing.runalg("saga:profilesfromlines",dem_pre_erupt,None,cross_profiles,"distance",False,0,proc_path+'/'+str(n)+'/x-section',proc_path+'/'+str(n)+'/x-sections')
    profiles = QgsVectorLayer(proc_path+'/'+str(n)+'/x-section.shp','x-sections',"ogr")
    
    ############################### Compute Channel Cross-Sectional Areas ##################################
    # create empty csv file to save results
    csvfile = (proc_path+"/results/results_x-sections_"+str(n)+"_"+year_dem+".csv")
    if os.path.exists(csvfile):
        os.remove(csvfile)
    #Run
    ComputeCrossSectionalAreas(csvfile,profiles)

    ################################ Compute Slope and Bearing Centerline ##################################    
    # create empty csv file to save results
    csvfile = (proc_path+"/results/results_slope_and_bearing_"+str(n)+"_"+year_dem+".csv")
    if os.path.exists(csvfile):
        os.remove(csvfile)
    #Run
    ComputeChannelSlopeAndBearing(csvfile,centerline)

    ############################### Create Polygons for Channel Capacity ###################################

    # Extract transect nodes to create polygons
    processing.runalg("qgis:extractnodes",transects_shp,proc_path+'/'+str(n)+'/nodes')
    nodes = QgsVectorLayer(proc_path+'/'+str(n)+'/nodes.shp','nodes',"ogr")

    #add XY coordinates to each feature in nodes.getFeatures
    processing.runalg("saga:addcoordinatestopoints",nodes,proc_path+'/'+str(n)+'/nodesXY')
    nodesXY = QgsVectorLayer(proc_path+'/'+str(n)+'/nodesXY.shp','nodesXY',"ogr")
    len_list = int(nodesXY.featureCount())
    uri = "MultiPoint?crs=epsg:" +str(epsg)+"&field=id:integer""&field=X:double(20,2)&field=Y:double(20,2)&index=yes"
    pts_W = QgsVectorLayer(uri,
                           "pts_W",
                           "memory")
    prov = pts_W.dataProvider()
    pnts = [ pts for pts in nodesXY.getFeatures() ]
    feats = []

    for i in range(0,(len_list-1),2):
        fet = QgsFeature(pts_W.fields())
        gPnt = QgsGeometry.fromPoint(QgsPoint(pnts[i].geometry().asPoint()))
        fet.setGeometry(gPnt)
        fet.setAttributes([pnts[i].id(),fet.geometry().asPoint()[0],fet.geometry().asPoint()[1]])
        feats.append(fet)

    prov.addFeatures(feats)

    uri = "MultiPoint?crs=epsg:" +str(epsg)+"&field=id:integer""&field=X:double(20,2)&field=Y:double(20,2)&index=yes"
    pts_E = QgsVectorLayer(uri,
                           "pts_E",
                           "memory")
    prov = pts_E.dataProvider()
    pnts = [ pts for pts in nodesXY.getFeatures() ]
    feats = []

    for i in range(1,(len_list-1),2):
        fet = QgsFeature(pts_E.fields())
        gPnt = QgsGeometry.fromPoint(QgsPoint(pnts[i].geometry().asPoint()))
        fet.setGeometry(gPnt)
        fet.setAttributes([pnts[i].id(),fet.geometry().asPoint()[0],fet.geometry().asPoint()[1]])
        feats.append(fet)

    prov.addFeatures(feats)

    # Create polygons
    uri = "MultiPolygon?crs=epsg:" +str(epsg)+"&field=id:integer""&field=area:double(20,2)&index=yes"
    polygons = QgsVectorLayer(uri,
                              "polygons",
                              "memory")
    prov = polygons.dataProvider()
    pnts_W = [ pt_feature for pt_feature in pts_W.getFeatures() ]
    pnts_E = [ pt_feature for pt_feature in pts_E.getFeatures() ]
    feats = []
    count = pts_W.featureCount()
    len_list = int(count)

    for i in range(0,(len_list-2),1):
        fet = QgsFeature(polygons.fields())
        gPnt1 = (QgsPoint(pnts_W[i].geometry().asPoint()))
        gPnt2 = (QgsPoint(pnts_E[i].geometry().asPoint()))
        gPnt3 = (QgsPoint(pnts_E[i+1].geometry().asPoint()))
        gPnt4 = (QgsPoint(pnts_W[i+1].geometry().asPoint()))
        gPolygon = QgsGeometry.fromPolygon([[gPnt1, gPnt2, gPnt3, gPnt4]])
        fet.setGeometry(gPolygon)
        fet.setAttributes([fet.id(),fet.geometry().area()])
        feats.append(fet)

    prov.addFeatures(feats)

    polygons.startEditing()
    fid = polygons.fieldNameIndex('id')
    for feat in polygons.getFeatures():
        polygons.changeAttributeValue(feat.id(), fid, feat.id())

    polygons.commitChanges()
    crs = QgsCoordinateReferenceSystem("epsg:"+str(epsg))
    polygon_shpfile = QgsVectorFileWriter.writeAsVectorFormat(polygons, proc_path+'/'+str(n)+'/polygons.shp', "UTF-8", crs , "ESRI Shapefile")
    polygons_shp = QgsVectorLayer(proc_path+'/'+str(n)+'/polygons.shp','polygons_'+str(n),"ogr")

    # Convert polygon vertices to points
    processing.runalg("saga:convertpolygonlineverticestopoints",polygons_shp,proc_path+'/'+str(n)+'/vertices')
    vertices = QgsVectorLayer(proc_path+'/'+str(n)+'/vertices.shp','vertices',"ogr")

    # Add z values to nodes
    processing.runalg("saga:addgridvaluestopoints",vertices,proc_path+dem_name,0,0,proc_path+'/'+str(n)+'/nodes_xyz')
    nodesXYZ = QgsVectorLayer(proc_path+'/'+str(n)+'/nodes_xyz.shp','nodes_xyz',"ogr")

    ######################## Create synthetic surface at the top of the channel ############################

    # IDW interpolation of a synthetic surface using nodes
    fields = nodesXYZ.pendingFields()
    field_name = str(fields.field(6).name())
    processing.runalg("grass7:v.surf.idw",nodesXYZ,12,2,field_name,True,"%f,%f,%f,%f"%(nodesXYZ.extent().xMinimum(),nodesXYZ.extent().xMaximum(),nodesXYZ.extent().yMinimum(),nodesXYZ.extent().yMaximum()),dem_res,-1,0.0001,proc_path+'/'+str(n)+'/nodes2idw')
    nodes2idw = QgsRasterLayer(proc_path+'/'+str(n)+'/nodes2idw.tif','nodes2idw')

    # Clip the synthetic surface and input DEM to channel extent
    processing.runalg("gdalogr:cliprasterbyextent",nodes2idw,-9999,"%f,%f,%f,%f"%(nodesXYZ.extent().xMinimum(),nodesXYZ.extent().xMaximum(),nodesXYZ.extent().yMinimum(),nodesXYZ.extent().yMaximum()),5,0,75,6,1,False,0,True,None,proc_path+'/'+str(n)+'/nodes2idw_clip_extent')
    nodes2idw_clip = QgsRasterLayer(proc_path+'/'+str(n)+'/nodes2idw_clip_extent.tif','nodes2idw_clip')
    processing.runalg("gdalogr:cliprasterbyextent",dem_pre_erupt,-9999,"%f,%f,%f,%f"%(nodesXYZ.extent().xMinimum(),nodesXYZ.extent().xMaximum(),nodesXYZ.extent().yMinimum(),nodesXYZ.extent().yMaximum()),5,0,75,6,1,False,0,True,None,proc_path+'/'+str(n)+'/orig_dem_clip_extent')
    dem_pre_clip = QgsRasterLayer(proc_path+'/'+str(n)+'/orig_dem_clip_extent.tif','orig_dem_clip')

    # Resample interpolated surface to the same resolution as input DEM
    processing.runalg("grass7:r.resample",nodes2idw_clip,"%f,%f,%f,%f"%(nodes2idw_clip.extent().xMinimum(),nodes2idw_clip.extent().xMaximum(),nodes2idw_clip.extent().yMinimum(),nodes2idw_clip.extent().yMaximum()),dem_res,proc_path+'/'+str(n)+'/nodes2idw_clip_extent_resamp')

    ####################### Compute raster difference b/w top and channel surface ##########################
    
    # Compute raster difference between clipped input DEM and synthetic surface
    real_dem = gdal.Open(proc_path+'/'+str(n)+'/orig_dem_clip_extent.tif')
    synthetic_surf = gdal.Open(proc_path+'/'+str(n)+'/nodes2idw_clip_extent_resamp.tif')
    # Get array dimensions
    pixels = real_dem.RasterXSize
    lines = real_dem.RasterYSize
    # Get band GDAL
    real_dem_band = real_dem.GetRasterBand(1)
    synthetic_surf_band = synthetic_surf.GetRasterBand(1)
    # Get band GDAL datatype and map to Numpy Array
    gdal_band_type = real_dem_band.DataType
    np_band_type = gdal_array.GDALTypeCodeToNumericTypeCode(gdal_band_type)
    # Read the DEM rasters into arrays
    real_dem_array = real_dem_band.ReadAsArray(0,0,pixels,lines).astype(np_band_type)
    synthetic_surf_array = synthetic_surf_band.ReadAsArray(0,0,pixels,lines).astype(np_band_type)
    # Get Geotransform infos
    geotrans = real_dem.GetGeoTransform()
    projInfo = real_dem.GetProjection()

    # Create output and apply projection and transform info
    driver = gdal.GetDriverByName("GTiff")
    dst_ds = driver.Create(proc_path+'/'+str(n)+'/diff_dem.tif',pixels,lines,1,gdal.GDT_Float32)
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
            diff_dem[i,j] = synthetic_surf_array[i,j] - real_dem_array[i,j]

    # Write diff_dem array to band
    dst_ds.GetRasterBand(1).WriteArray(diff_dem)
    dst_ds = driver.Create(proc_path+'/'+str(n)+'/diff_dem.tif',pixels,lines,1,gdal.GDT_Float32)
    dst_ds.SetGeoTransform(geotrans)
    sp_ref = osr.SpatialReference()
    sp_ref.ImportFromWkt(projInfo)
    dst_ds.SetProjection(sp_ref.ExportToWkt())
    QgsRasterLayer(proc_path+'/'+str(n)+'/diff_dem.tif','diff_dem')

    # Close band and datasets
    real_dem = None
    real_dem_band = None
    synthetic_surf = None
    synthetic_surf_band = None

    # Fill data gaps
    os.system('gdal_fillnodata.bat -md 100 -b 1 -of GTiff ' + proc_path+'/'+str(n)+'/diff_dem.tif'+' '+ proc_path+'/'+str(n)+'/diff_dem_fill_nodata.tif')
    QgsRasterLayer(proc_path+'/'+str(n)+'/diff_dem_fill_nodata.tif','diff_dem_fillnodata')
    
    # Compute volume for each cell
    os.system('gdal_calc --calc "A*'+str(dem_res)+'" --format GTiff --type Float32 -A '+proc_path+'/'+str(n)+'/diff_dem_fill_nodata.tif --A_band 1 --outfile '+proc_path+'/'+str(n)+'/diff_dem_multiply.tif')
    QgsRasterLayer(proc_path+'/'+str(n)+'/diff_dem_multiply.tif','diff_dem_multiply')
    
    # Get total volume for each polygon
    processing.runalg("qgis:zonalstatistics",proc_path+'/'+str(n)+'/diff_dem_multiply.tif',1,proc_path+'/'+str(n)+'/polygons.shp','_',True,proc_path+'/'+str(n)+'/zonal_statistics')
    zonal_stats = QgsVectorLayer(proc_path+'/'+str(n)+'/zonal_statistics.shp','zonal_stats',"ogr")

    # create empty csv file to store results
    csvfile = proc_path+"/results/results_capacity_"+str(n)+"_"+year_dem+".csv"
    if os.path.exists(csvfile):
        os.remove(csvfile)

    #Run
    ComputeChannelCapacity(csvfile, zonal_stats)
    ChangeColorPolygons(proc_path+"/results/results_capacity_"+str(n)+"_"+year_dem+".csv",proc_path+'/'+str(n)+'/polygons.shp')
    ChangeColorCenterline(proc_path+"/results/results_slope_and_bearing_"+str(n)+"_"+year_dem+".csv",proc_path+'/'+str(n)+'/centerline.shp')

# ------------------------------------------- Compute Channel Geomorphology Analysis -----------------------------------------
# Get number of valleys in folder
tot_files = len(glob.glob1(proc_path+'/valley_edges',"*.shp"))
total_valleys = tot_files/2

start_time = time.time()
for n in range (1,total_valleys+1):
    valley_edge_left = proc_path+'/valley_edges/valley_edge_W_'+str(n)+'.shp'
    valley_edge_right = proc_path+'/valley_edges/valley_edge_E_'+str(n)+'.shp'
    #Run
    ChannelMorphologyAnalysis(proc_path,dem_name,dem_source,year_dem,valley_edge_left,valley_edge_right,equidistance,display,n)

elapsed_time = time.time() - start_time
minutes = (np.round(elapsed_time))/60
print 'Elapsed time after processing is %s minutes: ' % (minutes)
