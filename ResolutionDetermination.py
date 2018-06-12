# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 15:20:01 2018

@author: J Garrett

This code is a way to determine the resolution of a (KPFM) image.

"""

import CypherKPFM as ck
import matplotlib
from matplotlib import pyplot as plt
from scipy import signal
from scipy import ndimage
from scipy.optimize import curve_fit
import numpy as np
from skimage.morphology import watershed
from skimage import filters
from skimage import segmentation
from skimage import measure
from skimage.feature import register_translation
from scipy.ndimage import fourier_shift
import csv
import os

	
filepath1 = "L:\\JGarrett\\AsylumResearchData\\171229\\FLG1H_3V0000.ibw"
filepath2 = "L:\\JGarrett\\AsylumResearchData\\171229\\FLG0H_3V0000.ibw"

baseotherfilepath = "L:\\JGarrett\\AsylumResearchData\\160721_FLG_HO_data"

def populateDict( alist ):
    myDict = {}    
    for i in alist:
        im = ck.load_and_prep(i)
        ck.avg_kpfms(im)
        ck.clean_kpfms(im)
        lab = os.path.splitext(i)[0]
        
        print(lab)        
        myDict[lab] = {'KPFM':'Cleaned_KPFM','data':im }
    return myDict

def plotFromDict( myDict, key):
    ck.quickplot(ck.get_image(myDict[key]['data'],myDict[key]['KPFM']))
    print(key)
    return myDict[key]['data']

def load2params( path, b_pts = False):
    scan = ck.load_and_prep(path)
    kpfm = ck.avg_kpfms(scan)
    b_list, b_dist, segments = getBoundaries(kpfm)
    if b_pts == False:
        b_pts = plotwboundary(kpfm, b_list[0], num = 50)
    
    vCuts, vGroup, vSelect = manyVoltages(b_pts, kpfm, b_dist, 35, segments)
    all_params = parametersManyPoints(b_pts, vCuts)
    return all_params, float(scan['cleannote']['dx'])
    
def prepped2Params( data ):
    kpfm = ck.get_image(data,"Cleaned_KPFM")
    b_list, b_dist, segments = getBoundaries(kpfm)
    b_pts = plotwboundary(kpfm, b_list[0], num = 10)
    
    vCuts, vGroup, vSelect = manyVoltages(b_pts, kpfm, b_dist, 50, segments)
    all_params, all_errs = parametersManyPoints(b_pts, vCuts)
    return all_params, all_errs, vCuts
    
def shiftCorrection( original, new, toReturn = 'Averaged_KPFM'):
    oTopo = ck.basicTopography(original)
    nTopo = ck.basicTopography(new)
    
    shift, error, diffphase = register_translation(oTopo, nTopo)
    shiftedData = ck.get_image(new, toReturn)
    shiftedData = fourier_shift(np.fft.fftn(shiftedData), shift)
    shiftedData = np.fft.ifftn(shiftedData)
    return shiftedData.real, shift

def allParams( dataDict, points, dist_image, max_dist, segments ):
    keys = list(dataDict.keys())
    base = 'FLG1H_3V0000'
    for i in keys:
        print('running {}'.format(i))
        shiftedKPFM, shift = shiftCorrection(dataDict[base]['data'],
                                             dataDict[i]['data'],
                                             toReturn = dataDict[i]['KPFM'])
        print('finished shifting {}'.format(i))
        vCuts, vGroup, sel = manyVoltages(points, shiftedKPFM,
                                          dist_image, max_dist, segments)
        print('cuts in {} determined'.format(i))
        all_params, all_errs = parametersManyPoints(points, vCuts)
        print('finished fitting {}'.format(i))
        dataDict[i]['params'] = all_params
        dataDict[i]['uncertainty'] = all_errs
        dataDict[i]['shift'] = shift
        dataDict[i]['cuts'] = vCuts
    return 0
    
def allParamsSeparate( dataDict ):
    keys = list(dataDict.keys())
    for i in keys:
        print('running {}'.format(i))
        all_params, all_errs, vCuts = prepped2Params(dataDict[i]['data'])
        print('finished fitting {}'.format(i))
        dataDict[i]['params'] = all_params
        dataDict[i]['uncertainty'] = all_errs
        dataDict[i]['cuts'] = vCuts
    return 0
        
def saveAllParams(dataDict):
    keys = list(dataDict.keys())
    keys.sort()
    base = 'FLG1H_3V0000'
    try:     
        points = list(dataDict[base]['params'].keys())
        points.sort()
        with open('resPoints.csv', 'w') as myfile:
            wr = csv.writer(myfile)
            for row in points:
                print(row)
                wr.writerow(row)
        dim2 = len(points)
    except KeyError:
        print("inconsistent points")
        k = list(dataDict.keys())[0]
        dim2 = len(list(dataDict[k]['params'].keys()))        
        
    all_params = np.zeros((len(keys),dim2,4))  
    all_errs = np.zeros((len(keys),dim2,4)) 
    all_chisq = np.zeros((len(keys),dim2))      
    for i in range(len(keys)):
        k = keys[i]
        for j in range(dim2):
            #print(i,j,k,points[j])
            points = list(dataDict[k]['params'].keys())
            all_params[i][j][:] = dataDict[k]['params'][points[j]][0]
            all_chisq[i][j] = dataDict[k]['params'][points[j]][1]
            all_errs[i][j][:] = dataDict[k]['uncertainty'][points[j]]
        
    names = ['DeltaV', 'a_fit','offset','x0']    
    for i in range(4):
        np.savetxt(names[i]+'.csv', all_params[:,:,i], delimiter = ',')
        np.savetxt(names[i]+'_uncertainty.csv', all_errs[:,:,i], delimiter = ',')
        
    np.savetxt("reducedChiSq.csv", all_chisq, delimiter = ',')
    return all_params

def prepForWatershed( ibw ):
    try: 
        kpfm = ck.get_image(ibw, 'Averaged_KPFM')
    except KeyError:
        try:
            kpfm = ck.get_image(ibw, 'UserIn1Retrace')
        except KeyError:
            return -1
    kpfm = signal.medfilt2d( kpfm, 9) 
    return kpfm

def paddedmedfilt( image, dim ):
    temp = np.pad(image, dim, 'edge')
    temp = signal.medfilt2d(temp, 2*dim+1)
    filtered = np.array(temp[dim:-dim,dim:-dim])
    return filtered

def getBoundaries( kpfm_image , lower = 0.7, upper = 1.4):
    kpfm_mask = ck.autoMask(kpfm_image, toMask = 'low')
    kpfm_flat = ck.flatten(kpfm_image, mask = kpfm_mask)
    padded = paddedmedfilt(kpfm_flat,15)
    kpfmn = ck.normalize(kpfm_image)
    val = filters.threshold_otsu(kpfmn)
    lines5 = filters.sobel(padded)
    
    markers = np.zeros_like(kpfm_image)
    markers[kpfmn < lower*val] = 1
    markers[kpfmn > upper*val] = 2
    
    segments = watershed(lines5, markers)
    boundary = segmentation.find_boundaries(segments)
    blist = measure.find_contours(segments, 1.5)
    boundary_dist = ndimage.distance_transform_edt(1-boundary)
    return blist, boundary_dist, segments

def plotwboundary( image, boundary, num = 0):
    plt.imshow(np.transpose(image))
    plt.plot(np.transpose(boundary)[0], np.transpose(boundary)[1], 
                 c='r')
    if num > 0:
        spacing = len(boundary)/num
        points = np.array([ boundary[int(spacing*(i+1/2))] \
            for i in range(int(num))])
        plt.scatter(np.transpose(points)[0],
                    np.transpose(points)[1],zorder=3, c = 'black')
        #plt.scatter(boundary[0][0])
        plt.show()
        return points
        
    return 0

def grabNearbyBoundary( point, image, dist_image,
                       max_dist, segmentations, lateral = 20):
    selected_points = {i:[] for i in range(int(max_dist))}
    voltages = {i:[] for i in range(int(2*max_dist))}
    mixedimage = np.zeros_like(image)
    seps = separations((point[1],point[0]),image)
    #plt.imshow(seps)
    #plt.plot(point[1],point[0], marker='o', markersize=3, color="red")
    #plt.show()
    #plt.imshow(image)
    #plt.show()
    
    for i in range(int(max_dist)):
        dlim = np.sqrt((i+1)**2+lateral**2)
        selected_points[i] = \
            (np.abs(dist_image-i)<=(0.5))&(seps<dlim)&(segmentations==1)
        voltages[int(max_dist)+i] = image[selected_points[i]]
                #      &(seps<dlim)
            #[image[(np.abs(dist_image-0.5)<=(i+0.5))\
             #      &(seps<dlim)&segmentations==1],
            #image[(np.abs(dist_image-0.5)<=(i+0.5))\
             #     &(seps<dlim)&segmentations==2]]
        mixedimage[selected_points[i]] = i 
    for i in range(int(max_dist)):
        dlim = np.sqrt((i+1)**2+lateral**2)
        selected_points[i] = \
            (np.abs(dist_image-i)<=(0.5))&(seps<dlim)&(segmentations==2)
        voltages[int(max_dist)-i] = image[selected_points[i]]
                #      &(seps<dlim)
            #[image[(np.abs(dist_image-0.5)<=(i+0.5))\
             #      &(seps<dlim)&segmentations==1],
            #image[(np.abs(dist_image-0.5)<=(i+0.5))\
             #     &(seps<dlim)&segmentations==2]]
        mixedimage[selected_points[i]] = -i 
        
        
    
    voltage_points = np.concatenate([lineScatter(voltages[i], i) \
                      for i in voltages.keys()])
    
    voltage_avgs = np.array([[np.nanmean(voltages[i]),i] \
                      for i in voltages.keys()])
    
    return mixedimage, voltage_points, np.transpose(voltage_avgs)

def separations(point, image):
    #seps = np.array([[np.linalg.norm(point-np.array([j,i])) \
                    #for i in range(image.shape[0])] \
                    #for j in range(image.shape[1])])
    locs = np.array([[[i,j] for i in range(image.shape[1])]\
             for j in range(image.shape[0])])
    
    seps = np.array([[np.linalg.norm(point-i) for i in j]\
                     for j in locs])
    return seps
        
def voltagesPlot( voltages ):
    return 0

def lineScatter(ys, x):
    points = np.array([ys,np.ones_like(ys)*x])
    return np.transpose(points)

def lineAverage():
    return 0

def manyVoltages(points, image, dist_image, max_dist, segmentations):
    '''return voltage_cuts, voltage_pts_grouped, selected_pts'''
    voltage_cuts = {tuple(i):[] for i in points}
    voltage_pts_grouped = {tuple(i):[] for i in points}
    selected_pts = np.zeros_like(image)
    
    for i in points:
        print('Voltage cut with point {}'.format(str(i)))
        t0, voltage_pts_grouped[tuple(i)], voltage_cuts[tuple(i)] = \
            grabNearbyBoundary( i, image,
            dist_image, max_dist, segmentations)
        selected_pts = selected_pts + t0
            
    return voltage_cuts, voltage_pts_grouped, selected_pts

def qpCuts( cuts, points, i):
    try:
        for j in i:
            p = cuts[tuple(points[j])]
            plt.plot(p[1],p[0])
    except TypeError:
        p = cuts[tuple(points[i])]
        plt.plot(p[1],p[0])
    
    plt.show()
    return 0

def tfunc(x, Vb, a, b, x0):
        return Vb*(np.tanh(2*a*(x-x0)))/2+b
    
def dropnans(points):
    '''Drops any 'nan' y-values from a list of points'''
    mask = np.ones_like(points, dtype = bool)
    for i in points:
        #print(i)
        if np.isnan(i[0]):
            mask[int(i[1])] = 0
            
    remaining = int(np.sum(mask)/2)
    newpoints = points[mask] 
    return np.transpose(newpoints.reshape((remaining,2)))
    
def edgeFit( cut ):
    cleancut = dropnans(np.transpose(cut))
    L = len(cleancut)
    try:    
        popt, pcov = curve_fit(tfunc, cleancut[1], cleancut[0],
                                p0 = [-.3,2e-2,0.5,35])    
        redChiSq = np.sum((tfunc(cleancut[1],*popt)-cleancut[0])**2)/L
        return popt, redChiSq, pcov
    except RuntimeError:
        return [0,0,0,0], 0, [0]

def edgePlot( cut, popt ):
    plt.plot(cut[1], cut[0], 'bo')    
    plt.plot(cut[1], tfunc(cut[1],*popt), 'r-')    
    return 0
    
def parametersManyPoints(points, cuts):
    all_params = {tuple(i):[] for i in points}    
    all_errs = {tuple(i):[] for i in points}  
    
    for i in points:
        tofit = cuts[tuple(i)]
        popt, redChiSq, pcov = edgeFit(tofit)
        all_params[tuple(i)] = [popt, redChiSq]
        all_errs[tuple(i)] = np.sqrt(np.diag(pcov))
        
    return all_params, all_errs
    
def plotHist( arrayDict, i ):
    length = len(arrayDict)
    toPlot = np.zeros((length,1))
    keys = list(arrayDict.keys())
    for j in range(length):
        toPlot[j] = arrayDict[keys[j]][i]
        
    plt.hist(toPlot)
    
    return toPlot
    
def saveAll( dataDict ):
    keys = list(dataDict.keys())
    for i in keys:
        ck.quickSave(dataDict[i]['data'],dataDict[i]['KPFM'],
                     i)
    
    return 0
#dataList = { name:{'path':'', 'KPFM':'', } for i in range(8)}
def plotFD( dictionary, k1, k2):
    edgePlot(cut,potp)
    return 0