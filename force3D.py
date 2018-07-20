# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 10:08:32 2018

@author: joe
"""
import numpy as np
import tables


def expand_filter( filt ):
    halfside = (filt.shape[0]-1)//2
    print(halfside)
    exp_filt = filt
    for i in range(halfside, 2*halfside+1):
        for j in range(halfside,i):
            exp_filt[i][j] = exp_filt[j][i]
    
    exp_filt[halfside:,:halfside] = exp_filt[halfside:,:halfside:-1]
    exp_filt[:halfside] = exp_filt[:halfside:-1]
            
    return exp_filt
    
class interpolate_packaging:
    def __init__(self, filt):
        #print(start,end)
        self.pfilter = filt
        
    def out(self, loc, x):
        output = self.pfilter[loc](x)
        #print(output, output == np.nan)
        if np.isfinite(output):
            return output
        else:
            return -1e12   
    
def make_filter( h, filter_interp):
    '''Makes the filter from the interpolated filter data'''
    nside = max([i for i in filter_interp.keys() if type(i) is tuple])[0]
    hlog = np.log(h)
    print(hlog)
    lf = interpolate_packaging( filter_interp )
    #print(list(filter_interp.keys()))
    hfilter = ([[np.exp(lf.out(racs((i,j)),hlog)) 
                 for j in range(-nside, nside+1)] for i in range(-nside, nside+1)])
    hfilter = hfilter / np.sum(hfilter)
    return hfilter
  
def make_calc_filters(filter_raw):
    '''Makes the filter from the interpolated filter data'''
    nside = max([i for i in filter_raw.keys() if type(i) is tuple])[0]
    #hs = filter_raw['separation']
    #lf = interpolate_packaging( filter_interp )
    #print(list(filter_interp.keys()))
    hfilters = {}
    for hx, h in enumerate(filter_raw['separation']):
        hfilter = np.array([[filter_raw[racs((j,i))][hx] 
                 for j in range(-nside, nside+1)] for i in range(-nside, nside+1)])
        print(np.sum(hfilter))
        hfilter = hfilter/np.sum(hfilter)
        hfilters[h] = hfilter
    #for i, h in enumerate(filter_raw['separation']):
    #  hfilter[:,:,i] = hfilter[:,:,i] / np.sum(hfilter[:,:,i])
      
    return hfilters

def racs( tup ):
    '''Reverse Absolute Coordinate Sort
    
    I wrote this function to simply make line 3 of "make_filter" shorter'''
    absx = abs(tup[0])
    absy = abs(tup[1])
    return (min(absx,absy),max(absx,absy))
  
def initialize_filtered_images( h_vals, raw_image, filter_int, savingfile = 0):
    '''Converts a cleaned KPFM voltages image into the image charge voltage at 
    certain separations.
    
    h_vals is list of separations.
    raw_image should be a np. array.
    filter_int is a dictionary describing a filter.
    the function returns a dictionary of the voltage from the image charges as several heights, 
    as a dictionary with the heights as the keys.'''
    filtered_images = {}
    starttime = time.time()
    for i in h_vals:
        print(time.time()-starttime)
        loc_filt = make_filter(i,filter_int)
        filtered_images[i] = signal.convolve2d(raw_image, loc_filt,
                                               boundary='symm',mode='same') 
        if savingfile is not 0:
          print(189)
          #savingfile(i) = filtered_images[i]
    
    return filtered_images