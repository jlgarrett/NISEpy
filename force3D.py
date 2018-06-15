# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 10:08:32 2018

@author: joe
"""
import numpy as np


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
            return -1e12class interpolate_packaging:
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
    nside = max(list(filter_interp.keys()))[0]
    hlog = np.log(h)
    print(hlog)
    lf = interpolate_packaging( filter_interp )
    #print(list(filter_interp.keys()))
    hfilter = ([[np.exp(lf.out(racs((i,j)),hlog)) 
                 for j in range(-nside, nside+1)] for i in range(-nside, nside+1)])
    hfilter = hfilter / np.sum(hfilter)
    return hfilter

def racs( tup ):
    '''Reverse Absolute Coordinate Sort
    
    I wrote this function to simply make line 3 of "make_filter" shorter'''
    absx = abs(tup[0])
    absy = abs(tup[1])
    return (min(absx,absy),max(absx,absy))
  
  