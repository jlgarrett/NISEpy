# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 10:07:52 2018

@author: joe
"""

import numpy as np
import scipy
from scipy.integrate import dblquad
import scipy.interpolate
import time
import force3D as f3
from scipy import signal
import mpmath as mp

'''This file is for calculating the patch filters to process an image. Both
the actual processing of KPFM data, and the calculation of the force from the 
images can be found in 'force3D.py'. The function 'completefiltering' creates
the filter, for a given number of pixes, pixel size, heights, when the function
is specified. It acts using the functions preceding it.'''

def interpolate_from_pix( pix ):
    pix2 = [i for i in pix.keys() if not isinstance(i, str)]
    hs = pix['separation']
    interpolations = {j:scipy.interpolate.interp1d(np.log(hs), np.log(pix[j]), kind = 'linear') for j in pix2}
    #nterpolations = {j:scipy.interpolate.interp1d(hs, pix[j], kind = 'cubic') for j in pix2}
    return interpolations    

def patchIntegrand( label = ""):
    '''Calculates the radial integral of the filter to apply to the patch 
    potential image for a particular force setting. x is k normalized by 
    the radius. dor is d normalized by the radius. For some reason, all the
    'pi's are gone.
    
    The label can be:
    
    'fs' for force from same surface
    'fo' for force from opposite surfaces
    'fds' for force derivative from same surface
    'fdo' for force derivative from other surface'''
     
    #Return a function based on the input 
    if label == 'fs':
        def intfunc(x, dor):
            return dor**4 * np.float128(x**3 / (np.sinh(x*dor))**2) * np.float128(scipy.special.j0(x))
        return intfunc
    elif label == 'fo': 
        def intfunc(x, dor):
            return dor**4 * x**3 / (np.sinh(x*dor))**2 * np.cosh(x*dor) \
            * scipy.special.j0(x)
        return intfunc
    elif label == 'fds':
        def intfunc(x, dor):
            return dor**5 * x**4 / (np.sinh(x*dor))**3 * np.cosh(x*dor) \
            * scipy.special.j0(x)
        return intfunc
    elif label == 'fdo':
        def intfunc(x, dor):
            return dor**5 * x**4 / (np.sinh(x*dor))**3 * (np.cosh(x*dor)**2 
            + 1) * scipy.special.j0(x)
        return intfunc
    else:
        return 0
      
def mpPatchIntegrand( label = ""):
    '''Calculates the radial integral of the filter to apply to the patch 
    potential image for a particular force setting. x is k normalized by 
    the radius. dor is d normalized by the radius. For some reason, all the
    'pi's are gone.
    
    The label can be:
    
    'fs' for force from same surface
    'fo' for force from opposite surfaces
    'fds' for force derivative from same surface
    'fdo' for force derivative from other surface
    
    unlike the original PatchIntegrand function, these are funcion of functions'''
     
    #Return a function based on the input 
    if label == 'fs':
        def intfunc(dor):
          def interfunc(x):
            return scipy.constants.pi/2*dor**4 * x**3 / mp.sinh(x*dor)**2 * mp.j0(x)
          return interfunc
        return intfunc
    elif label == 'fo': 
        def intfunc(dor):
          def interfunc(x):
            return scipy.constants.pi/2*dor**4 * x**3 / (mp.sinh(x*dor))**2 * mp.cosh(x*dor) \
            * mp.j0(x)
          return interfunc
        return intfunc
    elif label == 'fds':
      def intfunc(dor):
        def interfunc(x):
            return scipy.constants.pi/2*dor**5 * x**4 / (mp.sinh(x*dor))**3 * mp.cosh(x*dor) \
            * mp.j0(x)
        return interfunc
      return intfunc
    elif label == 'fdo':
      def intfunc(dor):
        def interfunc(x):
            return scipy.constants.pi/4*dor**5 * x**4 / (mp.sinh(x*dor))**3 * (mp.cosh(x*dor)**2 
            + 1) * mp.j0(x)
        return interfunc
      return intfunc
    else:
        return 0

def mpIntegralfunc( dor, force_integ):
    '''The mpmath package is used to arrive to integrate to lower dor than 
    possible with scipy.quad'''
   
    if dor > 100:
        starttime = time.time()
        inte = mp.quad(force_integ(dor),[0,mp.inf])
        print('the duration for dor = {0} is: {1}'.format(dor, time.time()-starttime))
        return float(inte)
    elif dor < 0.05:
        mp.mp.prec = 250
    elif dor < 0.1:
        mp.mp.prec = 100
    #else:
        #mp.mp.prec = 5
    
    j0zero = lambda n: mp.findroot(mp.j0, mp.pi*(n-0.25))
    
    starttime = time.time()
    inte = mp.quadosc(force_integ(dor),[0,mp.inf], zeros = j0zero)
    print('the duration for dor = {0} is: {1}'.format(dor, time.time()-starttime))
    return float(inte)
    
      
def integrate_rad(rfunc, minpoint, maxpoint):
        '''Integrates the function rfunc from minpoint to maxpoint (tuples).
        Includes code to prevent integral from going wild at small values'''
        p1 = minpoint
        p2 = maxpoint
        #If we let the function get too small, the integral never converges
        if abs(rfunc( np.sqrt(p1[0]**2+p2[0]**2) )) < 1e-12: 
            return 0
        else:
            return dblquad(lambda x,y: rfunc( np.sqrt(x**2+y**2)), p1[0], p2[0],
                           lambda x: p1[1], lambda x: p2[1], epsabs=1.49e-08)[0]
        
def define_filter( rfunc , side_n , spacing ):
    '''Calculates the patch potential filter from a radially-symmetric function.
    rfunc is the function in terms of $r$. side_n is the number of pixels per
    side, which must be odd, for symmetry reasons, and spacing is the spacing 
    between the centers of the pixels'''
    if side_n %2 == 0:
        print('Side must have odd number of pixels')
        return 0
    else:
        #We define indices to give the coordinates of all the points
        indices = [[(i,j) for i in range(-(side_n-1)//2,(side_n+1)//2)] \
        for j in range(-(side_n-1)//2,(side_n+1)//2)]
        #only one eighth of the filter is filled out because the rest is 
        #implied by radial symmetry
        rfilter = np.array([[0.0 for i in range(-(side_n-1)//2,(side_n+1)//2)]\
            for j in range(-(side_n-1)//2,(side_n+1)//2)])     
        
        #We fill in the filter by integrating over the radial function
        for i in range((side_n-1)//2,(side_n)):
            for j in range(i,(side_n)):
                p1 = tuple(spacing*indices[i][j][k]-spacing/2 for k in range(2))
                p2 = tuple(spacing*indices[i][j][k]+spacing/2 for k in range(2))
                inte = integrate_rad(rfunc, p1 ,p2)
                rfilter[i][j] = inte
        
        return rfilter

def internal_integral( dor , force_inte = 0, intlimit = 0):
    '''Integrates over the function force_inte. Assumes a force of the type
    generated by patchIntegrand()'''
    if force_inte == 0:
        force_inte = patchIntegrand('fs')
    if intlimit == 0:
      intlimit = 100/dor
    return scipy.integrate.quad(lambda x: force_inte(x,np.float128(dor)),0,intlimit, 
                                epsrel = 1e-12, epsabs = 1e-12, limit = 5000)[0]

def ii_array( dors , force_integ = 0, integral_func = 0):
    '''Creates an array of integrals of a function for many separations. 
    Assumes a function of the form of patchIntegrand()'''
    if force_integ == 0:
        force_integ = patchIntegrand('fs')
    output = dors
    
    
    for i in range(len(dors)):
        output[i] = integral_func( dors[i] , force_integ)
    
    return output

class InterpRes:
    '''This class stores the interpolation of a radial function describing
    the correlation between different points on a patch potential image. It is 
    useful because it stores the data on a log-log scale, if it can, but stores
    it on a log-linear scale if necessary. The keywords are:
    start -- is the smallest x value
    end -- is the largest x value to be interpolated
    func_2_interpolate -- is the function to be interpolated.'''
    def __init__(self, start=0.09, end=10**5,
                 func_2_interpolate=0, integral_func = 0, xs = 0, ys = 0):
        #print(start,end)
        if func_2_interpolate == 0:
            func_2_interpolate = patchIntegrand('fs')
            
        if integral_func == 0:
            integral_func = internal_integral
            
        try:
            if xs == 0 | ys == 0:
                self.startpoint = start
                self.endpoint = end
                lnx = np.linspace(np.log(start), np.log(end), num=200, endpoint=True)
                ys = ii_array(np.exp(lnx), func_2_interpolate, integral_func = integral_func)
            else:
                lnx = np.log(xs)
        except TypeError:
            lnx = np.log(xs)
            self.startpoint = xs.min()
            self.endpoint = xs.max()
        
        if min(ys) > 0:
            lny = np.log(ys)
            self.ispos = True
            self.loginterp = scipy.interpolate.interp1d(lnx, lny, kind = 'cubic')
        else:
            self.loginterp = scipy.interpolate.interp1d(lnx, ys, kind = 'cubic')
            self.ispos = False
            def basic(self, x):
                '''
                basic returns the part of the function independent of height
                '''
                if (x < self.endpoint )& (x > self.startpoint):
                    return self.loginterp(np.log(x))
                elif x < self.endpoint:
                    return 0
                else:
                    return self.loginterp(np.log(self.endpoint))    
   
    def basic(self, x):
        '''
        basic returns the part of the function independent of height
        '''
        if (x < self.endpoint )& (x > self.startpoint):
            if self.ispos:
                return np.exp(self.loginterp(np.log(x)))
            else:
                return self.loginterp(np.log(x))
        elif x < self.endpoint:
            return 0
        else:
            if self.ispos:
                return np.exp(self.loginterp(np.log(self.endpoint)))
            else: 
                return self.loginterp(np.log(self.endpoint))
                
    def transformed(self, h, r):
        '''
        returns the interpolated function once the height has been included
        '''
        if r>0:
            return self.basic(h/r)/h**2/(np.pi**2)
        else:
            return self.basic(self.endpoint)/h**2/(np.pi**2)
    
    def singleh(self, h):
        '''singleh returns a function to give the value of the interpolation
        at a single height, for all possible radii
        '''
        def rfunc( r ):
            return self.transformed(h, r)
        return rfunc

def filter_h_pixels(hs, nside, rinterp, dx):
    '''This function generates the patch potential filter as a function of 
    position. A dictionary which gives the filter as a function of height
    for each pixel is output. The inputs are:
    hs -- a list of all the heights
    nside -- the number of pixels in the filter
    rinterp -- the function interpolated by radius
    dx -- the spacing between pixels'''
    fullfilt = [[[0 for i in range(1,j+1)] for j in range(1,(nside+3)//2)] \
    for k in range(len(hs))]
    locs = [(i-1,j-1) for j in range(1,(nside+3)//2) for i in range(1,j+1)] 
    
    for i in range(len(hs)):
        rfunc = rinterp.singleh(hs[i])
        locfilt = define_filter(rfunc,nside, dx)
        print(i)
        
        for j in range((nside+1)//2):
            for k in range(j+1):
                fullfilt[i][j][k] = locfilt[(nside-1)//2+k,(nside-1)//2+j]

    pixels = {j:[fullfilt[i][j[1]][j[0]] for i in range(len(hs))] for j in locs}            
    pixels['separation'] = hs
    
    return pixels
    
def complete_filtering(hs, nside, rlabel, dx):
    '''This function generates a filter of a given pixel size and list of 
    heights for a given radial function. The inputs are:
    hs -- heights
    nside -- the number of pixels on each side of the filter
    rlabel -- the type of function to use (see patchIntegrand())
    dx -- the spacing between pixel centers'''
    rfunc = InterpRes(func_2_interpolate = patchIntegrand(rlabel))
    mainfilter = filter_h_pixels(hs, nside, rfunc, dx)
    return mainfilter
    
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
  
def initialize_filtered_images( h_vals, raw_image, filter_int):
    '''Converts a cleaned KPFM voltages image into the image charge voltage at a certain separations.
    
    h_vals is list of separations.
    raw_image should be a np. array.
    filter_int is a dictionary describing a filter.
    the function returns a dictionary of the voltage from the image charges as several heights, 
    as a dictionary with the heights as the keys.'''
    filtered_images = {}
    starttime = time.time()
    for i in h_vals:
        print(time.time()-starttime)
        loc_filt = f3.make_filter(i,filter_int)
        filtered_images[i] = signal.convolve2d(raw_image, loc_filt, boundary='symm', mode='same') 
        
    return filtered_images

def KPFMimagefilter( raw_image, loc_filter):
    '''Converts a cleaned KPFM voltages image into the image charge voltage at a certain separations.
    
    raw_image should be a np.array.
    filter_int is a np.array also.
    the function returns a dictionary of the voltage from the image charges as several heights, 
    as a dictionary with the heights as the keys.'''
    med = np.median(raw_image)
    total_lf = np.sum(loc_filter)
    resid = 1 - total_lf
    #starttime = time.time()
    #print(time.time()-starttime)
    filtered_image = signal.convolve2d(raw_image, loc_filter, boundary='fill', fillvalue = med, mode='same') 
    filtered_image = filtered_image + med*resid*np.ones_like(filtered_image)    
    return filtered_image
  
def allfilteredImages(raw_image, filters, h_values):
    allfiltered = {}
    for x, i in enumerate(h_values):
        allfiltered[i] = KPFMimagefilter(raw_image, filters[x])
        
    return allfiltered
  
def thePhiller( submatrix ):
    a = len(submatrix)-1
    nside = 2*a + 1
    print(a,nside)
    array2fill = np.zeros((nside, nside))
    array2fill[a:,a:] = submatrix[:,:]
    array2fill[0:a,0:a] = submatrix[a:0:-1,a:0:-1]
    array2fill[0:a,a:] = submatrix[a:0:-1,:]
    array2fill[a:,0:a] = submatrix[:,a:0:-1]
    
    return array2fill