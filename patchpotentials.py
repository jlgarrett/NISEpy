# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 17:20:31 2017

@author: JLG

This .py file is dedicated to the calculation of patch potential forces from KPFM images.

The first section concerns estimating the force from a single image, while the  
later part of the file is dedicated to calculation of patches from curved
surfaces, in which the patch potentials are resolved on both surfaces
"""

#It is important to note that the igor package is required for this analysis
import numpy as np
import pickle
from scipy import signal
from scipy.spatial.distance import euclidean
from scipy.spatial import distance
from scipy import constants
from scipy import special
import scipy
import time
from scipy.integrate import dblquad



#We are going to base all of our image analysis on the original .ibw, at least initially,
#so that we do not lose access to any of the data contained therin
#these first few functions are written to clear up the 








def ES_force(v1,v2,d):
    #constants.epsilon_0
    try:
        len(d)
    except TypeError:
        d = np.zeros(v1.shape)+d
        
    try:
        len(v2)
    except TypeError:
        v2 = np.zeros(v1.shape)+v2
    
    force2sum = constants.epsilon_0*(v1-v2)**2/d**2/2
    
    return sum(sum(force2sum))

def ES_force_1term(v1,v2,d):
    #constants.epsilon_0
    try:
        len(d)
    except TypeError:
        d = np.zeros(v1.shape)+d
        
    try:
        len(v2)
    except TypeError:
        v2 = np.zeros(v1.shape)+v2
    
    force2sum = constants.epsilon_0*v1*v2/d**2/2
    
    return np.sum(force2sum)

def force_with_filter(v1,d,filt):
    filtered = signal.convolve2d(v1, filt, boundary='symm', mode='same')    
    return ES_force_1term(v1,filtered,d)

def self_force_filt_interp_plates(v1, ds, filt_interpolated):
    force = np.zeros(len(ds))
    badforce = np.zeros(len(ds))
    starttime = time.time()
    for i in range(len(ds)):
        print(time.time()-starttime)
        filt_d = make_filter(ds[i],filt_interpolated)
        force[i] = force_with_filter(v1, ds[i], filt_d)
        badforce[i] = ES_force_1term(v1, v1, ds[i])
        
    return force, badforce

def initialize_filtered_images( h_vals, raw_image, filter_int):
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
        filtered_images[i] = signal.convolve2d(raw_image, loc_filt, boundary='symm', mode='same') 
        
    return filtered_images

def minimizing_voltage_force_image( voltage_map, heights , inner = False):
    '''calculates the minimizing voltage'''
    
    v0_unnorm = -0.5*np.sum(voltage_map/heights**2)
    v0_normalizer = np.sum(1/heights**2)
    
    if inner:
        return v0_unnorm, v0_normalizer
    else:
        return v0_unnorm/v0_normalizer

def minimizing_voltage_Force_PFA( voltage_map, heights, dx, R = 4e-5, ext_voltage = 0 ):
    v0_unnorm, v0_normalizer = minimizing_voltage_force_image( voltage_map, heights, inner = True)
    d = heights.min()
    dx = float(dx)
    
    
    return v0_unnorm/(2*np.pi*R/d)
    
def centered_spheremap(radius, n_pixels, closest, dx):
    center = (n_pixels//2, n_pixels//2)
    sphere = np.array([[np.sqrt(radius**2 - (dx*(i-center[1]))**2 - (dx*(j-center[0]))**2) for i in range(n_pixels)] for 
                                j in range(n_pixels)])
    sphere = - sphere + radius + closest
    return sphere

def patches_on_sphere(sphere, init_patches, num_heights, filter_interpolation):
    maxh, minh = np.amax(sphere), np.amin(sphere)
    h_values = np.logspace(np.log10(minh),np.log10(maxh),num_heights+1)
    patches_out = np.zeros(init_patches.shape)
    
    for i in range(len(h_values)-1):
        hloc = h_values[i]
        print(hloc)
        filt = make_filter(h_values[i],ffilt_int)
        convolved = signal.convolve2d(init_patches, filt, boundary='symm', mode='same')
        patches_out = np.where(sphere > hloc,convolved,patches_out)
        
    return patches_out

def patches_on_sphere2(sphere, filtered_patch_dict):
    heights = list(filtered_patch_dict.keys())
    heights.sort()
    patches_out = np.zeros(filtered_patch_dict[heights[0]].shape)
    
    for i in range(len(heights)-1):
        hloc = heights[i]
        patches_out = np.where(sphere > hloc,filtered_patch_dict[hloc],patches_out)
    
    return patches_out
        
def spheremap(radius, n_pixels, dx, center = 0, closest = 0):
    '''This function creates a map of the heights of a sphere above a surface'''
    
    
    
    if center == 0:
        center = (n_pixels//2, n_pixels//2)
    # print(center)
    dx = float(dx)
    try:
        sphere = np.array([[np.sqrt(radius**2 - (dx*(i-center[0]))**2 - (dx*(j-center[1]))**2) for j in range(n_pixels)] for i in range(n_pixels)])
    except TypeError:
        sphere = np.array([[np.sqrt(radius**2 - (dx*(i-center[0]))**2 - (dx*(j-center[1]))**2) for j in range(n_pixels[1])] for i in range(n_pixels[0])])
      
    sphere = - sphere + radius + closest
    return sphere

def shifted_sphere(orig_pix, shift, orig_center,  dx, radius = 4e-5):
    try:
        newx,newy = orig_pix[0]+abs(shift[0]), orig_pix[1]+abs(shift[1])
    except TypeError:
        newx,newy = orig_pix+abs(shift[0]), orig_pix+abs(shift[1])
    
    new_center = list(orig_center)
    for i in range(2):    
        if shift[i] < 0:
            #print(i)
            new_center[i] = new_center[i] - shift[i]
            
    #print(new_center, newx, newy)
            
    #print('radius is type {}'.format(type(radius)))
    #print('newx is type {}'.format(type(newx)))
    #print('dx is type {}'.format(type(dx)))
    #print('center is type {}'.format(type(new_center)))
    sphere = spheremap(radius, (newx,newy), dx=dx, center = new_center)
    
    return sphere

def subselect(orig_shape, shift, image, bottom = False):
    if bottom:
        shift = (-shift[0],-shift[1])
    
    start_x, start_y = 0,0
    if shift[0] > 0:
        start_x = shift[0]
    if shift[1] > 0:
        start_y = shift[1]
    
    try:
        return image[start_x:start_x+orig_shape[0],start_y:start_y+orig_shape[1]]
    except TypeError:
        return image[start_x:start_x+orig_shape,start_y:start_y+orig_shape]
      
def shiftFill(shift, image, bottom = False, fill = -1):
    shiftIm = np.zeros((image.shape[0]+abs(shift[0]),image.shape[1]+abs(shift[1])))
    if fill == -1:
        shiftIm[:,:] = np.median(image)
    else:
        shiftIm[:,:] = fill
    
    if bottom:
        shift = (-shift[0],-shift[1])
    
    start_x, start_y = 0,0
    if shift[0] > 0:
        start_x = shift[0]
    if shift[1] > 0:
        start_y = shift[1]
        
    shiftIm[start_x:start_x+image.shape[0], start_y:start_y + image.shape[1]] = image[:,:]
    
    return shiftIm
    
def invert( center, shift):
    orig = shifted_sphere(512, (0,0), center,  dx =  1e-5/512, radius = 4e-5)
    shifted = shifted_sphere(512, shift, center,  dx =  1e-5/512, radius = 4e-5)
    topTopo = subselect(512, shift, shifted, bottom = True)
    return orig, topTopo, shifted
  
def shift4heights( orig_shape, orig_center , shift, dx , radius = 4e-5):
    '''orig_center - where the bottom sphere is centered
    shift - how the top sphere is shifted relative to the bottom sphere'''
    shifted = shifted_sphere(orig_shape, shift, orig_center,  dx, radius = radius)
    topTopo = subselect(orig_shape, shift, shifted)
    botTopo = subselect(orig_shape, (-shift[0],-shift[1]), shifted)
    return topTopo, botTopo, shifted

def shiftByCenters( orig_shape, top_center, bot_center, dx, radius = 4e-5):
    shift = (bot_center[0]-top_center[0], bot_center[1]-top_center[1])
    return shift4heights(orig_shape, bot_center, shift, dx, radius = radius), shift