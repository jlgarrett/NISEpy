# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 10:08:32 2018

@author: joe
"""
import numpy as np
import tables as tb
from scipy import constants
import patchpotentials as pt
import tqdm

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

def expand_from_dict(filt_dict):
    for i in filt_dict.keys():
        try:
            if i[0] > maxi:
                maxi = i[0]
        except TypeError:
            next
    
    edge = 2*maxi+1
    to_fill = np.zeros((edge, edge))
    
    for i in filt_dict.keys():
      #print(i,i[0],i[1])
        try:
            to_fill[maxi+i[0], maxi+i[1]] = filt_dict[i][0]
        except TypeError:
            next  
            
    to_fill = expand_filter(to_fill)
    return to_fill
    
  
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
    #print(hlog)
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
  
def make_grid( center, lside, nside):
    '''Returns a set of points as a list of tuples, centered 
    with side length lside and nside points'''
    dist_btw_pts = lside//nside
    while lside % nside > 0: #we include this just incase stuff doesn't quite line up
        lside = lside + 1
    minx, miny = tuple([center[i]-dist_btw_pts*(nside//2) for i in range(2)])
    points = [ (minx+i*dist_btw_pts, miny+j*dist_btw_pts) for i in range(nside) for j in range(nside)]
    return points
  
def V0_for_df( v1, v1_tilde, sphere_heights, dx, v2 = 0 , v2_hat = 0, R = 4e-5, ext_voltage = 0 ):
    dx = float(dx)
    
    voltage2min = v2 + v2_hat - (v1 + v1_tilde)
    v0_unnorm, v0_normalizer = minimizing_voltage_fd_image( voltage2min, sphere_heights, dx, inner = True)
    d = sphere_heights.min()
    
    
    v0 = (v0_unnorm-ext_voltage*(spherePFAdF(d,R) - v0_normalizer))/spherePFAdF(d, R)
    
    return v0

def V0_for_force( v1, v1_tilde, sphere_heights, dx, v2 = 0 , v2_hat = 0, R = 4e-5, ext_voltage = 0 ):
    '''This is like the minimizing_voltage_PFA function, except that it takes the actual maps rather
    than dictionaries'''
    dx = float(dx)
    
    voltage2min = v2 + v2_hat - (v1 + v1_tilde)
    v0_unnorm, v0_normalizer = minimizing_voltage_force_image( voltage2min, sphere_heights, dx, inner = True)
    d = sphere_heights.min()
    
    #should probably include the full PFA solution here
    v0 = (v0_unnorm-ext_voltage*(spherePFA(d,R)- v0_normalizer))/spherePFA(d,R)
    
    return v0
  
def spherePFA( d, R ):
    #I've set it up to use the full rather than the leading-order PFA in order to better incorporate the image
    #The full PFA actually is no better at predicting the exact solution than the partial PFA
    PFA = 2*np.pi*(R/d-np.log(R/d +1))
    return PFA
  
def spherePFAdF( d, R ):
    #I've set it up to use the full rather than the leading-order PFA in order to better incorporate the image
    #The full PFA actually is no better at predicting the exact 
    PFA = 2*np.pi*R/d**2*(1-1/(R/d +1))
    return PFA
  
def minimizing_voltage_fd_image( voltage_map, heights , dx,  inner = False):
    '''calculates the minimizing voltage'''
    
    v0_unnorm = 2*0.5*np.sum(voltage_map/heights**3)*dx**2
    v0_normalizer = 2*np.sum(1/heights**3)*dx**2
    #print(v0_unnorm,v0_normalizer)
    
    if inner:
        return v0_unnorm, v0_normalizer
    else:
        return v0_unnorm/v0_normalizer

def minimizing_voltage_force_image( voltage_map, heights , dx,  inner = False):
    '''calculates the minimizing voltage'''
    
    v0_unnorm = 0.5*np.sum(voltage_map/heights**2)*dx**2
    v0_normalizer = np.sum(1/heights**2)*dx**2
    #print(v0_unnorm,v0_normalizer)
    
    if inner:
        return v0_unnorm, v0_normalizer
    else:
        return v0_unnorm/v0_normalizer

def minimizing_voltage_PFA( v1, v1_tilde_d, heights, dx,v2=0 , v2_hat_d=0, R = 4e-5, ext_voltage = 0 ):
    if v2 == 0:
        v2 = np.zeros(v1.shape)
        v2_hat_d = {0:np.zeros(v1.shape)}
        
    v1tilde = patches_on_sphere2(heights, v1_tilde_d)
    v2_hat = patches_on_sphere2(heights, v2_hat_d)
    
    voltage2min= v2+v2_hat-(v1+v1tilde)
    v0_unnorm, v0_normalizer = minimizing_voltage_force_image( voltage2min, heights, dx, inner = True)
    d = heights.min()
    #print(d,v0_unnorm,v0_normalizer,R,ext_voltage)
    dx = float(dx)
    
    v0 = (v0_unnorm-ext_voltage*(2*np.pi*R/d - v0_normalizer))/(2*np.pi*R/d)
    #print(v0)
    
    return v0
  
def force_from_heights( sphere_image, dx, v1, v1_tilde, v1_hat, v2, v2_tilde, v2_hat,\
                      ext_voltage = 0, app_voltage = 0, radius = 4e-5):
    
    d = np.amin(sphere_image)
    area2sum = 1/sphere_image**2*dx**2
    force2sum = constants.epsilon_0/2*((v1+app_voltage)*(v1_tilde+app_voltage)- \
                                       (v1+app_voltage)*v2_hat - v2*(v1_hat+app_voltage) + v2*v2_tilde)*area2sum
    total_force = np.sum(force2sum) + \
        constants.epsilon_0/2*(ext_voltage+app_voltage)**2*(spherePFA( d, radius ) - np.sum(area2sum))
    
    return total_force
  
def force_expanded( sphere_image, dx, v1, v1_tilde, v1_hat, v2, v2_tilde, v2_hat,\
                      ext_voltage = 0, app_voltage = 0, radius = 4e-5):
    
    d = np.amin(sphere_image)
    area2sum = 1/sphere_image**2*dx**2
    force2sum = constants.epsilon_0/2*np.sum(area2sum)*np.array([np.sum((v1+app_voltage)*(v1_tilde+app_voltage)),
                                                np.sum(-(v1+app_voltage)*v2_hat),np.sum( - v2*(v1_hat+app_voltage)), np.sum(v2*v2_tilde)])
    residual_force = constants.epsilon_0/2*(ext_voltage+app_voltage)**2*(spherePFA( d, radius ) - np.sum(area2sum))
    
    return force2sum, residual_force
  
def df_from_heights(sphere_image, dx, v1, v1_tilde, v1_hat, v2, v2_tilde, v2_hat,\
                      ext_voltage = 0, app_voltage = 0, radius = 4e-5):
    
    d = np.amin(sphere_image)
    area2sum = 1/sphere_image**3*dx**2
    force2sum = constants.epsilon_0*((v1+app_voltage)*(v1_tilde+app_voltage)- \
                                       (v1+app_voltage)*v2_hat - v2*(v1_hat+app_voltage) + v2*v2_tilde)*area2sum
    total_force = np.sum(force2sum) + \
        constants.epsilon_0/2*(ext_voltage+app_voltage)**2*(spherePFAdF(d, radius) - np.sum(area2sum))
    
    return total_force
  
def quickforcesave(fileName, pixels, dataList, description = '' ):
    with tb.open_file(fileName, mode = 'w', title = description) as h5file:
        pix = np.array([[i[0],i[1]] for i in pixels])
        parray = h5file.create_array('/', 'center_pix', pix, 'the pixels on which sphere is centered')
        parray.attrs.plotting_note = 'remember x and y coordinates are switched when plotting on an image'
        
        dgroup = h5file.create_group(h5file.root, 'calcData', 'numerical data calculated from measurements at each separation')
        
        for i in dataList:
            tablestyle = {}
            for j in dataList[i]:
                if j == 'Title':
                    Title = dataList[i][j]
                    continue

                tablestyle[j] = tb.FloatCol()
                
            newtable = h5file.create_table(dgroup, i, tablestyle, title = Title )
                
            aRow = newtable.row            
            for m, j in enumerate(dataList[i]['separations']):
                for k in tablestyle:
                    #print(j)
                    aRow[k] = dataList[i][k][m]
                aRow.append()
            newtable.flush()
    return 0
  
def calcForceData( sphere, plate, grid, seps, justOne = False):
    all_data = {}
    l = len(seps)

    
    if justOne:
        grid2 = [grid[0]]
    else:
        grid2 = grid

    if not sphere.KPFM.shape == plate.KPFM.shape:
        print('The sphere and plate images (currently) need to be the same size')
        return -1
    if check(sphere, plate):
        all_forces = {}
        all_v0s = {}
        print('forces')
        for u,i in enumerate(grid2):
            for w, j in tqdm.tqdm(enumerate(grid)):
                labelstr  = 's'+str(u).zfill(2) + 'p' + str(w).zfill(2)
                
                (toptopo, bottopo, wholesphere), shift = pt.shiftByCenters(sphere.KPFM.shape, i, j, sphere.dx)
                vPkpfm = pt.shiftFill(shift, plate.KPFM, bottom = True)
                vSkpfm = pt.shiftFill(shift, sphere.KPFM)
        
                v0s = np.zeros(l)
                forces = np.zeros(l)
                for v, h in enumerate(seps):
                    vSs1 = pt.patches_on_sphere2(toptopo + h, sphere.fs)
                    vSo1 = pt.patches_on_sphere2(toptopo + h, sphere.fo)
                    vPs1 = pt.patches_on_sphere2(bottopo + h, plate.fs)
                    vPo1 = pt.patches_on_sphere2(bottopo + h, plate.fo)
                    vPs = pt.shiftFill( shift, vPs1, bottom = True)
                    vPo = pt.shiftFill( shift, vPo1, bottom = True)
                    vSs = pt.shiftFill( shift, vSs1)
                    vSo = pt.shiftFill( shift, vSo1)
                    v0 = V0_for_force(vSkpfm, vSs, wholesphere+h, dx=sphere.dx, v2 = vPkpfm, v2_hat = vPo, R = sphere.R)
                    force = force_from_heights( wholesphere+h, sphere.dx, vSkpfm, vSs, vSo, vPkpfm, vPs, vPo,
                      ext_voltage = 0, app_voltage = v0, radius = sphere.R)
                    v0s[v] = v0
                    forces[v] = force
            
                all_forces[labelstr] = forces
                all_v0s[labelstr] = v0s
        all_forces['Title'] = 'The forces calculated form the sphere ' + sphere.name
        all_v0s['Title'] = 'The force-minimizing voltages calculated form the sphere ' + sphere.name
        all_data['f'] = all_forces
        all_data['v0f'] = all_v0s
        
    if check(sphere, plate, f = 'fd'):
        all_df = {}
        all_v0df = {}
        print('force gradients')
        for u,i in enumerate(grid2):
            for w, j in tqdm.tqdm(enumerate(grid)):
                labelstr  = 's'+str(u).zfill(2) + 'p' + str(w).zfill(2)
                
                (toptopo, bottopo, wholesphere), shift = pt.shiftByCenters(sphere.KPFM.shape, i, j, sphere.dx)
                vPkpfm = pt.shiftFill(shift, plate.KPFM, bottom = True)
                vSkpfm = pt.shiftFill(shift, sphere.KPFM)
        
                v0s = np.zeros(l)
                forces = np.zeros(l)
                for v, h in enumerate(seps):
                    vSs1 = pt.patches_on_sphere2(toptopo + h, sphere.fds)
                    vSo1 = pt.patches_on_sphere2(toptopo + h, sphere.fdo)
                    vPs1 = pt.patches_on_sphere2(bottopo + h, plate.fds)
                    vPo1 = pt.patches_on_sphere2(bottopo + h, plate.fdo)
                    vPs = pt.shiftFill( shift, vPs1, bottom = True)
                    vPo = pt.shiftFill( shift, vPo1, bottom = True)
                    vSs = pt.shiftFill( shift, vSs1)
                    vSo = pt.shiftFill( shift, vSo1)
                    v0 = V0_for_df(vSkpfm, vSs, wholesphere+h, dx=sphere.dx, v2 = vPkpfm, v2_hat = vPo, R = sphere.R)
                    force = df_from_heights( wholesphere+h, sphere.dx, vSkpfm, vSs, vSo, vPkpfm, vPs, vPo,
                      ext_voltage = 0, app_voltage = v0, radius = sphere.R)
                    v0s[v] = v0
                    forces[v] = force
            
                all_df[labelstr] = forces
                all_v0df[labelstr] = v0s

        all_df['Title'] = 'The force derivatives calculated form the sphere ' + sphere.name
        all_v0df['Title'] = 'The force derivative-minimizing voltages calculated form the sphere ' + sphere.name
        all_data['df'] = all_df
        all_data['v0df'] = all_v0df        
    return all_data
    
def check(sphere, plate, f = 'f'):
    if f == 'f':
        return not ((sphere.fs == 0) | (sphere.fo == 0) | (plate.fo == 0) | (plate.fs == 0))
    else:
        return not ((sphere.fds == 0) | (sphere.fdo == 0) | (plate.fdo == 0) | (plate.fds == 0))
    
def quick_f_V0(sph, pl, center, shift):
    hmap = 1e-7
    vSs = pt.patches_on_sphere2(hmap, sph.fs)
    vPo = pt.patches_on_sphere2(hmap, pl.fo)
    v0 = V0_for_force(sph.KPFM, vSs, hmap, sph.dx, v2 = pl.KPFM, v2_hat = vPo, R = sph.R)
    return v0

def quick_force(sph, pl, hmap):
    vSs = pt.patches_on_sphere2(hmap, sph.fs)
    vSo = pt.patches_on_sphere2(hmap, sph.fo)
    vPs = pt.patches_on_sphere2(hmap, pl.fs)
    vPo = pt.patches_on_sphere2(hmap, pl.fo)
    force = force_from_heights( hmap, sph.dx, sph.KPFM, vSs, vSo, pl.KPFM, vPs, vPo, radius = sph.R)
    return force
  
class surface:
    def __init__( self, KPFM, fs = 0 , fo = 0, fds = 0, fdo = 0):
        self.KPFM = KPFM
        self.fs = fs
        self.fo = fo
        self.fds = fds
        self.fdo = fdo

class sphere(surface):
    def __init__(self, KPFM, fs = 0 , fo = 0, fds = 0, fdo = 0, R = 0, dx = 0, name = ''):
        super(sphere, self).__init__(KPFM, fs, fo, fds, fdo)
        self.R = R
        self.dx = dx
        self.name = name
        
def to_dictionaries( h5):
    with tb.open_file(h5, mode = 'r') as h5fil:
        separations = h5fil.root.separations.read()
        try:
            mo = h5fil.root.fo.read()
            ms = h5fil.root.fs.read()
        except Exception:
            mo = h5fil.root.fdo.read()
            ms = h5fil.root.fds.read()
        
        o = {}    
        s = {}    
        for i, x in enumerate(separations):
            o[x] = mo[:,:,i]
            s[x] = ms[:,:,i]
   
    return s, o  
  
def self_force_filt_interp_plates(v1, ds, filt_interpolated):
    force = np.zeros(len(ds))
    badforce = np.zeros(len(ds))
    #starttime = time.time()
    for i in range(len(ds)):
        #print(time.time()-starttime)
        filt_d = make_filter(ds[i],ffilt_int)
        force[i] = force_with_filter(v1, ds[i], filt_d)
        badforce[i] = ES_force_1term(v1, v1, ds[i])
        
    return force, badforce
