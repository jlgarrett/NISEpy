# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 17:03:25 2018

@author: joe
"""
import matplotlib.pyplot as plt
import numpy as np
import pickle
from scipy import signal
from scipy.spatial.distance import euclidean
from scipy.spatial import distance
from scipy import constants
import scipy

'''The next several functions are designed for calculating the force from a 
single KPFM image in the parallel-plate formulation'''


def autocorrelate_linesfast(data):
    '''Gives the autocorrelation of a 2D array of numerical as a function of 
    the number of pixels based on the auto correlation of individual vertical 
    and horizontal lines.
    
    The input data is anticipated as a 2D numpy array. This function is faster
    than the full autocorrelation of autocorrelate_2d_norm()'''
    
    #processing a null signal relieves us of the necessity of combinitorics
    null_data = np.ones((data.shape[0],data.shape[1])) 
              
    #We calculate the autocorrelation function as a line    
    autoc_line = np.zeros((2*data.shape[0]-1))                    
    null_autoc = np.zeros(2*data.shape[1]-1)
    
    #We calculate the autocorrelation function in both the x and y directions
    for x in range(data.shape[0]): 
        autoc_line += signal.fftconvolve(data[x,:], data[x,::-1], mode='full')
        null_autoc += signal.fftconvolve(null_data[x,:], null_data[x,::-1], mode='full')
                    
    for x in range(data.shape[1]):
        autoc_line += signal.fftconvolve(data[:,x], data[::-1,x], mode='full')
        null_autoc += signal.fftconvolve(null_data[:,x], null_data[::-1,x], mode='full')
      
    #This prevents dividing by zero  
    for x in np.nditer(null_autoc, op_flags=['readwrite']):
        x[...] = max(x,1)
        
    #The autocorrelation is normalized by the size of the image    
    autoc_line /= null_autoc
    
    autocorrelation = autoc_line[(data.shape[0]-1):]
    
    return autocorrelation

def autocorrelate_2d_norm(data):
    '''This function computes the 2D autocorrelation function from 2D data as 
    an image, and normalizes it'''
    null_data = np.ones(data.shape)             
    
    
    autoc_2d = signal.fftconvolve(data, data[::-1,::-1], mode='full')
    null_autoc_2d = signal.fftconvolve(null_data, null_data[::-1,::-1], mode='full')
    
    for x in np.nditer(null_autoc_2d, op_flags=['readwrite']):
        x[...] = max(x,1)
    
    autoc_2d /= null_autoc_2d

    return autoc_2d

def autocorrelate_2Dac_toline(autocor_2d):
    '''This function adds up all the points in a 2D autocorrelation function to
    make a 1D autocorrelation function'''
    length = int(np.floor(np.sqrt(((autocor_2d.shape[0]+1)/2)**2 \
        + ((autocor_2d.shape[1]+1)/2)**2)))
    autocor = np.zeros(length)
    null_autocor = np.zeros(length)
    center = np.array([(autocor_2d.shape[0]+1)/2-1,
                       (autocor_2d.shape[1]+1)/2-1])
    
    for x in range(autocor_2d.shape[0]):
        for y in range(autocor_2d.shape[1]):
            #we need to find which box we're looking in!
            box = int(np.floor(euclidean(np.array([x,y]),center))) 
            autocor[int(box)] += autocor_2d[x,y]
            null_autocor[int(box)] += 1
    
    for x in np.nditer(null_autocor, op_flags=['readwrite']):
        x[...] = max(x,1)
        
    autocor /= null_autocor
    
    return autocor#autocor

def coord_mat(dimx, dimy):
    '''Gives a matrix listing the coordinates'''
    xs = np.arange(dimx)
    ys = np.arange(dimy)
    coords = np.empty((dimx, dimy,2), dtype=np.intp)
    coords[..., 0] = xs[:, None]
    coords[..., 1] = ys
    return coords

def autocorrelate_2Dac_toline_v2(autocor_2d):
    '''This function adds up all the points in a 2D autocorrelation function to
    make a 1D autocorrelation function'''
    center = np.array([(autocor_2d.shape[0]+1)/2-1,(autocor_2d.shape[1]+1)/2-1])
    coords = coord_mat(autocor_2d.shape[0],autocor_2d.shape[1])
    dist = distance.cdist(coords.reshape((autocor_2d.shape[0]*autocor_2d.shape[1],2)),
                          [center], 'euclidean')
    dist = dist.reshape((autocor_2d.shape[0],autocor_2d.shape[1]))
    
    length = int(np.floor(dist.max()-0.5))
    pointcounter = np.ones(autocor_2d.shape)
    autocor = np.zeros(length)
    null_autocor = np.zeros(length)
    
    for i in range(length): 
        acsum  = sum(autocor_2d[abs(dist-i)<(0.5)])
        autocor[i] = acsum
        nullsum = sum(pointcounter[abs(dist-i)<0.5])
        null_autocor[i] = nullsum
    
    autocor /= null_autocor
    return autocor

def autocorrelate_2D_toline(data):
    '''Takes 2D data and outputs the radial autocorrelation function for that data'''
    autoC_2d = autocorrelate_2d_norm(data)
    autoC_line = autocorrelate_2Dac_toline(autoC_2d)
    return autoC_line

def autocorrelate_2D_toline_v2(data):
    '''Takes 2D data and outputs the radial autocorrelation function for that data'''
    autoC_2d = autocorrelate_2d_norm(data)
    autoC_line = autocorrelate_2Dac_toline_v2(autoC_2d)
    return autoC_line

class MultipleAutocorrelations:
    '''this class is a way to store our calculated autocorrelation functions and plot them easily'''
    
    def __init__(self, num):
        self.num = num
        self.ACdata = [{'data' : [],
                      'description' : '',
                      'scale' : 1 } for j in range(num)]
        self.current = 0
        self.MainDescription = ''
            
    def __iter__(self):
        self.current = 0
        return self
        
    def __next__(self):
        trial = self.current + 1
        self.current = trial% (self.num)
        return self.current    

    def add_2_MA(self, lista):
        self.ACdata[self.current]['data'] = lista
        return (self.current)
    
    def change_scale(self, scale):
        self.ACdata[self.current]['scale'] = scale
        return (self.current)
    
    def add_2_descriptor(self, descriptor_string):
        self.ACdata[self.current]['description'] = descriptor_string
        return (self.current)
    
    def plot_autocorrs(self, mylist ):
        #self.xaxes = [len(self.ACdata[i]) for i in range(self.num)]
        #fig, arr = plt.subplots(1,len(mylist))
        for i in range(len(mylist)):
            plt.plot(self.ACdata[mylist[i]]['data'])
        
        plt.show()
        return 0
    
    def ploglog(self, mylist):
        for i in range(len(mylist)):
            xs = [ j*self.ACdata[mylist[i]]['scale'] for j in range(len(self.ACdata[mylist[i]]['data']))]
            plt.loglog(xs,abs(self.ACdata[mylist[i]]['data']))
        
        plt.show()
        return 0
    
    def add_main_description(self, description):
        self.MainDescription = description
        return 0
    
    def save(self, name):
        combineddict = { 'Main Description' : self.MainDescription,
                       'ACdata' : self.ACdata,
                       'ACversion' : 1}
        with open( name, 'wb') as f:
            pickle.dump(combineddict,f)
        return 0

class OneSizeAutocorrelations(MultipleAutocorrelations):
    '''this class is is like MultipleAutocorrelations, but forces all the autocorrelation functions to be the same size, and 
    so allow for computation of the average value, and the standard deviation '''
    
    def __init__(self, num):
        self.num = num
        self.ACdata = [{'data' : [],
                      'description' : '',
                      'scale' : 1 } for j in range(num)]
        self.current = 0
        self.setsize = 0
    
    def add_2_MA(self, lista):  
        #print(np.size(np.array(self.ACdata)))
        if(np.size(np.array([self.ACdata[i]['data'] for i in range(len(self.ACdata))])) == 0):
            self.setsize = len(lista)
            #self.ACdata=np.zeros(self.num,self.setsize)
            self.ACdata[self.current]['data'] = lista
            added = 1
        else:
            if(len(lista)==self.setsize):
                self.ACdata[self.current]['data'] = lista
                added = 1
            else:
                print("wrong size!")
                added = 0
        return (self.current, added)
    
    
    def is_not_empty(self):
        return np.array([len(self.ACdata[i]['data'])>0 for i in range(self.num)])
    
    def averages(self):
        to_average = np.array([self.ACdata[i]['data'] for i in range(self.num)])[self.is_not_empty()]
        #print(to_average)
        avg = np.array([sum([to_average[i][j] for i in range(len(to_average))]) for j in range(self.setsize)])
        avg =avg/len(to_average)
        return avg
    
    def sdevs(self):
        averages = self.averages()
        to_sdev = np.array([self.ACdata[i]['data'] for i in range(self.num)])[self.is_not_empty()]
        #print(averages)
        sdev = np.array([sum([(to_sdev[i][j]-averages[j])**2 for i in range(len(to_sdev))]) for j in range(self.setsize)])
        sdev = sdev/len(to_sdev)
        sdev = np.sqrt(sdev)
        return sdev
    
def loadAC( filename ):
    '''Loads autocorrelation functions stored as described above'''
    with open( filename, 'rb') as f:
        combineddict = pickle.load(f)
    
    if(len(combineddict)!=3):
        return 1
    
    if(rectangular([combineddict['ACdata'][i]['data'] for i in range(
                                            len(combineddict['ACdata']))])):
        loadedAC = OneSizeAutocorrelations(len(combineddict['ACdata']))
        loadedAC.setsize = len(combineddict['ACdata'][0]['data'])
    else:
        loadedAC = MultipleAutocorrelations(len(combineddict['ACdata']))
        
    loadedAC.ACdata = combineddict['ACdata']
    loadedAC.MainDescription = combineddict['Main Description']
    
    return loadedAC


def rectangular(n):
    '''Checks to see if data are rectangular'''
    lengths = {len(i) for i in n}
    return len(lengths) == 1

def Ck_integrand( x , dor):
    '''integrand of the broadening calculation function
    
    Here, x has been normalized by the radius'''
    
    return dor**4 * x**3 / (np.sinh(x*dor))**2 * scipy.special.j0(x)

def convert_Cr_to_Ck( Crs, dx , L=0):
    '''Converts autocorrelation function in r-space to k-space to convert 
    to a force'''
    if L == 0:
        L = Crs.shape[0]//3
    rs = [dx*i for i in range(L)]
    k2sum = np.zeros(L)
    ks = np.array([ 10**(i*0.001) for i in range(int(1e4))])
    Cks = np.zeros(ks.shape)
    
    for j in range(len(ks)):
        k2sum = np.array([scipy.special.j0(rs[i]*ks[j])*i*dx*Crs[i] for i in range(L)])
        Cks[j] = 2*np.pi*np.trapz(k2sum,rs)
        
    return Cks, ks

def pressure_from_C11(C11, ks, ds):
    '''Converts k-space autocorrelation data into a pressure for separations
    ds, using k values ks'''
    np.seterr(all = 'ignore')
    e0 = constants.epsilon_0
    klen = len(ks)
    dlen = len(ds)
    pressures = np.zeros(dlen)
    
    for i in range(dlen):
        p2sum = [C11[j]*ks[j]**3/np.sinh(ks[j]*ds[i])**2 for j in range(klen)]
        pressures[i] = e0/(4*np.pi)*np.trapz(p2sum, ks)
        
    return pressures
  
def dpressure_from_C11(C11, ks, ds):
    '''Converts k-space autocorrelation data into a pressure for separations
    ds, using k values ks'''
    np.seterr(all = 'ignore')
    e0 = constants.epsilon_0
    klen = len(ks)
    dlen = len(ds)
    pressures = np.zeros(dlen)
    
    for i in range(dlen):
        p2sum = [C11[j]*ks[j]**4/np.sinh(ks[j]*ds[i])**3*np.cosh(ks[i]*ds[i]) for j in range(klen)]
        pressures[i] = e0/(2*np.pi)*np.trapz(p2sum, ks)
        
    return pressures