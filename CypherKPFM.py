# -*- coding: utf-8 -*-
"""
load_and_prep() loads an .ibw image file stored according to the traditions of 
Asylum research and stores it as a dictionary (because that is what the 'igor' 
package does by default)

to access metadata about the file, x['cleannote']
to list the data in the file, x['labels']
to extract a particular image, get_image(loadedibw, nameOfImage)

Created on Mon Jan 22 17:00:33 2018

@author: joe
"""

'''This first section of the code pertains to loading .ibw files and cleaning
the data therein for further processing'''

import igor.binarywave as igoribw #The package igor is not too hard to find 
import numpy as np
#import matplotlib  #  A powerful graphics package.
from matplotlib import cm
import matplotlib.pyplot as plt
from scipy import signal
from skimage import filters



def load_and_prep( ibw_location ):
    '''Loads an ibw file
    
    Given the location of an ibw file, this function loads the file and cleans
    the data as specified in 'cleanup_ibw_labels()' and 'cleanup_note().' The
    returned ibw is stored as a dictionary.'''
    
    #Code from the igor package loads the .ibw file
    loadedibw = igoribw.load(ibw_location)
    
    #the labels are rearranged to be more amenable to our analysis
    cleanup_ibw_labels( loadedibw )
    
    #The note, which stores all of the metadata, is rearranged to be more searchable
    cleanup_note( loadedibw )
    
    return loadedibw
    
def cleanup_ibw_labels( ibw ):
    '''This function translates the loaded IBW labels into a more readable format.
    
    The labels are stored in a dictionary called 'cleanlabels' in the form name:index
    The current version sorts through Trace and Retrace and should also separate out
    modified versions, but that feature has not been tested yet
    '''
    
    #The original labels are stored in an odd location
    originalLabels = ibw['wave']['labels'][2]
    newlabels = {}
    
    #Iterate through all the labels to copy them to the new location
    for i in range(len(originalLabels)):
        istr = ibw['wave']['labels'][2][i].decode("utf-8")
        
        #Basing the name on the location of 'trace' cuts of unwanted chars
        traceloc = (istr.lower()).find('trace')
        if traceloc > 0:
            scanName = istr[:(traceloc+5)]
            
            #We can also accommodate data that has been modified in Igor
            modloc = istr.find('Mod')
            if modloc > 0:
                scanName = scanName + istr[modloc:]
            
            #All the names are stored in newlabels
            newlabels[scanName] = i
    
    #the cleanlabels are stored in the original file 
    ibw['cleanlabels'] = newlabels        
    return newlabels
    
def cleanup_note( ibw ):
    '''Turns the note of the ibw file into a dictionary,
    which can be accessed as 'cleannote'.
    
    The function is based on the Asylum Research note notation'''
    
    
    originalnote = ibw['wave']['note'].decode('utf-8','ignore')


    newnote = {}
    dictlabel = ''
    dictvalue = ''
    tolabel = 1 #1 if characters are being added to label, 0 otherwise
    
    #Iterate through the note and create new dictionary entries whenever the
    #return character is found
    for i in originalnote:
        #print(i)
        if tolabel:
            if i ==':':
                tolabel = 0
            else:
                dictlabel= dictlabel + i
                
        else:
            if i == '\r':
                tolabel = 1
                newnote[dictlabel] = dictvalue
                dictlabel = ''
                dictvalue = ''
            else:
                dictvalue = dictvalue + i
    
    #add the size of a pixel to the note (this will be useful later!)
    try:
        newnote['dx'] = str(float(newnote['SlowScanSize'])/int(newnote['ScanLines']))
    except KeyError:
        newnote['dx'] =  None
        
    #add the note to the Scanning probe data
    ibw['cleannote'] = newnote
    
    return newnote
    
def check_labels( ibw ):
    '''Checks whether or not the labels are updated'''
    if 'cleanlabels' in ibw:
        return 1
    else:
        print('Error! You must update the labels first')
        return 0
    
def get_image( ibw , name ):
    '''Extracts one image from the .ibw file'''
    
    if check_labels( ibw ):
        #the -1 comes from an oddity of how wData is stored in the ibw file
        index = ibw['cleanlabels'][name] - 1  
        return ibw['wave']['wData'][:,:,index]
    
    else:
        print('Error! No such image')
        return None
    

def add_image( ibw, image ,name):
    '''Adds an image to an uploaded ibw
    
    ibw must be a loaded .ibw AR (as dictionary) file with cleanlabels
    image should be an ndarray, name should be a string'''
    
    #First check to see if the an image with the name has already been loaded 
    #into the .ibw, so that nothing is overwritten            
    if check_labels( ibw ):
        if name in ibw['cleanlabels']:
            print('image already loaded')
            return None
        else:
            #Ensure that the new data has the right shape!
            ibw_dims = (int(ibw['cleannote']['Initial ScanPoints']),
                int(ibw['cleannote']['Initial PointsLines']))
            if image.shape == ibw_dims:
                #to be consistent, it is necessary to add 1 to the index                
                indextoadd = max(ibw['cleanlabels'].values()) + 1
                ibw['cleanlabels'][name] = indextoadd
                #We concatenate the data to add new layers
                newtotal = np.concatenate((ibw['wave']['wData'],
                                           image[...,None]), axis=2)
                ibw['wave']['wData'] = newtotal
                return 1
            else:
                print('wrong size')
                return None
            
def rename_image(ibw, oldname, newname):
    '''Changes the label on an image'''
    if oldname in ibw['cleanlabels']:
        ibw['cleanlabels'][newname] = ibw['cleanlabels'].pop(oldname)
       
        

def avg_kpfms( ibw ):
    '''This function averages the KPFM trace and retrace data together.
    
    It should only be used if the gain was sufficiently high during the 
    collection of the data. In the future, I hope to replace it with more
    comprehensive preprocessing. But until then...'''
    
    #Our H-KPFM code stores the KPFM data in UserIn1
    #Check to ensure that both the trace and retrace are present in the data
    tracepresent = 'UserIn1Trace' in ibw['cleanlabels']
    retracepresent = 'UserIn1Retrace' in ibw['cleanlabels']
    if tracepresent & retracepresent:
        KPFM = 0.5*get_image( ibw, 'UserIn1Trace')+\
            0.5*get_image( ibw, 'UserIn1Retrace')
        add_image( ibw, KPFM ,'Averaged_KPFM')
        return KPFM
    else:
        return None
    
def clean_kpfms( ibw, setpoint = 0, tolerance = 0, consistencyTolerance = 0):
    '''This function combines the KPFM trace and retrace data together 
    in a way that considers the noise in the system.
    
    It is necessary to have UserIn2 data, as well as UserIn1 data'''
    
    print(ibw['cleanlabels'])    
    
    #First, check to see if we have all the necessary data
    avg_present = 'Averaged_KPFM' in ibw['cleanlabels']
    traceerrpresent = 'UserIn2Trace' in ibw['cleanlabels']
    retraceerrpresent = 'UserIn2Retrace' in ibw['cleanlabels']
    
    insufficientUserIn2 = ~(traceerrpresent & retraceerrpresent)    
    
    if insufficientUserIn2==True:
        print("insufficientUserIn2",traceerrpresent,
              retraceerrpresent,insufficientUserIn2)
        return None
    
    if avg_present:
        KPFM_avg = get_image(ibw, 'Averaged_KPFM')
        KPFM_best = np.empty_like(KPFM_avg)
        KPFM_manynans = np.empty_like(KPFM_avg)
        KPFM_best[:] = KPFM_avg
        KPFM_manynans[:] = KPFM_avg
    
    
    
    #We load all the data that are to be processed
    #Both the KPFM signal
    
    ui1r = get_image(ibw, 'UserIn1Retrace')   
    ui1t = get_image(ibw, 'UserIn1Trace')
    
    #And the error signal
    ui2t = get_image(ibw, 'UserIn2Trace')
    ui2r = get_image(ibw, 'UserIn2Retrace')
    
    #We setup the defaults of the sepoint and tolerance
    if setpoint == 0:
        setpoint = 0.5*(np.mean(ui2t)+np.mean(ui2r))
    
    if tolerance == 0:
        tolerance = .7e-3
    
    if consistencyTolerance == 0:
        consistencyTolerance = .1
    
    #Check and find which points suffice with a simple average and which 
    #require more sophisticated treatment
    consistent_bound = abs(ui1t-ui1r) < consistencyTolerance
    ui2t_inbound = (abs(ui2t - setpoint) < tolerance) & consistent_bound
    ui2r_inbound = (abs(ui2r - setpoint) < tolerance) & consistent_bound
    bothbound = (ui2t_inbound) & (ui2r_inbound)

    #These help to keep tabs on the opperation
    nans = []
    u2ts = []
    u2rs = []
    Lx, Ly=KPFM_best.shape
    count = 0
    countt = 0
    countr = 0
    
    #The program gives the best KPFM estimate at each point
    for i in range(Lx):
        for j in range(Ly):
            if bothbound[i][j]:
                count = count +1
                pass
            elif ui2t_inbound[i][j] == 1:
                KPFM_best[i][j] = ui1t[i][j]
                KPFM_manynans[i][j] = np.nan
                countt = countt + 1
                u2ts.append((ui1t[i][j],ui1r[i][j]))
            elif ui2r_inbound[i][j] == 1:
                KPFM_best[i][j] = ui1r[i][j]
                KPFM_manynans[i][j] = np.nan
                countr = countr + 1
                u2rs.append((ui1r[i][j],ui1t[i][j]))
            else: 
                KPFM_best[i][j] = np.nan
                KPFM_manynans[i][j] = np.nan
                nans.append((i,j))#adding 1 helps in the next section
    
    #use median filter (neglecting nans) to fill in the points where the KPFM
    #data we bad in both the trace and retrace channels    
    #KPFM_bounded = np.pad(KPFM_out,1, mode= 'median')
    nanout = []
    print(('nans length is {}, but {} were '
           'passed, {} in trace, {} in retrace').format(len(nans),
                    count,countt, countr))
    count = 0
    while (np.sum(np.isnan(KPFM_best))) > 0 & (count < 10):
        for i in nans:
            box_min_x = max(i[0]-3,0)
            box_min_y = max(i[1]-3,0)
            box_max_x = min(i[0]+4, Lx)
            box_max_y = min(i[1]+4, Ly)
            box = KPFM_manynans[box_min_x:box_max_x,box_min_y:box_max_y]
            #print(box)
            box = box.flatten('F')
            box = box[np.nonzero(box)]
            #print(box)
            nanval = float(np.nanmedian(box))
            #print(nanval)
            KPFM_best[i[0],i[1]] = nanval
            nanout.append(nanval)
        if np.sum(np.isnan(KPFM_best)) > 0:
            KPFM_manynans = KPFM_best
            KPFM_best = np.empty_like(KPFM_manynans)
            KPFM_best = KPFM_manynans[:][:]
        count = count + 1
        
    #print(nanout)
    
    #one more median filter cleans the image up a little bit more    
    KPFM_out = signal.medfilt2d(flatten(KPFM_best),5)
    
    #flattening helps eliminate line noise
    KPFM_out = flatten(KPFM_out)
    
    #Add the image to the ibw
    add_image( ibw, KPFM_best ,'Cleaned_KPFM')
    return KPFM_best
    
    
def flatten( image, axis = 0, mask = '0', order = 0):
    '''Flattens the individual lines of an image.

    Right now it supports only 0th order flattening, which is most common
    for KPFM image analysis. The option 'axis' permits one to change the 
    direction of the flattening'''
    #transpose for other axis    
    if axis == 1:
        image = image.transpose()
    
    
    flattened_image = np.array(image)
    flattening = np.array(image)
    
    try:
        flattening[mask] = np.nan
    except IndexError:
        print('no mask')
        
    #we create another value to act on    

    for i in range(image.shape[0]):
        pfit = np.polyfit(np.indices(flattened_image[:,i].shape).flatten(),flattened_image[:,i], order)
        p = np.poly1d(pfit)
        flattened_image[:,i]=flattened_image[:,i]-p(np.indices(flattened_image[:,i].shape))
        
    #transpose back    
    if axis == 1:
        flattened_image = flattened_image.transpose() 
        
    return flattened_image

def quickplot( kpfm, clim = 0, boundary = 0):
    '''Quickly plot a ndarray. Default settings are for KPFM image'''
    if ~clim:
        clim = (np.min(kpfm), np.max(kpfm))
    fig = plt.imshow(np.transpose(kpfm), interpolation='nearest',  clim=clim,
                     cmap=cm.get_cmap('Greens_r', 50))
    #if boundary > 0:
    #    fig = ax.plot(boundary)
    plt.colorbar(label = 'Potential (V)')
    plt.xlabel('pixels')
    plt.ylabel('pixels')
    plt.show()
    
def basicTopography( ibw ):
    '''Extracts topography and 0th order flattening'''
    #want different cases for trace, retrace
    tracepresent = 'ZSensorTrace' in ibw['cleanlabels']
    retracepresent = 'ZSensorRetrace' in ibw['cleanlabels']
    #print(tracepresent, retracepresent)
    if tracepresent:
        hT = get_image( ibw, 'ZSensorTrace')
        hT = flatten(hT)
        hTmask = autoMask(hT)
        hT = flatten(hT, mask = hTmask)
        add_image( ibw, hT ,'ZST_flattened0')
    if retracepresent:
        hR = get_image( ibw, 'ZSensorRetrace')
        hR = flatten(hR)
        hRmask = autoMask(hR)
        hR = flatten(hR, mask = hRmask)
        add_image( ibw, hR ,'ZSR_flattened0')
        return hR
    else:
        return None
    
def autoMask( image , toMask = 'high'):
    '''Still under construction'''
    mask = np.array(image)
    val = filters.threshold_otsu(image)
    if 'h' in toMask:
        return mask > val
    else:
        return mask < val

    return 0
        
def normalize(image):
    return (image-np.min(image))/(np.max(image)-np.min(image))

def quickSave(image, name, prefix = ''):
    np.savetxt(prefix+name+'.csv', get_image(image,name), delimiter = ',')
    return 0
    
