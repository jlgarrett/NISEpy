import CypherKPFM as ck
import numpy as np
import filters3D as filt3
import force3D as f3
import patchpotentials as pp
import patches2D as p2d
import pickle

print('loading files....')
n8B_plate = ck.load_and_prep('n8B_kpfm_surfloc0002.ibw')
n8B_sphere = ck.load_and_prep('n8B_kpfm0013.ibw')

print('cleaning KPFM....')
plate_avg = ck.avg_kpfms(n8B_plate)
plate_clean = ck.clean_kpfms(n8B_plate)
avg = ck.avg_kpfms(n8B_sphere)
sphere_ckpfm = ck.clean_kpfms(n8B_sphere)

print('calculating filters...')
h_values = np.logspace(-7.7,-5,21)

print('first filter')
try:
  gii_packed = pickle.load(open("gii_packed.pkl", "rb"))
except FileNotFoundError:
  gii_packed = filt3.complete_filtering(h_values,243, 'fs', float(n8B_plate['cleannote']['dx']))
  pickle.dump(gii_packed,open('gii_packed.pkl','wb') )

print('second filter')
try:
  gpii_packed = pickle.load(open("gpii_packed.pkl", "rb"))
except FileNotFoundError:
  gpii_packed = filt3.complete_filtering(h_values,243, 'fds', float(n8B_plate['cleannote']['dx']))
  pickle.dump(gpii_packed,open('gii_packed.pkl','wb') )

print('third filter')
try:
  gij_packed = pickle.load(open("gij_packed.pkl", "rb"))
except FileNotFoundError:
  gij_packed = filt3.complete_filtering(h_values,243, 'fo', float(n8B_plate['cleannote']['dx']))
  pickle.dump(gij_packed,open('gij_packed.pkl','wb') )

print('fourth filter')
try:
  gpij_packed = pickle.load(open("gpij_packed.pkl", "rb"))
except FileNotFoundError:
  gpij_packed = filt3.complete_filtering(h_values,243, 'fdo', float(n8B_plate['cleannote']['dx']))
  pickle.dump(gpii_packed,open('gij_packed.pkl','wb') )

print('processing plate...')
try:
  all_filtered = pickle.load(open('plate_filtered.pkl', 'rb'))
except FileNotFoundError:
  all_filtered = filt3.allfilteredImages(plate_clean, gii_packed)
  pickle.dump(all_filtered,open('plate_filtered.pkl','wb'))

try:
  fds_filtered = pickle.load(open('plate_filtered_fds.pkl', 'rb'))
except FileNotFoundError:
  fds_filtered = filt3.allfilteredImages(plate_clean, gpii_packed)
  pickle.dump(fds_filtered,open('plate_filtered_fds.pkl','wb') )

try:
  plate_fo = pickle.load(open('plate_fo.pkl', 'rb'))
except FileNotFoundError:
  plate_fo = filt3.allfilteredImages(plate_clean, gij_packed)
  pickle.dump(plate_fo,open('plate_fo.pkl','wb'))

try:
  plate_fdo = pickle.load(open('plate_fdo.pkl', 'rb'))
except FileNotFoundError:  
  plate_fdo = filt3.allfilteredImages(plate_clean, gpij_packed)
  pickle.dump(plate_fdo,open('plate_fdo.pkl','wb'))

print('processing sphere...')
try:
  sph_fs = pickle.load(open('sphere_fs.pkl', 'rb'))
except FileNotFoundError:  
  sph_fs = filt3.allfilteredImages(sphere_ckpfm, gii_packed)
  pickle.dump(sph_fs,open('sphere_fs.pkl','wb'))

try:
  sph_fo = pickle.load(open('sphere_fo.pkl', 'rb'))
except FileNotFoundError:  
  sph_fo = filt3.allfilteredImages(sphere_ckpfm, gij_packed)
  pickle.dump(sph_fo,open('sphere_fo.pkl','wb'))

try:
  sph_fds = pickle.load(open('sphere_fds.pkl', 'rb'))
except FileNotFoundError:    
  sph_fds = filt3.allfilteredImages(sphere_ckpfm, gpii_packed)
  pickle.dump(sph_fds,open('sphere_fds.pkl','wb'))

try:
  sph_fdo = pickle.load(open('sphere_fdo.pkl', 'rb'))
except FileNotFoundError:  
  sph_fdo = filt3.allfilteredImages(sphere_ckpfm, gpij_packed)
  pickle.dump(sph_fdo,open('sphere_fdo.pkl','wb'))

try:
  fgrad_sp = f3.h5_to_dictionaries('savedforcedata.h5')
except OSError:
  print('calculating sphere-plate interaction...')
  plate =f3.surface(plate_clean, fs = all_filtered, fo = plate_fo, fds = fds_filtered, fdo = plate_fdo)
  sphere = f3.sphere(sphere_ckpfm, fs = sph_fs, fo = sph_fo,
                        fds = sph_fds, fdo =  sph_fdo, R = 2.46e-5, dx = 1e-5/512)
  mygrid = f3.make_grid( (256,256), 342, 3)
  fgrad_sp = f3.calcForceData( sphere, plate, mygrid, h_values)

  print('saving sphere-plate interaction...')
  #save the data once your calculations are complete

  f3.fullforcesave('savedforcedata.hdf5', mygrid ,fgrad_sp, h_values)

print('yer done')
