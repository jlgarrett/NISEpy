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
gii_packed = filt3.complete_filtering(h_values,243, 'fs', float(n8B_plate['cleannote']['dx']))
gpii_packed = filt3.complete_filtering(h_values,243, 'fds', float(n8B_plate['cleannote']['dx']))
gij_packed = filt3.complete_filtering(h_values,243, 'fo', float(n8B_plate['cleannote']['dx']))
gpij_packed = filt3.complete_filtering(h_values,243, 'fdo', float(n8B_plate['cleannote']['dx']))


pickle.dump(gii_packed,open('gii_packed.pkl','wb') )
pickle.dump(gpii_packed,open('gpii_packed.pkl','wb') )
pickle.dump(gij_packed,open('gij_packed.pkl','wb') )
pickle.dump(gpij_packed,open('gpij_packed.pkl','wb') )

print('processing plate...')
all_filtered = filt3.allfilteredImages(plate_clean, gii_packed)
fds_filtered = filt3.allfilteredImages(plate_clean, gpii_packed)
plate_fo = filt3.allfilteredImages(plate_clean, gij_packed)
plate_fdo = filt3.allfilteredImages(plate_clean, gpij_packed)

pickle.dump(all_filtered,open('plate_filtered.pkl','wb'))
pickle.dump(fds_filtered,open('plate_filtered_fds.pkl','wb') )
pickle.dump(plate_fo,open('plate_fo.pkl','wb'))
pickle.dump(plate_fdo,open('plate_fdo.pkl','wb'))

print('processing sphere...')
sph_fs = filt3.allfilteredImages(sphere_ckpfm, gii_packed)
sph_fo = filt3.allfilteredImages(sphere_ckpfm, gij_packed)
sph_fds = filt3.allfilteredImages(sphere_ckpfm, gpii_packed)
sph_fdo = filt3.allfilteredImages(sphere_ckpfm, gpij_packed)


pickle.dump(sphere_fs,open('sphere_fs.pkl','wb'))
pickle.dump(sphere_fo,open('sphere_fo.pkl','wb'))
pickle.dump(sphere_fds,open('sphere_fds.pkl','wb'))
pickle.dump(sphere_fdo,open('sphere_fdo.pkl','wb'))

print('calculating sphere-plate interaction...')
plate =f3.surface(plate_clean, fs = all_filtered, fo = plate_fo, fds = fds_filtered, fdo = plate_fdo)
sphere = f3.sphere(sphere_ckpfm, fs = sph_fs, fo = sph_fo,
                        fds = sph_fds, fdo =  sph_fdo, R = 3.2e-5, dx = 1e-5/512)
mygrid = f3.make_grid( (256,256), 342, 3)
fgrad_sp = (f3.calcForceData( sphere, plate, mygrid, h_values))


print('saving sphere-plate interaction...')
import savingdict
#save the data once your calculations are complete
f3.quickforcesave('savedforcedata.hdf5', mygrid ,fgrad_sp)

print('yer done')