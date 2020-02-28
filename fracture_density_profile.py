#Fracture density function
import numpy as np
import os
from skimage import io
def sys():
    import sys
    sys.exit()
def indx_nearest_neighbours(point,arr, n):
    """Return the index of the closest element in arr"""
    import scipy
    tree=scipy.spatial.cKDTree(arr)
    dist,ind=tree.query(point, k=n)
    return ind, dist
def mm2inch(*tupl):
    inch = 25.4
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)
#
path_data = "path to the data"
arr_h0 = np.load(path_data+"mfps_H.npy")
arr_h = arr_h0[arr_h0[:,3]>=0]; arr_h=arr_h[:,:3]
arr_v0 = np.load(path_data+"mfps_V.npy")
arr_v = arr_v0[arr_v0[:,3]>=0]; arr_v=arr_v[:,:3]
print ("Array H shape: {}".format(arr_h.shape))
print ("Array V shape: {}\n".format(arr_v.shape))
#############################################################################################
path_img = "path to the image"
os.chdir(path_img)
img = io.imread ("fracture_separated.tif", plugin='tifffile')
img = np.array(img)
if img.ndim!=3:
    raise ValueError('Image is wrong, use another plugin for skimage.io.imread such as imageio or pil')
from glob import glob
data0 = sorted (glob("*.txt"))[:-1]
data0 = np.array([[int(el[:2]), int(el[-7:-4])] for el in data0], dtype=int)
data0 = data0[data0[:,1]<999] #0 is the label value and 1 is the dip direction
del img
#
img = io.imread (path_data+"img.tif", plugin='tifffile')
img = np.array(img)
if img.ndim!=3:
    raise ValueError('Image is wrong, use another plugin for skimage.io.imread such as imageio or pil')
res = 0.09
slc_area = float(len(img[0].nonzero()[0]))*res
vol_area = float(slc_area*img.shape[0])*res
del img
#
path = "path to the image"
os.chdir(path)
#
img = io.imread ("fracture_separated.tif", plugin='tifffile')
img = np.array(img)
if img.ndim!=3:
    raise ValueError('Image is wrong, use another plugin for skimage.io.imread such as imageio or pil')
lbls = np.unique(img[img>0])
#
r=6
data=[]
v_check=[]
h_check=[]
lbl0=[]
for i in data0:
    indx = np.column_stack((np.where(img==i[0])))
    i_ftp = []
    for k in indx:
            if 45<=i[1]<=135 or 225<=i[1]<=315:
                slc = np.where(arr_h[:,0]==k[0])[0]
                ind, d = indx_nearest_neighbours(k, arr_h[slc,:3], 1)
                ind = ind+np.amin(slc)
                if ind not in h_check and d<=r:
                    h_check.append(ind)
                    i_ftp.append(arr_h[ind].astype(int))
                else:
                    slc = np.where(arr_v[:,0]==k[0])[0]
                    ind, d = indx_nearest_neighbours(k, arr_v[slc,:3], 1)
                    ind = ind+np.amin(slc)
                    if ind not in v_check and d<=r:
                        v_check.append(ind)
                        i_ftp.append(arr_v[ind].astype(int))
            else:
                slc = np.where(arr_v[:,0]==k[0])[0]
                ind, d = indx_nearest_neighbours(k, arr_v[slc,:3], 1)
                ind = ind+np.amin(slc)
                if ind not in v_check and d<=r:
                    v_check.append(ind)
                    i_ftp.append(arr_v[ind].astype(int))
                else:
                    slc = np.where(arr_h[:,0]==k[0])[0]
                    ind, d = indx_nearest_neighbours(k, arr_h[slc,:3], 1)
                    ind = ind+np.amin(slc)
                    if ind not in h_check and d<=r:
                        h_check.append(ind)
                        i_ftp.append(arr_h[ind].astype(int))
    i_ftp = np.array(i_ftp)
    if len(i_ftp)>0:
        lbl0.append ([i[0], len(i_ftp)])
        data.append(i_ftp)
    print ("Analyzing {} with {} elements:\t\t{} found.".format(i[0], len(indx), len(i_ftp)))
shape = img.shape
del img
img = np.zeros(shape)
lbl0 = np.array(lbl0)
print ("New image created.")
for n, i in enumerate(data):
    crds = np.array(i)
    if len(crds)>1:
        img[crds[:,0],crds[:,1],crds[:,2]]=lbl0[n,0]
print ("Imaged filled with FTPs.")
########################################################################################################################
lbls = np.unique(img[img>0])
prof_all = ([(float(len(el.nonzero()[0]))/slc_area)*res for el in img])
print ("P32 is {}".format(np.sum(prof_all)))
prof_ind = []
for n, lbl in enumerate(lbls):
    prof_lbl = []
    for slc in img:
        prof_lbl.append((float(len(np.where(slc==lbl)[0]))/slc_area)*res)
    prof_lbl = np.array(prof_lbl)
    prof_lbl = np.array([lbls[n],prof_lbl])
    prof_ind.append(prof_lbl)
prof_ind = np.array(prof_ind)
#PLOTS
import matplotlib.pyplot as plt
plt.subplots(figsize=mm2inch(194.865, 92.4))#; plt.grid()
mrksp='s'
mrksz=2.2 #MarkerSize
lw=0.8
plt.plot(prof_all, 'k', marker=mrksp, markersize=mrksz, linewidth=lw)
#Density profile plot
col = ['red', 'red', 'red', 'red', 'red', 'blue', 'red', 'red', 'red', 'red', 'red', 'red', 'deepskyblue', 'red', 'red', 'purple', 'red', 'red', 'red', 'red', 'red', 'red', 'green', 'red', 'red', 'gold', 'red', 'red', 'red',]
col = [w.replace('red', 'grey') for w in col]
for n, el in enumerate(prof_ind):
    print (el[0],np.sum(el[1]))
    prof = el[1].copy()
    prof = prof.astype(float)
    elements = prof.nonzero()[0]
    if elements[0]>0:
        prof[:elements[0]-1]=np.nan
    if elements[-1]<len(prof_all)-1:
        prof[elements[-1]+2:]=np.nan
    plt.plot(prof, marker=mrksp, markersize=mrksz, c=col[n], linewidth=lw)
plt.xlabel('Image slice (z-axis)', fontsize=8); plt.ylabel('Fracture trace (mm/mm^2)', fontsize=8)
plt.show()
del img