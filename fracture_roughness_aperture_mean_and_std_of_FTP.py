#Connect fracture points to MFPs
import numpy as np
import os
from skimage import io
#
def indx_nearest_neighbours(point,arr, n):
    """Return the index of the closest element in arr"""
    import scipy
    tree=scipy.spatial.cKDTree(arr)
    dist,ind=tree.query(point, k=n)
    return ind, dist
#
path_data = "path to data"
arr_h = np.load(path_data+"mfps_H.npy")
arr_h = arr_h[arr_h[:,3]>=0]
arr_v = np.load(path_data+"mfps_V.npy")
arr_v = arr_v[arr_v[:,3]>=0]
#############################################################################################
path_img = "path to image"
os.chdir(path_img)
img = io.imread ("fracture_separated.tif", plugin='tifffile')
img = np.array(img)
if img.ndim!=3:
    raise ValueError('Image is wrong, use another plugin for skimage.io.imread such as imageio or pil')
from glob import glob
data = sorted (glob("*.txt"))[:-1]
data = np.array([[int(el[:2]), int(el[-7:-4])] for el in data], dtype=int)
data = data[data[:,1]<999] #0 is the label value and 1 is the dip direction
#Start the analysis
thr_d = 6
res = 12
h_check = []
v_check = []
ap_all = []
for i in data:
    indx = np.column_stack((np.where(img==i[0])))
    if len(indx)>7000: #Analyze only the 5 big fractures
        print ("Analyzing label {} with {} elements.".format(i[0], len(indx)))
        ap = []
        ornt = []
        n_h=0
        n_v=0
        for k in indx:
            if 45<=i[1]<=135 or 225<=i[1]<=315:
                slc = np.where(arr_h[:,0]==k[0])[0]
                ind, d = indx_nearest_neighbours(k, arr_h[slc,:3], 1)
                ind = ind+np.amin(slc)
                if ind not in h_check and d<=thr_d:
                    dip_angle = arr_h[ind,3].copy()
                    ma = arr_h[ind,9].copy()/10
                    ma = np.sin(np.deg2rad(dip_angle))*ma
                    ap_ma = ap_ma = round(0.84*ma,2)
                    ap.append([ap_ma])
                    ornt.append([dip_angle, arr_h[ind,4]])
                    h_check.append(ind)
                    n_h+=1
                else:
                    slc = np.where(arr_v[:,0]==k[0])[0]
                    ind, d = indx_nearest_neighbours(k, arr_v[slc,:3], 1)
                    ind = ind+np.amin(slc)
                    if ind not in v_check and d<=thr_d:
                        dip_angle = arr_v[ind,3].copy()
                        ma = arr_v[ind,9].copy()/10
                        ma = np.cos(np.deg2rad(dip_angle))*ma
                        ap_ma = round(0.84*ma,2)
                        ap.append([ap_ma])
                        ornt.append([dip_angle, arr_v[ind,4]])
                        v_check.append(ind)
                        n_v+=1
            else:
                slc = np.where(arr_v[:,0]==k[0])[0]
                ind, d = indx_nearest_neighbours(k, arr_v[slc,:3], 1)
                ind = ind+np.amin(slc)
                if ind not in v_check and d<=thr_d:
                    dip_angle = arr_v[ind,3].copy()
                    ma = arr_v[ind,9].copy()/10
                    ph = arr_v[ind,8].copy()/10
                    ma = np.cos(np.deg2rad(dip_angle))*ma
                    ap_ma = round(0.84*ma,2)
                    ap.append([ap_ma])
                    ornt.append([dip_angle, arr_v[ind,4]])
                    v_check.append(ind)
                    n_v+=1
                else:
                    slc = np.where(arr_h[:,0]==k[0])[0]
                    ind, d = indx_nearest_neighbours(k, arr_h[slc,:3], 1)
                    ind = ind+np.amin(slc)
                    if ind not in h_check and d<=thr_d:
                        dip_angle = arr_h[ind,3].copy()
                        ma = arr_h[ind,9].copy()/10
                        ma = np.sin(np.deg2rad(dip_angle))*ma
                        ap_ma = round(0.84*ma,2)
                        ap.append([ap_ma])
                        ornt.append([dip_angle, arr_h[ind,4]])
                        h_check.append(ind)
                        n_h+=1
        ap = np.array(ap, dtype=float, copy=False)
        ma = np.ravel(ap[:]); ap_ma_all = ma[ma>0.].copy(); ap_ma_mean = np.mean(ap_ma_all); s_ma = round(np.std(ap_ma_all),3)
        ap_ma_all_mm = ap_ma_all*0.097
        ap_ma_mean_mm = np.mean(ap_ma_all_mm); s_ma_mm = round(np.std(ap_ma_all_mm),3)
        ornt = np.array(ornt, dtype=int); da = ornt[:,0].copy(); dd = ornt[:,1].copy()
        print ("{} measurements made ({} H and {} V):\n[{} data points]\nMean MA:\t{}\tsd={}\nMean MA:\t{}\tsd={}\nrugosity={}\n\n"
               .format(len(ap), n_h, n_v, len(ap_ma_all), round(ap_ma_mean,3), s_ma, round(ap_ma_mean_mm,3), s_ma_mm, round(s_ma_mm/ap_ma_mean_mm,3)))
        ap_all.append(ap)