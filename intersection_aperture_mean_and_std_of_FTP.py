#Interconnection aperture
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
#
path_data = "path to the data"
arr_h = np.load(path_data+"mfps_H.npy")
arr_h = arr_h[arr_h[:,3]>=0]
arr_v = np.load(path_data+"mfps_V.npy")
arr_v = arr_v[arr_v[:,3]>=0]
print ("Array H shape: {}".format(arr_h.shape))
print ("Array V shape: {}\n".format(arr_v.shape))
#############################################################################################
path_img = "path to the image"
os.chdir(path_img)
img = io.imread ("fracture_intersections_labelled.tif", plugin='tifffile')
img = np.array(img)
print (img.shape)
if img.ndim!=3:
    raise ValueError('Image is wrong, use another plugin for skimage.io.imread such as imageio or pil')
res = 0.097
lbls = np.unique(img[img>0])
print (lbls)
#
r=6
data=[]
v_check=[]
h_check=[]
lbl0=[]
for i in lbls:
    indx = np.column_stack((np.where(img==i)))
    ap = []
    for k in indx:
        #H
        slc_h = np.where(arr_h[:,0]==k[0])[0]
        ind_h, d_h = indx_nearest_neighbours(k, arr_h[slc_h,:3], 1)
        ind_h = ind_h+np.amin(slc_h)
        #V
        slc_v = np.where(arr_v[:,0]==k[0])[0]
        ind_v, d_v = indx_nearest_neighbours(k, arr_v[slc_v,:3], 1)
        ind_v = ind_v+np.amin(slc_v)
        #Aperture
        H = d_h<=r and ind_h not in h_check
        V =  d_v<=r and ind_v not in v_check
        if H or V:
            if H and V:
                if d_h<d_v:
                    dip_angle = arr_h[ind_h,3].copy()
                    ma = arr_h[ind_h,9].copy()/10
                    ma = np.sin(np.deg2rad(dip_angle))*ma
                    ap_ma = ap_ma = round(0.84*ma,2)
                    ap.append([ap_ma, ap_ma*res])
                    h_check.append(ind_h)
                else:
                    dip_angle = arr_v[ind_v,3].copy()
                    ma = arr_v[ind_v,9].copy()/10
                    ma = np.cos(np.deg2rad(dip_angle))*ma
                    ap_ma = round(0.84*ma,2)
                    ap.append([ap_ma, ap_ma*res])
                    v_check.append(ind_v)
            else:
                if H:
                    dip_angle = arr_h[ind_h,3].copy()
                    ma = arr_h[ind_h,9].copy()/10
                    ma = np.sin(np.deg2rad(dip_angle))*ma
                    ap_ma = ap_ma = round(0.84*ma,2)
                    ap.append([ap_ma, ap_ma*res])
                    h_check.append(ind_h)
                else:
                    dip_angle = arr_v[ind_v,3].copy()
                    ma = arr_v[ind_v,9].copy()/10
                    ma = np.cos(np.deg2rad(dip_angle))*ma
                    ap_ma = round(0.84*ma,2)
                    ap.append([ap_ma, ap_ma*res])
                    v_check.append(ind_v)
    ap = np.array(ap)
    ap_px = round(np.mean(ap[:,0]),2); ap_mm = round(np.mean(ap[:,1]),2)
    std_px = round(np.std(ap[:,0]),2); std_mm = round(np.std(ap[:,1]),2)
    print ("Label {}\thas {}/{} FTPs.\tAperture is:\t{} [{}] px\t {} [{}] mm".format(i, len(ap), len(indx), ap_px, std_px, ap_mm, std_mm))