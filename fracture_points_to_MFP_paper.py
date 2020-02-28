#Connect fracture points to MFPs
import numpy as np
import os
from skimage import io
def sys():
    import sys
    sys.exit()
def myround(x, base, c=None):
    n = base * round(float(x)/base)
    if c is None:
        n = n if (n-base)>0 else n+base
    return n
def mm2inch(*tupl):
    inch = 25.4
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)
#
def indx_nearest_neighbours(point,arr, n):
    """Return the index of the closest element in arr"""
    import scipy
    tree=scipy.spatial.cKDTree(arr)
    dist,ind=tree.query(point, k=n)
    return ind, dist
#
path_data = "C:\\Users\\Franco\\Desktop\\PhD_Paper\\test_auckland\\Minima_Analysis\\"
arr_h = np.load(path_data+"mfps_H.npy")
arr_h = arr_h[arr_h[:,3]>=0]
arr_v = np.load(path_data+"mfps_V.npy")
arr_v = arr_v[arr_v[:,3]>=0]
print (arr_h.shape)
print (arr_v.shape)
#############################################################################################
path_img = "C:\\Users\\Franco\\Desktop\\PhD_Paper\\test_auckland\\Data_fractures_separated\\"
os.chdir(path_img)
img = io.imread ("fracture_separated.tif", plugin='tifffile')
img = np.array(img)
if img.ndim!=3:
    raise ValueError('Image is wrong, use another plugin for skimage.io.imread such as imageio or pil')
#import matplotlib.pyplot as plt
#img = img.astype(float)
#img[img==0]=np.nan
#plt.subplot(111)
#for i in np.unique(img[img>0]):
#    img[img==i]=len(np.where(img==i)[0])
#print (np.unique(img[img>0]))
#plt.imshow(img[51], 'jet')
#sys()
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
                    ph = arr_h[ind,8].copy()/10
                    ma = np.sin(np.deg2rad(dip_angle))*ma
                    ap_ma = ap_ma = round(0.84*ma,2)
                    ap_ph = round(0.1*ph-4.3,2)
                    ap.append([ap_ma, ap_ph])
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
                        ph = arr_v[ind,8].copy()/10
                        ma = np.cos(np.deg2rad(dip_angle))*ma
                        ap_ma = round(0.84*ma,2)
                        ap_ph = round(0.1*ph-4.3,2)
                        ap.append([ap_ma, ap_ph])
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
                    ap_ph = round(0.1*ph-4.3,2)
                    ap.append([ap_ma, ap_ph])
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
                        ph = arr_h[ind,8].copy()/10
                        ma = np.sin(np.deg2rad(dip_angle))*ma
                        ap_ma = round(0.84*ma,2)
                        ap_ph = round(0.1*ph-4.3,2)
                        ap.append([ap_ma, ap_ph])
                        ornt.append([dip_angle, arr_h[ind,4]])
                        h_check.append(ind)
                        n_h+=1
        ap = np.array(ap, dtype=float, copy=False)
        ma = np.ravel(ap[:,0]); ap_ma_all = ma[ma>0.].copy(); ap_ma_mean = np.mean(ap_ma_all); s_ma = round(np.std(ap_ma_all),2)
        ph = np.ravel(ap[:,1]); ap_ph_all = ph[ph>0.].copy(); ap_ph_mean = np.mean(ap_ph_all); s_ph = round(np.std(ap_ph_all),2)
        ornt = np.array(ornt, dtype=int); da = ornt[:,0].copy(); dd = ornt[:,1].copy()
        print ("{} measurements made ({} H and {} V):\nMean MA:\t{} [{} data points]\tsd={}\nMean PH:\t{} [{} data points]\tsd={}\n"
               .format(len(ap), n_h, n_v, round(ap_ma_mean,2), len(ap_ma_all), s_ma, round(np.mean(ap_ph_mean),2), len(ap_ph_all), s_ph))
#        import matplotlib.pyplot as plt
#        plt.figure(figsize=mm2inch(92.4/2, 92.4/2))
#        plt.subplot(111); b=15
#        n2, b, _ = plt.hist(ph, b, facecolor='r', edgecolor='k', linewidth=0.75, density=True)
#        n, b, _ = plt.hist(ma, b, facecolor='g', edgecolor='k', linewidth=0.75, density=True, alpha=0.65)
#        plt.xticks(np.arange(myround(np.amin(np.hstack((ma,ph))),5,c=1), myround(np.amax(np.hstack((ma,ph))),5)+1, 5), fontsize=6)
#        plt.ylim(0,myround(np.amax(n),0.05))
#        plt.yticks(np.arange(0., myround(np.amax(np.hstack((n,n2))),0.05)+0.01, 0.05), fontsize=6)
#        plt.subplot(133)#; plt.grid()
#        plt.hist(da, bins=30, density=True)
#        plt.hist(da, bins=30, density=True)
#        plt.xlim(0,90); plt.xticks(range(0,91,30), fontsize=6)
#        plt.ylim(0,360); plt.yticks(range(0,361,45), fontsize=6)
        #
#        plt.tight_layout()
#        plt.show()
        ap_all.append(ap)