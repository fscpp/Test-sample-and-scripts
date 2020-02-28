import numpy as np
def sys():
    import sys
    sys.exit()
###########################################################################################################################################
def mm2inch(*tupl):
    inch = 25.4
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)
##############################################################################################################################################################
arr_h = np.genfromtxt ("path to \\mfps_H.txt", delimiter=',', skip_header=1)
arr_v = np.genfromtxt ("path to \\mfps_V.txt", delimiter=',', skip_header=1)
matrix=4770
air=2200
fwhm_val = ((matrix-air)/2)+air
def params (arr, fwhm_val):
    """Return parameters and correct relative value to the absolute one by using the dip angle"""
    import numpy as np
    arr = arr[np.where((arr[:,7]>0) & (arr[:,5]<fwhm_val))] #7 is FWHM and 5 is the MFP value
    print ("{} elements for the calibration.\n".format(arr.shape[0]))
    da = arr[:,3].astype(int)
    dd = arr[:,4].astype(int)
    val = arr[:,5].astype(int)
    res = arr[:,6].astype(int)
    fwhm = arr[:,7].astype(float)/10
    ph = arr[:,8].astype(float)/10
    ma = arr[:,9].astype(float)/10
    return fwhm, ma, ph, da, dd, val, res
fwhm_h, ma_h, ph_h, da_h, dd_h, val_h, res_h = params (arr_h, fwhm_val)
fwhm_v, ma_v, ph_v, da_v, dd_v, val_v, res_v = params (arr_v, fwhm_val)
###############################################################################
def true_width (fwhm, ma, da, side='h'):
    """From relative to true widths using dip angle and dip direction."""
    import numpy as np
    fwhm_corr = fwhm.copy()
    ma_corr = ma.copy()
    for n, i in enumerate(da):
        if side=='h':
            fwhm_corr[n] = np.sin(np.deg2rad(i))*fwhm[n]
            ma_corr[n] = np.sin(np.deg2rad(i))*ma[n]
        else:
            fwhm_corr[n] = np.cos(np.deg2rad(i))*fwhm[n]
            ma_corr[n] = np.cos(np.deg2rad(i))*ma[n]
    return fwhm_corr, ma_corr
fwhm_h, ma_h = true_width (fwhm_h, ma_h, da_h, side='h')
fwhm_v, ma_v = true_width (fwhm_v, ma_v, da_v, side='v')
###############################################################################
def ransac_2d(x, y):  
    from skimage.measure import LineModelND, ransac
    from sklearn.metrics import mean_squared_error, r2_score
    import numpy as np
    coords = np.column_stack([x, y])
    model, inliers = ransac(coords, LineModelND, min_samples=5, residual_threshold=1, max_trials=500)
    line_x = np.linspace(0, x.max()*10,len(x))
    line_y_robust = model.predict_y(line_x)
    outliers = inliers == False
    line_y_robust =line_y_robust-line_y_robust[0]
    #Equation of the line
    m = round((line_y_robust[-1]-line_y_robust[0])/(line_x[-1]-line_x[0]),2)
    b = round(line_y_robust[0],2)
    print ("Line equation:\ty={}*x+{}".format(m,b))
    #RMSE and R2
    line_y_R2 = model.predict_y(x)
    rmse = np.sqrt(mean_squared_error(y,line_y_R2))
    r2 = r2_score(y,line_y_R2)
    print ("RMSE:\t{}\nR2:\t{}".format(round(rmse,2), round(r2,2)))
    return line_x, line_y_robust, outliers
###############################################################################
fwhm = np.hstack((fwhm_h, fwhm_v)); del fwhm_h; del fwhm_v
ma = np.hstack((ma_h, ma_v)); del ma_h; del ma_v
ph = np.hstack((ph_h, ph_v)); del ph_h; del ph_v
val = np.hstack((val_h, val_v)); del val_h; del val_v
import matplotlib.pyplot as plt
mrk_sp = 's'; mrk_sz = 5; lw=0.7; cl = 'lightgrey'
f_lim=0,30,5; m_lim=0,30,5; p_lim=45,165,20
plt.figure(figsize=mm2inch(92.4, 72.4))
plt.subplot(111)
plt.xlabel ('MA', fontsize=8); plt.ylabel ('FWHM (px)', fontsize=8)
plt.xlim(m_lim[:-1]); plt.xticks(range(m_lim[0],m_lim[1]+1,m_lim[2]), fontsize=8)
plt.ylim(f_lim[:-1]); plt.yticks(range(f_lim[0],f_lim[1]+1,f_lim[2]), fontsize=8)
n = 1780
x = ma[val<n].copy()
y = fwhm[val<n].copy()
plt.scatter(x, y
            , c=cl, marker=mrk_sp, s=mrk_sz, edgecolor='k', linewidth=lw); plt.grid()
x_fit, y_fit, out = ransac_2d(x, y)
plt.plot(x_fit, y_fit, '--g', linewidth=1.2)
plt.tight_layout()
plt.show()