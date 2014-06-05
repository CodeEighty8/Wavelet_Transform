# This script applies Haar Transform on Daily Spatial means for various
# regions derived from TRMM's 3b42 data set.
# This script performs multi resolution analysis.
#
# INPUT = 1) s = Input time series of length 2^n, where n is an integer
#         2) LEVELS = Depth of MRA
#
# OUTPUT = 2D array of dimention (LEVELS, length_of_time_series)
#          Here, the n-1th row corresponds to time series representation of
#          nth level scaling coefficients and subsequent rows (i.e. n-2, n-3
#          , ..., 0) gives time series representation wavelet coefficients
#
# Ths structure of the output is as follows:
#                                    0 1 2 ... (2^n)-1
#          LEVEL-k     scaling coeff . . . ... .
#          LEVEL-k     wavelet coeff . . . ... .
#          LEVEL-(k-1) wavelet coeff . . . ... .
#          LEVEL-(k-2) wavelet coeff . . . ... .
#          .
#          .
#          .
#          LEVEL-1     wavelet coeff . . . ... .
#
#          Adding column-wise will give back the reconstructed signal

import numpy as np
import pylab as pl

# Import the Haar Transform function
from Haar_Transform import haar_transform as ht

def MRA(s, LEVELS):
    level = range(1, LEVELS+1)
    scale = [2**l for l in level]
    
    s_coef, w_coef, n = ht(s, LEVELS)

    ts_sc = np.zeros(s.size, dtype=float)
    temp = 0

    while temp < scale[LEVELS - 1]:
        ts_sc[temp::scale[LEVELS-1]] = s_coef[0][:]/np.power(np.sqrt(2), LEVELS)
        temp += 1
    
    ts_wc = np.zeros((LEVELS, s.size), dtype=float)

    for l in level:
        temp = 0
        while temp < scale[LEVELS-l]:
            if temp < scale[LEVELS-l]/2:
                ts_wc[LEVELS - l][temp::scale[LEVELS-l]] = w_coef[l-1][:]/np.power(np.sqrt(2), l)
                #print temp, scale[LEVELS-l], ts_wc[LEVELS - l][temp::scale[LEVELS-l]].shape, w_coef[l-1].shape
                temp += 1
            else:
                ts_wc[LEVELS - l][temp::scale[LEVELS-l]] = (-w_coef[l-1][:])/np.power(np.sqrt(2), l)
                #print temp, scale[LEVELS-l], ts_wc[LEVELS - l][temp::scale[LEVELS-l]].shape, w_coef[l-1].shape
                temp += 1
    
    return np.vstack((ts_wc, ts_sc))



if __name__ == "__main__":
    # Use your time series here
    s = np.load('/home/venu/programming/CAOS/Data/3b42_v7_daily_quarter_degree/trmm3b42v7_daily/Python_data_saves/daily.spatial.mean.Central.India.npy')['mean'][0:4096]

    s = (s - s.min())/(s.max() - s.min())
    
    LEVELS = 3
    
    if LEVELS == 0:
        LEVELS = np.int(np.log2(s.size))

    k = MRA(s, LEVELS)
    
    fig, ax = pl.subplots(LEVELS+1, sharex=True)
    
    for l in range(1, LEVELS+1):
        ax[l].plot(k[LEVELS-l])
        ax[l].set_ylabel('D'+str(LEVELS+1-l))
        ax[l].grid()
        ax[l].set_xlim(0, s.size-1)
        ax[l].set_ylim(-1, 1)

        ax[0].plot(k[LEVELS], 'b')
        ax[0].set_ylabel('A'+str(LEVELS))
        ax[0].grid()

        #ax[1].set_xlim(0, s.size-1)
        #ax[1].grid()

        #ax[0].plot(s, 'b')
        #ax[0].set_ylabel('s ($mm\ day^{-1}$)')
        #ax[0].set_xlim(0, s.size-1)
        #ax[0].grid()

    #pl.tight_layout()
    pl.show()
