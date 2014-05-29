# This script computes Haar Transform of 1D signal (eg. Time series)
# i.e. input, upto n levels where, at the n'th level there is only one Scaling
# and Wavelet coeficient, each.
#
# INPUT: 1) 1D array of size 2^n (n is an integer).
#	 2) Number of levels. (an integer greater than 1)
#
# OUTPUT: 2D list containing Scalig Cef.s and Wavelet Coef.s of INPUT in the
# 	    following form:
#		Scaling Coef: a1, a2, a3, ...
#		Wavelet Coef: d1, d2, d3, ...
#		n: Number of levels upto which Haar transform is calculated
#		   starting from 0 level (at 0-level, we have the actual series)
#		Each a's and d's are individual lists which are contained in
#		Scaling Coef and Wavelet Coef respectively, which are lists in
#		themselvels. (List of lists)

import numpy as np

def calc_wavelet_coef(scaling_coef):
	return (scaling_coef[0::2]-scaling_coef[1::2])/np.sqrt(2)

def calc_scaling_coef(scaling_coef):
	return (scaling_coef[0::2]+scaling_coef[1::2])/np.sqrt(2)

def haar_transform(series, levels=False):
	s_coef = []
	w_coef = []
	w_coef.insert(0, calc_wavelet_coef(series))
	s_coef.insert(0, calc_scaling_coef(series))

	n = 1

	if(levels == False):
		while(len(s_coef[0] != 1)):
			w_coef.insert(0, calc_wavelet_coef(s_coef[0]))
			s_coef.insert(0, calc_scaling_coef(s_coef[0]))
			n = n+1

	else:
		while(n < levels and len(s_coef[0]) != 1):
			w_coef.insert(0, calc_wavelet_coef(s_coef[0]))
			s_coef.insert(0, calc_scaling_coef(s_coef[0]))
			n = n+1

	return s_coef, w_coef, n-1
