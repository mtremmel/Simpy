import numpy as np
import pynbody


def specAM(beta,type='stars',mrange=None, munits='Msol'):
	if mrange is None:
		mrange = [1e9,1e11]
	#D. Obreschkow and K. Glazebrook 2014
	if type not in ['stars','baryons']:
		print 'WARNING do not understand type given. must be either stars or gas'
		return
	if type == 'stars':
		k = 0.89
		a = 0.94
		g = 7.03
	else:
		k = 0.77
		a = 0.98
		g = 6.65
	dlm = (np.log10(mrange[1]) - np.log10(mrange[0])) / 1000.
	marray = pynbody.array.SimArray(10**np.arange(np.log10(mrange[0]),np.log10(mrange[1]),dlm), munits)
	j = k * 1e3 * np.exp(-1.0*g*beta) * (marray.in_units('Msol')/1e10)**a
	return marray, pynbody.array.SimArray(j,'kpc km s**-1')



def betaSpecAM(type='stars',sAMrange=None, sAMunits='kpc km s**-1 Msol**-1'):
	#D. Obreschkow and K. Glazebrook 2014
	if sAMrange is None:
		sAMrange = [1e-9,1e-6]
	if type not in ['stars','baryons']:
		print 'WARNING do not understand type given. must be either stars or gas'
		return

	if type == 'stars':
		k1 = -0.34
		k2 = -0.04
	else:
		k1 = -0.3
		k2 = -0.01

	dlsAM = (np.log10(sAMrange[1]) - np.log10(sAMrange[0])) / 1000.
	sAMarray =  pynbody.array.SimArray(10**np.arange(np.log10(sAMrange[0]),np.log10(sAMrange[1]),dlsAM), sAMunits)
	beta = k1 * np.log10(sAMarray/1e-7) + k2
	return sAMarray, beta

