import numpy as np

#useful constants
c = 2.99792458e10
lbol_sun = 3.9e33
loglbol_sun = np.log10(lbol_sun)

def histogram(a,inbins,weights=None):
		if (np.size(np.shape(inbins))!=2) | (np.shape(inbins)[1] != 2):
				print "bins have wrong shape! shoudl be: [[a,b],[b,c],[c,d]...]"
				return
		edges = np.array(inbins.ravel())
		edges.sort()
		edges = np.unique(edges)
		n, edges= np.histogram(a,bins=edges,weights=weights)
		hist = np.zeros(len(inbins))
		binInds = np.searchsorted(edges,inbins)

		for i in range(len(inbins)):
				hist[i] = np.sum(n[binInds[i,0]:binInds[i,1]])

		return hist

def timeweightedAve(x,dt):
		T = np.sum(dt)
		mindt = dt[(dt>0)].min()
		xsum = np.sum(x*dt)
		ave = xsum/T
		variance = np.average((x-ave)**2, weights=dt/mindt)
		return ave, np.sqrt(variance), xsum

#convert monochromatic AB absolute magnitude to log(luminosity [ergs/s])
def mcABconv(mag,nu):
        C = 20.638
        return -(2./5.)*mag + C + np.log10(nu)

