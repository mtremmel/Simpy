import numpy as np

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
	mindt = dt.min()
	ave = np.sum(x*dt/T)
	variance = np.average((x-ave)**2, weights=dt/mindt)
	return (ave, np.sqrt(variance))
