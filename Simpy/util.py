import numpy as np
import scipy

#useful constants
c = 2.99792458e10
lbol_sun = 3.9e33
loglbol_sun = np.log10(lbol_sun)

def cutdict(target, goodinds):
    for key in target.keys():
        target[key] = target[key][goodinds]
    return

def partial_derivative(func, var=0, point=[]):
    args = point[:]
    def wraps(x):
        args[var] = x
        return func(*args)
    return scipy.misc.derivative(wraps, point[var], dx = 1e-8)

def init_iord_to_fpos(snap):
    iord_to_fpos = np.empty(snap['iord'].max()+1,dtype=np.int64)
    iord_to_fpos[snap['iord']] = np.arange(len(snap))
    return iord_to_fpos


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

def smoothdata(rawdat,nsteps=20,ret_std=False):
    nind = len(rawdat) - len(rawdat)%nsteps
    use = np.arange(nind)
    rawdat = rawdat[use].reshape((nind/nsteps,nsteps))
    meandat = rawdat.mean(axis=1)
    if ret_std is False:
        return meandat
    else:
        std = rawdat.std(axis=1)
        return meandat, std

