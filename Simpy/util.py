import numpy as np
import scipy
import pynbody

#useful constants
c = pynbody.array.SimArray(2.99792458e10,'cm s**-1')
lbol_sun = 3.9e33
G = pynbody.array.SimArray(6.67259e-8,'cm**3 s**-2 g**-1')
M_sun_g = 1.988547e33
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
    iord_to_fpos = np.empty(snap['iord'].max()+1,dtype=np.dtype('int64'))
    iord_to_fpos[snap['iord']] = np.linspace(0,len(snap)-1,len(snap)).astype(np.int64)
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

def get_rates_by_z(zin,s,range=[0,20],nbins=40,units='yr**-1',weights=None):
    import numpy as np
    from . import cosmology
    n,zbins = np.histogram(zin,weights=weights,range=range,bins=nbins)
    tedges = pynbody.array.SimArray([cosmology.getTime(z,s) for z in zbins],'Gyr')
    tedges = tedges.in_units(units)
    dt = np.abs(tedges[0:-1]-tedges[1:])
    n = n/dt
    return n, zbins


