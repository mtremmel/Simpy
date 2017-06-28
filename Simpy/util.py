import numpy as np
import scipy
import pynbody

#useful constants
c = pynbody.array.SimArray(2.99792458e10,'cm s**-1')
lbol_sun = 3.9e33
G = pynbody.array.SimArray(6.67259e-8,'cm**3 s**-2 g**-1')
M_sun_g = 1.988547e33
loglbol_sun = np.log10(lbol_sun)
mh = pynbody.array.SimArray(1.6726219e-24, 'g')
kb = pynbody.array.SimArray(1.380658e-16, 'erg K**-1')

def L_edd(mass):
    return mass*lbol_sun * 3.2e4

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

def wrap(relpos,scale,boxsize=25):
	bphys = boxsize*scale
	bad = np.where(np.abs(relpos) > bphys/2.)
	relpos[bad] = -1.0 * (relpos[bad]/np.abs(relpos[bad])) * np.abs(bphys/2. - np.abs(relpos[bad]))
	return


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
    newdat = rawdat[use].reshape((nind/nsteps,nsteps))
    meandat = newdat.mean(axis=1)
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
    dt = np.abs(tedges[0:-1]-tedges[1:])
    n = n/dt
    return n.in_units(units), zbins

def find_sats(Mvir, cen, Rvir, strict = False):
    sub = np.zeros(len(Mvir))
    for i in range(len(Mvir)):
        if strict is False:
            int = np.where((np.sqrt(np.sum((cen[i] - cen)**2,axis=1))<Rvir)&(np.sum((cen[i] - cen)**2,axis=1)>0)&(Mvir[i]<Mvir))
        else:
            int = np.where((np.sqrt(np.sum((cen[i] - cen)**2,axis=1))<Rvir)&(np.sum((cen[i] - cen)**2,axis=1)>0))
        if len(int[0])>0:
            sub[i] = 1
    return sub



