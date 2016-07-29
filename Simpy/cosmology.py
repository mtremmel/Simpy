from pynbody.analysis import pkdgrav_cosmo as cosmo
import pynbody.analysis.cosmology as other_cosmo
import pynbody.array as arr
import numpy as np
from scipy import optimize as opt
from scipy import integrate
from . import util

def getTime(z,sim):
        c = cosmo.Cosmology(sim=sim)
        return other_cosmo.age(sim,z=0)*c.Exp2Time(1.0 / (1+z))/c.Exp2Time(1)

def getScaleFactor(times,s):
        redshift = np.zeros(np.size(times))
        ntimes = np.size(times)
        for tt in range(np.size(times)):
                if tt%100==0:
                        print tt/np.float(np.size(times)) * 100, '% done'
                def func(z):
                        return getTime(z,s) - times.in_units('Gyr')[tt]
                try: redshift[tt] = opt.newton(func,0)
                except:
                        try:
                                print "trying again"
                                redshift[tt] = opt.newton(func,20)
                        except:
                                print "ERROR did not converge", times[tt],tt
                                redshift[tt] = -1
        scaleFac = 1./(1+redshift)
        return scaleFac, redshift

def HoverH0(z, omegaM, omegaL):
        return np.sqrt(omegaM*(1+z)**3 + omegaL)

def comoving_dist(z, omegaM, omegaL, h):
        dH = util.c/(1.0e5*h*100.)

        def func(z):
                return 1./HoverH0(z,omegaM, omegaL)

        def get_d(x):
                return dH * integrate.quad(func,0,x)[0]

        if isinstance(z, np.ndarray) or isinstance(z, list):
                return np.array(map(get_d, z))
        else:
                return get_d(z)

def lum_distance(z, omegaM, omegaL, h):
        dc = comoving_dist(z,omegaM, omegaL, h)
        return dc*(1+z)

def event_count(N, z, omegaM, omegaL, h):
        dc = comoving_dist(z, omegaM, omegaL, h)
        return 4*np.pi*util.c*1.02269032e-17*N*dc**2

def get_hmf_data(simpath,**kwargs):
        import pynbody
        sim = pynbody.load(simpath)
        sim.properties['sigma8'] = 0.77
        mass, sig, phi = pynbody.analysis.halo_mass_function(sim,pspec=pynbody.analysis.hmf.PowerSpectrumCAMBLive,**kwargs)
        return np.log10(mass),phi

def get_hmf_data_all(**kwargs):
        f = open('files.list','r')
        hmf = {}
        for l in f:
                name = l.strip('\n')
                print name
                hmf[name] = make_hmf(name,**kwargs)
        return hmf

class HMF(object):
        def __init__(self,**kwargs):
                self.data = get_hmf_data_all(**kwargs)
                self.hmf = {}

        def __getitem__(self,item):
                if len(self.hmf.keys())==0:
                        print "initializing functions..."
                        self.create_functions()
                return self.hmf[item]

        def create_functions(self):
                import scipy
                for k in self.data.keys():
                        self.hmf[key] = scipy.interpolate.interp1d(self.data[k][0],self.data[k][1])





