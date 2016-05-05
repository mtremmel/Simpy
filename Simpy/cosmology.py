from pynbody.analysis import pkdgrav_cosmo as cosmo
import pynbody.analysis.cosmology as other_cosmo
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

def HoverH0(z, sim):
        return np.sqrt(sim.properties['omegaM0']*(1+z)**3+sim.properties['omegaL0'])

def comoving_dist(z, sim):
        dH = util.c*1e5/(sim.properties['h']*100)

        def func(z):
                return HoverH0(z,sim)

        def get_d(x):
                return dH * integrate.quad(func,0,x)

        if isinstance(z, np.ndarray) or isinstance(z, list):
                return np.SimArray((map(get_d, z)),'Mpc')
        else:
                return get_d(z)

def lum_distance(z, sim):
        dc = comoving_dist(z,sim)
        if isinstance(z, np.ndarray) or isinstance(z, list):
                return np.SimArray(dc*(1+z),'Mpc')
        else:
                return dc*(1+z)



