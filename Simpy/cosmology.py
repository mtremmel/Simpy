from pynbody.analysis import pkdgrav_cosmo as cosmo
import pynbody.analysis.cosmology as other_cosmo
import numpy as np
from scipy import optimize as opt

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
