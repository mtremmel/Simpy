import numpy as np
from pynbody.analysis import pkdgrav_cosmo as cosmo

def getTime(z,sim):
        c = cosmo.Cosmology(sim=sim)
        return 13.7*c.Exp2Time(1.0 / (1+z))/c.Exp2Time(1)

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
                        print "ERROR did not converge", times[tt],tt
                        redshift[tt] = -1
        scaleFac = 1./(1+redshift)
        return scaleFac, redshift
