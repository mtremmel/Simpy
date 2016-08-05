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

def get_hmf_data(simpath,log_M_min=8., log_M_max=15.,**kwargs):
        import pynbody
        sim = pynbody.load(simpath)
        sim.properties['sigma8'] = 0.77
        log_M_min += np.log10(sim.properties['h'])
        log_M_max += np.log10(sim.properties['h'])
        mass, sig, phi = pynbody.analysis.halo_mass_function(sim,pspec=pynbody.analysis.hmf.PowerSpectrumCAMBLive,
                                                             log_M_max=log_M_max,log_M_min=log_M_min,**kwargs)
        return np.log10(mass/sim.properties['h']),phi*sim.properties['h']**3, sim.properties['z']

def get_hmf_data_all(log_M_min=8., log_M_max=15.,**kwargs):
        f = open('files.list','r')
        hmf = {'z':[], 'mass':[], 'phi':[]}
        for l in f:
                name = l.strip('\n')
                print name
                lm, phi, z = get_hmf_data(name,log_M_min=log_M_min, log_M_max=log_M_max,**kwargs)
                hmf['z'].append(z)
                hmf['mass'].append(lm)
                hmf['phi'].append(phi)
        hmf['z'] = np.array(hmf['z'])
        return hmf

class HMF(object):
        def __init__(self,log_M_min=8.0, log_M_max=15.0,delta_log_M=0.1,**kwargs):
                self.hmf = get_hmf_data_all(log_M_min=log_M_min, log_M_max=log_M_max,delta_log_M=delta_log_M,**kwargs)
                self.minlm = log_M_min
                self.maxlm = log_M_max
                self.delta = delta_log_M

        def calc_rho(self,logm,z):
                j, = np.where(z>self.hmf['z'])
                if logm > self.maxlm or logm < self.minlm:
                        print "value for logM excedes maximum or is less than minimum value calculated"
                        raise ValueError
                i = int((logm - self.minlm)/self.delta)
                if len(j) ==0:
                        return self.hmf['phi'][-1][i] * self.delta
                else:
                        if j > 0:
                                return (self.hmf['phi'][j-1][i] + \
                                        (self.hmf['phi'][j][i]-self.hmf['phi'][j-1][i])/(self.z[j]-self.z[j-1]) * (z - self.z[j-1]))*self.delta
                        else:
                                return self.hmf['phi'][j]*self.delta

        def get_N_halos(self,dbsim):
                self.z = []
                self.mbins = np.arange(self.minlm, self.maxlm+self.delta, self.delta)
                self.nhalos = np.zeros((len(dbsim.timesteps),len(self.mbins)-1))
                cnt = 0
                for step in dbsim.timesteps:
                        print step.extension
                        Mvir, = step.gather_property('Mvir')
                        N, bins = np.histogram(np.log10(Mvir),bins=self.mbins)
                        self.nhalos[cnt,:] = N
                        self.z.append(step.redshift)
                        cnt += 1
                self.z = np.array(self.z)












