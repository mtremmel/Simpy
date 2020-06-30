from pynbody.analysis import pkdgrav_cosmo as cosmo
import pynbody.analysis.cosmology as other_cosmo
import pynbody.array as arr
import numpy as np
from scipy import optimize as opt
from scipy import integrate
from scipy import interpolate
from . import util
import math

_interp_points = 1000

def _a_dot(a, h0, om_m, om_l):
    om_k = 1.0 - om_m - om_l
    return h0 * a * np.sqrt(om_m * (a ** -3) + om_k * (a ** -2) + om_l)

def _a_dot_recip(*args):
    return 1. / _a_dot(*args)

def tdyne(time, sim=None, OmegaM=0.3086, Lambda=0.6914, Omegab=0.456, h0=0.67769, overden=200):
        a, z = getScaleFactor(arr.SimArray([time], 'Gyr'), s=sim, OmegaM=OmegaM, h0=h0, Omegab=Omegab, Lambda=Lambda)
        rho = overden * rho_crit(z, OmegaM=OmegaM, h0=h0, Omegab=Omegab, Lambda=Lambda, unit='Msol kpc**-3')
        return np.sqrt(3*np.pi/(16*util.G.in_units('Msol**-1 kpc**3 Gyr**-2')*rho))

def ntdyne(tref, t, **kwargs):
        def func(x):
                return 1 / tdyne(x, **kwargs)
        if isinstance(t, np.ndarray) or isinstance(t, list):
                results = []
                for tt in t:
                        results.append(integrate.quad(func, tref, tt)[0])
        else:
                results = integrate.quad(func, tref, t)[0]
        return np.array(results)

def interp_ntdyne(tref=None, times = None, **kwargs):
        if tref is None and times is None and len(list(kwargs.keys()))==0:
                import pickle
                import os
                default_file = os.path.join(os.path.dirname(__file__), "default_tdyne_data.pkl")
                f = open(default_file, 'rb')
                tdyne_interp_data = pickle.load(f)
                f.close()
                td_ex = tdyne_interp_data['tdyne']
                times = tdyne_interp_data['times']
                print((tdyne_interp_data['notes']))
        else:
                if not tref:
                        tref = 0.2
                if not times:
                        times = np.arange(0.3, 13.9, 0.1)
                if not sim:
                        td_ex = ntdyne(tref, times, **kwargs)
                else:
                        td_ex = ntdyne(tref, times, sim=sim)
        return interpolate.interp1d(times, td_ex)

def rho_crit(z, sim=None, OmegaM=0.3086, Lambda=0.6914, Omegab=0.456, h0=0.67769, unit=None):
        if z is None:
                z = sim.properties['z']

        if unit is None:
                try:
                        unit = sim.dm["mass"].units / sim.dm["pos"].units ** 3
                except arr.units.UnitsException:
                        unit = arr.units.NoUnit()

        if hasattr(unit, "_no_unit"):
                unit = units.Unit("Msol kpc^-3 a^-3")

        if sim:
                omM = sim.properties['omegaM0']
                omL = sim.properties['omegaL0']
                h0 = sim.properties['h']
        else:
                omM = OmegaM
                omL = Lambda
                h0 = h0
        a = 1.0 / (1.0 + z)
        H_z = _a_dot(a, h0, omM, omL) / a
        H_z = arr.units.Unit("100 km s^-1 Mpc^-1") * H_z

        rho_crit = (3 * H_z ** 2) / (8 * math.pi * arr.units.G)

        return arr.SimArray(rho_crit.in_units(unit), unit)

def getTime(z,sim=None, OmegaM=0.3086, Lambda=0.6914, Omegab=0.456, h0=0.67769, unit='Gyr'):
        c = cosmo.Cosmology(sim=sim, Ob = Omegab, Om=OmegaM, L = Lambda)
        import scipy
        import scipy.integrate
        from scipy.interpolate import interp1d

        if z is None:
                z = sim.properties['z']

        if sim is not None:
                h0 = sim.properties['h']
                omM = sim.properties['omegaM0']
                omL = sim.properties['omegaL0']
        else:
                omM = OmegaM
                omL = Lambda

        conv = arr.units.Unit("0.01 s Mpc km^-1").ratio(unit)

        def get_age(x):
                x = 1.0 / (1.0 + x)
                return scipy.integrate.quad(_a_dot_recip, 0, x, (h0, omM, omL))[0] * conv

        if isinstance(z, np.ndarray) or isinstance(z, list):
                if len(z) > _interp_points:
                        a_vals = np.logspace(-3, 0, _interp_points)
                        z_vals = 1. / a_vals - 1.
                        log_age_vals = np.log(getTime(z_vals, sim))
                        interp = interp1d(np.log(a_vals), log_age_vals, bounds_error=False)
                        log_a_input = np.log(1. / (1. + z))
                        results = np.exp(interp(log_a_input))
                else:
                        results = np.array(list(map(get_age, z)))
                results = results.view(SimArray)
                results.units = unit
        else:
                results = get_age(z)
        return results
        #return other_cosmo.age(sim,z=z)*c.Exp2Time(1.0 / (1+z))/c.Exp2Time(1)

def getScaleFactor(times,s=None, OmegaM=0.3086, Lambda=0.6914, Omegab=0.456, h0=0.67769, verbose=False):
        redshift = np.zeros(np.size(times))
        ntimes = np.size(times)
        for tt in range(np.size(times)):
                if tt%100==0 and verbose is True:
                        print((tt/np.float(np.size(times)) * 100, '% done'))
                def func(z):
                        return getTime(z, s, OmegaM=OmegaM, Lambda=Lambda, Omegab=Omegab, h0=h0) - times.in_units('Gyr')[tt]
                try: redshift[tt] = opt.newton(func, 0)
                except:
                        try:
                                print("trying again")
                                redshift[tt] = opt.newton(func, 20)
                        except:
                                print(("ERROR did not converge", times[tt], tt))
                                redshift[tt] = -1
        scaleFac = 1./(1+redshift)
        return scaleFac, redshift

def HoverH0(z, omegaM, omegaL):
        return np.sqrt(omegaM*(1+z)**3 + omegaL)

def comoving_dist(z, omegaM, omegaL, h):
        dH = util.c/(1.0e5*h*100.)

        def func(z):
                return 1./HoverH0(z, omegaM, omegaL)

        def get_d(x):
                return dH * integrate.quad(func, 0, x)[0]

        if isinstance(z, np.ndarray) or isinstance(z, list):
                return np.array(list(map(get_d, z)))
        else:
                return get_d(z)

def lum_distance(z, omegaM, omegaL, h):
        dc = comoving_dist(z, omegaM, omegaL, h)
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
        mass, sig, phi = pynbody.analysis.halo_mass_function(sim, pspec=pynbody.analysis.hmf.PowerSpectrumCAMBLive,
                                                             log_M_max=log_M_max, log_M_min=log_M_min, **kwargs)
        return np.log10(mass/sim.properties['h']), phi*sim.properties['h']**3, sim.properties['z']

def get_hmf_data_all(log_M_min=8., log_M_max=15.,**kwargs):
        f = open('files.list', 'r')
        hmf = {'z':[], 'mass':[], 'phi':[]}
        for l in f:
                name = l.strip('\n')
                print(name)
                lm, phi, z = get_hmf_data(name, log_M_min=log_M_min, log_M_max=log_M_max, **kwargs)
                hmf['z'].append(z)
                hmf['mass'].append(lm)
                hmf['phi'].append(phi)
        hmf['z'] = np.array(hmf['z'])
        return hmf

class HMF(object):
        def __init__(self,log_M_min=8.0, log_M_max=15.0,delta_log_M=0.1,**kwargs):
                self.hmf = get_hmf_data_all(log_M_min=log_M_min, log_M_max=log_M_max, delta_log_M=delta_log_M, **kwargs)
                self.minlm = log_M_min
                self.maxlm = log_M_max
                self.delta = delta_log_M

        def calc_rho(self, logm, z):
                j, = np.where(z>self.hmf['z'])
                if logm > self.maxlm or logm < self.minlm:
                        print("value for logM excedes maximum or is less than minimum value calculated")
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

        def get_N_halos(self,dbsim, check_contam=True):
                self.z = []
                self.mbins = np.arange(self.minlm, self.maxlm+self.delta, self.delta)
                self.nhalos = np.zeros((len(dbsim.timesteps), len(self.mbins)-1))
                cnt = 0
                for step in dbsim.timesteps:
                        print((step.extension))
                        if not check_contam:
                                Mvir, = step.gather_property('Mvir')
                        else:
                                Mvir, contam =  step.gather_property('Mvir', 'contamination_fraction')
                                Mvir = Mvir[(contam<0.1)]
                        N, bins = np.histogram(np.log10(Mvir), bins=self.mbins)
                        self.nhalos[cnt,:] = N
                        self.z.append(step.redshift)
                        cnt += 1
                self.z = np.array(self.z)












