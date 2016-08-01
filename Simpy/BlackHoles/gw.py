from .. import util, cosmology
import numpy as np

class GWemit(object):
#GW emission approx. without spin given by Ajith+ 2008
# (http://journals.aps.org/prd/pdf/10.1103/PhysRevD.77.104017)
	def __init__(self, mass1, mass2, z,omegaM, omegaL, h):
		self.masses = (mass1, mass2)
		self.z = z
		self.cosmology = {'omegaM': omegaM, 'omegaL': omegaL, 'h':h}



	def freq(self, a, b, c):
		M = self.masses[0]+self.masses[1]
		n = self.masses[0]*self.masses[1]/M**2

		return (a*n**2 + b*n + c)/(np.pi*M) * util.c.in_units('cm s**-1')**3/util.G.in_units('Msol**-1 cm**3 s**-2')

	def freq_merge(self):
		a = 2.9740e-1
		b = 4.4810e-2
		c = 9.5560e-2
		return self.freq(a,b,c)

	def freq_ring(self):
		a = 5.94110e-1
		b = 8.97940e-2
		c = 1.91110e-1
		return self.freq(a,b,c)

	def freq_cut(self):
		a = 8.49450e-1
		b = 1.28480e-1
		c = 2.72990e-1
		return self.freq(a,b,c)

	def sigfreq(self):
		a = 5.08010e-1
		b = 7.75150e-2
		c = 2.23690e-1
		return self.freq(a,b,c)

	def ampGW(self, f):
		M = self.masses[0]+self.masses[1]
		n = self.masses[0]*self.masses[1]/M**2
		dL = cosmology.lum_distance(self.z,self.cosmolggy['omegaM'], self.cosmolgoy['omegaL'],self.cosmology['h'])

		#dL *= 3.086e+24
		C = (util.G.in_units('Mpc**3 Msol**-1 s**-2')**(5./6.)/util.c.in_untits('Mpc s**-1')**(3./2.)) * \
			M**(5./6.) * self.freq_merge()**(-7./6.) * (5.*n/24.)**0.5 / (dL*np.pi**(2./3.))

		fm = self.freq_merge()
		fr = self.freq_ring()
		fc = self.freq_cut()
		sig = self.sigfreq()

		def LL(sig,f,fr):
			return (1/2.*np.pi)*sig/((f-fr)**2 + sig**2/4.)

		if f < fm:
			return C * (f/fm)**(-7./6.)
		if f >= fm and f < fr:
			return C * (f/fm)**(-2./3.)
		if f >= fr and f < fc:
			return C * np.pi*sig/2. * (fr/fm)**(-2./3.) * LL(sig,f,fr)
		if f >=fc:
			print "ERROR frequency > f_cut"
			raise RuntimeError

	def ampGW_merger(self):
		return self.ampGW(self.freq_merge())

	def ampGW_ring(self):
		return self.ampGW(self.freq_ring())


class SalcidoGW(object):
	# GW calculations based on Salcido+ 2015, Flanagen + Hughes 1998, and others
	def __init__(self, mass1, mass2, z, omegaM, omegaL, h):
		self.masses = (mass1, mass2)
		self.z = z
		self.cosmology = {'omegaM': omegaM, 'omegaL': omegaL, 'h':h}

	def freq_qnr(self, a=0.95):
		Mtot = self.masses[0] + self.masses[1]
		Mtot *= util.M_sun_g
		return util.c**3 * (1 - 0.63*(1-a)**(3./10.)) / (2 * np.pi * util.G * Mtot)

	def freq_merger(self):
		Mtot = self.masses[0] + self.masses[1]
		return util.c**3 * 0.02 / (util.G * Mtot)

	def dEdf(self, eps=0.1, a=0.95):
		fqnr = self.freq_qnr(a)
		fm = self.freq_merger()
		M1 = self.masses[0]
		M2 = self.masses[1]
		M1 *= util.M_sun_g
		M2 *= util.M_sun_g
		F = (4.*M1*M2)**2/(M1+M2)**4
		return util.c**2*(M1+M2)*F*eps/(fqnr-fm)

	def strain(self, eps=0.1, a=0.95):
		dL = cosmology.lum_distance(self.z,self.cosmology['omegaM'], self.cosmology['omegaL'], self.cosmology['h'])
		dL *= 3.086e+24
		Ef = self.dEdf(eps, a)
		return np.sqrt(2.*util.G/util.c**3) * ((1+self.z)/(np.pi*dL)) * np.sqrt(Ef)

def eLisaLimitLPF(fobs, config=1):
	def sacc(fobs):
		return 1.551245e-29*(1 + (4.48833e-4/fobs)**2 + (fobs/0.0636734)**4) / (2*np.pi*fobs)**4

	somn = 8.65e-24

	if config == 1:
		ssn = 4.39e-23
	if config == 2:
		ssn =  3.90e-23
	if config == 5:
		ssn = 5.86e-23

	L = config * 1.0e6

	def sn(fobs):
		return 20./3. * (4.*sacc(fobs)+ssn+somn)/L**2 * (1+(fobs/(0.41*(util.c/(1.0e2*2.0*L))))**2)

	return np.sqrt(sn(fobs)*fobs)

#eLisa limits from Salcido+2015, taken from Amaro-Seoane+ 2013
def eLisaLimit_old(fobs):
	def sacc(fobs):
		return 1.37e-32 * (1+1e-4/fobs)*(fobs)**-4
	ssn = 5.25e-23
	somn = 6.28e-23
	L = 1.0e9
	def sn(fobs):
		return 20./3. * (4.*sacc(fobs)+ssn+somn)/L**2 * (1+(fobs/(0.41*(util.c/(1.0e2*2.0*L))))**2)

	return np.sqrt(sn(fobs)*fobs)