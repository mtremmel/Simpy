import numpy as np
from .. import readcol
from .. import util
from .. import plotting
import pynbody
import os
import gc


#Moster relations
def moster13allvar(logM,z,M10,M11,N10,N11,b10,b11,g10,g11):
	r = z/(z+1)
	logM1 = M10 + M11 * r
	N = N10 + N11 * r
	b = b10 + b11 * r
	g = g10 + g11 * r
	x = logM - logM1
	ratio = 2 * N / ( 10**(-b * x) + 10**(g * x) )
	return ratio

def moster13(logM, z):
    M10 = 11.590
    M11 =  1.195
    N10 =  0.0351
    N11 = -0.0247
    b10 =  1.376
    b11 = -0.826
    g10 =  0.608
    g11 =  0.329
    return moster13allvar(logM,z,M10,M11,N10,N11,b10,b11,g10,g11)

def errmoster13(logM,z):
	M10 = 11.590
	M11 =  1.195
	N10 =  0.0351
	N11 = -0.0247
	b10 =  1.376
	b11 = -0.826
	g10 =  0.608
	g11 =  0.329
	sigM10 = 0.236
	sigM11 = 0.353
	sigN10 = 0.0058
	sigN11 = 0.0069
	sigb10 = 0.153
	sigb11 = 0.225
	sigg10 = 0.059
	sigg11 = 0.173
	sigvar = [sigM10,sigM11,sigN10,sigN11,sigb10,sigb11,sigg10,sigg11]
	sigma = np.zeros(len(logM))
	for i in range(len(logM)):
		point = [logM[i],z,M10,M11,N10,N11,b10,b11,g10,g11]
		for j in range(8):
			sigma[i] += util.partial_derivative(moster13allvar,var=j+2,point=point)**2 * sigvar[j]**2
	sigma = np.sqrt(sigma)
	return sigma

#Behroozi data
def behroozi13(logM,a):
	z = a**-1 - 1
	v = np.exp(-4.*a**2)
	le = -1.777 + (-0.006*(a-1)*v) - 0.119*(a-1)
	lM1 = 11.514 + (-1.793*(a-1)-0.251*z)*v
	A = -1.412+0.731*(a-1)*v
	d = 3.508+(2.608*(a-1)-0.043*z)*v
	g = 0.316+(1.319*(a-1)+0.279*z)*v
	def f(x):
		return -1.0*np.log10(10**(A*x)+1) + d*((np.log10(1+np.exp(x)))**g)/(1+np.exp(10**(-x)))
	lMstar = le +lM1 + f(logM-lM1) - f(0)
	ratio = 10**(lMstar - logM)
	return ratio

#Kravstov data
def kravstov14(logM, z):
	loge = -1.685
	logM1 = 11.39
	alpha = -1.74
	delta = 4.335
	gamma = 0.531

	def f(x):
		return -1.0 * np.log10(10.0**(alpha * x) + 1.0) + (delta * np.log10(1.0 + np.exp(x))**gamma) / (1.0 + np.exp(10.**(-1.0 * x)))

	logMstar = loge + logM1 + f(logM - logM1) - f(0)
	return 10**(logMstar - logM)


def getstats(simname, step):
	print "getting halo stats..."
	if os.path.exists(simname + '.' + step + '.amiga.stat'):
		amigastat = readcol.readcol(simname + '.' + step + '.amiga.stat', asdict=True)
		type = 'amiga'
	else:
			print "amiga file failed, looking for rockstar"
			if os.path.exists(simname + '.' + step + '.rockstar.stat'):
				amigastat = readcol.readcol(simname + '.' + step + '.rockstar.stat', asdict=True)
				type = 'rockstar'
			else:
				print "ERROR cannot find recognized stat file (amiga.stat or rockstar.stat)"
	s = pynbody.load(simname+'.'+step)
	a = s.properties['a']
	del(s)
	gc.collect()
	return amigastat, a**-1 -1, type


def SMHM_db(sim, step, style, fitstyle=['k-','k--'], fit=['Mos', 'Krav'], minmass=None, maxmass=None, markersize=5,label=None, correct=True):
	import halo_db as db
	step = db.get_timestep(sim+'/%'+str(step))
	Mvir, Mstar = step.gather_property('Mvir', 'Mstar')

	if correct is True:
		Mstar *= 0.6
		Mvir /= 0.8

	plotting.plt.plot(np.log10(Mvir), np.log10(Mstar/Mvir), style, markersize=markersize, label=label)

	if minmass is None:
		minmass = Mvir.min()/2.
	if maxmass is None:
		maxmass = Mvir.max()*2.

	cnt = 0
	logmv_fit = np.arange(np.log10(minmass),np.log10(maxmass),0.1)
	if fit is not None:
		for ff in fit:
			if ff not in ['Mos','Beh', 'Krav', 'Moster', 'Behroozi', 'Kravtsov']:
				print "fit request ", ff, " not understood... should be in list", \
					['Mos','Beh', 'Krav', 'Moster', 'Behroozi', 'Kravtsov']
				continue
			print "found fit relation for", ff

			if ff in ['Mos', 'Moster']:
				fitfunc = moster13
				flabel = 'Moster+ 13'
			if ff in ['Beh', 'Behroozi']:
				fitfunc = behroozi13
				flabel = 'Behroozi+ 13'
			if ff in ['Krav', 'Kravtsov']:
				fitfunc = kravstov14
				flabel = 'Kravtsov+ 14'

			ratio_fit = fitfunc(logmv_fit, step.redshift)
			plotting.plt.plot(logmv_fit, np.log10(ratio_fit), fitstyle[cnt], label=flabel)

			cnt += 1

	#plotting.plt.plot(np.log10(Mvir), np.log10(Mstar/Mvir), style, markersize=markersize, label=label)
	plotting.plt.ylabel(r'M$_{*,central}$/M$_{vir}$',fontsize=30)
	plotting.plt.xlabel(r'M$_{vir}$ [M$_{\odot}$]',fontsize=30)
	#plotting.plt.xscale('log',base=10)
	#plotting.plt.yscale('log',base=10)
	plotting.plt.legend(loc='lower right',fontsize=20)





def SMHM(sim, step, style, fitstyle=['k-','k--'], fit=['Mos', 'Krav'], minmass=None, maxmass=None, markersize=5,label=None, correct=True, usedb=False, remove_sats=True, only_sats=False):
	if usedb is True:
		import halo_db as db
		dbstep = db.get_timestep(sim+'/%'+step)
		if remove_sats or only_sats:
			Mvir, Mstar, sub = dbstep.gather_property('Mvir', 'Mstar', 'Sub')
			if len(Mvir)==0:
				print "Warning Sub halos not implemented yet in your database! " \
					  "please turn off sats options or re-run with usedb False"
				return
		else:
			Mvir, Mstar = dbstep.gather_property('Mvir', 'Mstar')
		redshift = dbstep.redshift
	else:
		amigastat, redshift, type = getstats(sim, step)
		Mvir = amigastat['Mvir(M_sol)']
		Mstar = amigastat['StarMass(M_sol)']
		if remove_sats:
			if type == 'rockstar':
				ok, = np.where(amigastat['Satellite?'] == -1)
			if type == 'amiga':
				ok, = np.where(amigastat['Satellite?'] == 'no')
			Mvir = Mvir[ok]
			Mstar = Mstar[ok]
		if only_sats:
			if type == 'rockstar':
				ok, = np.where(amigastat['Satellite?'] != -1)
			if type == 'amiga':
				ok, = np.where(amigastat['Satellite?'] == 'yes')
		if only_sats or remove_sats:
			Mvir = Mvir[ok]
			Mstar = Mstar[ok]

	if correct is True:
		Mstar *= 0.6
		Mvir /= 0.8

	plotting.plt.plot(Mvir, Mstar/Mvir, style, markersize=markersize, label=label)

	if minmass is None:
		minmass = Mvir.min()/2.
	if maxmass is None:
		maxmass = Mvir.max()*2.

	logmv_fit = np.arange(np.log10(minmass),np.log10(maxmass),0.1)
	if fit is not None:
		for ff in fit:
			if ff not in ['Mos','Beh', 'Krav', 'Moster', 'Behroozi', 'Kravtsov']:
				print "fit request ", ff, " not understood... should be in list", \
					['Mos','Beh', 'Krav', 'Moster', 'Behroozi', 'Kravtsov']
				continue
			print "found fit relation for", ff

			if ff in ['Mos', 'Moster']:
				fitfunc = moster13
				flabel = 'Moster+ 13'
			if ff in ['Beh', 'Behroozi']:
				fitfunc = behroozi13
				flabel = 'Behroozi+ 13'
			if ff in ['Krav', 'Kravtsov']:
				fitfunc = kravstov14
				flabel = 'Kravtsov+ 14'

			ratio_fit = fitfunc(logmv_fit, step.redshift)
			plotting.plt.plot(10**logmv_fit, ratio_fit, fitstyle[cnt], label=flabel)

			cnt += 1

	plotting.plt.ylabel(r'M$_{*,central}$/M$_{vir}$',fontsize=30)
	plotting.plt.xlabel(r'M$_{vir}$ [M$_{\odot}$]',fontsize=30)
	plotting.plt.legend(loc='lower right',fontsize=20)
	plotting.plt.yscale('log',base=10)
	plotting.plt.xscale('log',base=10)
















