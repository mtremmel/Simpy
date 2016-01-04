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
	alpha = 1.74
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


def SMHM(sim, step, style, fitstyle=['k-','k--'], fit=['Mos', 'Krav'], minmass=None, maxmass=None, markersize=10,label=None, correct=True):
	import halo_db as db
	step = db.get_timestep(sim+'/%'+str(step))
	Mvir, Mstar = step.gather_property('Mvir', 'Mstar')

	if correct is True:
		Mstar *= 0.6
		Mvir /= 0.8

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

	plotting.plt.plot(np.log10(Mvir), np.log10(Mstar/Mvir), style, markersize=markersize, label=label)
	plotting.plt.ylabel(r'M$_{*,central}$/M$_{vir}$',fontsize=30)
	plotting.plt.xlabel(r'M$_{vir}$ [M$_{\odot}$]',fontsize=30)
	#plotting.plt.xscale('log',base=10)
	#plotting.plt.yscale('log',base=10)
	plotting.plt.legend(loc='lower right',fontsize=20)





def SMHM_old(simname, step, style, fitstyle='k-',skipgrp=None, maxgrp = None, minmass = None, plotratio=True, plotfit='mos', ploterr = True, nodata=False, lw=3, correct = True, marksize=10, label=None, overplot=False, sats='Remove'):
	amigastat, redshift, type = getstats(simname, step)
	print "z = ", redshift
	if minmass is None:
		minmass = 1e8
	minm = np.log10(minmass)
	maxm = np.log10(amigastat['Mvir(M_sol)'].max())

	if nodata is False:
		if skipgrp is not None:
			skipgrp = np.array(skipgrp)
			ok, = np.where(np.in1d(amigastat['Grp'], skipgrp) == False)
			util.cutdict(amigastat,ok)
		if maxgrp is not None:
			ok, = np.where(amigastat['Grp']<maxgrp)
			util.cutdict(amigastat,ok)
		if sats == 'Remove':
			if type == 'rockstar':
				ok, = np.where(amigastat['Satellite?'] == -1)
			if type == 'amiga':
				ok, = np.where(amigastat['Satellite?'] == 'no')
			util.cutdict(amigastat,ok)
		if sats == 'Only':
			if type == 'rockstar':
				ok, = np.where(amigastat['Satellite?'] != -1)
			if type == 'amiga':
				ok, = np.where(amigastat['Satellite?'] == 'yes')
		ok, = np.where(amigastat['Mvir(M_sol)']>minmass)
		util.cutdict(amigastat,ok)

		xdata = amigastat['Mvir(M_sol)']
		ydata = amigastat['StarMass(M_sol)']
		if correct is True:
			ydata *= 0.6
			xdata /= 0.8
		if plotratio is True:
			ydata = ydata/xdata
		ydata = np.log10(ydata)
		xdata = np.log10(xdata)
		plotting.plt.plot(xdata, ydata, style, markersize=marksize, label=label)

	if plotfit is not None:
		zlabel = '%0.4f'%(redshift)
		lmhaloline = np.arange(minm-1.0,maxm+0.5,0.01)
		if plotfit != 'mos' and plotfit != 'beh':
			print "WARNING! Input string for plotfit type is not recognized... defaulting to Moster+ 13"
			plotfit = 'mos'
		if plotfit == 'mos':
			ratiofit = moster13(lmhaloline, redshift)
			fitlabel = "Moster+ 13, z = "+zlabel
			if ploterr is True:
				sigma = errmoster13(lmhaloline, redshift)
		if plotfit == 'beh':
			ratiofit = behroozi13(lmhaloline, redshift)
			fitlabel = "Behroozi+ 13, z = "+zlabel
			if ploterr is True:
				print "WARNING ploterr is set to True but behroozi was chosen. Errors not implemented in this yet."
		plotting.plt.plot(lmhaloline, np.log10(ratiofit), fitstyle, linewidth=lw, label=fitlabel)
		if ploterr is True:
			ok, = np.where(np.isnan(np.log10(ratiofit-sigma))==False)
			plotting.plt.fill_between(lmhaloline[ok],np.log10(ratiofit[ok]-sigma[ok]),np.log10(ratiofit[ok]+sigma[ok]),facecolor='grey',alpha=0.4)


	plotting.plt.legend(loc='lower right', fontsize=25)
	plotting.plt.ylabel(r'M$_{*}$/M$_{vir}$',fontsize=30)
	plotting.plt.xlabel(r'M$_{vir}$',fontsize=30)














