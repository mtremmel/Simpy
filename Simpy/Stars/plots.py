import numpy as np
from .. import plotting, cosmology
from . import get

def SF_MS_hiz(Mstar,z):
	#Tomczak+ 2016
	so = 0.448 + 1.22*z - 0.174*z**2
	logMo = 9.458 + 0.865*z - 0.132*z**2
	gamma = 1.091

	Mo = 10**logMo

	return so - np.log10(1+(Mstar/Mo)**(-1.0*gamma))

def SF_MS_z0_SDSS(Mstar,z):
	#Elbaz+ 2007
	return 8.7 * (Mstar/1e11)**0.77

def dbSFH(halo, return_data=True, **kwargs):
	from .. import dbanalysis
	ssfr,t, mstar = dbanalysis.formation_history.sSFR(halo, return_mstar=True)
	plotting.plt.plot(t,ssfr,**kwargs)
	if return_data is True:
		return ssfr, t, mstar


def plot_sfms_t(mstar,t,z,scatter=0.5, alpha = 0.5, color='grey'):
	from .. import cosmology
	#a,z = cosmology.getScaleFactor(t,sim)
	logSFR = SF_MS(mstar,z)
	up = logSFR - np.log10(mstar) + scatter
	down = logSFR - np.log10(mstar) - scatter

	plotting.plt.fill_between(t,10**up,10**down,color=color,alpha=alpha)



def SFH(sim, color, style, lw=3,bins=50, trange=None, label=None, ret_hist=True, type='time',xlog=False,ylog=False, overplot=False):
	tform = get.tform_fromsnap(sim,'Gyr')
	massform = sim.stars['massform'].in_units('Msol')

	if type =='time':
		if trange is None:
			trange = [0,tform.max()]
		dt = (trange[1] - trange[0])*1e9/float(bins)
		data = np.histogram(tform, range=trange, bins=bins,weights=massform/dt)
		sfr = data[0]
		tedges = data[1]
		plotting.plt.step(tedges[0:-1],sfr,color, linestyle=style, label=label, linewidth=lw, where='post')
		if overplot is False:
			plotting.plt.xlabel(r'Time [Gyr]', fontsize=30)


	if type == 'redshift':
		if trange is None:
			trange = [0,25]
		if xlog is False:
			dz = (trange[1]-trange[0])/float(bins)
			zbins = np.arange(trange[0],trange[1]+dz,dz)
		else:
			dz = (np.log10(trange[1]+1)-np.log10(trange[0]+1))/float(bins)
			zbins = 10**np.arange(np.log10(trange[0]+1),np.log10(trange[1]+1)+dz,dz)
		zbins = zbins[::-1]  #flip array so times will be monotonically increasing
		tedges = np.array([cosmology.getTime(z,sim) for z in zbins])
		dt = (tedges[1:] - tedges[0:-1]) * 1e9
		data = np.histogram(tform, range=trange, bins=tedges,weights=massform)
		sfr = data[0]/dt
		if xlog is False:
			plotting.plt.step(zbins[0:-1],sfr,color, linestyle=style, label=label, linewidth=lw, where='post')
		if xlog is True:
			plotting.plt.step(zbins[0:-1]+1,sfr,color, linestyle=style, label=label, linewidth=lw, where='post')
			if overplot is False:
				plotting.plt.xticks([1,2,3,4,5,6,7,8,9,10],['0','1','2','3','4','5','6','7','8','9'])
		if overplot is False:
			plotting.plt.xlabel('z',fontsize=30)
	if overplot is False:
		plotting.plt.ylabel('SFR [M$_{\odot}$ yr$^{-1}$]', fontsize=30)
		if xlog:
			plotting.plt.xscale('log', base=10)
		if ylog:
			plotting.plt.yscale('log', base=10)
	if label is not None:
		plotting.plt.legend()
	if ret_hist is True:
		if type == 'time':
			return sfr, tedges
		if type == 'redshift':
			return sfr, tedges, zbins
	else:
		return


def cosmicSFH_by_step(dbsim, style, lw=3, range=None, label=None, volume=25**3, retdata=False, plot=True, plotobs=True,xlog=True,ylog=True):
	sfrdata = np.zeros(len(dbsim.timesteps))
	redshift = np.zeros(len(dbsim.timesteps))
	proplist = ['SFR']
	min = []
	max = []
	if range is not None:
		rangespl = range.split(' ')
		if len(rangespl)%3 != 0:
			print "range string given has incorrect syntax, " \
					"must go as property min max property min max, etc"
		proplist.extend(rangespl[0::3])
		min = np.array(rangespl[1::3]).astype(np.float)
		max = np.array(rangespl[2::3]).astype(np.float)
	proplist = tuple(proplist)
	cnt = 0
	for step in dbsim.timesteps:
		print "getting data for redshift", step.redshift
		redshift[cnt] = step.redshift
		data = step.gather_property(*proplist)
		SFR = data[0]
		ok = np.arange(len(SFR))
		if len(data)>1:
			rcnt = 0
			for dd in data[1::]:
				oknew = np.where((dd[ok] > min[rcnt])&(dd[ok] < max[rcnt]))[0]
				ok = ok[oknew]
				rcnt += 1
		SFR = SFR[ok]
		sfrdata[cnt] = np.sum(SFR)/volume
		cnt += 1

	if plot:
		if xlog is True:
			plotting.plt.plot(redshift+1,sfrdata,style,label=label,lw=lw)
			plotting.plt.xscale('log',base=10)
		else:
			plotting.plt.plot(redshift,sfrdata,style,label=label)
		if ylog is True:
			plotting.plt.yscale('log',base=10)
		if plotobs is True:
			cosmicSFH(None,None,xlog=xlog,ylog=ylog,zrange=[redshift.min(),redshift.max()])

	if retdata:
		return redshift,sfrdata




def cosmicSFH(sim, style, lw=3, bins=50, zrange=None, label=None, ret_hist=True,xlog=False,ylog=False,
			  plotdata=True,volume=25**3):
	if zrange is None:
			zrange = [0,25]
	if sim is not None:
		tform = get.tform_fromsnap(sim,'Gyr')
		massform = sim.stars['massform'].in_units('Msol')
		if xlog is True:
			dz = (np.log10(zrange[1]+1)-np.log10(zrange[0]+1))/float(bins)
			zbins = 10**np.arange(np.log10(zrange[0]+1),np.log10(zrange[1]+1)+dz,dz)
			tedges = np.array([cosmology.getTime(z - 1,sim) for z in zbins])
		else:
			dz = (zrange[1] - zrange[0])/float(bins)
			zbins = np.arange(zrange[0],zrange[1]+dz,dz)
			tedges = np.array([cosmology.getTime(z,sim) for z in zbins])

		tsorted = np.argsort(tedges)
		tedges = tedges[tsorted]

		dt = np.abs((tedges[0:-1] - tedges[1:]) * 1e9)
		data = np.histogram(tform, bins=tedges, weights=massform)
		sfr = data[0]/dt
		zbins = zbins[tsorted]
		#sfr = np.append(sfr,sfr[-1])
		plotting.plt.plot(zbins[1:],sfr/volume,style, label=label, linewidth=lw)


	if plotdata is True:
		import colldata as dat

		zfits = np.arange(zrange[0],zrange[1],0.01)
		fitB,sigB = dat.CSFRFit(zfits,type='beh')
		fitH,sigH = dat.CSFRFit(zfits,type='hop')
		if xlog is True:
			plotting.plt.plot(zfits+1,fitB,'k-',label='Behroozi+ 13',lw=2)
			#plotting.plt.plot(zfits+1,fitH,'k--',label='Hopkins+ 06')
			plotting.plt.fill_between(zfits+1,10**(np.log10(fitB)-sigB),10**(np.log10(fitB)+sigB),linewidth=1.5,facecolor='grey',alpha=0.2)
			#plotting.plt.fill_between(zfits+1,10**(np.log10(fitH)-sigH),10**(np.log10(fitH)+sigH),linewidth=1.5,facecolor='grey',alpha=0.2)

			plotting.plt.scatter(np.array(dat.duncan_z)+1, 10**dat.duncan_sf, color='r', marker='D', label='Duncan+ 14', s=100)
			plotting.plt.errorbar(np.array(dat.kist_zhigh)+1, dat.kist_sfr, fmt='o',markersize=10, yerr=[dat.kist_sfrminus, dat.kist_sfrplus], xerr=[dat.kist_zhighminus, dat.kist_zhighplus], ls='None', linewidth=1.5, color='k', label='Kistler+ 13')

		else:
			plotting.plt.plot(zfits,fitB,'k-',label='Behroozi+ 13',lw=2)
			#plotting.plt.plot(zfits,fitH,'k--',label='Hopkins+ 06')
			plotting.plt.fill_between(zfits,10**(np.log10(fitB)-sigB),10**(np.log10(fitB)+sigB),linewidth=1.5,facecolor='grey',alpha=0.2)
			#plotting.plt.fill_between(zfits,10**(np.log10(fitH)-sigH),10**(np.log10(fitH)+sigH),linewidth=1.5,facecolor='grey',alpha=0.2)

			plotting.plt.scatter(dat.duncan_z, 10**dat.duncan_sf, color='r', marker='D', label='Duncan+ 14', s=100)
			plotting.plt.errorbar(dat.kist_zhigh, dat.kist_sfr, fmt='o',markersize=10, yerr=[dat.kist_sfrminus, dat.kist_sfrplus], xerr=[dat.kist_zhighminus, dat.kist_zhighplus], ls='None', linewidth=1.5, color='k', label='Kistler+ 13')


	plotting.plt.ylabel(r'$\rho_{SFR}$ [M$_{\odot}$ yr$^{-1}$ Mpc$^{-1}$]',fontsize=30)
	plotting.plt.legend(loc='lower left',fontsize=35)
	if ylog is True:
		plotting.plt.yscale('log',base=10)
	if xlog is False:
		plotting.plt.xlabel(r'Redshift',fontsize=30)
		plotting.plt.xlim(zrange[0],zrange[1])
	else:
		plotting.plt.xscale('log',base=10)
		plotting.plt.xticks([1,2,3,4,5,6,7,8,9,10,11],['0','1','2','3','4','5','6','7','8','9','10'])
		plotting.plt.xlim(zrange[0]+1,zrange[1]+1)
		plotting.plt.xlabel(r'Redshift',fontsize=30)
