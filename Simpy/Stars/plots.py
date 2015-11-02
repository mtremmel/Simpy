import numpy as np
from .. import plotting, cosmology
from . import get

def SFH(sim, style, lw=3,bins=50, trange=None, label=None, ret_hist=True, type='time',xlog=False,ylog=False, overplot=False):
	tform = get.tform_fromsnap(sim,'Gyr')
	massform = sim.stars['massform'].in_units('Msol')

	if type =='time':
		if trange is None:
			trange = [0,tform.max()]
		dt = (trange[1] - trange[0])*1e9/float(bins)
		data = np.histogram(tform, range=trange, bins=bins,weights=massform/dt)
		sfr = data[0]
		tedges = data[1]
		plotting.plt.step(tedges[0:-1],sfr,style, label=label, linewidth=lw, where='post')
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
			plotting.plt.step(zbins[0:-1],sfr,style, label=label, linewidth=lw, where='post')
		if xlog is True:
			plotting.plt.step(zbins[0:-1]+1,sfr,style, label=label, linewidth=lw, where='post')
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

def cosmicSFH(sim, style, lw=3, bins=50, zrange=None, label=None, ret_hist=True,xlog=False,ylog=False, plotdata=True,volume=25**3, overplot=False):
	tform = get.tform_fromsnap(sim,'Gyr')
	massform = sim.stars['massform'].in_units('Msol')
	if zrange is None:
		zrange = [0,25]
	if xlog is False:
		dz = (zrange[1]-zrange[0])/float(bins)
		zbins = np.arange(zrange[0],zrange[1]+dz,dz)
	else:
		dz = (np.log10(zrange[1]+1)-np.log10(zrange[0]+1))/float(bins)
		zbins = 10**np.arange(np.log10(zrange[0]+1),np.log10(zrange[1]+1)+dz,dz)
	tedges = np.array([cosmology.getTime(z,sim) for z in zbins])
	tsorted = np.argsort(tedges)
	dt = (tedges[1:] - tedges[0:-1]) * 1e9
	tedges = tedges[tsorted]
	data = np.histogram(tform, bins=tedges, weights=massform)
	sfr = data[0]/dt
	zbins = zbins[tsorted]

	if xlog is False:
		plotting.plt.step(zbins[0:-1],sfr/volume,style, label=label, linewidth=lw, where='post')
	else:
		plotting.plt.step(zbins[0:-1]+1,sfr/volume,style, label=label, linewidth=lw, where='post')
		plotting.plt.xticks([1,2,3,4,5,6,7,8,9,10],['0','1','2','3','4','5','6','7','8','9'])
		plotting.plt.xscale('log',base=10)
	if ylog is True:
		plotting.plt.yscale('log',base=10)

	if plotdata is True:
		import colldata as dat

		zfits = np.arange(1,11,0.01)
		fitB,sigB = dat.CSFRFit(zfits,type='beh')
		fitH,sigH = dat.CSFRFit(zfits,type='hop')
		plotting.plt.plot(zfits,fitB,'k-',label='Behroozi+ 13')
		plotting.plt.plot(zfits,fitH,'k--',label='Hopkins+ 06')
		plotting.plt.fill_between(zfits,10**(np.log10(fitB)-sigB),10**(np.log10(fitB)+sigB),linewidth=1.5,facecolor='grey',alpha=0.2)
		plotting.plt.fill_between(zfits,10**(np.log10(fitH)-sigH),10**(np.log10(fitH)+sigH),linewidth=1.5,facecolor='grey',alpha=0.2)

		plotting.plt.scatter(dat.duncan_z, 10**dat.duncan_sf, color='k', marker='D', label='Duncan+ 14', s=40)
		plotting.plt.errorbar(dat.kist_zhigh, dat.kist_sfr, fmt='o',yerr=[dat.kist_sfrminus, dat.kist_sfrplus], xerr=[dat.kist_zhighminus, dat.kist_zhighplus], ls='None', linewidth=1.5, color='k', label='Kistler+ 13')

	if ylog is True:
		plotting.plt.yscale('log',base=10)
	if not overplot:
		plotting.plt.xticks([0,1,2,3,4,5,6,7,8,9,10],['0','1','2','3','4','5','6','7','8','9','10'])
	if overplot is False and plotdata is True:
		plotting.plt.ylabel(r'$\rho_{SFR}$ [M$_{\odot}$ yr$^{-1}$ Mpc$^{-1}$]',fontsize=30)
		plotting.plt.xlabel(r'Redshift',fontsize=30)
	plotting.plt.legend(loc='lower left',fontsize=20)