import numpy as np
from .. import plotting, cosmology
from . import get

def plt_SFH(sim, color='b', linestyle='solid', lw=3,bins=50, trange=None, label=None, ret_hist=True, type='time',xlog=False,ylog=False):
	tform = get.tform_fromsnap(sim,'Gyr')
	massform = sim.stars['massform'].in_units('Msol')

	if type =='time':
		if trange is None:
			trange = [0,tform.max()]
		dt = (trange[1] - trange[0])*1e9/float(bins)
		data = plotting.plt.hist(tform, range=trange, bins=bins, linestyle=linestyle, color=color,label=label,histtype='step',weights=massform/dt, linewidth=lw)
		sfr = data[0]
		tedges = data[1]
		plotting.plt.xlabel(r'Time [Gyr]', fontsize=30)


	if type == 'redshift':
		if trange is None:
			trange = [0,25]
		dz = (trange[1]-trange[0])/float(bins)
		zbins = np.arange(trange[0],trange[1]+dt,dt)
		tedges = np.array([cosmology.getTime(z,sim) for z in zbins])
		dt = (tedges[1:] - tedges[0:-1]) * 1e9
		data = np.histogram(tform, range=trange, bins=tedges,weights=massform)
		sfr = data[0]/dt
		plotting.plt.step(zbins[0:-1],sfr,linestyle=linestyle, color=color, label=label, linewidth=lw, where='post')

	plotting.plt.ylabel('SFR [M$_{\odot}$ yr$^{-1}$]')
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