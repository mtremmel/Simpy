import numpy as np
from .. import plotting, cosmology
from . import get

def plt_SFH(sim, color='b', linestyle='solid', lw=3,bins=50, trange=None, label=None, ret_hist=True, type='time',xlog=False,ylog=False):
	tform = get.tform_fromsnap(sim,'Gyr')
	massform = sim.stars['massform'].in_units('Msol')

	if type =='time':
		if trange is None:
			trange = [0,tform.max()]
		dt = (trange[1] - trange[0])*1e9/bins
		data = plotting.plt.hist(tform, range=trange, bins=bins, linestyle=linestyle, color=color,label=label,histtype='step',weights=massform/dt, linewidth=lw)
		sfr = data[0]
		tedges = data[1]
		plotting.plt.xlabel(r'Time [Gyr]', fontsize=30)


	if type == 'redshift':
		aform, zform = cosmology.getScaleFactor(tform,sim)
		if trange is None:
			trange = [0,zform.max()]
		n, zbins = np.histogram(zform, range=trange, bins=bins)
		tedges = np.array([cosmology.getTime(z, sim) for z in zbins])
		dt = (tedges[1:] - tedges[0:-1]) * 1e9
		sfr = n*massform/dt
		if xlog is False:
			plotting.plt.step(zbins[0:-1], sfr, linestyle=linestyle, color=color, label=label, linewidth=lw, where='post')
			plotting.plt.xlabel('z',fontsize=30)
		else:
			plotting.plt.step(zbins[0:-1]+1, sfr, linestyle=linestyle, color=color, label=label, linewidth=lw, where='post')
			plotting.plt.xlabel('z + 1',fontsize=30)


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