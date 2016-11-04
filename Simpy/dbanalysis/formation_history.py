import tangos as db
import tangos.relation_finding
import numpy as np


def get_merger_ratios(merger_list, property):
	ratio = []
	for mm in merger_list:
		try:
			ratio.append(min(mm[0][property],mm[1][property])/max(mm[0][property],mm[1][property]))
		except KeyError:
			break
	return np.array(ratio)


def sSFR(halo):
	SFH = halo.calculate("reassemble(SFR_histogram)")
	Mstar, time = halo.reverse_property_cascade("Mstar","t()")
	t_sfh = np.arange((len(SFH)))*0.01
	Mstar_interp = np.interp(t_sfh,time[::-1],Mstar[::-1])
	sSFR = SFH/Mstar_interp

	return sSFR, t_sfh


def plot_merger_times(halo, y_range, sim = None, ratio_property=None, convert_to_time=True, ratio_range=None, alpha=None, **kwargs):
	import tangos.examples.mergers as mergers
	from .. import plotting
	z, ratio, merger_list = mergers.get_mergers_of_major_progenitor(halo)
	ratio = ratio**-1
	if ratio_property is not None:
		ratio = get_merger_ratios(merger_list, ratio_property)
	if convert_to_time is True:
		from .. import cosmology
		x = np.array([cosmology.getTime(zz, sim) for zz in z])
	else:
		x = z

	if ratio_range is None:
			ratio_range=[0.1,1]

	for i in range(len(ratio)):
		if ratio[i] < ratio_range[0] or ratio[i] > ratio_range[1]:
			continue
		if alpha is None:
			alpha = 0.9/(ratio_range[1]-ratio_range[0]) * (ratio[i] - ratio_range[0]) + 0.1
		if alpha > 1: alpha = 1
		if alpha < 0: alpha = 0
		plotting.plt.plot([x[i],x[i]],y_range,alpha=alpha,**kwargs)
