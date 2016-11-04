import tangos as db
import tangos.relation_finding
import numpy as np


def get_merger_ratios(merger_list, property):
	ratio = [min(mm[0][property],mm[1][property])/max(mm[0][property],mm[1][property]) for mm in merger_list]
	return np.array(ratio)


def sSFR(halo):
	SFH = halo.calculate("reassemble(SFR_histogram)")
	Mstar, time = halo.reverse_property_cascade("Mstar","t()")
	t_sfh = np.arange((len(SFH)))*0.01
	Mstar_interp = np.interp(t_sfh,time[::-1],Mstar[::-1])
	sSFR = SFH/Mstar_interp

	return sSFR, t_sfh


def plot_merger_times(halo, y_range, sim = None, ratio_property=None, convert_to_time=True, ratio_min=0.1, ratio_norm=1.0, **kwargs):
	import tangos.examples.mergers as mergers
	from .. import plotting
	z, ratio, merger_list = mergers.get_mergers_of_major_progenitor(halo)
	if ratio_property is not None:
		ratio = get_merger_ratios(merger_list, ratio_property)
	if convert_to_time is True:
		from .. import cosmology
		x = np.array([cosmology.getTime(zz, sim) for zz in z])
	else:
		x = z

	for i in range(len(ratio)):
		if ratio[i] < ratio_min:
			continue
		alpha = ratio[i]/ratio_norm
		if alpha > 1: alpha = 1
		print alpha
		plotting.plt.plot([x,x],y_range,alpha=alpha,**kwargs)
