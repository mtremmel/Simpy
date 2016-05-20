import numpy as np
from .plotting import plt

def hist_by_step_intervals(dbsim, data, type='time',weight=None, **kwargs):
	t = []
	z = []
	for step in dbsim.timesteps:
		t.append(step.time_gyr)
		z.append(step.redshift)
	t = np.array(t)
	z = np.array(z[::-1])
	if weight is not None:
		if type(weight)!=list and type(weight)!=np.ndarray:
			weight = np.ones(len(data))*weight
	if type=='time':
		binedge = t
	if type=='redshift':
		binedge=z

	plt.hist(data,bins=binedge,weight=weight,**kwargs)


