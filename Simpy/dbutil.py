import numpy as np
from .plotting import plt

def list_db_time(dbsim):
	t = []
	for step in dbsim.timesteps:
		t.append(step.time_gyr)
	return np.array(t)

def list_db_redshift(dbsim):
	z = []
	for step in dbsim.timesteps:
		z.append(step.redshift)
	return np.array(z)


def hist_by_step_intervals(dbsim, data, type='time',weight=None, dxnorm=True, ret_data=False, plot=True, **kwargs):

	if weight is not None:
		if type(weight)!=list and type(weight)!=np.ndarray:
			weight = np.ones(len(data))*weight

	if type=='time':
		binedge = list_db_time(dbsim)
	if type=='redshift':
		binedge = list_db_redshift(dbsim)[::-1]


	n, bins = np.histogram(data,bins=binedge,weights=weight,**kwargs)

	if dxnorm is True:
		n = n/(bins[1:]-bins[:-1])

	if plot is True:
		plt.step(bins,np.append(n,n[-1]),**kwargs)

	if ret_data is True:
		return n, bins
	else:
		return


class PopEvol(object):
	def __init__(self, step, hidlist, proplist, type='backward'):
		if 'z()' not in proplist:
			proplist.append('z()')
		if 't()' not in proplist:
			proplist.append('t()')

		self.data = {'hid':[]}
		for p in proplist:
			self.data[p] = []

		for hid in hidlist:
			print "collecting data for halo ", hid
			self.data['hid'].append(hid)
			if type == 'backward':
				rawdat =  step.halos[hid-1].reverse_property_cascade(*proplist)
			if type == 'forward':
				rawdat =  step.halos[hid-1].property_cascae(*proplist)

			cnt = 0
			for p in proplist:
				self.data[p].append(rawdat[cnt])
				cnt += 1

	def get_data_zrange(self,zmin,zmax, key):
		output = []
		for i in range(len(self.data['hid'])):
			ztarget = np.where((self.data['z()'][i]>=zmin)&(self.data['z()'][i]<zmax))[0]
			output.append(self.data[key][i][ztarget])

		return np.array(output)



