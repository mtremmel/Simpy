import numpy as np


def property_array(simname, step, haloids, keylist):
	import halo_db as db
	uhids, inv = np.unique(haloids, return_inverse=True)
	outputlist = {}
	udatalist = {}
	for key in keylist:
		udatalist[key] = []

	for ii in range(len(uhids)):
		id = uhids[ii]
		if id < 0:
			for key in keylist:
				udatalist[key].append(np.nan)
			continue
		h = db.get_halo(simname + '/%' + str(step) + '/' + str(id))
		for key in keylist:
			try:
				udatalist[key].append(h[key])
			except KeyError:
				udatalist[key].append(np.nan)
	for key in keylist:
		if type(udatalist[key][0]) == np.ndarray:
			outputlist[key] = np.vstack(udatalist[key])[inv]
		else:
			outputlist[key] = np.array(udatalist[key])[inv]

	return tuple([outputlist[key] for key in outputlist.keys()])
