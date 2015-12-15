import numpy as np

def property_array(simname, step, haloids, keylist):
	import halo_db as db
	idsort = np.argsort(haloids)
	uhids, ind = np.unique(haloids[idsort],return_index=True)
	outputlist = {}
	for key in keylist:
		outputlist[key] = np.ones(len(haloids)) * -1
	for ii in range(len(uhids)):
		id = uhids[ii]
		h = db.get_halo(simname+'/%'+str(step)+'/'+str(id))
		if ii +1 < len(uhids):
			ss = idsort[ind[ii]:ind[ii+1]]
		else:
			ss = idsort[ind[ii]:]
		for key in keylist:
			if key not in h.keys():
				continue
			outputlist[key][ss] = h[key]

	return tuple([outputlist[key] for key in outputlist.keys()])