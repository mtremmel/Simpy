import numpy as np
from . import util
import os
import gc


def rockstar_iord_to_fpos(snap):
	iord_to_fpos = util.init_iord_to_fpos(snap)
	f = open(snap+'rockstar.halo_particles','rb')
	orig = np.fromfile(f,dtype=np.dtype('int64'),count=-1)
	f.close()
	return iord_to_fpos[orig]

def convert_all_rockstar():
	f = open('files.list','r')
	files = f.readlines()
	f.close()
	for ff in files:
		snap = ff.strip('\n')
		print snap
		if not os.path.exists(snap+'.rockstar.halo_particles'):
			print "snapshot does not have halo_particles file... ignoring this snapshot..."
			continue
		fpos = rockstar_iord_to_fpos(snap)
		if type(fpos[-1]) != np.int64:
			fpos = fpos.astype(np.int64)
		outf = open(snap+'.rockstar.halo_particles_fpos','wb')
		fpos.tofile(outf, format='%l')
		del(fpos)
		gc.collect()
		outf.close()

