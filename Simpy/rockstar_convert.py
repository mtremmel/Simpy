import numpy as np
from . import util
import os


def rockstar_iord_to_fpos(snap):
	iord_to_fpos = util.init_iord_to_fpos(snap)
	f = open(snap+'rockstar.halo_partiels',rb)
	orig = np.fromfile(f,dtype=np.dtype('int64'),count=-1)
	f.close()
	return iord_to_fpos[orig]

def convert_all_rockstar():
	f = open('files.list','r')
	files = f.readlines()
	f.close()
	for ff in files:
		print ff
		snap = ff.strip('\n')
		if not os.path.exists(snap+'rockstar.halo_partiels'):
			print "snapshot does not have halo_particles file... ignoring this snapshot..."
			continue
		fpos = rockstar_iord_to_fpos(snap)
