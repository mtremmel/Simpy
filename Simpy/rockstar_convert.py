import numpy as np
import pynbody
from . import util, Files
import os
import gc


def rockstar_iord_to_fpos(snapn):
	snap = pynbody.load(snapn)
	iord_to_fpos = util.init_iord_to_fpos(snap)
	del(snap)
	gc.collect()
	f = open(snapn+'.rockstar.halo_particles', 'rb')
	orig = np.fromfile(f, dtype=np.dtype('int64'), count=-1)
	f.close()
	return iord_to_fpos[orig]

def convert_all_rockstar(simname):
	Files.cklists(simname)
	f = open('steps.list', 'r')
	steps = f.readlines()
	f.close()
	for ff in steps:
		snapn = ff.strip('\n')
		snap = simname+'.'+snapn
		print(snap)
		if not os.path.exists(snap+'.rockstar.halo_particles'):
			print("snapshot does not have halo_particles file... ignoring this snapshot...")
			continue
		if os.path.exists(snap+'.rockstar.halo_particles_fpos'):
			print("fpos file already found! skipping...")
			continue
		fpos = rockstar_iord_to_fpos(snap)
		if not isinstance(fpos[-1], np.int64):
			fpos = fpos.astype(np.int64)
		outf = open(snap+'.rockstar.halo_particles_fpos', 'wb')
		fpos.tofile(outf, format='%l')
		del(fpos)
		gc.collect()
		outf.close()

