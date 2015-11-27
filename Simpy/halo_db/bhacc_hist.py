import numpy as np
import halo_db as db

def gen_bh_acc(simname, endstep, halonum, bhorbit)
	hcur = db.get_halo(simname+'/%'+endstep+'/'+str(halonum))
	h = hcur
	time = np.zeros(endstep)
	lummax = np.zeros(endstep)
	lumtot = np.zeros(endstep)
	bhorbit.getprogbhs()
	while h.previous:
		tcur = h.timestep.time_gyr
		tnext = h.previous.time_gyr
		bhids = np.array([bh.halo_number for bh in h['BH']])
		bhids_prev = np.array([bh.halo_number for bh in h.previous['BH']])
		orbitmatch, = np.where(np.in1d(bhorbit.bhiords, bhids))
		prevmatch, = np.where(np.in1d(bhids, bhids_prev))
		nomatch, = np.where(!np.in1d(bhids, bhids_prev))
		omatch, = np.where(np.in1d(bhorbit.bhiords,bhids[prevmatch]))
		onomatch, = np.where(np.in1d(bhorbit.bhiords,bhids[nomatch]))

		luminosity = np.array([])
		time = np.array([])
		step = np.array([])
		for id in bhids[prevmatch]:
			lum = bhorbit.single_BH_data(id, 'lum')
			t = bhorbit.single_BH_data(id, 'time')
			s = bhorbit.single_BH_data(id, 'step')
			np.append(luminosity, lum[((t>tnext)&(t<=tcur))])
			np.append(time, t[((t>tnext)&(t<=tcur))])
			np.append(step, s[((t>tnext)&(t<=tcur))])

		progs = []
		progt = []
		for i in onomatch:
			if len(bhorbit.prog['time'][i])==0: continue
			ok = np.where((np.array(bhorbit.prog['time'][i]) >= tnext)&(np.array(bhorbit.prog['time'][i]) < tnext))
			if len(ok) > 0:
				progmatch, = np.where(np.inqd(bhids_prev,bhorbit.prog['iord'][i][ok]))
				for i in progmatch:








