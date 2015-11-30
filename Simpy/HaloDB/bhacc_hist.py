import numpy as np
import pynbody


def track_halo_bh_acc(simname, endstep, halonum, bhorbit, active=1e42):
	import halo_db as db
	hcur = db.get_halo(simname+'/%'+str(endstep)+'/'+str(halonum))
	h = hcur
	if 'BH' not in h.keys():
		print "ERROR No BHs found in halo "+str(halonum)+" at step "+endstep+"\n"
		return {}

	bhorbit.getprogbhs()
	luminosity = np.array([])
	time = np.array([])
	step = np.array([])
	iord = np.array([])
	while h.previous:
		tcur = h.timestep.time_gyr
		tnext = h.previous.timestep.time_gyr
		print h.previous, h, tcur, tnext

		nbhs_prev = len(np.where(np.array(h.previous.keys())=='BH')[0])
		if nbhs_prev > 1:
			bhids_prev = np.array([bh.halo_number for bh in h.previous['BH']])
		if nbhs_prev == 1:
			bhids_prev = np.array([h.previous['BH'].halo_number])
		if nbhs_prev == 0:
			print "No BHs in previous timestep!"
			for id in bhids:
				lum = bhorbit.single_BH_data(id, 'lum')
				t = bhorbit.single_BH_data(id, 'time')
				s = bhorbit.single_BH_data(id, 'step')
				ids = np.ones(len(lum)).astype(np.int64)*id
				luminosity = np.append(luminosity, lum[(t<=tcur)])
				time = np.append(time, t[(t<=tcur)])
				step = np.append(step, s[(t<=tcur)])
				iord = np.append(iord,ids[(t<=tcur)])
			break

		prevmatch, = np.where(np.in1d(bhids, bhids_prev))
		nomatch, = np.where(np.in1d(bhids, bhids_prev)==False)
		onomatch, = np.where(np.in1d(bhorbit.bhiords,bhids[nomatch]))

		for id in bhids[prevmatch]:
			lum = bhorbit.single_BH_data(id, 'lum')
			t = bhorbit.single_BH_data(id, 'time')
			s = bhorbit.single_BH_data(id, 'step')
			ids = np.ones(len(lum)).astype(np.int64)*id
			luminosity = np.append(luminosity, lum[((t>tnext)&(t<=tcur))])
			time = np.append(time, t[((t>tnext)&(t<=tcur))])
			step = np.append(step, s[((t>tnext)&(t<=tcur))])
			iord = np.append(iord,ids[((t>tnext)&(t<=tcur))])

		for i in onomatch:
			if len(bhorbit.prog['time'][i])==0: continue
			ok, = np.where((np.array(bhorbit.prog['time'][i]) >= tnext)&(np.array(bhorbit.prog['time'][i]) < tnext))
			if len(ok) > 0:
				tfirst = tcur
				progmatch, = np.where(np.in1d(bhids_prev,bhorbit.prog['iord'][i][ok]))
				if len(progmatch)==0: continue
				for id in bhids_prev[progmatch]:
					lum = bhorbit.single_BH_data(id, 'lum')
					t = bhorbit.single_BH_data(id, 'time')
					s = bhorbit.single_BH_data(id, 'step')
					ids = np.ones(len(lum)).astype(np.int64)*id
					if t.max() > tcur:
						print "WARNING! time found for eated BH that excedes current time... should not be possible"
					if t.max() <= tfirst:
						tfirst = t.max()
					luminosity = np.append(luminosity, lum[((t>tnext)&(t<=tcur))])
					time = np.append(time, t[((t>tnext)&(t<=tcur))])
					step = np.append(step, s[((t>tnext)&(t<=tcur))])
					iord = np.append(iord,ids[((t>tnext)&(t<=tcur))])
				lum = bhorbit.single_BH_data(bhorbit.bhiords[i], 'lum')
				t = bhorbit.single_BH_data(bhorbit.bhiords[i], 'time')
				s = bhorbit.single_BH_data(bhorbit.bhiords[i], 'step')
				ids = np.ones(len(lum)).astype(np.int64)*id
				luminosity = np.append(luminosity, lum[((t>tfirst)&(t<=tcur))])
				time = np.append(time, t[((t>tfirst)&(t<=tcur))])
				step = np.append(step, s[((t>tfirst)&(t<=tcur))])
				iord = np.append(iord,ids[((t>tfirst)&(t<=tcur))])
		h = h.previous

	outstep = np.unique(step)
	outtime = np.zeros(len(outstep))
	maxlum = np.zeros(len(outstep))
	totlum = np.zeros(len(outstep))
	nactive = np.zeros(len(outstep))
	idmax = np.zeros(len(outstep))
	for i in range(len(outstep)):
		curt, = np.where(step == outstep[i])
		outtime[i] = np.mean(time[curt])
		maxlum[i] = luminosity[curt].max()
		totlum[i] = luminosity[curt].sum()
		act, = np.where(luminosity[curt]>active)
		nactive[i] = len(act)
		idmax[i] = iord[curt][np.argmax(luminosity[curt])]

	AccHist = {'maxlum':pynbody.array.SimArray(maxlum,'erg s**-1'),
			'totlum':pynbody.array.SimArray(totlum,'erg s**-1'),
			'step':outstep,
			'time':pynbody.array.SimArray(outtime,'Gyr'),
			'nactive':nactive,
			'idmaxlum': idmax,
			}

	return AccHist

def total_halo_bh_acc(simname, endstep, halonum, bhorbit, active=1e42):
	import halo_db as db
	from .. import BlackHoles
	h = db.get_halo(simname+'/%'+str(endstep)+'/'+str(halonum))
	tcur = h.timestep.time_gyr
	nbhs = len(np.where(np.array(h.keys())=='BH')[0])
	if nbhs > 1:
		bhids = np.array([bh.halo_number for bh in h['BH']])
	if nbhs == 1:
		bhids = np.array([h['BH'].halo_number])
	if nbhs == 0:
		return
	omatch, = np.where(np.in1d(bhorbit.bhiords, bhids))
	bhorbit.getprogbhs()
	cnt = 0

	for i in omatch:
		cnt += 1
		print "finding total progenitor list for BH ", bhorbit.bhiords[i], "("+str(cnt)+'/'+str(len(omatch))+")"
		proglist = BlackHoles.mergers.get_complete_prog_list(bhorbit,bhorbit.bhiords[i], tcur)
		bhids = bhids.append(np.array(proglist))

	print "getting total accretion history..."
	luminosity = np.array([])
	time = np.array([])
	step = np.array([])
	iord = np.array([])
	for bid in bhids:
		lumpart = bhorbit.single_BH_data(bid, 'lum')
		timepart = bhorbit.single_BH_data(bid, 'time')
		steppart = bhorbit.single_BH_data(bid, 'step')
		luminosity = np.append(luminosity, lumpart)
		step = np.append(step, steppart)
		time = np.append(time, timepart)
		ids = np.ones(len(lumpart)).astype(np.int64)*bid
		iord = np.append(iord, ids)

	outstep = np.unique(step)
	outtime = np.zeros(len(outstep))
	maxlum = np.zeros(len(outstep))
	totlum = np.zeros(len(outstep))
	nactive = np.zeros(len(outstep))
	idmax = np.zeros(len(outstep))

	for i in range(len(outstep)):
		curt, = np.where(step == outstep[i])
		outtime[i] = np.mean(time[curt])
		maxlum[i] = luminosity[curt].max()
		totlum[i] = luminosity[curt].sum()
		act, = np.where(luminosity[curt]>active)
		nactive[i] = len(act)
		idmax[i] = iord[curt][np.argmax(luminosity[curt])]

	AccHist = {'maxlum':pynbody.array.SimArray(maxlum,'erg s**-1'),
			'totlum':pynbody.array.SimArray(totlum,'erg s**-1'),
			'step':outstep,
			'time':pynbody.array.SimArray(outtime,'Gyr'),
			'nactive':nactive,
			'idmaxlum': idmax,
			}

	return AccHist













