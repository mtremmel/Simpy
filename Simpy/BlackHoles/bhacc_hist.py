import numpy as np
import pynbody


def track_halo_bh_acc(simname, endstep, halonum, bhorbit, active=1e42, type='Central'):
	import halo_db as db
	if type not in ['Central', 'All']:
		print "WARNING type argument not understood (Central or All). Assuming type to be All"
		type = 'All'
	hcur = db.get_halo(simname+'/%'+str(endstep)+'/'+str(halonum))
	h = hcur
	if type == 'Central':
		dictstr = 'BH_central'
	else:
		dictstr = 'BH'
	nbhs = len(np.where(np.array(h.keys())==dictstr)[0])
	if nbhs > 1:
		bhids = np.array([bh.halo_number for bh in h[dictstr]])
	if nbhs == 1:
		bhids = np.array([h[dictstr].halo_number])
	if nbhs == 0:
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

		nbhs_prev = len(np.where(np.array(h.previous.keys())==dictstr)[0])
		if nbhs_prev > 1:
			bhids_prev = np.array([bh.halo_number for bh in h.previous[dictstr]])
		if nbhs_prev == 1:
			bhids_prev = np.array([h.previous[dictstr].halo_number])
		if nbhs_prev == 0:
			print "No BHs in previous timestep!"
			for bid in bhids:
				lum = bhorbit.single_BH_data(bid, 'lum')
				t = bhorbit.single_BH_data(bid, 'time')
				s = bhorbit.single_BH_data(bid, 'step')
				ids = np.ones(len(lum)).astype(np.int64)*bid
				luminosity = np.append(luminosity, lum[(t<=tcur)])
				time = np.append(time, t[(t<=tcur)])
				step = np.append(step, s[(t<=tcur)])
				iord = np.append(iord,ids[(t<=tcur)])
			break

		prevmatch, = np.where(np.in1d(bhids, bhids_prev))
		noprevmatch, = np.where(np.in1d(bhids, bhids_prev)==False)
		omatch, = np.where(np.in1d(bhorbit.bhiords, bhids))

		for bid in bhids[prevmatch]:
			lum = bhorbit.single_BH_data(bid, 'lum')
			t = bhorbit.single_BH_data(bid, 'time')
			s = bhorbit.single_BH_data(bid, 'step')
			ids = np.ones(len(lum)).astype(np.int64)*bid
			luminosity = np.append(luminosity, lum[((t>tnext)&(t<=tcur))])
			time = np.append(time, t[((t>tnext)&(t<=tcur))])
			step = np.append(step, s[((t>tnext)&(t<=tcur))])
			iord = np.append(iord,ids[((t>tnext)&(t<=tcur))])

		for i in omatch:
			if len(bhorbit.prog['time'][i])==0: continue
			ok, = np.where((np.array(bhorbit.prog['time'][i]) >= tnext)&(np.array(bhorbit.prog['time'][i]) <= tcur))
			if len(ok) > 0:
				tfirst = tcur
				progmatch, = np.where(np.in1d(bhids_prev,np.array(bhorbit.prog['iord'][i])[ok]))
				if len(progmatch) == 0: continue

				for bid in bhids_prev[progmatch]:
					lum = bhorbit.single_BH_data(bid, 'lum')
					t = bhorbit.single_BH_data(bid, 'time')
					s = bhorbit.single_BH_data(bid, 'step')
					ids = np.ones(len(lum)).astype(np.int64)*bid
					if t.max() > tcur:
						print "WARNING! time found for eaten BH that excedes current time... should not be possible"
					if t.max() <= tfirst:
						tfirst = t.max()
					luminosity = np.append(luminosity, lum[((t>tnext)&(t<=tcur))])
					time = np.append(time, t[((t>tnext)&(t<=tcur))])
					step = np.append(step, s[((t>tnext)&(t<=tcur))])
					iord = np.append(iord,ids[((t>tnext)&(t<=tcur))])

				if bhorbit.bhiords[i] not in bhids_prev:
					lum = bhorbit.single_BH_data(bhorbit.bhiords[i], 'lum')
					t = bhorbit.single_BH_data(bhorbit.bhiords[i], 'time')
					s = bhorbit.single_BH_data(bhorbit.bhiords[i], 'step')
					ids = np.ones(len(lum)).astype(np.int64)*bhorbit.bhiords[i]
					luminosity = np.append(luminosity, lum[((t>tfirst)&(t<=tcur))])
					time = np.append(time, t[((t>tfirst)&(t<=tcur))])
					step = np.append(step, s[((t>tfirst)&(t<=tcur))])
					iord = np.append(iord,ids[((t>tfirst)&(t<=tcur))])

		bhids = bhids_prev
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

def total_halo_bh_acc(simname, endstep, halonum, bhorbit, active=1e42, type='Central'):
	import halo_db as db
	from .. import BlackHoles
	if type not in ['Central', 'All']:
		print "WARNING type argument not understood (Central or All). Assuming type to be All"
		type = 'All'
	h = db.get_halo(simname+'/%'+str(endstep)+'/'+str(halonum))
	tcur = h.timestep.time_gyr
	if type == 'Central':
		dictstr = 'BH_central'
	else:
		dictstr = 'BH'
	nbhs = len(np.where(np.array(h.keys())==dictstr)[0])
	if nbhs > 1:
		bhids = np.array([bh.halo_number for bh in h[dictstr]])
	if nbhs == 1:
		bhids = np.array([h[dictstr].halo_number])
	if nbhs == 0:
		return
	omatch, = np.where(np.in1d(bhorbit.bhiords, bhids))
	bhorbit.getprogbhs()
	cnt = 0

	for i in omatch:
		cnt += 1
		print "finding total progenitor list for BH ", bhorbit.bhiords[i], "("+str(cnt)+'/'+str(len(omatch))+")"
		proglist = BlackHoles.mergers.get_complete_prog_list(bhorbit,bhorbit.bhiords[i], tcur)
		bhids = np.append(bhids,proglist)

	print "getting total accretion history..."
	luminosity = np.array([])
	time = np.array([])
	step = np.array([])
	iord = np.array([])
	for bid in bhids:
		if len(np.where(np.in1d(bhorbit.bhiords,bid))[0]) == 0:
			continue
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

class BHAccHist(object):
	def __init__(self, simname, endstep, halonum, bhorbit, active=1e42):
		self.track_central = track_halo_bh_acc(simname, endstep, halonum, bhorbit, active=active, type='Central')
		self.track_all = track_halo_bh_acc(simname, endstep, halonum, bhorbit, active=active, type='All')
		self.total_central = total_halo_bh_acc(simname, endstep, halonum, bhorbit, active=active, type='Central')
		self.total_all = total_halo_bh_acc(simname, endstep, halonum, bhorbit, active=active, type='All')

	def plot_bh_hist(self, track=True, halo='Central', qty = 'totlum', linestyle='-', lw=3, label=None, color='blue'):
		from .. import plotting
		if halo not in ['Central','All']:
			print "WARNING halo input not understood (Central or All). Assuming it is Central"
			halo = 'Central'
		if qty not in ['totlum','maxlum', 'nactive']:
			print "WARNING type input not understood(totlum or maxlum or nactive). Assuming totlum."
			type = 'totlum'

		if track is True:
			if halo=='Central':
				plotting.plt.plot(self.track_central['time'], self.track_central[qty],
								  linestyle=linestyle, lw=lw, label=label, color=color)
			if halo=='All':
				plotting.plt.plot(self.track_all['time'], self.track_all[qty],
								  linestyle=linestyle, lw=lw, label=label, color=color)

		if track is False:
			if halo=='Central':
				plotting.plt.plot(self.total_central['time'], self.total_central[qty],
								  linestyle=linestyle, lw=lw, label=label, color=color)
			if halo=='All':
				plotting.plt.plot(self.total_all['time'], self.total_all[qty],
								  linestyle=linestyle, lw=lw, label=label, color=color)

		plotting.plt.legend()
		if qty == 'nactive':
			plotting.plt.ylabel('N(>1e43 ergs/s)')
		else:
			plotting.plt.ylabel('Luminosity (ergs/s)')

		plotting.plt.xlabel('Time (Gyr)')












