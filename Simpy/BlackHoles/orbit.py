import numpy as np
import pynbody

import os
import pickle
from .. import readcol
import gc
from .. import Files
from .. import util
from .. import cosmology


def mkAbridgedOrbit(simname, sim, endname, lmin=1e43, mmin=1e6):
	munits = sim.infer_original_units('Msol')
	tunits = sim.infer_original_units('Gyr')
	mdotunits = munits / tunits

	Mlimitsim = mmin / munits.in_units('Msol')
	mdotlimit = lmin / (0.1 * 3e10 * 3e10)
	mdotlimit /= mdotunits.in_units('g s**-1')
	cstr = """ awk '{if ($4 - $13> """ + str(Mlimitsim) + """ && $12 > """ + str(
		mdotlimit) + """) print  $4 " " $12 " " $13 " " $15 " " $16}' """ + simname + ".orbit > " + simname + ".orbit.abridged." + endname
	os.system(cstr)
	return


def sepOrbitbyStep(simname, minstep=0, maxstep=1000000000, MBHinit=1e6, NCHILADA=True):
	Files.cklists(simname, NCHILADA=NCHILADA)
	f = open('files.list', 'r')
	sim = pynbody.load(f.readline().strip('\n'))
	f.close()
	munits = sim.infer_original_units('Msol')
	MBHinit = MBHinit / float(munits)
	if not os.path.exists('orbitsteps'): os.system('mkdir orbitsteps')
	os.chdir('orbitsteps')
	print "separating orbit file by step..."
	os.system(
		"awk -F ' ' '{if($4 - $13 > " + str(MBHinit) + " && $3!=int($3) && $3 >= " + str(minstep) + " && $3 <= " + str(
			maxstep) + ") print>(int($3)+1); if($4 - $13 > " + str(MBHinit) + " && $3==int($3) && $3 >= " + str(
			minstep) + " && $3 <= " + str(maxstep) + ") print>$3}' " + '../' + simname + '.orbit')
	os.chdir('../')
	return


def extract_single_BH(simname, bhiord):
	os.system("awk -F ' ' '{if($1 == " + str(bhiord) + ") print>$1}' " + simname + ".orbit")
	os.system("mv " + str(bhiord) + "orbit."+str(bhiord))
	return


def read_full_orbit_file(filename, simname):
	f = open('files.list','r')
	files = f.readlines()
	s = pynbody.load(files[0].strip('\n'))
	f.close()
	a = None
	try:
		bhid, time, step, mass, x, y, z, vx, vy, vz, pot, mdot, dm, E, dt, a = readcol(filename, twod=False)
	except:
		bhid, time, step, mass, x, y, z, vx, vy, vz, pot, mdot, dm, E, dt = readcol(filename, twod=False)
	output = {'iord':bhid, 'time':time, 'step':step, 'mass':mass, 'x':x, 'y':y,
			  'z':z, 'vx':vx, 'vy':vy, 'vz':vz, 'pot':pot, 'mdot':mdot, 'dm':dm, 'E':E, 'dt':dt, 'a':a}
	if a is None:
		a, = cosmology.getScaleFactor(pynbody.array.SimArray(time,s.infer_original_units('Gyr')),s)
		output['a'] = a
	units = {'x':'kpc', 'y':'kpc', 'z':'kpc', 'vx':'km s**-1', 'vy':'km s**-1', 'vz':'km s**-1',
			 'mdot':'Msol yr**-1', 'dm':'Msol', 'dt':'Gyr', 'time':'Gyr', 'mass':'Msol'}

	for key in ['x', 'y', 'z', 'vx', 'vy', 'vz']:
		print "fixing scale factor for", key
		output[key] *= output['a']
	for key in units.keys():
		print "converting units for ", key
		origunit = s.infer_original_units(units[key])
		if key in ['x', 'y', 'z', 'vx', 'vy', 'vz']: origunit = origunit / pynbody.units.Unit('a')
		output[key] = pynbody.array.SimArray(output[key],origunit)
		output[key] = output[key].in_units(units[key])
	return output

def smooth_raw_orbit_data(output, key, nsteps, maxstep=4096, minstep=0):
	ok = np.where((output['step']>=minstep)&(output['step']<maxstep))[0]
	output = util.cutdict(output,ok)
	ord = np.argsort(output['step'])
	output = util.cutdict(output,ord)
	ustep, ind = np.unique(output['step'].astype(np.int),return_index=True)
	ss = np.where(ustep%nsteps==0)
	smoothed_dat = []
	stdev = []
	time = []
	for i in range(len(ss)):
		data = output[key][ind[ss[i]]:ind[ss[i+1]]]
		dt = output['dt'][ind[ss[i]]:ind[ss[i+1]]]
		xmean, xsd, xsum = util.timeweightedAve(data, dt)
		smoothed_dat.append(xmean)
		stdev.append(xsd)
		time.append((output['time'][ind[ss[i]]]+output['time'][ind[ss[i+1]]])/2.)
	return np.array(smoothed_dat), np.array(stdev), np.array(time)





def getOrbitValbyStep(minstep=1, maxstep=4096, clean=False, filename=None, ret_output=False):
	'''
	simname = name of simulation (i.e. simname.004096)
	minstep, maxstep = min, max steps the code will expect to have in the orbit file.
	clean = True/False, clean up the extra files made from the oritinal orbit file (this is done to spare memory, but they take up lots of space...
		'''
	output = {'iord': [], 'time': [], 'step': [], 'mass': [], 'x': [], 'y': [], 'z': [], 'vx': [], 'vy': [], 'vz': [],
			  'mdot': [], 'mdotmean': [], 'mdotsig': [], 'a': [], 'dM': []}
	oldform = False
	f = open('files.list','r')
	files = f.readlines()
	s = pynbody.load(files[0].strip('\n'))
	f.close()
	if not os.path.exists('orbitsteps/'):
		print "ERROR! can't find orbitsteps... run sepOrbitbyStep first!"
		return
	for i in range(minstep, maxstep + 1):
		if os.path.exists('orbitsteps/' + str(i)):
			print "getting data for step ", i
			try:
				bhid, time, step, mass, x, y, z, vx, vy, vz, pot, mdot, dm, E, dt, a = readcol.readcol(
				'orbitsteps/' + str(i), twod=False)
			except:
				oldform = True
				bhid, time, step, mass, x, y, z, vx, vy, vz, pot, mdot, dm, E, dt = readcol.readcol(
				'orbitsteps/' + str(i), twod=False)
			if clean == True: os.system('rm orbitsteps/' + str(i))
		else:
			continue
		ord = np.argsort(bhid)
		bhid = bhid[ord].astype(np.int64)
		time = time[ord]
		step = step[ord]
		mass = mass[ord]
		x = x[ord]
		y = y[ord]
		z = z[ord]
		vx = vx[ord]
		vy = vy[ord]
		vz = vz[ord]
		mdot = mdot[ord]
		dt = dt[ord]
		if oldform is False:
			a = a[ord]
		else:
			a,redsh= cosmology.getScaleFactor(pynbody.array.SimArray(time,s.infer_original_units('Gyr')),s)
			del(redsh)
			gc.collect()
		# bad, = np.where((mass - mdot*dt < MBHinit)|(mass<MBHinit))
		# mdot[bad] = 0
		ubhid, uind = np.unique(bhid, return_index=True)
		for ii in range(len(ubhid)):
			if ii < len(ubhid) - 1:
				idind = range(uind[ii], uind[ii + 1])
			else:
				idind = range(uind[ii], len(step))
			curstep, = np.where(step[idind] == i)
			if len(curstep) == 0: continue
			curstep = curstep[0]
			output['iord'].append(ubhid[ii])
			output['step'].append(step[idind][curstep])
			output['mass'].append(mass[idind][curstep])
			output['time'].append(time[idind][curstep])
			output['x'].append(x[idind][curstep])
			output['y'].append(y[idind][curstep])
			output['z'].append(z[idind][curstep])
			output['vx'].append(vx[idind][curstep])
			output['vy'].append(vy[idind][curstep])
			output['vz'].append(vz[idind][curstep])
			output['mdot'].append(mdot[idind][curstep])
			output['a'].append(a[idind][curstep])
			mean, std, dM = util.timeweightedAve(mdot[idind], dt[idind])
			output['dM'].append(dM)
			output['mdotmean'].append(mean)
			output['mdotsig'].append(std)
	if filename:
		f = open(filename, 'wb')
		pickle.dump(output, f)
		f.close()
	if ret_output:
		return output


def truncOrbitFile(simname, minstep=1, maxstep=4096, ret_output=False):
	output = getOrbitValbyStep(minstep=minstep, maxstep=maxstep, ret_output=True)
	outorder = ['iord', 'time', 'step', 'mass', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'mdot', 'mdotmean', 'mdotsig', 'a', 'dM']
	tofile = []
	for key in outorder:
		tofile.append(output[key])
	if ret_output == False:
		del (output)
		gc.collect()
	tofile = tuple(tofile)
	outputname = simname + '.shortened.orbit'
	print "saving to file..."
	np.savetxt(outputname, np.column_stack(tofile),
			   fmt=['%d', '%f', '%d', '%e', '%f', '%f', '%f', '%f', '%f', '%f', '%e', '%e', '%e', '%f','%e'])
	del (tofile)
	gc.collect()
	if ret_output:
		return output
	return


def sticthOrbitSteps(simname, nfiles, ret_output=False, overwrite=False, nstart=1):

	outorder = ['iord', 'time', 'step', 'mass', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'mdot', 'mdotmean', 'mdotsig', 'a','dM']
	outputname = simname + '.shortened.orbit'

	if os.path.exists(outputname) and overwrite == False:
		outf = open(outputname,'a')
	else:
		outf = open(outputname,'w')
	for i in range(nfiles):
		num = i + nstart
		print "reading in data from ", simname + '.shortened.orbit' + str(num)
		f = open(simname + '.shortened.orbit' + str(num), 'rb')
		output_part = pickle.load(f)
		f.close()
		tofile = []
		for key in outorder:
			tofile.append(output_part[key])
			del(output_part[key])
			gc.collect()
		del(output_part)
		gc.collect()
		tofile = tuple(tofile)
		np.savetxt(outf, np.column_stack(tofile),
			   fmt=['%d', '%f', '%d', '%e', '%f', '%f', '%f', '%f', '%f', '%f', '%e', '%e', '%e', '%f','%e'])
		del(tofile)
		gc.collect()
	outf.close()
	return


class Orbit(object):
	def __init__(self, simname, savefile=None, NCHILADA=True):
		self.simname = simname
		self.nchil = NCHILADA
		self.filename=savefile
		ofile = simname + ".shortened.orbit"
		if not os.path.exists(ofile):
			print "ERROR shortened orbit file not found! Exiting..."
			return
		Files.cklists(simname, NCHILADA=NCHILADA)

		# read raw data from shortened orbit file
		print "reading in data. . ."
		bhid, time, step, mass, x, y, z, vx, vy, vz, mdot, mdotmean, mdotsig, scalefac, dM = readcol.readcol(ofile,
																										 twod=False,
																										 nanval=0.0)
		self.data = {'iord': bhid, 'time': time, 'step': step, 'mass': mass, 'x': x, 'y': y, 'z': z, 'vx': vx, 'vy': vy,
					 'vz': vz, 'mdot': mdot, 'mdotmean': mdotmean, 'mdotsig': mdotsig, 'scalefac': scalefac, 'dM': dM}
		del (bhid, time, step, mass, x, y, z, vx, vy, vz, mdot, mdotmean, mdotsig, scalefac, dM)
		gc.collect()

		# make each data object a simulation array with correct units
		defunits = {'iord': None, 'time': 'Gyr', 'step': None, 'mass': 'Msol', 'mdot': 'Msol yr**-1',
					'mdotmean': 'Msol yr**-1', 'mdotsig': 'Msol yr**-1', 'x': 'kpc', 'y': 'kpc', 'z': 'kpc',
					'vx': 'km s**-1', 'vy': 'km s**-1', 'vz': 'km s**-1', 'scalefac': None, 'dM': 'Msol'}
		f = open('files.list', 'r')
		sim = f.readline()
		s = pynbody.load(sim.strip('\n'))
		f.close()

		# control for "fake" bhs caused by restarts from outputs
		print "checking for fake BHs. . ."
		sl = pynbody.tipsy.StarLog(self.simname + '.starlog')
		slbhiords = sl['iord'][(sl['tform'] < 0)]
		ok, = np.where(np.in1d(self.data['iord'], slbhiords))
		for key in self.data.keys():
			self.data[key] = self.data[key][ok]
		self.data['iord'][(self.data['iord']<0)] = 2*2147483648 + self.data['iord'][(self.data['iord']<0)]
		# convert comoving quantities to physical units
		for key in ['x', 'y', 'z', 'vx', 'vy', 'vz']:
			self.data[key] *= self.data['scalefac']
		for key in self.data.keys():
			unit = None
			if defunits[key] is not None:
				unit = s.infer_original_units(defunits[key])
				if key in ['x', 'y', 'z', 'vx', 'vy', 'vz']: unit = unit / pynbody.units.Unit('a')
			self.data[key] = pynbody.array.SimArray(self.data[key], unit)
			if defunits[key] is not None:
				self.data[key].convert_units(defunits[key])

		# get information on iord,step data for easy future data recovery
		print "reducing data. . ."
		self.bhiords, self.id_slice = self._get_slice_ind('iord')
		self.steps, self.step_slice = self._get_slice_ind('step')
		self.times = np.unique(self.data['time'])
		self._calc_lum(er=0.1)
		self.get_all_BH_tform(sl)

		if savefile:
			self.save(filename=savefile)
		return


	def _calc_lum(self, er=0.1):
		csq = pynbody.array.SimArray((2.998e10) ** 2, 'erg g**-1')
		self.data['lum'] = self.data['mdotmean'].in_units('g s**-1') * csq * er
		return


	def _get_slice_ind(self, key, orderby='time'):
		'''
		return unique values as well as a list of indices of data with each value of self.data[key]
		key = (string) target key to do operation on
		'''
		ord_ = np.argsort(self.data[key])
		uvalues, ind = np.unique(self.data[key][ord_], return_index=True)
		slice_ = []
		for i in range(len(uvalues) - 1):
			ss = ord_[ind[i]:ind[i + 1]]
			if orderby:
				sort_ = np.argsort(self.data[orderby][ss])
				ss = ss[sort_]
			slice_.append(ss)
		ss = ord_[ind[i + 1]:]
		if orderby:
			sort_ = np.argsort(self.data[orderby][ss])
			ss = ss[sort_]
		slice_.append(ss)
		return uvalues, slice_

	def save(self,filename=None):
		if filename is None:
			if self.filename is not None: ff = self.filename
			else:
				print "WARNING! Filename not found. Using default"
				ff = 'BHorbit.pkl'
		else: ff = filename
		f = open(ff,'wb')
		pickle.dump(self,f)
		f.close()
		return


	def single_BH_data(self, iord, key):
		o, = np.where(self.bhiords == iord)
		slice_ = self.id_slice[o[0]]
		return self.data[key][slice_]

	def single_BH_data_smooth(self,iord,key,nsteps=10):
		rawdat = self.single_BH_data(iord, key)
		#time = self.single_BH_data(iord, 'time')
		nind = len(rawdat) - len(rawdat)%nsteps
		use = np.arange(nind)
		rawdat = rawdat[use].reshape((nind/nsteps,nsteps))
		#time = time[use].reshape((nind/nsteps,nsteps))
		meandat = rawdat.mean(axis=1)
		#meantime = time.mean(axis=1)
		return meandat



	def single_step_data(self, step, key):
		o, = np.where(self.steps == step)
		slice_ = self.step_slice[o[0]]
		return self.data[key][slice_]


	def get_all_BH_tform(self, sl):
		# sl = pynbody.tipsy.StarLog(self.simname+'.starlog')
		ord = np.argsort(sl['iord'])
		bhind, = np.where(np.in1d(sl['iord'][ord], self.bhiords))
		self.tform = sl['tform'][ord][bhind] * -1
		if self.tform.min() < 0: print "WARNING! Positive tforms were found for BHs!"
		del sl
		gc.collect()


	def getprogbhs(self):
		time, step, ID, IDeat, ratio, kick = readcol.readcol(self.simname + '.mergers', twod=False)
		self.prog = {'iord': [[] for i in range(len(self.bhiords))], 'kick': [[] for i in range(len(self.bhiords))],
					 'ratio': [[] for i in range(len(self.bhiords))], 'step': [[] for i in range(len(self.bhiords))],
					 'time': [[] for i in range(len(self.bhiords))]}
		o = np.argsort(ID)
		ID = ID[o]
		step = step[o]
		IDeat = IDeat[o]
		time = time[o]
		ratio = ratio[o]
		kick = kick[o]
		uID, ind = np.unique(ID, return_index=True)
		eaterind, = np.where(np.in1d(self.bhiords, uID))
		for i in range(len(uID) - 1):
			self.prog['iord'][eaterind[i]].extend(IDeat[ind[i]:ind[i + 1]])
			self.prog['step'][eaterind[i]].extend(step[ind[i]:ind[i + 1]])
			self.prog['time'][eaterind[i]].extend(time[ind[i]:ind[i + 1]])
			self.prog['ratio'][eaterind[i]].extend(ratio[ind[i]:ind[i + 1]])
			if kick.max() > 0:
				self.prog['kick'][eaterind[i]].extend(kick[ind[i]:ind[i + 1]])


	def gethalos(self):
		f = open('files.list', 'r')
		simfiles = f.readlines()
		nsnaps = len(simfiles)
		f.close()
		f = open('steps', 'r')
		snaps = f.readlines()
		f.close()
		initarr = np.ones((len(self.bhiords), len(snaps))) * -1
		self.bhhalos = {'Grp': initarr}
		for i in range(nsnaps):
			if os.path.exists(simfiles[i].strip('\n') + '.rockstar.grp'):
				grpfile = simfiles[i].strip('\n') + '.rockstar.grp'
			else:
				if os.path.exists(simfiles[i].strip('\n') + '.amiga.grp'):
					grpfile = simfiles[i].strip('\n') + '.amiga.grp'
				else:
					print "ERROR there is no grp file for this step!"
					continue
			sim = pynbody.load(simfiles[i].strip('\n'))
			simind = np.where(np.in1d(sim.stars['iord'], self.bhiords))
			orbind = np.where(np.inqd(self.bhiords, sim.stars['iord']))
			simind += len(sim.dm) + len(sim.gas)
			del sim['iord']
			gc.collect()
			grp = readcol.readcol(grpfile, skipline=1)
			self.bhhalos['Grp'][orbind][:, i] = grp[simind]
			del grp, simind, orbind
			gc.collect()


#Plotting functions


	def plt_single_BH_data(self, iord, keyx, unitx, keyy, unity, style, lw=1, msize=10, ylog=True, xlog=False,
						   label=None, overplot=False, smooth=False, nsmooth=10):
		from .. import plotting
		if smooth is False:
			ydat = self.single_BH_data(iord, keyy)
			xdat = self.single_BH_data(iord, keyx)
		else:
			ydat = self.single_BH_data_smooth(iord, keyy, nsteps=nsmooth)
			xdat = self.single_BH_data_smooth(iord, keyx, nsteps=nsmooth)
		plotting.plt.plot(xdat.in_units(unitx), ydat.in_units(unity), style, label=label, linewidth=lw,
						  markersize=msize)
		if xlog is True and overplot is False:
			plotting.plt.xscale('log', base=10)
		if ylog is True and overplot is False:
			plotting.plt.yscale('log', base=10)
		if label:
			plotting.plt.legend()
		return


	def plt_acc_hist(self, style, minM = 1e6, maxM = None, minL = 1e42, maxL = None, type='redshift',xlog=False,ylog=False, label=None, lw=1.5, volume=25**3, plotdata=True, overplot=False):
		from .. import plotting
		ord = np.argsort(self.data['scalefac'])
		if maxL is None: maxL = self.data['lum'].max()*10
		if minL is None: minL = 0
		if maxM is not None or minM is not None:
			if maxM is None: maxM = self.data['mass'].max()*10
			if minM is None: minM = 0
			okM, = np.where((self.data['mass'][ord] > minM)&(self.data['mass'][ord] <= maxM))
			ord = ord[okM]
		if minL is None or maxL is None:
			if maxL is None: maxL = self.data['lum'].max()*10
			if minL is None: minL = 0
			okL, = np.where((self.data['lum'][ord] > minL)&(self.data['lum'][ord] <= maxL))
			ord = ord[okL]

		macc = np.cumsum(self.data['dM'][ord])

		if type == 'redshift':
			z = self.data['scalefac']**-1 -1
			if xlog is False:
				plotting.plt.plot(z[ord],macc/volume,style, linewidth=lw, label=label)
			else:
				plotting.plt.plot(z[ord]+1,macc/volume,style, linewidth=lw, label=label)
				if overplot is False:
					plotting.plt.xticks([1,2,3,4,5,6,7,8,9,10],['0','1','2','3','4','5','6','7','8','9'])

			if plotdata is True:
				import colldata as dat
				err = dat.shankar09H - dat.shankar09
				plotting.plt.errorbar([1.03],[dat.shankar09],yerr=[err],color='black',fmt='D',label="Shankar+ 09")
				plotting.plt.errorbar([dat.Salvaterra12z+1],[dat.Salvaterra12],color='black',fmt='x',xerr=[dat.Salvaterra12zH-dat.Salvaterra12z],yerr=0.5*dat.Salvaterra12,uplims=[True],label='Salvaterra+ 12')
				plotting.plt.errorbar(dat.Treister13z,dat.Treister13,color='black',fmt='o',xerr=dat.Treister13zErr,yerr=0.5*dat.Treister13,uplims=[True,True,True], label='Treister+ 13')
				plotting.plt.errorbar(dat.Hopkins07zp1,10**dat.Hopkins07,color='grey',fmt='o',yerr=(dat.Hopkins07merr,dat.Hopkins07perr),label='Hopkins+ 07')

		if type == 'time' and plotdata is True:
			print "WARNING! Data only valid for redshift plotting. Ignoring keyword for time plot"
		if type== 'time':
			plotting.plt.plot(self.data['time'].in_units('Gyr'),macc/volume,style, linewidth=lw, label=label)
		if overplot is False:
			if xlog is True: plotting.plt.xscale('log',base=10)
			if ylog is True: plotting.plt.yscale('log',base=10)
			if type == 'redshift':
				plotting.plt.xlabel('z',fontsize=30)
			if type == 'time':
				plotting.plt.xlabel('Time (Gyr)',fontsize=30)
			plotting.plt.ylabel(r'log($\rho_{acc}$ [M$_{\odot}$ Mpc$^{-3}$])',fontsize=30)
		if label is not None or plotdata is True: plotting.plt.legend()
		return

	def plt_lumfun(self, style, minM = 1e6, maxM = None, minL = 1e42, maxL = None, volume=25**3,overplot=False,label=None,bins=50, redshift=1, plotdata=True, dotitle=True,lw=2):
		from .. import plotting
		from .. import cosmology
		import colldata as dat
		f = open('files.list', 'r')
		simfiles = f.readlines()
		f.close()
		sim = pynbody.load(simfiles[0].strip('\n'))
		tt = self.single_BH_data(self.bhiords[0],'time')
		dt = tt[1] - tt[0]
		zz = np.where((redshift > dat.hop_bhlf_zbinsl)&(redshift < dat.hop_bhlf_zbinsh))
		print "single timestep time: ", dt
		if maxM is None: maxM = self.data['mass'].max()*2
		if minM is None: minM = 0
		if maxL is None: maxL = self.data['lum'].max()*2
		if minL is None: minL = 0
		ok, = np.where((self.data['lum'] > minL)&(self.data['lum'] <= maxL)&(self.data['mass']>minM)&(self.data['mass']<maxM)&(self.data['scalefac']**-1 - 1 > dat.hop_bhlf_zbinsl[zz])&(self.data['scalefac']**-1 - 1 < dat.hop_bhlf_zbinsh[zz]))
		tlz = cosmology.getTime(dat.hop_bhlf_zbinsl[zz],sim)
		thz = cosmology.getTime(dat.hop_bhlf_zbinsh[zz],sim)
		T = tlz - thz
		if minL >0: lrange = [np.log10(minL), np.log10(maxL)]
		else: lrange = [np.log10(self.data['lum'].min()), np.log10(maxL)]
		dlogl = (lrange[1] - lrange[0])/float(bins)
		data = np.histogram(np.log10(self.data['lum'][ok]),bins=bins,range=lrange)
		phi = data[0]*(dt/(T*dlogl*volume))
		lbins = data[1]
		plotting.plt.step(lbins[0:-1],np.log10(phi),style, label=label, linewidth=lw, where='post')
		if plotdata is True:
			#Hopkins 07 data
			tardat, = np.where(dat.hop_bhlf_obs['redshift']==dat.hop_bhlf_z[zz])
			plotting.plt.errorbar(dat.hop_bhlf_obs['lbol'][tardat] + util.loglbol_sun, dat.hop_bhlf_obs['dphi'][tardat],yerr=dat.hop_bhlf_obs['sig'][tardat],fmt='o',color='grey',ecolor='grey',label='Hopkins+ 2007 (Compilation)')
			if dat.hop_bhlf_z[zz]==6:
				#For z = 6, Barger+ 03 data
				plotting.plt.errorbar([dat.bar_bhlfz6_L],[dat.bar_bhlfz6_phi],xerr=dat.bar_bhlfz6_Lerr,yerr=dat.bar_bhlfz6_phierr,fmt='^',color='k',label='Barger+2003')
				#...and Fiore+ 12
				plotting.plt.errorbar([dat.fio_bhlf_L6F],[dat.fio_bhlf_phi6F],xerr=[[dat.fio_bhlf_L6Fm],[dat.fio_bhlf_L6Fp]],yerr=[[dat.fio_bhlf_errphi6Fm],[dat.fio_bhlf_errphi6Fp]],fmt='s',color='k',label='Fiore+ 2012')
			if dat.hop_bhlf_z[zz]==5:
				l1450 = np.log10(4.4)+util.mcABconv(dat.mcg_bhlf_obs['M1450'],util.c/(0.145e-4))
				dphi = 10**dat.mcg_bhlf_obs['logphi']
				dphip = (2./5.) * (dphi + dat.mcg_bhlf_obs['sig'])
				dphim = (2./5.) * (dphi - dat.mcg_bhlf_obs['sig'])
				dphi = np.log10((2./5.)*dphi)
				dphierr = [dphi-np.log10(dphim),np.log10(dphip)-dphi]
				plotting.plt.errorbar(l1450,dphi,yerr=dphierr,fmt='D',color='k',label='McGreer+ 2013')
		if overplot is False:
			if dotitle is True:
				plotting.plt.title(str(dat.hop_bhlf_zbinsl[zz[0]])+' < z < '+str(dat.hop_bhlf_zbinsh[zz[0]]))
			plotting.plt.xlabel(r'log$_{10}$($L_{bol}$ [ergs/s]))',fontsize=30)
			plotting.plt.ylabel(r'log$_{10}$($\phi$ [Mpc$^{-3}$ dex$^{-1}$])',fontsize=30)
		if label is not None or plotdata is True:
			plotting.plt.legend(loc='lower left',fontsize=20)
		return
