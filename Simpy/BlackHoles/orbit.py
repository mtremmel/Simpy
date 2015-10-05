import numpy as np
import pynbody

import os
import pickle
from .. import readcol
import gc
from .. import Files
from .. import util

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
	Files.cklists(simname,NCHILADA=NCHILADA)
	f = open('files.list','r')
	sim = pynbody.load(f.readline().strip('\n'))
	f.close()
	munits = sim.infer_original_units('Msol')
	MBHinit = MBHinit/float(munits)
	if not os.path.exists('orbitsteps'): os.system('mkdir orbitsteps')
	os.chdir('orbitsteps')
	print "separating orbit file by step..."
	os.system("awk -F ' ' '{if($4 - $13 > "+str(MBHinit)+" && $3!=int($3) && $3 >= " + str(minstep) + " && $3 <= " + str(maxstep) + ") print>(int($3)+1); if($4 - $13 > "+str(MBHinit)+" && $3==int($3) && $3 >= " + str(minstep) + " && $3 <= " + str(maxstep) + ") print>$3}' " + '../' + simname+'.orbit')
	os.chdir('../')
	return


def getOrbitValbyStep(minstep=1, maxstep=4096, clean=False, filename=None, ret_output=False):
	'''
	simname = name of simulation (i.e. simname.004096)
	minstep, maxstep = min, max steps the code will expect to have in the orbit file. 
	clean = True/False, clean up the extra files made from the oritinal orbit file (this is done to spare memory, but they take up lots of space...
		'''
	output = {'iord': [], 'time': [], 'step': [], 'mass': [], 'x': [], 'y': [], 'z': [], 'vx': [], 'vy': [], 'vz': [],
			  'mdot': [], 'mdotmean': [], 'mdotsig': [], 'a': []}

	if not os.path.exists('orbitsteps/'):
		print "ERROR! can't find orbitsteps... run sepOrbitbyStep first!"
		return
	for i in range(minstep, maxstep + 1):
		if os.path.exists('orbitsteps/' + str(i)):
			print "getting data for step ", i
			bhid, time, step, mass, x, y, z, vx, vy, vz, pot, mdot, dm, E, dt, a = readcol.readcol('orbitsteps/' + str(i), twod=False)
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
		a = a[ord]
		#bad, = np.where((mass - mdot*dt < MBHinit)|(mass<MBHinit))
		#mdot[bad] = 0
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
#			utimes, ind = np.unique(time[idind], return_index=True)
			mean, std = util.timeweightedAve(mdot[idind], dt[idind])
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
	outorder = ['iord', 'time', 'step', 'mass', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'mdot', 'mdotmean', 'mdotsig', 'a']
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
			   fmt=['%d', '%f', '%d', '%e', '%f', '%f', '%f', '%f', '%f', '%f', '%e', '%e', '%e', '%f'])
	del (tofile)
	gc.collect()
	if ret_output:
		return output
	return


def sticthOrbitSteps(simname, nfiles, ret_output=False, overwrite=True, nstart=1):
	ea = np.array([])
	output = {'iord': ea, 'time': ea, 'step': ea, 'mass': ea, 'x': ea, 'y': ea, 'z': ea, 'vx': ea, 'vy': ea, 'vz': ea,
			  'mdot': ea, 'mdotmean': ea, 'mdotsig': ea, 'a': ea}
	outorder = ['iord', 'time', 'step', 'mass', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'mdot', 'mdotmean', 'mdotsig', 'a']
	outputname = simname + '.shortened.orbit'

	if os.path.exists(outputname) and overwrite == False:
		cols = {'iord': 1, 'time': 2, 'step': 3, 'mass': 4, 'x': 5, 'y': 6, 'z': 7, 'vx': 8, 'vy': 9, 'vz': 10,
				'mdot': 11, 'mdotmean': 12, 'mdotsig': 13, 'a': 14}
		data = readcol.readcol(simname + '.shortened.orbit')
		for key in output.keys():
			output[key] = data[:, cols[key] - 1]

	for i in range(nfiles):
		num = i + nstart
		f = open(simname + '.shortened.orbit' + str(num), 'rb')
		output_part = pickle.load(f)
		f.close()
		for key in output.keys():
			output[key] = np.append(output[key], output_part[key])
		del (output_part)
		gc.collect()
	tofile = []
	for key in outorder:
		tofile.append(output[key])
	tofile = tuple(tofile)

	print "saving to file..."
	np.savetxt(outputname, np.column_stack(tofile),
			   fmt=['%d', '%f', '%d', '%e', '%f', '%f', '%f', '%f', '%f', '%f', '%e', '%e', '%e', '%f'])
	del (tofile)
	gc.collect()
	if ret_output == True:
		return output
	else:
		del (output)
		gc.collect()
		return


class Orbit(object):
	def __init__(self, simname, savefile=None, NCHILADA=True):
		self.simname = simname
		self.nchil = NCHILADA
		ofile = simname + ".shortened.orbit"
		if not os.path.exists(ofile):
			print "ERROR shortened orbit file not found! Exiting..."
			return
		Files.cklists(simname, NCHILADA=NCHILADA)

		# read raw data from shortened orbit file
		print "reading in data. . ."
		bhid, time, step, mass, x, y, z, vx, vy, vz, mdot, mdotmean, mdotsig, scalefac = readcol.readcol(ofile, twod=False,nanval=0.0)
		self.data = {'iord':bhid, 'time':time, 'step':step, 'mass':mass, 'x':x, 'y':y, 'z':z, 'vx':vx, 'vy':vy, 'vz':vz, 'mdot':mdot, 'mdotmean':mdotmean, 'mdotsig':mdotsig, 'scalefac':scalefac}
		del(bhid, time, step, mass, x, y, z, vx, vy, vz, mdot, mdotmean, mdotsig, scalefac)
		gc.collect()

		#make each data object a simulation array with correct units
		defunits = {'iord':None, 'time':'Gyr', 'step':None, 'mass':'Msol', 'mdot':'Msol yr**-1', 'mdotmean':'Msol yr**-1', 'mdotsig':'Msol yr**-1', 'x':'kpc', 'y':'kpc', 'z':'kpc', 'vx':'km s**-1', 'vy':'km s**-1', 'vz':'km s**-1','scalefac':None}
		f = open('files.list','r')
		sim = f.readline()
		s = pynbody.load(sim.strip('\n'))
		f.close()

		#control for "fake" bhs caused by restarts from outputs
		print "checking for fake BHs. . ."
		sl = pynbody.tipsy.StarLog(self.simname+'.starlog')
		slbhiords = sl['iord'][(sl['tform']<0)]
		ok, = np.where(np.in1d(self.data['iord'],slbhiords))
		for key in self.data.keys():
			self.data[key] = self.data[key][ok]

		#convert comoving quantities to physical units
		for key in ['x','y','z','vx','vy','vz']: self.data[key] *= self.data['scalefac']
		for key in self.data.keys():
			unit = None
			if defunits[key] is not None:
				unit = s.infer_original_units(defunits[key])
				if key in ['x','y','z','vx','vy','vz']: unit = unit/pynbody.units.Unit('a')
			self.data[key] = pynbody.array.SimArray(self.data[key],unit)
			if defunits[key] is not None:
				print unit, defunits[key]
				self.data[key].convert_units(defunits[key])

		# get information on iord,step data for easy future data recovery
		print "reducing data. . ."
		self.bhiords, self.id_slice = self._get_slice_ind('iord')
		self.steps, self.step_slice = self._get_slice_ind('step')
		self.times = np.unique(self.data['time'])
		self._calc_lum(er=0.1)
		self.get_all_BH_tform(sl)

		if savefile:
			f = open(savefile, 'wb')
			pickle.dump(self, f)
			f.close()
		return


	def _calc_lum(self, er=0.1):
		csq = pynbody.array.SimArray((2.998e10)**2,'erg g**-1')
		self.data['lum'] = self.data['mdot'].in_units('g s**-1')*csq*er
		return



	def _get_slice_ind(self, key,orderby='time'):
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


	def single_BH_data(self, iord, key):
		o, = np.where(self.bhiords == iord)
		slice_ = self.id_slice[o[0]]
		return self.data[key][slice_]


	def single_step_data(self, step, key):
		o, = np.where(self.steps == step)
		slice_ = self.step_slice[o[0]]
		return self.data[key][slice_]


	def get_all_BH_tform(self,sl):
		#sl = pynbody.tipsy.StarLog(self.simname+'.starlog')
		ord = np.argsort(sl['iord'])
		bhind, = np.where(np.in1d(sl['iord'][ord],self.bhiords))
		self.tform = sl['tform'][ord][bhind] * -1
		if self.tform.min() < 0: print "WARNING! Positive tforms were found for BHs!"
		del sl
		gc.collect()


	def get_single_BH_dm(self, iord):
		if hasattr(self,'tform') is False:
			self.get_all_BH_tform()
		time = self.single_BH_data(iord,'time')
		mdot = self.single_BH_data(iord,'mdotmean')
		time.convert_units('yr')
		timel = np.append(self.tform[(self.bhiords==iord)].in_units('yr'), time[:-1])
		assert timel[0] < time[0]
		dt = pynbody.array.SimArray(time - timel, 'yr')
		return dt*mdot.in_units('Msol yr**-1')


	def get_acc_hist(self):
		self.acchist = {'dm':np.array([]), 'time':np.array([]), 'mass':np.array([]), 'scalefac':np.array([]), 'lum':np.array([])}
		for iord in self.bhiords:
			dmpart = self.get_single_BH_dm(iord)
			mpart = self.single_BH_data(iord,'mass')
			lumpart = self.single_BH_data(iord,'lum')
			tpart = self.single_BH_data(iord,'time')
			scalepart = self.single_BH_data(iord,'scalefac')
			np.append(self.acchist['dm'],dmpart.in_units('Msol'))
			np.append(self.acchist['mass'], mpart.in_units('Msol'))
			np.append(self.acchist['lum'],lumpart.in_units('ergs s**-1'))
			np.append(self.acchist['time'], tpart.in_units('Gyr'))
			np.append(self.acchist['scalefac'], scalepart)

		ord = np.sort(self.acchist['scalefac'])
		for key in self.acchist.keys():
			self.acchist[key] = self.acchist[key][ord]


	def getprogbhs(self):
		time, step, ID, IDeat, ratio, kick = readcol.readcol(self.simname+'.mergers',twod=False)
		initlist = [[] for i in range(len(self.bhiords))]
		self.prog = {'iord':[[] for i in range(len(self.bhiords))],'kick':[[] for i in range(len(self.bhiords))], 'ratio':[[] for i in range(len(self.bhiords))], 'step':[[] for i in range(len(self.bhiords))], 'time':[[] for i in range(len(self.bhiords))]}
		o = np.argsort(ID)
		ID = ID[o]
		step = step[o]
		IDeat = IDeat[o]
		time = time[o]
		ratio = ratio[o]
		kick = kick[o]
		uID, ind = np.unique(ID,return_index=True)
		eaterind, = np.where(np.in1d(self.bhiords, uID))
		for i in range(len(uID)-1):
			self.prog['iord'][eaterind[i]].extend(IDeat[ind[i]:ind[i+1]])
			self.prog['step'][eaterind[i]].extend(step[ind[i]:ind[i+1]])
			self.prog['time'][eaterind[i]].extend(time[ind[i]:ind[i+1]])
			self.prog['ratio'][eaterind[i]].extend(ratio[ind[i]:ind[i+1]])
			if kick.max() > 0:
				self.prog['kick'][eaterind[i]].extend(kick[ind[i]:ind[i+1]])


	def gethalos(self):
		f = open('files.list','r')
		simfiles = f.readlines()
		nsnaps = len(simfiles)
		f.close()
		f = open('steps','r')
		snaps = f.readlines()
		f.close()
		initarr = np.ones((len(self.bhiords),len(snaps))) * -1
		self.bhhalos = {'Grp':initarr}
		for i in range(nsnaps):
			if os.path.exists(simfiles[i].strip('\n')+'.rockstar.grp'):
				grpfile = simfiles[i].strip('\n')+'.rockstar.grp'
			else:
				if os.path.exists(simfiles[i].strip('\n')+'.amiga.grp'):
					grpfile = simfiles[i].strip('\n')+'.amiga.grp'
				else:
					print "ERROR there is no grp file for this step!"
					continue
			sim = pynbody.load(simfiles[i].strip('\n'))
			simind = np.where(np.in1d(sim.stars['iord'],self.bhiords))
			orbind = np.where(np.inqd(self.bhiords,sim.stars['iord']))
			simind += len(sim.dm)+len(sim.gas)
			del sim['iord']
			gc.collect()
			grp = readcol.readcol(grpfile,skipline=1)
			self.bhhalos['Grp'][orbind][:,i] = grp[simind]
			del grp, simind, orbind
			gc.collect()


	def plt_single_BH_data(self, iord, keyx, unitx, keyy, unity, style, lw=1, msize=10, ylog=True, xlog=False, label=None):
		from .. import plotting

		ydat = self.single_BH_data(iord,keyy)
		xdat = self.single_BH_data(iord,keyx)
		plotting.plt.plot(xdat.in_units(unitx), ydat.in_units(unity), style, label=label, linewidth=lw, markersize=msize)
		if xlog:
			plotting.plt.xscale('log', base=10)
		if ylog:
			plotting.plt.yscale('log', base=10)
		if label:
			plotting.plt.legend()
		return
