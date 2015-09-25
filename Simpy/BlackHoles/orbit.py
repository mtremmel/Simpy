import numpy as np
import pynbody

import os
import pickle
import readcol
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


def sepOrbitbyStep(file, min=0, max=1000000000):
	if not os.path.exists('orbitsteps'): os.system('mkdir orbitsteps')
	os.chdir('orbitsteps')
	print "separating orbit file by step..."
	os.system("awk -F ' ' '{if($3!=int($3) && $3 >= " + str(min) + " && $3 <= " + str(
		max) + ") print>(int($3)+1); if($3==int($3) && $3 >= " + str(min) + " && $3 <= " + str(
		max) + ") print>$3}' " + '../' + file)
	os.chdir('../')
	return


def getOrbitValbyStep(simname, minstep=1, maxstep=4096, clean=False, filename=None, ret_output=False):
	'''
	simname = name of simulation (i.e. simname.004096)
	minstep, maxstep = min, max steps the code will expect to have in the orbit file. 
	clean = True/False, clean up the extra files made from the oritinal orbit file (this is done to spare memory, but they take up lots of space...
		'''
	output = {'iord': [], 'time': [], 'step': [], 'mass': [], 'x': [], 'y': [], 'z': [], 'vx': [], 'vy': [], 'vz': [],
			  'mdot': [], 'mdotmean': [], 'mdotsig': [], 'a': []}
	for i in range(minstep, maxstep + 1):
		if os.path.exists('orbitsteps/' + str(i)):
			print "getting data for step ", i
			bhid, time, step, mass, x, y, z, vx, vy, vz, pot, mdot, dm, E, dt, a = readcol.readcol(
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
		a = a[ord]
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
			utimes, ind = np.unique(time[idind], return_index=True)
			mean, std = util.timeweightedAve(mdot[idind][ind], dt[idind][ind])
			output['mdotmean'].append(mean)
			output['mdotsig'].append(std)
	if filename:
		f = open(filename, 'wb')
		pickle.dump(output, f)
		f.close()
	if ret_output:
		return output


def truncOrbitFile(simname, minstep=1, maxstep=4096, ret_output=False):
	output = getOrbitValbyStep(simname, minstep=minstep, maxstep=maxstep)
	outorder = ['iord', 'time', 'step', 'mass', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'mdot', 'mdotmean', 'mdotsig', 'a']
	tofile = []
	for key in outputorder:
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
	@classmethod
	def _get_slice_ind(self, key, timesort=True):
		'''
		return unique values as well as a list of indices of data with each value of self.data[key]
		key = (string) target key to do operation on
		timesort = True/False whether you want to return data sorted by time. Doesn't do anything for single time slices
		'''
		if key != 'time' and key != 'step' and timesort:
			tord = np.argsort(self.data['time'])
		ord_ = np.argsort(self.data[key][tord])
		uvalues, ind = np.unique(self.data[key][tord][ord_], return_index=True)
		slice_ = []
		for i in range(len(uvalues) - 1):
			slice_.append(tord[ord_[ind[i]:ind[i + 1]]])
		slice_.append(ord_[ind[i + 1]:])
		return uvalues, slice_

	def __init__(self, simname, savefile=None):
		ofile = simname + ".shortened.orbit"
		if not os.path.exists(ofile):
			print "ERROR shortened orbit file not found! Exiting..."
			return
		if not os.path.exists('files.list'):
			Files.getFileLists(simname)

		self.data = {}

		# read raw data from shortened orbit file
		print "reading in data. . ."
		self.data['iord'], self.data['time'], self.data['step'], self.data['mass'], self.data['x'], self.data['y'], \
		self.data['z'], self.data['vx'], self.data['vy'], self.data['vz'], self.data['mdot'], self.data['mdotmean'], \
		self.data['mdotsig'], self.data['scalefac'] = readcol.readcol(ofile, twod=False)

		# get information on iord,step data for easy future data recovery
		print "reducing data. . ."
		self.bhiords, self.id_slice = self._get_slice_ind('iord')
		self.steps, self.step_slice = self._get_slice_ind('step')
		self.times = np.unique(self.data['time'])

		if savefile:
			f = open(savefile, 'wb')
			pickle.dump(self, f)
			f.close()

	def single_BH_data(self, iord, key):
		o, = np.where(self.bhiords == iord)
		slice_ = self.id_slice[o[0]]
		return self.data[key][slice_]

	def single_step_data(self, iord, key):
		o, = np.where(self.bhiords == iord)
		slice_ = self.step_slice[o[0]]
		return self.data[key][slice_]
