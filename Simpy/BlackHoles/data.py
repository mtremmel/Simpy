import numpy as np
#import matplotlib.pyplot as plt
import pynbody
#from pynbody.analysis import profile, angmom, halo
#from pynbody import snapshot, filt, units, config, array
#from scipy import optimize as opt
import warnings
import math
import os
#from pynbody.analysis import pkdgrav_cosmo as cosmo
import pickle
import readcol
import gc
from .. import Files
from .. import util

#orbitcolumns = {'iord':1,'time':2,'step':3,'mass':4,'pos':[5,6,7],'vel':[8,9,10],'mdot':12,'dm':13,'dt':15,'a':16}

def sepOrbitbyStep(file):
	if not os.path.exists('orbitsteps'): os.system('mkdir orbitsteps')
	os.chdir('orbitsteps')
	print "separating orbit file by step..."
	os.system("awk -F ' ' '{if($3!=int($3)) print>(int($3)+1); else print>$3}' "+'../'+file)
	os.chdir('../')
	return

def getOrbitValbyStep(simname,minstep=1,maxstep=4096,clean=True,filename=None,ret_output=False):
        '''
	simname = name of simulation (i.e. simname.004096)
	minstep, maxstep = min, max steps the code will expect to have in the orbit file. 
	clean = True/False, clean up the extra files made from the oritinal orbit file (this is done to spare memory, but they take up lots of space...
        '''
	ea = []
	output = {'iord':[],'time':[],'step':[],'mass':[],'x':[],'y':[],'z':[],'vx':[],'vy':[],'vz':[],'mdot':[],'mdotmean':[],'mdotsig':[],'a':[]}
	for i in range(minstep,maxstep+1):
		if os.path.exists('orbitsteps/'+str(i)):
			print "getting data for step ", i
			bhid,time,step,mass,x,y,z,vx,vy,vz,pot,mdot,dm,E,dt,a = readcol.readcol('orbitsteps/'+str(i),twod=False)
			if clean==True: os.system('rm orbitsteps/'+str(i))
		else: continue
		print "done reading"
		ord = np.argsort(bhid)
		bhid = bhid[ord]
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
		ubhid,uind= np.unique(bhid,return_index=True)	
		print "here"
		for ii in range(len(ubhid)):
			if ii < len(ubhid)-1: idind = range(uind[ii],uind[ii+1])
			else: idind = range(uind[ii],len(step))
			curstep, = np.where(step[idind]==i)
			if len(curstep)==0:continue
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
			utimes, ind = np.unique(time[idind],return_index=True)
			mean,std = util.timeweightedAve(mdot[idind][ind],dt[idind][ind])
			output['mdotmean'].append(mean)
                        output['mdotsig'].append(std)
		print "here 2"
	if filename:
		f = open(filename,'wb')
		pickle.dump(output, f)
		f.close()
	if ret_output:
		return output

def truncOrbitFile(simname,ret_output=False):
	output = getOrbitValbyStep(simname)
	tofile = []
	for key in output:
		tofile.append(output[key])
	if ret_output==False:
		del(output)
		gc.collect()
	tofile = tuple(tofile)
	outputname = simname+'.shortened.orbit'
	print "saving to file..."
	np.savetxt(outputname,np.column_stack(tofile))
	del(tofile)
	gc.collect()
	if ret_output==True:
		return output
	return

def sticthOrbitSteps(simname,nfiles, ret_output=False,overwrite=True, nstart=1):
	ea = np.array([])
	output = {'iord':ea,'time':ea,'step':ea,'mass':ea,'x':ea,'y':ea,'z':ea,'vx':ea,'vy':ea,'vz':ea,'mdot':ea,'mdotmean':ea,'mdotsig':ea,'a':ea}
	outputname = simname+'.shortened.orbit'

	if os.path.exists(outputname) and overwrite==False:
		cols = {'iord':1,'time':2,'step':3,'mass':4,'x':5,'y':6,'z':7,'vx':8,'vy':9,'vz':10,'mdot':11,'mdotmean':12,'mdotsig':13,'a':14}
		data = readcol.readcol(simname+'.shortened.orbit')
		for key in output.keys():
			output[key] = data[:,cols[key]-1]

	for i in range(nfiles):
		num = i+nstart
		f = open(simname+'.shortened.orbit'+str(num),'rb')
		output_part = pickle.load(f)
		f.close()
		for key in output.keys():
			print len(output_part[key]), key
			output[key] = np.append(output[key],output_part[key])
		del(output_part)
		gc.collect()
	tofile = []
        for key in output:
		print key, len(output[key])
                tofile.append(output[key])
	tofile = tuple(tofile)

	print "saving to file..."
        np.savetxt(outputname,np.column_stack(tofile))
	del(tofile)
	gc.collect()
	if ret_output==True:
                return output
	else:
		del(output)
                gc.collect()
		return
