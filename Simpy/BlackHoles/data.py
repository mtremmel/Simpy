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

def sepOrbitbyStep(file,min=0,max=1000000000):
	if not os.path.exists('orbitsteps'): os.system('mkdir orbitsteps')
	os.chdir('orbitsteps')
	print "separating orbit file by step..."
	os.system("awk -F ' ' '{if($3!=int($3) && $3 >= "+str(min)+" && $3 <= "+str(max)+") print>(int($3)+1); if($3==int($3) && $3 >= "+str(min)+" && $3 <= "+str(max)+") print>$3}' "+'../'+file)
	os.chdir('../')
	return

def getOrbitValbyStep(simname,minstep=1,maxstep=4096,clean=True,filename=None,ret_output=False):
        '''
	simname = name of simulation (i.e. simname.004096)
	minstep, maxstep = min, max steps the code will expect to have in the orbit file. 
	clean = True/False, clean up the extra files made from the oritinal orbit file (this is done to spare memory, but they take up lots of space...
        '''
	ea = np.array([])
	output = {'iord':ea,'time':ea,'step':ea,'mass':ea,'x':ea,'y':ea,'z':ea,'vx':ea,'vy':ea,'vz':ea,'mdot':ea,'mdotmean':ea,'mdotsig':ea,'a':ea}
	for i in range(minstep,maxstep+1):
		if os.path.exists('orbitsteps/'+str(i)):
			print "getting data for step ", i
			bhid,time,step,mass,x,y,z,vx,vy,vz,pot,mdot,dm,E,dt,a = readcol.readcol('orbitsteps/'+str(i),twod=False)
			if clean==True: os.system('rm orbitsteps/'+str(i))
		else: continue
		del(E)
		del(dm)
		del(pot)
		gc.collect()
		curstep, = np.where(step==i)
		if len(curstep)==0: continue
		ubhid,cnt = np.unique(bhid[curstep],return_index=True) #only add BHs that are actually present at the current timestep
		output['iord'] = np.append(output['iord'],ubhid)
		output['step'] = np.append(output['step'],step[curstep][cnt])
		output['mass'] = np.append(output['mass'],mass[curstep][cnt])
		output['time'] = np.append(output['time'],time[curstep][cnt])
		output['x'] = np.append(output['x'],x[curstep][cnt])
		output['y'] = np.append(output['y'],y[curstep][cnt])
		output['z'] = np.append(output['z'],z[curstep][cnt])
		output['vx'] = np.append(output['vx'],vx[curstep][cnt])
                output['vy'] = np.append(output['vy'],vy[curstep][cnt])
                output['vz'] = np.append(output['vz'],vz[curstep][cnt])
		output['mdot'] = np.append(output['mdot'],mdot[curstep][cnt])
		output['a'] = np.append(output['a'],a[curstep][cnt])
		#nbh = len(ubhid)
		#cc = 0
		for id in ubhid:
			#cc+=1
			#if cc%100==0: print float(cc)*100/nbh, " % done getting ave mdot properties"
			curbh, = np.where(bhid==id)
			utimes, ind = np.unique(time[curbh],return_index=True)
			mean,std = util.timeweightedAve(mdot[curbh[ind]],dt[curbh[ind]])	
			output['mdotmean'] = np.append(output['mdotmean'],mean)
			output['mdotsig'] = np.append(output['mdotsig'],std)
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
