import os


def getFileLists(simname, NCHILADA=True):
	simname_split = simname.split('.')
	num = len(simname_split)
	os.system('ls ' + simname + '.00???? --directory > files.list')
	os.system('ls ' + simname + '.00???? --directory | cut -d"." -f' + str(num + 1) + ' > steps.list')

#	if not NCHILADA:
#		os.system('ls *.iord | cut -d"." -f1-' + str(num + 1) + ' > files.list')
#		os.system('ls *.iord | cut -d"." -f' + str(num + 1) + ' > steps.list')
#	else:
#		os.system('ls ' + simname + '.00???? --directory > files.list')
#		os.system('ls ' + simname + '.00???? --directory | cut -d"." -f' + str(num + 1) +' > steps.list')
	##os.system('ls ' + parentdir + simname +'.00*.grp | cut -d "." -f1-' + str(num + 1) + '> grpfiles.list' )
	return


def listfilebytype(simname, ending):
	simname_split = simname.split('.')
	num = len(simname_split)
	os.system('ls ' + simname +'.00*.'+ending+' | cut -d "." -f1-' + str(num + 1) + '> grpfiles.list')
	return


def cklists(simname,NCHILADA=True):
	if not os.path.exists('files.list'):
		print "files,steps lists not found, running getFileLists..."
		getFileLists(simname, NCHILADA=NCHILADA)
	return


def mkAbridge(filename,columns=[1,2,3,4],condition=None,suffix=None):
	if suffix==None:
		outfile = filename+'.abridged'
	else:
		outfile = filename + suffix
	if not os.path.exists(outfile):
		print "making abridged file for ", filename
		cstr = 'awk {'
		if condition:
			cstr += condition
		cstr += ' print '
		for c in columns:
			cstr = cstr + ' $'+cstr(c)
		cstr += '} '
		cstr += filename
		cstr += ' > '
		cstr =+ outfile
		print cstr
		os.system(cstr)
	return


def makelinks(parentdir, simname):
	cklists(parentdir+simname, NCHILADA=True)
	f = open('files.list', 'r')
	files = f.readlines()
	f.close()
	f = open('steps.list', 'r')
	steps = f.readlines()
	f.close()
	for i in range(len(files)):
		if not os.path.exists(simname + '.' + steps[i].strip('\n')):
			os.system('ln -s ' + files[i].strip('\n') + ' ' + simname + '.' + steps[i].strip('\n'))
		else:
			print "link already found. skipping. . ."
			continue
	return


def mkXML(startstep=None, endstep=None, path='~trq/bin/'):
	wdir = os.getcwd()
	f = open('files.list','r')
	files = f.readlines()
	f.close()
	f = open('steps.list', 'r')
	steps = f.readlines()
	f.close()
	for i in range(len(files)):
		stint = int(steps[i].strip('\n'))
		print stint
		if startstep and stint < startstep: continue
		if endstep and stint > endstep: break
		os.chdir(files[i].strip('\n'))
		if os.path.exists('description.xml'):
			print "xml file found!"
			os.chdir(wdir)
			continue
		os.system(path+'/make_xml.pl')
		os.chdir(wdir)
	return

def run_NCHIL_to_TIPSY(target_dir='TipsyFiles', exec_path='~/utility/TreeDataFormat/', overwrite=False, nstart=0, nend=-1):
	if not os.path.exists(target_dir):
		os.system('mkdir ' + target_dir)
	f = open('files.list')
	lines = f.readlines()
	ntodo = len(lines)
	if nend>0:
		ntodo = nend-nstart
	f.close()
	cnt = 0
	for ll in lines:
		if cnt < nstart:
			cnt +=1
			continue
		if cnt > nend:
			break
		print ll
		print cnt/float(ntodo)*100, '% done'
		if os.path.exists(target_dir+'/'+ll.strip('\n')):
			if overwrite is False:
				print "file already found! skipping"
				cnt += 1
				continue
			else:
				print "file already found! overwrite flag is on, so proceeding anyway"
		os.system(exec_path+'salsa2tipsy '+ll.strip('\n') + ' ' + target_dir+'/'+ll.strip('\n'))
		cnt += 1

def create_halo_tree_files(dbsim, h=0.6776931508813172, nmax=None):
	import numpy as np
	cnt = 0
	desc = open('DescScales.txt','w')
	desc.close()
	import pynbody
	for step in dbsim.timesteps:
		halodat = step.gather_property('halo_number()', 'Mvir', 'vmax_dm_global()', 'Rvir', 'NDM()', 'NStar()', 'NGas()', 'SSC', 'Vcom')
		halo_linkdat = step.gather_property('halo_number()','vmax_dm_global()','later(1).halo_number()','later(1).vmax_dm_global()')

		desc = np.ones(len(halodat[0]))*-1
		desc[(np.in1d(halodat[0],halo_linkdat[0]))] = halo_linkdat[2]
		ID = halodat[0]
		Mvir = halodat[1] * h
		Vmax = halodat[2]
		Rvir = halodat[3] * h
		Np = halodat[4]+halodat[5]+halodat[6]
		Xc = pynbody.array.SimArray(halodat[7][:,0], 'kpc')
		Yc = pynbody.array.SimArray(halodat[7][:,1], 'kpc')
		Zc = pynbody.array.SimArray(halodat[7][:,2], 'kpc')
		VXc = halodat[8][:,0]
		VYc = halodat[8][:,1]
		VZc = halodat[8][:,2]
		Jx = np.zeros(len(halodat[0]))
		Jy = np.zeros(len(halodat[0]))
		Jz = np.zeros(len(halodat[0]))
		Spin = np.zeros(len(halodat[0]))
		Rs = np.zeros(len(halodat[0]))
		Vrms = np.zeros(len(halodat[0]))

		outf = open('out_'+str(cnt)+'.list','w')

		if nmax:
			ok = np.where((ID<nmax)&(desc<nmax))[0]

			tofile = [ID[ok], desc[ok], Mvir[ok], Vmax[ok], Vrms[ok], Rvir[ok], Rs[ok], Np[ok],
				  Xc.in_units('Mpc')[ok]*h, Yc.in_units('Mpc')[ok] * h, Zc.in_units('Mpc')[ok] * h,
				  VXc[ok], VYc[ok], VZc[ok], Jx[ok], Jy[ok], Jz[ok], Spin[ok]]

		else:
			tofile = [ID, desc, Mvir, Vmax, Vrms, Rvir, Rs, Np,
				  Xc.in_units('Mpc')*h, Yc.in_units('Mpc') * h, Zc.in_units('Mpc') * h,
				  VXc, VYc, VZc, Jx, Jy, Jz, Spin]

		np.savetxt(outf, np.column_stack(tofile),
				   header="#ID DescID Mvir Vmax Vrms Rvir Rs Np X Y Z VX VY VZ JX JY JZ Spin",
				   fmt=['%d', '%d', '%e', '%f', '%f', '%f', '%f', '%d', '%f', '%f', '%f', '%f', '%f', '%f','%e', '%e','%e','%f'])
		outf.close()

		desc = open('DescScales.txt','a')
		desc.write(str(cnt)+' '+str(1./(1+step.redshift))+'\n')
		desc.close()

		cnt += 1






