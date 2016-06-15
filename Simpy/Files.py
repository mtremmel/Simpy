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


def mkXML(startstep=None, endstep=None):
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
		os.system('~trq/bin/make_xml.pl')
		os.chdir(wdir)
	return

def create_halo_tree_files(dbsim, add_props=[]):
	import numpy as np
	properties = ['halo_number()','Mvir', 'Vmax', 'Rvir', 'NDM()', 'NStar()', 'NGas()', 'SSC', 'Vcom']
	if len(add_props)>0:
		properties = properties.extend(add_props)

	cnt = 0
	for step in dbsim.timesteps:
		outf = open('out_'+str(cnt)+'.list','w')

		data = step.gather_property(*properties)
		desc = np.ones(len(data[0]))*-1
		n, nn = step.gather_property('halo_number()','later(1).halo_number()')
		desc[(np.in1d(data[0],n))] = nn
		tofile = [data[0],desc,data[1], data[2], np.zeros(len(data[0])), data[3], np.zeros(len(data[0])), data[4]+data[5],+data[6],
				  data[7][:,0], data[7][:,1], data[7][:,2], data[8][:,0], data[8][:,1], data[8][:,2],
				  np.zeros(len(data[0])), np.zeros(len(data[0])), np.zeros(len(data[0])), np.zeros(len(data[0]))]

		np.savetxt(outf, np.column_stack(tofile),
				   fmt=['%d', '%d', '%f', '%f', '%f', '%f', '%f', '%f', '%f', '%f', '%e', '%e', '%e', '%f','%e'])






