import os


def getFileLists(simname, NCHILADA=True):
	simname_split = simname.split('.')
	num = len(simname_split)
	if not NCHILADA:
		os.system('ls *.iord | cut -d"." -f1-' + str(num + 1) + ' > files.list')
		os.system('ls *.iord | cut -d"." -f' + str(num + 1) + ' > steps.list')
	else:
		os.system('ls ' + simname + '.00* --directory > files.list')
		os.system('ls ' + simname + '.00* --directory | cut -d"." -f' + str(num + 1) +' > steps.list')
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
	if sufffix==None:
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
			sos.system('ln -s ' + files[i].strip('\n') + ' ' + simname + '.' + steps[i].strip('\n'))
		else:
			print "link already found. skipping. . ."
			continue
	os.system('ln -s ' + pathdir + simname + '.orbit ' + simname + '.orbit')
	os.system('ln -s ' + pathdir + simname + '.starlog ' + simname + '.starlog')
	os.system('ln -s ' + pathdir + simname + '.BHAccLog ' + simname + '.BHAccLog')
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
		if startstep and stint < statstep: continue
		if endstep and stint > endstep: break
		os.chdir(files[i].strip('\n'))
		os.system('~trq/bin/make_xml.pl')
	os.chdir(wdir)
	return

