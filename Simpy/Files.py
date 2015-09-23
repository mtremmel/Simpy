import os


def getFileLists(parentdir, simname, NCHILADA=True):
	totalpath = parentdir + simname
	simname_split = totalpath.split('.')
	num = len(simname_split)
	if not NCHILADA:
		os.system('ls ' + parentdir + '*.iord | cut -d"." -f1-' + str(num + 1) + ' > files.list')
		os.system('ls ' + parentdir + '*.iord | cut -d"." -f' + str(num + 1) + ' > steps.list')
	else:
		os.system('ls ' + parentdir + simname+'.00* --directory > files.list')
		os.system('ls ' + parentdir + simname + '.00* --directory | cut -d"." -f' + str(num + 1) +' > steps.list')
	##os.system('ls ' + parentdir + simname +'.00*.grp | cut -d "." -f1-' + str(num + 1) + '> grpfiles.list' )
	return


def listfilebytype(parentdir, simname, ending):
	simname_split = totalpath.split('.')
	num = len(simname_split)
	os.system('ls ' + parentdir + simname +'.00*.'+ending+' | cut -d "." -f1-' + str(num + 1) + '> grpfiles.list')
	return


def cklists(parentdir,simname,NCHILADA=True):
	if not os.path.exists('files.list'):
		print "files,steps lists not found, running getFileLists..."
		getFileLists(parentdir, simname, NCHILADA=NCHILADA)
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
	cklists(parentdir, simname, NCHILADA=True)
	f = open('files.list', 'r')
	files = f.readlines()
	f.close()
	f = open('steps.list', 'r')
	steps = f.readlines()
	f.close()
	for i in range(len(files)):
		os.system('ln -s ' + files[i].strip('\n') + ' ' + simname + '.' + steps[i].strip('\n'))
	os.system('ln -s ' + pathdir + simname + '.orbit ' + simname + '.orbit')
	os.system('ln -s ' + pathdir + simname + '.starlog ' + simname + '.starlog')
	os.system('ln -s ' + pathdir + simname + '.BHAccLog ' + simname + '.BHAccLog')
	return


def mkXML(do='All',start=0):
	wdir = os.getcwd()
	cklists(parentdir, simname, NCHILADA=True)
	f = open('files.list','r')
	files = f.readlines()
	if do == 'All': num = len(files)
	else: num = do
	for i in range(num)+start:
		os.chdir(files[i].strip('\n'))
		os.system('~trq/bin/make_xml.pl')
	os.chdir(wdir)
	return

