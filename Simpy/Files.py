import os


def getFileLists(simname):
        simname_split = simname.split('.')
        num = len(simname_split)
        os.system('ls  *.iord | cut -d"." -f1-'+str(num+1)+' > files.list')
        os.system('ls *.iord | cut -d"." -f'+str(num+1)+' > steps.list')
	os.system('ls '+simname+'.00*.grp | cut -d "." -f1-'+str(num+1)+ '> grpfiles.list' )
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
