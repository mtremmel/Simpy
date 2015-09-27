import gc
import matplotlib.pyplot as plt
from .. import readcol
import pickle
import numpy as np
import pynbody


def getRhoAccHist(inputname, sim, savefile='RhoAcc.pkl', volume=25**3, Mmin=1e6, ):
	mass, mdot, dm, dt, scale = readcol.readcol(inputname,twod=False)
	del(mass)
	del(dm)
	gc.collect()
	o = np.argsort(scale)

	munits = sim.infer_original_units('Msol')
	tunits = sim.infer_original_units('Gyr')
	dM = pynbody.array.SimArray(mdot[o]*dt[o], munits)
	scale = scale[o]
	del(o)
	del(mdot)
	del(dt)
	gc.collect()

	print "summing..."
	rhoBH = np.cumsum(dM)

	if savefile:
		print "saving data..."
		f = open(savefile,'wb')
		pickle.dump([rhoBH,scale],f)
		f.close()
	return rhoBH,scale


def printRhoAccHist(rhoBH, scale, data=True, style='b-', ylog=True, xlog=True, overplot=False, lw=2, label=False):
	shankar09L = 3.2e5
	shankar09H = 5.4e5
	Salvaterra12 = 0.66e4
	Salvaterra12zH = 9
	Salvaterra12zL = 5
	Treister13 = np.array([851.,666.,674.])
	Treister13z = np.array([6.5,7.5,8.5])
	Treister13zErr = np.array([.5,.5,.5])
	Hopkins07zp1, Hopkins07 = readcol.readcol("/nobackupp8/mtremmel/DATA/QSOdata/RhoAccZ.csv",twod=False)
	Hopkins07zp1H, Hopkins07H = readcol.readcol("/nobackupp8/mtremmel/DATA/QSOdata/RhoAccZPLUS.csv",twod=False)
	Hopkins07zp1L, Hopkins07L = readcol.readcol("/nobackupp8/mtremmel/DATA/QSOdata/RhoAccZMINUS.csv",twod=False)
	Hopkins07perr = 10**Hopkins07H - 10**Hopkins07
	Hopkins07merr = 10**Hopkins07 - 10**Hopkins07L
	plt.plot(scale**-1,rhoBH,style,linewidth=lw,label=label)
	if data:
		shankar09 = (shankar09H +shankar09L) / 2.
		err = shankar09H - shankar09
		plt.errorbar([1.03],[shankar09],yerr=[err],color='black',fmt='D',label="Shankar+ 09")
		Salvaterra12z = (Salvaterra12zH + Salvaterra12zL)/2.
		plt.errorbar([Salvaterra12z+1],[Salvaterra12],color='black',fmt='x',xerr=[Salvaterra12zH-Salvaterra12z],yerr=0.5*Salvaterra12,uplims=[True],label='Salvaterra+ 12')
		plt.errorbar(Treister13z,Treister13,color='black',fmt='o',xerr=Treister13zErr,yerr=0.5*Treister13,uplims=[True,True,True], label='Treister+ 13')
		plt.errorbar(Hopkins07zp1,10**Hopkins07,color='grey',fmt='o',yerr=(Hopkins07merr,Hopkins07perr),label='Hopkins+ 07')
	if not overplot:
		if ylog: plt.yscale('log',base=10)
		if xlog: plt.xscale('log',base=10)
		plt.xticks([1,2,3,4,5,6,7,8,9,10],['0','1','2','3','4','5','6','7','8','9'])
		plt.ylabel(r'log($\rho_{acc}$ [M$_{\odot}$ Mpc$^{-3}$])',fontsize=30)
		plt.xlabel('Redshift',fontsize=30)
	return