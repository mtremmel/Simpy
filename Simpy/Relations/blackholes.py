import numpy as np
from .. import readcol
from .. import util
from .. import plotting
import pynbody
import os
import gc


def BHMstar(logMstar):
    '''
    BH mass given a stellar mass
    based on analysis by Haring and Rix 2004 and Schramm + Silverman 2013
    '''
    c = 8.31
    a = 1.12
    b = 11.
    predlogMBH = c + a * (logMstar-b)
    return predlogMBH


def BHMBulge(logMbulge):
    '''
    BH mass given a bulge mass
    based on analysis by Haring and Rix 2004 and Kormendy and Ho 2013
    '''
    c = 8.69
    a = 1.16
    b = 11.
    predlogMBH = c + a * (logMbulge-b)
    return predlogMBH


def plt_BHMStar(simname, step, marker='o', size = 100, color='blue', label=None, fit=True, fiterr=True, remove_sats=True):
    from .. import plotting
    import tangos as db
    print "getting data from database..."
    bhdata, Mstar, cen, Rvir, Mvir = db.get_timestep(simname+'/%'+step).gather_property('bh().BH_mass', 'Mstar', 'SSC', 'Rvir', 'Mvir')
    sats = np.zeros(len(Mstar))
    if remove_sats is True:
        for i in range(len(Mstar)):
            D = np.sqrt(np.sum((cen[i] - cen)**2,axis=1))
            close = np.where((D>0)&(D<Rvir)&(Mvir[i]<Mvir))[0]
            if len(close)>0:
                sats[i] = 1
        ok = np.where(sats==0)
        Mstar = Mstar[ok]
        bhdata = bhdata[ok]
    plotting.plt.scatter(Mstar*0.6, bhdata, marker=marker, s=size, color=color, label=label)
    if fit is True:
        lmstar = np.arange(np.log10(Mstar[(Mstar > 0)]).min() - 1., np.log10(Mstar[(Mstar > 0)]).max() + 1., 0.1)
        bhmass = BHMstar(lmstar)
        plotting.plt.plot(10**lmstar, 10**bhmass,'k-', label='Schramm+Silverman 13',lw=2)
        if fiterr is True:
            plotting.plt.fill_between(10**lmstar,10**(bhmass+0.3),10**(bhmass-0.3),color='Grey',alpha=0.5)

    plotting.plt.yscale('log',base=10)
    plotting.plt.xscale('log',base=10)
    plotting.plt.legend(loc='upper left',fontsize=25)
    plotting.plt.ylabel(r"M$_{BH}$ [M$_{\odot}$]")
    plotting.plt.xlabel(r"M$_{\star}$ [M$_{\odot}$]")