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


def plt_BHMStar(simname, step, marker='o', size = 100, color='blue', label=None, fit=True, fiterr=True):
    from .. import plotting
    import tangos as db
    print "getting data from database..."
    bhdata, Mstar = db.get_timestep(simname+'/%'+step).gather_property('bh().BH_mass', 'Mstar')
    plotting.plt.scatter(Mstar, bhdata, marker=marker, s=size, color=color, label=label)
    if fit is True:
        lmstar = np.arange(np.log10(Mstar[(Mstar > 0)]).min() - 1., np.log10(Mstar[(Mstar > 0)]).max() + 1., 0.1)
        bhmass = BHMstar(lmstar)
        plotting.plt.plot(10**lmstar, 10**bhmass,'k-', label='Schramm+Silverman 13')
        if fiterr is True:
            plotting.plt.fill_between(10**lmstar,10**(bhmass+0.3),10**(bhmass-0.3),color='Grey',alpha=0.5)

    plotting.plt.yscale('log',base=10)
    plotting.plt.xscale('log',base=10)
    plt.legend(loc='upper left')