import pynbody
import numpy as np
import os
from .. import readcol

def writebhmark(simname, step, Name=None, iord=False, massrange=False):
    if not Name:
        f = open('BH.' + step + '.mark', 'w')
    else:
        f = open(Name, 'w')
    s = pynbody.load(simname + '.' + step)
    f.write(str(len(s)) + ' ' + str(len(s.gas)) + ' ' + str(len(s.star)) + '\n')
    if not iord:
        if not massrange:
            bhind, = np.where(s.stars['tform'] < 0)
        else:
            if len(massrange) != 2:
                print "error massrange must be a length 2 tuple!"
                return
            bhind, = np.where((s.stars['tform'] < 0) & (s.stars['mass'].in_units('Msol') < massrange[1]) & (
            s.stars['mass'].in_units('Msol') > massrange[0]))
    else:
        bhind = np.array([])
        for ii in range(len(iord)):
            tmpind, = np.where(s.stars['iord'] == iord[ii])
            if len(tmpind) == 0: print "uh oh... iord ", iord[ii], " not found!"
            bhind = np.append(bhind, tmpind)
    bhindreal = bhind + len(s.dark) + len(s.gas) + 1
    for ii in range(len(bhindreal)):
        f.write(str(bhindreal[ii]) + '\n')
    f.close()
    del (s)
    return


def getBHiords(simname):
    if not os.path.exists("BHid.list"):
        print "finding IDs for all BHs that ever existed..."
        os.system("awk '{print $1}' " + simname + ".orbit > BHid.list")
        f = open("BHid.list", 'r')
        id = f.readlines()
        id = np.array(id)
        id = id.astype('int')
        id = np.unique(id)
        f.close()
        os.system("rm BHid.list")
        np.savetxt("BHid.list", id)
    else:
        print "previous BHid.list file found! reading it..."
        id, = readcol.readcol("BHid.list", twod=False)

    return id

def get_new_masses(bhid_list, time_list, mass_list, bhmergers, orig_seed = 1e6, new_seed = 1e3, useonly=None):
    from . import mergers
    new_masses = np.copy(mass_list)
    cnt = 0
    for id, t, m in zip(bhid_list,time_list,mass_list):
        print id, t, m
        massdiff = 0
        prog_list, masses, times = mergers.get_complete_prog_list(bhmergers, id, t, useonly=useonly, return_details=True)
        small = np.where((masses < orig_seed)&(masses > new_seed))
        big = np.where(masses >= orig_seed)
        if len(small > 0):
            massdiff += np.sum(masses[small] - new_seed)
        if len(big) > 0:
            massdiff += len(big) * (orig_seed - new_seed)
        new_masses[cnt] -= massdiff
        if m >= orig_seed:
            new_masses[cnt] -= (orig_seed - new_seed)
        if m < orig_seed and m > new_seed:
            new_masses[cnt] -= (m - new_seed)
        cnt += 1
    return new_masses

