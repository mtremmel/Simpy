from .. import Files, readcol, dbutil
import numpy as np
import pynbody
import os


def calc_lum(mdot,er=0.1):
        csq = pynbody.array.SimArray((2.998e10) ** 2, 'erg g**-1')
        return pynbody.array.SimArray(mdot,'Msol yr**-1').in_units('g s**-1') * csq * er


class StepList(object):
    def __init__(self, steplist, db, boxsize, simname):
        self._steplist = steplist
        self.data = {}
        for step in steplist:
            print "gathering data for step ", step
            dbstep = db.get_timestep(simname+'/%'+step)
            self.data[step] = StepData(step,dbstep,boxsize)

    def __getitem__(self,step):
        if int(step) not in self._steplist.astype(np.int):
            raise KeyError
        else:
            target = np.where(self._steplist.astype(np.int)==int(step))[0]
        return self.data[self._steplist[target]]


class StepData(object):
    def __init__(self, step, dbstep, boxsize):
        self.redshift = dbstep.redshift
        self.time = dbstep.time_gyr
        self.a = 1./(1.+self.redshift)
        self.id = step
        self.boxsize = boxsize

        print "Gathering BH data..."

        bhids, bhmass, bhmdot, offset, dist, hostnum = \
                dbstep.gather_property('N', 'BH_mass', 'BH_mdot_ave', 'BH_central_offset', 'BH_central_distance','host')

        nbhs = len(bhids)

        if nbhs==0:
            print "No BHs Found in This Step"
            self.bh = {'lum':[], 'dist':[], 'pos':[], 'halo':[], 'mass':[],
                   'bhid':[], 'nearhalo':[],'mdot':[], 'neardist':[]}
            self.halo_properties = {'Mvir':[], 'Mstar':[], 'Rvir':[], 'Mgas':[], 'SSC':[], 'N':[]}
            return


        self.bh = {'lum':calc_lum(bhmdot), 'dist':dist, 'pos':offset, 'host':hostnum, 'mass':bhmass,
                   'bhid':bhids, 'mdot':bhmdot}

        print "slicing data..."
        self.host_ids, self._halo_slices,self._host_indices = self._get_halo_slices()

        nhalos = len(self.host_ids)

        print "Gathering halo data"
        Mvir, Mstar, Rvir, Mgas, SSC, hid = dbstep.gather_property('Mvir', 'Mstar', 'Rvir', 'Mgas', 'SSC', 'N')

        self.halo_properties = {'Mvir':np.zeros(nbhs), 'Mstar':np.zeros(nbhs), 'Rvir':np.zeros(nbhs), 'Mgas':np.zeros(nbhs),
                                'SSC':np.zeros((nbhs,3)), 'N':np.zeros(nbhs)}

        print "matching halo data to BHs"
        for i in range(nhalos):
            target = np.where(hid==self.host_ids[i])[0]
            self.halo_properties['Mvir'][self._halo_slices[i]] = Mvir[target]
            self.halo_properties['Mstar'][self._halo_slices[i]] = Mstar[target]
            self.halo_properties['Mgas'][self._halo_slices[i]] = Mgas[target]
            self.halo_properties['Rvir'][self._halo_slices[i]] = Rvir[target]
            self.halo_properties['N'][self._halo_slices[i]] = hid[target]
            self.halo_properties['SSC'][self._halo_slices[i]] = SSC[target]


        print "finding nearby halos"
        self.nearby_halo_properties = \
            {'Mvir':np.zeros(nbhs), 'Mstar':np.zeros(nbhs), 'Rvir':np.zeros(nbhs),
             'Mgas':np.zeros(nbhs), 'SSC':np.zeros(nbhs), 'N':np.zeros(nbhs)}
        self.bh['nearhalo'] = np.zeros(nbhs)
        self.bh['nearpos'] = np.zeros((nbhs,3))

        hpos = self.host_prop('SSC')
        hN = self.host_prop('N')

        for i in range(nhalos):
            relpos = SSC - hpos[i]
            bad = np.where(np.abs(relpos) > self.boxsize.in_units('kpc', a=self.a)/2.)
            relpos[bad] = -1.0 * (relpos[bad]/np.abs(relpos[bad])) * \
                              (self.boxsize.in_units('kpc', a=self.a) - np.abs(relpos[bad]))
            reldist = np.sqrt(np.sum(relpos**2, axis=1))

            ok, = np.where(hid < hN[i])
            if len(ok) > 0:
                amin = ok[np.argmin(reldist[ok])]

                self.nearby_halo_properties['Mvir'][self._halo_slices[i]] = Mvir[amin]
                self.nearby_halo_properties['Mstar'][self._halo_slices[i]] = Mstar[amin]
                self.nearby_halo_properties['Rvir'][self._halo_slices[i]] = Rvir[amin]
                self.nearby_halo_properties['N'][self._halo_slices[i]] = hid[amin]
                self.nearby_halo_properties['Mgas'][self._halo_slices[i]] = Mgas[amin]
                self.nearby_halo_properties['SSC'][self._halo_slices[i]] = SSC[amin]
                self.bh['nearhalo'][self._halo_slices[i]] = hid[amin]
                for j in range(3):
                    self.bh['nearpos'][self._halo_slices[i]][:,j] = \
                        self.bh['pos'][self._halo_slices[i]][:,j] + self.halo_properties['SSC'][self._halo_slices[i]][:,j] - SSC[amin][j]

        self.bh['neardist'] = np.sqrt(np.sum(self.bh['nearpos']**2,axis=1))

    def host_prop(self,key):
        return self.halo_properties[key][self._host_indices]

    def host_bhs(self, N, key):
        target = np.where(self.host_ids==N)[0]
        return self.bh[key][self._halo_slices[target]]

    def _get_halo_slices(self):
        ord_ = np.argsort(self.bh['host'])
        uvalues, ind = np.unique(self.bh['host'][ord_], return_index=True)
        slice_ = []
        for i in range(len(uvalues) - 1):
            ss = ord_[ind[i]:ind[i + 1]]
            slice_.append(ss)
        ss = ord_[ind[i + 1]:]
        slice_.append(ss)
        return uvalues, slice_, ord_[ind]

    #def get_mergers(self, ID1, ID2, ratio, time, step):




class BHhalocat(object):

    def __init__(self, simname, boxsize='25 Mpc'):
        self.simname = simname
        if not os.path.exists('steps.list'):
            print "ERROR cannot find steps.list file"
            return
        f = open('steps.list', 'r')
        steps_str = np.array(f.readlines())
        for ii in range(len(steps_str)):
            steps_str[ii] = steps_str[ii].strip('\n')
        steps = steps_str.astype(np.int64)

        self.steps = steps_str
        self.time = np.zeros(len(steps))
        self.boxsize = pynbody.units.Unit(boxsize+' a')

        import halo_db as db

        self.data = StepList(self.steps, db, self.boxsize, self.simname)

    def __getitem__(self,N):
        if type(N)==int:
            return self.data[self.steps[N]]
        if type(N)==str:
            return self.data[N]

    def addsteps(self):
        import halo_db as db
        f = open('steps.list', 'r')
        steps_str = np.array(f.readlines())
        self._steplist = steplist
        for step in steplist:
            print "gathering data for step ", step
            dbstep = db.get_timestep(simname+'/%'+step)
            self.data.data[step] = StepData(step,dbstep,boxsize)
