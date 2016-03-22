from .. import readcol,util
import numpy as np
import pynbody
import os
import pickle


def calc_lum(mdot,er=0.1):
    csq = pynbody.array.SimArray((2.998e10) ** 2, 'erg g**-1')
    return pynbody.array.SimArray(mdot,'Msol yr**-1').in_units('g s**-1') * csq * er


class StepList(object):
    def __init__(self, steplist, db, boxsize, simname):
        self._steplist = steplist
        self.data = {}
        self.simname = simname
        for step in steplist:
            print "gathering data for step ", step
            dbstep = db.get_timestep(simname+'/%'+step)
            self.data[step] = StepData(step,dbstep,boxsize)

    def __getitem__(self,step):
        if int(step) not in self._steplist.astype(np.int):
            raise KeyError
        else:
            target = np.where(self._steplist.astype(np.int)==int(step))[0]
        return self.data[self._steplist[target[0]]]

    def add_host_property(self, db, *plist):
        for step in self._steplist:
            dbstep = db.get_timestep(self.simname+'/%'+step)
            self.data[step].add_host_property(dbstep, *plist)

    def add_nearby_property(self, db, *plist):
        for step in self._steplist:
            dbstep = db.get_timestep(self.simname+'/%'+step)
            self.data[step].add_nearby_property(dbstep, *plist)

    def addstep(self,step, db, boxsize):
        dbstep = db.get_timestep(self.simname+'/%'+step)
        self.data[step] = StepData(step,dbstep,boxsize)
        self._steplist = np.append(self._steplist,step)

    def get_BH_mergers(self, time, step, ID, IDeat, ratio, kick):
        for i in range(len(self._steplist)-1):
            curstep = int(self._steplist[i])
            nextstep = int(self._steplist[i+1])
            target = np.where((step > curstep)&(step<nextstep))
            self.data[self._steplist[i]].get_BH_mergers(time[target],step[target],ID[target], IDeat[target],
                                                        ratio[target], kick[target])


class StepData(object):
    def __init__(self, step, dbstep, boxsize):
        self.redshift = dbstep.redshift
        self.time = dbstep.time_gyr
        self.a = 1./(1.+self.redshift)
        self.id = step
        self.boxsize = boxsize
        self.mergers = {}
        self._merger_halo_indices = []
        self._eaten_halo_indices = []

        print "Gathering BH data..."

        bhids, bhmass, bhmdot, offset, dist, hostnum, Mvir, Mstar, Rvir, Mgas, SSC = \
                dbstep.gather_property('halo_number()', 'BH_mass', 'BH_mdot_ave', 'BH_central_offset', 'BH_central_distance',
                                       'host_halo.halo_number()', 'host_halo.Mvir', 'host_halo.Mstar', 'host_halo.Rvir',
                                       'host_halo.Mgas', 'host_halo.SSC')

        nbhs = len(bhids)

        if nbhs==0:
            print "No BHs Found in This Step"
            self.bh = {'lum':[], 'dist':[], 'pos':[], 'host':[], 'mass':[],
                   'bhid':[], 'nearhalo':[],'mdot':[], 'neardist':[]}
            self.halo_properties = {'Mvir':[], 'Mstar':[], 'Rvir':[], 'Mgas':[], 'SSC':[]}
            return

        self.bh = {'lum':calc_lum(bhmdot), 'dist':dist, 'pos':offset, 'host':hostnum, 'mass':bhmass,
                   'bhid':bhids, 'mdot':bhmdot}
        self.halo_properties = {'Mvir':Mvir, 'Mstar':Mstar, 'Rvir':Rvir, 'Mgas':Mgas, 'SSC':SSC, 'N':hostnum}

        print "ordering data..."
        ord = np.argsort(bhids)
        util.cutdict(self.bh,ord)
        util.cutdict(self.halo_properties,ord)

        print "slicing data..."
        self.host_ids, self._halo_slices,self._host_indices = self._get_halo_slices()

        nhalos = len(self.host_ids)

        print "Finding nearby halos..."
        Mvir, Mstar, Rvir, Mgas, SSC, hid = dbstep.gather_property('Mvir', 'Mstar', 'Rvir', 'Mgas', 'SSC', 'halo_number()')

        self.nearby_halo_properties = \
            {'Mvir':np.zeros(nbhs), 'Mstar':np.zeros(nbhs), 'Rvir':np.zeros(nbhs),
             'Mgas':np.zeros(nbhs), 'SSC':np.zeros((nbhs,3)), 'N':np.zeros(nbhs)}
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
                    #print self.bh['pos'][self._halo_slices[i]][:,j], self.halo_properties['SSC'][self._halo_slices[i]][:,j], SSC[amin][j]
                    self.bh['nearpos'][self._halo_slices[i],j] = \
                        self.bh['pos'][self._halo_slices[i],j] + hpos[i][j] - SSC[amin][j]

        self.bh['neardist'] = np.sqrt(np.sum(self.bh['nearpos']**2,axis=1))

        self.near_ids, self._near_slices,self._near_indices = self._get_halo_slices(near=True)

    def host_prop(self,key):
        return self.halo_properties[key][self._host_indices]

    def host_bhs(self, N, key):
        target = np.where(self.host_ids==N)[0]
        return self.bh[key][self._halo_slices[target]]

    def nearby_prop(self,key):
        return self.nearby_halo_properties[key][self._near_indices]

    def _get_halo_slices(self, near=False):
        if near is True:
            key = 'nearhalo'
        else:
            key = 'host'
        ord_ = np.argsort(self.bh[key])
        uvalues, ind = np.unique(self.bh[key][ord_], return_index=True)
        slice_ = []
        for i in range(len(uvalues) - 1):
            ss = ord_[ind[i]:ind[i + 1]]
            slice_.append(ss)
        ss = ord_[ind[i + 1]:]
        slice_.append(ss)
        return uvalues, slice_, ord_[ind]

    def add_host_property(self, dbstep, *plist):
        nbh = len(self.bh['bhid'])
        plist = list(plist)
        finallist = ['halo_number()']
        for key in plist:
            if key not in self.halo_properties.keys():
                self.halo_properties[key] = np.ones(nbh) * -1
                finallist.append('host_halo.'+key)
        if len(finallist)==0:
            return

        finallist = tuple(finallist)
        print "gathering data..."
        data = dbstep.gather_property(*finallist)
        bhid_new = data[0]
        ord_ = np.argsort(bhid_new)
        match = np.in1d(self.bh['bhid'],bhid_new)
        match2 = np.in1d(bhid_new[ord_],self.bh['bhid'])

        if not np.array_equal(self.bh['bhid'][match],bhid_new[ord_][match2]):
            print "ERROR with match"
            return

        cnt = 1
        for key in plist:
            self.halo_properties[key][match] = data[cnt][ord_][match2]
            cnt += 1

    def add_nearby_property(self,dbstep,*plist):
        nbh = len(self.bh['bhid'])
        plist = list(plist)
        finallist = ['halo_number()']
        for key in plist:
            if key not in self.nearby_halo_properties.keys():
                self.nearby_halo_properties[key] = np.ones(nbh) * -1
                finallist.append(key)
        if len(finallist)==0:
            return
        finallist = tuple(finallist)
        data = dbstep.gather_property(*finallist)
        hid = data[0]
        ord_ = np.argsort(self.nearby_halo_properties['N'])
        ord2_ = np.argsort(hid)
        match = np.in1d(self.nearby_halo_properties['N'][ord_],hid)
        match2 = np.in1d(hid[ord2_], self.nearby_halo_properties['N'])
        umatch, uinv = np.unique(self.nearby_halo_properties['N'][ord_][match],return_inverse=True)
        if not np.array_equal(self.nearby_halo_properties['N'][ord_][match],hid[ord2_][match2][uinv]):
            print "ERROR in matching"
            return
        cnt = 1
        for key in plist:
            self.nearby_halo_properties[key][ord_][match] = data[cnt][ord2_][match2][uinv]
            cnt += 1

    def get_BH_mergers(self, time, step, ID1, ID2, ratio, kick):
        self.mergers = {'bhid': ID1, 'eaten_bhid': ID2, 'ratio': ratio, 'time': time,
                        'step': step, 'halo': -1*np.ones(len(ID1)), 'eaten_halo': -1*np.ones(len(ID1)),
                        'mass1': np.zeros(len(ID1)), 'mass2': np.zeros(len(ID1)),
                        'lum1': np.zeros(len(ID1)), 'lum2': np.zeros(len(ID1)), 'kick':kick}

        self._merger_halo_indices = []
        self._eaten_halo_indices = []

        for i in range(len(ID1)):
            o1, = np.where(self.bh['bhid']==ID1[i])
            o2, = np.where(self.bh['bhid']==ID2[i])
            if len(o1) > 0:
                o1 = o1[0]
                h1 = self.bh['host'][o1]
                self.mergers['halo'][i] = h1
                m1 = self.bh['mass'][o1]
                self.mergers['mass1'][i] = m1
                md1 = self.bh['mdot'][o1]
                self.mergers['lum1'][i] = calc_lum(md1)
                h1index, = np.where(self.host_ids==h1)
                if len(h1index)>0:
                    h1index = h1index[0]
                    self._merger_halo_indices.append(h1index)
            if len(o2) > 0:
                h2 = self.bh['host'][o2]
                self.mergers['eaten_halo'][i] = h2
                m2 = self.bh['mass'][o2]
                self.mergers['mass2'][i] = m2
                md2 = self.bh['mdot'][o2]
                self.mergers['lum2'][i] = calc_lum(md2)
                h2index, = np.where(self.host_ids==h2)
                if len(h2index)>0:
                    h2index = h2index[0]
                    self._eaten_halo_indices.append(h2index)

    def BH_merger_halo_props(self,key):
        main = self.host_prop(key)[self._merger_halo_indices]
        other = self.host_prop(key)[self._eaten_halo_indices]
        return np.concatenate([[main],[other]]).T



class BHhalocat(object):

    def __init__(self, simname, boxsize='25 Mpc', filename=None):
        self.filename = filename
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

        if filename:
            if os.path.exists(filename):
                print "ERROR ", filename, " already exists. NOT SAVING OBJECT TO FILE"
            else:
                fout = open(filename,'wb')
                pickle.dump(self,fout)
                fout.close()

    def __getitem__(self,N):
        if type(N)==int:
            return self.data[self.steps[N]]
        if type(N)==str:
            return self.data[N]

    def save(self, newname=None):
        if newname:
            if os.path.exists(newname):
                print "ERROR ", newname, " already exists. NOT SAVING OBJECT TO FILE"
                return
            else:
                self.filename=newname
                if os.path.exists(self.filename):
                    os.system('rm self.filename')
        fout = open(self.filename,'wb')
        pickle.dump(self, fout)
        fout.close()

    def addsteps(self, *newsteps):
        import halo_db as db
        if len(newsteps) == 0:
            print "checking for new steps in file steps.list"
            f = open('steps.list', 'r')
            steps_str = np.array(f.readlines())
        else:
            steps_str = newsteps
        for step in steps_str:
            if step in self.steps:
                continue
            self.steps = np.append(self.steps,step)
            print "gathering data for step ", step
            self.data.addstep(step, db, self.boxsize)

    def add_host_property(self,*plist):
        import halo_db as db
        self.data.add_host_property(db, *plist)

    def add_nearby_property(self,*plist):
        import halo_db as db
        self.data.add_nearby_property(db, *plist)


    def get_mergers(self,simname):
        time, step, ID, IDeat, ratio, kick = readcol.readcol(simname+'.mergers',twod=False)
        self.data.get_BH_mergers(time, step, ID, IDeat, ratio, kick)