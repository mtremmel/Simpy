from .. import Files, readcol
import numpy as np
import pynbody

default_prop_list = ['Mvir', 'Mstar', 'Rvir', 'Mgas', 'MHIGas', 'MColdGas']

class BHhalocat(object):

    def __init__(self, simname, bhorbit, boxsize='25 Mpc', hostproperties=default_prop_list):
        self.simname = simname
        Files.cklists(simname)
        f = open('steps.list', 'r')
        steps = np.array(f.readlines()).astype(np.int64)
        self.steps = steps
        self.time = np.zeros(len(steps))
        self.host_properties = {}
        self.boxsize = pynbody.units.Unit(boxsize+' a')

        import halo_db as db

        self.bh = {'lum':[], 'dist':[], 'pos':[], 'halo':[], 'mass':[], 'bhid':[], 'tophalo':[],'mdot':[], 'other_dist':[], 'other_halo':[]}
        self.halo_properties = {}
        self.other_halo_properties = {}
        self.bhmergers = {}

        for step in self.steps:
            bhids = bhorbit.single_step_data(step, 'iord')
            if len(bhids) == 0:
                for key in self.bh:
                    self.bh['key'] = np.array([])
                continue
            self.bh['bhid'].append(bhids)
            self.bh['mass'].append(bhorbit.single_step_data(step,'mass'))
            self.bh['mdot'].append(bhorbit.single_step_data(step,'mdotmean'))
            self.bh['lum'].append(bhorbit.single_step_data(step,'lum'))

            hostnum = []
            pos = []
            dist = []

            for id in bhids:
                bh_db = db.get_halo(self.simname+'/%'+str(step)+'/1.'+str(id))
                hostnum.append(bh_db.host_halo.halo_number)
                pos.append(bh_db['BH_central_offset'] + bh_db.host_halo['SSC'])
                dist.append(bh_db['BH_central_distance'])

            hostnum = np.array(hostnum)
            pos = np.vstack(tuple(pos))
            dist = np.array(dist)
            self.bh['halo'].append(hostnum)
            self.bh['pos'].append(pynbody.array.SimArray(pos,'kpc'))
            self.bh['dist'].append(pynbody.array.SimArray(dist,'kpc'))

            dbstep = db.get_step(self.simname+'/%'+str(step))
            nearhalo = np.ones(len(bhids)*-1)
            distnear = np.ones(len(bhids)*-1)
            for halo in dbstep.halos:
                if 'Rvir' not in halo.keys() or 'SSC' not in halo.keys():
                    continue
                if halo.halo_number > hostnum.max():
                    break
                relpos = pos - halo['SSC']
                bad, = np.where(np.abs(relpos) > self.boxsize/2.)
                relpos[bad] = -1.0 * (relpos[bad]/np.abs(relpos[bad])) * (self.boxsize - np.abs(relpos[bad]))
                reldist = np.sqrt(np.sum(relpos**2, axis=1))
                near = np.where((reldist < halo['Rvir']) & (hostnum > halo.halo_num))[0]
                closer = np.where((distnear[near] > 0) & (distnear[near] < reldist[near]))[0]
                np.delete(near, closer)
                nearhalo[near] = halo.halo_number
                distnear[near] = reldist[near]
            self.bh['other_halo'].append(nearhalo)
            self.bh['other_dist'].append(distnear)

        if hostproperties and len(hostproperties)>0:
            self.add_host_property(hostproperties)

    def add_host_property(self,keylist):
        import halo_db as db
        for key in keylist:
            self.halo_properties[key] = []
            self.other_halo_properties[key] = []
        for ii in range(len(self.steps)):
            for key in keylist:
                self.halo_properties[key].append(-1*np.ones(len(self.bh['halo'][ii])))
                self.other_halo_properties[key].append(-1*np.ones(len(self.bh['other_halo'][ii])))
            ord_ = np.argsort(self.bh['halo'][ii])
            o_ord_ = np.argsort(self.bh['other_halo'][ii])

            uhalos, ind = np.unique(self.bh['halo'][ii][ord_], return_index=True)
            o_uhalos, o_ind = np.unique(self.bh['other_halo'][ii][o_ord_], return_index=True)

            for jj in range(len(uhalos)):
                hid = uhalos[jj]
                h = db.get_halo(self.simname+'/%'+str(self.steps[ii])+'/'+str(hid))
                for key in keylist:
                    if key not in h.keys():
                        continue
                    if jj + 1 < len(ind):
                        ss = ord_[ind[jj]:ind[jj + 1]]
                    else:
                        ss = ord_[ind[jj]:]
                    self.halo_properties[key][ii][ss] = h[key]

            for jj in range(len(o_uhalos)):
                hid = o_uhalos[jj]
                h = db.get_halo(self.simname+'/%'+str(self.steps[ii])+'/'+str(hid))
                for key in keylist:
                    if key not in h.keys():
                        continue
                    if jj + 1 < len(ind):
                        ss = o_ord_[o_ind[jj]:o_ind[jj + 1]]
                    else:
                        ss = o_ord_[o_ind[jj]:]
                    self.other_halo_properties[key][ii][ss] = h[key]

    def get_mergers(self):
        time, step, ID, IDeat, ratio, kick = readcol(self.simname + '.mergers', twod=False)

        self.bhmergers['time'] = pynbody.array.SimArray(time,'Gyr')
        self.bhmergers['ratio'] = ratio
        self.bhmergers['step'] = step
        self.bhmergers['ID'] = ID
        self.bhmergers['EatenID'] = IDeat

        if kick.max() > 0:
            self.bhmergers['kick'] = kick

        self.bhmergers['step_before'] = np.array([])
        self.bhmergers['step_after'] = np.array([])

        for i in range(len(self.mergers['time'])):
            self.bhmergers['step_before'] = np.append(self.mergers['step_before'],
                                                    self.step[(self.time <= self.mergers['time'][i])][-1])
            self.bhmergers['step_after'] = np.append(self.mergers['step_after'],
                                                   self.step[(self.time >= self.mergers['time'][i])][0])


        sord = np.argsort(self.mergers['step_after'])
        steps, ind = np.unique(self.meregers['step_after'][sord],return_index=True)
        self.mergers['hosts'] = np.ones(len(self.mergers['ID']))*-1
        for ii in range(len(steps)):
            step = steps[ii]
            ss = sord[ind[ii]:ind[ii]+1]
            idsort = np.argsort(self.mergers['ID'][ss])
            stepind = np.where(self.steps==step)[0]
            stepind = stepind[0]
            match = np.where(np.in1d(self.bh['bhid'][stepind],self.mergers['ID'][ss][idsort]))[0]

            self.bhmergers['host'][ss][idsort] = self.bh['halo'][stepind][match]












