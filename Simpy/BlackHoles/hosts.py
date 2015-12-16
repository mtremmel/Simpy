from .. import Files, readcol, dbutil
import numpy as np
import pynbody

default_prop_list = ['Mvir', 'Mstar', 'Rvir', 'Mgas', 'MHIGas', 'MColdGas', 'SFR']

class BHhalocat(object):

    def __init__(self, bhorbit, boxsize='25 Mpc', hostproperties=default_prop_list):
        self.simname = bhorbit.simname
        Files.cklists(self.simname)
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
            print "getting black hole data for step ", step
            try:
                bhids = bhorbit.single_step_data(step, 'iord')
            except IndexError:
                bhids = np.array([])
            if len(bhids) == 0:
                for key in self.bh:
                    self.bh[key].append(np.array([]))
                continue
            self.bh['bhid'].append(bhids)
            self.bh['mass'].append(bhorbit.single_step_data(step, 'mass'))
            self.bh['mdot'].append(bhorbit.single_step_data(step, 'mdotmean'))
            self.bh['lum'].append(bhorbit.single_step_data(step, 'lum'))

            hostnum = []
            pos = []
            dist = []

            print "querying database..."
            for id in bhids:
                bh_db = db.get_halo(self.simname+'/%'+str(step)+'/1.'+str(id))
                hostok = True
                try:
                    hid = bh_db.host_halo.halo_number
                except:
                    hid = -1
                    hostok = False
                hostnum.append(hid)
                if hostok:
                    pos.append(bh_db['BH_central_offset'] + bh_db.host_halo['SSC'])
                    dist.append(bh_db['BH_central_distance'])
                else:
                    pos.append(np.array([0,0,0]))
                    dist.append(0)

            hostnum = np.array(hostnum)
            pos = np.vstack(tuple(pos))
            dist = np.array(dist)
            self.bh['halo'].append(hostnum)
            self.bh['pos'].append(pynbody.array.SimArray(pos, 'kpc'))
            self.bh['dist'].append(pynbody.array.SimArray(dist, 'kpc'))

            dbstep = db.get_timestep(self.simname+'/%'+str(step))
            nearhalo = np.ones(len(bhids))*-1
            distnear = np.ones(len(bhids))*-1
            print "finding nearby halos..."
            for halo in dbstep.halos:
                if 'Rvir' not in halo.keys() or 'SSC' not in halo.keys():
                    continue
                if halo.halo_number > hostnum.max():
                    break
                relpos = pos - halo['SSC']
                bad, = np.where(np.abs(relpos) > self.boxsize/2.)
                relpos[bad] = -1.0 * (relpos[bad]/np.abs(relpos[bad])) * (self.boxsize - np.abs(relpos[bad]))
                reldist = np.sqrt(np.sum(relpos**2, axis=1))
                near = np.where((reldist < halo['Rvir']) & ((hostnum > halo.halo_num) | (hostnum == -1)))[0]
                closer = np.where((distnear[near] > 0) & (distnear[near] < reldist[near]))[0]
                np.delete(near, closer)
                nearhalo[near] = halo.halo_number
                distnear[near] = reldist[near]
            self.bh['other_halo'].append(nearhalo)
            self.bh['other_dist'].append(distnear)

        print "finding host galaxy properties"
        if hostproperties and len(hostproperties)>0:
            self.add_host_property(hostproperties)

    def add_host_property(self,keylist):
        for key in keylist:
            self.halo_properties[key] = []
            self.other_halo_properties[key] = []

        for ii in range(len(self.steps)):
            print "finding host galaxy properties for step ", self.steps[ii]
            step = self.steps[ii]
            data = dbutil.property_array(self.simname, step, self.bh['halos'][ii], keylist)
            data_other = dbutil.property_array(self.simname, step, self.bh['other_halos'][ii], keylist)

            for jj in range(len(data)):
                self.halo_properties[keylist[jj]].append(data[jj])

            for jj in range(len(data_other)):
                self.other_halo_properties[keylist[jj]].append(data_other[jj])

    def trace_bh(self, bhid, key, return_steps=False):
        alldata = []
        id = []
        step = []
        if key not in self.halo_properties.keys() and key not in self.bh.keys():
            print "Cannot find given data key in either bh or halo_properties"
            return
        for i in range(len(self.steps)):
            if key in self.bh.keys():
                alldata.extend(self.bh[key][i])
            if key in self.halo_properties.keys():
                alldata.extend(self.halo_properties[key])
            id.extend(self.bh['bhid'])
            step.extend(self.bh['step'])
        alldata = np.array(alldata)
        id = np.array(id)
        step = np.array(step)
        target = np.where(id == bhid)[0]
        target_data = alldata[target]
        ord = np.argsort(step)
        if return_steps is True:
            return target_data[ord[::-1]], step[ord[::-1]]
        else:
            return target_data[ord[::-1]]

    def get_mergers(self):
        time, step, ID, IDeat, ratio, kick = readcol(self.simname + '.mergers', twod=False)

        self.bhmergers['time'] = pynbody.array.SimArray(time, 'Gyr')
        self.bhmergers['ratio'] = ratio
        self.bhmergers['step'] = step
        self.bhmergers['ID'] = ID
        self.bhmergers['EatenID'] = IDeat

        if kick.max() > 0:
            self.bhmergers['kick'] = kick

        self.bhmergers['step_before'] = np.array([])
        self.bhmergers['step_after'] = np.array([])

        for i in range(len(self.bhmergers['time'])):
            self.bhmergers['step_before'] = np.append(self.bhmergers['step_before'],
                                                    self.steps[(self.time <= self.bhmergers['time'][i])][-1])
            self.bhmergers['step_after'] = np.append(self.bhmergers['step_after'],
                                                   self.steps[(self.time >= self.bhmergers['time'][i])][0])

        sord = np.argsort(self.bhmergers['step_after'])
        steps, ind = np.unique(self.bhmergers['step_after'][sord],return_index=True)
        self.bhmergers['hosts'] = np.ones(len(self.bhmergers['ID']))*-1
        for ii in range(len(steps)):
            step = steps[ii]
            ss = sord[ind[ii]:ind[ii]+1]
            idsort = np.argsort(self.bhmergers['ID'][ss])
            stepind = np.where(self.steps==step)[0]
            stepind = stepind[0]
            match = np.where(np.in1d(self.bh['bhid'][stepind],self.bhmergers['ID'][ss][idsort]))[0]

            self.bhmergers['host'][ss][idsort] = self.bh['halo'][stepind][match]











