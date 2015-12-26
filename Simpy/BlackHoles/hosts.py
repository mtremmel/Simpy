from .. import Files, readcol, dbutil
import numpy as np
import pynbody

default_prop_list = ['Mvir', 'Mstar', 'Rvir', 'Mgas', 'MHIGas', 'MColdGas', 'SFR']


class BHhalocat(object):

    def __init__(self, bhorbit, boxsize='25 Mpc', hostproperties=default_prop_list):
        self.simname = bhorbit.simname
        Files.cklists(self.simname)
        f = open('steps.list', 'r')
        steps_str = np.array(f.readlines())
        for ii in range(len(steps_str)):
            steps_str[ii] = steps_str[ii].strip('\n')
        steps = steps_str.astype(np.int64)

        self.steps = steps_str
        self.time = np.zeros(len(steps))
        self.host_properties = {}
        self.boxsize = pynbody.units.Unit(boxsize+' a')

        import halo_db as db

        self.bh = {'lum':[], 'dist':[], 'pos':[], 'halo':[], 'mass':[],
                   'bhid':[], 'nearhalo':[],'mdot':[], 'neardist':[]}
        self.halo_properties = {}
        self.other_halo_properties = {}
        self.bhmergers = {}

        if hostproperties and type(hostproperties) == list and len(hostproperties)>0:
            for prop in hostproperties:
                self.halo_properties[prop] = []
                self.other_halo_properties[prop] = []

        for step in self.steps:
            print "getting black hole data for step ", step
            try:
                bhids = bhorbit.single_step_data(step.astype(np.int64), 'iord')
            except IndexError:
                bhids = np.array([])
            if len(bhids) == 0:
                for key in self.bh:
                    self.bh[key].append(np.array([]))
                continue
            self.bh['bhid'].append(bhids)
            self.bh['mass'].append(bhorbit.single_step_data(step.astype(np.int64), 'mass'))
            self.bh['mdot'].append(bhorbit.single_step_data(step.astype(np.int64), 'mdotmean'))
            self.bh['lum'].append(bhorbit.single_step_data(step.astype(np.int64), 'lum'))

            hostnum = []
            pos = []
            dist = []
            hprops = {}
            other_hprops = {}
            if len(self.halo_properties.keys())>0:
                for key in hostproperties:
                    hprops[key] = []
                    other_hprops[key] = []

            print "querying database for BH data, host properties..."
            for id in bhids:
                bh_db = db.get_halo(self.simname+'/%.'+step+'/1.'+str(id))
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
                    for key in hprops.keys():
                        if key in bh_db.host_halo.keys():
                            hprops[key].append(bh_db.host_halo[key])
                        else:
                            hprops[key].append(np.nan)
                else:
                    pos.append(np.array([0,0,0]))
                    dist.append(0)
                    for key in hprops.keys():
                        hprops[key].append(np.nan)

            hostnum = np.array(hostnum)
            pos = np.vstack(tuple(pos))
            dist = np.array(dist)
            self.bh['halo'].append(hostnum)
            self.bh['pos'].append(pynbody.array.SimArray(pos, 'kpc'))
            self.bh['dist'].append(pynbody.array.SimArray(dist, 'kpc'))

            for key in hprops.keys():
                self.halo_properties[key].append(np.array(hprops[key]))

            dbstep = db.get_timestep(self.simname+'/%'+step)
            a = 1/(1+dbstep.redshift)
            haloids = []
            halossc = []

            print "getting step halo data..."
            for halo in dbstep.halos:
                if 'Rvir' not in halo.keys() or 'SSC' not in halo.keys():
                    continue
                if halo.halo_number > hostnum.max():
                    break
                if halo.halo_type == 1:
                    break
                halossc.append(halo['SSC'])
                haloids.append(halo.halo_number)
                for key in other_hprops.keys():
                    if key in halo.keys():
                        other_hprops[key].append(halo[key])
                    else:
                        other_hprops[key].append(np.nan)

            if len(halossc) == 0:
                print "No halo data found in this snapshot"
                continue
            halossc = pynbody.array.SimArray(np.vstack(halossc),'kpc')
            haloids = np.array(haloids)
            for key in other_hprops.keys():
                other_hprops[key] = np.array(other_hprops[key])

            print "Searching for closest nearby halos..."

            sscl = []
            for i in range(len(bhids)):
                sscl.append(halossc)
            ssca = np.vstack(sscl).reshape((len(bhids), len(haloids), 3))

            posl = []
            for i in range(len(haloids)):
                posl.append(pos)
            posa = np.vstack(posl).reshape((len(bhids), len(haloids), 3))

            relpos = posa - ssca
            bad = np.where(np.abs(relpos) > self.boxsize.in_units('kpc', a=a)/2.)
            relpos[bad] = -1.0 * (relpos[bad]/np.abs(relpos[bad])) * (self.boxsize.in_units('kpc', a=a) - np.abs(relpos[bad]))
            reldist = np.sum(relpos**2, axis=2)
            dsort = reldist.argsort(axis=1)

            halonear = []
            distnear = []
            indnear = []
            for ii in range(len(bhids)):
                if hostnum[ii] > 0:
                    ok = np.where(haloids[dsort[ii]] < hostnum[ii])[0]
                    if len(ok)>0:
                        halonear.append(haloids[dsort[ii]][ok[0]])
                        distnear.append(reldist[ii][dsort[ii]][ok[0]])
                        indnear.append(dsort[ii][ok[0]])
                    else:
                        halonear.append(-1)
                        distnear.append(-1)
                else:
                    halonear.append(haloids[dsort[ii]][0])
                    distnear.append(reldist[ii][dsort[ii]][0])
                    indnear.append(dsort[ii][0])
            self.bh['nearhalo'].append(np.array(halonear))
            self.bh['neardist'].append(np.array(distnear))
            indnear = np.array(indnear)

            for key in other_hprops.keys():
                self.other_halo_properties[key].append(other_hprops[key][indnear])

    def add_host_property(self,keylist):
        for key in keylist:
            self.halo_properties[key] = []
            self.other_halo_properties[key] = []

        for ii in range(len(self.steps)):
            print "finding host galaxy properties for step ", self.steps[ii]
            step = self.steps[ii]
            data = dbutil.property_array(self.simname, step, np.append(self.bh['halos'][ii],self.bh['other_halos'][ii]), keylist)

            for jj in range(len(data)):
                self.halo_properties[keylist[jj]].append(data[jj][:len(self.bh['halos'][ii])])
                self.other_halo_properties[keylist[jj]].append(data[jj][len(self.bh['halos'][ii]):])

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
            else:
                alldata.extend(self.halo_properties[key])
            id.extend(self.bh['bhid'][i])
            step.extend(np.ones(len(self.bh['bhid'])) * self.steps[i])
        alldata = np.array(alldata)
        id = np.array(id)
        step = np.array(step)
        target = np.where(id == bhid)[0]
        target_data = alldata[target]
        step_target = step[target]
        ord = np.argsort(step_target)
        if return_steps is True:
            return target_data[ord], step_target[ord]
        else:
            return target_data[ord]

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











