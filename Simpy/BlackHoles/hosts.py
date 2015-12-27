from .. import Files, readcol, dbutil
import numpy as np
import pynbody

default_prop_list = ['Mvir', 'Mstar', 'Rvir', 'Mgas', 'MHIGas', 'MColdGas', 'SFR']


class BHhalocat(object):

    def __init__(self, simname, boxsize='25 Mpc', hostproperties=default_prop_list):
        self.simname = simname
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

        if 'SSC' not in hostproperties:
            hostproperties.append('SSC')
        if 'N' not in hostproperties:
            hostproperties.append('N')

        if hostproperties and type(hostproperties) == list and len(hostproperties)>0:
            for prop in hostproperties:
                self.halo_properties[prop] = []
                self.other_halo_properties[prop] = []

        for step in self.steps:
            print "getting black hole data for step ", step

            dbstep = db.get_timestep(self.simname+'/%'+step)

            hprops = {}
            if len(self.halo_properties.keys())>0:
                for key in hostproperties:
                    hprops[key] = np.array([])

            print "querying database for BH data..."
            bhids, bhmass, bhmdot, offset, dist, hostnum = \
                dbstep.gather_property(
                    'N', 'BH_mass', 'BH_mdot_ave', 'BH_central_offset', 'BH_central_distance', 'host')

            if len(bhids)==0:
                print "No BHs in halos found in this step"
                continue

            print "quering database for halo properties"
            rawdat = dbstep.gather_property(*hostproperties)
            cnt = 0
            for key in hprops.keys():
                hprops[key] = rawdat[cnt]
                cnt += 1

            self.bh['bhid'].append(bhids)
            self.bh['mass'].append(bhmass)
            self.bh['mdot'].append(bhmdot)
            self.bh['dist'].append(dist)
            self.bh['halo'].append(hostnum)

            uhost, ind = np.unique(hostnum, return_inverse=True)
            match = np.where(np.in1d(hprops['N'],uhost))[0]
            match2 = np.where(np.in1d(uhost,hprops['N']))[0]
            nomatch2 = np.where(np.in1d(uhost,hprops['N'])==False)[0]
            print uhost
            print match2
            for key in hprops.keys():
                d1 = hprops[key][match]
                if key == 'SSC':
                    d2 = np.zeros((len(uhost),3))
                else:
                    d2 = np.zeros(len(uhost))
                d2[match2] = d1
                d2[nomatch2] = np.nan
                self.halo_properties[key].append(d2[ind])

            pos = offset + self.halo_properties['SSC'][-1]
            self.bh['pos'].append(pos)

            a = 1/(1+dbstep.redshift)

            print "looking for closest other halos..."
            hnear = np.ones(len(bhids))*-1
            distnear = np.ones(len(bhids))*-1
            for key in hprops.keys():
                if key == 'SSC':
                    self.other_halo_properties[key].append(np.zeros((len(bhids),3)))
                else:
                    self.other_halo_properties[key].append(np.zeros(len(bhids)))

            for i in range(len(bhids)):
                ok = np.where(hprops['N'] != hostnum[i])[0]
                if len(ok) == 0:
                    hnear[i] = np.nan
                    distnear[i] = np.nan
                    for key in hprops.keys():
                        self.other_halo_properties[key][i] = np.nan
                relpos = hprops['SSC'] - pos[i]
                bad = np.where(np.abs(relpos) > self.boxsize.in_units('kpc', a=a)/2.)
                relpos[bad] = -1.0 * (relpos[bad]/np.abs(relpos[bad])) * \
                              (self.boxsize.in_units('kpc', a=a) - np.abs(relpos[bad]))
                reldist = np.sum(relpos**2, axis=2)
                amin = np.argmin(reldist[(hprops['N']<hostnum[i])])
                hnear[i] = hprops['N'][(hprops['N']<hostnum[i])][amin]
                distnear[i] = reldist[(hprops['N']<hostnum[i])][amin]
                for key in hprops.keys():
                    self.other_halo_properties[key][i] = hprops[key][(hprops['N']<hostnum[i])][amin]


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











