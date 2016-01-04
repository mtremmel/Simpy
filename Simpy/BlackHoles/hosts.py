from .. import Files, readcol, dbutil
import numpy as np
import pynbody

default_prop_list = ['Mvir', 'Mstar', 'Rvir', 'Mgas']


class BHhalocat(object):

    def __init__(self, simname, boxsize='25 Mpc'):
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

        hostproperties = ['Mvir', 'Mstar', 'Rvir', 'Mgas', 'SSC', 'N']

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
                for key in self.bh.keys():
                    self.bh[key].append(np.array([]))
                for key in self.halo_properties.keys():
                    self.halo_properties[key].append(np.array([]))
                    self.other_halo_properties[key].append(np.array([]))
                continue

            print "quering database for halo properties"
            rawdat = dbstep.gather_property(*hostproperties)
            cnt = 0
            for key in hostproperties:
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
                        self.other_halo_properties[key][-1][i] = np.nan
                    continue
                relpos = hprops['SSC'][ok] - pos[i]
                bad = np.where(np.abs(relpos) > self.boxsize.in_units('kpc', a=a)/2.)
                relpos[bad] = -1.0 * (relpos[bad]/np.abs(relpos[bad])) * \
                              (self.boxsize.in_units('kpc', a=a) - np.abs(relpos[bad]))
                reldist = np.sum(relpos**2, axis=1)
                amin = np.argmin(reldist)
                hnear[i] = hprops['N'][ok][amin]
                distnear[i] = reldist[amin]
                for key in hprops.keys():
                    self.other_halo_properties[key][-1][i] = hprops[key][ok][amin]
            self.bh['nearhalo'].append(hnear)
            self.bh['neardist'].append(distnear)


    def add_host_property(self,keylist):
        import halo_db as db
        prevkeys = self.halo_properties.keys()
        newkeys = []
        for key in keylist:
            if key in prevkeys:
                print key, "already in halo_properties... will not update this"
                continue
            newkeys.append(key)
            self.halo_properties[key] = []
            self.other_halo_properties[key] = []

        for i in range(len(self.steps)):
            step = self.steps[i]
            if len(self.bh['bhid'][i]) == 0:
                for key in newkeys:
                    self.halo_properties[key].append(np.array([]))
                    self.other_halo_properties[key].append(np.array([]))
                continue
            print "getting galaxy data for step ", step
            dbstep = db.get_timestep(self.simname+'/%'+step)

            hprops = {}
            keylist.append('N')
            rawdat = dbstep.gather_property(*keylist)
            cnt = 0
            for key in keylist:
                hprops[key] = rawdat[cnt]
                cnt += 1

            print "connecting properties to hosts"
            uhost, ind = np.unique(self.bh['halo'][i], return_inverse=True)
            match = np.where(np.in1d(hprops['N'],uhost))[0]
            match2 = np.where(np.in1d(uhost,hprops['N']))[0]
            nomatch2 = np.where(np.in1d(uhost,hprops['N'])==False)[0]

            for key in newkeys:
                d1 = hprops[key][match]
                d2 = np.zeros(len(uhost))
                d2[match2] = d1
                d2[nomatch2] = np.nan
                self.halo_properties[key].append(d2[ind])

            print "connecting properties to nearby halos"
            if len(np.where(self.bh['nearhalo']>=0)[0]) == 0:
                for key in newkeys:
                    self.other_halo_properties[key].append(np.array([]))
                continue
            uotherhost, otherind = np.unique(self.bh['nearhalo'][i], return_inverse=True)
            othermatch = np.where(np.in1d(hprops['N'],uotherhost))[0]
            othermatch2 = np.where(np.in1d(uotherhost,hprops['N']))[0]
            othernomatch2 = np.where(np.in1d(uotherhost,hprops['N'])==False)[0]
            if len(othermatch) == 0:
                for key in newkeys:
                    self.other_halo_properties[key].append(np.array([]))
                continue
            for key in newkeys:
                d1 = hprops[key][othermatch]
                d2 = np.zeros(len(uotherhost))
                d2[othermatch2] = d1
                d2[othernomatch2] = np.nan
                self.other_halo_properties[key].append(d2[otherind])

    def calc_lum(self,er=0.1):
        csq = pynbody.array.SimArray((2.998e10) ** 2, 'erg g**-1')
        self.bh['lum'] = []
        for i in range(len(self.steps)):
            self.bh['lum'].append(pynbody.array.SimArray(self.bh['mdot'][i],'Msol yr**-1').in_units('g s**-1') * csq * er)
        return

    def gen_arrays(self):
        bhdata = {}
        hostdata = {}
        neardata = {}
        steps = []
        for i in range(len(self.steps)):
            steps.append(np.ones(len(self.bh['bhid'][i]))*self.steps[i])
        steps = np.concatenate(steps)
        for key in self.bh.keys():
            if key == 'pos':
                continue
            bhdata[key] = np.concatenate(self.bh[key])
        for key in self.halo_properties.keys():
            if key == 'SSC':
                continue
            hostdata = np.concatenate(self.halo_properties[key])
            neardata = np.concatenate(self.other_halo_properties[key])
        return bhdata, hostdata, neardata, steps

    def extract_AGN(self,threshold = 1e42, distmin = 1, distmax = 50):
        if not self.bh['lum']:
            self.calc_lum()
        act_halos = {'uhid':[], 'step':[], 'nbh':[], 'type':[]}

        bhdata, hostdata, neardata, steps = self.gen_arrays()
        bright = np.where(bhdata['lum']>threshold)[0]

        act_halos['step'] = steps[bright]

        otherh - np.where(bhdata['nearhalo'][bright])

        off = np.where(bhdata['neardist'][bright])

        bright = np.where(bhdata['lum']>threshold)
        bhalo = bhdata['halo'][bright]




        data_multi = {'step':[], 'nbhs':[], 'halo':[], 'lummain':[], 'lumother':[],
                'dist':[], 'distother':[], 'mass':[], 'massother':[]}
        for i in range(len(self.steps)):
            step = self.steps[i]
            if len(self.bh['bhid'][i]) == 0:
                continue
            bright = np.where(self.bh['lum'][i] > threshold)

            hosts = self.bh['halo'][i][bright]
            hord = np.argsort(hosts)
            hosts = hosts[ord]
            uhosts, ind, inv, cnt = np.unique(hosts,return_count=True, return_inverse = True, return_index = True)

            multi = np.where(cnt[inv] > 1)[0]
            data_multi['step'].append(np.ones(len(multi)*step))








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











