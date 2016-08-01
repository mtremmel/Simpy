import numpy as np
from .. import Files, cosmology, plotting, util
from ..readcol import readcol
import os
import pynbody
import re

def findmergers(outname, diagname='*out*'):
    os.system("awk '/BHSink/ && /Merge/ && /eating/' " + diagname + " > " + outname)
    return


def reducedata(simname, RetData=False, outname='*out*', mergename='BHmerge.txt', NCHILADA=True):
    if not os.path.exists(mergename):
        print "didn't find merger file... getting mergers from outputs"
        findmergers(mergename, diagname=outname)
    if not os.path.exists('files.list'):
        Files.getFileLists(simname, NCHILADA=NCHILADA)
    f = open('files.list', 'r')
    farr = f.readlines()
    f.close()
    f = open('steps.list', 'r')
    sarr = f.readlines()
    f.close()

    s = pynbody.load(farr[-1].strip('\n'))
    lstep = float(sarr[-1].strip('\n'))

    simtime = s.properties['time'].in_units('Gyr')#cosmology.getTime(s.properties['a'] ** -1 - 1, s)
    dt = simtime / lstep
    tunit = s.infer_original_units('Gyr')

    a, b, ID, c, IDeat, d, time, e, f, kick, g, h, mratio = readcol(mergename, twod=False)
    output = {}
    time = pynbody.array.SimArray(time, tunit)
    output['time'] = time.in_units('Gyr')

    output['ID'] = ID
    output['IDeat'] = IDeat
    output['step'] = time.in_units('Gyr') / dt
    output['ratio'] = mratio
    output['kick'] = kick

    output['ID'][(output['ID']<0)] = 2*2147483648 + output['ID'][(output['ID']<0)]
    output['IDeat'][(output['IDeat']<0)] = 2*2147483648 + output['IDeat'][(output['IDeat']<0)]

    outkeys = ['time', 'step', 'ID', 'IDeat', 'ratio', 'kick']
    tofile = []
    for key in outkeys:
        tofile.append(output[key])
    tofile = tuple(tofile)
    print "saving to file..."
    np.savetxt(simname + '.mergers', np.column_stack(tofile), fmt=['%f', '%f', '%d', '%d', '%f', '%f'])
    if RetData == True:
        return output
    else:
        return

def get_complete_prog_list(bhorbit, bhid, tmax):
    target, = np.where(bhorbit.bhiords == bhid)
    bhorbit.getprogbhs()
    idlist = np.array([])
    if len(bhorbit.prog['iord'][target])==0:
        return np.array([])
    idnew = np.array(bhorbit.prog['iord'][target])[(np.array(bhorbit.prog['time'][target])<tmax)]
    idlist = np.append(idlist,idnew)

    deep = 0
    while len(idnew) > 0:
        deep += 1
        idnext = np.array([])
        for bid in idnew:
            newtarget, = np.where(bhorbit.bhiords==bid)
            if len(newtarget)>1:
                print "Warning! multiple matches in bhiords found for ", bid
            if len(newtarget)==0:
                print str(bid)+" not found in orbit object! moving on..."
                continue
            newtarget = newtarget[0]
            newpart = np.array(bhorbit.prog['iord'][newtarget])
            idnext = np.append(idnext,newpart)
        idnew = idnext
        idlist = np.append(idlist,idnew)
    print "finished with ", deep, "steps\n"
    return idlist

def plt_merger_rates(time,sim, color='b',linestyle='-', vol_weights=1./25.**3, bins=50,
                     tzrange=[0,25], xlog=True, ylog=True, lw=3, label=None, ret_data=False,
                     pltredshift=True):
    if type(vol_weights)==list or type(vol_weights)==np.ndarray:
        if len(vol_weights)!=1 and len(vol_weights) != len(time):
            print "ERROR do not understand vol_weights format... aborting"
            return
    else:
        print "here"
        vol_weights = np.ones(len(time))*vol_weights
    if type(bins)!= list and type(bins)!=np.ndarray:
        dtz = (tzrange[1]-tzrange[0])/float(bins)
        tzbins = np.arange(tzrange[0],tzrange[1]+dtz,dtz)
    else:
        tzbins = bins
    if pltredshift:
        tedges = np.array([cosmology.getTime(z,sim) for z in tzbins])
    else:
        tedges = tzbins

    tsorted = np.argsort(tedges)
    tedges = tedges[tsorted]

    dt = np.abs((tedges[0:-1] - tedges[1:]) * 1e9)

    data = np.histogram(time, bins=tedges, weights=vol_weights)
    rate = data[0]/dt
    tzbins = tzbins[tsorted]
    if pltredshift is False:
        plotting.plt.step(tzbins[0:],np.append(rate,rate[-1]), color=color,linestyle=linestyle, label=label, linewidth=lw, where='post')
        plotting.plt.xlabel('Time (Gyr)')
    else:
        if xlog is False:
            plotting.plt.step(tzbins[0:],np.append(rate,rate[-1]), color=color,linestyle=linestyle, label=label, linewidth=lw, where='post')

        else:
            plotting.plt.step(tzbins[0:]+1,np.append(rate,rate[-1]), color=color,linestyle=linestyle, label=label, linewidth=lw, where='post')

    plotting.plt.xlabel('Redshift')

    if xlog is True:
        plotting.plt.xscale('log',base=10)
        if pltredshift is True:
            plotting.plt.xticks([1,2,3,4,5,6,9,11,16,21],['0','1','2','3','4','5','8','10','15','20'])
    if ylog is True:
        plotting.plt.yscale('log',base=10)

    if ret_data is True:
        return rate, tzbins,tedges


class mergerCat(object):
    def __init__(self,simname):
        self.simname = simname
        self.db_mergers = {}
        mergerfile = simname+'.mergers'
        print "reading .mergers file..."
        time, step, ID, IDeat, ratio, kick = readcol(mergerfile,twod=False)

        print "checking for bad IDs..."
        bad = np.where(ID<0)[0]
        if len(bad)>0:
            ID[bad] = 2*2147483648 + ID[bad]
        bad2 = np.where(IDeat<0)[0]
        if len(bad2)>0:
            IDeat[bad2] = 2*2147483648 + IDeat[bad2]

        uIDeat, indices = np.unique(IDeat, return_index=True)

        self.rawdat = {'time':time, 'ID1':ID, 'ID2':IDeat, 'ratio':ratio, 'kick':kick, 'step':step}
        util.cutdict(self.rawdat,indices)
        ordr = np.argsort(self.rawdat['ID2'])
        util.cutdict(self.rawdat,ordr)

        f = open('files.list','r')
        s = pynbody.load(f.readline().strip('\n'))
        f.close()
        scale, red = cosmology.getScaleFactor(pynbody.array.SimArray(self.rawdat['time'],'Gyr'),s)
        self.rawdat['redshift'] = red

        uIDeat, cnt = np.unique(self.rawdat['ID2'], return_counts=True)
        if len(np.where(cnt>1)[0])>0:
            print "SHIIIIIIT"

    def get_halo_info(self,dbsim,halo_props = ['halo_number()','Mvir', 'Mstar', 'Mgas']):
        steps, = readcol('steps.list',twod=False)
        import tangos as db
        self.rawdat['snap_after'] = np.zeros(len(self.rawdat['ID1'])).astype('S100')
        self.rawdat['snap_before'] = np.zeros(len(self.rawdat['ID1'])).astype('S100')
        for p in halo_props:
            self.rawdat[p] = np.ones(len(self.rawdat['ID1'])) * -1
        for i in range(len(self.rawdat['ID1'])):
            if i%100 == 0: print float(i)/len(self.rawdat['ID1']) * 100, '% done'
            id1 = self.rawdat['ID1'][i]
            id2 = self.rawdat['ID2'][i]
            if self.rawdat['step']>steps.max():
                continue
            ind = np.where(self.rawdat['step'][i]<=steps)[0][0]
            self.rawdat['snap_after'][i] = dbsim.timesteps[ind].extension
            if ind > 0:
                self.rawdat['snap_before'][i] = dbsim.timesteps[ind-1].extension
            bh = db.get_halo(str(dbsim.timesteps[ind].path)+'/'+str(id1))
            if bh is None:
                bh = db.get_halo(str(dbsim.timesteps[ind-1].path)+'/'+str(id1))
            if bh is None:
                bh = db.get_halo(str(dbsim.timesteps[ind-1].path)+'/'+str(id2))
            if bh is None:
                continue
            for p in halo_props:
                self.rawdat[p][i] = bh.calculate('host_halo.'+p)


    def get_db_data(self,dbsim,properties=['host_halo.Mvir', 'host_halo.Mstar', 'host_halo.Mgas']):
        proplist = ['halo_number()', 'BH_merger_next.halo_number()', 'host_halo.halo_number()',
                    'BH_merger_next.host_halo.halo_number()', 'BH_merger_next.earlier(1).host_halo.halo_number()']
        for prop in properties:
            proplist.append(prop)
            proplist.append('BH_merger_next.'+prop)
            proplist.append('BH_merger_next.earlier(1).'+prop)

        self.db_mergers = {'ID1':[], 'ID2':[], 'ratio':[], 'kick':[], 'step':[],
                    'time':[], 'tsnap_prev':[], 'tsnap_after':[], 'snap_prev':[], 'snap_after':[],
                    'host_N_1':[], 'host_N_2':[], 'host_N_f':[], 'redshift':[]}

        for p in properties:
            self.db_mergers[p+'_1'] = []
            self.db_mergers[p+'_2'] = []
            self.db_mergers[p+'_f'] = []


        self.nmergers = []
        self.steptimes = []

        for step in dbsim.timesteps:
            print step
            try:
                data = step.gather_property(*proplist)
            except:
                print "Nothing found in this step"
                self.nmergers.append(0)
                self.steptimes.append(step.time_gyr)
                continue

            stepnum = re.match("^(.*)\.(0[0-9]*)$",step.filename).groups()[1]
            stepnumA = re.match("^(.*)\.(0[0-9]*)$",step.next.filename).groups()[1]

            self.nmergers.append(len(data[0]))
            self.steptimes.append(step.time_gyr)

            self.db_mergers['ID1'].extend(data[1])
            self.db_mergers['ID2'].extend(data[0])

            tt1, tt2 = np.unique(data[0],return_counts=True)
            if len(np.where(tt2>1)[0])>0:
                print "Double counted IDs: ", tt1[(tt2>1)]
                raise RuntimeError("ERROR double counted IDeat in database analysis")

            self.db_mergers['host_N_1'].extend(data[4])
            self.db_mergers['host_N_2'].extend(data[2])
            self.dadb_mergersta['host_N_f'].extend(data[3])

            self.db_mergers['snap_after'].extend(int(stepnumA) * np.ones(len(data[1])).astype(np.int))
            self.db_mergers['snap_prev'].extend(int(stepnum) * np.ones(len(data[1])).astype(np.int))
            self.db_mergers['tsnap_after'].extend(step.next.time_gyr * np.ones(len(data[1])))
            self.db_mergers['tsnap_prev'].extend(step.time_gyr * np.ones(len(data[1])))

            index = 5
            for i in range(len(properties)):
                self.db_mergers[properties[i]+'_1'].extend(data[index+2])
                self.db_mergers[properties[i]+'_2'].extend(data[index])
                self.db_mergers[properties[i]+'_f'].extend(data[index+1])
                index += 3

        for key in self.db_mergers.keys():
            self.db_mergers[key] = np.array(self.db_mergers[key])

        self._match_data_to_raw('ratio', 'kick', 'time','step', 'redshift')

    def __getitem__(self, item):
        return self.rawdat[item]

    def keys(self):
        return self.rawdat.keys()

    def get_snap_name(self):
        if 'step' not in self.data.keys():
            print "ERROR 'step' is not a property that has been calculated"
            raise RuntimeError
        small = np.where(self.data['step_after']<100)[0]
        med = np.where(self.data['step_after']<1000)[0]
        big = np.where(self.data['step_after']>=1000)[0]
        self.data['stepname'] = np.zeros(len(self.data['ID1'])).astype('S100')
        if len(small)>0:
            self.data['stepname'][small] = self.simname+'.0000'+self.data['step_after'][small]
        if len(med)>0:
            self.data['stepname'][med] = self.simname+'.000'+self.data['step_after'][med]
        if len(big)>0:
            self.data['stepname'][big] = self.simname+'.00'+self.data['step_after'][big]

    def _match_data_to_raw(self,*properties):
        ordee = np.argsort(self.db_mergers['ID2'])
        match = np.where(np.in1d(self.db_mergers['ID2'][ordee],self.rawdat['ID2']))[0]
        match2 = np.where(np.in1d(self.rawdat['ID2'],self.data['ID2'][ordee]))[0]

        if len(match) != len(match2):
            print self.db_mergers['ID2'][ordee], self.rawdat['ID2']
            print match, match2
            raise RuntimeError("ERROR match not returning same number of elements")

        if len(np.where(np.equal(self.db_mergers['ID2'][ordee[match]],self.rawdat['ID2'][match2]) is False)[0])>0:
            raise RuntimeError("ERROR something wrong with array matching")


        for p in properties:
            self.db_mergers[p] = np.ones(len(self.db_mergers['ID2'])) * -1
            self.db_mergers[p][ordee[match]] = self.rawdat[p][match2]


    def calc_GW_emission(self, sim):
        from . import gw

        if 'merge_mass_1' not in self.rawdat.keys():
            print "ERROR need to find final masses for BHs"
            return

        ok = np.where((self.rawdat['merge_mass_1'] > 0)&(self.rawdat['merge_mass_2']>0))[0]

        self.gwemit = gw.GWemit(self.rawdat['merge_mass_1'][ok], self.rawdat['merge_mass_2'][ok],self.rawdat['redshift'])

        self.rawdat['GW_freq_merge'] = np.ones(len(self.rawdat['ID1']))*-1
        self.rawdat['GW_strain_merge'] = np.ones(len(self.rawdat['ID1']))*-1
        self.rawdat['GW_freq_ring'] = np.ones(len(self.rawdat['ID1']))*-1
        self.rawdat['GW_strain_ring'] = np.ones(len(self.rawdat['ID1']))*-1

        self.rawdat['GW_freq_merge'][ok] = self.gwemit.freq_merge()
        self.rawdat['GW_freq_ring'][ok] = self.gwemit.freq_ring()
        self.rawdat['GW_strain_merge'][ok] = self.gwemit.ampGW_merger()
        self.rawdat['GW_strain_ring'][ok] = self.gwemit.ampGW_ring()

        self._match_data_to_raw('GW_freq_merge', 'GW_freq_ring', 'GW_strain_merge', 'GW_strain_ring')

    def get_final_values(self,bhorbit):
        self.rawdat['merge_mass_2'] = np.ones(len(self.rawdat['ID1']))*-1
        self.rawdat['merge_mass_1'] = np.ones(len(self.rawdat['ID1']))*-1
        self.rawdat['merge_mdot_2'] = np.ones(len(self.rawdat['ID1']))*-1
        self.rawdat['merge_mdot_1'] = np.ones(len(self.rawdat['ID1']))*-1
        self.rawdat['merge_lum_2'] = np.ones(len(self.rawdat['ID1']))*-1
        self.rawdat['merge_lum_1'] = np.ones(len(self.rawdat['ID1']))*-1
        self.rawdat['tform1'] = np.ones(len(self.rawdat['ID1']))*-1
        self.rawdat['tform2'] = np.ones(len(self.rawdat['ID1']))*-1
        for i in range(len(self.rawdat['ID1'])):
            mass1 = bhorbit.single_BH_data(self.rawdat['ID1'][i],'mass')
            mass2 = bhorbit.single_BH_data(self.rawdat['ID2'][i],'mass')

            mdot1 = bhorbit.single_BH_data(self.rawdat['ID1'][i],'mdotmean')
            mdot2 = bhorbit.single_BH_data(self.rawdat['ID2'][i],'mdotmean')

            lum1 = bhorbit.single_BH_data(self.rawdat['ID1'][i],'lum')
            lum2 = bhorbit.single_BH_data(self.rawdat['ID2'][i],'lum')

            time1 = bhorbit.single_BH_data(self.rawdat['ID1'][i],'time')
            time2 = bhorbit.single_BH_data(self.rawdat['ID2'][i],'time')

            if len(time1)>0:
                oo = np.where(bhorbit.bhiords == self.rawdat['ID1'][i])[0]
                oo = oo[0]
                self.rawdat['tform1'][i] = bhorbit.tform[oo]
            if len(time2)>0:
                oo = np.where(bhorbit.bhiords == self.rawdat['ID2'][i])[0]
                oo = oo[0]
                self.rawdat['tform2'][i] = bhorbit.tform[oo]


            if len(mass2)>0:
                self.rawdat['merge_mass_2'][i] = mass2[-1]
                self.rawdat['merge_mdot_2'][i] = mdot2[-1]
                self.rawdat['merge_lum_2'][i] = lum2[-1]

            if len(mass1)>0:
                if len(time2)>0:
                    argm = np.argmin(np.abs(time2[-1]-time1))
                else:
                    argm = np.argmin(np.abs(self.rawdat['time'][i]-time1))
                    if argm > 0:
                        argm = argm - 1
                self.rawdat['merge_mass_1'][i] = mass1[argm]
                self.rawdat['merge_mdot_1'][i] = mdot1[argm]
                self.rawdat['merge_lum_1'][i] = lum1[argm]
        if len(self.db_mergers.keys()) > 0:
            self._match_data_to_raw('merge_mass_1', 'merge_mass_2', 'merge_mdot_1', 'merge_mdot_2',
                                'merge_lum_1', 'merge_lum_2', 'tform1','tform2')

    def get_dual_frac(self,bhorbit,minL=1e43,maxD=10,boxsize=25,comove=True):
        if comove:
            tstr = 't_'+str(maxD)+'c'
            fstr = 'frdual_'+str(minL)+'_'+str(maxD)+'c'

        self.rawdat[fstr] = np.ones(len(self.rawdat['ID1']))*-1
        self.rawdat[tstr] = np.ones(len(self.rawdat['ID1']))*-1
        for i in range(len(self.rawdat['ID1'])):
            if i %100==0:
                print float(i)/float(len(self.rawdat['ID1']))*100, '% done'

            x1 = bhorbit.single_BH_data(self.rawdat['ID1'][i],'x')
            x2 = bhorbit.single_BH_data(self.rawdat['ID2'][i],'x')
            y1 = bhorbit.single_BH_data(self.rawdat['ID1'][i],'y')
            y2 = bhorbit.single_BH_data(self.rawdat['ID2'][i],'y')
            z1 = bhorbit.single_BH_data(self.rawdat['ID1'][i],'z')
            z2 = bhorbit.single_BH_data(self.rawdat['ID2'][i],'z')

            lum1 = bhorbit.single_BH_data(self.rawdat['ID1'][i],'lum')
            lum2 = bhorbit.single_BH_data(self.rawdat['ID2'][i],'lum')

            time1 = bhorbit.single_BH_data(self.rawdat['ID1'][i],'time')
            time2 = bhorbit.single_BH_data(self.rawdat['ID2'][i],'time')

            scale1 = bhorbit.single_BH_data(self.rawdat['ID1'][i],'scalefac')
            scale2 = bhorbit.single_BH_data(self.rawdat['ID2'][i],'scalefac')
            if len(time1) == 0 or len(time2) == 0:
                continue

            mint = np.max([time1.min(),time2.min()])
            maxt = time2.max()

            use1 = np.where((time1<=maxt)&(time1>=mint))[0]
            use2 = np.where((time2<=maxt)&(time2>=mint))[0]

            if len(use1) == 0 or len(use2) == 0:
                continue

            if len(use1) != len(use2):
                if len(use1)< len(use2):
                    use1 = np.append(use1,use1[-1])
                    if len(use1) != len(use2):
                        print "SHIIIIIIT"
                else:
                    print "SHIIIIIIIT"
                    print self.rawdat['ID1'][i], self.rawdat['ID2'][i]

            xd = x1[use1]-x2[use2]
            yd = y1[use1]-y2[use2]
            zd = z1[use1]-z2[use2]

            bphys = boxsize*scale1[use1]*1e3
            badx = np.where(xd > bphys/2)[0]
            xd[badx] = -1.0 * (xd[badx]/np.abs(xd[badx])) * \
                              (bphys[badx] - np.abs(xd[badx]))

            bady = np.where(yd > bphys/2)[0]
            yd[bady] = -1.0 * (yd[bady]/np.abs(yd[bady])) * \
                              (bphys[bady] - np.abs(yd[bady]))
            badz = np.where(zd > bphys/2)[0]
            zd[badz] = -1.0 * (zd[badz]/np.abs(zd[badz])) * \
                              (bphys[badz] - np.abs(zd[badz]))

            dist = np.sqrt(xd**2 + yd**2 + zd**2)
            if len(use1) != len(dist) or len(use2) != len(dist):
                print len(use1), len(use2), len(dist)
            if comove:
                dist /= scale1[use1]
            close = np.where(dist<maxD)[0]
            if len(close) > 0:
                try:
                    dt = np.sum(time1[use1[close]] - time1[use1[close]-1])
                except:
                    dt = np.sum(time1[use1[close]] - time1[use1[close]+1])
                self.rawdat[tstr][i] = dt
                dual = np.where((lum1[use1[close]]>minL)&(lum2[use2[close]]>minL))[0]
                self.rawdat[fstr][i] = float(len(dual))/float(len(close))

        self._match_data_to_raw(tstr, fstr)

    def get_halo_merger(self,dbsim,overwrite=False):
        import tangos as db
        if 'dt_hmerger' not in self.data.keys() or overwrite==True:
            self.rawdat['dt_hmerger'] = np.ones(len(self.rawdat['ID1']))*-1
            self.rawdat['dt_hmerger_min'] = np.ones(len(self.rawdat['ID1']))*-1
        nodiff = 0
        badmatch = 0
        for i in range(len(self.rawdat['ID1'])):
            if i%30 == 0:
                print float(i)/float(len(self.rawdat['ID1']))*100, '% done'
            if self.rawdat['dt_hmerger'][i] >= 0 and overwrite == False:
                continue
            bh1 = db.get_halo(str(dbsim.path)+'/%'+str(self.rawdat['step_before'])+'/'+str(self.rawdat['ID1']))
            bh2 = db.get_halo(str(dbsim.path)+'/%'+str(self.rawdat['step_before'])+'/'+str(self.rawdat['ID2']))

            if bh1 is None or bh2 is None:
                self.data['dt_hmerger'][i] = self.data['time'][i] - min(self.rawdat['tform1'][i],self.rawdat['tform2'][i])
                self.data['dt_hmerger_min'][i] = 0
                continue

            time1, hn1 = bh1.reverse_property_cascade('t()', 'host_halo.halo_number()')
            time2, hn2 = bh2.reverse_property_cascade('t()', 'host_halo.halo_number()')
            match1 = np.where(np.in1d(time1,time2))[0]
            match2 = np.where(np.in1d(time2,time1))[0]
            if not np.array_equal(time1[match1],time2[match2]):
                print "WARNING time arrays don't match!"
            if len(match1)==0 or len(match2)==0:
                badmatch += 1
                continue
            diff = np.where(hn1[match1]!=hn2[match2])[0]
            if len(diff)==0:
                nodiff += 1
                self.data['dt_hmerger'][i] = self.data['time'][i] - min(self.rawdat['tform1'][i],self.rawdat['tform2'][i])
                self.data['dt_hmerger_min'][i] = self.data['time'][i] - max(time1.min(),time2.min())
                continue
            th1 = time1[match1[diff[0]]]
            th2 = time2[match2[diff[0]]]
            if diff[0] != 0:
                th1p = time1[match1[diff[0]-1]]
            else:
                th1p = self.data['time'][i]
            if th1 != th2:
                print "WARNING halo merge times not correct"

            self.data['dt_hmerger'][i] = self.data['time'][i] - th1
            if self.data['time'][i] - th1p > 0:
                self.data['dt_hmerger_min'][i] = self.data['time'][i] - th1p
            else:
                self.data['dt_hmerger_min'][i] = 0
        print "finished with ", nodiff, "BHs having never been in different halos"











