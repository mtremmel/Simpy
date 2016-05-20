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
        vol_weights = np.ones(len(time))*vol_weights
    print type(bins)
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

# GW calculations based on Salcido+ 2015, Flanagen + Hughes 1998, and others
def freq_qnr(Mtot, a=0.95):
    Mtot *= util.M_sun_g
    return util.c**3 * (1 - 0.63*(1-a)**(3./10.)) / (2 * np.pi * util.G * Mtot)

def freq_merger(Mtot):
    Mtot *= util.M_sun_g
    return util.c**3 * 0.02 / (util.G * Mtot)

def dEdf(M1, M2, eps=0.1, a=0.95):
    fqnr = freq_qnr(M1+M2,a)
    fm = freq_merger(M1+M2)
    M1 *= util.M_sun_g
    M2 *= util.M_sun_g
    F = (4.*M1*M2)**2/(M1+M2)**4
    return util.c**2*(M1+M2)*F*eps/(fqnr-fm)

def strain(M1, M2, z, omegaM, omegaL, h0, eps=0.1, a=0.95):
    dL = cosmology.lum_distance(z,omegaM, omegaL, h0)
    dL *= 3.086e+24
    Ef = dEdf(M1, M2, eps, a)
    return np.sqrt(2.*util.G/util.c**3) * ((1+z)/(np.pi*dL)) * np.sqrt(Ef)

#eLisa limits from Salcido+2015, taken from Amaro-Seoane+ 2013
def eLisaLimit(fobs):
    def sacc(fobs):
        return 1.37e-32 * (1+1e-4/fobs)*(fobs)**-4
    ssn = 5.25e-23
    somn = 6.28e-23
    L = 1.0e9
    def sn(fobs):
        return 20./3. * (4.*sacc(fobs)+ssn+somn)/L**2 * (1+(fobs/(0.41*(util.c/(1.0e2*2.0*L))))**2)

    return np.sqrt(sn(fobs)*fobs)




class mergerCat(object):
    def __init__(self, dbsim, simname, properties=[], bhorbit=None):
        proplist = ['halo_number()', 'BH_merger_next.halo_number()', 'host_halo.halo_number()',
                    'BH_merger_next.host_halo.halo_number()', 'BH_merger_next.earlier(1).host_halo.halo_number()']
        for prop in properties:
            proplist.append(prop)
            proplist.append('BH_merger_next.'+prop)
            proplist.append('BH_merger_next.earlier(1).'+prop)
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

        uIDeat, cnt = np.unique(self.rawdat['ID2'], return_counts=True)
        if len(np.where(cnt>1)[0])>0:
            print "SHIIIIIIT"

        self.data = {'ID1':[], 'ID2':[], 'ratio':[], 'kick':[], 'step':[],
                    'time':[], 'tsnap_prev':[], 'tsnap_after':[], 'snap_prev':[], 'snap_after':[],
                    'host_N_pre_1':[], 'host_N_pre_2':[], 'host_N_post':[]}

        for p in properties:
            self.data[p+'_pre_1'] = []
            self.data[p+'_pre_2'] = []
            self.data[p+'_post'] = []


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

            self.data['ID1'].extend(data[1])
            self.data['ID2'].extend(data[0])

            tt1, tt2 = np.unique(data[0],return_counts=True)
            if len(np.where(tt2>1)[0])>0:
                print "Double counted IDs: ", tt1[(tt2>1)]
                raise RuntimeError("ERROR double counted IDeat in database analysis")

            self.data['host_N_pre_1'].extend(data[4])
            self.data['host_N_pre_2'].extend(data[2])
            self.data['host_N_post'].extend(data[3])

            self.data['snap_after'].extend(int(stepnumA) * np.ones(len(data[1])).astype(np.int))
            self.data['snap_prev'].extend(int(stepnum) * np.ones(len(data[1])).astype(np.int))
            self.data['tsnap_after'].extend(step.next.time_gyr * np.ones(len(data[1])))
            self.data['tsnap_prev'].extend(step.time_gyr * np.ones(len(data[1])))

            index = 5
            for i in range(len(properties)):
                self.data[properties[i]+'_pre_1'].extend(data[index+2])
                self.data[properties[i]+'_pre_2'].extend(data[index])
                self.data[properties[i]+'_post'].extend(data[index+1])
                index += 3

        for key in self.data.keys():
            self.data[key] = np.array(self.data[key])

        self._match_data_to_raw('ratio', 'kick', 'time','step')

    def __getitem__(self, item):
        return self.data[item]

    def keys(self):
        return self.data.keys()

    def _match_data_to_raw(self,*properties):
        ordee = np.argsort(self.data['ID2'])
        match = np.where(np.in1d(self.data['ID2'][ordee],self.rawdat['ID2']))[0]
        match2 = np.where(np.in1d(self.rawdat['ID2'],self.data['ID2'][ordee]))[0]

        if len(match) != len(match2):
            print self.data['ID2'][ordee], self.rawdat['ID2']
            print match, match2
            raise RuntimeError("ERROR match not returning same number of elements")

        if len(np.where(np.equal(self.data['ID2'][ordee[match]],self.rawdat['ID2'][match2]) is False)[0])>0:
            raise RuntimeError("ERROR something wrong with array matching")


        for p in properties:
            self.data[p] = np.ones(len(self.data['ID2'])) * -1
            self.data[p][ordee[match]] = self.rawdat[p][match2]


    def calc_GW_emission(self, sim, eps=0.1, a=0.95):
        self.rawdat['GW_freq_merge'] = np.ones(len(self.rawdat['ID1']))*-1
        self.rawdat['GW_strain_merge'] = np.ones(len(self.rawdat['ID1']))*-1

        if 'merge_mass_1' not in self.rawdat.keys():
            print "ERROR need to find final masses for BHs"
            return

        ok = np.where((self.rawdat['merge_mass_1'] > 0)&(self.rawdat['merge_mass_2']>0))[0]
        if "redshift" in self.rawdat.keys():
            z = self.rawdat['redshift']
        else:
            z = pynbody.analysis.cosmology.redshift(sim, self.rawdat['time'])
            self.rawdat['redshift'] = z
        hc = strain(self.rawdat['merge_mass_1'][ok], self.rawdat['merge_mass_2'][ok], self.rawdat['redshift'][ok],
                    sim.properties['omegaM0'], sim.properties['omegaL0'], sim.properties['h'], eps=eps, a=a)
        fm = freq_merger(self.rawdat['merge_mass_1']+self.rawdat['merge_mass_2'])
        self.rawdat['GW_freq_merge'][ok] = fm
        self.rawdat['GW_strain_merge'][ok] = hc

        self._match_data_to_raw('GW_freq_merge', 'GW_strain_merge')



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
                self.rawdat['tform1'][i] = bhorbit.tform.in_units('Gyr')[oo]
            if len(time2)>0:
                oo = np.where(bhorbit.bhiords == self.rawdat['ID2'][i])[0]
                oo = oo[0]
                self.rawdat['tform2'][i] = bhorbit.tform.in_units('Gyr')[oo]


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

        self._match_data_to_raw('merge_mass_1', 'merge_mass_2', 'merge_mdot_1', 'merge_mdot_2',
                                'merge_lum_1', 'merge_lum_2', 'tform1','tform2')

    def get_dual_frac(self,bhorbit,minL=1e43,maxD=10):
        self.rawdat['frdual_'+str(minL)+'_'+str(maxD)] = np.ones(len(self.rawdat['ID1']))*-1
        self.rawdat['t_'+str(maxD)] = np.ones(len(self.rawdat['ID1']))*-1
        for i in range(len(self.rawdat['ID1'])):
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

            mint = np.min([time1.min(),time2.min()])
            maxt = self.rawdat['time'][i]

            use1 = np.where((time1<maxt)&(time1>mint))[0]
            use2 = np.where((time2<maxt)&(time2>mint))[0]

            if len(use1) != len(use2):
                print "Shit!!!!!"

            dist = np.sqrt((x1[use1]-x2[use2])**2 + (y1[use1]-y2[use2])**2 + (z1[use1]-z2[use2])**2)
            close = np.where(dist<maxD)[0]
            dt = time1[use1[close[-1]]] - time1[use1[close[0]]]
            self.rawdat['t_'+str(maxD)][i] = dt
            dual = np.where((lum1[use1[close]]>minL)&(lum2[use2[close]]>minL))
            self.rawdat['frdual_'+str(minL)+'_'+str(maxD)][i] = float(len(dual))/float(len(close))

        self._match_data_to_raw('t_'+str(maxD), 'frdual_'+str(minL)+'_'+str(maxD))












