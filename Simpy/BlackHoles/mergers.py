import numpy as np
from .. import Files, cosmology, plotting, util
from ..readcol import readcol
import os
import pynbody


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

    a, b, ID, c, IDeat, d, time, e, f, kick, g, h, mratio = readcol.readcol(mergename, twod=False)
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

    dtz = (tzrange[1]-tzrange[0])/float(bins)
    tzbins = np.arange(tzrange[0],tzrange[1]+dtz,dtz)
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
            plotting.plt.step(tzbins[0:-1],rate, color=color,linestyle=linestyle, label=label, linewidth=lw, where='post')

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
    def __init__(self, dbsim, simname, properties=[]):
        proplist = ['halo_number()', 'later(1).halo_number()', 'host_halo.halo_number()', 'later(1).halo_number()']
        for prop in properties:
            proplist.append(prop)
            proplist.append('later(1).'+prop)
        mergerfile = simname+'.mergers'
        print "reading .mergers file..."
        time, step, ID, IDeat, ratio, kick = readcol(mergerfile,twod=False)
        ID = ID.astype(np.int64)
        IDeat = IDeat.astype(np.int64)
        print "checking for bad IDs..."
        bad = np.where(ID<0)
        if len(bad)>0:
            ID[bad] = 2*2147483648 + ID[bad]
        bad2 = np.where(IDeat<0)
        if len(bad2)>0:
            IDeat[bad2] = 2*2147483648 + IDeat[bad2]

        self.rawdat = {'time':time, 'ID1':ID, 'ID2':IDeat, 'ratio':ratio, 'kick':kick}
        ordr = np.argsort(self.rawdat['ID2'])
        util.cutdict(self.rawdat,ordr)

        self.data = {'ID1':[], 'ID2':[], 'ratio':[], 'kick':[], 'step':[],
                    'tmerge':[], 'tstep_prev':[], 'tstep_after':[],
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
                continue
            bhid = data[0]
            bhid_next = data[1]
            good = np.where(np.in1d(bhid_next,bhid))[0]

            bhid = bhid[good]
            bhid_next = bhid_next[good]
            for i in range(len(data[2:])):
                data[i+2] = data[i+2][good]

            host_n = data[2]
            host_n_next = data[3]

            ordd = np.argsort(bhid)

            ubhid, inv, cnt = np.unique(bhid_next,return_counts=True, return_inverse=True)
            mm = np.where(cnt>1)
            if len(mm)==0:
                print "No BH mergers found this step"
                continue
            self.nmergers.append((cnt[mm]-1).sum())
            self.steptimes.append(step.time_gyr)
            print self.nmergers, self.steptimes

            eat = np.where((bhid_next != bhid)&(cnt[inv]>1))
            bheat = bhid[eat]
            bhmain, inv = np.unique(bhid_next[eat],return_inverse=True)
            match = np.where(np.in1d(bhid[ordd],bhmain))[0]
            print np.where(np.equal(bhid[ordd[match]],bhmain)==False)
            print np.where(np.equal(bhid[ordd[match[inv]]],bhmain[inv]))

            self.data['host_N_pre_1'].append(host_n[match[inv]])
            self.data['host_N_pre_2'].append(host_n[eat])
            self.data['host_N_post'].append(host_n_next[eat])

            index = 4
            for i in range(len(properties)):
                self.data[properties[i]+'_pre_1'].append(data[index][match[inv]])
                self.data[properties[i]+'_pre_2'].append(data[index][eat])
                self.data[properties[i]+'_post'].append(data[index+1][eat])
                index += 2

        ordee = np.argsort(self.data['ID2'])
        match = np.where(np.in1d(self.data['ID2'][ordee],self.rawdat['ID2']))[0]
        match2 = np.where(np.in1d(self.rawdat['ID2'],self.data['ID2'][ordee]))[0]

        self.data['ratio'] = np.zeros(len(self.data['ID2']))
        self.data['kick'] = np.zeros(len(self.data['ID2']))
        self.data['time'] = np.zeros(len(self.data['ID2']))
        self.data['step'] = np.zeros(len(self.data['ID2']))

        self.data['ratio'][ordee[match]] = self.rawdat['ratio'][match2]
        self.data['kick'][ordee[match]] = self.rawdat['kick'][match2]
        self.data['time'][ordee[match]] = self.rawdat['time'][match2]
        self.data['step'][ordee[match]] = self.rawdat['step'][match2]

