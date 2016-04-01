import numpy as np
from .. import Files, cosmology, readcol, plotting
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
        plotting.plt.step(tzbins[0:],np.append(rate,rate[-1], color=color,linestyle=linestyle, label=label, linewidth=lw, where='post')
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
            plotting.plt.xticks([1,2,3,4,5,6,9,11,16,21],['0','1','2','3','4','5','6','8','10','15','20'])
    if ylog is True:
        plotting.plt.yscale('log',base=10)

    if ret_data is True:
        return rate, tzbins,tedges




class mergerCat(object):
    def __init__(self,bhhalocat,bhorbit,mergerfile):
        time, step, ID, IDeat, ratio, kick = readcol(mergerfile,twod=False)
        self.data = {'time':time, 'step':np.ones(len(ID))*-1, 'ID1':ID, 'ID2':IDeat, 'ratio':ratio, 'kick':kick,
                     'mass1':np.ones(len(ID))*-1, 'mass2':np.ones(len(ID))*-1,
                     'mdot1':np.ones(len(ID))*-1,'mdot2':np.ones(len(ID))*-1,
                     'lum1':np.ones(len(ID))*-1,'lum2':np.ones(len(ID))*-1,
                     'snap_prev':np.ones(len(ID))*-1, 'snap_post':np.ones(len(ID))*-1,
                     'tform1':np.ones(len(ID))*-1,'tform2':np.ones(len(ID))*-1}

        for i in range(len(ID)):
            mass1 = bhorbit.single_BH_data(ID[i],'mass')
            mass2 = bhorbit.single_BH_data(IDeat[i],'mass')
            mdot1 = bhorbit.single_BH_data(ID[i],'mdot')
            mdot2 = bhorbit.single_BH_data(IDeat[i],'mdot')
            lum1 = bhorbit.single_BH_data(ID[i],'lum')
            lum2 = bhorbit.single_BH_data(IDeat[i],'lum')
            time1 = bhorbit.single_BH_data(ID[i],'time')
            time2 = bhorbit.single_BH_data(IDeat[i],'time')
            if len(mass2)>0:
                self.data['step'][i] = bhorbit.single_BH_data(IDeat[i],'step')[-1]
                timemerge = time2[-1]
                self.data['mass2'][i] = mass2[-1]
                self.data['lum2'][i] = lum2[-1]
                self.data['mdot2'][i] = mdot2[-1]
            else:
                timemerge = time[i]

            if len(mass1)>0:
                argmerge1 = np.argmin(np.abs(time1-timemerge))
                self.data['mass1'][i] = mass1[argmerge1]
                self.data['mdot1'][i] = mdot1[argmerge1]
                self.data['lum1'][i] = lum1[argmerge1]

            if (i+1)/float(len(ID)) % 0.1 == 0:
                print i/float(len(ID)) * 100, '% done'

        ord_ = np.argsort(self.data['ID1'])
        match1 = np.where(np.in1d(bhorbit.bhiords,self.data['ID1'][ord_]))[0]
        match2 = np.where(np.in1d(self.data['ID1'][ord_],bhorbit.bhiords))[0]
        umatch1, uinv1 = np.unique(self.data['ID1'][ord_[match2]],return_inverse=True)
        print len(match1), len(match2), len(ord_), len(self.data['tform1']), len(uinv1)
        print match1.max(), match2.max(), len(bhorbit.tform), uinv1.max(), ord_.max()
        self.data['tform1'][ord_[match2]] = bhorbit.tform.in_units('Gyr')[match1[uinv1]]

        ord_ = np.argsort(self.data['ID2'])
        match1 = np.where(np.in1d(bhorbit.bhiords,self.data['ID2'][ord_]))[0]
        match2 = np.where(np.in1d(self.data['ID2'][ord_],bhorbit.bhiords))[0]
        umatch2, uinv2 = np.unique(self.data['ID2'][ord_[match2]],return_inverse=True)
        self.data['tform2'][ord_[match2]] = bhorbit.tform.in_units('Gyr')[match1[uinv2]]

        self._prev_snap_slice_1 = {}
        self._prev_snap_slice_inv_1 = {}
        self._prev_snap_slice_2 = {}
        self._prev_snap_slice_inv_2 = {}

        self._post_snap_slice = {}
        self._post_snap_slice_inv = {}


        for ii in range(len(bhhalocat.steps)-1):
            curstep = int(bhhalocat.steps[ii])
            nextstep = int(bhhalocat.steps[ii+1])
            print curstep,nextstep
            curdata = np.where((self.data['step']>curstep)&(self.data['step']<nextstep))[0]

            self.data['snap_prev'][curdata] = curstep
            self.data['snap_post'][curdata] = nextstep

            ord1 = np.argsort(self.data['ID1'][curdata])
            ord2 = np.argsort(self.data['ID2'][curdata])

            match1 = np.where(np.in1d(bhhalocat[str(curstep)].bh['bhid'],self.data['ID1'][curdata[ord1]]))[0]
            match2 = np.where(np.in1d(self.data['ID1'][curdata[ord1]],bhhalocat[str(curstep)].bh['bhid']))[0]
            umatch, uinv = np.unique(self.data['ID1'][curdata[ord1[match2]]],return_inverse=True)
            self._prev_snap_slice_1[str(curstep)] = match1[uinv]
            self._prev_snap_slice_inv_1[str(curstep)] = curdata[ord1[match2]]

            match1 = np.where(np.in1d(bhhalocat[str(curstep)].bh['bhid'],self.data['ID2'][curdata[ord2]]))[0]
            match2 = np.where(np.in1d(self.data['ID2'][curdata[ord2]],bhhalocat[str(curstep)].bh['bhid']))[0]
            umatch, uinv = np.unique(self.data['ID2'][curdata[ord2[match2]]],return_inverse=True)
            self._prev_snap_slice_2[str(curstep)] = match1[uinv]
            self._prev_snap_slice_inv_2[str(curstep)] = curdata[ord2[match2]]

            match1 = np.where(np.in1d(bhhalocat[str(nextstep)].bh['bhid'],self.data['ID1'][curdata[ord1]]))[0]
            match2 = np.where(np.in1d(self.data['ID1'][curdata[ord1]],bhhalocat[str(nextstep)].bh['bhid']))[0]
            umatch, uinv = np.unique(self.data['ID1'][curdata[ord1[match2]]],return_inverse=True)
            self._post_snap_slice[str(nextstep)] = match1[uinv]
            self._post_snap_slice_inv[str(nextstep)] = curdata[ord1[match2]]

    def get_snap_info(self,key,bhhalocat):
        self.data['prev_'+key+'1'] = np.ones(len(self.data['ID1']))*-1
        self.data['prev_'+key+'2'] = np.ones(len(self.data['ID1']))*-1
        self.data['post_'+key] = np.ones(len(self.data['ID1']))*-1
        for ii in range(len(bhhalocat.steps)-1):
            print bhhalocat.steps[ii]
            curstep = int(bhhalocat.steps[ii])
            nextstep = int(bhhalocat.steps[ii+1])
            if key in bhhalocat[str(curstep)].bh.keys():
                data = bhhalocat[str(curstep)].bh[key]
                datanext = bhhalocat[str(nextstep)].bh[key]
            else:
                data = bhhalocat[str(curstep)].halo_properties[key]
                datanext = bhhalocat[str(nextstep)].halo_properties[key]
            if len(self._prev_snap_slice_1[str(curstep)])>0:
                self.data['prev_'+key+'1'][self._prev_snap_slice_inv_1[str(curstep)]] = data[self._prev_snap_slice_1[str(curstep)]]
                self.data['prev_'+key+'2'][self._prev_snap_slice_inv_2[str(curstep)]] = data[self._prev_snap_slice_2[str(curstep)]]
            if len(self._post_snap_slice_inv[str(nextstep)])>0:
                self.data['post_'+key][self._post_snap_slice_inv[str(nextstep)]] = datanext[self._post_snap_slice[str(nextstep)]]


    def __getitem__(self,item):
        return self.data[item]

    def keys(self):
        return self.data.keys()
