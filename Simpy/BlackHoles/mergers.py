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

    simtime = cosmology.getTime(s.properties['a'] ** -1 - 1, s)
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

def plt_merger_rates(time,sim, style, vol_weights=1./25.**3, bins=50,
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
        plotting.plt.step(tzbins[0:-1],rate, style, label=label, linewidth=lw, where='post')
        plotting.plt.xlabel('Time (Gyr)')
    else:
        if xlog is False:
            plotting.plt.step(tzbins[0:-1],rate, style, label=label, linewidth=lw, where='post')
            plotting.plt.xlabel('Redshift')

        else:
            plotting.plt.step(tzbins[0:-1]+1,rate, style, label=label, linewidth=lw, where='post')
            plotting.plt.xlabel('z + 1')


    if xlog is True:
        plotting.plt.xscale('log',base=10)
    if ylog is True:
        plotting.plt.yscale('log',base=10)

    if ret_data is True:
        return rate, tzbins




class mergerCat(object):
    def __init(self,bhhalocat,bhorbit,mergerfile):
        time, step, ID, IDeat, ratio, kick = readcol.readcol(mergerfile,twod=False)
        self.data = {'time':time, 'step':step, 'ID1':ID, 'ID2':IDeat, 'ratio':ratio, 'kcik':kick,
                     'mass1':np.ones(len(ID))*-1, 'mass2':np.ones(len(ID))*-1,
                     'mdot1':np.ones(len(ID))*-1,'mdot2':np.ones(len(ID))*-1,
                     'lum1':np.ones(len(ID))*-1,'lum2':np.ones(len(ID))*-1,
                     'snap_prev':np.ones(len(ID))*-1, 'snap_post':np.ones(len(ID))*-1,
                     'tform1':np.ones(len(ID))*-1,'tform2':np.ones(len(ID))*-1}

        stepfl = np.floor(step)
        sord = np.argsort(stepfl)
        ustepfl, ind = np.unique(stepfl[sord],return_index=True)
        self.step_slice = {}
        for i in range(len(ustepfl) - 1):
            ss = sord[ind[i]:ind[i + 1]]
            self.step_slice[str(ustepfl[i])] = ss
        ss = sord[ind[i + 1]:]
        self.step_slice[str(ustepfl[i+1])] = ss

        print "gathering black hole data from orbit file..."

        for st in self.step_slice.keys():
            bhids = bhorbit.single_step_data(int(st),'iord')
            ord_orbit = np.argsort(bhids)
            masses = bhorbit.single_step_data(int(st),'mass')
            masses = masses[ord_orbit]
            mdot = bhorbit.single_step_data(int(st),'mdotmean')
            mdot = mdot[ord_orbit]
            lum = bhorbit.single_step_data(int(st),'lum')
            lum = lum[ord_orbit]
            stepid1 = self.data['ID1'][self.step_slice[st]]
            stepid2 = self.data['ID2'][self.step_slice[st]]
            ord_ = np.argsort(stepid1)
            ord2_ = np.argsort(stepid2)
            omatch = np.where(np.in1d(bhids,stepid1[ord_]))[0]
            omatch2 = np.where(np.in1d(stepid1[ord_],bhids))[0]
            self.data['mass1'][self.step_slice[st]][ord_][omatch] = masses[omatch2]
            self.data['mdot1'][self.step_slice[st]][ord_][omatch] = mdot[omatch2]
            self.data['lum1'][self.step_slice[st]][ord_][omatch] = lum[omatch2]
            omatch = np.where(np.in1d(bhids,stepid2[ord2_]))[0]
            omatch2 = np.where(np.in1d(stepid2[ord2_],bhids))[0]
            self.data['mass2'][self.step_slice[st]][ord2_][omatch] = masses[omatch2]
            self.data['mdot2'][self.step_slice[st]][ord2_][omatch] = mdot[omatch2]
            self.data['lum2'][self.step_slice[st]][ord2_][omatch] = lum[omatch2]

        ord_ = np.argsort(self.data['ID1'])
        match1 = np.where(np.in1d(bhorbit.bhiords,self.data['ID1'][ord_]))[0]
        match2 = np.where(np.in1d(self.data['ID1'][ord_],bhorbit.bhiords))[0]
        self.data['tform1'][match2] = bhorbit.tform[match1]

        ord_ = np.argsort(self.data['ID2'])
        match1 = np.where(np.in1d(bhorbit.bhiords,self.data['ID2'][ord_]))[0]
        match2 = np.where(np.in1d(self.data['ID2'][ord_],bhorbit.bhiords))[0]
        self.data['tform2'][match2] = bhorbit.tform[match1]

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

            self._prev_snap_slice_1[str(curstep)] = \
                np.where(np.in1d(bhhalocat[str(curstep)].bh['bhid'],self.data['ID1'][curdata[ord1]]))[0]
            self._prev_snap_slice_inv_1[str(curstep)] = \
                curdata[ord1[np.where(np.in1d(self.data['ID1'][curdata[ord1]],bhhalocat[str(curstep)].bh['bhid'])[0])]]

            self._prev_snap_slice_2[str(curstep)] = \
                np.where(np.in1d(bhhalocat[str(curstep)].bh['bhid'],self.data['ID2'][curdata[ord2]]))[0]
            self._prev_snap_slice_inv_2[str(curstep)] = \
                curdata[ord2[np.where(np.in1d(self.data['ID2'][curdata[ord2]],bhhalocat[str(curstep)].bh['bhid'])[0])]]

            self._post_snap_slice[str(nextstep)] = \
                np.where(np.in1d(bhhalocat[str(nextstep)].bh['bhid'],self.data['ID1'][curdata[ord1]]))[0]
            self._post_snap_slice_inv[str(nextstep)] = \
                curdata[ord1[np.where(np.in1d(self.data['ID1'][curdata[ord1]],bhhalocat[str(nextstep)].bh['bhid'])[0])]]



    def get_snap_info(self,key,bhhalocat):
        self.data['prev_'+key+'1'] = np.ones(len(self.data['ID1']))*-1
        self.data['prev_'+key+'2'] = np.ones(len(self.data['ID1']))*-1
        self.data['post_'+key] = np.ones(len(self.data['ID1']))*-1
        for ii in range(len(bhhalocat.steps)-1):
            print bhhalocat.steps[ii]
            curstep = int(bhhalocat.steps[ii])
            nextstep = int(bhhalocat.steps[ii+1])
            if key in bhhalocat[bhhalocat.steps[ii]].bh.keys():
                data = bhhalocat[curstep].bh[key]
                datanext = bhhalocat[nextstep].bh[key]
            else:
                data = bhhalocat[curstep].halo_properties[key]
                datanext = bhhalocat[nextstep].halo_properties[key]
            self.data['prev_'+key+'1'][self._prev_snap_slice_inv_1[curstep]] = data[self._prev_snap_slice_1[curstep]]
            self.data['prev_'+key+'2'][self._prev_snap_slice_inv_2[curstep]] = data[self._prev_snap_slice_2[curstep]]
            self.data['post_'+key][self._post_snap_slice_inv[nextstep]] = datanext[self._post_snap_slice[nextstep]]




