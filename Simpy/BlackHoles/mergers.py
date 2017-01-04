import numpy as np
from .. import Files, cosmology, plotting, util
from ..readcol import readcol
import os
import pynbody
import re

def findmergers(outname, diagname='*out*'):
    os.system("awk '/BHSink/ && /Merge/ && /eating/' " + diagname + " > " + outname)
    return


def reducedata(simname, RetData=False, outname='*out*', mergename='BHmerge.txt', NCHILADA=True, out_ext='.mergers'):
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
    np.savetxt(simname + out_ext, np.column_stack(tofile), fmt=['%f', '%f', '%d', '%d', '%f', '%f'])
    if RetData == True:
        return output
    else:
        return

# def get_complete_prog_list(bhorbit, bhid, tmax):
#     target, = np.where(bhorbit.bhiords == bhid)
#     bhorbit.getprogbhs()
#     idlist = np.array([])
#     if len(bhorbit.prog['iord'][target])==0:
#         return np.array([])
#     idnew = np.array(bhorbit.prog['iord'][target])[(np.array(bhorbit.prog['time'][target])<tmax)]
#     idlist = np.append(idlist,idnew)
#
#     deep = 0
#     while len(idnew) > 0:
#         deep += 1
#         idnext = np.array([])
#         for bid in idnew:
#             newtarget, = np.where(bhorbit.bhiords==bid)
#             if len(newtarget)>1:
#                 print "Warning! multiple matches in bhiords found for ", bid
#             if len(newtarget)==0:
#                 print str(bid)+" not found in orbit object! moving on..."
#                 continue
#             newtarget = newtarget[0]
#             newpart = np.array(bhorbit.prog['iord'][newtarget])
#             idnext = np.append(idnext,newpart)
#         idnew = idnext
#         idlist = np.append(idlist,idnew)
#     print "finished with ", deep, "steps\n"
#     return idlist

def get_complete_prog_list(bhmergers,bhid,tmax,useonly=None, return_details=False):
    if useonly is None:
        useonly = np.arange(len(bhmergers['ID1']))
    match, = np.where((bhmergers['ID1'][useonly]==bhid)&(bhmergers['time'][useonly]<=tmax))
    if len(match)==0:
        if return_details is False:
            return np.array([])
        else:
            return np.array([]), np.array([]), np.array([])
    idnew = np.copy(bhmergers['ID2'][useonly[match]])
    idlist = np.copy(idnew)
    deep = 0
    if return_details is True:
        massnew = np.copy(bhmergers['merge_mass_2'][useonly[match]])
        masslist = np.copy(massnew)
        timenew = np.copy(bhmergers['time'][useonly[match]])
        timelist = np.copy(timenew)
    while len(idnew)>0:
        deep +=1
        idnext = np.array([])
        massnext = np.array([])
        timenext = np.array([])
        for eid in idnew:
            match, = np.where(bhmergers['ID1'][useonly]==eid)
            if len(match)>0:
                idnext = np.append(idnext,bhmergers['ID2'][useonly[match]])
                if return_details is True:
                    massnext = np.append(massnext,bhmergers['merge_mass_2'][useonly[match]])
                    timenext = np.append(timenext,bhmergers['time'][useonly[match]])
        idnew = idnext
        idlist = np.append(idlist,idnew)
        if return_details is True:
            masslist = np.append(masslist,massnext)
            timelist = np.append(timelist,timenext)
    if return_details is True:
        return idlist, masslist, timelist
    else:
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


def mass_binned_counts(redshift,Mvir,hmf,s,weights=None,zrange=[0,10],dz=0.5,tnorm=True, logz=True):
    if logz is True:
        lzbins = np.arange(zrange[0],zrange[1]+dz,dz)
        zbins = 10**lzbins
        zmid = zbins[0:-1] + (zbins[1:]-zbins[0:-1])*0.5
    else:
        zbins = np.arange(zrange[0],zrange[1]+dz,dz)
        zmid = zbins[0:-1]+dz/2.
    Mbins = np.arange(hmf.minlm,hmf.maxlm+hmf.delta,hmf.delta)

    mbin_counts = np.zeros((len(zbins)-1,len(Mbins)-1))
    mbin_total = np.zeros((len(zbins)-1,len(Mbins)-1))

    for i in range(len(zbins)-1):
        o, = np.where((redshift>=zbins[i])&(redshift<zbins[i+1]))
        n,b = np.histogram(np.log10(Mvir[o]),bins=Mbins,weights=weights[o])
        mbin_counts[i,:] = n
        zind = np.where((hmf.z >= zbins[i])&(hmf.z < zbins[i+1]))[0]
        if len(zind) > 0:
            mbin_total[i,:] = np.sum(hmf.nhalos[zind,:],axis=0)/float(len(zind))
        if len(zind) == 0:
            z1 = np.where(hmf.z > zmid[i])[0]
            z2 = np.where(hmf.z < zmid[i])[0]
            if len(z1)== 0 or len(z2) == 0:
                print "ERROR make sure the z bins are within the available snapshot z range"
                raise RuntimeError
            z1 = z1[-1]
            z2 = z2[0]

            mbin_total[i,:] = hmf.nhalos[z1] + \
                (hmf.nhalos[z2]-hmf.nhalos[z1])/(hmf.z[z2]-hmf.z[z1]) * (zmid[i]-hmf.z[z1])
        bad = np.where((mbin_total[i,:]==0)&(mbin_counts[i,:]>0))[0]
        while len(bad)>0:
            for bi in bad:
                mbin_counts[i,bi-1] += mbin_counts[i,bi]
                mbin_counts[i,bi] = 0
            bad = np.where((mbin_total[i,:]==0)&(mbin_counts[i,:]>0))[0]

    if tnorm is True:
        tedges = pynbody.array.SimArray([cosmology.getTime(z,s) for z in zbins],'Gyr')
        dt = np.abs(tedges[0:-1]-tedges[1:])
    else:
        dt = None

    return mbin_counts, mbin_total, dt, dz, zmid


def raw_mass_bin_to_rates(cnt, tot, z, norm, hmf, ptype='obs'):
    f = open('files.list','r')
    s = pynbody.load(f.readline().strip('\n'))
    frac = np.nan_to_num(cnt/tot)
    n = np.zeros_like(frac)
    for i in range(len(z)):
        ind = np.argmin(np.abs(z[i]-hmf.z))
        n[i,0:-1] = frac[i,0:-1]*hmf.hmf['phi'][ind]*0.3

    if ptype == 'obs':
        nall = n.sum(axis=1)/norm
        nobs = cosmology.event_count(nall,z,s.properties['omegaM0'], s.properties['omegaL0'], s.properties['h'])
        return nobs
    if ptype == 'rate':
        return n.sum(axis=1)/norm


def combine_merger_data(z,mh,hmf,s, weights=None,zrange=[0,10],dz=0.5,tnorm=True,rel_weights=[8,1], logz=True, ptype='obs'):
    todo = len(z)
    if weights is None:
        weights = []
        for i in range(todo):
            weights.append(np.ones(len(z[i])))
    cnt, tot, dt, dz, zmid = mass_binned_counts(z[0], mh[0], hmf, s, weights[0], zrange, dz, tnorm, logz)
    cnt *= rel_weights[0]
    tot *= rel_weights[0]
    for i in range(todo-1):
        cn, tn, dt, dz, zmid = mass_binned_counts(z[i+1], mh[i+1], hmf, s, weights[i+1], zrange, dz, tnorm, logz)
        cnt += cn*rel_weights[i+1]
        tot += tn*rel_weights[i+1]
    if ptype=='obs':
        nobs = raw_mass_bin_to_rates(cnt, tot, zmid, dz, hmf,ptype=ptype)
    if ptype=='rate':
        nobs = raw_mass_bin_to_rates(cnt, tot, zmid, dt*1e9, hmf,ptype=ptype)
    if ptype=='zdist':
        nobs = raw_mass_bin_to_rates(cnt, tot, zmid, dz, hmf,ptype='rate')
    return nobs, zmid

def cal_weights(z, mh, hmf, s, hmf_comp = None, rel_weight=None):
    weights = np.ones_like(z)*-1
    for i in range(len(hmf.z)-1):
        for j in range(len(hmf.mbins)-2):
            target = np.where((z<hmf.z[i])&(z>=hmf.z[i+1])&(np.log10(mh)>hmf.mbins[j])&(np.log10(mh)<=hmf.mbins[j+1]))[0]
            #if hmf_comp is None:
            if len(target)==0: continue
            nh_hmf = hmf.nhalos[i+1,j]
            jj = j
            while nh_hmf == 0:
                jj = jj -1
                nh_hmf = hmf.nhalos[i+1,jj]
            weights[target] = hmf.hmf['phi'][i+1][jj]*0.3/nh_hmf
    return weights
  #          else:
   #             norm  = hmf.nhalos[i+1,j]
    #            cnt = 0
     #           for hmfc in hmf_comp:
     #               dir = np.where(hmfc.z==hmf.z[i+1])[0]
     ###               if len(dir)>0:
      #                  norm += hmfc.nhalos[dir[0],j]*rel_weight[cnt]
      #              else:
      #                  ocomp_l = np.where(hmf_comp.z<=hmf.z[i+1])[0]
      ##                  ocomp_h = np.where(hmf_comp.z>hmf.z[i+1])[0]
       #                 if len(ocomp_h)==0:
       #                     norm += hmfc.nhalos[ocomp_l[0],j]*rel_weight[cnt]
       #                 else:
       #                     norm += hmfc.nhalos[ocomp_l[0],j] + \
       #                             (hmfc.nhalos[ocomp_h[0],j]-hmfc.nhalos[ocomp_l[0],j])/((hmfc.z[ocomp_l[0]]-hmfc.z[ocomp_h[0]])) *

def calc_nobs(z, m1, m2, mh, hmf, s, weights=None, rel_weights=[1],
              msum_range=None, mmin_range=None, mh_range=None, ratio_range=None, zrange=[0,10], dz = 0.5, logz=True, ptype='obs'):
    if type(z) != list or type(m1) != list or type(m2) != list or type(mh) != list:
        print "expecting lists of arrays!"
        raise ValueError
    zuse_l = []
    m1use_l = []
    m2use_l = []
    mhuse_l = []
    w_l = []
    for i in range(len(z)):
        zuse = np.copy(z[i])
        m1use = np.copy(m1[i])
        m2use = np.copy(m2[i])
        mhuse = np.copy(mh[i])
        if weights is not None:
            wuse = np.copy(weights[i])
        if msum_range is not None:
            ok = np.where((m1use+m2use>=msum_range[0])&(m1use+m2use<msum_range[1]))[0]
            zuse = zuse[ok]
            m1use = m1use[ok]
            m2use = m2use[ok]
            mhuse = mhuse[ok]
            if weights:
                wuse = wuse[ok]
        if mmin_range is not None:
            ok = np.where((np.minimum(m1use,m2use)>=mmin_range[0])&(np.minimum(m1use,m2use)<mmin_range[1]))[0]
            zuse = zuse[ok]
            m1use = m1use[ok]
            m2use = m2use[ok]
            mhuse = mhuse[ok]
            if weights:
                wuse = wuse[ok]
        if mh_range is not None:
            ok = np.where((mhuse>=mh_range[0])&(mhuse<mh_range[1]))[0]
            zuse = zuse[ok]
            m1use = m1use[ok]
            m2use = m2use[ok]
            mhuse = mhuse[ok]
            if weights:
                wuse = wuse[ok]
        if ratio_range is not None:
            ratio = np.minimum(m1use,m2use)/np.maximum(m1use,m2use)
            ok = np.where((ratio>=ratio_range[0])&(ratio<ratio_range[1]))[0]
            zuse = zuse[ok]
            m1use = m1use[ok]
            m2use = m2use[ok]
            mhuse = mhuse[ok]
            if weights:
                wuse = wuse[ok]
        if weights:
            w_l.append(wuse)
        zuse_l.append(zuse)
        m1use_l.append(m1use)
        m2use_l.append(m2use)
        mhuse_l.append(mhuse)
    if weights is None:
        if len(z)==1:
            cnt, tot, dt, dz, zmid = mass_binned_counts(zuse_l[0], mhuse_l[0], hmf, s, np.ones(len(zuse_l[0])), zrange, dz, True, logz)
            if ptype=='obs':
                return raw_mass_bin_to_rates(cnt, tot, zmid, dz, hmf,ptype=ptype), zmid
            if ptype=='rate':
                return raw_mass_bin_to_rates(cnt, tot, zmid, dt*1e9, hmf,ptype=ptype), zmid
            if ptype=='zdist':
                return raw_mass_bin_to_rates(cnt, tot, zmid, dz, hmf,ptype='rate'), zmid
        else:
            return combine_merger_data(z,mh,hmf,s, weights=None,zrange=zrange,dz=dz,rel_weights=rel_weights, logz=logz, ptype=ptype)
    else:
        if logz is True:
            zbins = 10**np.arange(zrange[0],zrange[1]+dz,dz)
        else:
            zbins = np.arange(zrange[0],zrange[1]+dz,dz)
        zmid = zbins[0:-1]+0.5*(zbins[1:]-zbins[0:-1])
        n, zbins = np.histogram(zuse_l[0],weights=w_l[0],bins=zbins)
        n *= rel_weights[0]
        for i in range(len(zuse_l)-1):
            nn, zbinsn = np.histogram(zuse_l[i+1],weights=w_l[i+1],bins=zbins)
            n += nn*rel_weights[i+1]
        if ptype=='obs':
            return cosmology.event_count(n,zmid,s.properties['omegaM0'], s.properties['omegaL0'], s.properties['h'])/dz, zmid
        if ptype=='rate':
            tedges = pynbody.array.SimArray([cosmology.getTime(z,s) for z in zbins],'Gyr')
            dt = np.abs(tedges[0:-1]-tedges[1:]) * 1e9
            return n/dt,zmid
        if ptype=='zdist':
            return n/dz, zmid



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
        import tangos as db
        self.rawdat['snap_after'] = np.zeros(len(self.rawdat['ID1'])).astype('S100')
        self.rawdat['snap_before'] = np.zeros(len(self.rawdat['ID1'])).astype('S100')
        self.rawdat['sat_flag'] = np.zeros(len(self.rawdat['ID1']))
        zsteps = []
        for step in dbsim.timesteps:
            zsteps.append(step.redshift)
        zsteps = np.array(zsteps)
        for p in halo_props:
            self.rawdat[p] = np.ones(len(self.rawdat['ID1'])) * -1
        for i in range(len(self.rawdat['ID1'])):
            if i%100 == 0: print float(i)/len(self.rawdat['ID1']) * 100, '% done'
            id1 = self.rawdat['ID1'][i]
            id2 = self.rawdat['ID2'][i]
            if self.rawdat['redshift'][i]<zsteps.min():
                continue
            ind = np.where(self.rawdat['redshift'][i]>=zsteps)[0][0]
            self.rawdat['snap_after'][i] = dbsim.timesteps[ind].extension
            if ind > 0:
                self.rawdat['snap_before'][i] = dbsim.timesteps[ind-1].extension
            bh = db.get_halo(str(dbsim.timesteps[ind].path)+'/1.'+str(id1))
            if bh is None:
                if ind > 0:
                    bh = db.get_halo(str(dbsim.timesteps[ind-1].path)+'/1.'+str(id1))
                    if bh is not None:
                        bh = bh.next
            if bh is None:
                if ind > 0:
                    bh = db.get_halo(str(dbsim.timesteps[ind-1].path)+'/1.'+str(id2))
                    if bh is not None:
                        bh = bh.next
            if bh is None:
                ii = np.where(self.rawdat['ID2']==self.rawdat['ID1'][i])[0]
                if len(ii)>0:
                    id3 = self.rawdat['ID1'][ii[0]]
                    bh = db.get_halo(str(dbsim.timesteps[ind].path)+'/1.'+str(id3))
            if bh is None:
                continue
            if 'host_halo' in bh.keys():
                hosth = bh['host_halo']
            else:
                hosth = None
            if hosth is not None:
                if "BH_central" in hosth.keys():
                    if type(hosth['BH_central'])==list:
                        cenbhs = np.array([halobh.halo_number for halobh in hosth['BH_central']])
                    else:
                        cenbhs = np.array([hosth['BH_central']])
                    if bh.halo_number not in cenbhs:
                        self.rawdat['sat_flag'][i] = 1
                else:
                    self.rawdat['sat_flag'][i] = 1
            else:
                self.rawdat['sat_flag'][i] = -1
            for p in halo_props:
                try:
                    self.rawdat[p][i] = bh.calculate('host_halo.'+p)
                except:
                    continue

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

        self.gwemit = gw.GWemit(self.rawdat['merge_mass_1'][ok], self.rawdat['merge_mass_2'][ok],self.rawdat['redshift'][ok],
                                sim.properties['omegaM0'], sim.properties['omegaL0'], sim.properties['h'])

        self.rawdat['GW_freq_merge'] = np.ones(len(self.rawdat['ID1']))*-1
        self.rawdat['GW_strain_merge'] = np.ones(len(self.rawdat['ID1']))*-1
        self.rawdat['GW_freq_ring'] = np.ones(len(self.rawdat['ID1']))*-1
        self.rawdat['GW_strain_ring'] = np.ones(len(self.rawdat['ID1']))*-1

        self.rawdat['GW_freq_merge'][ok] = self.gwemit.freq_merge()
        self.rawdat['GW_freq_ring'][ok] = self.gwemit.freq_ring()
        self.rawdat['GW_strain_merge'][ok] = self.gwemit.ampGW_merger() * 2.0 * self.rawdat['GW_freq_merge'][ok]
        self.rawdat['GW_strain_ring'][ok] = self.gwemit.ampGW_ring() * 2.0 * self.rawdat['GW_freq_ring'][ok]

    def get_vol_weight_hmf(self,hmf):
        self.rawdat['volweight'] = np.ones(len(self.rawdat['ID1'])) * -1
        if 'Mvir' not in self.rawdat.keys():
            print "halo information not found!"
            raise RuntimeError
        for i in range(len(self.rawdat['ID1'])):
            if self.rawdat['Mvir'][i] < 0:
                continue
            self.rawdat['volweight'][i] = hmf.calc_rho(np.log10(self.rawdat['Mvir'][i]),self.rawdat['snap_after'][i])

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

    def get_dual_frac(self,bhorbit,minL=1e43,maxD=10,boxsize=25,comove=True, gather_array=False, timestep=None):
        if gather_array is True:
            self.closeBHevol = {'dist':[], 'lum1':[], 'lum2':[],'ID2':[], 'ID1':[], 'z':[], 'time':[],'mass1':[],'mass2':[]}
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

            mass1 = bhorbit.single_BH_data(self.rawdat['ID1'][i],'mass')
            mass2 = bhorbit.single_BH_data(self.rawdat['ID2'][i],'mass')
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
                              np.abs(bphys[badx] - np.abs(xd[badx]))

            bady = np.where(yd > bphys/2)[0]
            yd[bady] = -1.0 * (yd[bady]/np.abs(yd[bady])) * \
                              np.abs(bphys[bady] - np.abs(yd[bady]))
            badz = np.where(zd > bphys/2)[0]
            zd[badz] = -1.0 * (zd[badz]/np.abs(zd[badz])) * \
                              np.abs(bphys[badz] - np.abs(zd[badz]))

            dist = np.sqrt(xd**2 + yd**2 + zd**2)
            if len(use1) != len(dist) or len(use2) != len(dist):
                print len(use1), len(use2), len(dist)
            if comove:
                dist /= scale1[use1]
            close = np.where(dist<maxD)[0]
            if len(close) > 0:
                if timestep:
                    dt = len(close)*timestep
                else:
                    try:
                        dt = np.sum(np.abs(time1[use1[close]] - time1[use1[close]-1]))
                    except:
                        dt = np.sum(np.abs(time1[use1[close]] - time1[use1[close]+1]))
                self.rawdat[tstr][i] = dt
                dual = np.where((lum1[use1[close]]>minL)&(lum2[use2[close]]>minL))[0]
                self.rawdat[fstr][i] = float(len(dual))/float(len(close))
                if gather_array is True:
                    self.closeBHevol['dist'].extend(dist[close])
                    self.closeBHevol['lum1'].extend(lum1[use1[close]])
                    self.closeBHevol['lum2'].extend(lum2[use2[close]])
                    self.closeBHevol['ID1'].extend(np.ones(len(close))*self.rawdat['ID1'][i])
                    self.closeBHevol['ID2'].extend(np.ones(len(close))*self.rawdat['ID2'][i])
                    self.closeBHevol['z'].extend(scale1[use1[close]]**-1 - 1)
                    self.closeBHevol['time'].extend(time1[use1[close]])
                    self.closeBHevol['mass1'].extend(mass1[use1[close]])
                    self.closeBHevol['mass2'].extend(mass2[use2[close]])
        if gather_array is True:
            for key in self.closeBHevol.keys():
                self.closeBHevol[key] = np.array(self.closeBHevol[key])

    def get_halo_interaction(self, dbsim):
        import tangos as db
        self.rawdat['hinteract_mvir_1'] = np.ones(len(self.rawdat['ID1']))*-1
        self.rawdat['hinteract_mvir_2'] = np.ones(len(self.rawdat['ID1']))*-1
        self.rawdat['hinteract_mstar_1'] = np.ones(len(self.rawdat['ID1']))*-1
        self.rawdat['hinteract_mstar_2'] = np.ones(len(self.rawdat['ID1']))*-1
        self.rawdat['hinteract_mgas_1'] = np.ones(len(self.rawdat['ID1']))*-1
        self.rawdat['hinteract_mgas_2'] = np.ones(len(self.rawdat['ID1']))*-1
        self.rawdat['hinteract_mbh_1'] = np.ones(len(self.rawdat['ID1']))*-1
        self.rawdat['hinteract_mbh_2'] = np.ones(len(self.rawdat['ID1']))*-1
        self.rawdat['dt_hinteract'] = np.ones(len(self.rawdat['ID1']))*-1

        nodiff = 0
        badmatch = 0
        for i in range(len(self.rawdat['ID1'])):
            if i%30 == 0:
                print float(i)/float(len(self.rawdat['ID1']))*100, '% done'
            try:
                bh1 = db.get_halo(str(dbsim.path)+'/%'+str(self.rawdat['snap_before'][i])+'/1.'+str(self.rawdat['ID1'][i]))
                bh2 = db.get_halo(str(dbsim.path)+'/%'+str(self.rawdat['snap_before'][i])+'/1.'+str(self.rawdat['ID2'][i]))
            except:
                continue

            try:
                time1, hn1, mv1, mg1, ms1, mbh1, dbh1, ssc, rvir = bh1.reverse_property_cascade('t()', 'host_halo.halo_number()',
                                                          'host_halo.Mvir', 'host_halo.Mgas','host_halo.Mstar', 'BH_mass','BH_central_distance',
                                                            'host_halo.SSC', 'host_halo.Rvir')
                time2, hn2, mv2, mg2, ms2, mbh2, dbh2, ssc2, rvir2 = bh2.reverse_property_cascade('t()', 'host_halo.halo_number()',
                                                          'host_halo.Mvir', 'host_halo.Mgas','host_halo.Mstar', 'BH_mass','BH_central_distance',
                                                            'host_halo.SSC', 'host_halo.Rvir')
            except:
                continue

            match1 = np.where(np.in1d(time1,time2))[0]
            match2 = np.where(np.in1d(time2,time1))[0]
            if not np.array_equal(time1[match1],time2[match2]):
                print "WARNING time arrays don't match!"
            if len(match1)==0 or len(match2)==0:
                continue
            ssc = ssc[match1]
            ssc2 = ssc2[match2]
            rvir = rvir[match1]
            rvir2 = ssc2[match2]

            dist = np.sqrt(np.sum((ssc-ssc2)**2,axis=1))
            good = np.where(dist > np.maximum(rvir,rvir2))[0]
            if len(good)==0: continue
                #good = np.where(dist >0)[0]
                #if len(good) == 0:
                #    continue
            self.rawdat['hinteract_mvir_1'][i] = mv1[match1[good[0]]]
            self.rawdat['hinteract_mvir_2'][i] = mv2[match2[good[0]]]
            self.rawdat['hinteract_mstar_1'][i] = ms1[match1[good[0]]]
            self.rawdat['hinteract_mstar_2'][i] = ms2[match2[good[0]]]
            self.rawdat['hinteract_mgas_1'][i] = mv1[match1[good[0]]]
            self.rawdat['hinteract_mgas_2'][i] = mv2[match2[good[0]]]
            self.rawdat['hinteract_mbh_1'][i] = mbh1[match1[good[0]]]
            self.rawdat['hinteract_mbh_2'][i] = mbh2[match2[good[0]]]

            self.rawdat['dt_hinteract'][i] = self.rawdat['time'][i] = time1[match2[good[0]]]



    def get_halo_merger(self,dbsim,overwrite=False, detail=False):
        import tangos as db
        if 'dt_hmerger' not in self.rawdat.keys() or overwrite==True:
            self.rawdat['dt_hmerger'] = np.ones(len(self.rawdat['ID1']))*-1
            self.rawdat['dt_hmerger_min'] = np.ones(len(self.rawdat['ID1']))*-1
            if detail == True:
                self.rawdat['hmerger_mvir_1'] = np.ones(len(self.rawdat['ID1']))*-1
                self.rawdat['hmerger_mvir_2'] = np.ones(len(self.rawdat['ID1']))*-1
                self.rawdat['hmerger_mgas_1'] = np.ones(len(self.rawdat['ID1']))*-1
                self.rawdat['hmerger_mgas_2'] = np.ones(len(self.rawdat['ID1']))*-1
                self.rawdat['hmerger_mstar_1'] = np.ones(len(self.rawdat['ID1']))*-1
                self.rawdat['hmerger_mstar_2'] = np.ones(len(self.rawdat['ID1']))*-1
                self.rawdat['hmerger_mbh_1'] = np.ones(len(self.rawdat['ID1']))*-1
                self.rawdat['hmerger_mbh_2'] = np.ones(len(self.rawdat['ID1']))*-1
                self.rawdat['hmerger_dbh_1'] = np.ones(len(self.rawdat['ID1']))*-1
                self.rawdat['hmerger_dbh_2'] = np.ones(len(self.rawdat['ID1']))*-1
            else:
                self.rawdat['hmerger_ndm_1'] = np.ones(len(self.rawdat['ID1']))*-1
                self.rawdat['hmerger_ndm_2'] = np.ones(len(self.rawdat['ID1']))*-1
        nodiff = 0
        badmatch = 0
        for i in range(len(self.rawdat['ID1'])):
            if i%30 == 0:
                print float(i)/float(len(self.rawdat['ID1']))*100, '% done'
            if self.rawdat['dt_hmerger'][i] >= 0 and overwrite == False:
                continue
            try:
                bh1 = db.get_halo(str(dbsim.path)+'/%'+str(self.rawdat['snap_before'][i])+'/1.'+str(self.rawdat['ID1'][i]))
                bh2 = db.get_halo(str(dbsim.path)+'/%'+str(self.rawdat['snap_before'][i])+'/1.'+str(self.rawdat['ID2'][i]))
            except:
                continue

            if bh1 is None or bh2 is None:
                self.rawdat['dt_hmerger'][i] = self.rawdat['time'][i] - min(self.rawdat['tform1'][i],self.rawdat['tform2'][i])
                self.rawdat['dt_hmerger_min'][i] = 0
                continue

            try:
                if detail==True:
                    time1, hn1, mv1, mg1, ms1, mbh1, dbh1 = bh1.reverse_property_cascade('t()', 'host_halo.halo_number()',
                                                          'host_halo.Mvir', 'host_halo.Mgas','host_halo.Mstar', 'BH_mass','BH_central_distance')
                    time2, hn2, mv2, mg2, ms2, mbh2, dbh2 = bh2.reverse_property_cascade('t()', 'host_halo.halo_number()',
                                                          'host_halo.Mvir', 'host_halo.Mgas','host_halo.Mstar', 'BH_mass','BH_central_distance')
                else:
                    time1, hn1, ndm1 = bh1.reverse_property_cascade('t()', 'host_halo.halo_number()', 'host_halo.NDM()')
                    time2, hn2, ndm2 = bh2.reverse_property_cascade('t()', 'host_halo.halo_number()', 'host_halo.NDM()')
            except:
                badmatch +=1
                continue
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
                self.rawdat['dt_hmerger'][i] = self.rawdat['time'][i] - max(self.rawdat['tform1'][i],self.rawdat['tform2'][i])
                self.rawdat['dt_hmerger_min'][i] = self.rawdat['time'][i] - max(time1.min(),time2.min())
                if detail == True:
                    self.rawdat['hmerger_mvir_1'][i] = mv1[-1]
                    self.rawdat['hmerger_mgas_1'][i] = mg1[-1]
                    self.rawdat['hmerger_mstar_1'][i] = ms1[-1]
                    self.rawdat['hmerger_mbh_1'][i] = mbh1[-1]
                    self.rawdat['hmerger_mbh_2'][i] = mbh2[-1]
                    self.rawdat['hmerger_dbh_1'][i] = dbh1[-1]
                    self.rawdat['hmerger_dbh_2'][i] = dbh2[-1]
                else:
                    self.rawdat['hmerger_ndm_1'][i] = ndm1[-1]
                continue

            th1 = time1[match1[diff[0]]]
            th2 = time2[match2[diff[0]]]
            if diff[0] != 0:
                th1p = time1[match1[diff[0]-1]]
            else:
                th1p = self.rawdat['time'][i]
            if th1 != th2:
                print "WARNING halo merge times not correct"

            self.rawdat['dt_hmerger'][i] = self.rawdat['time'][i] - th1
            if self.rawdat['time'][i] - th1p > 0:
                self.rawdat['dt_hmerger_min'][i] = self.rawdat['time'][i] - th1p
            else:
                self.rawdat['dt_hmerger_min'][i] = 0

            if detail == True:
                self.rawdat['hmerger_mvir_1'][i] = mv1[match1[diff[0]]]
                self.rawdat['hmerger_mvir_2'][i] = mv2[match2[diff[0]]]
                self.rawdat['hmerger_mgas_1'][i] = mg1[match1[diff[0]]]
                self.rawdat['hmerger_mgas_2'][i] = mg2[match2[diff[0]]]
                self.rawdat['hmerger_mstar_1'][i] = ms1[match1[diff[0]]]
                self.rawdat['hmerger_mstar_2'][i] = ms2[match2[diff[0]]]
                self.rawdat['hmerger_mbh_1'][i] = mbh1[match1[diff[0]]]
                self.rawdat['hmerger_mbh_2'][i] = mbh2[match2[diff[0]]]
                self.rawdat['hmerger_dbh_1'][i] = dbh1[match1[diff[0]]-1]
                self.rawdat['hmerger_dbh_2'][i] = dbh2[match2[diff[0]]-1]
            else:
                self.rawdat['hmerger_ndm_1'][i] = ndm1[match1[diff[0]]]
                self.rawdat['hmerger_ndm_2'][i] = ndm2[match2[diff[0]]]
        print "finished with ", nodiff, "BHs having never been in different halos and ", badmatch, "bad matches"

    def get_close_time(self,bhorbit, maxD, boxsize=25,comove=False, timestep=None):
        tstr = 'dt_hmerger_'+str(maxD)
        if comove:
            tstr = tstr+'c'
        for i in range(len(self.rawdat['ID1'])):
            if i %100==0:
                print float(i)/float(len(self.rawdat['ID1']))*100, '% done'

            x1 = bhorbit.single_BH_data(self.rawdat['ID1'][i],'x')
            x2 = bhorbit.single_BH_data(self.rawdat['ID2'][i],'x')
            y1 = bhorbit.single_BH_data(self.rawdat['ID1'][i],'y')
            y2 = bhorbit.single_BH_data(self.rawdat['ID2'][i],'y')
            z1 = bhorbit.single_BH_data(self.rawdat['ID1'][i],'z')
            z2 = bhorbit.single_BH_data(self.rawdat['ID2'][i],'z')

            time1 = bhorbit.single_BH_data(self.rawdat['ID1'][i],'time')
            time2 = bhorbit.single_BH_data(self.rawdat['ID2'][i],'time')

            scale1 = bhorbit.single_BH_data(self.rawdat['ID1'][i],'scalefac')
            scale2 = bhorbit.single_BH_data(self.rawdat['ID2'][i],'scalefac')

            if len(time1) == 0 or len(time2) == 0:
                continue

            tsame = self.rawdat['time'][i]-self.rawdat['dt_hmerger_min'][i]
            mint = np.max([time1[(time1>=tsame)].min(),time2[(time2>=tsame)].min()])
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
                              np.abs(bphys[badx] - np.abs(xd[badx]))

            bady = np.where(yd > bphys/2)[0]
            yd[bady] = -1.0 * (yd[bady]/np.abs(yd[bady])) * \
                              np.abs(bphys[bady] - np.abs(yd[bady]))
            badz = np.where(zd > bphys/2)[0]
            zd[badz] = -1.0 * (zd[badz]/np.abs(zd[badz])) * \
                              np.abs(bphys[badz] - np.abs(zd[badz]))

            dist = np.sqrt(xd**2 + yd**2 + zd**2)
            if len(use1) != len(dist) or len(use2) != len(dist):
                print len(use1), len(use2), len(dist)
            if comove:
                dist /= scale1[use1]

            close = np.where(dist < maxD)
            if len(close) > 0:
                if timestep:
                    dt = len(close)*timestep
                else:
                    try:
                        dt = np.sum(np.abs(time1[use1[close]] - time1[use1[close]-1]))
                    except:
                        dt = np.sum(np.abs(time1[use1[close]] - time1[use1[close]+1]))

                self.rawdat[tstr] = dt

    def get_new_masses(self,orig_seed = 1e6, new_seed = 1e3, useonly=None, doonly=None):
        from . import util
        self.rawdat['newmass1'] = util.get_new_masses(self.rawdat['ID1'],self.rawdat['time']-0.001,self.rawdat['merge_mass_1'],
                                                      self, orig_seed=orig_seed, new_seed = new_seed, useonly=useonly)
        self.rawdat['newmass2'] = util.get_new_masses(self.rawdat['ID2'],self.rawdat['time']-0.001,self.rawdat['merge_mass_2'],
                                                      self, orig_seed=orig_seed, new_seed = new_seed, useonly=useonly)













