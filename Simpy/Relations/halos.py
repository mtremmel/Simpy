import numpy as np
from ..readcol import readcol
from .. import util
from .. import plotting
import pynbody
import os
import gc


#Moster relations
def moster13allvar(logM,z,M10,M11,N10,N11,b10,b11,g10,g11):
    r = z/(z+1)
    logM1 = M10 + M11 * r
    N = N10 + N11 * r
    b = b10 + b11 * r
    g = g10 + g11 * r
    x = logM - logM1
    ratio = 2 * N / ( 10**(-b * x) + 10**(g * x) )
    return ratio

def moster13(logM, z):
    M10 = 11.590
    M11 =  1.195
    N10 =  0.0351
    N11 = -0.0247
    b10 =  1.376
    b11 = -0.826
    g10 =  0.608
    g11 =  0.329
    return moster13allvar(logM,z,M10,M11,N10,N11,b10,b11,g10,g11)

def errmoster13(logM,z):
    M10 = 11.590
    M11 =  1.195
    N10 =  0.0351
    N11 = -0.0247
    b10 =  1.376
    b11 = -0.826
    g10 =  0.608
    g11 =  0.329
    sigM10 = 0.236
    sigM11 = 0.353
    sigN10 = 0.0058
    sigN11 = 0.0069
    sigb10 = 0.153
    sigb11 = 0.225
    sigg10 = 0.059
    sigg11 = 0.173
    sigvar = [sigM10,sigM11,sigN10,sigN11,sigb10,sigb11,sigg10,sigg11]
    sigma = np.zeros(len(logM))
    for i in range(len(logM)):
        point = [logM[i],z,M10,M11,N10,N11,b10,b11,g10,g11]
        for j in range(8):
            sigma[i] += util.partial_derivative(moster13allvar,var=j+2,point=point)**2 * sigvar[j]**2
    sigma = np.sqrt(sigma)
    return sigma

#Behroozi data
def behroozi13(logM,a):
    z = a**-1 - 1
    v = np.exp(-4.*a**2)
    le = -1.777 + (-0.006*(a-1)*v) - 0.119*(a-1)
    lM1 = 11.514 + (-1.793*(a-1)-0.251*z)*v
    A = -1.412+0.731*(a-1)*v
    d = 3.508+(2.608*(a-1)-0.043*z)*v
    g = 0.316+(1.319*(a-1)+0.279*z)*v
    def f(x):
        return -1.0*np.log10(10**(A*x)+1) + d*((np.log10(1+np.exp(x)))**g)/(1+np.exp(10**(-x)))
    lMstar = le +lM1 + f(logM-lM1) - f(0)
    ratio = 10**(lMstar - logM)
    return ratio

#Kravstov data
def kravstov14(logM, z = 0.1):
    loge = -1.685
    logM1 = 11.39
    alpha = -1.74
    delta = 4.335
    gamma = 0.531

    def f(x):
        return -1.0 * np.log10(10.0**(alpha * x) + 1.0) + (delta * np.log10(1.0 + np.exp(x))**gamma) / (1.0 + np.exp(10.**(-1.0 * x)))

    logMstar = loge + logM1 + f(logM - logM1) - f(0)
    return 10**(logMstar - logM)


def getstats(simname, step):
    print "getting halo stats..."
    if os.path.exists(simname + '.' + step + '.amiga.stat'):
        amigastat = readcol(simname + '.' + step + '.amiga.stat', asdict=True)
        type = 'amiga'
    else:
            print "amiga file failed, looking for rockstar"
            if os.path.exists(simname + '.' + step + '.rockstar.stat'):
                amigastat = readcol(simname + '.' + step + '.rockstar.stat', asdict=True)
                type = 'rockstar'
            else:
                print "ERROR cannot find recognized stat file (amiga.stat or rockstar.stat)"
    s = pynbody.load(simname+'.'+step)
    a = s.properties['a']
    del(s)
    gc.collect()
    return amigastat, a**-1 -1, type


def SMHM_db(sim, step, style, fitstyle=['k-','k--'], fit=['Mos', 'Krav'], minmass=None,
            maxmass=None, markersize=5,label=None, correct=True, error=True):
    import tangos as db
    step = db.get_timestep(sim+'/%'+str(step))
    Mvir, Mstar = step.gather_property('Mvir', 'Mstar')

    if correct is True:
        Mstar *= 0.6
        Mvir /= 0.8

    plotting.plt.plot(np.log10(Mvir), np.log10(Mstar/Mvir), style, markersize=markersize, label=label)

    if minmass is None:
        minmass = Mvir.min()/2.
    if maxmass is None:
        maxmass = Mvir.max()*2.

    cnt = 0
    logmv_fit = np.arange(np.log10(minmass),np.log10(maxmass),0.1)
    if fit is not None:
        for ff in fit:
            if ff not in ['Mos','Beh', 'Krav', 'Moster', 'Behroozi', 'Kravtsov']:
                print "fit request ", ff, " not understood... should be in list", \
                    ['Mos','Beh', 'Krav', 'Moster', 'Behroozi', 'Kravtsov']
                continue
            print "found fit relation for", ff

            if ff in ['Mos', 'Moster']:
                fitfunc = moster13
                flabel = 'Moster+ 13'
            if ff in ['Beh', 'Behroozi']:
                fitfunc = behroozi13
                flabel = 'Behroozi+ 13'
            if ff in ['Krav', 'Kravtsov']:
                fitfunc = kravstov14
                flabel = 'Kravtsov+ 14'

            ratio_fit = fitfunc(logmv_fit, step.redshift)
            plotting.plt.plot(logmv_fit, np.log10(ratio_fit), fitstyle[cnt], label=flabel)
            if ff in ['Mos','Moster'] and error is True:
                sigma = errmoster13(logmv_fit, step.redshift)
                plotting.plt.fill_between(logmv_fit,np.log10(ratio_fit-sigma),np.log10(ratio_fit+sigma),
                                        facecolor='grey',alpha=0.5)

            cnt += 1

    plotting.plt.ylabel(r'M$_{*,cen}$/M$_{vir}$',fontsize=40)
    plotting.plt.xlabel(r'M$_{vir}$ [M$_{\odot}$]',fontsize=40)
    plotting.plt.legend(loc='lower right',fontsize=30)


def SMHM(sim, step, style, color, fitstyle=['k-','k--'], fit=['Mos', 'Krav'], minmass=None, maxmass=None,
         markersize=5,label=None, correct=True, usedb=False, remove_sats=True, only_sats=False, error=True, alpha=1.0,
         remove_sats_hard=False, boxsize=25, ret_data=False):
    if usedb is True:
        import tangos as db
        dbstep = db.get_timestep(sim+'/%'+step)
        if remove_sats or only_sats:
            Mvir, Mstar, sub = dbstep.gather_property('Mvir', 'Mstar', 'Sub')
            if len(Mvir)==0:
                print "Warning Sub halos not implemented yet in your database! " \
                      "please turn off sats options or re-run with usedb False"
                return
        else:
            Mvir, Mstar = dbstep.gather_property('Mvir', 'Mstar')
        redshift = dbstep.redshift
    else:
        amigastat, redshift, type = getstats(sim, step)
        Mvir = amigastat['Mvir(M_sol)'].astype(np.float)
        Mstar = amigastat['StarMass(M_sol)'].astype(np.float)
        Mgas = amigastat['GasMass(M_sol)'].astype(np.float)
        xC = amigastat['Xc'].astype(np.float)
        yC = amigastat['Yc'].astype(np.float)
        zC = amigastat['Zc'].astype(np.float)
        Rvir = amigastat['Rvir(kpc)'].astype(np.float)
        grp = amigastat['Grp'].astype(np.int64)

        if minmass is None:
            minmass = Mvir.min()/2.
        if maxmass is None:
            maxmass = Mvir.max()*2.

        if remove_sats is True:
            if type == 'rockstar':
                ok, = np.where(amigastat['Satellite?'] == -1)
            if type == 'amiga':
                ok, = np.where(amigastat['Satellite?'] == 'no')
        if remove_sats_hard is True:
            satsarray = np.zeros(len(amigastat['Grp']))
            for i in range(len(amigastat['Grp'])):
                if Mvir[i] < minmass or Mvir[i] > maxmass:
                    continue
                xd = xC[i] - xC
                yd = yC[i] - yC
                zd = zC[i] - zC
                oo = np.where(np.abs(xd) > boxsize*1e3/2.)[0]
                xd[oo] = -1.0 * (xd[oo]/np.abs(xd[oo])) * (boxsize*1e3 - np.abs(xd[oo]))
                oo = np.where(np.abs(yd) > boxsize*1e3/2.)[0]
                yd[oo] = -1.0 * (yd[oo]/np.abs(yd[oo])) * (boxsize*1e3 - np.abs(yd[oo]))
                oo = np.where(np.abs(zd) > boxsize*1e3/2.)[0]
                zd[oo] = -1.0 * (zd[oo]/np.abs(zd[oo])) * (boxsize*1e3 - np.abs(zd[oo]))

                dist = np.sqrt(xd**2 + yd**2 + zd**2)
                bad = np.where((dist*1e3<2*(Rvir+Rvir[i]))&(0.5*Mvir[i]<Mvir)&(grp[i] != grp))[0]
                if len(bad)>0 or Mgas[i]==0:
                    satsarray[i] = 1
            if type == 'rockstar':
                ok = np.where((satsarray==0)&(amigastat['Satellite?']==-1))[0]
            if type == 'amiga':
                ok = np.where((satsarray==0)&(amigastat['Satellite?']=='no'))[0]
        if only_sats is True:
            if type == 'rockstar':
                ok, = np.where(amigastat['Satellite?'] != -1)
            if type == 'amiga':
                ok, = np.where(amigastat['Satellite?'] == 'yes')
        if only_sats is True or remove_sats is True or remove_sats_hard:
            Mvir = Mvir[ok]
            Mstar = Mstar[ok]
            grp = grp[ok]

    if correct is True:
        Mstar *= 0.6
        Mvir /= 0.8

    if minmass is None:
        minmass = Mvir.min()/2.
    if maxmass is None:
        maxmass = Mvir.max()*2.

    logmv_fit = np.arange(np.log10(minmass),np.log10(maxmass),0.1)
    cnt = 0
    plotting.plt.plot(Mvir, Mstar/Mvir, style, color=color, markersize=markersize, label=label, alpha=alpha,zorder=1)
    if fit is not None:
        for ff in fit:
            if ff not in ['Mos','Beh', 'Krav', 'Moster', 'Behroozi', 'Kravtsov']:
                print "fit request ", ff, " not understood... should be in list", \
                    ['Mos','Beh', 'Krav', 'Moster', 'Behroozi', 'Kravtsov']
                continue
            print "found fit relation for", ff

            if ff in ['Mos', 'Moster']:
                fitfunc = moster13
                flabel = 'Moster+ 13, z = %0.2f' % (redshift)
            if ff in ['Beh', 'Behroozi']:
                fitfunc = behroozi13
                flabel = 'Behroozi+ 13, z = %0.2f' % (redshift)
            if ff in ['Krav', 'Kravtsov']:
                fitfunc = kravstov14
                flabel = 'Kravtsov+ 14, z < 0.1'

            ratio_fit = fitfunc(logmv_fit, redshift)
            plotting.plt.plot(10**logmv_fit, ratio_fit, fitstyle[cnt], label=flabel, lw=5,alpha=0.75,zorder=10)
            if ff in ['Mos', 'Moster'] and error is True:
                sigma = errmoster13(logmv_fit, redshift)
                plotting.plt.fill_between(10**logmv_fit,ratio_fit-sigma,ratio_fit+sigma, facecolor='grey',
                                          edgecolor='k', lw=1.5, alpha=0.5,zorder=10)

            cnt += 1

    plotting.plt.ylabel(r'M$_{*}$/M$_{vir}$')
    plotting.plt.xlabel(r'M$_{vir}$ [M$_{\odot}$]')
    plotting.plt.legend(loc='lower right',fontsize=25)
    plotting.plt.yscale('log',base=10)
    plotting.plt.xscale('log',base=10)
    plotting.plt.xlim(minmass, maxmass)

    if ret_data:
        return Mvir, Mstar, grp

def mergerhist(dbsim, moreprops=None, names=None, hmax=3000,nstart=0,nend=-1):
    data = {'time': [], 'Mstar1':[], 'Mstar2':[], 'Mstarf':[],'Mvir1':[], 'Mvir2':[], 'Mvirf':[],
            'Mgas1':[], 'Mgas2':[], 'Mgasf':[], 'nMerge':[], 'dbstep':[],'redshift':[],
            'N1':[], 'N2':[], 'Nf':[]}
    for step in dbsim.timesteps[nstart:nend]:
        print step
        N = []
        Nf = []
        Mvir = []
        Mstar = []
        Mgas = []
        Mvirf = []
        Mgasf = []
        Mstarf = []
        cnt = 0
        for h in step.halos:
            if h.halo_number > hmax: break
            if h.next is None: continue
            #if 'Mvir' not in h.keys() or 'Mvir' not in h.next.keys(): continue
            if 'BH_central' not in h.keys() or 'BH_central' not in h.next.keys(): continue
            N.append(h.halo_number)
            Nf.append(h.next.halo_number)
            Mvir.append(h.NDM*3.4e5)
            Mvirf.append(h.next.NDM*3.4e5)
            Mstar.append(h.NStar)
            Mstarf.append(h.next.NStar)
            Mgas.append(h.NGas)
            Mgasf.append(h.next.NGas)
            if cnt%200 == 0:
                print cnt/2000.
            cnt += 1
        N = np.array(N)
        Nf = np.array(Nf)
        Mvir = np.array(Mvir)
        Mvirf = np.array(Mvirf)
        Mstar = np.array(Mstar)
        Mstarf = np.array(Mstarf)
        Mgas = np.array(Mgas)
        Mgasf = np.array(Mgasf)


        ziparr = np.array(zip(Nf,1./Mvir),dtype=[('nf','int64'),('mv','float')])
        order = np.argsort(ziparr,order=('nf','mv'))
        Nf = Nf[order]
        N = N[order]
        Mvir = Mvir[order]
        Mvirf = Mvirf[order]
        Mstar = Mstar[order]
        Mstarf = Mstarf[order]
        Mgas = Mgas[order]
        Mgasf = Mgasf[order]

        uNf, ind, inv, cnt = np.unique(Nf,return_index=True, return_inverse=True,return_counts=True)

        mm = np.where(cnt > 1)[0]

        if len(mm)==0:
            print "No Mergers This Step"
        indm = ind[mm]

        data['Mvirf'].extend(Mvirf[indm])
        data['Mstarf'].extend(Mstarf[indm])
        data['Mgasf'].extend(Mgasf[indm])
        data['Nf'].extend(Nf[indm])

        data['Mvir1'].extend(Mvir[indm])
        data['Mvir2'].extend(Mvir[indm+1])
        data['Mstar1'].extend(Mstar[indm])
        data['Mstar2'].extend(Mstar[indm+1])
        data['Mgas1'].extend(Mgas[indm])
        data['Mgas2'].extend(Mgas[indm+1])
        data['N1'].extend(N[indm])
        data['N2'].extend(N[indm+1])
        data['redshift'].extend(np.ones(len(indm))*step.redshift)
        fextent = np.zeros(len(indm)).astype(dtype='S64')
        fextent[:] = step.relative_filename
        data['dbstep'].extend(fextent)

        data['nMerge'].extend(cnt[mm])
    for key in data.keys():
        data[key] = np.array(data[key])

    return data


def mergerhist_slow(dbsim, moreprops=None, names=None, ret_totals=False, nsteps = None):
    data = {'time': [], 'Mstar1':[], 'Mstar2':[], 'Mstarf':[],'Mvir1':[], 'Mvir2':[], 'Mvirf':[],
            'Mgas1':[], 'Mgas2':[], 'Mgasf':[], 'nMerge':[], 'dbstep':[],'redshift':[],
            'N1':[], 'N2':[], 'Nf':[]}
    allprops = ['t()', 'halo_number()','later(1).halo_number()', 'Mvir', 'later(1).Mvir',
                'Mstar', 'later(1).Mstar', 'Mgas', 'later(1).Mgas']

    if names is None or len(names)!=len(moreprops):
        names = moreprops

    if names is not None:
        for n in names:
            data[n+'1'] = []
            data[n+'2'] = []
            data[n+'f'] = []
    else:
        names = []

    for p in moreprops:
        allprops.append(p)
        allprops.append('later(1).'+p)

    print allprops

    nmergers = []
    redshifts = []
    times = []

    loopcnt = 0

    for step in dbsim.timesteps:
        if nsteps is not None:
            if loopcnt >= nsteps: break
        loopcnt +=1
        print step
        try:
            dbdata = step.gather_property(*allprops)

        except:
            print "No halo data this step"
            nmergers.append(0)
            times.append(step.time_gyr)
            redshifts.append(step.redshift)
            continue
        time = dbdata[0]
        N = dbdata[1]
        Nf = dbdata[2]
        Mvir = dbdata[3]
        Mvirf = dbdata[4]
        Mstar = dbdata[5]
        Mstarf = dbdata[6]
        Mgas = dbdata[7]
        Mgasf = dbdata[8]

        ziparr = np.array(zip(Nf,1./Mvir),dtype=[('nf','int64'),('mv','float')])
        order = np.argsort(ziparr,order=('nf','mv'))
        Nf = Nf[order]
        time = time[order]
        N = N[order]
        Mvir = Mvir[order]
        Mvirf = Mvirf[order]
        Mstar = Mstar[order]
        Mstarf = Mstarf[order]
        Mgas = Mgas[order]
        Mgasf = Mgasf[order]

        for i in range(len(dbdata)-9):
            dbdata[9+i] = dbdata[9+i][order]

        uNf, ind, inv, cnt = np.unique(Nf,return_index=True, return_inverse=True,return_counts=True)

        mm = np.where(cnt > 1)[0]

        if len(mm)==0:
            print "No Mergers This Step"
            nmergers.append(0)
            times.append(step.time_gyr)
            redshifts.append(step.redshift)
            continue
        indm = ind[mm]

        nmergers.append(np.sum(cnt[mm]-1))
        times.append(step.time_gyr)
        redshifts.append(step.redshift)

        data['time'].extend(time[indm])
        data['Mvirf'].extend(Mvirf[indm])
        data['Mstarf'].extend(Mstarf[indm])
        data['Mgasf'].extend(Mgasf[indm])
        data['Nf'].extend(Nf[indm])

        data['Mvir1'].extend(Mvir[indm])
        data['Mvir2'].extend(Mvir[indm+1])
        data['Mstar1'].extend(Mstar[indm])
        data['Mstar2'].extend(Mstar[indm+1])
        data['Mgas1'].extend(Mgas[indm])
        data['Mgas2'].extend(Mgas[indm+1])
        data['N1'].extend(N[indm])
        data['N2'].extend(N[indm+1])
        data['redshift'].extend(np.ones(len(indm))*step.redshift)
        fextent = np.zeros(len(indm)).astype(dtype='S64')
        fextent[:] = step.relative_filename
        data['dbstep'].extend(fextent)

        indtrack = 9

        for n in names:
            data[n+'1'].extend(dbdata[indtrack][indm])
            data[n+'2'].extend(dbdata[indtrack][indm+1])
            data[n+'f'].extend(dbdata[indtrack+1][indm])
            indtrack = indtrack+2

        data['nMerge'].extend(cnt[mm])
    for key in data.keys():
        data[key] = np.array(data[key])
    if ret_totals:
        return data, nmergers, times, redshifts
    else:
        return data

class HaloMergers(object):
    def __init__(self, dbsim, mhist=None, moreprops=None, names=None):
        if dbsim and not mhist:
            self.data = mergerhist(dbsim, moreprops, names)
        if mhist:
            self.data = mhist
        if not mhist and not dbsim:
            print "ERROR requires either a dbsim object or an already made merger data structure"

    def __getitem__(self, item):
        return self.data[item]

    def keys(self):
        return self.data.keys()

    def dtBHmerge(self):
        import tangos as db
        self.data['dtBHmerge'] = np.ones(len(self.data['time']))*-1
        self.data['dtBHmerge_min'] = np.ones(len(self.data['time']))*-1
        if 'time_next' not in self.data.keys():
            self.get_tnext(self)
        for ii in range(len(self.data['time'])):
            if ii % 100 == 0:
                print float(ii)/float(len(self.data['time'])) * 100, "% done"
            bh1 = db.get_halo(self.data['dbstep'][ii]+'/1.'+str(self.data['bh().halo_number()1'][ii]))
            bh2 = db.get_halo(self.data['dbstep'][ii]+'/1.'+str(self.data['bh().halo_number()2'][ii]))
            bhc1, time1 = bh1.property_cascade('halo_number()', 't()')
            bhc2, time2 = bh2.property_cascade('halo_number()', 't()')
            im = np.where(bhc1==bhc2)[0]
            if len(im)>0:
                tmerge = time1[im[0]]
                if tmerge != time2[im[0]]:
                    print "WEIRD ONE"
                self.data['dtBHmerge'][ii] = tmerge - self.data['time'][ii]
                self.data['dtBHmerge_min'][ii] = self.data['dtBHmerge'][ii] - self.data['time_next'][ii]


    def get_tnext(self):
        import tangos as db
        self.data['time_next'] = np.ones(len(self.data['time']))*-1
        self.data['redshift_next'] = np.ones(len(self.data['time']))*-1
        for ii in range(len(self.data['time'])):
            if ii % 100 == 0:
                print float(ii)/float(len(self.data['time'])) * 100, "% done"
            step = db.get_timestep(self.data['dbstep'][ii])
            self.data['time_next'][ii] = step.next.time_gyr
            self.data['redshift_next'][ii] = step.next.redshift

    def get_tstart(self):
        import tangos as db
        self.data['tstart'] = np.ones(len(self.data['time']))*-1
        self.data['Mvir1_start'] = np.ones(len(self.data['time']))*-1
        self.data['Mvir2_start'] = np.ones(len(self.data['time']))*-1
        self.data['Mstar1_start'] = np.ones(len(self.data['time']))*-1
        self.data['Mstar2_start'] = np.ones(len(self.data['time']))*-1
        self.data['Mgas1_start'] = np.ones(len(self.data['time']))*-1
        self.data['Mgas2_start'] = np.ones(len(self.data['time']))*-1
        for ii in range(len(self.data['time'])):
            if ii % 100 == 0:
                print float(ii)/float(len(self.data['time'])) * 100, "% done"
            h1 = db.get_halo(self.data['dbstep'][ii]+'/'+str(self.data['N1'][ii]))
            h2 = db.get_halo(self.data['dbstep'][ii]+'/'+str(self.data['N2'][ii]))

            pos1, Rvir1, time1, Mvir1, Mstar1, Mgas1 = h1.reverse_property_cascade('SSC', 'Rvir', 't()', 'Mvir', 'Mstar', 'Mgas')
            pos2, Rvir2, time2, Mvir2, Mstar2, Mgas2 = h2.reverse_property_cascade('SSC', 'Rvir', 't()', 'Mvir', 'Mstar', 'Mgas')
            ll = min(len(time1), len(time2))
            dist = np.sqrt(np.sum((pos1[0:ll] - pos2[0:ll])**2,axis=1))
            sep = np.where(dist>Rvir1[0:ll]+Rvir2[0:ll])[0]
            if len(sep)==0:
                #print "cannot find output when objects were not close"
                continue
            self.data['tstart'][ii] = time1[sep[0]]
            self.data['Mvir1_start'][ii] = Mvir1[sep[0]]
            self.data['Mvir2_start'][ii] = Mvir2[sep[0]]
            self.data['Mstar1_start'][ii] = Mstar1[sep[0]]
            self.data['Mstar2_start'][ii] = Mstar2[sep[0]]
            self.data['Mgas1_start'][ii] = Mgas1[sep[0]]
            self.data['Mgas2_start'][ii] = Mgas2[sep[0]]











