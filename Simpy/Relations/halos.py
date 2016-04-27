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
    import halo_db as db
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
         remove_sats_hard=False, boxsize=25):
    if usedb is True:
        import halo_db as db
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

        if remove_sats is True:
            if type == 'rockstar':
                ok, = np.where(amigastat['Satellite?'] == -1)
            if type == 'amiga':
                ok, = np.where(amigastat['Satellite?'] == 'no')
        if remove_sats_hard is True:
            satsarray = np.zeros(len(amigastat['Grp']))
            for i in range(len(amigastat['Grp'])):
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


def mergerhist(dbsim,volume=25**3,mfconv=False, ret_totals=False):
    data = {'time': [], 'Mstar1':[], 'Mstar2':[], 'Mstarf':[],'Mvir1':[], 'Mvir2':[], 'Mvirf':[],
            'Mgas1':[], 'Mgas2':[], 'Mgasf':[], 'nMerge':[],
            'N1':[], 'N2':[], 'Nf':[]}

    nmergers = []
    redshifts = []
    times = []

    for step in dbsim.timesteps:
        print step
        try:
            time, N, Nf, Mvir, Mvirf, Mstar, Mstarf, Mgas, Mgasf = \
                step.gather_property('t()', 'halo_number()','later(1).halo_number()', 'Mvir', 'later(1).Mvir',
                                    'Mstar', 'later(1).Mstar', 'Mgas', 'later(1).Mgas')
        except:
            print "No halo data this step"
            nmergers.append(0)
            times.append(step.time_gyr)
            redshifts.append(step.redshift)
            continue


        ordM = np.argsort(Mvir)
        Nf = Nf[ordM]
        N = N[ordM]
        Mvir = Mvir[ordM]
        Mvirf = Mvirf[ordM]
        Mstar = Mstar[ordM]
        Mstarf = Mstarf[ordM]
        Mgas = Mgas[ordM]
        Mgasf = Mgas[ordM]

        ordNf = np.argsort(Nf)
        Nf = Nf[ordNf]
        time = time[ordNf]
        N = N[ordNf]
        Mvir = Mvir[ordNf]
        Mvirf = Mvirf[ordNf]
        Mstar = Mstar[ordNf]
        Mstarf = Mstarf[ordNf]
        Mgas = Mgas[ordNf]
        Mgasf = Mgas[ordNf]

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

        data['nMerge'].extend(cnt[mm])
    for key in data.keys():
        data[key] = np.array(data[key])
    if ret_totals:
        return data, nmergers, times, redshifts
    else:
        return data









