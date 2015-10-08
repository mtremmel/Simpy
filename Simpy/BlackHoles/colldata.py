from .. import readcol
import numpy as np

#Accretion history data
#Shankar 09
shankar09L = 3.2e5
shankar09H = 5.4e5
shankar09 = (shankar09H +shankar09L) / 2.
#Salvaterra 12
Salvaterra12 = 0.66e4
Salvaterra12zH = 9
Salvaterra12zL = 5
Salvaterra12z = (Salvaterra12zH + Salvaterra12zL)/2.
#Treister 13
Treister13 = np.array([851.,666.,674.])
Treister13z = np.array([6.5,7.5,8.5])
Treister13zErr = np.array([.5,.5,.5])
#Hopkins 07
Hopkins07zp1,Hopkins07 = readcol("../DATA/QSOdata/RhoAccZ.csv",twod=False)
Hopkins07zp1H,Hopkins07H = readcol("../DATA/QSOdata/RhoAccZPLUS.csv",twod=False)
Hopkins07zp1L,Hopkins07L = readcol("../DATA/QSOdata/RhoAccZMINUS.csv",twod=False)
Hopkins07perr = 10**Hopkins07H - 10**Hopkins07
Hopkins07merr = 10**Hopkins07 - 10**Hopkins07L

#BH Luminosity Function Data
#Hopkins 07
hop_bhlf_z =      np.array([1.  ,1.5 ,2.  ,2.5 ,3.  ,3.5 ,4.  ,4.5 ,5.  ,5.5 ,6.])
hop_bhlf_zbinsl = np.array([0.75,1.25,1.75,2.25,2.75,3.25,3.75,4.25,4.75,5.25,5.75])
hop_bhlf_zbinsh = np.array([1.25,1.75,2.25,2.75,3.25,3.75,4.25,4.75,5.25,5.75,6.25])
hop_bhlf_obs = readcol("../DATA/QSOdata/bol_lf_point_dump.dat",twod=False,asdict=True,skipline=38)
#McGreer 13
mcg_bhlf_obs = readcol("../DATA/QSOdata/M1450z5_McGreer13.dat",twod=False,asdict=True,skipline=1)
#Barger 03, z = 6
bar_bhlfz6_phi = np.log10(1e-6)
bar_bhlfz6_phierr = 0.6
bar_bhlfz6_L = 45.05
bar_bhlfz6_Lerr = 0.5
#Fiore+ 2012 number for z>5.8 quasars
fio_bhlf_phi6F = np.log10(0.66e-5)
fio_bhlf_errphi6Fp = np.log10(0.66e-5+1.1e-5) - logphi6F
fio_bhlf_errphi6Fm = logphi6F - np.log10(0.66e-5-0.5e-5)
fio_bhlf_L6F = 45.55
fio_bhlf_L6Fm = 45.55 - 44.9
fio_bhlf_L6Fp = 46.2 - 45.55