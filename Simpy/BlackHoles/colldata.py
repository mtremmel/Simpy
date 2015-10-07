from .. import readcol
import numpy as np

shankar09L = 3.2e5
shankar09H = 5.4e5
Salvaterra12 = 0.66e4
Salvaterra12zH = 9
Salvaterra12zL = 5
Salvaterra12z = (Salvaterra12zH + Salvaterra12zL)/2.
Treister13 = np.array([851.,666.,674.])
Treister13z = np.array([6.5,7.5,8.5])
Treister13zErr = np.array([.5,.5,.5])
Hopkins07zp1,Hopkins07 = readcol.readcol("/nobackupp8/mtremmel/DATA/QSOdata/RhoAccZ.csv",twod=False)
Hopkins07zp1H,Hopkins07H = readcol.readcol("/nobackupp8/mtremmel/DATA/QSOdata/RhoAccZPLUS.csv",twod=False)
Hopkins07zp1L,Hopkins07L = readcol.readcol("/nobackupp8/mtremmel/DATA/QSOdata/RhoAccZMINUS.csv",twod=False)
Hopkins07perr = 10**Hopkins07H - 10**Hopkins07
Hopkins07merr = 10**Hopkins07 - 10**Hopkins07L