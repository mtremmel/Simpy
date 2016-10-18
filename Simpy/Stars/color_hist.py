import numpy as np
from Simpy import plotting

CANDELS_M31 = {'redshift':[0.45,  0.8, 1.0,1.25, 1.55, 1.85,  2.1,  2.5,3.15],
        'lMstar':[10.85,10.81,10.8,10.7,10.62,10.48,10.36,10.15, 9.8],
        'n': [4.2,3.6,3.0,2.5,2.2,1.8,1.0,1.1,1.3],
        'n+':[1.3,1.3,1.2,2.5,1.7,1.7,1.5,2.0,1.9],
        'n-':[1.5,1.0,1.3,1.3,1.3,1.1,0.5,0.6,0.7],
        'UV':[2.0,1.9,1.7,1.7,1.6,1.5,1.2,0.9,0.6],
        'UV+':[0.2,0.2,0.2,0.2,0.3,0.3,0.5,0.6,0.4],
        'UV-':[0.2,0.3,0.3,0.4,0.3,0.3,0.4,0.3,0.3],
        'VJ': [1.3,1.3,1.4,1.2,1.3,1.3,1.1,0.8,0.3],
        'VJ+':[0.1,0.2,0.2,0.3,0.4,0.4,0.5,0.5,0.9],
        'VJ-':[0.1,0.2,0.2,0.2,0.2,0.3,0.4,0.4,0.6]
        }
CANDELS_MW = { 'redshift':[0.45,  0.8,  1.0, 1.25, 1.55,1.85, 2.1, 2.5],
        'lMstar': [10.6,10.47,10.35,10.21,10.06,9.88, 9.7,9.48],
        'n':  [3.4,2.7,2.1,1.5,1.2,1.1,1.3,1.3],
        'n+': [1.7,1.4,1.5,2.1,1.5,1.1,1.4,1.4],
        'n-': [1.7,1.4,1.2,0.8,0.6,0.6,0.6,0.7],
        'UV': [1.9,1.7,1.6,1.4,1.1,0.9,0.6,0.6],
        'UV+':[0.2,0.2,0.3,0.4,0.5,0.4,0.4,0.3],
        'UV-':[0.3,0.3,0.3,0.4,0.3,0.3,0.3,0.3],
        'VJ': [1.3,1.3,1.2,1.2,1.0,0.8,0.5,0.3],
        'VJ+':[0.2,0.3,0.4,0.4,0.4,0.5,0.4,0.5],
        'VJ-':[0.1,0.2,0.2,0.3,0.3,0.4,0.3,0.4]
        }

def plt_UVJ_quench_region(lw=1.5, shade=True):
        plotting.plt.plot([1.5,1.5], [1.9,2.7], 'k--', linewidth=lw)
        plotting.plt.plot([1.0,1.5], [1.3,1.9], 'k--', linewidth=lw)
        plotting.plt.plot([-0.5,1.0], [1.3,1.3], 'k--', linewidth=lw)
        if shade is True:
            plotting.plt.fill_between([-0.5,1.0], [1.3,1.3], [2.5,2.5], color='red', alpha=0.1)
            plotting.plt.fill_between([1.005,1.5], [1.3,1.9], [2.5,2.5], color='red', alpha=0.1)

        return


def plt_colorcolor_multi(chlist,c1, c2, c3, c4, dust=True, data=True, cbar=True,
                         cmap="Blues", color='blue',marksize=100, mark='o', label=None, overplot=False, zcolor_range=None):
    if zcolor_range is None:
        zcolor_range = [4, 0]
    nh = len(chlist)
    cnt = 0
    overplot = overplot
    for ch in chlist:
        if cnt > 0 and overplot is False:
            overplot = True
        if cnt < nh - 1:
            ch.plt_colorcolor(c1, c2, c3, c4,
                              dust=dust, data=False, cbar=False, cmap=cmap, color=color,
                              marksize=marksize, mark=mark, label=None, overplot=overplot, zcolor_range=zcolor_range)
        else:
            ch.plt_colorcolor(c1, c2, c3, c4,
                              dust=dust, data=data, cbar=cbar, cmap=cmap, color=color,
                              marksize=marksize, mark=mark, label=label, overplot=overplot, zcolor_range=zcolor_range)
        if overplot is False:
            overplot = True
        cnt += 1


class ColorHist(object):
    def __init__(self, h):
        print "getting color history from database..."
        self.U, self.V, self.B, self.K, self.J, self.I, \
            self.dustU, self.dustV, self.dustB, self.dustK, self.dustJ, self.dustI,\
            self.z, self.t =\
            h.earliest.property_cascade('AB_U', 'AB_V', 'AB_B', 'AB_K', 'AB_J', 'AB_I',
                                        'dustExt_U', 'dustExt_V', 'dustExt_B', 'dustExt_K', 'dustExt_J', 'dustExt_I',
                                        'z','t')
        self.colors = {}

    def get_color(self, c1, c2, dust = True):
        mags = {'U': self.U, 'V': self.V, 'J': self.J, 'K': self.K, 'I': self.I, 'B': self.B}
        dustmags = {'U': self.dustU, 'V': self.dustV, 'J': self.dustJ, 'K': self.dustK, 'I': self.dustI, 'B': self.dustB}
        color = mags[c1]-mags[c2]
        if dust is True:
            print "calculation dust extinction..."
            dustcolor = dustmags[c1]-dustmags[c2]
            color += dustcolor
        self.colors[c1+c2] = color

    def plt_colorcolor(self, c1, c2, c3, c4, dust=True, data=True, cbar=True,
                       cmap="Blues", color='blue', marksize=100, mark='o', label=None, overplot=False, zcolor_range=None):
        if c1+c2 not in self.colors.keys():
            self.get_color(c1, c2, dust=dust)
        if c3+c4 not in self.colors.keys():
            self.get_color(c3, c4, dust=dust)


        import matplotlib.colors as pltcolors

        if zcolor_range is None:
            redcNorm = pltcolors.Normalize(1./(self.z.max()+1), 1./(self.z.min()+1))
        else:
            redcNorm = pltcolors.Normalize(1./(1+zcolor_range[0]), 1./(1+zcolor_range[1]))

        if data==True:
            plotting.plt.errorbar(CANDELS_MW['VJ'], CANDELS_MW['UV'],
                                  xerr=[CANDELS_MW['VJ-'],CANDELS_MW['VJ+']], yerr=[CANDELS_MW['UV-'],CANDELS_MW['UV+']],
                                  fmt='o-', color='grey', markersize=0, elinewidth=.75, linewidth=2)
            plotting.plt.errorbar(CANDELS_M31['VJ'], CANDELS_M31['UV'],
                                  xerr=[CANDELS_M31['VJ-'],CANDELS_M31['VJ+']], yerr=[CANDELS_M31['UV-'],CANDELS_M31['UV+']],
                                  fmt='o-', color='k', markersize=0, elinewidth=.75, linewidth=2)

            plotting.plt.scatter(CANDELS_MW['VJ'], CANDELS_MW['UV'], c=1./(1+np.array(CANDELS_MW['redshift'])),
                                 norm=redcNorm, cmap='Greys', s=150, marker='D',label='Papovich+ 14 MW', linewidth=1.5, color='k')
            plotting.plt.scatter(CANDELS_M31['VJ'], CANDELS_M31['UV'], c=1./(1+np.array(CANDELS_M31['redshift'])),
                                 norm=redcNorm, cmap='Greys', s=150, marker='^',label='M31', linewidth=1.5, color='k')

        plotting.plt.scatter(self.colors[c1+c2], self.colors[c3+c4],
                             c=1./(1+self.z), norm=redcNorm, cmap=cmap, s=marksize, marker=mark, label=label, color=color)


        if cbar is True:
            cbar = plotting.plt.colorbar(ticks=[0.25, 1./3., 0.5, 1./1.5, 1])
            cbar.set_label('Redshift', fontsize=30, labelpad=2)
            cbar.set_ticklabels(['3','2','1','0.5','0'])

        if overplot is False:
            plotting.plt.xlabel(c1+'-'+c2+' color')
            plotting.plt.ylabel(c3+'-'+c4+' color')

        plotting.plt.legend(loc='upper left',fontsize=30)
        return

    def plt_UVJ_quench_region(self, lw=1.5, shade=True):
        plotting.plt.plot([1.5,1.5], [1.9,2.7], 'k--', linewidth=lw)
        plotting.plt.plot([1.0,1.5], [1.3,1.9], 'k--', linewidth=lw)
        plotting.plt.plot([-0.5,1.0], [1.3,1.3], 'k--', linewidth=lw)
        if shade is True:
            plotting.plt.fill_between([-0.5,1.0], [1.3,1.3], [2.5,2.5], color='red', alpha=0.1)
            plotting.plt.fill_between([1.005,1.5], [1.3,1.9], [2.5,2.5], color='red', alpha=0.1)

        return
