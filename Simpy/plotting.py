import matplotlib.pyplot as plt
import numpy as np
plt.ion()
plt.rc('xtick',labelsize=15)
plt.rc('ytick',labelsize=15)
plt.rc('font', weight='medium', size=15)
plt.rc('axes', linewidth=2)
plt.rc('xtick.major',width=2, pad=10)
plt.rc('ytick.major',width=2, pad=5)
plt.rc('text', usetex=True)


def two_d_grid(data,tick_array1,tick_array2,color_norm,cmap,noticks=True,cbar=True,cbarlabel=None):
	plt.imshow(data,norm=color_norm,cmap=cmap)
	if noticks is True:
		plt.xticks([])
		plt.yticks([])
	else:
		plt.xticks(np.arange(len(data[0])),tick_array1)
		plt.yticks(np.arange(len(data[:,0])),tick_array2)
	if cbar is True:
		cbar = plt.colorbar()
		if cbarlabel:
			cbar.set_label(cbarlabel, fontsize=30, labelpad=2)
	return

def name_ticks_2d(xnames,ynames):
	plt.xticks(np.arange(len(xnames)),xnames)
	plt.yticks(np.arange(len(xnames)),ynames[::-1])

def make_two_d_plot(data1, bins1, data2, bins2, zdata, cmap, zdata_function = np.mean,
					zdata_range=None, makeplot=True, return_data=False, cbar=True, cbarlabel=None):

	image = np.zeros((len(bins1)-1,len(bins2)-1))

	for i in range(len(bins1[:-1])):
		for j in range(len(bins2[:-1])):
			obin = np.where((data1 < bins1[i+1]) & (data1 > bins1[i]) & (data2 < bins2[j+1]) & (data2 > bins2[j]))
			image[i,j] = zdata_function(zdata[obin])

	if makeplot is True:
		import matplotlib.colors as colors
		if zdata_range:
			norm = colors.Normalize(zdata_range[0],zdata_range[1])
		else:
			norm = colors.Normalize(image.min(),image.max())

		two_d_grid(image,bins1[:-1],bins2[:-1],norm,cmap,noticks=True,cbar=cbar,cbarlabel=cbarlabel)

	if return_data is True:
		return image

	else:
		return
