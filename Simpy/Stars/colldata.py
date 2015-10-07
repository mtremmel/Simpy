import numpy as np


# Duncan 14
duncan_z = np.array([4, 5, 6, 7])
duncan_sf = np.array([-1.14, -1.33, -1.58, -1.78])

# Kistler 13:
kist_sfr = np.array([0.0653, 0.03, 0.041, 0.0276, 0.025])
kist_zhigh = [4.5, 5.5, 6.75, 8, 9.4]
kist_zhighplus = np.array([0.5, 0.5, 0.75, 0.5, 1.025]) + 1
kist_zhighminus = np.array([0.5, 0.5, 0.75, 0.5, 1.025]) + 1
kist_sfrplus = np.array([0.0653, 0.03, 0.0405, 0.0647, 0.058])
kist_sfrminus = np.array([0.0326, 0.015, 0.023, 0.023, 0.021])
kist_logsfrplus = np.log10(kist_sfrplus + kist_sfr) - np.log10(kist_sfr)
kist_logsfrminus = np.abs(np.log10(kist_sfr - kist_sfrminus) - np.log10(kist_sfr))


def genCSFRfit(z, z0, A, B, C):
	return C / (10 ** (A * (z - z0)) + 10 ** (B * (z - z0)))


def CSFRFit(z, type='beh'):
	if type != 'beh' and type != 'hop':
		print "WARNING, type of SFR fit not understood (must be either hop or beh). Assuming beh(roozi)"
		type = 'beh'
	if type == 'beh':
		# Behroozi 13
		z0 = 1.243
		C = 0.18
		A = -0.997
		B = 0.241
	if type == 'hop':
		# Hopkins 06
		z0 = 0.840
	C = 0.143
	A = -1.311
	B = 0.085
	sigma = np.zeros(len(z))
	sigma[(z <= 0.9)] = 0.13
	sigma[((0.9 < z) & (z <= 1.5))] = 0.17
	sigma[((1.5 < z) & (z <= 3))] = 0.19
	sigma[(3 < z)] = 0.27
	return genCSFRfit(z, z0, A, B, C), sigma
