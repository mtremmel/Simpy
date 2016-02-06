# Reading Black Hole data from ChaNGa simulations

**THE ORBIT FILE**

Simulations run with the newest black hole implementation in ChaNGa output data on *every* black hole at *every* small timestep they experience. This is a lot of data, particularly for big simulations.

To make the data easier to handle, Simpy has scripts that take the raw orbit file data and output a new file with data for every "big" timestep of the simulation (usually 4096 - 8192 steps). To create this file, run the following script

```
Simpy.BlackHoles.orbit.truncOrbitFile(simname, minstep=1, maxstep=4096, ret_output=False, MBHinit=1e6, overwrite=False, newdata=False):
```
Usage:

	Extract raw data from the orbit files, after separating the data out by step to save memory
	Average BH data over each simulation step to make it easier to handle, create new "shortened" orbit file
	Required to create the BHOrbit class object

	simname = name of simulation (i.e. simname.004096)
	minstep, maxstep = min, max steps the code will expect to have in the orbit file.
	ret_output = True/False return the extracted data
	MBHinit = the initial BH mass (as specified in param file)
	overwrite = True/False overwrite shortened.orbit file if it exists. Useful for adding timesteps to an already existing file

This will create an ascii file with the following collumns, one row for each (big) step of the siulation (all values will still be in simulation units):

  BH ID number, time, step, BH mass, x, y, z, vel_x, vel_y, vel_z, mdot, average mdot during step, std dev of mdot during step, scalefactor of step, mass accreted during the step

**Creating an Orbit object**

In order to work with and manipulate the black hole data, create a BHOrbit object

```
bhorbit = Simpy.BlackHoles.orbit.Orbit(simname, savefile='bhorbit.pkl')
```
Usage:
	simname = name of the simulation
	savefile = where to save the object after creation. Pickle formatted

**Utilizing the Orbit class data structure**

The orbit object created above is very useful for analyzing black hole data in a variety of ways. The following are a few basic aspects that are useful.

	bhorbit.data
basically a dictionary with all of the raw data from the shortened orbit file with the addition of a “lum” column, which is just calculated from mdot assuming some constant conversion (0.1 * c^2). In ergs/s. All columns should be arrays that have units associated with it (pynbody.array.SimArray objects)

	bhorbit.bhiords
the ids of all BHs that ever existed in the simulation

	bhorbit.times
the time of all the steps in the simulation

	bhorbit.steps
the steps in the simulation

	bhorbit.single_BH_data(iord, key)
returns an array with the given data (key) for bh with id equal to iord

	bhorbit.single_step_data(step,key)
returns an array with all of the given data (key) for bhs at a given step

	bhorbit.plt_acc_hist(self, style, minM = 1e6, maxM = None, minL = 1e42, maxL = None, type='redshift',xlog=False,ylog=False, label=None, lw=1.5, volume=25**3, plotdata=True, overplot=False)
plots the cumulative accreted mass density over cosmic time for all the black holes in the simulation with luminosity and mass fitting between the max/minM and max/minL values given. Will also plot observed data points if requested.

	bhorbit.plt_lumfun(self, style, minM = 1e6, maxM = None, minL = 1e42, maxL = None, volume=25**3, overplot=False,label=None, bins=50, redshift=1, plotdata=True, dotitle=True,lw=2)
plots the luminosity function of black holes in a given redshift and with masses and luminosities between the given max/min values. 
