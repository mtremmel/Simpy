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
