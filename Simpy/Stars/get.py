def tform_fromsnap(sim, units):
	try: tform = sim.stars['tform'].in_units(units)
	except:
		try: tform = sim.stars['timeform']
		except:
			print "ERROR! Cannot find formation times for this format!"
			return
	return tform

