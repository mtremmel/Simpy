import numpy as np
import pynbody

csq = pynbody.array.SimArray((2.998e10) ** 2, 'erg g**-1')

class BHhalocat(object):

    def __init__(self, simname, bhorbit, er = 0.1):

        Files.cklists(simname)
        f = open('steps.list', 'r')
        steps = np.array(f.readlines).astype(np.int64)
        self.steps = steps
        self.time = np.zeros(len(steps))

        import halo_db as db
        sim = db.get_simulation(simname)
        dbstep = sim.timesteps[0]

        self.data = {'lum':[], 'dist':[], 'x':[], 'y':[], 'z':[], 'halo':[], 'mass':[]}
        cnt = 0
        while dbstep:
            self.time[cnt] = dbstep.time_gyr
            bhh = dbstep.gather_property('BH')
            luminosity = np.array([])
            x = np.array([])
            y = np.array([])
            z = np.array([])
            mass = np.array([])
            dist = np.array([])
            halonum = np.array([])
            for bhl in bhh:
                if type(bhl)!=list:
                    bhl = [bhl]
                mdotarr = pynbody.array.SimArray([bh['BH_mdot'] for bh in bhl],'Msol yr**-1')
                distarr = pynbody.array.SimArray([bh['BH_central_distance'] for bh in bhl],'kpc')
                iord = np.array([bh.halo_number for bh in bhl])
                hnumarr = np.array([bh.host_halo.halo_number for bh in bhl])
                massarr = np.array([bh['BH_mass'] for bh in bhl])

                luminosity = np.append(luminosity, mdotarr.in_units('g s**-1')*er*csq)



