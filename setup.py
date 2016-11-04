from distutils.core import setup


DESCRIPTION = 'Tools to help analyze black holes in large nbody simulations (in conjunction with pynbody)'
LONG_DESCRIPTION = open('README.md').read()
NAME = 'Simpy'
VERSION = '1.0'
AUTHOR = 'Michael Tremmel'
AUTHOR_EMAIL = 'mjt29@uw.edu'
MAINTAINER = 'Michael Tremmel'
MAINTAINER_EMAIL = 'mjt29@uw.edu'
URL = ''
DOWNLOAD_URL = ''
LICENSE = 'BSD'



setup(name=NAME,
      version=VERSION,
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      author=AUTHOR,
      author_email= AUTHOR_EMAIL,
      maintainer=MAINTAINER,
      maintainer_email=MAINTAINER_EMAIL,
      url=URL,
      download_url=DOWNLOAD_URL,
      license=LICENSE,
      package_dir={'/':''},
      packages=['Simpy',
                'Simpy/BlackHoles',
                'Simpy/Stars',
                'Simpy/Relations',
                'Simpy/dbanalysis'],
#      package_data={'nbdpt':['data/mf_planck13.dat',
 #                            'data/mf_wmap3.dat',
  #                           'data/mf_wmap7.dat']},
      package_data={'Simpy/BlackHoles':['data/QSOdata/M1450z5_McGreer13.dat',
					'data/QSOdata/RhoAccZ.csv',
					'data/QSOdata/RhoAccZMINUS.csv',
					'data/QSOdata/RhoAccZPLUS.csv',
					'data/QSOdata/UnedaRhoBH.csv',
					'data/QSOdata/bol_lf_point_dump.dat',
                    'data/QSOdata/Lacy15_bestfit_rhoBHz.csv',
                    'data/QSOdata/Lacy15_constFaintSlope_rhoBHz.csv',
                    'data/QSOdata/Lacy15_bestfit_epsp06_rhoBHz.csv',
                    'data/Klein16_BHmerge_MBHd.csv',
                    'data/Klein16_BHmerge_MBHnod.csv',
                    'data/SAM_mergers_data/mergers_Q3_delays.dat',
                    'data/SAM_mergers_data/mergers_Q3_short_delays.dat'],
		    'Simpy/Stars':['data/GasFrac/M1450z5_McGreer13.dat',
					'data/GasFrac/RhoAccZ.csv',
					'data/GasFrac/RhoAccZMINUS.csv',
					'data/GasFrac/RhoAccZPLUS.csv',
					'data/GasFrac/UnedaRhoBH.csv,'
					'data/GasFrac/bol_lf_point_dump.dat'
					'data/Morishita15/Morishita15.pkl',
					'data/behrooziData/sfh.pkl']},
      classifiers = ["Development Status :: 3 - Alpha",
                     "Intended Audience :: Developers",
                     "Intended Audience :: Science/Research",
                     "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
                     "Programming Language :: Python :: 2",
                     "Topic :: Scientific/Engineering :: Astronomy"])

