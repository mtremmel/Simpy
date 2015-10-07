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
                'Simpy/Stars'],
#      package_data={'nbdpt':['data/mf_planck13.dat',
 #                            'data/mf_wmap3.dat',
  #                           'data/mf_wmap7.dat']},
      package_data={'Simpy':['DATA/QSOData/M1450z5_McGreer13.dat',
					'DATA/QSOData/RhoAccZ.csv',
					'DATA/QSOData/RhoAccZMINUS.csv',
					'DATA/QSOData/RhoAccZPLUS.csv',
					'DATA/QSOData/UnedaRhoBH.csv',
					'DATA/QSOData/bol_lf_point_dump.dat',
					'DATA/GasFrac/M1450z5_McGreer13.dat',
					'DATA/GasFrac/RhoAccZ.csv',
					'DATA/GasFrac/RhoAccZMINUS.csv',
					'DATA/GasFrac/RhoAccZPLUS.csv',
					'DATA/GasFrac/UnedaRhoBH.csv,'
					'DATA/GasFrac/bol_lf_point_dump.dat'
					'DATA/Morishita15/Morishita15.pkl',
					'DATA/behrooziData/sfh.pkl']},
      classifiers = ["Development Status :: 3 - Alpha",
                     "Intended Audience :: Developers",
                     "Intended Audience :: Science/Research",
                     "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
                     "Programming Language :: Python :: 2",
                     "Topic :: Scientific/Engineering :: Astronomy"])

