#!/usr/bin/env python3

import os
from setuptools import setup, Extension


# Should fail if the readme is missing
long_des = open('README.rst', 'r').read()

setup(name='pynmrstar',
      version='alpha',
      packages=['pybmrb2'],
      install_requires=['requests>=2.21.0,<=3','plotly>=4.1.0','pynmrstar>=3.0.4','numpy>1.15'],
      python_requires='>=3.6',
      author='Kumaran Baskaran',
      author_email='baskaran@uchc.edu',
      description='PyBMRB provides tools for visializing BMRB entries and  NMR-STAR files. '
                  'Maintained by the BMRB.',
      long_description=long_des,
      long_description_content_type='text/x-rst',
      keywords=['bmrb', 'nmr', 'nmrstar', 'biomagresbank',
                'biological magnetic resonance bank',
                'HSQC','TOCSY','Chemical shifts',
                'Histogram'],
      url='https://github.com/uwbmrb/PyBMRB',
      license='MIT',
      package_data={'pynmrstar': ['reference_files/schema.csv',
                                  'reference_files/comments.str',
                                  'reference_files/data_types.csv']},
      classifiers=[
          'Development Status :: 3 - Alpha',
          'Environment :: Console',
          'Programming Language :: Python :: 3 :: Only',
          'Programming Language :: Python :: 3.6',
          'Programming Language :: Python :: 3.7',
          'Programming Language :: Python :: 3.8',
          'Programming Language :: Python :: 3.9',
          'Intended Audience :: Education',
          'License :: OSI Approved :: MIT License',
          'Natural Language :: English',
          'Operating System :: MacOS',
          'Operating System :: POSIX :: Linux',
          'Operating System :: Microsoft :: Windows',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'Topic :: Software Development :: Libraries',
          'Topic :: Software Development :: Libraries :: Python Modules'
      ]
      )