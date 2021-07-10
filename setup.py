#!/usr/bin/env python3

import os
from setuptools import setup, Extension


def get_version():
    internal_file_location = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'pybmrb', '__init__.py')

    with open(internal_file_location, 'r') as internal_file:
        for line in internal_file:
            if line.startswith('__version__'):
                delim = '"' if '"' in line else "'"
                return line.split(delim)[1]
        else:
            raise RuntimeError("Unable to find version string.")

# Should fail if the readme is missing
long_des = open('README.rst', 'r').read()

setup(name='pybmrb',
      version=get_version(),
      packages=['pybmrb'],
      install_requires=['pandas','requests>=2.21.0,<=3','plotly>=4.1.0','pynmrstar>=3.0.4','numpy>1.15'],
      python_requires='>=3.6',
      author='Kumaran Baskaran',
      author_email='baskaran@uchc.edu',
      description='PyBMRB provides tools for visualizing BMRB entries and  NMR-STAR files. '
                  'Maintained by the BMRB.',
      long_description=long_des,
      long_description_content_type='text/x-rst',
      keywords=['bmrb', 'nmr', 'nmrstar', 'biomagresbank',
                'biological magnetic resonance bank',
                'HSQC','TOCSY','Chemical shifts',
                'Histogram'],
      url='https://github.com/uwbmrb/PyBMRB',
      license='MIT',
      classifiers=[
          'Development Status :: 5 - Production/Stable',
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