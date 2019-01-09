from setuptools import setup, find_packages

from pybmrb.pybmrb import __version__

setup(name='pybmrb',
      version=__version__,
      packages = ['pybmrb'],
      author='Kumaran Baskaran',
      author_email='kbaskaran@bmrb.wisc.edu',
      description='PyBMRB provides tools to visualize chemical shift data in BMRB',
      keywords=['bmrb', 'hsqc', 'chemical shift', 'nmrstar', 'biomagresbank', 'biological magnetic resonance bank'],
      url='https://github.com/uwbmrb/PyBMRB',
      package_data = {'pybmrb':['data/*','examples/*']},
      license='MIT')
