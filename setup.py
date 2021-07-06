from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='pybmrb',
      version='2.0.1',
      packages=['pybmrb'],
      author='Kumaran Baskaran',
      author_email='baskaran@uchc.edu',
      description='PyBMRB provides tools to visualize chemical shift data in BMRB',
      long_description=long_description,
      long_description_content_type='text/markdown',
      keywords=['bmrb', 'hsqc', 'chemical shift', 'nmrstar', 'biomagresbank', 'biological magnetic resonance bank'],
      url='https://github.com/uwbmrb/PyBMRB',
      package_data={'pybmrb': ['data/*', 'examples/*']},
      install_requires=[
          'pynmrstar>=3.0.4',
          'plotly>=4.5.4',
          'numpy>=1.15.0'],
      license='MIT')
