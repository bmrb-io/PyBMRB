from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='pybmrb',
      version='1.2.97',
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
          'pynmrstar',
          'plotly',
          'numpy'],
      license='MIT')
