from setuptools import setup, find_packages, Extension
import numpy

"""
This module is part of *Dandi (custom, pbmanis)*.

Support::

    NIH grants: 1R01NS128873-01  Kato and Manis
    
Paul B. Manis, 2022
"""

version = '0.1.1'

setup(name='dandi_tools',
      version=version,
      description='Dandi/NWB tools ',
      url='http://github.com/pbmanis/dandi_tools',
      author='Paul B. Manis',
      author_email='pmanis@med.unc.edu',
      license='MIT',
      packages=find_packages(include=['dandi*']),
      # ext_modules=[
      #     Extension(
      #               include_dirs=[numpy.get_include()]),
      #     ],
      python_requires='>=3.9',
      zip_safe=False,
      entry_points={
          'console_scripts': [
               'dispnwb=src.display_nwb:main',
               'acq4tonwb=src.acq4tonwb:main'k
               ],
      },
      classifiers = [
             "Programming Language :: Python :: 3.9+",
             "Development Status ::  Beta",
             "Environment :: Console",
             "Intended Audience :: Manis/Kato Labs",
             "License :: MIT",
             "Operating System :: OS Independent",
             "Topic :: Software Development :: Tools :: Python Modules",
             "Topic :: Brain Initiative :: Neuroscience",
             ],
    )
