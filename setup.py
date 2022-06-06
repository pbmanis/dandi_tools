from setuptools import setup, find_packages, Extension
import numpy

"""
This module is part of *Dandi (custom, pbmanis)*.

Support::

    NIH grants: 1R01NS128873-01  Kato and Manis
    

Paul B. Manis, 2022
"""

version = '0.9.3'

setup(name='vcnmodel',
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
               # 'model_run=vcnmodel.model_run2:main',
               # 'allgbcivs=vcnmodel.all_gbc_iv:main',
               # 'show_swc=vcnmodel.util.show_swc:main',
               # 'render=vcnmodel.util.render:main',
               # 'plot_sims=vcnmodel.plotters.plot_sims:main',
               # 'datatable=vcnmodel.DataTablesVCN:main',
               # 'hocswcmap = vcnmodel.util.hoc_swc_sectionmap:main',
               ],
      },
      classifiers = [
             "Programming Language :: Python :: 3.6+",
             "Development Status ::  Beta",
             "Environment :: Console",
             "Intended Audience :: Manis Lab",
             "License :: MIT",
             "Operating System :: OS Independent",
             "Topic :: Software Development :: Tools :: Python Modules",
             "Topic :: Brain Initiative :: Neuroscience",
             ],
    )
