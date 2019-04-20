#!/usr/bin/env python2

from distutils.core import setup

setup(name='cischecker',
      version='0.1.0',
      description='Python tools to check if two mutations in a same gene are located in the same allele (in cis) or not (in trans) using next-generation sequencing data.',
      author='Yuki Saito',
      author_email='js3050ys1990@gmail.com',
      url='---',
      package_dir = {'': 'lib'},
      packages=['cischecker'],
      scripts=['cischecker'],
      license='GPL-3'
     )
