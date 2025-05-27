#!/usr/bin/env python3

from setuptools import setup

setup(
    name='cischecker',
    version='0.1.0',
    description='Python tools to check if two mutations in the same gene are located on the same allele (in cis) or not (in trans) using next-generation sequencing data.',
    author='Yuki Saito',
    author_email='js3050ys1990@gmail.com',
    url='---',
    package_dir={'': 'lib'},
    packages=['cischecker'],
    scripts=['cischecker'],
    license='GPL-3.0'
)
