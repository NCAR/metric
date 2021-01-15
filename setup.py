#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from os.path import exists

from setuptools import find_packages, setup

if exists('README.rst'):
    with open('README.rst') as f:
        long_description = f.read()
else:
    long_description = ''


with open('requirements.txt') as f:
    install_requires = f.read().strip().split('\n')

CLASSIFIERS = [
    'Development Status :: 4 - Beta',
    'License :: OSI Approved :: Apache Software License',
    'Operating System :: OS Independent',
    'Intended Audience :: Science/Research',
    'Programming Language :: Python',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7',
    'Topic :: Scientific/Engineering',
]



setup(
    name='metric',
    packages= find_packages(),
    version = '0.1.0',
    description='Calculates observational-style decomposition of AMOC using '
                'output from an ocean general circulation model.',
    author="Fred Castruccio",
    author_email="fredc@ucar.edu",
    install_requires=install_requires,
    python_requires='>3.5',
    license='GNU Lesser General Public License',
    long_description=long_description,
    classifiers=CLASSIFIERS,
    url="",
    package_dir={'metric': 'metric'},
    package_data={'metric': ['etc/*']},
    zip_safe=False,
)
