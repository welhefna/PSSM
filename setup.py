# -*- coding: utf-8 -*-

from setuptools import setup, find_packages


with open('README.rst') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='modb',
    version='0.0.1',
    description='Protein PSSM Parser',
    long_description=readme,
    author='Wessam Elhefnawy',
    author_email='welhe001@odu.edu',
    url='http://cs.odu.edu/~welhefna',
    license=license,
    packages=find_packages(exclude=('tests', 'docs'))
)