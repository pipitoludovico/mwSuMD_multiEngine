from setuptools import find_packages
from setuptools import setup

setup(
    name='mwSuMD',
    version='1.1',
    description='mwSuMD dependencies',
    author='Giuseppe Deganutti, Ludovico Pipitò',
    packages=find_packages(),
    install_requires=[
        'MDAnalysis',
        'GPUtil',
        'numpy',
        'pandas',
        'setuptools'
    ],
)
