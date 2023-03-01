from setuptools import find_packages
from setuptools import setup

setup(
    name='mwSuMD',
    version='1.1',
    description='mwSuMD dependencies',
    author='Giuseppe Deganutti',
    packages=find_packages(),
    install_requires=[
        'MDAnalysis',
        'mdtraj',
        'moleculekit',
        'numpy',
        'pandas',
        'setuptools'
    ],
)
