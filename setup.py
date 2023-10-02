from setuptools import find_packages
from setuptools import setup

setup(
    name='mwSuMD',
    version='1.18.9',
    description='multiple walker Supervised Molecular Dynamics',
    author='Giuseppe Deganutti, Ludovico Pipitò',
    author_email='pipitol@uni.coventry.ac.uk',
    python_requires=">=3.6.6",
    packages=find_packages(),
    install_requires=[
        'MDAnalysis',
        'GPUtil',
        'numpy',
        'pandas',
        'setuptools'
    ],
)

