from setuptools import find_packages
from setuptools import setup

setup(
    name='mwSuMD',
    version='1.23.6',
    description='multiple walker Supervised Molecular Dynamics',
    author='Giuseppe Deganutti, Ludovico Pipitò',
    author_email='pipitol@uni.coventry.ac.uk',
    python_requires=">=3.6.6",
    packages=find_packages(),
    # the CHARMM36 topologies/parameters live inside the package: without this they are
    # left out of the wheel and an installed mwSuMD finds an empty parameters folder
    package_data={'mwSuMD_lib.parameters': ['*.prm', '*.rtf', '*.str', '*.par', '*.top', '*.param']},
    install_requires=[
        'MDAnalysis',
        'GPUtil',
        'numpy',
        'pandas',
        'setuptools'
    ],
)

