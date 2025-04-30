from setuptools import setup

setup(
    name='memelite',
    version='0.1.0',
    author='Jacob Schreiber',
    author_email='jmschreiber91@gmail.com',
    packages=['memelite'],
    #scripts=[''],
    url='https://github.com/jmschrei/memesuite-lite',
    license='LICENSE.txt',
    description='Biological sequence analysis for the modern age.',
    install_requires=[
        "numpy >= 1.14.2, <= 2.0.1",
		"numba >= 0.55.1",
        "pandas >= 1.3.3",
        "pyfaidx >= 0.7.2.1",
        "tqdm >= 4.64.1"
    ],
)