#! /usr/bin/env python
from setuptools import setup, find_packages

version = '0.2.0'
dl_version = 'master' if 'dev' in version else '{}'.format(version)

setup(
    name='WisecondorX',
    version=version,
    author='Matthias De Smet, Lennart Raman',
    author_email='Lennart.raman@ugent.be',
    description="WisecondorX -- an evolved WISECONDOR",
    long_description=__doc__,
    keywords=['bioinformatics', 'biology', 'sequencing', 'NGS', 'next generation sequencing',
              'CNV', 'SWGS', 'Shallow Whole Genome Sequencing'],
    download_url='https://github.com/CenterForMedicalGeneticsGhent/WisecondorX/archive/v{}.tar.gz'.format(dl_version),
    license='Attribution-NonCommercial-ShareAlike CC BY-NC-SA',
    packages=find_packages('.'),
    python_requires='==2.7',
    include_package_data=True,
    zip_safe=False,
    install_requires=[
        'scipy',
        'scikit-learn',
        'pysam',
        'numpy'
    ],
    entry_points={
        'console_scripts': ['WisecondorX = wisecondorX.main:main']
    },
    scripts=['scripts/fix_convert_npz.py'],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ]
)
