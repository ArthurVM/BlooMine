import subprocess
from glob import glob
from os import listdir
from os.path import splitext, basename, dirname, realpath, exists, join
from setuptools import setup, find_packages
from setuptools.command.install import install as DistutilsInstall


class Install(DistutilsInstall):

    def run(self):
        self._compile_bloomine()

    def _compile_bloomine(self):
        """ build and compile BlooMine from source
        """
        currant_dir = dirname(realpath(__file__))

        subprocess.call("mkdir build; cd build; cmake ..; make; make install", cwd=currant_dir, shell=True)

        DistutilsInstall.run(self)


with open("README", 'r') as f:
    long_description = f.read()

exec(open('./bloomine/_version.py').read())

setup(
    name='BlooMine',
    version=__version__,
    description='Read mining and allele heterogeneity analysis tool.',
    license="MIT",
    long_description=long_description,
    author='Arthur V. Morris',
    author_email='morrisa28@cardiff.ac.uk',
    url="https://github.com/ArthurVM/BlooMine",
    packages=find_packages(''),
    include_package_data=True,
    classifiers=[
        "Intended Audience :: Developers",
        "Operating System :: Unix",
        "Operating System :: POSIX",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: Implementation :: CPython",
        "Programming Language :: Python :: Implementation :: PyPy",
        "Topic :: Utilities",
    ],
    install_requires=[
    'Biopython==1.79',
    'numpy==1.20.3',
    'pysam==0.18.0',
    'pandas==1.3.4',
    'scipy==1.7.1',
    'rich==13.7.1'
    ],
    entry_points={
        'console_scripts': [
            'BlooMine = bloomine.cli:main',
        ],
    },
    # cmdclass={'install': Install},
)
