from setuptools import setup, find_packages
from setuptools.command.install import install
from pathlib import Path
import os

# read the contents of your README file
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

dependencies = [
    'scipy',
    'pandas',
    'biopython',
    'numpy',
    'primer3-py',
    'networkx',
    'scikit-learn',
    'requests',
    'openpyxl'
]

class CustomInstallCommand(install):
    def run(self):
        install.run(self)
        r=os.system("wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.14.0/ncbi-blast-2.14.0+-x64-linux.tar.gz")
        if r: raise Exception("wget failed")
        r=os.system("tar -xzvf ncbi-blast-2.14.0+-x64-linux.tar.gz")
        if r: raise Exception("tar failed")
        # os.system("export PATH=$PATH:%s" % os.path.abspath("$(pwd)/ncbi-blast-2.14.0+/bin"))
        # - this command has no effect
        # (the 'export' lasts only until system() returned)
        # os.environ["PATH"] += ":"+os.path.abspath("$(pwd)/ncbi-blast-2.14.0+/bin")
        # ^^ this is what you meant if you
        # wanted the new PATH to take effect
        # in any future programs we call here

setup(
    name='primerJinn',
    version='1.0.3',
    url='https://github.com/SemiQuant/primerJinn',
    install_requires=dependencies,
    description='In silico PCR tool',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Jason D Limberis',
    author_email='Jason.Limberis@ucsf.edu',
    keywords=['PCR', 'in silico PCR'],
    license='MIT',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'getMultiPrimerSet=primerJinn.getMultiPrimerSet:main',
            'PCRinSilico=primerJinn.PCRinSilico:main',
            'primerJinn=primerJinn.primerJinn:main'
        ],
    },
    cmdclass={
        'install': CustomInstallCommand,
    }
)
