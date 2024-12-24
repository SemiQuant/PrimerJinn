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
    'primer3-py==2.0.3',
    'networkx',
    'scikit-learn',
    'requests',
    'openpyxl'
]

class CustomInstallCommand(install):
    def run(self):
        install.run(self)
        
        # First check if BLAST is already functional
        print("\nChecking for existing BLAST+ installation...")
        try:
            # Try to run blastn -version
            if os.system("blastn -version > /dev/null 2>&1") == 0:
                print("Found working BLAST+ installation. Skipping download.")
                return
        except:
            pass
            
        # If we get here, we need to install BLAST
        print("No working BLAST+ installation found.")
        
        # Determine platform-specific download URL
        platform = os.uname().sysname.lower()
        if platform == 'darwin':
            blast_url = "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.14.0/ncbi-blast-2.14.0+-x64-macosx.tar.gz"
            dirname = "ncbi-blast-2.14.0+"
        elif platform == 'linux':
            blast_url = "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.14.0/ncbi-blast-2.14.0+-x64-linux.tar.gz"
            dirname = "ncbi-blast-2.14.0+"
        else:
            raise Exception(f"Unsupported platform: {platform}")
            
        # Get the installation directory
        install_dir = os.path.join(self.install_lib, 'primerJinn', 'blast')
        os.makedirs(install_dir, exist_ok=True)
        
        print(f"Downloading BLAST+ for {platform}...")
        tarfile = os.path.join(install_dir, "blast.tar.gz")
        r = os.system(f"wget {blast_url} -O {tarfile}")
        if r: raise Exception("wget failed")
        
        # Extract to the correct location
        current_dir = os.getcwd()
        os.chdir(install_dir)
        r = os.system(f"tar -xzf {tarfile}")
        if r: raise Exception("tar failed")
        os.remove(tarfile)
        os.chdir(current_dir)
        
        # Add BLAST to PATH for this session
        blast_bin = os.path.join(install_dir, dirname, "bin")
        os.environ["PATH"] = blast_bin + os.pathsep + os.environ["PATH"]
        
        # Verify installation
        if os.system("blastn -version > /dev/null 2>&1") != 0:
            raise Exception("BLAST+ installation failed verification")
        print("BLAST+ installation successful!")

setup(
    name='primerJinn',
    version='1.1.1',
    url='https://github.com/SemiQuant/primerJinn',
    install_requires=dependencies,
    description='In silico PCR tool',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Jason D Limberis',
    author_email='Jason.Limberis@ucsf.edu',
    keywords=['PCR', 'in silico PCR'],
    license='GPL',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
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
