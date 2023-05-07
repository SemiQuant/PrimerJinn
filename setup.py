from setuptools import setup

setup(
    name='getMultiPrimerSet',
    version='1.0.1',
    description='A package for generating multiplex PCR primer sets',
    author='Jason D Limberis',
    author_email='Jason.Limberis@ucsf.edu',
    long_description='''
This program designs PCR primers for a multiplex of given target regions of a DNA sequence in a FASTA file.

The script takes several arguments:

| Parameter | Description |
| --- | --- |
| region_file (required) | The path to the primer design region file. This file should have three columns: name, start position, and end position (1-based). It can be in either tsv or xlxs format. |
| input_file (required) | The path to the input FASTA file. |
| target_tm | The desired melting temperature (Tm) for the primers. Default is 60. |
| primer_len | The desired length of the primers. Default is 20. |
| product_size_min | The desired min size for the PCR product. Default is 400. |
| product_size_max | The desired max size for the PCR product. Default is 800. |
| ret | The maximum number of primer pairs to return. Default is 100. |
| Q5 | A boolean indicating whether to use NEB Q5 hotstart polymerase settings for primer3. Default is True. |
| background | The path to the mispriming library FASTA file. Default is an empty string. |
| output | The name of the output file. Default is 'MultiPlexPrimerSet'. |
| eval | The maximum number of primer sets to evaluate. Default is 10000. |

The main function of the script is 'design_primers', which takes the input FASTA file, start and end positions of a target region, and the arguments specified using argparse, and returns the best primer set for that region as a list of dictionaries containing information about the primer pairs.

The 'design_primers' function performs the following steps:

-   Parses the input FASTA file and extracts the target region sequence.
-   Sets up the primer3 input parameters based on the specified arguments.
-   Runs primer3 to design the primers.
-   Extracts information about the primer pairs.
-   Determines the best set of multiplex primers.


### Installation
```
pip install getMultiPrimerSet
```

### Example usage
```
getMultiPrimerSet \
--region_file "./example/primer_regions.tsv" \
--input_file "./example/ref.fasta" \
--target_tm 65 \
--primer_len 20 \
--product_size_min 400 \
--product_size_max 800 \
--ret 100 \
--eval 1000 \
--Q5 \
--background "" \
--output "example"
```
    ''',
    long_description_content_type='text/markdown',
    packages=['getMultiPrimerSet'],
    install_requires=[
        'biopython',
        'pandas',
        'numpy',
        'primer3-py'
    ],
    entry_points={
        'console_scripts': [
            'getMultiPrimerSet=getMultiPrimerSet.getMultiPrimerSet:main',
        ],
    }
)