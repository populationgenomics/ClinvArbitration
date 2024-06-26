"""
automated installation instructions
"""

from setuptools import find_packages, setup

with open('README.md', encoding='utf-8') as handle:
    readme = handle.read()


def read_reqs(filename: str) -> list[str]:
    """
    Read requirements from a file, return as a list

    Args:
        filename (str): the requirements file to parse

    Returns:
        list[str]: the requirements
    """
    with open(filename, encoding='utf-8') as filehandler:
        return [line.strip() for line in filehandler if line.strip() and not line.startswith('#')]


setup(
    name='clinvarbitration',
    description='CPG ClinVar Re-interpretation',
    long_description=readme,
    version='1.1.0',
    author='Matthew Welland, CPG',
    author_email='matthew.welland@populationgenomics.org.au, cas.simons@populationgenomics.org.au',
    url='https://github.com/populationgenomics/ClinvArbitration',
    license='MIT',
    classifiers=[
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    packages=find_packages(),
    include_package_data=True,
    install_requires=read_reqs('requirements.txt'),
    extras_require={'test': read_reqs('requirements-dev.txt')},
    entry_points={
        'console_scripts': [
            # Step 1; re-summarise ClinVar using altered conflict resolution
            'resummary = clinvarbitration.resummarise_clinvar:cli_main',
            # Step 2, post-annotation; obtain PM5 annotations from VEP annotated clinvar
            'pm5_table = clinvarbitration.clinvar_by_codon:cli_main',
        ]
    }
)
