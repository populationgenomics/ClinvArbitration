"""
automated installation instructions
"""

from setuptools import find_packages, setup

with open('README.md', encoding='utf-8') as handle:
    readme = handle.read()


setup(
    name='clinvarbitration',
    description='CPG ClinVar Re-interpretation',
    long_description=readme,
    version='2.0.0',
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
    packages=find_packages('src'),
    package_dir={'': 'src'},
    package_data={'clinvarbitration/cpg_flow': ['config_template.toml']},
    include_package_data=True,
    install_requires=[
        'cpg-flow>=v0.1.2',
        'cpg-utils',
        'hail>=0.2.133',
        'pandas>=2.0.3',
        'pyspark>=3.3.3',
    ],
    extras_require={
        'test': [
            'bump2version',
            'black',
            'pre-commit',
            'pytest',
            'pytest-xdist>=3.6.0',
        ],
    },
    entry_points={
        'console_scripts': [
            # sets off the whole workflow in cpg-flow orchestrated Stages
            'run_workflow = clinvarbitration.cpg_flow.run_workflow:cli_main',
            # Step 1; re-summarise ClinVar using altered conflict resolution
            'resummary = clinvarbitration.scripts.resummarise_clinvar:cli_main',
            # Step 2, post-annotation; obtain PM5 annotations from VEP annotated clinvar
            'pm5_table = clinvarbitration.scripts.clinvar_by_codon:cli_main',
        ],
    },
)
