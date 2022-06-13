from setuptools import setup

setup(
    name='SirIsaac',
    setup_requires=[
        'numpy',
    ],
    install_requires=[
        'matplotlib',
        'numpy',
        'scipy',
        'SloppyCell',
    ],
    extras_require={
        'multi-processor': ['mpi4py'],
        'graphs': ['pygraphviz'],
        'SBML': ['python-libsbml'],
    },
    dependency_links=[
        'https://github.com/GutenkunstLab/SloppyCell/tarball/master#egg=SloppyCell',
    ]
)
