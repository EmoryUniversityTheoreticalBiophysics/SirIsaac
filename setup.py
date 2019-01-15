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
        'multi-processor': ['pypar'],
        'graphs': ['pygraphviz'],
        'SBML': ['python-libsbml'],
    },
    dependency_links=[
        'https://github.com/daleroberts/pypar/tarball/master#egg=pypar',
        'https://github.com/GutenkunstLab/SloppyCell/tarball/master#egg=SloppyCell',
    ]
)
