Installing SirIsaac
===================

## Install dependencies

SirIsaac depends on the following packages:

- Python 2.6 or later (not Python 3)
- Scipy
- Matplotlib
- SloppyCell (https://github.com/GutenkunstLab/SloppyCell)

There exist several relatively simple ways to install
the first three packages above at once, including

- Anaconda: http://www.anaconda.com
- Sage: http://www.sagemath.org

These systems also have the added benefit of coming
prepackaged with other useful software such as
iPython and Jupyter.

Note: As of June 2020, SirIsaac does not support Python 3 (we are working on it).  Installing Python 2 using the above package managers is likely no longer the default option as it is now officially out of date.

### Install SloppyCell

SirIsaac depends on the latest
version of SloppyCell, available on GitHub.
In the simplest case, you can install with two steps.  First download the git repository by running

    git clone https://github.com/GutenkunstLab/SloppyCell.git
    
which will create a folder named `SloppyCell` in the current directory.  Next, install SloppyCell by running setup.py.  The easiest way to do this is using `pip`:

    pip install -e SloppyCell/

## Install SirIsaac

SirIsaac is similarly available as a git repository on GitHub.   To install, first download by running

    git clone https://github.com/EmoryUniversityTheoreticalBiophysics/SirIsaac.git

which will create a folder named `SirIsaac` in the current directory.  Next, install SirIsaac by running setup.py.   The easiest way to do this is using `pip`: 

	pip install -e SirIsaac/

## Test installation

A basic test of the SirIsaac and SloppyCell installation can be
run by descending into the `SirIsaac/SirIsaac/` directory and running

    python SloppyCellTest.py

More comprehensive tests are found in the `test` subfolder, and can be run using, e.g., `nosetests`.

To further help you get up and running, 
code to fit and analyze a simple example dataset 
using SirIsaac is provided in two formats: 
Python (simpleExample.py) and a Jupyter iPython 
notebook (simpleExample.ipynb).  The 
iPython notebook opens in a web browser and 
includes plots in an interactive format.  To 
open the .ipynb file, run:
    
    jupyter notebook simpleExample.ipynb

To run the .py file in iPython at the command line, run:

    ipython --pylab
    %run simpleExample.py
    show()

