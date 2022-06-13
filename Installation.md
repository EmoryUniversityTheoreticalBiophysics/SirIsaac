Installing SirIsaac
===================

## Install dependencies

SirIsaac depends on the following packages:

- Python 3
- Scipy
- Numpy
- Matplotlib
- SloppyCell (https://github.com/GutenkunstLab/SloppyCell)

There exist several relatively simple ways to install
the first four packages above at once, including

- Anaconda: http://www.anaconda.com
- Sage: http://www.sagemath.org

These systems also have the added benefit of coming
prepackaged with other useful software such as
iPython and Jupyter.

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

A basic suite of unit tests can be run by moving into the SirIsaac directory and running

    python -m unittest

This should take about a minute to run.  Note that there may be lots of warnings and compiler optimization messages (we're working on it...), but if you see something like the following then all tests have passed: 

	---------------------------------------------
	Ran 9 tests in 56.687s
	
	OK

To further help you get up and running, 
code to fit and analyze a simple example dataset 
using SirIsaac is provided in two formats: 
a Jupyter iPython 
notebook (simpleExample.ipynb) and a simple Python script (simpleExample.py).  The 
iPython notebook opens in a web browser and 
includes plots in an interactive format.  To 
open the .ipynb file, run:
    
    jupyter notebook simpleExample.ipynb

To run the .py file in iPython at the command line, run:

    ipython --pylab
    %run simpleExample.py
    plt.show()

## Running in parallel

To run parameter fitting in parallel, you will need to install mpi4py, which further depends on an installation of mpi (often OpenMPI, depending on your operating system).  For more information about installing mpi4py, see here: https://pypi.org/project/mpi4py/

Following are notes from one user detailing OpenMPI installation—these are somewhat outdated (using Python 2.7) but could potentially be useful:

	OpenMPI Installation
	
	Download the latest version(4.1.1) from https://www.open-mpi.org/software/ompi/v4.1/
	Refer the Building MPI from sources section for the installation.
	
	Another source for OpenMPI Installation - https://gist.github.com/mrosemeier/088115b2e34f319b913a
	
	Other Installations that were done for Ubuntu
	
	pip install mpi4py
	
	sudo apt-get install python-dev  \
	     build-essential libssl-dev libffi-dev \
	     libxml2-dev libxslt1-dev zlib1g-dev \
	     
	
	sudo apt install libopenmpi-dev
	
	sudo apt-get install python2.7-dev
	sudo apt-get install build-essential
	sudo apt-get install gcc
	sudo apt-get install python-dev gcc
	sudo apt-get install python2-dev build-essential gcc libpq-dev
	sudo apt-get install libblas-dev libatlas-base-dev
	sudo apt-get install build-essential gcc gfortran git
	sudo apt install gfortran


