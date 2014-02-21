Installing SirIsaac
===================

## Install dependencies

SirIsaac depends on the following packages:

- Python 2.6 or later (not Python 3)
- Scipy
- Matplotlib
- SloppyCell (http://sloppycell.sourceforge.net)

One relatively simple way to install the first three packages 
above is using Sage: http://www.sagemath.org.

As of February 2014, SirIsaac depends on the latest
developer version of SloppyCell.  This software is
available through a CVS repository.
We describe here one typical way to download and 
install SloppyCell.
More information is available here: https://sourceforge.net/p/sloppycell/code/

Typically you will want to install SloppyCell in the Python's
site-packages folder, usually found at

* /lib/python2.7/site-packages/ for Linux
(* /Library/Python/2.7/site-packages for Mac?)
* /sw/lib/python2.7/site-packages/

In the simplest case, you can anonymously 'check out' the latest version 
by changing to the above directory and running

> cd /sw/lib/python2.7/site-packages/
> sudo cvs -z3 -d:pserver:anonymous@sloppycell.cvs.sourceforge.net:/cvsroot/sloppycell co -P SloppyCell

You should see a bunch of files being downloaded.  Next,
build the Fortran libraries by changing into the new 'SloppyCell' 
directory and running setup.py:

> cd SloppyCell
> sudo python setup.py build install --install-lib=..

(prefacing by sudo if necessary for write privileges).  (You will
otherwise get the error 'No module named _daskr' when trying to
import SirIsaac.FittingProblem.)  If this is successful, you'll 
eventually see something like 

> Installed /sw/lib/python2.7/site-packages/SloppyCell-CVS-py2.7-macosx-10.7-x86_64.egg
> Processing dependencies for SloppyCell==CVS
> Finished processing dependencies for SloppyCell==CVS


## Install SirIsaac

SirIsaac is available as a git repository hosted on GitHub.  
To download the code:

> cd /sw/lib/python2.7/site-packages/
> git clone https://github.com/EmoryUniversityTheoreticalBiophysics/SirIsaac.git

This will create a 'SirIsaac' directory in the current
location containing the SirIsaac software.

## Test installation

A basic test of the SirIsaac and SloppyCell installation can be
run using

	python SloppyCellTest.py

Code to fit and analyze a simple example dataset 
using SirIsaac is provided in two formats: 
Python (simpleExample.py) and an iPython 
notebook (simpleExample.ipynb).  The 
iPython notebook opens in a web browser and 
includes plots in an interactive format.  To 
open the .ipynb file, run:
    
    ipython notebook --pylab=inline simpleExample.ipynb

To run the .py file in iPython at the command line, run:

    ipython --pylab
    %run simpleExample.py
    show()
