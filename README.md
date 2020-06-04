SirIsaac
========

Automated dynamical systems inference.

Main goal:
Given experimental dynamical systems trajectories, find a dynamical system that can predict future trajectories.


References
==========

An example of the SirIsaac algorithm applied to experimental data appears in the following publication:

* Daniels, B. C., Ryu, W. S., & Nemenman, I. (2019).  Automated, predictive, and interpretable inference of _Caenorhabditis elegans_ escape dynamics.  Proc. Natl. Acad. Sci. USA.  
https://doi.org/10.1073/pnas.1816531116

Details of the theory and rationale behind the SirIsaac approach are described here:  

* Daniels, B. C., & Nemenman, I. (2015). Automated adaptive inference of phenomenological dynamical models. Nature Communications, 6, 8133.  
https://doi.org/10.1038/ncomms9133

* Daniels, B. C., & Nemenman, I. (2015). Efficient Inference of Parsimonious Phenomenological Models of Cellular Dynamics Using S-Systems and Alternating Regression. Plos One, 10(3), e0119821.  
https://doi.org/10.1371/journal.pone.0119821



Dependencies
============

Python 2.6 or later (not Python 3)  
Scipy  
Matplotlib  
(One way to install the above is with Anaconda or Sage.  See Installation.md.)

SloppyCell (https://github.com/GutenkunstLab/SloppyCell)  


Optional dependencies
=====================

mpi4py (for running on multiple processors)  
SBML (systems biology markup language)  
BioNetGen  
Pygraphviz (for creating network diagrams)  
ipython (for reading ipython notebook file describing example usage)  


Contributors
============

Bryan Daniels, Ilya Nemenman





