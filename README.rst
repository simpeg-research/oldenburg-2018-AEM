oldenburg-2018-AEM
===================

.. image:: https://travis-ci.org/simpeg-research/oldenburg-2018-AEM.svg?branch=master
    :target: https://travis-ci.org/simpeg-research/oldenburg-2018-AEM

.. image:: https://mybinder.org/badge.svg
    :target: https://mybinder.org/v2/gh/simpeg-research/oldenburg-2018-AEM/master
    
.. image:: https://zenodo.org/badge/131676406.svg
   :target: https://zenodo.org/badge/latestdoi/131676406


Notebooks and python scripts to reproduce the figures shown in
"3D electromagnetic modelling and inversion: A case for open-source,"
submitted to the AEM 2018 workshop.

.. image:: currents.png
    :width: 70%

**Summary**

Electromagnetics has an important role to play in solving the next generation of geoscience problems. These problems are multidisciplinary, complex, and require collaboration. This is especially true at the base scientific level where the underlying physical equations need to be solved, and data, associated with physical experiments, need to be inverted. In this paper, we present arguments for adopting an open-source methodology for geophysics and provide some background about open-source software for electromagnetics. Immediate benefits are the reduced time required to carry out research, being able to collaborate, having reproducible results, and being able to disseminate results quickly. To illustrate the use of an open-source methodology in electromagnetics, we present two challenges. The first is to simulate data from a time domain airborne system over a conductive plate buried in a more resistive earth. The second is to jointly invert airborne TDEM and FDEM data with ground TDEM. SimPEG, Simulation and Parameter Estimation in Geophysics, (https://simpeg.xyz/) is used for the open-source software. The figures in this paper can be reproduced by downloading the Jupyter Notebooks we provide with this paper (https://github.com/simpeg-research/oldenburg-2018-AEM). Access to the source code allows the researcher to explore the simulations and inversions by changing model and inversion parameters, plot fields and fluxes to gain further insight about the EM phenomena, and solve a new research problem by using open-source software as a base. By providing results in a manner that allows others to reproduce, further explore, and even extend them, we hope to demonstrate that an open-source paradigm has the potential to enable more rapid progress of the geophysics community as a whole.


**Notebooks**

 - `3D_TDEM_simulation_sphere_movie.ipynb <notebooks/3D_TDEM_simulation_sphere_movie.ipynb>`_
 - `3D_TDEM_simulation_topography_movie.ipynb <notebooks/3D_TDEM_simulation_topography_movie.ipynb>`_
 - `Halfspace.ipynb <notebooks/Halfspace.ipynb>`_
 - `Joint EM inversion.ipynb <notebooks/Joint%20EM%20inversion.ipynb>`_

**Usage**

Dependencies are specified in `requirements.txt <https://github.com/simpeg-research/oldenburg-2018-AEM/blob/master/requirements.txt>`_

.. code::

    pip install -r requirements.txt

Please `make an issue <https://github.com/simpeg-research/oldenburg-2018-AEM/issues>`_ if you encounter any problems while trying to run the notebooks.
