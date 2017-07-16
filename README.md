# GalIMF version 1.0.0

## Contents

 - [Overview](https://github.com/Azeret/galIMF#overview)
   - [Scientific motivation](https://github.com/Azeret/galIMF#scientific-motivation)
   - [Main features of the module](https://github.com/Azeret/galIMF#main-features-of-the-module)
 - [Deployment](https://github.com/Azeret/galIMF#deployment)
 - [Getting Started](https://github.com/Azeret/galIMF#getting-started)
   - [Prerequisites](https://github.com/Azeret/galIMF#prerequisites)
   - [Running the test](https://github.com/Azeret/galIMF#running-the-test)
 - [Employ GalIMF for your own program](https://github.com/Azeret/galIMF#employ-galimf-for-your-own-program)
    - [For Python programs](https://github.com/Azeret/galIMF#for-the-python-program)
    - [For non-Python programs](https://github.com/Azeret/galIMF#for-a-non-python-program)
 - [Inputs/parameters](https://github.com/Azeret/galIMF#inputsparameters)
   - [Basic inputs](https://github.com/Azeret/galIMF#basic-inputs)
   - [Other adjustable parameters](https://github.com/Azeret/galIMF#other-adjustable-parameters)
   - [Internal parameters of the theory](https://github.com/Azeret/galIMF#internal-parameters-of-the-theory)
 - [Versioning](https://github.com/Azeret/galIMF#versioning)
 - [Authors](https://github.com/Azeret/galIMF#authors)
 - [License](https://github.com/Azeret/galIMF#license)
 - [Acknowledgment](https://github.com/Azeret/galIMF#acknowledgment)

## Overview

GalIMF stands for the Galaxy-wide Initial Mass Function. It is a Python 3 module that allows users to compute galaxy-wide initial stellar mass functions based on locally derived empirical constraints following the IGIMF theory (see, e.g., [Weidner et al. 2013](http://adsabs.harvard.edu/abs/2013MNRAS.436.3309W); [Kroupa et al. 2013](http://adsabs.harvard.edu/abs/2013pss5.book..115K)). The module is described here, but more detailed comments can be found in the source code 
together with the support PDF file ([supplementary-document-galimf.pdf](https://github.com/Azeret/galIMF/blob/master/supplementary-document-galimf.pdf)), where all equations are derived in detail and labeled in a consistent way with the source code.


### Scientific motivation

The initial stellar mass function (IMF) can be defined as a mass distribution of stars formed during one star formation event in a region of approximately the size of 1 pc. The stellar IMF, therefore, dictates the number of supernova explosions, the chemical enrichment, how bright the stellar population and unresolved objects are and many other issues which affect directly or indirectly a vast majority of astrophysical fields. 

Despite ongoing research, we still cannot formulate and predict the IMF self-consistently. Therefore it is necessary to look at the empirical evidence which might help us to understand the IMF. 

One possible way of doing so is to look at Galactic star forming regions and try to estimate the shape of the IMF. Then we can look at whole galaxies which can have very different chemical compositions and very different densities and other physical parameters. For a galaxy, we can measure the emitted light and we can try to constrain the total mass and deduce the composite galaxy-wide IMF. 

And now here comes the question: if we take the locally constrained empirical laws and integrate them so we create this galaxy-wide IMF, will we get the same as the locally constrained stellar IMF? If yes, well it would be great and if not we can learn something more about the local IMF based on other galaxies. To help with exactly this problem we present this Python module GalIMF.   

GalIMF represents a Python 3 module which allows computing galaxy-wide IMFs under various assumptions. With the module, we distribute an example script where we use the invariant two-part power-law canonical IMF (Kroupa 2001) as a benchmark and the grid of Salpeter slopes (2.3) drawn into the figures is used for demonstration of IMF variations. 

For the computational details, please, look at Yan, Jerabkova & Kroupa (2017, A&A, in press) and [Schulz, Pflamm-Altenburg & Kroupa (2015)](http://adsabs.harvard.edu/abs/2015A%26A...582A..93S).



### Main features of the module

* IGIMF in its integrated form

Based on a local IMF (can be the fixed universal Kroupa IMF or the systematically varying IMF based on Marks et al. 2012), GalIMF will produce the galaxy-wide IMF in a data file with contents: stellar mass [Msun] vs. IGIMF values [number of stars Msun^(-1)] normalized to the total mass of a stellar population (see Yan et al. 2017 for details).

* OSGIMF 

The optimally sampled galaxy-wide IMF is based on the same local assumptions as the IGIMF above, however, it represents a discrete formulation of the IMF. That is, its output is in the form of the number of stars in mass bins. The optimal sampling methodology implemented here is from Schulz et al. (2015) which generates stars from the IMF without Poison noise (therefore "optimal sampling"). This method presents a probe which is capable of testing the nature of star formation as an alternative to random / stochastic sampling of the IMF.

* Additional functions

To compute the IGIMF or the OSGIMF, the GalIMF module contains all local IMF properties (e.g. the dependence of the stellar IMF on the metallicity, on the density of the star-cluster forming molecular cloud cores), and this software module can, therefore, be also used to obtain only the stellar IMF with various prescriptions, or to investigate other features of the stellar population such as what is the most massive star that can be formed in a star cluster.

* Features of the code: IGIMF and OSGIMF

This script needs two input values: SFR (the star formation rate, in Msun/yr) and the [Fe/H] value. The output is the IGIMF and OSGIMF as a function of stellar mass, normalized to the total stellar mass. The IGIMF and OSGIMF values are also written into the output file.  To demonstrate the shape of the generated IMF, we include a grid of Salpeter power-law indices into the figures.



## Deployment

For users without any experience with Python, we recommend using [Anaconda](https://www.continuum.io/) to install Python 3 and all required packages. 

For users who already have installed Python 3, we recommend using [pip](https://pip.pypa.io/en/stable/) to install all required packages. An instruction can be found [here](https://www.scipy.org/install.html).

For users already having both Python 2 and Python 3 installed, you have two "Python" interpreter coexisting in your computer. See [this page](https://stackoverflow.com/questions/341184/can-i-install-python-3-x-and-2-x-on-the-same-computer).



## Getting Started

In the following subsections, we describe how to install and set up the module. We tested this on MACOSX, Linux and Windows platforms. 

### Prerequisites

The GalIMF module is written in Python 3, therefore you need to install Python 3 and the following packages: 

For analyzing and visualize the results as our example script [test_gimf.py](https://github.com/Azeret/galIMF/blob/master/test_gimf.py) does, one needs [numpy](http://www.numpy.org/), [scipy](https://www.scipy.org/), [matplotlib](https://matplotlib.org/). 



### Running the test

To learn how to use the code and to present its main features also to researchers not familiar with Python, we prepared an example implementation of the GalIMF module. This example implementation is called [test_gimf.py](https://github.com/Azeret/galIMF/blob/master/test_gimf.py) and is included together with the module.

First, make sure you are using Python 3, then write:
```
python directory_of_test_gimf/test_gimf.py
```
into a terminal to run our example program. Further instructions will show up in the terminal.

The [test_gimf.py](https://github.com/Azeret/galIMF/blob/master/test_gimf.py) will generate a TXT file and a PDF file in the same directory as the basic output of GalIMF.



## Employ GalIMF for your own program

### For Python programs

You can download the GalIMF repository and call the [galIMF.py](https://github.com/Azeret/galIMF/blob/master/galIMF.py) module based on the placement in your computer.

If the Module directory is in the same directory as your own Python script (this is the case of the presented examples) you will import [galIMF.py](https://github.com/Azeret/galIMF/blob/master/galIMF.py) as: 
```python
import galIMF
```
If it is in a different directory, it is also possible to call it from its directory using:
```python
import directory.galIMF
```
The third option is to put [galIMF.py](https://github.com/Azeret/galIMF/blob/master/galIMF.py) into the Python directory structure so that you can easily deploy [galIMF.py](https://github.com/Azeret/galIMF/blob/master/galIMF.py) for your Python 3 project in any directory all as:
```python
import galIMF
```

To do that open Python interpreter and run: 
```python
import sys
sys.path
```
This will locate Python libraries on your computer (usually there is something similar to "...\lib\site-packages"). If GalIMF is placed in this directory the [galIMF.py](https://github.com/Azeret/galIMF/blob/master/galIMF.py) module can be called from the Python script located anywhere simply as:
```python
import galIMF
```
Then you can treat GalIMF as the same as any other packages.

### For non-Python programs

Use a pipeline to first run, e.g., the [test_gimf.py](https://github.com/Azeret/galIMF/blob/master/test_gimf.py). Then read in the output files of GalIMF, e.g., GalIMF_IGIMF.txt for your own program.



## Inputs/parameters

### Basic inputs

To apply the IGIMF theory on different galaxies, the following parameters should be changed. They are also the required input of our example code [test_gimf.py](https://github.com/Azeret/galIMF/blob/master/test_gimf.py):

* galaxy-wide star formation rate: SFR

* galaxy-wide metallicity [Fe/H]: Fe_over_H

* OSGIMF stellar mass resolution: resolution


### Other adjustable parameters

The following inputs are not essential for the IGIMF theory and can be changed according to the research context:

* Resolution parameter

bindw = defined in the module
will automatically change the resolution of histograms for optimal sampling

* IMF model parameters

M_str_L=0.08
stellar lower mass limit [Solar mass]

M_str_U=150
stellar upper mass limit [Solar mass]

M_turn=0.5
first mass at which the power-law index of the stellar IMF changes [Solar mass] (i.e. in the canonical IMF, the IMF power-law index changes from alpha_1=1.3 to alpha_2=2.3 at a stellar mass of 0.5 Msun)

M_turn2=1.
second mass at which the power-law index of the stellar IMF changes [Solar mass] (i.e. in the canonical IMF, the IMF power-law index changes from alpha_2=2.3 to alpha_3=2.3 at a stellar mass of 1.0 Msun, i.e., the canonical IMF has a Salpeter index above 0.5 Msun)

alpha3_model=1
IMF high-mass-end power-law index model, see Function_alpha_3_change

alpha2_model = 1 
see Function_alpha_2_change

alpha1_model = 1 
see Function_alpha_1_change


* ECMF model parameters

beta_model=1
see Function_beta_change

M_ecl_U=10^9 
upper limit of the embedded cluster mass in stars [Solar mass] 

M_ecl_L=5.
lower limit of the embedded cluster mass in stars [Solar mass] 


### Internal parameters of the theory

The following parameters should not be changed as they are part of the IGIMF theory. Read Yan et al. (2017) carefully if you do intend to change them.

* galaxy star formation assumption

delta_t=10
duration of star formation epoch [Myr] (the time-scale on which the ISM forms molecular clouds and from them a new stellar population in embedded clusters)

* Optimal sampling rule

I_ecl = 1.
normalization factor in the optimal sampling condition equation 

I_str=1.
normalization factor in the optimal sampling condition equation 



## Versioning

This is GalIMF version 1.0.0. We use [SemVer](http://semver.org/) for versioning.

For the versions available, see the [GalIMF homepage](https://sites.google.com/view/galimf/home). 



## Authors

The main author of this module is:

* Zhiqiang Yan,
University of Bonn,
beijingjuzikong(at)gmail.com

The other members of this project are:

* Tereza Jerabkova,
University of Bonn, Charles University,
tereza(at)sirrah.troja.mff.cuni.cz

* Prof. Dr. Pavel Kroupa,
Universty of Bonn, Charles University,
pavel(at)astro.uni-bonn.de



## License

This program is free software. You can redistribute it and/or modify it. However, you must state all modifications carefully and contain the License file when doing so.

When publishing results based on this software or parts of it (the executable and / or the source code) please cite: 
Yan, Z., Jerabkova, T., Kroupa, P. 2017, A&A, in press

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the [LICENSE](https://github.com/Azeret/galIMF/blob/master/LICENSE) file for details.



## Acknowledgment

We thank Vaclav Pavlik for kindly testing the module GalIMF.
