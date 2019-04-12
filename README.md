# GalIMF version 1.1.3

last update: 12.04.2019

## Contents
 - [Overview](https://github.com/Azeret/galIMF#overview)
 - [Scientific motivation](https://github.com/Azeret/galIMF#scientific-motivation)
 - [The IGIMF theory](https://github.com/Azeret/galIMF#the-igimf-theory)
   - [Development](https://github.com/Azeret/galIMF#development)
   - [Validity](https://github.com/Azeret/galIMF#validity)
   - [Controversy](https://github.com/Azeret/galIMF#controversy)
 - [The IGIMF model](https://github.com/Azeret/galIMF#the-igimf-model)
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
 - [The galaxy evolution model](https://github.com/Azeret/galIMF#the-galaxy-evolution-model)
   - [To be published](https://github.com/Azeret/galIMF#to-be-published)
 - [Other information](https://github.com/Azeret/galIMF#other-information)
   - [Versioning](https://github.com/Azeret/galIMF#versioning)
   - [Updates](https://github.com/Azeret/galIMF#updates)
   - [Authors](https://github.com/Azeret/galIMF#authors)
   - [License](https://github.com/Azeret/galIMF#license)
   - [Acknowledgment](https://github.com/Azeret/galIMF#acknowledgment)




 
## Overview

GalIMF stands for the Galaxy-wide Initial Mass Function originated from the Integrated-Galactic-IMF (IGIMF) theory [Kroupa & Weidner (2003)](http://adsabs.harvard.edu/abs/2003ApJ...598.1076K).

GalIMF version 1.0.0 is the original version developed for [Yan, Jerabkova, Kroupa (2017)](http://adsabs.harvard.edu/abs/2017A%26A...607A.126Y). It is a Python 3 module that allows users to compute galaxy-wide initial stellar mass functions based on locally derived empirical constraints following the IGIMF theory (see, e.g., [Weidner et al. 2013](http://adsabs.harvard.edu/abs/2013MNRAS.436.3309W); [Kroupa et al. 2013](http://adsabs.harvard.edu/abs/2013pss5.book..115K)).

GalIMF version 1.1.3 is the version developed for Yan et al. (2019, in prepare). It is a Python 3 module that applies the IGIMF theory on a single-zone closed-box galaxy chemical evolution simulation.

Please check [Versioning](https://github.com/Azeret/galIMF#versioning) below.


## Scientific motivation

The initial stellar mass function (IMF) can be defined as a mass distribution of stars formed during one star formation event in a region of approximately the size of 1 pc. The stellar IMF, therefore, dictates the number of supernova explosions, the chemical enrichment, how bright the stellar population and unresolved objects are and many other issues which affect directly or indirectly a vast majority of astrophysical fields. 

Despite ongoing research, we still cannot formulate and predict the IMF self-consistently. Therefore it is necessary to look at the empirical evidence which might help us to understand the IMF. 

One possible way of doing so is to look at Galactic star forming regions and try to estimate the shape of the IMF as is done by [Marks et al. (2012)](http://adsabs.harvard.edu/abs/2012MNRAS.422.2246M). Then we can look at whole galaxies which can have very different chemical compositions and very different densities and other physical parameters. For a galaxy, we can measure the emitted light and we can try to constrain the total mass and deduce the composite galaxy-wide IMF. For example, one could use the Hα EW and galaxy color as in [Hoversten & Glazebrook (2008)](http://adsabs.harvard.edu/abs/2008ApJ...675..163H) or the ratio between Hα and UV luminosity as in [Lee et al. (2009)](http://adsabs.harvard.edu/abs/2009ApJ...706..599L).

And now here comes the question: 

1. If we take the locally constrained empirical laws and integrate them so we create this galaxy-wide IMF, will we get the same as the locally constrained stellar IMF? If yes, well it would be great and if not we can learn something more about the local IMF based on other galaxies. To help with exactly this problem we present this Python module GalIMF. And indeed, [Yan, Jerabkova, Kroupa (2017)](http://adsabs.harvard.edu/abs/2017A%26A...607A.126Y) demonstrated that the observations are consistent with prediction given by the IGIMF theory.
2. Since the galaxy-wide IMF systematical vary with the galactic properties (see [Development](https://github.com/Azeret/galIMF#development) below),  the galaxy evolution history should be different from the estimates applying the canonical invariant IMF. What would be the influence of such a modification and what are the new implications? This is exactly what we are trying to answer in our new Yan et al. (2019, in prepare) paper.


## The IGIMF theory

### Development

[Kroupa & Weidner (2003)](http://adsabs.harvard.edu/abs/2003ApJ...598.1076K) first point out that almost all star forms in a star cluster, thus *PREDICT*s that the galactic-field IMF is the summing up the stellar IMFs of all the star clusters.
It also *SUGGEST*s that the most massive stellar mass in a star cluster, m<sub>max</sub>, is determined by the total mass of the embedded star cluster, M<sub>ecl</sub>, following a m<sub>max</sub>—M<sub>ecl</sub> relation.

[Weidner et al. (2004)](http://adsabs.harvard.edu/abs/2004MNRAS.350.1503W) *DETERMINE*s a galaxy-wide star formation epoch of about 10 Myr and a star cluster population formed in the epoch with an embedded cluster mass funtion (ECMF) power-law-index, ß, of about 2.

[Weidner & Kroupa (2005)](http://adsabs.harvard.edu/abs/2005ApJ...625..754W) point out that the IGIMF theory is related to the observed M<sub>ecl,max</sub>—galaxy-wide-SFR relation 
and the reduced number of supernovae per star observed in dwarf galaxies.

[Weidner & Kroupa (2006)](http://adsabs.harvard.edu/abs/2006MNRAS.365.1333W); [Weidner et al. (2010)](http://adsabs.harvard.edu/abs/2010MNRAS.401..275W); [(2013)](http://adsabs.harvard.edu/abs/2013MNRAS.434...84W); [(2014)](http://adsabs.harvard.edu/abs/2014MNRAS.441.3348W) **VERIFIED** the m<sub>max</sub>—M<sub>ecl</sub> relation;

[Pflamm-Altenburg+ (2007)](http://adsabs.harvard.edu/abs/2007ApJ...671.1550P) *PREDICT*s the Hα depletion of low-SFR galaxies using the IGIMF theory.

[Hoversten & Glazebrook (2008)](http://adsabs.harvard.edu/abs/2008ApJ...675..163H) **VERIFIED** the top-light gwIMF of low-SFR galaxies using Hα—color signal;

[Pflamm-Altenburg et al. (2009a)](http://adsabs.harvard.edu/abs/2009MNRAS.395..394P); [(2009b)](http://adsabs.harvard.edu/abs/2009ApJ...706..516P); [Lee et al. (2009)](http://adsabs.harvard.edu/abs/2009ApJ...706..599L) **VERIFIED** the top-light gwIMF of dwarf galaxies using Hα/UV signal;

[Pflamm-Altenburg & Kroupa (2010)](http://adsabs.harvard.edu/abs/2010MNRAS.404.1564P); [Gvaramadze et al. (2012)](http://adsabs.harvard.edu/abs/2012MNRAS.424.3037G); [Stephens et al. (2017)](http://adsabs.harvard.edu/abs/2017ApJ...834...94S) **VERIFIED** that the apparent isolated massive stars were born in star clusters.

[Dabringhausen et al. (2009)](http://adsabs.harvard.edu/abs/2009MNRAS.394.1529D); [(2012)](http://adsabs.harvard.edu/abs/2012ApJ...747...72D); [Marks et al. (2012)](http://adsabs.harvard.edu/abs/2012MNRAS.422.2246M) *DETERMINE*s the α<sub>1 & 2 & 3</sub>—M<sub>ecl</sub> relation from an analysis of observed open and globular cluster and ultra-compact dwarf galaxies.

[Weidner et al. (2011)](http://adsabs.harvard.edu/abs/2011MNRAS.412..979W) *INCORPORATE* an observational based α<sub>3</sub>—M<sub>ecl</sub> relation given in [Marks et al. (2012, in prepare at the time)](http://adsabs.harvard.edu/abs/2012MNRAS.422.2246M) and first suggest a possible variation of ß.

[Gunawardhana et al. (2011)](http://adsabs.harvard.edu/abs/2011MNRAS.415.1647G); [Romano et al. (2017)](http://adsabs.harvard.edu/abs/2017MNRAS.470..401R) **VERIFIED** the top-heavy gwIMF of high-SFR galaxies;

[Kroupa et al. (2013)](http://adsabs.harvard.edu/abs/2013pss5.book..115K) gives a comprehensive summary of the IGIMF theory.

[Weidner et al. (2013)](http://adsabs.harvard.edu/abs/2013MNRAS.436.3309W) *INCORPORATE* an observational based ß—SFR relation

[Kroupa et al. (2013)](http://adsabs.harvard.edu/abs/2013pss5.book..115K); [Schulz et al. (2015)](http://adsabs.harvard.edu/abs/2015A%26A...582A..93S) develop the optimal sampling formulation.

[Weidner et al. (2013)](http://adsabs.harvard.edu/abs/2013MNRAS.436.3309W); [Yan et al. (2017)](http://adsabs.harvard.edu/abs/2017A%26A...607A.126Y); [Zhang et al. (2018)](http://adsabs.harvard.edu/abs/2018Natur.558..260Z): **VERIFIED** the systematic variation of the high-mass gwIMF slope;

[Randriamanakoto et al. (2013)](http://adsabs.harvard.edu/abs/2013ApJ...775L..38R); [Yan et al. (2017)](http://adsabs.harvard.edu/abs/2017A%26A...607A.126Y) **VERIFIED** the M<sub>ecl,max</sub>—gwSFR relation .

[Yan et al. (2017)](http://adsabs.harvard.edu/abs/2017A%26A...607A.126Y) realise the optimal sampling of the entire galaxy; GalIMF open Python code;

[Jerábková et al. (2018)](http://adsabs.harvard.edu/abs/2018A%26A...620A..39J) *INCORPORATE* the α<sub>1 & 2</sub>—M<sub>ecl</sub> relation given by Marks+ (2012).

[Jerábková et al. (2018)](http://adsabs.harvard.edu/abs/2018A%26A...620A..39J) point out that the suggested low-mass gwIMF variation may be explained by the IGIMF theory.

Yan et al. (2019, in prepare) develop a galaxy evolution model coupling the IGIMF theory.



### Validity

The important insight given by the IGIMF theory is to link and test statistically the consistency between the independent observations of the star cluster scale and galaxy-scale IMF variations. 

It is conceivable that star cluster IMF should depend on the temperature, metallicity, and density of the pre-cluster gas cloud.
However, the measurement of this dependency is difficult as it involves dynamical simulation of the stellar system that constantly forming new stars, 
possibly with some level of initial mass segregation (with massive stars forming in the center of the cloud) and losing 
stars through dynamical ejection and cluster expansion that leads to a preferential lost of low-mass stars when they pass the virial radius. 
Thus the constrains on star-cluster IMF variation suffers this systematic uncertainties.

On the other hand, galaxy-wide IMF measurements are free from the dynamical evolution and unresolved multiplicity of single stars
as it fit the integrated emission from the entire galaxy. But the interpretation of the galactic spectrum and the IMF sensitive spectral features involves the complicated and uncertain stellar evolution model and stellar atmosphere model.
In addition, the spectra of the galaxy are influenced strongly by the star formation history, mean stellar metallicity/age, and the dust extinction from the galaxy thus the 
interpretation of the gwIMF variation cannot easily be claimed to be the last solution. 

Knowing the above, it is extremely encouraging to see that the independent set of star cluster scale and galaxy scale observations, 
involving entirely different method and performed on different scales can be linked and explained consistently by integrating the star-cluster IMF to form the galaxy-wide IMF as formulated 
by the IGIMF theory (see [Yan et al. 2017](http://adsabs.harvard.edu/abs/2017A%26A...607A.126Y)). 



### Controversy

There are several controversies of the IGIMF theory. For example, one of the main objection comes from the seemingly isolated massive stars which has now been proved to be 
dynamically ejected from the star clusters ([Pflamm-Altenburg & Kroupa 2010](http://adsabs.harvard.edu/abs/2010MNRAS.404.1564P); 
[Gvaramadze et al. 2012](http://adsabs.harvard.edu/abs/2012MNRAS.424.3037G); 
[Stephens et al. 2017](http://adsabs.harvard.edu/abs/2017ApJ...834...94S)). 
Here we mention two other apparent falsifications of the IGIMF theory in publications and explain why they are incorrect.
    
The current formulation of the IGIMF theory consider the observed m<sub>max</sub>—M<sub>ecl</sub> relation as one of the indications that the star formation process is highly self-regulated and 
the mass distribution of a stellar population closely follows the IMF with an intrinsic variation far smaller than a 
random sampling procedure would give [(Kroupa et al. 2013)](http://adsabs.harvard.edu/abs/2013pss5.book..115K). This suggest that all the less massive stars are 
likely to also follow closely with the IMF and thus the M<sub>ecl</sub> as a summation of them. Then the observed 
m<sub>max</sub>—M<sub>ecl</sub> relation is nothing but a representative of the second most massive star—M<sub>ecl</sub> 
relation, third most massive star—M<sub>ecl</sub> relation and so on and so forth to the lowest mass star—M<sub>ecl</sub> 
relation. Thus it is possible to assert the initial mass of every star in a star cluster according to the IMF without 
going into the random sampling procedure for the stellar masses and achieving a better result. This sampling method is 
the optimal sampling. With this interpretation, the stellar mass of every 
star is determined once we know the M<sub>ecl</sub>. It happens that some studies 
take the m<sub>max</sub>—M<sub>ecl</sub> relation as an upper limit of the stellar mass in a star cluster with given mass 
and then apply a m<sub>max</sub> limited random sampling to form the stellar masses of all the stars. This would only 
result in a lower m<sub>max</sub> for a given star cluster than the observed value and is inconsistent with the IGIMF theory.

Another test have been performed on the IGIMF theory is the [Fe/H] distribution of the stars in comparison with the 
observed distribution of giants in the Galactic bulge and solar neighbors. Apparent mismatch was found by 
[Ballero et al. (2007)](http://adsabs.harvard.edu/abs/2007A%26A...467..117B) but it does not stand as a falsify of the IGIMF theory for two reasons. Firstly, 
the employed IGIMF formulation is an early version where the top-heavy IMF was not in place. 
The massive stars from a top-heavy IMF could already contribute a decent portion of the iron production thus the 
development of the IGIMF formulation could not be ignored without a detailed examination. Second but more 
importantly, the simulation performed is not fully self-consistent. When consider 
different environments that are able to cause different integrated IMFs, the IMF variation is no longer the only 
variable in the chemical evolution model since the vast difference in the stellar density should also lead to a 
difference in the SNIa rate. This variation is omitted when applying, e.g., the description of the SNIa rate in 
[Matteucci & Recchi (2001)](http://adsabs.harvard.edu/abs/2001ApJ...558..351M) which depends only on the number of potential progenitors but not their spacial density. 
With the density taken into consideration, the [Fe/H] distribution of stars in Galactic bulge should peak at a higher 
value in comparison with the solar neighborhood stars. The exact effect of this correction is not to be declared 
without a reasonable N-body simulation but it certainly release the tension of the mismatch.



## The IGIMF model

The module calculating the galaxy-wide IMF is described here while more detailed comments can be found in the source code 
together with the support PDF file ([supplementary-document-galimf.pdf](https://github.com/Azeret/galIMF/blob/master/supplementary-document-galimf.pdf)), where all equations are derived in detail and labeled in a consistent way with the source code [galIMF.py](https://github.com/Azeret/galIMF/blob/master/galIMF.py).

An example file, [example_galaxy.py](https://github.com/Azeret/galIMF/blob/master/example_galaxy.py), is provided for a quick test and also serve as an easy entrance for the most basic usage of [galIMF.py](https://github.com/Azeret/galIMF/blob/master/galIMF.py).

GalIMF is also able to optimally sample not an entire galaxy, but only one embedded star cluster with given mass and metallicity. This is demonstrated in [example_star_cluster.py](https://github.com/Azeret/galIMF/blob/master/example_star_cluster.py).


GalIMF represents a Python 3 module which allows computing galaxy-wide IMFs under various assumptions. With the module, we distribute an example script where we use the invariant two-part power-law canonical IMF (Kroupa 2001) as a benchmark and the grid of Salpeter slopes (2.3) drawn into the figures is used for demonstration of IMF variations. 

For the computational details, please, look at [Yan, Jerabkova, Kroupa (2017)](http://adsabs.harvard.edu/abs/2017A%26A...607A.126Y) and [Schulz, Pflamm-Altenburg & Kroupa (2015)](http://adsabs.harvard.edu/abs/2015A%26A...582A..93S).


#### Main features of the module

The generated stellar mass distribution depends on the galaxy-wide star formation rate (**SFR**, which is related to the total mass of a galalxy) and the galaxy-wide **metallicity** ([M/H], see alpha1_model, alpha2_model, and alpha3_model in the code file [galIMF.py](https://github.com/Azeret/galIMF/blob/master/galIMF.py)).

The code can generate a galaxy-wide IMF, i.e., IGIMF. It can also generate all the stellar masses within a galaxy with optimal sampling, i.e., OSGIMF:

* IGIMF in its integrated form

Based on a local IMF (can be the fixed universal Kroupa IMF or the systematically varying IMF based on [Marks et al. 2012](http://adsabs.harvard.edu/abs/2012MNRAS.422.2246M)), GalIMF will produce the galaxy-wide IMF in a data file with contents: stellar mass [Msun] vs. IGIMF values [number of stars Msun^(-1)] normalized to the total mass of a stellar population (see [Yan, Jerabkova, Kroupa 2017](http://adsabs.harvard.edu/abs/2017A%26A...607A.126Y) for details).

* OSGIMF 

The optimally sampled galaxy-wide IMF is based on the same local assumptions as the IGIMF above, however, it represents a discrete formulation of the IMF. That is, its output is in the form of the number of stars in mass bins. The optimal sampling methodology implemented here is from [Schulz, Pflamm-Altenburg & Kroupa (2015)](http://adsabs.harvard.edu/abs/2015A%26A...582A..93S) which generates stars from the IMF without Poison noise (therefore "optimal sampling"). This method presents a probe which is capable of testing the nature of star formation as an alternative to random / stochastic sampling of the IMF.

* Additional functions

To compute the IGIMF or the OSGIMF, the GalIMF module contains all local IMF properties (e.g. the dependence of the stellar IMF on the metallicity, on the density of the star-cluster forming molecular cloud cores), and this software module can, therefore, be also used to obtain only the stellar IMF with various prescriptions, or to investigate other features of the stellar population such as what is the most massive star that can be formed in a star cluster.

* Features of the code: IGIMF and OSGIMF

This script needs two input values: SFR (the star formation rate, in Msun/yr) and the [M/H] value. The output is the IGIMF and OSGIMF as a function of stellar mass, normalized to the total stellar mass. The IGIMF and OSGIMF values are also written into the output file.  To demonstrate the shape of the generated IMF, we include a grid of Salpeter power-law indices into the figures.



### Deployment

For users without any experience with Python, we recommend using [Anaconda](https://www.continuum.io/) to install Python 3 and all required packages. 

For users who already have installed Python 3, we recommend using [pip](https://pip.pypa.io/en/stable/) to install all required packages. An instruction can be found [here](https://www.scipy.org/install.html).

For users having both Python 2 and Python 3 installed, your two Python interpreter may coexisting with the same name. See [this page](https://stackoverflow.com/questions/341184/can-i-install-python-3-x-and-2-x-on-the-same-computer) if it causes any problem.



### Getting Started

In the following subsections, we describe how to install and set up the module. We tested this on MACOSX, Linux and Windows platforms. 

#### Prerequisites

The GalIMF module is written in Python 3, therefore you need to install [Python 3](https://www.python.org/download/releases/3.0/) and the following packages: 

For analyzing and visualize the results as our example script [example_galaxy.py](https://github.com/Azeret/galIMF/blob/master/example_galaxy.py) does, one needs [numpy](http://www.numpy.org/), [scipy](https://www.scipy.org/), [matplotlib](https://matplotlib.org/). 



#### Running the test

To learn how to use the code and to present its main features also to researchers not familiar with Python, we prepared an example implementation of the GalIMF module. This example implementation is called [example_galaxy.py](https://github.com/Azeret/galIMF/blob/master/example_galaxy.py) and is included together with the module.

First, make sure you are using Python 3, then write:
```
python directory_of_example_galaxy/example_galaxy.py
```
into a terminal to run our example program. Further instructions will show up in the terminal.

The [example_galaxy.py](https://github.com/Azeret/galIMF/blob/master/example_galaxy.py) will generate a TXT file and a PDF file in the same directory as the basic output of GalIMF.



### Employ GalIMF for your own program

#### For Python programs

You can download the GalIMF repository and call the [galIMF.py](https://github.com/Azeret/galIMF/blob/master/galIMF.py) module based on the placement in your computer.

If the Module directory is in the same directory as your own Python script (this is the case of the presented example [example_galaxy.py](https://github.com/Azeret/galIMF/blob/master/example_galaxy.py)) you will import [galIMF.py](https://github.com/Azeret/galIMF/blob/master/galIMF.py) as: 
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

#### For non-Python programs

Use a pipeline to first run, e.g., the [example_galaxy.py](https://github.com/Azeret/galIMF/blob/master/example_galaxy.py). Then read in the output files of GalIMF, e.g., GalIMF_IGIMF.txt, for your own program.



### Inputs/parameters

#### Basic inputs

To apply the IGIMF theory on different galaxies, the following parameters should be changed. They are also the required input of our example code [example_galaxy.py](https://github.com/Azeret/galIMF/blob/master/example_galaxy.py):

* galaxy-wide star formation rate: SFR

* galaxy-wide metallicity [M/H]: M_over_H

* OSGIMF stellar mass resolution: resolution


#### Other adjustable parameters

The following inputs are not essential for the IGIMF theory and can be changed according to the research context:

* Resolution parameter

bindw: defined in [example_galaxy.py](https://github.com/Azeret/galIMF/blob/master/example_galaxy.py),
will automatically change the resolution of histograms for optimal sampling.

* IMF model parameters

M_str_L = 0.08:
stellar mass lower limit [Solar mass]

M_str_U = 150:
stellar mass upper limit [Solar mass]

M_turn = 0.5:
first mass at which the power-law index of the stellar IMF changes [Solar mass] (i.e. in the canonical IMF, the IMF power-law index changes from alpha_1=1.3 to alpha_2=2.3 at a stellar mass of 0.5 Msun)

M_turn2 = 1:
second mass at which the power-law index of the stellar IMF changes [Solar mass] (i.e. in the canonical IMF, the IMF power-law index changes from alpha_2=2.3 to alpha_3=2.3 at a stellar mass of 1.0 Msun, i.e., the canonical IMF has a Salpeter index above 0.5 Msun)

alpha3_model = 1:
IMF high-mass-end power-law index model, see Function_alpha_3_change

alpha2_model = 1:
see Function_alpha_2_change

alpha1_model = 1:
see Function_alpha_1_change


* ECMF model parameters

beta_model = 1:
see Function_beta_change

M_ecl_U = 10^9:
upper limit of the embedded cluster mass in stars [Solar mass] 

M_ecl_L = 5:
lower limit of the embedded cluster mass in stars [Solar mass] 


#### Internal parameters of the theory

The following parameters should not be changed as they are part of the IGIMF theory. Read [Yan, Jerabkova, Kroupa (2017)](http://adsabs.harvard.edu/abs/2017A%26A...607A.126Y) carefully if you do intend to change them.

* galaxy star formation assumption

delta_t = 10:
duration of star formation epoch [Myr] (the time-scale on which the ISM forms molecular clouds and from them a new stellar population in embedded clusters)

* Optimal sampling rule

I_ecl = 1:
normalization factor in the optimal sampling condition equation 

I_str = 1:
normalization factor in the optimal sampling condition equation 



## The galaxy evolution model

### To be published

The galaxy evolution model coupling the IGIMF theory will be published in Yan et al. (2019, in prepare).

Detailed explanation and application information of the new model will be released.


## Other information

### Versioning

This site always keep the newest GalIMF version and the old version used in our publications.

We use [SemVer](http://semver.org/) for versioning. That is, 

1. the first digits stand for MAJOR update that is incompatible with the earlier version;
2. the second digits indicate MINOR update when new features or function is added in a backwards-compatible manner;
3. the last digits note the backwards-compatible error corrections (e.g., adding comments, change parameter names, etc.).

See the [GalIMF homepage](https://sites.google.com/view/galimf/home) for an overlook of all the major versions, related scientific publications, and simulations demonstrations.



### Updates

The major updates include:
1. Add example files that demonstrates how to construct star cluster and galaxy-wide IMF as well as getting each stellar mass in the star cluster or the galaxy applying the IGIMF theory with the galIMF.py model.
2. Change the IMF metal dependence parameter from the iron abundance, [Fe/H], to the total metallicity, [M/H], indicating that the IMF variation depend on general or total metallicity instead of solely on the iron abundance. The old version applying [Fe/H] follows the formulation in Marks et al. 2012 (MNRAS.422.2246M) correctly but the author of this paper (through private communication) actually consider the [Fe/H] to be a representative of [M/H]. It makes more sense that all metal element should have a similar effect (if not an identical effect) to the IMF variation.
3. Some function names in the file galIMF.py are changed to lower case letters. This may cause incompatible issues. Please change the function names accordingly or use the galIMF_version_1.0.py instead of galIMF.py.
4. Add the galaxy evolution model, galaxy_evol.py, and corresponding supporting data files. The new model adopt the galaxy-wide IMF for a single 10 Myr star formation epoch calculated by galIMF.py to the galaxy formation and evolution in a 10 Gyr timescale. (01.01.2019)
5. The approximated stellar luminosity weighted results is now available. The "stellar luminosity" adopted are the luminosity of the star during its main-sequence stage and do not consider any stellar evolution, i.e., the luminosity is only a function the stellar initial mass but not its age or metallicity. (10.02.2019)
6. An uniform outflow (uniform in the sense that the element ratios are the same as the well-mixed gas phase) that is proportional to the stellar mass formed is added. It has a minor effect on the final total gas mass (roughly 0.3 dex) and the metal abundances (roughly 0.1 dex) and a negligible effect on the galaxy final metal abundance ratios (roughly 0.05 dex).



### Authors

The main author of the Python program is:

* [Yan Zhiqiang (闫智强)](https://orcid.org/0000-0001-7395-1198),
University of Bonn, Charles University,
yan(at)astro.uni-bonn.de

Other members of this project include:

* [Tereza Jeřábková](http://sirrah.troja.mff.cuni.cz/~tereza/index.html),
Charles University, European Southern Observatory
tereza(at)sirrah.troja.mff.cuni.cz

* Prof. Dr. [Pavel Kroupa](https://astro.uni-bonn.de/~pavel/),
Universty of Bonn, Charles University,
pavel(at)astro.uni-bonn.de



### License and Citation

This program is free software. You can redistribute it and/or modify it. However, you must state all modifications carefully and contain the License file when doing so.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the [LICENSE](https://github.com/Azeret/galIMF/blob/master/LICENSE) file for details.

GalIMF is developed by our group with **a large amount of effort**. If GalIMF contributes to a project that leads to a scientific publication, please acknowledge this work by citing the project. 

When publishing results based on this software or parts of it (the executable and / or the source code) please cite the relevant publications with this the ready-made citation entry: 
1. [Yan, Z., Jerabkova, T., Kroupa, P. (2017)](https://sites.google.com/view/galimf/citation) for the IGIMF model (galIMF.py);
2. and Yan et al. (2019, in prepare) for the galaxy evolution model (galaxy_evol.py).



### Acknowledgment

We thank Vaclav Pavlik for kindly testing the code GalIMF 1.0.
