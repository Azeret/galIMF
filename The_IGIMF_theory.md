## The IGIMF theory

Last update: 6.6.2019

The Integrated-Galactic-IMF (IGIMF) theory, originally formulated by [Pavel Kroupa](https://astro.uni-bonn.de/~pavel/) and Carsten Weidner, is the foundation of this work. Here we record the outline of the theory and related references for the reader.


#### Contents
- [Idea](https://github.com/Azeret/galIMF#idea)
- [Method](https://github.com/Azeret/galIMF#method)
- [Development](https://github.com/Azeret/galIMF#development)
- [Validity](https://github.com/Azeret/galIMF#validity)
- [Controversy](https://github.com/Azeret/galIMF#controversy)


### Idea

The fundamental insight underlying the IGIMF theory is that the empirical systematic variation of the IMF of galaxies (galaxy-wide IMF), which appears to correlate with galactic SFR and metallicity, has its origin from the variation of the IMF on a molecular-cloud core, i.e., embedded star cluster scale. In other words, there exists a universal law of the star-cluster-scale IMF shape which leads to the various IMF shapes of different composite systems.

Thus, the IGIMF theory is the link between the star-cluster-scale IMF and the galaxy-wide IMF.


### Method

1. Empirically determine the star-cluster-scale IMF and its variation (as a function of, e.g., star cluster mass and metallicity).

2. Empirically determine the initial star cluster mass distribution function (a.k.a., Embedded star Cluster Mass Function, ECMF) of a galaxy.

Add assumption: The star-cluster-scale IMF, as well as its variation law, is universal and invariant.

3. Integrate the star cluster IMF with ECMF to achieve the galaxy-wide IMF.

4. Compare the galaxy-wide with the observation of star-forming galaxies.

Note that one should compare the galaxy-wide IMF with the recently formed stellar population (luminosity dominated by massive stars). Old population (lower-mass stars) live for a long time and are a composition of multiple star formation event thus multiple galaxy-wide IMF. The comparison of low-mass stars can only be performed in galaxy evolution simulation taking into consideration of the star formation history.


### Development

[Kroupa & Weidner (2003)](https://ui.adsabs.harvard.edu/abs/2003ApJ...598.1076K/abstract) first point out that almost all star forms in a star cluster, thus *PREDICT*s that the galactic-field IMF is the summing up the stellar IMFs of all the star clusters.
It also *SUGGEST*s that the most massive stellar mass in a star cluster, m<sub>max</sub>, is determined by the total mass of the embedded star cluster, M<sub>ecl</sub>, following a m<sub>max</sub>—M<sub>ecl</sub> relation.

[Weidner et al. (2004)](http://adsabs.harvard.edu/abs/2004MNRAS.350.1503W) *DETERMINE*s a galaxy-wide star formation epoch of about 10 Myr and a star cluster population formed in the epoch with an embedded cluster mass function (ECMF) power-law-index, ß, of about 2.

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

[Weidner et al. (2013)](http://adsabs.harvard.edu/abs/2013MNRAS.436.3309W) *INCORPORATE* an empirical ß—SFR relation

[Kroupa et al. (2013)](http://adsabs.harvard.edu/abs/2013pss5.book..115K); [Schulz et al. (2015)](http://adsabs.harvard.edu/abs/2015A%26A...582A..93S) develop the optimal sampling formulation.

[Weidner et al. (2013)](http://adsabs.harvard.edu/abs/2013MNRAS.436.3309W); [Yan et al. (2017)](http://adsabs.harvard.edu/abs/2017A%26A...607A.126Y); [Zhang et al. (2018)](http://adsabs.harvard.edu/abs/2018Natur.558..260Z): **VERIFIED** the systematic variation of the high-mass gwIMF slope;

[Randriamanakoto et al. (2013)](http://adsabs.harvard.edu/abs/2013ApJ...775L..38R); [Yan et al. (2017)](http://adsabs.harvard.edu/abs/2017A%26A...607A.126Y) **VERIFIED** the M<sub>ecl,max</sub>—gwSFR relation .

[Yan et al. (2017)](http://adsabs.harvard.edu/abs/2017A%26A...607A.126Y) realise the optimal sampling of the entire galaxy; GalIMF open Python code;

[Jerábková et al. (2018)](http://adsabs.harvard.edu/abs/2018A%26A...620A..39J) *INCORPORATE* the α<sub>1 & 2</sub>—M<sub>ecl</sub> relation given by Marks+ (2012).

[Jerábková et al. (2018)](http://adsabs.harvard.edu/abs/2018A%26A...620A..39J) point out that the suggested low-mass gwIMF variation may be explained by the IGIMF theory.

Yan et al. (2019, in prepare) develop a galaxy evolution model coupling the IGIMF theory.



### Validity

The important insight given by the IGIMF theory is to link and test statistically the consistency between the independent observations of the star cluster scale and galaxy-scale IMF variations. 

It is conceivable that the star cluster IMF should depend on the temperature, metallicity, and density of the pre-cluster gas cloud.
However, the measurement of this dependency is difficult as it involves dynamical simulation of the stellar system that constantly forming new stars, 
possibly with some level of initial mass segregation (with massive stars forming in the center of the cloud) and losing 
stars through dynamical ejection and cluster expansion that leads to a preferential lost of low-mass stars when they pass the virial radius. 
Thus the constraints on star-cluster IMF variation suffers these systematic uncertainties.

On the other hand, galaxy-wide IMF measurements are free from the dynamical evolution and unresolved multiplicity of single stars
as they fit the integrated emission from the entire galaxy. But the interpretation of the galactic spectrum and the IMF 
sensitive spectral features involves the complicated and uncertain stellar evolution model and stellar atmosphere model.
In addition, the spectra of the galaxy are influenced strongly by the star formation history, mean stellar metallicity/age, and the dust extinction from the galaxy thus the interpretation of the gwIMF variation cannot easily be claimed to be the last solution. 

Knowing the above, it is extremely encouraging to see that the independent set of star cluster scale and galaxy-scale observations, 
involving entirely different method and performed on different scales can be linked and explained consistently by integrating the star-cluster IMF to form the galaxy-wide IMF as formulated 
by the IGIMF theory (see [Yan et al. 2017](http://adsabs.harvard.edu/abs/2017A%26A...607A.126Y)). 


### Controversy

There are several controversies of the IGIMF theory. For example, one of the main objection comes from the seemingly isolated massive stars which have now been proved to be 
dynamically ejected from the star clusters ([Pflamm-Altenburg & Kroupa 2010](http://adsabs.harvard.edu/abs/2010MNRAS.404.1564P); 
[Gvaramadze et al. 2012](http://adsabs.harvard.edu/abs/2012MNRAS.424.3037G); 
[Stephens et al. 2017](http://adsabs.harvard.edu/abs/2017ApJ...834...94S)). 
Here we mention two other apparent falsifications of the IGIMF theory in publications and explain why they are incorrect.
    
The current formulation of the IGIMF theory consider the observed m<sub>max</sub>—M<sub>ecl</sub> relation as one of the indications that the star formation process is highly self-regulated and 
the mass distribution of a stellar population closely follows the IMF with an intrinsic variation far smaller than a 
random sampling procedure would give [(Kroupa et al. 2013)](http://adsabs.harvard.edu/abs/2013pss5.book..115K). This suggests that all the less massive stars are likely to also follow closely with the IMF and thus the M<sub>ecl</sub> as a summation of them. Then the observed 
m<sub>max</sub>—M<sub>ecl</sub> relation is nothing but a representative of the second most massive star—M<sub>ecl</sub> 
relation, third most massive star—M<sub>ecl</sub> relation and so on and so forth to the lowest mass star—M<sub>ecl</sub> 
relation. Thus it is possible to assert the initial mass of every star in a star cluster according to the IMF without going into the random sampling procedure for the stellar masses and achieving a better result. This sampling method is the optimal sampling. With this interpretation, the stellar mass of every star is determined once we know the M<sub>ecl</sub>. It happens that some studies 
take the m<sub>max</sub>—M<sub>ecl</sub> relation as an upper limit of the stellar mass in a star cluster with given mass 
and then apply a m<sub>max</sub> limited random sampling to form the stellar masses of all the stars. This would only result in a lower m<sub>max</sub> for a given star cluster than the observed value and is inconsistent with the IGIMF theory.

Another test has been performed on the IGIMF theory is the [Fe/H] distribution of the stars in comparison with the observed distribution of giants in the Galactic bulge and solar neighbors. An apparent mismatch was found by 
[Ballero et al. (2007)](http://adsabs.harvard.edu/abs/2007A%26A...467..117B) but it does not stand as a falsify of the IGIMF theory for two reasons. Firstly, 
the employed IGIMF formulation is an early version where the top-heavy IMF was not in place. 
The massive stars from a top-heavy IMF could already contribute a decent portion of the iron production thus the development of the IGIMF formulation could not be ignored without a detailed examination. Second but more importantly, the simulation performed is not fully self-consistent. When considering different environments that are able to cause different integrated IMFs, the IMF variation is no longer the only variable in the chemical evolution model since the vast difference in the stellar density should also lead to a 
difference in the SNIa rate. This variation is omitted when applying, e.g., the description of the SNIa rate in 
[Matteucci & Recchi (2001)](http://adsabs.harvard.edu/abs/2001ApJ...558..351M) which depends only on the number of potential progenitors but not their spatial density. 
With the density taken into consideration, the [Fe/H] distribution of stars in Galactic bulge should peak at a higher value in comparison with the solar neighborhood stars. The exact effect of this correction is not to be declared without a reasonable N-body simulation but it certainly releases the tension of the mismatch.


