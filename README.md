ComGamHarm
================
Adam Lang

-   [Overview](#overview)
-   [Installation](#installation)
-   [ComGamHarm](#comgamharm)
    -   [Arguments](#arguments)
-   [Example](#example)
    -   [Cross-Sectional Data](#cross-sectional-data)
    -   [Harmonize](#harmonize)
-   [Applying Harmonization To New
    Data](#applying-harmonization-to-new-data)
    -   [Example](#example-1)
-   [Plots](#plots)
    -   [Plots Before Harmonization](#plots-before-harmonization)
        -   [ROI1](#roi1)
        -   [ROI2](#roi2)
        -   [ROI3](#roi3)
        -   [ROI4](#roi4)
    -   [Plots After Harmonization](#plots-after-harmonization)
        -   [ROI1](#roi1-1)
        -   [ROI2](#roi2-1)
        -   [ROI3](#roi3-1)
        -   [ROI4](#roi4-1)

Overview
========

Based on
[ComBatHarmonization](https://github.com/Jfortin1/ComBatHarmonization)
and [neuroHarmonize](https://github.com/rpomponio/neuroHarmonize)

ComGamHarm estimates and removes site-effect for each feature while
preserving biological variance. If Empirical Bayes is specified,
ComGamHarm can pool information across features and estimate site effect
with more accuracy. ComGamHarm utilizes
[gam](https://www.rdocumentation.org/packages/mgcv/versions/1.8-38/topics/gam)
models to allow for preservation of non-linear biological variance. This
is often important in AD biomarker data as there is a known non-linear
relationship between AD biomarkers and age.

<br>

ComGamHarm requires cross sectional data in order to estimate site
effect. If a user wants to harmonize longitudinal data, they can first
run ComGamHarm on the cross sectional data and use
<code>ApplyHarm</code> to then harmonize their longitudinal data. Both
the longitudinal data and cross sectional data must have the same batch
labels. This also works with any new data containing the same batch
labels

<br>

If you are using this package please cite the follwing papers:

|                                                                                                        |                                                                                                                                                                                          Cite                                                                                                                                                                                          |                                             Link                                             |
|--------------------------------------------------------------------------------------------------------|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|:--------------------------------------------------------------------------------------------:|
| ComBat for multi-site DTI data                                                                         |                                  Jean-Philippe Fortin, Drew Parker, Birkan Tunc, Takanori Watanabe, Mark A Elliott, Kosha Ruparel, David R Roalf, Theodore D Satterthwaite, Ruben C Gur, Raquel E Gur, Robert T Schultz, Ragini Verma, Russell T Shinohara. Harmonization Of Multi-Site Diffusion Tensor Imaging Data. NeuroImage, 161, 149-170, 2017                                  | [Here](https://www.sciencedirect.com/science/article/abs/pii/S1053811917306948?via%3Dihub#!) |
| ComBat for multi-site cortical thickness measurements                                                  | Jean-Philippe Fortin, Nicholas Cullen, Yvette I. Sheline, Warren D. Taylor, Irem Aselcioglu, Philip A. Cook, Phil Adams, Crystal Cooper, Maurizio Fava, Patrick J. McGrath, Melvin McInnis, Mary L. Phillips, Madhukar H. Trivedi, Myrna M. Weissman, Russell T. Shinohara. Harmonization of cortical thickness measurements across scanners and sites. NeuroImage, 167, 104-120, 2018 |       [Here](https://www.sciencedirect.com/science/article/abs/pii/S105381191730931X)        |
| Original ComBat paper for gene expression array                                                        |                                                                                                                 W. Evan Johnson and Cheng Li, Adjusting batch effects in microarray expression data using empirical Bayes methods. Biostatistics, 8(1):118-127, 2007.                                                                                                                  |            [Here](https://academic.oup.com/biostatistics/article/8/1/118/252073)             |
| Harmonization of large MRI datasets for the analysis of brain imaging patterns throughout the lifespan |                                                                                                   Pomponio, R., Shou, H., Davatzikos, C., et al., (2019). “Harmonization of large MRI datasets for the analysis of brain imaging patterns throughout the lifespan.” Neuroimage 208.                                                                                                    |                      [Here](https://pubmed.ncbi.nlm.nih.gov/31821869/)                       |

<br>

Installation
============

Both ***ComGamHarmFunction.R*** and ***ComGamHarmFunctionHelpers.R***
must be downloaded. The following packages must be installed:

``` r
install.packages("ggplot2")
install.packages("mgcv")
install.packages("dplyr")
install.packages("reshape2")
install.packages("matrixStats")
install.packages("ebbr")
install.packages("BiocParallel")
install.packages("furniture")
```

ComGamHarm
==========

Arguments
---------

**ComGamHarm** requires the following arguments:

-   **feature.data** <code> data.frame </code> data.frame with features
    for harmonization

-   **covar.data** <code> data.frame </code> data.frame with covariates
    for harmonization. All covariates included in
    <code>covar.data</code> will be preserved during harmonization.
    <code>covar.data</code> must contain a column <code>STUDY</code>
    corresponding to batch

-   **eb** <code> logical </code> whether or not to perform Empirical
    Bayesion estimation of site effect. <code> TRUE/FALSE </code> <code>
    defaults FALSE </code>

-   **parametric** <code>logical</code> parametric EB adjustments
    <code>defaults TRUE</code>

-   **smooth.terms** <code>character</code> vector of non-linear
    covariates for harmonization Ex: <code>c(“Age”)</code>

-   **k.val** <code>numeric</code> vector of smooth parameters for
    non-linear covariates for harmonization. Must be ordered in respect
    to <code>smooth.terms</code>

-   **verbose** <code>logical</code> print model fit progress
    <code>defaults TRUE</code>

-   **model.diagnostics** <code>logical</code> return gam fit
    diagnostics <code>defaults FALSE</code>

Example
=======

Cross-Sectional Data
--------------------

``` r
#data generated in DataSimulation.R
#not real medical data. 

#first 5 subjects in each Site
knitr::kable(full.simulated.data.cs[c(1:5, 301:305, 601:605),])
```

|     ROI1 |     ROI2 |     ROI3 |      ROI4 |      Age | Sex | ApoE | Education | STUDY |  id |
|---------:|---------:|---------:|----------:|---------:|:----|:-----|----------:|:------|----:|
| 6560.067 | 1031.122 | 36833.68 | 15209.550 | 72.19762 | 2   | 1    |  18.47753 | A     |   1 |
| 5845.162 | 1044.454 | 35690.31 | 12440.928 | 73.84911 | 1   | 1    |  17.81281 | A     |   2 |
| 4601.741 | 1009.818 | 33542.62 | 12279.165 | 82.79354 | 2   | 2    |  16.75394 | A     |   3 |
| 6633.555 | 1210.975 | 33618.14 | 11469.827 | 75.35254 | 2   | 2    |  17.23873 | A     |   4 |
| 4696.089 | 1412.610 | 33241.65 | 12656.505 | 75.64644 | 2   | 1    |  17.58859 | A     |   5 |
| 5340.044 | 1507.867 | 39532.46 | 15977.043 | 72.19762 | 2   | 1    |  18.47753 | B     | 301 |
| 7050.704 | 1477.679 | 32607.81 | 16950.297 | 73.84911 | 1   | 1    |  17.81281 | B     | 302 |
| 4768.223 | 1367.974 | 33102.17 | 15598.596 | 82.79354 | 2   | 2    |  16.75394 | B     | 303 |
| 5027.141 | 1409.213 | 34729.24 | 16162.799 | 75.35254 | 2   | 2    |  17.23873 | B     | 304 |
| 5470.245 | 1299.244 | 31844.18 | 15610.725 | 75.64644 | 2   | 1    |  17.58859 | B     | 305 |
| 4927.961 | 1112.786 | 33460.39 | 11613.298 | 72.19762 | 2   | 1    |  18.47753 | C     | 601 |
| 5268.513 | 1183.476 | 32266.08 | 13217.589 | 73.84911 | 1   | 1    |  17.81281 | C     | 602 |
| 3887.168 | 1079.390 | 29715.05 |  9979.098 | 82.79354 | 2   | 2    |  16.75394 | C     | 603 |
| 5162.835 | 1049.117 | 35851.78 | 12106.939 | 75.35254 | 2   | 2    |  17.23873 | C     | 604 |
| 4852.780 | 1208.911 | 33630.76 | 11994.121 | 75.64644 | 2   | 1    |  17.58859 | C     | 605 |

Harmonize
---------

``` r
#Feature Data
feature.data <- full.simulated.data.cs[,c("ROI1", "ROI2", "ROI3", "ROI4")]

#Covariate Data
covariate.data <- full.simulated.data.cs[,c("Age","Sex", "ApoE", "Education", "STUDY")]


#Harmonize
harmfeats <- ComGamHarm(feature.data = feature.data,
                        covar.data   = covariate.data,
                        eb           = TRUE,
                        parametric   = TRUE,
                        smooth.terms = c("Age"),
                        k.val        = c(5), 
                        verbose      = FALSE)
```

Harmonized features are contained in
<code>harmfeats\[\[“harm.results”\]\]</code>. They are in matrix format
with dimension <code>n\_features X n\_observations</code> and must be
transposed to match the original dimensons.

``` r
harmonized.features <- as.data.frame(t(harmfeats$harm.results))
```

We can then recombine our new harmonized features with the rest of our
data

``` r
harm.data <- cbind(harmonized.features, covariate.data)
```

Applying Harmonization To New Data
==================================

Using <code>ApplyHarm</code> we can apply pre-trained harmonization
models to new data under the condition the new data comes from the same
sites.

**ApplyHarm** requires the following arguments:

-   **feature.data** <code> data.frame </code> data.frame with features
    for harmonization

-   **covar.data** <code> data.frame </code> data.frame with covariates
    for harmonization. All covariates included in
    <code>covar.data</code> will be preserved during harmonization.
    <code>covar.data</code> must contain a column <code>STUDY</code>
    corresponding to batch

-   **comgam.out** <code>list</code> output from <code>ComGamHarm</code>
    function

Example
-------

``` r
harm.longitudinal <- ApplyHarm(feature.data   = feature.data.long,
                               covariate.data = covariate.data.long,
                               comgam.out     = harmfeats)

harm.longitudinal <- as.data.frame(t(harm.longitudinal))

harm.longitudinal <- cbind(harm.longitudinal, covariate.data.long)
```

Plots
=====

Plots Before Harmonization
--------------------------

### ROI1

[![roi1unharmed.png](https://i.postimg.cc/6qS9FB3G/roi1unharmed.png)](https://postimg.cc/fScn9nYw)

### ROI2

[![roi2unharmed.png](https://i.postimg.cc/q7Lp4sqr/roi2unharmed.png)](https://postimg.cc/kVB0ytWY)

### ROI3

[![roi3unharmed.png](https://i.postimg.cc/YqWkhK0b/roi3unharmed.png)](https://postimg.cc/CzS9Q2jq)

### ROI4

[![roi4unharmed.png](https://i.postimg.cc/L6LXGZt9/roi4unharmed.png)](https://postimg.cc/vc8bVT2k)

Plots After Harmonization
-------------------------

### ROI1

[![roi1harmed.png](https://i.postimg.cc/SQWbhLvv/roi1harmed.png)](https://postimg.cc/145j0FPG)

### ROI2

[![roi2harmed.png](https://i.postimg.cc/15sh5nx5/roi2harmed.png)](https://postimg.cc/DmjYxzyt)

### ROI3

[![roi3harmed.png](https://i.postimg.cc/xCMh9Ytz/roi3harmed.png)](https://postimg.cc/626Y0ssW)

### ROI4

[![roi4harmed.png](https://i.postimg.cc/m28qLmN8/roi4harmed.png)](https://postimg.cc/9rw11tk7)
