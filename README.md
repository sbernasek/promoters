Repository overview
-------------------

This repository and its dependencies contain all of the code necessary to reproduce our study of the relationship between Gene Regulation and Metabolism. The code offers two primary functions:

  1. Running the simulations described in our manuscript. In general, these simulations probe pulse response sensitivity to gene regulatory network (GRN) perturbations under a variety of different metabolic conditions. *promoters* provides the tools necessary to set up each of these simulations and analyze the resultant dynamics. The simulations themselves are performed by our stochastic simulation package, [GeneSSA](https://github.com/sebastianbernasek/genessa).

  2. Quantifying Yan protein expression dynamics in the Drosophila eye. Yan level measurements were extracted from confocal microscopy data using [FlyEye Silhouette](http://www.silhouette.amaral.northwestern.edu/), our macOS platform for eye cell segmentation and annotation. The annotated `.silhouette` files are available in our [data repository](https://arch.library.northwestern.edu/concern/generic_works/n296wz31t?locale=en). We analyzed these measurements using [FlyEye Analysis](https://github.com/sebastianbernasek/flyeye), our pipeline for analyzing FlyEye Silhouette data.


Note on reproducibility
-----------------------

This repository contains all of the code used to generate the content in our manuscript. Please note that the vast majority of our figures are based on large-scale simulations that were executed on a high performance computing cluster. These simulations are therefore not well suited for reproduction on a personal computer. We have provided the scripts necessary to reproduce all of our results, but we caution that their direct execution would take a very long time.

As an alternative, we have also provided the output from all of our [completed simulations](https://doi.org/10.21985/n2-j361-8e86) along with a series of [Jupyter notebooks](https://github.com/sebastianbernasek/promoters/tree/master/notebooks) that walk the user through the steps necessary to analyze these results and reproduce each of our figures. The notebooks also provide users with an opportunity to set up and execute each type of individual simulation that appears in our manuscript. We leave it to the user to design a means to execute these simulations en masse.


Supporting data
---------------

[data.zip](https://doi.org/10.21985/n2-j361-8e86) (~15 MB) contains the completed output from each of our simulations.


Download the file, then unzip its contents to a common directory. In order to successfully run the provided Jupyter notebooks you will need to point the ``../data`` filepath toward this directory.


Installation
============

Before attempting to install *promoters*, we suggest creating a clean virtual environment and installing all necessary dependencies first. If you intend to reproduce our numerical simulations, you will first need to compile and install our stochastic simulation package, [GeneSSA](https://github.com/sebastianbernasek/genessa).


System requirements
-------------------

 - Python 3.6+
 - [Scipy](https://www.scipy.org/)
 - [Pandas](https://pandas.pydata.org/)
 - [Matplotlib](https://matplotlib.org/)
 - [GeneSSA](https://github.com/sebastianbernasek/genessa)
 - [FlyEye Analysis](https://github.com/sebastianbernasek/flyeye) 


Install promoters
-----------------

Download the [latest distribution](https://github.com/sebastianbernasek/promoters/archive/v1.0.tar.gz).

The simplest method is to install it via ``pip``:

    pip install promoters-1.0.tar.gz



Package contents
================

The ``promoters`` package consists of a set of python modules, scripts, and notebooks that walk the user through reproducing all of our simulation results. Their contents are as follows.


Modules
-------

The ``promoters`` package consists of several modules:

  * ``promoters.models`` provides templates for each type of GRN discussed in our manuscript.

  * ``promoters.simulation`` provides methods for performing pulse response simulations.

  * ``promoters.analysis`` provides methods for comparing expression dynamics between simulations, e.g. evaluating "error frequency".

  * ``promoters.sweep`` provides methods for constructing and executing a parameter sweep of each type of model.

  * ``promoters.figures`` provides templates for each type of figure that appears in our manuscript.


Scripts
-------

The ``promoters`` package contains several python scripts in ``promoters/scripts``. Those that may prove helpful include:

  * ``build_sweep.py`` initializes a parameter sweep of a specified model.

  * ``run_simulation.py`` runs an individual ``ConditionSimulation``.

  * ``run_batch.py`` runs a batch of ``ConditionSimulation`` instances.


Jupyter notebooks
-----------------

  * [Example Simulation (Fulll vs Partial Synthesis).ipynb](https://github.com/sbernasek/promoters/blob/master/notebooks/Example%20Simulation%20(Full%20vs%20Partial%20Synthesis).ipynb) walks the user through conducting and analyzing simulations in which auxiliary promoters are lost.

  * [Example Simulation (Fulll vs Partial Repression).ipynb](https://github.com/sbernasek/promoters/blob/master/notebooks/Example%20Simulation%20(Full%20vs%20Partial%20Repression).ipynb) walks the user through conducting and analyzing simulations in which auxiliary repressors are lost.

  * [Parameter Sweeps.ipynb](https://github.com/sbernasek/promoters/blob/master/notebooks/Parameter%20Sweeps.ipynb) walks the user through conducting each of our parameter sweeps and visualizing the results.

  * [Simulating Promoter Loss at Multiple Levels of Output Regulation.ipynb](https://github.com/sbernasek/promoters/blob/master/notebooks/Simulating%20Promoter%20Loss%20at%20Multiple%20Levels%20of%20Output%20Regulation.ipynb) walks the user through conducting and analyzing simulations in which one of a pair of promoters is removed.

  * [Simulating Repressor Loss at Multiple Levels of Output Regulation.ipynb](https://github.com/sbernasek/promoters/blob/master/notebooks/Simulating%20Repressor%20Loss%20at%20Multiple%20Levels%20of%20Output%20Regulation.ipynb) walks the user through conducting and analyzing simulations in which one of a pair of repressors is removed


Additional resources
====================

For examples detailing the usage of our stochastic simulation software, please see [GeneSSA](https://github.com/sebastianbernasek/genessa). For examples demonstrating the analysis of protein expression in the Drosophila eye, please see [FlyEye Analysis](https://github.com/sebastianbernasek/flyeye).
