# Nanoscale Interactions from Surface Electrostatics in Python

### (pronounced "nice pie")


## Preface

For most electrostatic problems, conductors can be assumed to be at an equipotential across their surface. However, at the nanoscale this is no longer true. One of the consequences is the presence of [an electrostatic force between conductors](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.90.160403) even if they are at, on average, the same electrostatic potential. This patch force presents itself as an artifact in certain precision measurements, such as measurements of the Casimir force, which was how I got involved. 

In 2014, I was working on an Casimir force measurement that wasn't working. Based on [an earlier paper](https://arxiv.org/abs/1108.1761), I thought that patch potentials might be the issue. Thus, I learned how to perform Kelvin probe force microscopy (KPFM), and [imaged the patches on a few gold plates](https://arxiv.org/abs/1409.5012). The patches were not the issue [in the end](https://drum.lib.umd.edu/handle/1903/20311), but were interesting enough to catch my attention. 

Later, I managed to measure the patch potentials on the sphere as well as the plate, by making a few changes to how I implemented KFPM. However, the way that the earlier theory was phrased made it very difficult to calculate the patch potential force from curved plates. A few changes to the formulae, were sufficient to transform the calculation from *k*-space to real-space, which allowed for its application to different geometries. 

Applying the formulae was not straight-forward for a variety of reasons: the kernels involved with the convolution changed with separation (*h*), they were large (at least 81x81), and they required a long time to calculate (about half an hour on a laptop). I wrote the code in this repository as a way of organizing the formulae, and solving them in a straight-forward way. The unfortunate part of this is that, in order to save time, many intermediate steps must be saved, which takes 100s of MB of storage.

Although the code in this repository is quite old, dating back to 2017-2018, and does not adhere at all to the PEP8 standard that I now tend to follow, I think that it is still worth posting. First of all, it was the basis for [a paper that was recently published](https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.2.023355). Second of all, I hope that it encourages others currently involved in force measurements to calculate the effect of patch potentials in their own setups. My hope is that one of you, dear readers, decides to take up development of this code for your own projects. 

## About the repository

The repository is divided into several .py files:

CypherKPFM.py - code to import raw data from the Cypher AFM, when a particular setup is used for KPFM. It is included so that the example data can be imported. 

filters3D.py - functions to create the filtered versions of the patch potentials

force3D.py - functions and classes to calculate the force in 3D geometries once the filtered patch potentials are created

patches2D.py - contains the code to calculate the patch force in the original *k*-space formulation

patchpotentials.py - includes functions for generating the sphere and calculating the force from pre-filtered images for a generic height map. The functions is this are mostly more simple/less general than those found in the other .py files. Also, they are used in the other functions. 

patchPotentials_testing.py - contains a function for creating simulated data for testing. I didn't really use it

preprocessing.py - a script for automatically calculating the force from a particular set of images. I recommend not using this, because it is difficult to interpret

SavingDict.py - used for saving files

## How to use the code

Although the documentation leaves something to be desired (sorry about that), I think that the tutorial (the .ipynb file) is fairly straight-forward to follow. 
