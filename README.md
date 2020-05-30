# Nanoscale Interactions from Surface Electrostatics in Python

## Preface

For most electrostatic problems, conductors can be assumed to be at an equipotential across their surface. However, at the nanoscale this is no longer true. One of the consequences is the presence of [an electrostatic force between conductors](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.90.160403) even if they are at, on average, the same electrostatic potential. This patch force presents itself as an artifact in certain precision measurements, such as measurements of the Casimir force, which was how I got involved. 

In 2014, I was working on an Casimir force measurement that wasn't working. Based on [an earlier paper](https://arxiv.org/abs/1108.1761), I thought that patch potentials might be the issue. Thus, I learned how to perform Kelvin probe force microscopy (KPFM), and [imaged the patches on a few gold plates](https://arxiv.org/abs/1409.5012). The patches were not the issue in the end, but were interesting enough to catch my attention. 

Later, I managed to measure the patch potentials on the sphere as well as the plate, by making a few changes to how I implemented KFPM. However, the way that the earlier theory was phrased made it very difficult to calculate the patch potential force from curved plates. A few changes to the formulae, were sufficient to transform the calculation from k-space to real-space, which allowed for its application to different geometries. 

## About the repository

## How to use the code

## License

LGPL 

