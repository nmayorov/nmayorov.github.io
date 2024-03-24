---
title: Open source projects
HideMeta: true
ShowToc: false
---

Here are Python open source packages I've implemented with somewhat decent quality and usefulness.

### [pyins](https://github.com/nmayorov/pyins) 

Implementation of basic concepts and algorithms in aided Inertial Navigation Systems.
Can be useful for studying INS and estimation in general, algorithm prototyping or performance evaluation.

### [allan-variance](https://github.com/nmayorov/allan-variance) 

Tiny package for computation of Allan variance and automatic noise parameter estimation from it.
Can be used for real data analysis like evaluation of inertial sensors.

### [dyne](https://github.com/nmayorov/dyne)

A collection of estimation algorithms for discrete-time dynamic systems.
The idea is to have a playground library where estimation algorithms can be easily implemented and tested thanks to the simplicity of problem formulation and flexibility of Python.
Currently contains standard (like EKF) as well as more interesting and experimental algorithms (like nonlinear batch optimization).
