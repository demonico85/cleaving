---
title: 'Cleaving: a LAMMPS package to compute surface free energies'
tags:
  - C++
  - LAMMPS
  - Molecular simulations
  - Phase Transitions
authors:
  - name: Nicodemo Di Pasquale^[Corresponding author]
    orcid:  0000-0001-5676-8527
    affiliation: 1
  - name: Ruslan Davidchack
    orcid: 0000-0001-9418-5322
    affiliation: 4
  - name: Lorenzo Rovigatti
    orcid: 0000-0001-5017-2829
    affiliation: "2,3"
affiliations:
 - name: Department of Chemical Engineering, Brunel University London, United Kingdom
   index: 1
 - name: Department of Physics, Sapienza University of Rome, Italy
   index: 2
 - name: CNR-ISC UoS Sapienza, Rome, Italy
   index: 3
 - name: Department of Mathematics, University of Leicester, United Kingdom
   index: 4
date: 06 Apr 2023
bibliography: paper.bib
---

# Summary

The calculation of the properties at the interface of different phases has always been an important problem, with applications ranging from the understanding of the behaviour of materials at the interface to technological applications. The main driving in the careful study of interfaces is because important phenomena happen at the interface, or depend on the interface properties, such as the nucleation of the relative stability of the different directions of a crystal. The thermodynamic definition of the Surface Free Energy (SFE) requires the calculation of a free energy, which cannot be directly calculated with Molecular Dynamics (MD) simulations. Different models are available to access the free energy type of quantity of a system and one of the most used is the Thermodynamic Integration (TI). Here, we present a package for the calculation of SFE in MD that interfaces with widely used MD package LAMMPS (Large-Scale Atomic/Molecular Massively Parallel Simulator)[@Thompson2022].

# Statement of need

Interest in the modelling of the properties of interfaces arises from the fact that several phenomena (freezing, nucleation, confinement) and technoloical processes (casting, wielding, formulation) involving solid phases require the detailed knowledge of the structure and thermodynamics properties of the interface between solid and other phases. While the problem of the determination of such surface properties of materials was formalized by Gibbs for solid and liquid interface, their exact determination in experimenst is still very difficult due to the strict control necessary on the experimental conditions required for a precise measurements (e.g., ensuring a constant irregularity or porosity in all the sample on which the measure is carried on). 

Given these difficulties, the in silico experiments offers a viable way to derive interfacial properties directly from the thermodynamic definitions, giving access to a complete surface characterization. Among the different methods to calculate the SFE in MD simulations, those based on Thermodynamic Integration are gaining lots of attention for their ability to provide accurate and direct results. In TI methods a continuous thermodynamic paths is defined between two points in the space of thermodynamics parameters of the system. The free energy difference between these points is determined by the reversible work needed to transform the system from the initial point to the final point and is calculated by integrating the enesemble average of some configuration dependent function (e.g., the potential energy) with respect a parameter defining the path from an initial to a final thermodynamic point.
 
Among different models available to calculate the SFE, this package presents the cleaving method for solid-liquid interfaces. Cleaving was popularized in [@davidchack2003direct] and used for the calculations of the SFE in several systems, including molecular crystals (mannitol) [@DiPasquale2022cleaving] and Lennard-Jones liquid-solid systems [@DiPasquale2020shuttleworth]. A sketch of the cleaving path to obtain the SFE between two generic phases $\alpha$ and $\beta$ is shown in  XXX.

# Functionality

This package includes several new additional functions written in the LAMMPS language needed to calculate the SFE through the cleaving.



# Acknowledgements

We acknowledge contributions from XXX and support from the CECAM award XXX.

# References

