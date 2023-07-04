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

Calculation of the properties at the interface between different materials or different phases of the same
material has always been an important problem due to the challenges of measuring such properties
experimentally and insight into mechanisms underlying various interfacial phenomena. While the free energy of
liquid-liquid and liquid-gas interfaces is equal to the surface tension which can be calculated from the
anisotropy of the pressure tensor at the interface, the same approach cannot be used for the interfaces with
a solid phase, where the presence of an interface introduces residual stress in the nearby solid, which
changes the interfacial tension relative to the excess free energy [@DiPasquale2020shuttleworth].  This
requires direct calculation of the free energy change during the formation of an interface from two separate
bulk phases.  While free energy calculation methods have been widely available in various molecular
simulation packages, their use in conjunction with the complex geometry of an interfacial system is quite
challenging.  Here, we present a package for the calculation of the Interfacial Free Energy
(IFE) using Thermodynamic Integration (TI) of solid-liquid and solid-solid interfaces within a widely used MD package LAMMPS (Large-Scale
Atomic/Molecular Massively Parallel Simulator)[@Thompson2022].

# Statement of need

Interest in the study of the properties of interfaces with solid phases arises from the fact that many
physical phenomena (*e.g.*, freezing, nucleation, confinement) and technological processes (*e.g.*, casting,
welding, formulation) involving solid phases are governed by the structure and thermodynamics properties of
the interface between a solid and other phases.  Even though the definitions of various interfacial
properties have been established a long time ago, starting with the works by Willard Gibbs
[@GibbsCollectedWorks], their exact determination in experimensts remains very difficult due to the fact that
the interface is surrounded by dense bulk phases, as well as the need for a strict control of the
experimental conditions required for precise measurements (*e.g.*, ensuring constant irregularity or porosity
in all the samples on which the measurement is carried out).

Given these difficulties, *in silico* experiments offer a viable way to obtain interfacial properties
directly from their definitions in terms of intermolecular interactions, giving access to complete structural
and thermodynamic characterisation of the interface.  One of the fundamental thermodynamic properties of the
interface is its excess free energy, which is defined as the amount of reversible work required to create a
unit area of the interface between two bulk phases at coexistence conditions.  Several indirect methods for
the determination of this quantity exist, thanks to its link with the rate of nucleation (via, *e.g.*, the
classical nucleation theory), with contact angles between three or more phases (via, *e.g.*, the Young
equation), or with capillary fluctuations of a diffuse interface. However, the accuracy of these methods is
limited by the approximations inherent to the relevant theories.

Direct determination of the IFE from the reversible work of creating an interface requires careful design of
a thermodynamic transformation path that starts from separate bulk phases equilibrated at the coexistence
conditions and ends with the same phases in contact across an interface. The cleaving method, originally
proposed by Broughton and Gilmer [@Broughton1986cleaving] for the calculation of the solid-liquid IFE in a
Lennard-Jones system, constructs such a path from four basic steps: 1,2 - introduce external 'cleaving
potentials' in separate bulk solid and liquid systems, 3 - rearrange the boundary conditions to merge the
systems while maintaining the cleaving potentials, 4 - remove the cleaving potentials.

![The thermodynamic path used to compute the interfacial free energy between phases $\alpha$ and $\beta$ with
the cleaving method\label{fig:cleaving}](Fig/joss.png){ width=30% }

A sketch of the cleaving path to obtain the IFE between two generic phases $\alpha$ and $\beta$ is shown in
Fig.\autoref{fig:cleaving}. The reversible work in each step can be determined by different free energy
computation methods, but those based on Thermodynamic Integration are the most straightfoward and provide
accurate and direct results. With TI methods the free energy difference between two thermodynamic states
connected by a transformation path is calculated by integrating the ensemble average of some configuration
dependent function (*e.g.*, the potential energy) with respect to a parameter defining the path from the
initial to the final thermodynamic state.

To date, the cleaving method has been used to calculate the solid-liquid IFE in hard and soft spheres
[@Davidchack00prl,@Davidchack05prl], Lennard-Jones [@Davidchack03direct], TIP4P model water
[@Handel2008prl,@Davidchack12ice], Embedded-Atom Model (EAM) iron [@Liu13cms], silver-ethylene glycol
[@Qi2016jcp], orcinol-chloroform and orcinol-nitromethane [@addula2020computation], as well as for the
calculation of the surface free energy of molecular crystals, such as mannitol [@DiPasquale2022cleaving]. 
However, due to the rather involved multi-step calculations needed for the cleaving method, its use so far
has been limited to a few research groups.  Wider adoption of the cleaving approach for the calculation of
solid-liquid and solid-solid IFE for different forcefields requires a more convenient simulation interface
for the construction of the thermodynamic transfomation path and calculation of the reversible work along
this path. The developed CLEAVING package is designed to make it more straightforward for LAMMPS users to
access the calculation of solid-liquid and solid-solid IFE in a wide variety of systems.

# Functionality

Here, we are presenting a package designed to be integrated in LAMMPS for calculation of interface free
energy using Thermodynamic Integration. Therefore, this package includes several new additional functions
written in the style suited to be directly patched into LAMMPS. The package includes both a version of the
interactions between atoms needed to run the MD simulation with the cleaving model and some auxiliary
functionalities (in LAMMPS terms, new "fix" and "compute") needed to perform the operations required in the
different step of the calculation.

The version of the code we present here is hosted on GitHub, where a detailed instruction manual completed by
few step-by-step examples reproducing some of the published results on the topic. The packages includes new
definitions of pair potentials not already included in LAMMPS, such as the Broughton-Gilmer modified
Lennard-Jones potential [@Broughton1983] and modifications of some of the already existing pair_styles in
LAMMPS. In addition, we included a modified version of the lj/cut and coul/dsf potentials. These
modifications allow to run the step of the cleaving method (which we usually refers as step3) where the
interactions between the phases are switched-on and off. In this case, the work performed during this
particular step depend on the interactions between atoms involved in the switching and we need to be able to
keep track of them. These new pair_style includes several version to use different scaling (i.e., different
power of the coupling parameter $\lambda$, see e.g., [@DiPasquale2022cleaving]). The processing of the interactions to generate the work in step3 is then left to the new compute 'cleavpairs'
to generate the output file needed for the calculation of the IFE. 

The cleaving method owes its name to the fact that each bulk phase is "cleaved", i.e., cut at a certain position to create the new interface. The cut
is modelled with an external potential (where we define "external potential" a potential depending on the
absolute position of atoms within the simulation box). In the cleaving package there are two external
potentials available the "walls" and the "wells". These two external potential are defined as two additional fix styles:
"wallforce" and "wellPforce". We refer to the documentation for a detailed explanation of their use. A new
fix style "move/dupl" used in the step3 to switch the interactions between two phases in contact is also
included.

The new functionalities presented here are extremely general and can be applied to a variety of problem, a
fact that makes the package presented extremely flexibile. This package answer to one of the problem in the
community, the presence of a very few systematic and well maintained codes for these kind of calculations
which usually relies on in-house codes and extensions. 


# Acknowledgements

We acknowledge the support from the CECAM and CCP5 through the CECAM/CCP5 sandpits 2022.

# References

