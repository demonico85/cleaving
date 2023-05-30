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

The calculation of the properties at the interface between different materials or different phases of the same material has always been an important problem due to the challenges of measuring such properties experimentally and insight into mechanisms underlying various interfacial phenomena.
While the free energy of liquid-liquid and liquid-vapour interfaces is equal to the surface tension which can be calculated from the anisotropy of the pressure tensor at the interface, the same approach cannot be used for the interfaces with a solid phase, where the presence of an interface introduces residual stress in the nearby solid, which changes the interfacial tension relative to the excess free energy [@DiPasquale2020shuttleworth].  This requires direct calculation of the free energy change during the formation of an interface from two separate bulk phases.  While free energy calculation methods have been widely available in various molecular simulation packages, their use in conjunction with the complex geometry of an interfacial system is quite challenging.  Thermodynamic Integration (TI) is one of the most commonly used approaches for computing free energy differences.  Therefore, here we present a package for the calculation of the Interfacial Free Energy (IFE) of solid-liquid and solid-solid interfaces within a widely used MD package LAMMPS (Large-Scale Atomic/Molecular Massively Parallel Simulator)[@Thompson2022].

# Statement of need

Interest in the study of the properties of interfaces with solid phases arises from the fact that many physical phenomena (e.g., freezing, nucleation, confinement) and technological processes (e.g., casting, welding, formulation) involving solid phases are governed by the structure and thermodynamics properties of the interface between solid and other phases.  Even though the definitions of various interfacial properties have been established a long time ago, starting from the works by Willard Gibbs [@GibbsCollectedWorks], their exact determination in experimensts remains very difficult due to the fact that the interface is surrounded by dense bulk phases, as well as the need for a strict control of the experimental conditions required for precise measurements (e.g., ensuring constant irregularity or porosity in all the samples on which the measurement is carried out). 

Given these difficulties, the in silico experiments offers a viable way to obtain interfacial properties directly from their definitions in terms of intermolecular interactions, giving access to complete structural and thermodynamic characterisation of the interface.  One of the fundamental thermodynamic properties of the interface is its excess free energy, which is defined as the amount of reversible work required to create a unit area of the interface between two bulk phases at coexistence conditions.  Several indirect methods for determination of this quantity exist due to its link to the rate of nucleation (via, e.g., the classical nucleation theory), or contact angles between three or more phases (via, e.g., the Young equation), or capillary fluctuations of a diffuse interface.  However, the accuracy of these methods is limited by the approximations inherent in the relevant theories.  

Direct determination of the IFE from the reversible work of creating an interface requires careful design of a thermodynamic transformation path that starts from separate bulk phases equilibrated at the coexistence conditions and ends with the same phases in contact across the interface.  The cleaving method, originally proposed by Broughton and Gilmer [@Broughton1986cleaving] for the calculation of the solid-liquid IFE in the Lennard-Jones system, constructs such a path from four basic steps: 1,2 - introduce extenernal 'cleaving potentials' in separate bulk solid and liquid systems, 3 - rearrange the boundary conditions to merge the systems while maintaining the cleaving potentials, 4 - remove the cleaving potentials.  A sketch of the cleaving path to obtain the IFE between two generic phases $\alpha$ and $\beta$ is shown in XXX.  The reversible work in each step can be determined by different free energy computation methods, but those based on Thermodynamic Integration (TI) are the most straightfoward and provide accurate and direct results. In the TI methods the free energy difference between two thermodynamic states connected by a transformation path is calculated by integrating the enesemble average of some configuration dependent function (e.g., the potential energy) with respect to a parameter defining the path from the initial to the final thermodynamic state.

To date, the cleaving method has been used to calculate the solid-liquid IFE in hard and soft spheres [@Davidchack00prl,@Davidchack05prl], Lennard-Jones [@Davidchack03direct], TIP4P model water [@Handel2008prl,@Davidchack12ice], Embedded-Atom Model (EAM) iron [@Liu13cms], silver-ethylene glycol [@Qi2016jcp], orcinol-chloroform and orcinol-nitromethane [@addula2020computation], as well as for the calculation of the surface free energy of molecular crystals, such as mannitol [@DiPasquale2022cleaving].  However, due to the rather involved multi-step calculations needed for the cleaving method, its use so far has been limited to a few research groups.  Wider adoption of the cleaving approach for the calculation of solid-liquid and solid-solid IFE for different forcefields requires a more convenient simulation interface for the construction of the thermodynamic transfomation path and calculation of the reversible work along this path.  The developed CLEAVING package is designed to make it more straightforward for LAMMPS users to access the calculation of solid-liquid and solid-solid IFE in a wide variety of systems.

# Functionality

This package includes several new additional functions written in the LAMMPS language needed to calculate the SFE through the cleaving.



# Acknowledgements

We acknowledge contributions from XXX and support from the CECAM award XXX.

# References

