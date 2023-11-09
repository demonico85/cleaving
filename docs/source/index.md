# CLEAVING

The LAMMPS package CLEAVING contains all the modified potentials and additional computes/fixes needed by the cleaving model to use Thermodynamic Integration to calculate the Surface Free Energy of a given system.

## Background

Direct determination of the IFE from the reversible work of creating an interface requires the careful design of a thermodynamic transformation path that starts from separate bulk phases equilibrated at coexistence conditions and ends with the same phases in contact across an interface. The cleaving method, originally proposed by Broughton and Gilmer {footcite:p}`Broughton1986Cleaving` for the calculation of the solid-liquid IFE in a Lennard-Jones system, constructs such a path from four basic steps:

1. introduce an external 'cleaving potential' in the first phase (*e.g.* a solid);
2. same as 1., but for the second phase (*e.g*.) a liquid;
3. rearrange the boundary conditions to merge the systems while maintaining the cleaving potentials;
4. remove the cleaving potentials.

A sketch of a cleaving path to obtain the IFE between two generic phases $\alpha$ and $\beta$ can be represented as:

![cleaving](../figs/joss.png "cleaving")

The reversible work in each step can be determined by different free energy computation methods, but those based on Thermodynamic Integration are the most straightfoward and provide accurate and direct results. With TI methods the free energy difference between two thermodynamic states connected by a transformation path is calculated by integrating the ensemble average of some configuration dependent function (*e.g.*, the potential energy) with respect to a parameter defining the path from the initial to the final thermodynamic state.

To date, the cleaving method, which owes its name to the fact that each bulk phase is "cleaved", i.e., cut at a certain position to create the new interface, has been used to calculate the solid-liquid IFE in hard and soft spheres {footcite:p}`Davidchack00prl,Davidchack05prl`, Lennard-Jones {footcite:p}`Davidchack03direct`, TIP4P model water {footcite:p}`Handel08prl,Davidchack12ice`, Embedded-Atom Model (EAM) iron {footcite:p}`Liu2013csm`, silver-ethylene glycol {footcite:p}`Qi2016jcp`, orcinol-chloroform and orcinol-nitromethane {footcite:p}`addula2020computation`, as well as for the calculation of the surface free energy of molecular crystals, such as mannitol {footcite:p}`DiPasquale2022cleaving`. However, due to the rather involved multi-step calculations required by the cleaving method, its use so far has been limited to a few research groups. Wider adoption of the cleaving approach for the calculation of solid-liquid and solid-solid IFE for different forcefields requires a more convenient simulation interface for the construction of the thermodynamic transformation path and calculation of the reversible work along this path. Here we attempt to fill this need with the CLEAVING package, which is designed to make it more straightforward for LAMMPS users to access the calculation of solid-liquid and solid-solid IFE in a wide variety of systems.

```{eval-rst}
.. toctree::
   :caption: Table of Contents
   :maxdepth: 2
   
   install.md
   lammps.md
   example.md
   utils.md
   contributing.md
```

## Bibliography

```{footbibliography}

```
