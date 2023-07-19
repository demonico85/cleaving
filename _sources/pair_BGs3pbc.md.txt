# pair_style lj/BGcleavpbc

## Syntax

```text
pair_style lj/BGcleavpbc cutoff1 cutoff2 dir
```

* `cutoff1` = global internal cut-off
* `cutoff2` = global external cut-off
* `dir`     = direction normal to the cleaving plane

`pair_coeff` accepts the same arguments as the pair_style except for direction

```text
pair_coeff a b cutoff1 cutoff2
```

where

```text
a       = atom of type a [mandatory]
b       = atom of type b [mandatory]  
cutoff1 = global internal cut-off
cutoff2 = global external cut-off
```

the arguments in the `pair_coeff` are all optional except for the atom types. Note that the order is important the internal cut-off must be declared before the external one

## Examples

```
pair_style lj/BGcleavs3 2.3 2.5 
pair_coeff 1 1 
pair_coeff 1 2 
```

## Description

This pair style implements the Broughton and Gilmer modification to Lennard-Jones potential {footcite:t}`Broughton1983` (see [pair lj/BG]{pairBG.md}) to be used in the step3 of the cleaving algorithm. 
The switching off of the interaction in this case is not obtained by varying a $\lambda$ parameter in the interval $[0,1]$ but by moving away the images of the box from the original one (in a specified direction) until the box itself is not interacting with the images anymore, see next figure.

![definition](../figs/dupl1.png "Interactions switching-off process")


This pair style returns also an array with all the calculated interactions which are needed to calculate the work in the step3. This array can be accessed using the new compute [compute paircleav]}{compute_pcleav.md}. 


The potential implemented in this pair style is the same described in [pair lj/BG]{pairBG.md}. The only difference is the fact that this pair style keeps track of the interactions across the periodic images which are the ones needed to calculate the work in this step. The movement of the periodic images can be obtained using the LAMMPS command
[change_box](https://docs.lammps.org/change_box.html). An example of the application of this pair_style is reported in [example SV-wells]{../example_SV_wells.md}.

We refer to {footcite:t}`DiPasquale2020` for the details of the calculation of the work for this pair style.


```{footbibliography}

```
