# pair_style lj/cleavcutNlint



- args=list of the possible arguments
```
<name pair style> args =  cut-off lambda i_IDflag orientation alpha
    cutoff        = global internal cut-off
    lambda        = global scaling of the potential
    i_IDflag      = new (integer) property to be added to the atoms
    orientation   = direction perpendicular to the cleaving plane
    N             = power of the polynomial for $\lambda$
    rspace        = short range cut-off (distance units)
```

`pair_coeff` accepts only the following arguments

```
pair_coeff a b lambda 
```

where

```
    a       = atom of type a [mandatory]
    b       = atom of type b [mandatory]
    epsilon = energy constant (energy units)
    sigma   = VdW radius (distance units)
    lambda  = scaling of the potential
    cut-off  = long range cut-off for LJ potential
```


This pair style is derived from the `pair/ljcut` in LAMMPS and we refer to the [LAMMPS documentation](https://docs.lammps.org/pair_lj.html) for the
description of its main features. Here, only the modifications needed for the cleaving calculations
will be considered. 


This pair style modifies the short range portion of the Lennard-Jones potential in order to avoid the "Lennard-Jones catastrophe". The potential is modified by smoothly interpolating it to zero for distances smaller than $r_s$ as:
$$
	U(r_{ln}) =
		\begin{cases}
			4\epsilon\left(\left(\frac{\sigma}{r_{ln}}\right)^{12} -\left(\frac{\sigma}{r_{ln}}\right)^{6}  \right)+C_1,\;\mbox{if}\; r_{ln} \leq r_{s} \\
			4\epsilon\left(\left(\frac{\sigma}{r_{ln}}\right)^{12} -\left(\frac{\sigma}{r_{ln}}\right)^{6}  \right)+C_1,\;\mbox{if}\; r_{ln} \leq r_{cut}	
		\end{cases}
$$

where $r_{ln}=|\mathbf{r}_l-\mathbf{r}_n|$ for each couple of atoms $l,n$ in the system and $r_c$ is the cut-off for the potential.


The description of the cleaving methodology using this pair style is analogous to the one described for the pair style [coul/dsf[cleav]](./pair_coul_cleav.md) to which we refer.


```{footbibliography}

```
