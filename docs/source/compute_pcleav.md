# Compute Pcleav

```
compute_style cleavpairs args
```

- args=list of the possible arguments

```
    cleavpairs args = < pair_style name> norm N <direction>
    < pair_style name > = name of the pair_style from which interactions are collected
    norm N = number N of new surfaces created to normalize the interactions
    <direction> = possible values x, y, z. The plane perpendicular to the selected direction is used
```

This compute collects the scaled interactions which are used to calculate the work to create a new interface. Note that it cannot be used for any pair_style, since the calculation of the scaled interactions for the cleaving must be obtained within the pair_style itself. For this reason, modified pair_styles are available in the package (see ADD).

The norm is the surface of the new interface(s) created. Note that at the moment the normalization must be specified with at least one interface.

The direction specify the normal to the newly created surface. If e.g., z is specified, then the new interface $S$ will be equal to  $S = x_{edge} * y_{edge}$.


Note that at the moment this works only for non-triclinic box.
