# compute cleavpairs

## Syntax

```
compute ID group_ID cleavpairs <pair_style name> norm N <direction>
```

* `<pair_style name>` = name of the pair_style from which interactions are collected
* `norm` = include this keyword if you want to normalize the results in units of surface created 
* `N` = number N of new surfaces created to normalize the interactions (requires norm keyword)
* `<direction>` = possible values x, y, z. The plane perpendicular to the selected direction is used (requires norm keyword)

## Description

This compute collects the scaled interactions which are used to calculate the work to create a new interface. Note that it cannot be used for any pair_style, since the calculation of the scaled interactions for the cleaving must be obtained within the pair_style itself. For this reason, modified pair_styles are available in the package (see the documented pair styles).

The `norm` keyword normalize the interaction for units of interface created (i.e., the result will be in units of [Energy] / [lenght]$^2$). The keyword must be followed by:

1. a positive integer $N$ specifying the number of interfaces created, which is used to calculate the total area created to normalize the work
2. a direction (x,y,z). The direction specify the normal to the newly created surface. If e.g., z is specified, then the new interface $S$ will be equal to  $S = x_{\rm edge} \cdot y_{\rm edge}$. 

The total area used for the normalization will therefore be $S_{\rm TOT}=S \cdot N=x_{\rm edge} \cdot y_{\rm edge} \cdot N$, where $x_{\rm edge}$ and $y_{\rm edge}$ are the size of the box in the x and y directions.

Note: that at the moment this works only for non-triclinic box.

Note: If the keyword norm is not specified, then the compute will only report the energy taken from the specified pair_style (no normalization). 
