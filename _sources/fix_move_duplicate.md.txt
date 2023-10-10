# fix move/dupl


## Syntax 

```
fix ID group-ID move/dupl C_ID
```

* `C_ID` = ID of the compute displace/atom

## Description

This fix allows to translate the duplicate atoms of the same distance crossed by the corresponding real atoms. In step3 of the cleaving methodology the duplicate atoms allow to calculate the cross-cleaving plane interactions but they need to "follow" the real atoms when they move in a new position each time step. 


In the step3 the system (solid+liquid) is duplicated and `real` and `duplicate` atom types are defined (see [example LJ-SL](example_SL_walls.md) for more details) as showin in the next Figure.


![definition](../figs/dupl1.png "Definition of duplicate atoms")


The duplicate atoms are used to keep track of the interactions crossing the cleaving planes, which are the ones to be calculated in order to find the total work in step3. 
However, the equation of motion must be integrated _only_ for the real atoms. After the integration steo the _real_ atoms move by a distance $\Delta \mathbf{x}$ (see next Figure).

![move](../figs/dupl2.png "Move the real atoms")


This fix allows the duplicated atoms to follow the movement of the corresponding real atoms. By specifying this fix we are moving the _duplicated_ atoms by the same quantity calculated for the real atoms as shown in the next Figure.

![move](../figs/dupl3.png "Move the real atoms")


```{footbibliography}

```
