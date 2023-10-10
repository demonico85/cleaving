# compute displace/atom_cleav


## Syntax

```
compute ID group_ID displace/atom_cleav 
```

This compute does not need any arguments. It calculates the $\Delta x$, $\Delta y$, $\Delta z$ of all the atoms where $\Delta$ is the difference in the position of the atom between the current time-step (after the integration of the equation of motion), and the previous time-step.
The compute is needed in the step3 of the cleaving model to determine the evolution of the duplicate atoms (see [fix move/dupl](fix_move_duplicate.md))
