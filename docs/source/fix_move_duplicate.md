# Fix move duplicate


## Syntax 

```
fix_style ID group-ID move/dupl args
```

- args=list of the possible arguments

```
   C_ID = ID of the compute displace/atom
```

This fix allows to translate the duplicate atoms of the same distance crossed by the corresponding real atoms. In step3 of the cleaving methodology the duplicate atoms allow to calculate the cross-cleaving plane interactions but they need to "follow" the real atoms when they move in a new position each time step. 



```{footbibliography}

```
