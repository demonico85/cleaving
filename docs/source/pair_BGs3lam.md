# pair_style lj/BGNlcleavs3

```
pair_style lj/BGNlcleavs3args
```

- args=list of the possible arguments

```
lj/BGNlcleavs3 args = cutoff1 cutoff2 lambda Dfac 
    cutoff1 = global internal cut-off
    cutoff2 = global external cut-off
    lambda  = global scaling of the potential
    npow    = exponent of the function lambda^npow
    Dfac    = multiplicative factor to multiply the derivative of the interaction with respect to lambda 
```


This pair style implements the Broughton and Gilmer modification to Lennard-Jones potential {footcite:t}`Broughton1983` (see [pair lj/BG]{pairBG.md}) to be used in the step3 of the cleaving algorithm. This pair style returns also an array with all the calculated interactions which are needed to calculate the work in the step3. This array can be accessed using the new compute [compute paircleav]}{compute_paircleav.md}. 


The potential implemented in this pair style is 

$$
	U(r_{ln},\lambda) =
		\begin{cases}
			4\epsilon f(\lambda)^n\left(\left(\frac{\sigma}{r_{ln}}\right)^{12} -\left(\frac{\sigma}{r_{ln}}\right)^{6}  \right)+C_1,\;\mbox{if}\; r_{ln} \leq 2.3\sigma \\
							C_2\left(\frac{\sigma}{r_{ln}}\right)^{12} + C_3\left(\frac{\sigma}{r_{ln}}\right)^{6} + C_4\left(\frac{r_{ln}}{\sigma}\right)^2 + C_5,\;\mbox{if}\; r_{ln} \leq 2.5\sigma \\
				0, 		\; r_{ln} \geq 2.5\sigma		
		\end{cases}
$$

where $r_{ln}=|\mathbf{r}_l-\mathbf{r}_n|$ for each couple of atoms $l,n$ in the system, and $C_1, C_2, C_3, C_4, C_5$ are constants we used the values reported in {footcite:t}`davidchack2003direct`.

````{note}
* $f(\lambda)$ is a continuous function $lambda$ which is equal to 0 or 1 when $lambda$ is equal to 0 or 1
````
   


The constants are hardcorded within the pair style and they not need to be defined.

```{footbibliography}

```
