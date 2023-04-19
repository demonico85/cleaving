# pair_style lj/BGcleavs3

## Syntax

```
pair_style lj/BGcleavs3 args
```

- args=list of the possible arguments

```
lj/BGNlcleavs3 args = cutoff1 cutoff2 lambda Dfac 
    cutoff1 = global internal cut-off
    cutoff2 = global external cut-off
    lambda  = global scaling of the potential
    Dfac    = multiplicative factor for the derivative of the interaction with respect to lambda 
```

`pair_coeff` accept the same arguments of the pair_style. However, keep in mind that the order is not the same

```
pair_coeff a b lambda Dfac cutoff1 cutoff2
```

where

```
    a       = atom of type a [mandatory]
    b       = atom of type b [mandatory]
    lambda  = global scaling of the potential
    Dfac    = multiplicative factor for the derivative of the interaction with respect to lambda     
    cutoff1 = global internal cut-off
    cutoff2 = global external cut-off
```

the arguments in the `pair_coeff` are all optional except for the atom types. Note that the order is important. If we want to specify a new lambda for a different pair of types we can just write

```
pair_coeff 1 1 1.0
```

However, in order to specify a new (e.g., external) cutoff we need to rewrite all the preceeding arguments:

```
pair_coeff 1 1 lambda Dfac cutoff1 cutoff2
```

## Examples

```
pair_style lj/BGNlcleavs3 2.3 2.5 
pair_coeff 1 1 1.0
pair_coeff 1 2 1.0
pair_coeff 1 3 -1.0
```

## Description

This pair style implements the Broughton and Gilmer modification to Lennard-Jones potential {footcite:t}`Broughton1983` (see [pair lj/BG]{pairBG.md}) to be used in the step3 of the cleaving algorithm. This pair style returns also an array with all the calculated interactions which are needed to calculate the work in the step3. This array can be accessed using the new compute [compute paircleav]}{compute_pcleav.md}. 


The potential implemented in this pair style is 

$$
	U(r_{ln},\lambda) =
		\begin{cases}
			\lambda\left[ 4\epsilon \left(\left(\frac{\sigma}{r_{ln}}\right)^{12} -\left(\frac{\sigma}{r_{ln}}\right)^{6}  \right)+C_1 \right],\;\mbox{if}\; r_{ln} \leq 2.3\sigma \\
			\lambda\left[ C_2\left(\frac{\sigma}{r_{ln}}\right)^{12} + C_3\left(\frac{\sigma}{r_{ln}}\right)^{6} + C_4\left(\frac{r_{ln}}{\sigma}\right)^2 + C_5\right],\;\mbox{if}\; r_{ln} \leq 2.5\sigma \\
				0, 		\; r_{ln} \geq 2.5\sigma		
		\end{cases}
$$

where $r_{ln}=|\mathbf{r}_l-\mathbf{r}_n|$ for each couple of atoms $l,n$ in the system, and $C_1, C_2, C_3, C_4, C_5$ are constants we used the values reported in {footcite:t}`davidchack2003direct`.
The constants are hardcorded within the pair style and they not need to be defined.

````{note}
* $lambda$ must vary in the (closed) interval [0,1]
* There is no internal check that $\lambda$ is changing consistently (i.e., always decreasing or increasing)
````
   
The work performed on the system to switch from $\lambda=1$ to $\lambda=0$ is then given by

$$
	\frac{\partial U(r_{ln},\lambda)}{\partial \lambda} =
		\begin{cases}
			\left[ 4\epsilon \left(\left(\frac{\sigma}{r_{ln}}\right)^{12} -\left(\frac{\sigma}{r_{ln}}\right)^{6}  \right)+C_1 \right],\;\mbox{if}\; r_{ln} \leq 2.3\sigma \\
			\left[ C_2\left(\frac{\sigma}{r_{ln}}\right)^{12} + C_3\left(\frac{\sigma}{r_{ln}}\right)^{6} + C_4\left(\frac{r_{ln}}{\sigma}\right)^2 + C_5\right],\;\mbox{if}\; r_{ln} \leq 2.5\sigma \\
				0, 		\; r_{ln} \geq 2.5\sigma		
		\end{cases}
$$

In this case the the expression of the potential is multiplied by the derivatibe of $\lambda$ with respect itself which is equal to 1 (i.e., $\partial \lambda / \partial \lambda=1$), and therefore the value of Dfrac must be set equal to 1.

When there is a contemporary switching (i.e., some interactions are multiplied by $\lambda$ and others by $1-lambda$) then the derivative  $\partial (1-\lambda) / \partial \lambda=-1$). In this case we need to set Dfrac=-1.

```
variable lambda file lambda.dat
variable minl   equal 1-lambda

pair_style lj/BGNlcleavs3  cutoff1  cutoff2 lambda 1.0
pair_coeff 1 2 minl -1.0

```



```{footbibliography}

```
