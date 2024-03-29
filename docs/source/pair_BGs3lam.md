# pair_style lj/BGcleavs3

## Syntax

```text
pair_style lj/BGcleavs3 cutoff1 cutoff2 lambda Dfac
```

* `cutoff1` = global internal cut-off
* `cutoff2` = global external cut-off
* `lambda`  = global scaling of the potential
* `Dfac`    = multiplicative factor for the derivative of the interaction with respect to lambda 

`pair_coeff` accepts the same arguments as the pair_style. However, keep in mind that parameters in pair coeff have a specific order:

```text
pair_coeff a b lambda Dfac cutoff1 cutoff2
```

where

```text
a       = atom of type a [mandatory]
b       = atom of type b [mandatory]
lambda  = scaling of the potential
Dfac    = multiplicative factor for the derivative of the interaction with respect to lambda     
cutoff1 = global internal cut-off
cutoff2 = global external cut-off
```

The arguments in the `pair_coeff` are all optional except for the atom types. However, if we need to specify a particular argument, we need to specify also all the preceeding ones. 
If we want to specify a new lambda for a different pair of types we can just write

```
pair_coeff 1 1 1.0
```

Instead, in order to specify a new (e.g., external) cutoff we need to rewrite all the preceeding arguments:

```
pair_coeff 1 1 lambda Dfac cutoff1 cutoff2
```

## Examples

```
pair_style lj/BGcleavs3 2.3 2.5 
pair_coeff 1 1 1.0
pair_coeff 1 2 1.0
pair_coeff 1 3 -1.0
```

## Description

This pair style implements the {footcite:t}`Broughton1983` modification to Lennard-Jones potential  (see [pair lj/BG](pairBG.md)) to be used in the step3 of the cleaving algorithm. This pair style returns also an array with all the calculated interactions which are needed to calculate the work in the step3. This array can be accessed using the new compute [compute paircleav](compute_pcleav.md). 


The potential implemented in this pair style is 

$$
	U(r_{ln},\lambda) =
		\begin{cases}
			\lambda\left[ 4\epsilon \left(\left(\frac{\sigma}{r_{ln}}\right)^{12} -\left(\frac{\sigma}{r_{ln}}\right)^{6}  \right)+C_1 \right],\;\mbox{if}\; r_{ln} \leq 2.3\sigma \\
			\lambda\left[ C_2\left(\frac{\sigma}{r_{ln}}\right)^{12} + C_3\left(\frac{\sigma}{r_{ln}}\right)^{6} + C_4\left(\frac{r_{ln}}{\sigma}\right)^2 + C_5\right],\;\mbox{if}\; r_{ln} \leq 2.5\sigma \\
				0, 		\; r_{ln} \geq 2.5\sigma		
		\end{cases}
$$

where $r_{ln}=|\mathbf{r}_l-\mathbf{r}_n|$ for each couple of atoms $l,n$ in the system, and $C_1, C_2, C_3, C_4, C_5$ are constants hardcorded within the pair style, therefore they not need to be defined. For C_1, C_2, C_3, C_4, C_5$ we used the values reported in {footcite:t}`davidchack2003direct`. 

````{note}
* $\lambda$ must vary in the (closed) interval [0,1]
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

In this case, the the expression of the potential is multiplied by the derivative of $\lambda$ with respect itself, which is equal to 1 (i.e., $\partial \lambda / \partial \lambda=1$). Therefore, the value of Dfrac must be set equal to 1.

Using this style we can represent different switches between different states of the system:

* We want to describe the switching between two different states of the system, e.g., $a$ and $a+b$. 
Let us assume the total Hamiltonian of the system, $H^{T}(\lambda)$, is defined as 

$$
	H^{T}(\lambda)= \lambda H_{a} + H_{b} 
$$  

for $\lambda \in [0,1]$, where we dropped the dependence on $\mathbf{p}$ and $\mathbf{q}$ for simplicity.

Let's consider this example:

```
variable lambda file lambda.dat

pair_style lj/BGNlcleavs3  2.3 2.5 1.0 1.0
pair_coeff 1 1 1.0 1.0
pair_coeff 2 2 1.0 1.0

pair_coeff 1 2 ${lambda} 1.0
```

The "self" interactions between atoms of the same type remains, but the mixed interactions between atoms of different types are `switched-on` (if $\lambda$ goes from 0 to 1) or `switched-off` (if $\lambda$ goes from 1 to 0).

````{note}
If the switching required is from state $a$ to state $a+b$, then we can use this pair style with a polynomial function of $\lambda$. In this case we need to pass $\lambda^N$ as an argument to the pair style, e.g., with $N=3$:

```
variable lam file lambda.dat
variable lambda equal ${lam}*${lam}*${lam}

pair_style lj/BGlcleavs3  2.3 2.5 1.0 1.0
pair_coeff 1 1 1.0 1.0
pair_coeff 2 2 1.0 1.0

pair_coeff 1 2 ${lambda} 1.0
```
However, the work calculated  by using the compute [paircleav](./compute_pcleav.md) needs to be adjusted by the correct derivative with respect to $\lambda$. E.g., in the case just showed we need to multiply the work by $N \lambda^{N-1}$. For this reason, even if it possible in principle to use this pair style with a polynomial function of $\lambda$ with $N>1$, we suggest for $N>1$ to use the pair_style [lj/BGNlcleavs3](./pair_BGs3Nlam.md). 
````

* We want to describe the simultaneous switch between two different states of the system, e.g., $a$ and $b$. 
Let us assume the total Hamiltonian of the system, $H^{T}(\lambda)$, is defined as 

$$
	H^{T}(\lambda)= \lambda H_{a} + (1-\lambda) H_{b} 
$$  

for $\lambda \in [0,1]$, where we dropped the dependence on $\mathbf{p}$ and $\mathbf{q}$ for simplicity.
In this latter case, the command Dfrac needs to be used when we are passing to the pair_style an expression like $(1-\lambda)$. Here, we need to set Dfrac=-1 for some pair coefficients.

```
variable lambda file lambda.dat
variable minl   equal 1-lambda

pair_style lj/BGcleavs3  2.3 2.5 1.0 1.0
pair_coeff 1 1 1.0 1.0
pair_coeff 2 2 1.0 1.0
pair_coeff 3 3 1.0 1.0

pair_coeff 1 2 ${lambda} 1.0
pair_coeff 1 3 ${minl} -1.0
pair_coeff 2 3 1.0 1.0
```

The "self" interactions are not modified during the run, as well as the interactions between atoms of type 2 and 3. However, interactions between atoms of type 1 and 2 are `switched-on` (`switched-off`) and interactions between atoms of type 1 and 3 are  `switched-off` (`switched-on`) if $\lambda$ goes from 0 to 1 (from 1 to 0).



````{warning}
This pair style can be used to simulate a polynomial switching from state $a$ to state $b$ in the form given by:

$$
	H^{T}(\lambda)= \lambda^N H_{a} + (1-\lambda^N) H_{b} 
$$  

by defining the polynomial in $\lambda$ outside the pair style, e.g.:

```
variable lam file lambda.dat
variable N equal 4
variable lambda equal exp($N*log(${lam}))
variable minl   equal 1-${lambda}

pair_style lj/BGcleavs3  2.3 2.5 1.0 1.0
pair_coeff 1 1 1.0 1.0
pair_coeff 2 2 1.0 1.0
pair_coeff 3 3 1.0 1.0

pair_coeff 1 2 ${lambda} 1.0
pair_coeff 1 3 ${minl} -1.0
pair_coeff 2 3 1.0 1.0
```

However, the energy obtained through the compute [paircleav](./compute_pcleav.md) must be multiplied by $N\lambda^{N-1}$ as the pair style does not recognize that $\lambda$ is not a linear parameter. 
Note that, if a scaling of the form 

$$
	H^{T}(\lambda)= \lambda^N H_{a} + (1-\lambda)^N H_{b} 
$$  

is needed, while in principle can be also defined with this pair style, we suggest to use the pair_style [lj/BGNlcleavs3](./pair_BGs3Nlam.md) instead. The reason is that, with such a scaling the coefficients for $H_{a}$ interactions is going to be $N\lambda^{N-1}$, while the coefficients for $H_{b}$ interactions is given by $-N(1-\lambda)^{N-1}$.  The pair_style [lj/BGNlcleavs3](./pair_BGs3Nlam.md) is able to automatically account for this.
````



```{footbibliography}

```
