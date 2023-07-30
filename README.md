# Simulation of the accretion disk of the binary star system PHL 1445 using parallel version of SPH algorithm

![alt-text](https://github.com/iliazhav/astro/blob/main/ezgif.com-crop.gif)

<b>Ballistic model of an accretion disk</b>

Taking into account the classical concepts of the 3-body problem and using the theorems of theoretical mechanics, it is possible to create a system of differential equations of motion for each component of the system:

```math
\begin{cases}
\ddot{x} = 2\dot{y} + x - \displaystyle\frac{(1-\mu)(x+\mu)}{p_1^3} - \frac{\mu (x-(1-\mu)}{p_2^3},\\
\ddot{y} = -2\dot{x} + y - \displaystyle\frac{(1-\mu)y}{p_1^3} - \frac{\mu y}{p_2^3},\\
\ddot{z} = -\displaystyle\frac{(2\pi)^2(1-\mu)z}{p_1^3}-\displaystyle\frac{(2\pi)^2\mu z}{p_2^3},
\end{cases}
```
where $`x`$, $`y`$, $`z`$ - coordinates of each particle, $`\mu`$ - the ratio of stars' masses, $`p_1`$ and $`p_2`$ - distances to each star.
