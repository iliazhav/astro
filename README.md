# Simulation of the accretion disk of the binary star system PHL 1445 using parallel version of SPH algorithm

![alt-text](https://github.com/iliazhav/astro/blob/main/pictures/ezgif.com-crop.gif)


A binary star system is a system consisting of two gravitationally bound stars that revolve in closed orbits around a common center of mass. Сlose doubles systems - a kind of binary systems in which the components included in it can exchange mass. This project considers the stage of evolution, in which the flow of matter occurs from a red dwarf (donor) to the vicinity of a white dwarf (accretor). In this work, to simulate the transfer of matter, the method of smoothed particle hydrodynamics (SPH) is used, which allows you to get a more detailed picture, compared with the ballistic approximation. Taking into account the cooling function of the ionized gas and the features of the radiation of stellar surfaces, one can obtain a light curve, which is the main way to study the characteristics of such systems. An important component in the formation of data curves is the luminosity at the point of maximum accumulation of matter in the accretion disk. The ballistic approximation of such a problem makes it possible to estimate the approximate size and location of the region of maximum luminosity and to determine the phases of the eclipse.


<b>Ballistic model of an accretion disk</b>

<img src="https://github.com/iliazhav/astro/blob/main/pictures/Motion2D1.jpg" width="50%" height="50%">

Taking into account the classical concepts of the 3-body problem and using the theorems of theoretical mechanics (Kepler's laws, Coriolis force), it is possible to create a system of differential equations of motion for each component of the system:

```math
\begin{cases}
\ddot{x} = 2\dot{y} + x - \displaystyle\frac{(1-\mu)(x+\mu)}{p_1^3} - \frac{\mu (x-(1-\mu)}{p_2^3},\\
\ddot{y} = -2\dot{x} + y - \displaystyle\frac{(1-\mu)y}{p_1^3} - \frac{\mu y}{p_2^3},\\
\ddot{z} = -\displaystyle\frac{(2\pi)^2(1-\mu)z}{p_1^3}-\displaystyle\frac{(2\pi)^2\mu z}{p_2^3},
\end{cases}
```
where $`x`$, $`y`$, $`z`$ - coordinates of each particle, $`\mu`$ - the ratio of stars' masses, $`p_1`$ and $`p_2`$ - distances to each star.
