# Simulation of the accretion disk of the binary star system PHL 1445 using parallel version of SPH algorithm

![alt-text](https://github.com/iliazhav/astro/blob/main/pictures/ezgif.com-crop.gif)


A binary star system is a system consisting of two gravitationally bound stars that revolve in closed orbits around a common center of mass. Сlose doubles systems - a kind of binary systems in which the components included in it can exchange mass. This project considers the stage of evolution, in which the flow of matter occurs from a red dwarf (donor) to the vicinity of a white dwarf (accretor). In this work, to simulate the transfer of matter, the method of smoothed particle hydrodynamics (SPH) is used, which allows you to get a more detailed picture, compared with the ballistic approximation. Taking into account the cooling function of the ionized gas and the features of the radiation of stellar surfaces, one can obtain a light curve, which is the main way to study the characteristics of such systems. An important component in the formation of data curves is the luminosity at the point of maximum accumulation of matter in the accretion disk. The ballistic approximation of such a problem makes it possible to estimate the approximate size and location of the region of maximum luminosity and to determine the phases of the eclipse.


<b>Ballistic model of an accretion disk</b>

Let's introduce a coordinate system for components with masses $`m_1`$ and $`m_2`$. We will choose the center of reference at the center of mass of two stars. We write down the potential associated with the gravitational fields of the components and the centrifugal force:

```math
\begin{equation}
\Phi = -\frac{Gm_1}{p_1} - \frac{Gm_2}{p_2} - \frac12 \omega^2(x^2 + y^2),
\end{equation}
```
where $`p_1 = \sqrt{(x - a_1)^2 + y^2 + z^2}`$, $`p_1 = \sqrt{(x - a_2)^2+y^2+z^2}`$, $\omega$ - rotation frequency, $`a_1`$ and $`a_2`$ - distances from stars to the center of mass.

<img src="https://github.com/iliazhav/astro/blob/main/pictures/Motion2D1.jpg" width="50%" height="50%">

The Lagrangian of this system:
```math
L=\frac12 m_3 ((\dot{x})^2 + (\dot{y})^2 + 2x\dot{y}\omega - 2\dot{x}y\omega + \omega^2(x^2+y^2)) - G \frac{m_3 m_1}{p_1} - G \frac{m_3 m_2}{p_2}.
```
Lagrange equation of the second kind
```math
\frac{d}{dt}\left(\frac{\partial L}{\partial \dot{q_i}}\right) - \frac{\partial L}{\partial q_i} = 0.
```

Taking into account the classical concepts of the 3-body problem and using the theorems of theoretical mechanics (Kepler's laws, Coriolis force), it is possible to create a system of differential equations of motion for smallest component of the system:

```math
\begin{cases}
\ddot{x} = 2\dot{y} + x - \displaystyle\frac{(1-\mu)(x+\mu)}{p_1^3} - \frac{\mu (x-(1-\mu)}{p_2^3},\\
\ddot{y} = -2\dot{x} + y - \displaystyle\frac{(1-\mu)y}{p_1^3} - \frac{\mu y}{p_2^3},\\
\ddot{z} = -\displaystyle\frac{(2\pi)^2(1-\mu)z}{p_1^3}-\displaystyle\frac{(2\pi)^2\mu z}{p_2^3},
\end{cases}
```
where $`x`$, $`y`$, $`z`$ - coordinates of each particle, $`\mu`$ - the ratio of stars' masses, $`p_1`$ and $`p_2`$ - ratio of masses.

This method uses gas-dynamic equations in Lagrangian variables and does not require the introduction of a difference grid [6].
