# Simulation of the accretion disk of the binary star system PHL1445 using parallel version of SPH algorithm

![alt-text](https://github.com/iliazhav/astro/blob/main/pictures/ezgif.com-crop.gif)


A binary star system is a system consisting of two gravitationally bound stars that revolve in closed orbits around a common center of mass. Сlose doubles systems - a kind of binary systems in which the components included in it can exchange mass. This project considers the stage of evolution, in which the flow of matter occurs from a red dwarf (donor) to the vicinity of a white dwarf (accretor). In this work, to simulate the transfer of matter, the method of smoothed particle hydrodynamics (SPH) is used. Taking into account the cooling function of the ionized gas and the features of the radiation of stellar surfaces, it is possible to obtain a light curve, which is the main way to study the characteristics of such systems. An important component in the formation of data curves is the luminosity at the point of maximum accumulation of matter in the accretion disk.


<b>Ballistic model of an accretion disk</b>

Let's introduce a coordinate system for components with masses $`m_1`$ and $`m_2`$. We will choose the center of reference at the center of mass of two stars. We write down the potential associated with the gravitational fields of the components and the centrifugal force:

```math
Phi = -\frac{Gm_1}{p_1} - \frac{Gm_2}{p_2} - \frac12 \omega^2(x^2 + y^2),
```
where $`p_1 = \sqrt{(x - a_1)^2 + y^2 + z^2}`$, $`p_1 = \sqrt{(x - a_2)^2+y^2+z^2}`$, $`\omega`$ - rotation frequency, $`a_1`$ and $`a_2`$ - distances from stars to the center of mass.

<img src="https://github.com/iliazhav/astro/blob/main/pictures/Motion2D1.jpg" width="50%" height="50%">

The Lagrangian of this system:
```math
L=\frac12 m_3 ((\dot{x})^2 + (\dot{y})^2 + 2x\dot{y}\omega - 2\dot{x}y\omega + \omega^2(x^2+y^2)) - G \frac{m_3 m_1}{p_1} - G \frac{m_3 m_2}{p_2}.
```
Lagrange equation:
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
The ballistic approximation of such a problem makes it possible to estimate the approximate size and location of the region of maximum luminosity and to determine the phases of the eclipse.

<b>Basics of SPH Method</b>

Let's complicate the model - let's introduce the interaction between individual clumps of matter. We will emit particles from the Lagrange point $`L_1`$, which will replace the continuous field.

<img src="https://github.com/iliazhav/astro/blob/main/pictures/BigRoche.png" width="30%" height="30%">

This method uses gas-dynamic equations in Lagrangian variables and does not require the introduction of a difference grid.
Consider the identity
```math
f(\vec{r}) = \int\limits_V f(\vec{r}^,) \delta(\vec{r} - \vec{r}^,)d\vec{r}^,,
```
where $`f (\vec{r})`$ is some scalar function defined in a three-dimensional coordinate system, $`\delta(\vec{r})`$ is the Dirac delta function, and $`\vec{r} ^,`$ - dummy variable, which is defined in volume $`V`$, $`\vec{r}`$ - coordinate.

One can generalize the delta function to the smoothing kernel $W$ with
smoothing length $`h`$, such that

```math
\lim\limits_{h\rightarrow0} W(\vec{r}, h) = \delta(\vec{r}),
```

for normalization we require

```math
\int\limits_V W(\vec{r}, h) d\vec{r}^, = 1.
```

In passing from a continuous to a discrete representation, imagine that the computational domain is filled with particles with coordinates $`\boldsymbol{r}_{i}`$ and masses $`m = \rho (\vec{r}^,)d\vec{r}^,`$. Then for $`f`$ the integral can be replaced by the sum
```math
f(\mathbf{r})\approx\sum_{i} m_{i} \frac{f_{i}}{\rho_{i}} W\left(\mathbf{r}-\mathbf{r}_ {i}, h\right),
```
where the index $`i`$ determines the number of the particle.

The following relations were also obtained in the work:

Density equation:
```math
\rho_{j}=\sum_{i} m_{i} W\left(\mathbf{r}_{j}-\mathbf{r}_{i}, h\right) =\sum_{i} m_{i} W_{j i}.
```
Momentum equation: 
```math
\frac{d \boldsymbol{v_{i}}}{d t}=-\sum_{j} m_{j}\left(\frac{P_{j}}{\rho_{j}^{2}}+\frac{P_{i}}{\rho_{i}^{2}}\right) \boldsymbol{\nabla_{i}} W_{i j}.
```
The pressure force acting from one particle to another is expressed as
```math
\boldsymbol{F}_{i j}=m_{i} m_{j}\left(\frac{P_{i}}{\rho_{i}^{2}}+\frac{P_{j}}{ \rho_{j}^{2}}\right)\left(\mathbf{r}_{i}-\mathbf{r}_{j}\right) W^{\prime}(r, h) / r.
```
The energy equation in the SPH formalism:
```math
\begin{aligned}
\frac{\mathrm{d} \epsilon_{j}}{\mathrm{~d} t} =\frac{P_{j}}{\rho_{j}^{2}} \frac{\mathrm{d } \rho_{j}}{\mathrm{~d} t} =\frac{P_{j}}{\rho_{j}^{2}} \sum_{i} m_{i}\left(\mathbf {v}_{j}-\mathbf{v}_{i}\right) \cdot \nabla_{i} W_{j i}.
\end{aligned}
```
The mutual attraction of particles (self-gravity) in the SPH method can be taken into account. To do this, it is necessary to introduce an additional term $`\boldsymbol{f}^{\text {grav }}`$, which will describe the gravitational attraction of each particle to the rest:
```math
\frac{d \boldsymbol{v}_{a}}{d t}=\boldsymbol{f}^{\text {hydro }}+\boldsymbol{f}^{\text {grav }}.
```

<b>Computing</b>

During testing, various methods for the numerical solution of differential equations were tested. The use of the fourth-order Runge-Kutta method made it possible to avoid the rapid accumulation of errors in the region of self-intersection of trajectories. The iterative process for the method is as follows:
```math
\vec{y_{n+1}} = \vec{y_n} + \frac h6 (\vec{k_1} + \vec{k_2} + \vec{k_3} + \vec{k_4}),
```
where $`y_i`$, $`k_i`$ are vectors.

Let there be $n$ particles in the computational domain at the last cycle. The SPH method, which obliges to consider $`\displaystyle n(n-1)`$ combinations, leads to the complexity of the algorithm $`O(n^2)`$. In this work, calculations were carried out on one node in several threads, which significantly reduced the program execution time. For example, the processing time of a model with a maximum number of particles equal to 500 was 56 seconds on one thread, and 5 seconds on 18 threads.

<img src="https://github.com/iliazhav/astro/blob/main/pictures/proj1.png" width="30%" height="30%">

Ilia Z.
