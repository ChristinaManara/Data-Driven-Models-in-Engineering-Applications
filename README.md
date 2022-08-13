# Data-Driven-Models-in-Engineering-Applications


# Exercise 1

A) Consider the following stochastic field: 

$$
E(x)=10(1+f(x))
$$

where f(x) is a zero-mean stationary Gaussian field with unit variance and $x \in[0,5]$(m). The autocorrelation function for f is $R_{f}(\tau)=\exp (-|\tau| / 2)$.

1. Use the Karhunen-Loeve series expansion method to generate N = 5000 realizations of the field E(x). 
2. Justify the number of terms you retained in the KL-expansion.
3. Calculate the ensemble average and the ensemble variance from these realizations. To which values would they converge as we increase the number N of realizations?

B) Consider the zero-mean Gaussian process X(t), t∈[0,10] (sec),  which has the following one-sided power-spectrum

$$
G(\omega)= \begin{cases}\omega-1, & 1 \leq \omega \leq 2 \\ 3-\omega, & 2<\omega \leq 3 \\ 0, & \text { otherwise }\end{cases}
$$

1. Use the Spectral Representation method to generate N=5000 time-histories (realizations) of the process X(t).
2. Calculate the ensemble average and the ensemble variance from these time-histories. To which values would they converge as we increased the number N of realizations.
3. Calculate the temporal average and temporal variance from a single realization. What do you observe?




# Exercise 2

Consider a rectangular plate of dimensions 1 m × 1 m. There is a candle located at position (0.55, 0.45) below the plate which heats it. The steady state equation that describes the temperature field T(x,y) along the plate at thermal equilibrium is:

$$
-\left(\frac{\partial^{2} T}{\partial x^{2}}+\frac{\partial^{2} T}{\partial y^{2}}\right)=f(x, y)
$$

where f(x,y) is the external heat source due the candle, given by:

$$
f(x, y)=100 \cdot \exp \left(-\frac{(x-0.55)^{2}+(y-0.45)^{2}}{r}\right)
$$

with $r \sim N(0.05,0.005)$ being a normal random variable. 
At all edges of the plate we assign T=0 (Dirichlet Boundary conditions).

1. 	Discretize the plate using a 40×40 grid, with h=Δx=Δy=1/40. Using the following second-order central difference approximation to 2nd derivatives

- $\frac{\partial^{2} T}{\partial x^{2}} \approx \frac{T(x+h, y)-2 T(x, y)+T(x-h, y)}{h^{2}}$
- $\frac{\partial^{2} T}{\partial y^{2}} \approx \frac{T(x, y+h)-2 T(x, y)+T(x, y-h)}{h^{2}}$

you can get the finite difference scheme:

$$
-T(x+h, y)-T(x, y+h)+4 T(x, y)-T(x-h, y)-T(x, y-h) \approx h^{2} f(x, y)
$$

Derive the linear system of equations Kx=b for the problem. (Comment: Instead of this scheme, you can use quadrilateral finite elements to derive a linear system for this problem, or any other scheme of your choice.)

2. Perform Monte Carlo simulation to obtain the probability density function of the temperature at the midpoint of the plate (0.5, 0.5).

3. Perform a small number of deterministic simulations for different values of r, and use these solutions as your initial data set. Implement the PCA/POD method to reduce the dimensionality of the linear system that describes the problem and perform the Monte Carlo simulation on the reduced system. Compare the pdf of T at point (0.5,0.5) to the one from the previous question.


# Exercise 3

Extreme Driving Detection

You are given the records of the sensors of a driver’s smartphone (accelerometer, gyroscope and GPS) for a number of trips. The data are in a completely anonymized format.
You are asked to detect extreme turning movements. Compare the results of, at least, three different algorithms.

