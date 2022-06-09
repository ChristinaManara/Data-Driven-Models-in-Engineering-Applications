###############################  Numerical Laplace Equation Solution using Finite Difference Method  ###############################
import numpy as np
import matplotlib.pyplot as plt

# Part A 

# Maximum iteration
max_Iter = 100

# Boundaries 
T_left = 0
T_right = 0
T_bottom = 0
T_top = 100

# Dimensions of plate -> 1m x 1m and delta
dim_X = 40
dim_Y = 40
delta = 1

# Initial guess of interior grid
T_guess = 30

# Set array size and set the interior value with Tguess
T = np.empty((dim_X, dim_Y))
T.fill(T_guess)

# Set Boundary condition
T[(dim_Y-1):, :] = T_top
T[:1, :] = T_bottom
T[:, (dim_X-1):] = T_right
T[:, :1] = T_left

# Iteration (We assume that the iteration is convergence in maxIter = 500)
for iteration in range(0, max_Iter):
    for i in range(1, dim_X-1, delta):
        for j in range(1, dim_Y-1, delta):
            T[i, j] = 0.25 * (T[i+1][j] + T[i-1][j] + T[i][j+1] + T[i][j-1])
    if iteration == 1:
        print(T)
        break

print("Iteration finished")
