###############################  Numerical Laplace Equation Solution using Finite Difference Method  ###############################
from itertools import count
from scipy.linalg import block_diag
import numpy as np
import matplotlib.pyplot as plt

# Part A 

N = 40
Temp = 100
nblock = N - 2

def K(nblock):
    # Diagonal has -4 in every cell 
    K_diag = -4 * np.eye(nblock) 

    # Create a list with size that is defined by (nblock -1) and fill the diagonal with 1 -> above the main diagonal
    K_upper = np.diag([1] * (nblock - 1), 1) 

    # Create a list with size that is defined by (nblock -1) and fill the diagonal with -1 -> above the main diagonal
    K_lower = np.diag([1] * (nblock - 1), -1)

    # Creation of matrix B
    B = K_diag + K_upper + K_lower 

    # Creat a list [B,B,B,...,B] with nblock Bs
    blst = [B] * nblock

    # Unpack the list of diagonal blocks ’blst’
    A = block_diag(*blst)

    # Upper diagonal array offset by nblock: we’ve got (nblock-1) I blocks
    # each containing nblock ones
    Dupper = np.diag(np.ones(nblock * (nblock - 1)), nblock)

    # Lower diagonal array offset by -nblock
    Dlower = np.diag(np.ones(nblock * (nblock - 1)), -nblock)

    A += Dupper + Dlower
    
    return A


def b(nblock, Temp):
    b = np.zeros(nblock**2)
    b[:nblock] = -Temp
    print("b = \n", b)
    return b

def solve_equation():
    u = np.linalg.solve(K, b)
    T = u.reshape((nblock, nblock))
    Tfull = embed(T)
    print("Tfull = \n", Tfull)






# # Maximum iteration
# max_Iter = 100

# # Boundaries 
# T_left = 0
# T_right = 0
# T_bottom = 0
# T_top = 100

# # Dimensions of plate -> 1m x 1m and delta
# dim_X = 40
# dim_Y = 40
# delta = 1
# h = delta/dim_X # or dim_Y it does not matter, type -> h=1/40=Δx=Δy
# r = 0.05 # follow normal distribution r ~ N(0.05, 0.005)

# # Initial guess of interior grid
# T_guess = 30

# # Set array size and set the interior value with Tguess
# T = np.empty((dim_X, dim_Y))
# T.fill(T_guess)

# # Set Boundary condition
# T[(dim_Y-1):, :] = T_top
# T[:1, :] = T_bottom
# T[:, (dim_X-1):] = T_right
# T[:, :1] = T_left

# # Set array b 
# # b = np.empty(dim_X * dim_Y)
# # counter = 0

# # for i in range(1, dim_X-1, delta):
# #     for j in range(1, dim_Y-1, delta):
# #        b[counter]  = (h ** 2) * 100.0 * np.exp(-((x - 0.55) ** 2 + (y - 0.45) ** 2 ) / r)
# #        counter += 1


# # Iteration (We assume that the iteration is convergence in maxIter = 500)
# # for iteration in range(0, max_Iter):
# #     for i in range(1, dim_X-1, delta):
# #         for j in range(1, dim_Y-1, delta):
# #             T[i, j] = 0.25 * (T[i+1][j] + T[i-1][j] + T[i][j+1] + T[i][j-1])
# #     if iteration == 1:
# #         print(T)
# #         break

# # print("Iteration finished")






