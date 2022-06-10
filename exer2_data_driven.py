###############################  Numerical Laplace Equation Solution using Finite Difference Method  ###############################
from difflib import restore
from itertools import count
from unittest import result
from scipy.linalg import block_diag
import numpy as np
import matplotlib.pyplot as plt

# Part A 

N = 40
Temp = 100
nblock = N - 2
dim_X = 40
dim_Y = 40
delta = 1
max_Iter = 10000
h = delta/dim_X # or dim_Y it does not matter, type -> h=1/40=Δx=Δy

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

# Creation the b vector for the system of linear equations

def b(nblock, Temp):
    b = np.zeros(nblock**2)
    b[:nblock] = -Temp
    print("b = \n", b)
    return b

# We want to add the boundary conditions to the array, and we have again a 40x40 array instead of just 38x38 of unknown values. 
# This can be done by embedding the unknown array into a larger array. 

def embed(T, Temp):
    N = T.shape[0] + 2
    T_filled = np.zeros((N,N))
    T_filled[0] = Temp
    T_filled[1:-1, 1:-1] = T
    return T_filled

def solve_equation(K, b, nblock, Temp):
    u = np.linalg.solve(K, b)
    T = u.reshape((nblock, nblock))
    T_filled = embed(T, Temp)
    print("T_filled = \n", T_filled)
    return T


K_arr = K(nblock)
print("S", len(K_arr))
b_vec = b(nblock, Temp)
T_filled = solve_equation(K_arr, b_vec, nblock, Temp)



# Colorplot 
def colorplot(T_filled, N, h):
    x = y = np.linspace(0, 1, N-2)


    X, Y = np.meshgrid(x,y)
    h = plt.contourf(X, Y, T_filled)
    plt.axis('scaled')
    plt.colorbar()
    plt.show()
    
    # plt.figure(1)
    # plt.clf()
    # plt.pcolor(X, Y, T_filled)
    # plt.axis("scaled")
    # plt.colorbar()
    # plt.xlabel("x (m)")
    # plt.ylabel("y (m)")
    # plt.title("T(x,y) on %dx%d grid" % (N,N))
    # plt.show()

# colorplot(T_filled, N, h)


# Part B

def f_x_y(h, r):
    f = []
    for i in range(1, dim_X-1, delta):
        for j in range(1, dim_Y-1, delta):
            res = (h ** 2) * 100.0 * np.exp(-((i/N - 0.55) ** 2 + (j/N - 0.45) ** 2 ) / r)
            f.append(res)
    return f
       

mid_point = int((nblock**2 - 1) / 2)
mean = 0.05
std = 0.005
K_arr_inv = np.linalg.inv(K_arr)
results = []
for i in range(max_Iter):
    r = np.random.normal(mean, std)
    tmp = f_x_y(h, r)
    simulation = np.dot(K_arr_inv, tmp)
    results.append(simulation)

print(mid_point)


def plot_samples(results, ax, orientation='vertical'):
    
    ax.hist(results, bins=50, density=True, label="Monte Carlo Samples",
           orientation=orientation)

plt.figure(2)  
fig, ax = plt.subplots()
plot_samples(results, ax)
ax.set_yticks(())
plt.show()

# Part C









