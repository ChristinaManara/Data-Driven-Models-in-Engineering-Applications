###############################  Numerical Laplace Equation Solution using Finite Difference Method  ###############################
from cgi import print_directory
from difflib import restore
from itertools import count
from pipes import Template
from re import M
from unittest import result
from scipy.linalg import block_diag
import numpy as np
import matplotlib.pyplot as plt
import sys
import seaborn as sns
import math
import scipy.stats as stats


# Part A 

N = 40
Temp = 0
nblock = N - 1
dim_X = 40
dim_Y = 40
delta = 1
max_Iter = 5000
h = delta/dim_X # or dim_Y it does not matter, type -> h=1/40=Δx=Δy

def K(nblock):
    # Diagonal has -4 in every cell 
    K_diag = 4 * np.eye(nblock**2) 

    sh_diagonal = - 1 * np.ones(nblock**2 - 1)

    sh_diagonal_N = - 1 * np.ones(nblock**2 - nblock)

    # Create a list with size that is defined by (nblock - 1) and fill the diagonal with 1 -> above the main diagonal
    K_upper = np.diag(sh_diagonal, 1) 

    # Create a list with size that is defined by (nblock - 1) and fill the diagonal with -1 -> above the main diagonal
    K_lower = np.diag(sh_diagonal, -1)

    K_upper_N = np.diag(sh_diagonal_N, nblock) 

    K_lower_N = np.diag(sh_diagonal_N, -nblock) 

    # Creation of matrix B
    B = K_diag + K_upper + K_lower + K_upper_N + K_lower_N

    return B

# Creation the b vector for the system of linear equations

def b(nblock, Temp):
    b = np.zeros(nblock**2)
    b[:nblock] = -Temp
    print("b = \n", b)
    return b

def solve_equation(K, b, nblock, Temp):
    u = np.linalg.solve(K, b)
    T = u.reshape((nblock, nblock))
    return T


K_arr = K(nblock)
b_vec = b(nblock, Temp)
T_solved = solve_equation(K_arr, b_vec, nblock, Temp)

# Colorplot 
# def colorplot(T_solved, N, Temp):
#     x = y = np.linspace(0, 1, N-1)
#     X, Y = np.meshgrid(x,y)
#     # h = plt.contourf(X, Y, T_filled)
#     # plt.axis('scaled')
#     # plt.colorbar()
#     # plt.show()
    
#     plt.figure(1)
#     plt.clf()
#     plt.pcolor(X, Y, T_solved)
#     plt.axis("scaled")
#     plt.colorbar()
#     plt.xlabel("x (m)")
#     plt.ylabel("y (m)")
#     plt.title("T(x,y) on %dx%d grid" % (N,N))
#     plt.show()

# colorplot(T_solved, N, h)

# Part B

def f_x_y(h, r):
    f = []
    for i in range(1, dim_X, delta):
        for j in range(1, dim_Y, delta):
            res = (h ** 2) * 100.0 * np.exp(-((i/(N-1) - 0.55) ** 2 + (j/(N-1) - 0.45) ** 2 ) / r)
            f.append(res)
    return f


mid_point = int((nblock**2 - 1) / 2)
mean = 0.05
std = 0.005
K_arr_inv = np.linalg.inv(K_arr)       

def monte_carlo():
    results = []
    for i in range(max_Iter):
        r = np.random.normal(mean, std)
        tmp = f_x_y(h, r)
        simulation = np.dot(K_arr_inv, tmp)
        # Append the point we are intrested in
        results.append(simulation[mid_point])

def plot_samples(results):

    plt.figure(figsize = (10, 10))
    sns.distplot(results, label = "PDF of temperature in (0.5, 0.5)")
    plt.vlines(np.mean(results), ymin = 0, ymax = 1.6, colors = "green", label = "Mean of PDF", linestyles = "dotted")


    mu = np.mean(results)
    variance = np.var(results)
    sigma = math.sqrt(variance)

    x_1 = np.linspace(mu - 5*sigma, mu + 5*sigma, 100) # truncated normal

    plt.plot(x_1, stats.norm.pdf(x_1, mu, sigma), label = "Truncated normal to 5 SD")
    plt.legend()
    plt.show()
    print("mean: {}\nvariance: {}".format(np.mean(results), np.var(results)))


# Part C

pca_iter = 100
reduced_res = []
for i in range(max_Iter):
    r = np.random.normal(mean, std)
    tmp = f_x_y(h, r)
    simulation = np.dot(K_arr_inv, tmp)
    # Append the point we are intrested in
    reduced_res.append(simulation[mid_point])

pca()

def pca(X , num_components):
     
    # Step-1
    X_meaned = X - np.mean(X , axis = 0)
     
    # Step-2
    cov_mat = np.cov(X_meaned , rowvar = False)
     
    # Step-3
    eigen_values , eigen_vectors = np.linalg.eigh(cov_mat)
     
    # Step-4
    sorted_index = np.argsort(eigen_values)[::-1]
    sorted_eigenvalue = eigen_values[sorted_index]
    sorted_eigenvectors = eigen_vectors[:,sorted_index]
     
    # Step-5
    eigenvector_subset = sorted_eigenvectors[:,0:num_components]
     
    # Step-6
    X_reduced = np.dot(eigenvector_subset.transpose() , X_meaned.transpose() ).transpose()
     
    return X_reduced




