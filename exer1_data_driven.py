from tkinter import N
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy import optimize

# Christina Manara 03400142
# Part A
a = 2.5 
b = 2 
X = np.linspace(-a, a, 300)

def lamda_n(w):
    return (2 * b) /( 1 + w**2 * b**2)

def ln(w):
    return 1 / np.sqrt(a - np.sin(2 * w * a) / (2 * w))

def cn(w):
    return 1 / np.sqrt(a + np.sin(2 * w * a) / (2 * w))

def phi_odd(w, t):
    return cn(w) * np.cos(w * t)

def phi_even(w, t):
    return ln(w) * np.sin(w * t) 

def eq_odd(w):
    return (1 / b) - (w * np.tan(a * w))
 
def eq_even(w):
    return ((1 / b) * np.tan(a * w)) + w

def solve_omg_even_or_odd(n):
    if (n % 2) != 0:
        tmp = optimize.fsolve(eq_odd, (n-0.5) * np.pi/a - 0.1/n)
    else:
        tmp = optimize.fsolve(eq_even, (n-0.5) * np.pi/a + 0.1/n)
    return float(tmp)
     
def eigenvalues():
    eigenvalues = []
    sol_w = []
    for n in range(1, 50):
        sol = solve_omg_even_or_odd(n)
        sol_w.append(sol)
        l = lamda_n(sol)
        eigenvalues.append(l) 
        print(eigenvalues)
    return eigenvalues, sol_w

values, solutions = eigenvalues()
plt.plot(range(1, 50), values)
plt.xlabel('KL Order (M)')
plt.ylabel('Eigenvalues')
plt.show()

def eigenfunction(values, solutions):
    M = 10
    values = values[:M]
    solutions = solutions[:M]
    eigenfun = []
    for n in range(1, M + 1):
        w_n = solutions[n-1]
        if (n % 2) != 0: 
            eigenfun.append(phi_odd(w_n, X))
        else: 
            eigenfun.append(phi_even(w_n, X))
    return values, eigenfun


eigenval, eigenfunctions = eigenfunction(values, solutions)
#print(eigenfunctions)

def realization(M):
    r = 0
    z_n = np.random.normal(size = M)
    for m in range(0, M):
        tmp = math.sqrt(eigenval[m]) * eigenfunctions[m] * z_n[m]
        r = r + tmp 
    return r

def multiple_realizations(M, N, points = 300):
    mul_realizations = np.zeros((N, points))
    for i in range(N):
        mul_realizations[i, :] = realization(M)
    return mul_realizations

# print(multiple_realizations(M=10, N=5000, points = 300))

def plot_samples(M = 10):
    R_5000 = 10 * (1 + multiple_realizations(M, N=5000))
    R_15000 = 10 * (1 + multiple_realizations(M, N=15000))
    R_45000 = 10 * (1 + multiple_realizations(M, N=45000))
    R_75000 = 10 * (1 + multiple_realizations(M, N=75000))
    R_110000 = 10 * (1 + multiple_realizations(M, N=110000))

    mean_5k = np.mean(R_5000, axis = 0)
    var_5k = np.var(R_5000, axis = 0)
    # print(var_5k)

    mean_15k = np.mean(R_15000, axis = 0)
    var_15k = np.var(R_15000, axis = 0)

    mean_45k = np.mean(R_45000, axis = 0)
    var_45k = np.var(R_45000, axis = 0)

    mean_75k = np.mean(R_75000, axis = 0)
    var_75k = np.var(R_75000, axis = 0)

    mean_110k = np.mean(R_110000, axis = 0)
    var_110k = np.var(R_110000, axis = 0)

    plt.plot(X + 2.5, var_5k, label="N=5000, M=10")

    plt.plot(X + 2.5, var_15k, label="N=15000, M=10")

    plt.plot(X + 2.5, var_45k, label="N=45000, M=10")

    plt.plot(X + 2.5, var_75k, label="N=750000, M=10")

    plt.plot(X + 2.5, var_110k, label="N=110000, M=10")

    plt.legend()
    plt.title('Ensemble Variance')
    plt.xlabel('x')
    plt.ylabel('Variance')
    plt.show()
    
    plt.plot(X + 2.5, mean_5k, label="N=5000, M=10")

    plt.plot(X + 2.5, mean_15k, label="N=15000, M=10")

    plt.plot(X + 2.5, mean_45k, label="N=45000, M=10")

    plt.plot(X + 2.5, mean_75k, label="N=750000, M=10")

    plt.plot(X + 2.5, mean_110k, label="N=110000, M=10")

    plt.legend()
    plt.title('Ensemble Average')
    plt.xlabel('x')
    plt.ylabel('Average')
    plt.show()

plot_samples()

# Part B

N = 100 
wu = 3
dw = wu / N
t = np.linspace(0, 10, 300)

# Power spectrum G
def G(w):
    if (w >= 1 and w <= 2):
        return w - 1
    elif (w > 2 and w <= 3):
        return 3 - w
    else:
        return 0

wn = []
G_values = []
An = []

for n in range(N):
    wn.append(n * dw)
    G_values.append(G(n * dw))
    An.append(math.sqrt(G(n * dw) * dw))


def realization_b():
    phi_n = np.random.uniform(low = 0, high = 2 * np.pi, size = N)
    r = 0
    for m in range(N): 
        sum = An[m] * np.cos(wn[m] * t + phi_n[m])
        r = r + sum 
    return math.sqrt(2) * r

def multiple_realizations_b(N):
    multiple_r = np.zeros((N, len(t)))
    for i in range(N):
        multiple_r[i, :] = realization_b()
    return multiple_r


def plot_samples_b_part(M = 10):
    R_5000 = multiple_realizations_b(N=5000)
    R_15000 = multiple_realizations_b(N=15000)
    R_45000 = multiple_realizations_b(N=45000)
    R_75000 = multiple_realizations_b(N=75000)
    R_110000 = multiple_realizations_b(N=110000)

    mean_5k = np.mean(R_5000, axis = 0)
    var_5k = np.var(R_5000, axis = 0)
    # print(var_5k)

    mean_15k = np.mean(R_15000, axis = 0)
    var_15k = np.var(R_15000, axis = 0)

    mean_45k = np.mean(R_45000, axis = 0)
    var_45k = np.var(R_45000, axis = 0)

    mean_75k = np.mean(R_75000, axis = 0)
    var_75k = np.var(R_75000, axis = 0)

    mean_110k = np.mean(R_110000, axis = 0)
    var_110k = np.var(R_110000, axis = 0)

    plt.plot(X + 2.5, var_5k, label="N=5000, M=10")

    plt.plot(X + 2.5, var_15k, label="N=15000, M=10")

    plt.plot(X + 2.5, var_45k, label="N=45000, M=10")

    plt.plot(X + 2.5, var_75k, label="N=750000, M=10")

    plt.plot(X + 2.5, var_110k, label="N=110000, M=10")

    plt.legend()
    plt.title('Ensemble Variance')
    plt.xlabel('x')
    plt.ylabel('Variance')
    plt.show()
    
    plt.plot(X + 2.5, mean_5k, label="N=5000, M=10")

    plt.plot(X + 2.5, mean_15k, label="N=15000, M=10")

    plt.plot(X + 2.5, mean_45k, label="N=45000, M=10")

    plt.plot(X + 2.5, mean_75k, label="N=750000, M=10")

    plt.plot(X + 2.5, mean_110k, label="N=110000, M=10")

    plt.legend()
    plt.title('Ensemble Average')
    plt.xlabel('x')
    plt.ylabel('Average')
    plt.show()

plot_samples_b_part()