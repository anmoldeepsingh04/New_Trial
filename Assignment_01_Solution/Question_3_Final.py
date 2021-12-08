"""
Question3 of Assignment1.
Completed by:
Anmoldeep Singh 180030002
Radhika Gujar   180030022
"""

# Importing relevant libraries and modules
import math as m
import sys
import time
import matplotlib.pyplot as plt
import numpy as np


# Defining function for calculating the spectral radius, which is the largest eigenvalue of a matrix
def spec_radius(a):
    try:
        v = np.linalg.eigvals(a)
        for i in range(0, len(v)):
            v[i] = abs(v[i])
        return max(v)
    except:
        return 1


# defining the Jacobi Algorithm
def jacobi(a, b, n, k):
    u = np.triu(a, 1)  # Upper triangular matrix of Coefficient matrix
    l = np.tril(a, -1)  # Lower triangular matrix of Coefficient matrix
    d = a - u - l  # Diagonal matrix of Coefficient matrix
    z = np.linalg.inv(d)  # Inverting Diagonal matrix
    p = np.matmul(z, -1 * (u + l))  # Matrix multiplication of z and u+l
    q = np.matmul(z, b)  # Matrix multiplication of z and rhs matrix
    x = np.zeros((n, 1))
    error_jac = np.zeros(k)  # for storing error value at each step
    test = np.zeros(k)
    for i in range(0, k):
        temp = x.copy()
        x = np.matmul(p, x) + q
        for j in range(0, num_mesh):
            test[j] = abs(x[j] - temp[j])
        error_jac[i] = max(test)
    return error_jac


# defining the Successive Over-Relaxation Algorithm
def sor(a, b, n, k):
    u = np.triu(a, 1)  # Upper triangular matrix of Coefficient matrix
    l = np.tril(a, -1)  # Lower triangular matrix of Coefficient matrix
    d = a - u - l  # Diagonal matrix of Coefficient matrix
    z = np.linalg.inv(d + l)  # Inverting d+l matrix
    j = np.matmul(z, (-1 * u))  # Matrix multiplication of z and u
    w = 2 / (1 + np.sqrt(1 - spec_radius(j)))  # Calculating the relaxation parameter
    z = np.linalg.inv(d + w * l)
    p = np.matmul(z, ((w - 1) * l - u))
    q = np.matmul(z, b)
    x = np.zeros((n, 1))
    error_sor = np.zeros(k)  # for storing error value at each step
    test = np.zeros(k)
    for i in range(0, k):
        temp = x.copy()
        x = np.matmul(p, x) + q
        for j in range(0, num_mesh):
            test[j] = abs(x[j] - temp[j])
        error_sor[i] = max(test)
    return error_sor


# defining the Gauss-Seidel Algorithm
def seidel(a, b, n, k):
    u = np.triu(a, 1)  # Upper triangular matrix of Coefficient matrix
    l = np.tril(a, -1)  # Lower triangular matrix of Coefficient matrix
    d = a - u - l  # Diagonal matrix of Coefficient matrix
    z = np.linalg.inv(d + l)  # Inverting Diagonal matrix
    p = np.matmul(z, (-1 * u))  # Matrix multiplication of z and u
    q = np.matmul(z, b)  # Matrix multiplication of z and rhs matrix
    x = np.zeros((n, 1))
    error_sed = np.zeros(k)  # for storing error value at each step
    test = np.zeros(k)
    for i in range(0, k):
        temp = x.copy()
        x = np.matmul(p, x) + q
        for j in range(0, num_mesh):
            test[j] = abs(x[j] - temp[j])
        error_sed[i] = max(test)
    return error_sed


# start discritization
num_mesh = 51
del_x = 1.0 / (num_mesh - 1.0)
alpha = -2.0 - del_x ** 2
x_mesh = np.zeros(num_mesh)

# Number of iterations
iter_count = 200

if iter_count < num_mesh:
    exit("Number of mesh points must be less than or equal to iteration counts.")

# Constructing the coefficient and rhs matrix
a = np.zeros((num_mesh, num_mesh))
b = np.zeros((num_mesh, 1))
a[0][0] = 1.0
a[num_mesh - 1][num_mesh - 1] = 1.0

for i in range(1, num_mesh - 1):
    a[i][i] = -2.0 - (del_x ** 2)
    a[i][i - 1] = 1.0
    a[i][i + 1] = 1.0

# Constructing the grid point matrix
for i in range(0, num_mesh - 1):
    x_mesh[i] = i * del_x

for i in range(0, num_mesh):
    b[i][0] = -(del_x ** 2) * (m.pi ** 2 + 1.0) * m.sin(m.pi * x_mesh[i])

# Matrix of iteration counts
iteration = np.zeros(iter_count)
for i in range(0, iter_count):
    iteration[i] = i

# Calling different iterative methods
jac_ans = jacobi(a, b, num_mesh, iter_count)
sed_ans = seidel(a, b, num_mesh, iter_count)
sor_ans = sor(a, b, num_mesh, iter_count)

# Plotting the results
plt.plot(iteration, jac_ans, 'r', label='Jacobi Method')
plt.plot(iteration, sed_ans, 'g', label='Gauss-Seidel Method')
plt.plot(iteration, sor_ans, 'b', label='Successive Over-Relaxation Method')
plt.title('Convergence of Jacobi, Gauss-Seidel and Successive over-relaxation method')
plt.xlabel('Iteration', fontsize=14)
plt.ylabel('L∞-norm', fontsize=14)
plt.legend()
# plt.xscale('log')
# plt.yscale('log')
plt.show()
