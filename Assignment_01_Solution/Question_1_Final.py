"""
Question2 of Assignment1.
Completed by:
Anmoldeep Singh 180030002
Radhika Gujar   180030022
"""

# Importing relevant libraries and modules
import numpy as np
import math as m
import matplotlib.pyplot as plt
import sys


# defining the TDMA Algorithm
def tdma(num, low, up, dia, func, sol):
    d1 = np.zeros(num)
    rhs1 = np.zeros(num)
    d1[0] = dia[0]
    rhs1[0] = func[0]
    for j in range(1, len(dia)):
        d1[j] = dia[j] - (low[j] * up[j - 1] / d1[j - 1])
        rhs1[j] = func[j] - (low[j] * rhs1[j - 1] / d1[j - 1])
    sol[num - 1] = func[num - 1] / dia[num - 1]
    for k in range(len(dia) - 2, -1, -1):
        sol[k] = (rhs1[k] - up[k] * sol[k + 1]) / d1[k]


# Applying the boundary conditions
def apply_bc(u):
    u[0] = 0.0
    if bc == "Dirichlet":
        u[num_mesh - 1] = 0.0
    else:
        sys.exit('Set correct boundary condition at right boundary')


# Defining required quantities
num_mesh = int(input("Enter the number of mesh points: "))
del_x = 1.0 / (num_mesh - 1.0)
alpha = -2.0 - del_x ** 2
bc = "Dirichlet"
norm = 0.0

x_mesh = np.zeros(num_mesh)
exact_sol = np.zeros(num_mesh)
numerical_sol = np.zeros(num_mesh)
rhs = np.zeros(num_mesh)
norm_val = np.zeros(num_mesh)

main_dia = np.zeros(num_mesh)
lower_dia = np.zeros(num_mesh)
upper_dia = np.zeros(num_mesh)

# Inserting values in diagonals
for i in range(0, num_mesh):
    x_mesh[i] = i * del_x
    exact_sol[i] = m.sin(m.pi * x_mesh[i])
    rhs[i] = -(del_x ** 2) * (m.pi ** 2 + 1) * m.sin(m.pi * x_mesh[i])
apply_bc(exact_sol)
apply_bc(numerical_sol)
apply_bc(rhs)


# Applying the boundary conditions
def apply_bc_dia():
    lower_dia[0] = 0.0
    upper_dia[0] = 0.0
    main_dia[0] = 1.0
    lower_dia[num_mesh - 1] = 0.0
    upper_dia[num_mesh - 1] = 0.0
    main_dia[num_mesh - 1] = 1.0


apply_bc_dia()

for i in range(1, len(numerical_sol) - 1):
    lower_dia[i] = 1
    upper_dia[i] = 1
    main_dia[i] = alpha

# Calling the function
tdma(num_mesh, lower_dia, upper_dia, main_dia, rhs, numerical_sol)
apply_bc(numerical_sol)

# Plotting the results
plt.plot(x_mesh, numerical_sol, 'g^')
plt.plot(x_mesh, exact_sol, 'r')
plt.title('Computational solution and Exact solution for %d grid points using TDMA' % num_mesh)
plt.xlim(0, 1)
plt.ylim(0, )
plt.xlabel('$\Delta x$', fontsize=14)
plt.ylabel('Solution Value', fontsize=14)
plt.show()
