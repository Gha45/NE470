#1.Consider the differential equation analytical and numerical solution presented in class,
# in Lecture Packet 04c “Simple Discretization Example.”
#a.Reproduce both of these solutions on your own, use MATLAB or Python.
# [Graduate or honor students use FORTRAN, C++, Python, etc.].
# b.Modify the program to use a variable number of nodes “N” to solve
# the numerical problem more accurately
#c.Explore how many nodes you need for the solutions to agree within about 0.01%.
# part a: Reproduce 04c solutions
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
a = np.array([[-3, 1, 0], [1, -3, 1], [0, 1, -3]])
b = np.array([[-2], [0], [-54.61647]])
f = np.linalg.solve(a, b)

soln = np.vstack((-b[0], f, -b[2]))

x = np.linspace(0, 4, 5)

analytical = (np.exp(x) + np.exp(-x)).reshape(-1, 1)

plt.plot(x, soln, ':+', label='Numerical')
plt.plot(x, analytical, '-o', label='Analytical')

plt.xlabel('x')
plt.ylabel('f(x)')
plt.legend(loc=2, title='Legend', title_fontsize='13', fontsize='10')
plt.title('Part A')
plt.show()

print("part a):")
print("Analytical:")
print(analytical)
print("Numerical:")
print(soln)
# part b
analytical = (np.exp(x) + np.exp(-x)).reshape(-1, 1)
def solve_numerical(N, x0, xN, f0, fN):
    dx = (xN - x0) / (N - 1)  # Step size

    # Initialize matrix A
    A = np.zeros((N - 2, N - 2))

    # Fill the main diagonal
    for i in range(N - 2):
        A[i, i] = -(2.0 / dx ** 2) - 1.0

    # Fill the off diagonals
    for i in range(1, N - 2):
        A[i, i - 1] = 1.0 / dx ** 2
        A[i - 1, i] = 1.0 / dx ** 2

    # Initialize vector b
    b = np.zeros(N - 2)
    b[0] = -f0 / dx ** 2
    b[-1] = -fN / dx ** 2

    # Solve the system r
    numsol = np.linalg.solve(A, b)

    return numsol

N = 5  # Number of nodes including boundary points
x0 = 0.0  # Start of the domain
xN = 4.0  # End of the domain
f0 = 2.0  # Boundary condition at x0
fN = 54.61647  # Boundary condition at xN

numsol = solve_numerical(N, x0, xN, f0, fN)

# Generate the full solution including the boundary points
full_sol = np.hstack((f0, numsol, fN))

print("part b):")
print("Analytical:")
print(analytical)
print("Numerical:")
print(full_sol)

# Generate the x values including the boundary points
x_values = np.linspace(x0, xN, N)
plt.plot(x, analytical, '-o', label='Analytical')
plt.plot(x_values, full_sol, '-o', label='Numerical Solution with Variable Nodes')
plt.xlabel('x')
plt.ylabel('f(x)')
plt.title('Part B')
plt.legend()
plt.show()

#Part C
# Iterating over possible N values to find the minimum N that satisfies the condition
analytical = (np.exp(x) + np.exp(-x)).reshape(-1, 1)

def solve_numerical(N, x0, xN, f0, fN):
    dx = (xN - x0) / (N - 1)  # Step size
    x = np.linspace(x0, xN, N)

    # Initialize matrix A
    A = np.zeros((N - 2, N - 2))

    # Fill the main diagonal
    for i in range(N - 2):
        A[i, i] = -(2.0 / dx ** 2) - 1.0

    # Fill the off diagonals
    for i in range(1, N - 2):
        A[i, i - 1] = 1.0 / dx ** 2
        A[i - 1, i] = 1.0 / dx ** 2

    # Initialize vector b
    b = np.zeros(N - 2)
    b[0] = -f0 / dx ** 2
    b[-1] = -fN / dx ** 2

    # Solve the system
    numsol = np.linalg.solve(A, b)

    return np.hstack((f0, numsol, fN)), x

N = 5  # Number of nodes including boundary points
x0 = 0.0  # Start of the domain
xN = 4.0  # End of the domain
f0 = 2.0  # Boundary condition at x0
fN = 54.61647  # Boundary condition at xN

max_error = 1.0
x_high_res = np.linspace(x0, xN, 1000)
analytical_high_res = (np.exp(x_high_res) + np.exp(-x_high_res))

# Increase N until the maximum error is below the threshold
while max_error > 0.01:
    full_sol, x_values = solve_numerical(N, x0, xN, f0, fN)
    analytical_at_nodes = np.interp(x_values, x_high_res, analytical_high_res)
    error = np.abs(full_sol - analytical_at_nodes)
    max_error = np.max(error)
    if max_error > 0.01:
        N += 1  # Increase N for the next iteration

# Plotting the results
plt.figure(figsize=(12, 6))

plt.subplot(1, 2, 1)
plt.plot(x_high_res, analytical_high_res, '-o', label='Analytical')
plt.plot(x_values, full_sol, '-x', label=f'Numerical Solution (N={N})')
plt.xlabel('x')
plt.ylabel('f(x)')
plt.legend()
plt.title('Part C')
plt.subplot(1, 2, 2)
plt.plot(x_values, error, '-x', label='Error')
plt.xlabel('x')
plt.ylabel('Error')
plt.title(f'Max Error: {max_error:.2f}')
plt.legend()

plt.tight_layout()
plt.show()