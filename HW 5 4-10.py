import numpy as np
import matplotlib.pyplot as plt


N = 5  #  four segments
x = np.linspace(0, 4, N)  # This creates an array [0, 1, 2, 3, 4]
dx = x[1] - x[0]  # This should be 1

# Initialize the matrix A and vector b
A = np.zeros((N-2, N-2))  # We have N-2 unknowns because f0 and fN are known (boundary conditions)
b = np.zeros(N-2)

# Construct the matrix A using the given differential equation
for i in range(1, N-1):
    A[i-1, i-1] = -2/dx**2 + x[i]**2
    if i-1 > 0:
        A[i-1, i-2] = 1/dx**2
    if i+1 < N-1:
        A[i-1, i] = 1/dx**2
    b[i-1] = 2*x[i]*(4 - x[i])


# Since f(0) = f(4) = 0 and these values are at the edges of the matrix, they are not included in A or b.

# Solve the system
f_interior = np.linalg.solve(A, b)

# Including the boundary conditions in the full solution
f = np.zeros_like(x)
f[1:-1] = f_interior

plt.plot(x, f, 'o-', label='Numerical solution')
plt.xlabel('x')
plt.ylabel('f(x)')
plt.title('solution of the differential equation')
plt.legend()
plt.grid(True)
plt.show()
print(" the f(x) solution at each node is:")
print(f)