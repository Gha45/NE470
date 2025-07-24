import numpy as np
import matplotlib.pyplot as plt
import math
# Properties from D&H Table 5-3
Sigma_Tr = 0.0362
Sigma_A = 0.1532
Sigma_B = 0.1532
nusigfA = 0.1570
nusigfB = 0.1570

# Total width and sub-widths
W = 67.6
Wa = 33.8
Wb = 33.8

# Calculation of diffusion coefficients and critical dimension
Da = 1 / (3 * (Sigma_A + Sigma_Tr))
Db = 1 / (3 * (Sigma_B + Sigma_Tr))
xcrit = np.pi / np.sqrt((nusigfA - Sigma_A) / Da)

Na = 26
Nb = 25
N = Na + Nb  # Total nodes
# Mesh sizes
Delta_Xa = Wa / Na
Delta_Xb = Wb / Na



# Initialize the matrix
matrix_A = np.zeros((N, N))

# Loop to fill the matrix
for i in range(N):
    if i == 0:  # First node in set A
        matrix_A[i, i] = (2 * (Da / Delta_Xa)) + (Sigma_A * Delta_Xa)
        matrix_A[i, i + 1] = -Da / Delta_Xa
    elif i < Na - 1:  # Remaining nodes within set A, excluding the interface node
        matrix_A[i, i] = (2 * (Da / Delta_Xa)) + (Sigma_A * Delta_Xa)
        matrix_A[i, i + 1] = -Da / Delta_Xa
        matrix_A[i, i - 1] = -Da / Delta_Xa
    elif i == Na - 1:  # Interface node
        matrix_A[i, i - 1] = -Da / Delta_Xa
        matrix_A[i, i] = (Da / Delta_Xa) + (Db / Delta_Xb) + ((Sigma_A * Delta_Xa) / 2) + ((Sigma_B * Delta_Xb) / 2)
        matrix_A[i, i + 1] = -Db / Delta_Xb
    elif i < N - 1:  # Nodes within set B, excluding the last node
        matrix_A[i, i] = (2 * (Db / Delta_Xb)) + (Sigma_B * Delta_Xb)
        matrix_A[i, i + 1] = -Db / Delta_Xb
        matrix_A[i, i - 1] = -Db / Delta_Xb
    else:  # Last node within set B
        matrix_A[i, i] = (2 * (Db / Delta_Xb)) + (Sigma_B * Delta_Xb)
        matrix_A[i, i - 1] = -Db / Delta_Xb
print(matrix_A)
matrix_B = np.zeros((N, N))

# Loop to fill matrix B
for i in range(N):
    if i < Na - 1:  # Nodes within set A, excluding the interface node
        matrix_B[i, i] = nusigfA * Delta_Xa
    elif i == Na - 1:  # Interface node
        matrix_B[i, i] = (nusigfA * (Delta_Xa / 2)) + (nusigfB * (Delta_Xb / 2))
    else:  # Nodes within set B
        matrix_B[i, i] = nusigfB * Delta_Xb

flux = np.ones(N)
k = 1.00
S = matrix_B.dot(flux)

flxdiff = 1.0
iter = 0

plt.figure(1)
plt.ion()

while flxdiff > 0.0001:
    iter += 1

    oldflux = flux.copy()
    oldk = k
    oldS = S.copy()

    flux = np.linalg.inv(matrix_A).dot(oldS) / oldk
    S = matrix_B.dot(flux)
    k = oldk * (sum(S) / sum(oldS))

    flxdiff = np.sqrt(np.sum((flux - oldflux) ** 2))

    plt.plot(flux, label=f'Iteration {iter}')

plt.xlabel('Node Number')
plt.ylabel('Flux')
plt.title('Iterative Un-normalized Flux')
plt.legend()
plt.ioff()
plt.show()

# Normalize flux
total = sum(flux)
normflux = flux / total
W = Delta_Xa * (Na - 1) + Delta_Xb * Nb

def analytical_solution(x, W):
    return np.cos((math.pi * x / W)-(math.pi/2))


x_nodes = np.linspace(0, W, N + 1)

plt.figure(2)
plt.plot(x_nodes[:-1] + np.diff(x_nodes) / 2, normflux, ':+')
analytical = analytical_solution(x_nodes, W)
plt.plot(x_nodes, analytical / sum(analytical))
plt.text(W / 4, max(normflux) / 2, 'k-eff = ')
plt.text(W / 3, max(normflux) / 2, f'{k}')
plt.xlabel('Location, x')
plt.ylabel('Flux')
plt.title('Normalized Flux')
plt.show()

