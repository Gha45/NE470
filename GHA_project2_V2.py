import numpy as np
import matplotlib.pyplot as plt
import math


def setup_matrices(Na, Nb, Sigma_A, Sigma_B, Sigma_Tr, nusigfA, nusigfB, Wa, Wb):
    N = Na + Nb
    Delta_Xa = Wa / Na
    Delta_Xb = Wb / Nb
    Da = 1 / (3 * (Sigma_A + Sigma_Tr))
    Db = 1 / (3 * (Sigma_B + Sigma_Tr))
    xcrit = np.pi / np.sqrt((nusigfA - Sigma_A) / Da)
    print('X critical dimension:'+xcrit)
    matrix_A = np.zeros((N, N))
    matrix_B = np.zeros((N, N))

    for i in range(N):
        if i < Na - 1:
            matrix_A[i, i] = (2 * (Da / Delta_Xa)) + (Sigma_A * Delta_Xa)
            matrix_A[i, i + 1] = -Da / Delta_Xa
            matrix_A[i, i - 1] = -Da / Delta_Xa
        elif i == Na - 1:
            matrix_A[i, i - 1] = -Da / Delta_Xa
            matrix_A[i, i] = (Da / Delta_Xa) + (Db / Delta_Xb) + ((Sigma_A * Delta_Xa) / 2) + ((Sigma_B * Delta_Xb) / 2)
            matrix_A[i, i + 1] = -Db / Delta_Xb
        elif i < N - 1:
            matrix_A[i, i] = (2 * (Db / Delta_Xb)) + (Sigma_B * Delta_Xb)
            matrix_A[i, i + 1] = -Db / Delta_Xb
            matrix_A[i, i - 1] = -Db / Delta_Xb
        else:
            matrix_A[i, i] = (2 * (Db / Delta_Xb)) + (Sigma_B * Delta_Xb)
            matrix_A[i, i - 1] = -Db / Delta_Xb

        # Fill matrix_B with source terms
        if i < Na - 1:
            matrix_B[i, i] = nusigfA * Delta_Xa
        elif i == Na - 1:
            matrix_B[i, i] = (nusigfA * (Delta_Xa / 2)) + (nusigfB * (Delta_Xb / 2))
        else:
            matrix_B[i, i] = nusigfB * Delta_Xb

    return matrix_A, matrix_B, Delta_Xa, Delta_Xb


def analytical_solution(x, W):
    return np.cos((math.pi * x / W) - (math.pi / 2))


def compute_error(num_flux, x_nodes, W):
    analytical_flux = analytical_solution(x_nodes, W)
    analytical_flux /= sum(analytical_flux)
    num_flux /= sum(num_flux)
    return np.linalg.norm(analytical_flux - num_flux) / np.linalg.norm(analytical_flux)


# Initialize parameters
Sigma_Tr = 0.0362
Sigma_A = 0.1532
Sigma_B = 0.1532
nusigfA = 0.1570
nusigfB = 0.1570
W = 67.6
Wa = 33.8
Wb = 33.8
Na = 10  # Start with a guess
Nb = 10
flux = None
error_threshold = 0.01  # 1% error
max_iter = 1000
iter = 0

while True:
    matrix_A, matrix_B, Delta_Xa, Delta_Xb = setup_matrices(Na, Nb, Sigma_A, Sigma_B, Sigma_Tr, nusigfA, nusigfB, Wa,
                                                            Wb)
    flux = np.ones(Na + Nb)
    k = 1.00
    S = matrix_B.dot(flux)

    # Iteration for convergence
    while True:
        old_flux = flux.copy()
        old_k = k
        old_S = S.copy()

        flux = np.linalg.inv(matrix_A).dot(old_S) / old_k
        S = matrix_B.dot(flux)
        k = old_k * (sum(S) / sum(old_S))
        flxdiff = np.sqrt(np.sum((flux - old_flux) ** 2))

        if flxdiff < 0.0001:
            break

    x_nodes = np.linspace(0, W, Na + Nb + 1)
    error = compute_error(flux, x_nodes[:-1] + np.diff(x_nodes) / 2, W)

    if error < error_threshold:
        break
    else:
        Na += 2  # Increase the mesh density
        Nb += 2  # Keep Na and Nb equal

    if Na + Nb > max_iter:
        print("Max iterations exceeded without reaching desired accuracy")
        break

# Plot final results
plt.figure(figsize=(12, 6))
plt.plot(x_nodes[:-1] + np.diff(x_nodes) / 2, flux / sum(flux), ':+', label='Numerical Solution')
plt.plot(x_nodes, analytical_solution(x_nodes, W) / sum(analytical_solution(x_nodes, W)), label='Analytical Solution')
plt.xlabel('Location, x')
plt.ylabel('Flux')
plt.title('Normalized Flux Comparison')
plt.text(W / 4, max(flux / sum(flux)) / 2, 'k-eff = ')
plt.text(W / 3, max(flux / sum(flux)) / 2, f'{k}')
plt.legend()
plt.show()

print(f"Final k-eff: {k}")
print(f"Mesh Points: Na = {Na}, Nb = {Nb}")
