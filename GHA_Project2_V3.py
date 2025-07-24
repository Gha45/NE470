import numpy as np
import matplotlib.pyplot as plt
import math
# Deliverable 5&6
def analytical_solution(x, W):
    return np.cos((math.pi * x / W)-(math.pi/2))

def calculate_flux(Na, Nb):
    # Constants and material properties
    Sigma_Tr = 0.0362
    Sigma_A = 0.1532
    Sigma_B = 0.1532
    nusigfA = 0.1570
    nusigfB = 0.1570
    W = 67.6  # Total width
    Wa = 33.8  # Width of material A
    Wb = 33.8  # Width of material B
    Da = 1 / (3 * (Sigma_A + Sigma_Tr))
    Db = 1 / (3 * (Sigma_B + Sigma_Tr))
    xcrit = np.pi / np.sqrt((nusigfA - Sigma_A) / Da)
    print('X critical dimension:' + f'{xcrit}')
    N = Na + Nb
    Delta_Xa = Wa / Na
    Delta_Xb = Wb / Nb

    # Matrix A setup
    matrix_A = np.zeros((N, N))
    for i in range(N):
      if i< Na:
          matrix_A[i,i] = (Sigma_A * Delta_Xa) +((2*Da)/Delta_Xa)
          matrix_A[i+1,i] = -Da/Delta_Xa
          matrix_A[i,i+1] = -Da/Delta_Xa
      elif i == Na :
          matrix_A[i,i]= (Da/Delta_Xa) + (Db/Delta_Xb) +((Sigma_A * Delta_Xa)/2)+((Sigma_B *Delta_Xb)/2)
      else:
        matrix_A[i, i] = (Sigma_B * Delta_Xb) + ((2 * Db) / Delta_Xb)
        matrix_A[i - 1, i] = -Db / Delta_Xb
        matrix_A[i, i - 1] = -Db / Delta_Xb

    matrix_B = np.zeros((N, N))
    for i in range(N):
        if i < Na:
            matrix_B[i,i] = (nusigfA* Delta_Xa)
        elif i == Na:
            matrix_B[i,i] = (nusigfA*(Delta_Xa/2))+(nusigfB*(Delta_Xb/2))
        else:
            matrix_B[i,i] = nusigfB * Delta_Xb

    # Solve the eigenvalue problem
    flux = np.ones(N)
    k = 1.00
    S = matrix_B.dot(flux)
    flxdiff = 1.0
    iter_count = 0
    while flxdiff > 0.0001:
        iter_count += 1
        oldflux = flux.copy()
        oldk = k
        oldS = S.copy()

        flux = np.linalg.inv(matrix_A).dot(oldS) / oldk
        S = matrix_B.dot(flux)
        k = oldk * (sum(S) / sum(oldS))

        flxdiff = np.sqrt(np.sum((flux - oldflux) ** 2))

    # Normalize the flux
    total = np.max(flux)
    normflux = flux / total
    normflux = [0,*normflux, 0]
    print(normflux)
    return normflux, k, iter_count, W
def plot_flux(Na, Nb):
    normflux, k, iter_count, W = calculate_flux(Na, Nb)
    x_nodes = np.linspace(0, W, Na + Nb + 3)
    print(x_nodes)
    plt.figure(figsize=(12, 6))
    plt.plot(x_nodes[:-1] + np.diff(x_nodes) / 2, normflux, ':+', label='Numerical')

    # Analytical solution plotting
    analytical_flux = analytical_solution(x_nodes, W)
    total_analytical = np.max(analytical_flux)
    normalized_analytical_flux = analytical_flux / total_analytical

    plt.plot(x_nodes, normalized_analytical_flux, '--', label='Analytical')

    plt.xlabel('Location, x')
    plt.ylabel('Normalized Flux')
    plt.title(f'Normalized Flux Distribution (Na={Na}, Nb={Nb}, k={k:.5f}, Iterations={iter_count})')
    plt.legend()
    plt.show()

# Example configurations
configurations = [(10, 10), (25, 25), (5, 5),(50,50)]
for Na, Nb in configurations:
    plot_flux(Na,Nb)


def calculate_fluxB(Na,Nb,Wa,Wb):
    # Constants and material properties
    Sigma_TrA = 0.0362
    Sigma_TrB = 1.79*10**(-2)
    Sigma_A = 0.1532
    Sigma_B = 8.08*10**(-3)
    nusigfA = 0.1570
    nusigfB = 0
    W = Wa+Wb
    #W = 67.6  # Total width
    #Wa = 33.8  # Width of material A
    #Wb = 33.8  # Width of material B
    Da = 1 / (3 * (Sigma_A + Sigma_TrA))
    Db = 1 / (3 * (Sigma_B + Sigma_TrB))
    N = Na + Nb
    Delta_Xa = Wa / Na
    Delta_Xb = Wb / Nb
    xcritA = np.pi / np.sqrt((nusigfA - Sigma_A) / Da)
    print('X critical dimension for region A:'+ f'{xcritA}')
    xcritB = np.pi / np.sqrt((nusigfB - Sigma_B) / Db)
    print('X critical dimension for region B:'+f'{xcritB}')
    # Matrix A setup
    matrix_A = np.zeros((N, N))
    for i in range(N):
      if i< Na:
          matrix_A[i,i] = (Sigma_A * Delta_Xa) +((2*Da)/Delta_Xa)
          matrix_A[i+1,i] = -Da/Delta_Xa
          matrix_A[i,i+1] = -Da/Delta_Xa
      elif i == Na :
          matrix_A[i,i]= (Da/Delta_Xa) + (Db/Delta_Xb) +((Sigma_A * Delta_Xa)/2)+((Sigma_B *Delta_Xb)/2)
      else:
        matrix_A[i, i] = (Sigma_B * Delta_Xb) + ((2 * Db) / Delta_Xb)
        matrix_A[i - 1, i] = -Db / Delta_Xb
        matrix_A[i, i - 1] = -Db / Delta_Xb

    matrix_B = np.zeros((N, N))
    for i in range(N):
        if i < Na:
            matrix_B[i,i] = (nusigfA* Delta_Xa)
        elif i == Na:
            matrix_B[i,i] = (nusigfA*(Delta_Xa/2))+(nusigfB*(Delta_Xb/2))
        else:
            matrix_B[i,i] = nusigfB * Delta_Xb

    # Solve the eigenvalue problem
    flux = np.ones(N)
    k = 1.00
    S = matrix_B.dot(flux)
    flxdiff = 1.0
    iter_count = 0
    while flxdiff > 0.0001:
        iter_count += 1
        oldflux = flux.copy()
        oldk = k
        oldS = S.copy()

        flux = np.linalg.inv(matrix_A).dot(oldS) / oldk
        S = matrix_B.dot(flux)
        k = oldk * (sum(S) / sum(oldS))

        flxdiff = np.sqrt(np.sum((flux - oldflux) ** 2))

    # Normalize the flux
    total = np.max(flux)
    normflux = flux / total
    normflux = [0,*normflux,0]

    return normflux, k, iter_count, W
def plot_fluxb(Na, Nb,Wa,Wb):

    normflux, k, iter_count, W = calculate_fluxB(Na, Nb,Wa,Wb)
    x_nodes = np.linspace(0, W, Na + Nb + 3)
    print(W)
    plt.figure(figsize=(12, 6))
    plt.plot(x_nodes[:-1] + np.diff(x_nodes) / 2, normflux, ':+', label='Numerical')

    # Analytical solution plotting
    analytical_flux = analytical_solution(x_nodes, W)
    total_analytical = np.max(analytical_flux)
    normalized_analytical_flux = analytical_flux / total_analytical

   # plt.plot(x_nodes, normalized_analytical_flux, '--', label='Analytical')

    plt.xlabel('Location, x')
    plt.ylabel('Normalized Flux')
    plt.title(f'Normalized Flux Distribution (Wa={Wa},Wb={Wb},Na={Na}, Nb={Nb}, k={k:.5f}, Iterations={iter_count})')
    plt.legend()
    plt.show()

#configurations2 = [(10, 25,10,5),(10, 10,20,20), (25, 25,33.8,33.8), (50, 50,67.6,300)]
#for Na, Nb,Wa,Wb in configurations2:
#    plot_fluxb(Na, Nb,Wa,Wb)