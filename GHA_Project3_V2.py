import numpy as np
import matplotlib.pyplot as plt
import math

Na = 50
Nb = 50
N = Na+ Nb
# Group 1 of 2 terms(fast)
nusigmaf12A = 0.008476
nusigmaf12B = 0
SigmaA12 = 0.01207
SigmaB12 = 0.0494
Da12 = 1.2627
Db12 = 1.13
Removala12 = 0.02619
Removalb12 = 0.0494
SigmaS12A = Removala12 - SigmaA12
SigmaS12B = 0.0494 - 0.0004
Sigma_f1 = 0.003320

# Group 2 of 2 terms(thermal)
nusigmaf22A = 0.18514
nusigmaf22B = 0
SigmaA22 = 0.1210
SigmaB22 = 0.0197
Removal22 = 0.1210
Da22 = 0.3543
Db22 = 0.16
Sigma_f2 = 0.07537

 # Total width
Wa = 47.35 # Width of material A
Wb = 27.65# Width of material B
W = Wa+Wb
Delta_Xa = Wa / Na
Delta_Xb = Wb / Nb
matrix_A1 = np.zeros((N, N))
matrix_A2 = np.zeros((N, N))
matrix_B1 = np.zeros((N, N))
matrix_B2 = np.zeros((N, N))
matrix_B3 = np.zeros((N, N))
#(loop for filling out matrix A1)
for i in range(N):
    if i < Na:
        matrix_A1[i, i] = (Removala12 * Delta_Xa) + ((2 * Da12) / Delta_Xa)
        matrix_A1[i + 1, i] = -Da12 / Delta_Xa
        matrix_A1[i, i + 1] = -Da12 / Delta_Xa
    elif i == Na:
        matrix_A1[i, i] = (Da12 / Delta_Xa) + (Removala12*(Delta_Xa/2))+ (Db12 / Delta_Xb) + (Removalb12*(Delta_Xb/2))
    else:
        matrix_A1[i, i] = (Removalb12 * Delta_Xb) + ((2 * Db12) / Delta_Xb)#matrix_A[i, i] = (Sigma_B * Delta_Xb) + ((2 * Db) / Delta_Xb)
        matrix_A1[i - 1, i] = -Db12 / Delta_Xb
        matrix_A1[i, i - 1] = -Db12 / Delta_Xb
# (loop for filling out matrix A2)
for b in range(N):
    if b < Na:
        matrix_A2[b, b] = (SigmaA22 * Delta_Xa) + ((2 * Da22) / Delta_Xa)
        matrix_A2[b + 1, b] = -Da22 / Delta_Xa
        matrix_A2[b, b + 1] = -Da22 / Delta_Xa
    elif b == Na:
        matrix_A2[b, b] = (Da22 / Delta_Xa) + (Db22 / Delta_Xb) + ((SigmaA22 * Delta_Xa) / 2) + ((SigmaB22 * Delta_Xb) / 2)
    else:
        matrix_A2[b, b] = (SigmaB22 * Delta_Xb) + ((2 * Db22) / Delta_Xb)#matrix_A[i, i] = (Sigma_B * Delta_Xb) + ((2 * Db) / Delta_Xb)
        matrix_A2[b - 1, b] = -Db22 / Delta_Xb
        matrix_A2[b, b - 1] = -Db22 / Delta_Xb
#(Loop for filling out Matrix B1)
for j in range(N):
    if j < Na:
        matrix_B1[j, j] =Delta_Xa*(nusigmaf12A)
    elif j == Na:
        matrix_B1[j, j] =((Delta_Xa*(nusigmaf12A)/2)) + ((Delta_Xb*(nusigmaf12B))/2)
    else:
        matrix_B1[j, j] = Delta_Xb*(nusigmaf12B)#matrix_A[i, i] = (Sigma_B * Delta_Xb) + ((2 * Db) / Delta_Xb)
#(Loop for filling out Matrix B2)
for k in range(N):
    if k  < Na:
        matrix_B2[k , k ] = Delta_Xa*nusigmaf22A
    elif k  == Na:
        matrix_B2[k , k ] =((Delta_Xa*nusigmaf22A)/2) + ((Delta_Xb*nusigmaf22B)/2)
    else:
        matrix_B2[k , k ] = Delta_Xb*nusigmaf22B#matrix_A[i, i] = (Sig
#(Loop for filling out Matrix B3)
for l in range(N):
    if l < Na:
        matrix_B3[l, l] = Delta_Xa * SigmaS12A
    elif l == Na:
        matrix_B3[l, l] =((Delta_Xa* SigmaS12A)/2) + ((Delta_Xb* SigmaS12B)/2)
    else:
        matrix_B3[l, l] = Delta_Xb* SigmaS12B #matrix_A[i, i] = (Sigma_B * Delta_Xb) + ((2 * Db) / Delta_Xb)

# Section for calculation of K-effictive adn curve

flux1 = np.ones(N)
flux2 = np.ones(N)
k = 1.00
Source1 = matrix_B1.dot(flux1) + matrix_B2.dot(flux2)
Source2 = matrix_B3.dot(flux1)

flxdiff1 = 1.0
flxdiff2 = 1.0
g = 1

while flxdiff1 > 0.0001 and flxdiff2 > 0.0001:
    oldflux1 = flux1.copy()
    oldflux2 = flux2.copy()
    oldk = k
    oldsource1 = Source1.copy()
    oldsource2 = Source2.copy()

    flux1 = np.linalg.inv(matrix_A1).dot(oldsource1) / oldk
    flxdiff1 = np.sqrt(np.sum((flux1 - oldflux1) ** 2))


    flux2 = np.linalg.inv(matrix_A2).dot(oldsource2) / oldk
    flxdiff2 = np.sqrt(np.sum((flux2 - oldflux2) ** 2))
    Source1 = matrix_B1.dot(flux1) + matrix_B2.dot(flux2)
    Source2 = matrix_B3.dot(flux1)
    k = oldk * np.sum(Source1) / np.sum(oldsource1)
    g += 1

    # Normalize the flux

Pth = 3400 * (10 ** 6) * 1.602e19  # Convert MWth to eV/s
kappa = 200 * (10 ** 6)  # MeV per fission
dV =  W/2
S = np.sum(Sigma_f1 * flux1 + Sigma_f2 * flux2)
C = Pth / (kappa * S * np.sum(dV))

normflux1 = C * flux1
normflux1 = [0,*normflux1,0]
normflux2 = C * flux2
normflux2 = [0,*normflux2,0]

print(g)
xa_nodes = np.linspace(0, Wa, Na+2)
xb_nodes = np.linspace(Wa, W, Nb+2)
x_nodes = np.append(xa_nodes,xb_nodes[1:])
plt.figure(figsize=(12, 6))
plt.axvline(x=Wa)
plt.plot(x_nodes[:-1] + np.diff(x_nodes) / 2, normflux1, ':+', label='Fast')
plt.plot(x_nodes[:-1] + np.diff(x_nodes) / 2, normflux2, ':+', label='Thermal')
plt.xlabel('Location, x')
plt.ylabel('Normalized Flux')
plt.title(f'Normalized Flux Distribution (Na={Na}, Nb={Nb}, k={k:.5f}, Wa ={Wa},Wb ={Wb})')
plt.legend()
plt.show()