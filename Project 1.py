import numpy as np
import matplotlib.pyplot as plt
import math
So = 1e8
a = 10
N1 = 5
N2 = N1
x1 = np.linspace(0,a/2,N1)
x2 = np.linspace(a/2, 10*a, N2)
Sigma_a1_tr = 1.79e-2
Sigma_a1_abs = 8.08e-3
D1 = 1 / (3*(Sigma_a1_abs + Sigma_a1_tr ))
L1 = math.sqrt(D1 / Sigma_a1_abs)
Sigma_a2_tr = 8.77e-6
Sigma_a2_abs = 3.41e-2
D2 = 1 / (3*(Sigma_a2_abs + Sigma_a2_tr ))
L2 = math.sqrt(D2 / Sigma_a2_abs)
Anum = (D1 * L2 * math.cosh(a/(2*L1))) + (D2 * L1 * math.sinh(a/(2*L1)))
Adenom =  (D2 * L1 * math.cosh(a/(2*L1))) + (D1 * L2 * math.sinh(a/(2*L1)))
A = ((So * L1) / (2*D1)) * (Anum/Adenom)
B = (-So * L1) / (2 * D1)
Cnum = np.exp(a/(2*L2))
Cdenom = (D2 * L1 *math.cosh(a/(2*L1))) + (D1*L2*math.sinh(a/(2*L1)))
C = ((So*L1*L2) / 2 ) * (Cnum/Cdenom)
fanalytical_1 = (A* np.cosh( x1 / L1)) + (B* np.sinh( x1 / L1))
fanalytical_2 = C*np.exp(-x2/L2)
plt.title("Analytical Solution for a = %i cm" %(a))
plt.xlabel("Distance From Center")
plt.ylabel('Flux $(n/cm^2-s)$' )
plt.plot(x1, fanalytical_1)
plt.plot(x2, fanalytical_2)
plt.show()
analytical_flux = np.concatenate((fanalytical_1,fanalytical_2))
###
N = N1+N2+1
dx1 = ((a/2) / (N1))
dx2 = ((10*a/2) - (a/2)) / (N2)
x_nodes = np.concatenate((x1,x2))
b_0 = So/(2*dx1)
B = np.zeros(N-1)
B[0] = b_0
A = np.zeros((N-1,N-1))
A[0,0] = (D1/dx1**2) + (Sigma_a1_abs /2)
A[0,1] = -(D1/(dx1**2))
for i in range(1,N1-1):
    A[i,i-1] =  -(D1/(dx1**2))
for j in range(1,N1-1):
    A[j,j] = (Sigma_a1_abs + (2*D1/(dx1**2)))
for k in range(1,N1-1):
    A[k,k+1] =  -(D1/(dx1**2))
A[N1-1,N1-2] = -(D1/(dx1))
A[N1-1,N1-1] = (D1/dx1) + (D2/dx2) + ((dx1*Sigma_a1_abs)/2) +((dx2*Sigma_a2_abs)/2)
A[N1-1,N1] = -(D2/(dx2))
for i in range(N1,N-1):
    A[i,i-1] =  -(D2/(dx2**2))
for j in range(N1,N-1):
    A[j,j] = (Sigma_a2_abs + (2*D2/(dx2**2)))
for k in range(N1,N-2):
    A[k,k+1] =  -(D2/(dx2**2))
a_m = np.linalg.inv(A)
phi = a_m*B
numerical_flux = phi[:,0]
diff = abs(analytical_flux - numerical_flux)
fig, (ax1, ax2) = plt.subplots(2)
fig.suptitle('Neutron Flux Distribution (a=%i cm / N1=N2=%i)' %(a, N1))
ax1.plot(x_nodes,analytical_flux,marker='.', lw = 1,c ='g',label = 'Analytical')
ax1.plot(x_nodes,numerical_flux,marker='.', lw = 1,c ='b',label = 'Numerical')
ax1.set(xlabel=' ', ylabel='Flux $(n/cm^2-s)$')
ax1.legend()

ax2.plot(x_nodes,diff,marker='.', lw = 1,c ='k',label = 'Analytical')
ax2.set(xlabel='Distance From Center (cm)', ylabel='Difference')
plt.show()