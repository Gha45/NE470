import numpy as np
import matplotlib.pyplot as plt
neutron_speed = 2.43
micro_fission_U235 =524*(10**-24)
Micro_absorption_U235 =683*(10**-24)
Micro_absorption_U238 = 7.58*(10**-24)
weights = np.arange(0.007, 1, 0.001)
k_inf_values = []
for weight in weights:

    N_U235=((weight)*(9.255)*(6.022*10**23))/235
    N_U238 =((1-weight)*9.255*(6.022*10**23))/238

    Macro_fission_U235 = N_U235 * micro_fission_U235
    Macro_absorption_U235 = N_U235 * Micro_absorption_U235
    Macro_absorption_U238 = N_U238 * Micro_absorption_U238
    k_inf = neutron_speed*(Macro_fission_U235/(Macro_absorption_U235 + Macro_absorption_U238))
    k_inf_values.append(k_inf)
plt.figure(figsize=(10, 6))
plt.plot(weights, k_inf_values, marker='o')
plt.title('Infinite Multiplication Factor (k$_{inf}$) vs U-235 Weight Fraction')
plt.xlabel('Weight Fraction of U-235')
plt.ylabel('k$_{inf}$')
plt.grid(True)
plt.show()