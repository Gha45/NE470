from scipy.optimize import fsolve
K_inf_target = 1
neutrons = 2.43 #neutrons per fission
micro_fission_U235 = 524*(10**-24)
Micro_abs_U235 =683*(10**-24)
micro_abs_graph = 0.004 * (10**-24)
micro_abs_ber = 0.01 * (10**-24)
micro_abs_water= 0.66 * (10**-24)
micro_abs_heavy = 0.001 * (10**-24)

def NmodNfuel(neutrons, micro_fission_fuel, micro_abs_mod, micro_abs_fuel):
    return (neutrons * micro_fission_fuel - micro_abs_fuel) / micro_abs_mod

print(NmodNfuel(neutrons,micro_fission_U235,micro_abs_graph,Micro_abs_U235))
micro_abs_mod_values = {
    "Graphite": NmodNfuel(neutrons,micro_fission_U235,micro_abs_graph,Micro_abs_U235),
    "Beryllium": NmodNfuel(neutrons,micro_fission_U235,micro_abs_ber,Micro_abs_U235),
    "Water":NmodNfuel(neutrons,micro_fission_U235,micro_abs_water,Micro_abs_U235),
    "Heavy Water": NmodNfuel(neutrons,micro_fission_U235,micro_abs_heavy,Micro_abs_U235)}
for modulator,ratio in micro_abs_mod_values.items():
    print(f"The ration of moderator to feul density for {modulator} is:{ratio}")
