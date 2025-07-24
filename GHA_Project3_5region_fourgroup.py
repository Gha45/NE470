import numpy as np
import matplotlib.pyplot as plt
import math

### Values###

# Define material properties
material_properties = {
    "MOX": {
        "nusigmaF14": 0.01244,
        "nusigmaF24": 0.01891,
        "nusigmaF34": 0.06742,
        "nusigmaF44": 0.1989,
        'Chi14': 0.9953,
        'Chi24': 0.004683,
        'Chi34': 3.512E-7,
        'Chi44': 1.902E-9,
        "D14": 2.1623,
        "D24": 1.0867,
        "D34": 0.6318,
        "D44": 0.3543,
        "Removal14": 4.631e-2,
        "Removal24": 1.027e-1,
        "Removal34": 1.665e-1,
        "Removal44": 1.284e-1,
        "SigmaS11": 0.0,
        "SigmaS12": 4.121e-2,
        "SigmaS13": 9.151e-5,
        "SigmaS14": 4.688e-7,
        "SigmaS21": 0.0,
        "SigmaS22": 0.0,
        "SigmaS23": 8.485e-2,
        "SigmaS24": 4.49e-4,
        "SigmaS31": 0.0,
        "SigmaS32": 0.0,
        "SigmaS33": 0.0,
        "SigmaS34": 9.931e-2,
        "SigmaS41": 0.0,
        "SigmaS42": 0.0,
        "SigmaS43": 2.205e-3,
        "SigmaS44": 0.0

    },
    "Feul10": {
        "nusigmaF14": 4.4487e-3,
        "nusigmaF24": 3.971e-3,
        "nusigmaF34": 1.896e-2,
        "nusigmaF44": 1.324e-1,
        'Chi14': 0.9955E-1,
        'Chi24': 4.545E-3,
        'Chi34': 4.701E-7,
        'Chi44': 2.546E-9,
        "D14": 2.1623,
        "D24": 1.0867,
        "D34": 0.6318,
        "D44": 0.3543,
        "Removal14": 4.542e-2,
        "Removal24": 9.931E-2,
        "Removal34": 1.427e-1,
        "Removal44": 9.367e-2,
        "SigmaS11": 0.0,
        "SigmaS12": 4.286e-2,
        "SigmaS13": 9.518e-5,
        "SigmaS14": 4.882e-7,
        "SigmaS21": 0.0,
        "SigmaS22": 0.0,
        "SigmaS23": 8.87e-2,
        "SigmaS24": 4.69e-4,
        "SigmaS31": 0.0,
        "SigmaS32": 0.0,
        "SigmaS33": 0.0,
        "SigmaS34": 1.075e-1,
        "SigmaS41": 0.0,
        "SigmaS42": 0.0,
        "SigmaS43": 1.757e-3,
        "SigmaS44": 0.0
    },
"Reflector": {
        "nusigmaF14": 0.0,
        "nusigmaF24": 0.0,
        "nusigmaF34": 0.0,
        "nusigmaF44": 0.0,
        'Chi14': 0.0,
        'Chi24': 0.0,
        'Chi34': 0.0,
        'Chi44': 0.0,
        "D14": 1.13,
        "D24": 1.13,
        "D34": 1.13,
        "D44": 0.16,
        "Removal14": 3.096e-2,
        "Removal24": 7.707e-2,
        "Removal34": 1.033e-1,
        "Removal44": 9.351e-3,
        "SigmaS11": 0.0,
        "SigmaS12": 3.058e-2,
        "SigmaS13": 6.841e-5,
        "SigmaS14": 3.52e-7,
        "SigmaS21": 0.0,
        "SigmaS22": 0.0,
        "SigmaS23": 7.42e-2,
        "SigmaS24": 3.947e-4,
        "SigmaS31": 0.0,
        "SigmaS32": 0.0,
        "SigmaS33": 0.0,
        "SigmaS34": 1.018e-1,
        "SigmaS41": 0.0,
        "SigmaS42": 0.0,
        "SigmaS43": 2.373e-4,
        "SigmaS44": 0.0
    }
}


def get_input():
    num_regions = int(input("Enter the number of regions: "))
    regions = []
    for i in range(num_regions):
        length = float(input(f"Enter the length of region {i + 1}: "))
        nodes = int(input(f"Enter the number of nodes in region {i + 1}: "))
        material = input(f"Enter the material of region {i + 1}: ")
        regions.append((length, material, nodes))
    return regions,


# regions_info = get_input()

# print("Regions Info:", regions_info)
Z = 10
Y = 8
regions_info = [(Y, 'Reflector', Z),(Y, 'Feul10', Z) , (Y, 'MOX', Z),(Y, 'Feul10', Z),(Y, 'Reflector', Z)]
# regions_info = [(27.65, 'Reflector', 20),(18.43, 'MOX', 20) , (18.43, 'MOX', 20),(18.43, 'MOX', 20),(27.65, 'Reflector', 20)]
# regions_info = [(27.65, 'Reflector', 20),(18.43, 'FeulU', 20) , (18.43, 'FeulU', 20),(18.43, 'FeulU', 20),(7, 'Reflector', 20)]
#regions_info = [(23.675, 'Reflector', 20), (47.35, 'PWR', 60), (23.675, 'Reflector', 20)]
#regions_info = [(56.625, 'PWR', 40)]
# Define properties for each region
nodes_per_region = []  # Number of regions
materials = []
lengths = []
for region in regions_info:
    nodes_per_region.append(
        region[2])  # Number of nodes per region
    materials.append(region[1])
    lengths.append(region[0])

    R = len(regions_info)
import numpy as np

N = sum(region[2] for region in
        regions_info)  # Total number of nodes
matrix_A1 = np.zeros((N, N))
matrix_A2 = np.zeros((N, N))
matrix_A3 = np.zeros((N, N))
matrix_A4 = np.zeros((N, N))
matrix_B = {}
for i in ['1','2','3','4']:
    for j in ['1','2','3','4']:
        matrix_B[i+j]= np.zeros((N,N))

# Set up of non leakage proabability and diffusion lengths
L1 = []
L2 = []
PNL1 = []
PNL2 = []
# Matrix A1
start_index1 = 0
for idx, region in enumerate(regions_info):
    length, material, nodes = region
    properties = material_properties[material]

    Delta_X = length / nodes
    Removal = properties['Removal14']
    D = properties['D14']
    L1.append(np.sqrt(D / Removal))

    node_index_start = start_index1
    node_index_end = start_index1 + nodes - 1

    for i in range(nodes):
        node_index = start_index1 + i

        # Diagonal element for all nodes
        matrix_A1[node_index, node_index] = (Removal * Delta_X) + ((2 * D) / Delta_X)

        # Off-diagonal elements within the current region
        if i < nodes - 1:
            matrix_A1[node_index, node_index + 1] = -D / Delta_X
            matrix_A1[node_index + 1, node_index] = -D / Delta_X

            # Handling the boundary between two regions
    if idx < len(regions_info) - 1:
        next_region = regions_info[idx + 1]
        next_length, next_material, next_nodes = next_region
        next_properties = material_properties[next_material]
        next_D = next_properties['D14']
        next_Removal = next_properties['Removal14']
        next_Delta_X = next_length / next_nodes
        print(next_material)

        # Modify the last node of current region considering the interface
        matrix_A1[node_index_end, node_index_end] /= 2
        matrix_A1[node_index_end, node_index_end] += (next_D / next_Delta_X) + (next_Removal * next_Delta_X / 2)
        if node_index_end + 1 < N:
            matrix_A1[node_index_end, node_index_end + 1] = -next_D / next_Delta_X
            matrix_A1[node_index_end + 1, node_index_end] = -next_D / next_Delta_X

    start_index1 += nodes

# Matrix A2
start_index1 = 0
for idx, region in enumerate(regions_info):
    length, material, nodes = region
    properties = material_properties[material]
    Delta_X = length / nodes
    Removal = properties['Removal24']
    D = properties['D24']

    node_index_start = start_index1
    node_index_end = start_index1 + nodes - 1

    for i in range(nodes):
        node_index = start_index1 + i

        # Diagonal element for all nodes
        matrix_A2[node_index, node_index] = (Removal * Delta_X) + ((2 * D) / Delta_X)

        # Off-diagonal elements within the current region
        if i < nodes - 1:
            matrix_A2[node_index, node_index + 1] = -D / Delta_X
            matrix_A2[node_index + 1, node_index] = -D / Delta_X

            # Handling the boundary between two regions
    if idx < len(regions_info) - 1:
        next_region = regions_info[idx + 1]
        next_length, next_material, next_nodes = next_region
        next_properties = material_properties[next_material]
        next_D = next_properties['D24']
        next_Removal = next_properties['Removal24']
        next_Delta_X = next_length / next_nodes

        # Modify the last node of current region considering the interface
        matrix_A2[node_index_end, node_index_end] /= 2
        matrix_A2[node_index_end, node_index_end] += (next_D / next_Delta_X) + (next_Removal * next_Delta_X / 2)
        if node_index_end + 1 < N:
            matrix_A2[node_index_end, node_index_end + 1] = -next_D / next_Delta_X
            matrix_A2[node_index_end + 1, node_index_end] = -next_D / next_Delta_X

    start_index1 += nodes
# Matrix 3
start_index1 = 0
for idx, region in enumerate(regions_info):
    length, material, nodes = region
    properties = material_properties[material]

    Delta_X = length / nodes
    Removal = properties['Removal34']
    D = properties['D34']
    L1.append(np.sqrt(D / Removal))

    node_index_start = start_index1
    node_index_end = start_index1 + nodes - 1

    for i in range(nodes):
        node_index = start_index1 + i

        # Diagonal element for all nodes
        matrix_A3[node_index, node_index] = (Removal * Delta_X) + ((2 * D) / Delta_X)

        # Off-diagonal elements within the current region
        if i < nodes - 1:
            matrix_A3[node_index, node_index + 1] = -D / Delta_X
            matrix_A3[node_index + 1, node_index] = -D / Delta_X

            # Handling the boundary between two regions
    if idx < len(regions_info) - 1:
        next_region = regions_info[idx + 1]
        next_length, next_material, next_nodes = next_region
        next_properties = material_properties[next_material]
        next_D = next_properties['D34']
        next_Removal = next_properties['Removal34']
        next_Delta_X = next_length / next_nodes
        print(next_material)

        # Modify the last node of current region considering the interface
        matrix_A3[node_index_end, node_index_end] /= 2
        matrix_A3[node_index_end, node_index_end] += (next_D / next_Delta_X) + (next_Removal * next_Delta_X / 2)
        if node_index_end + 1 < N:
            matrix_A3[node_index_end, node_index_end + 1] = -next_D / next_Delta_X
            matrix_A3[node_index_end + 1, node_index_end] = -next_D / next_Delta_X

    start_index1 += nodes

# Matrix A4
start_index1 = 0
for idx, region in enumerate(regions_info):
    length, material, nodes = region
    properties = material_properties[material]
    Delta_X = length / nodes
    Removal = properties['Removal44']
    D = properties['D44']

    node_index_start = start_index1
    node_index_end = start_index1 + nodes - 1

    for i in range(nodes):
        node_index = start_index1 + i

        # Diagonal element for all nodes
        matrix_A4[node_index, node_index] = (Removal * Delta_X) + ((2 * D) / Delta_X)

        # Off-diagonal elements within the current region
        if i < nodes - 1:
            matrix_A4[node_index, node_index + 1] = -D / Delta_X
            matrix_A4[node_index + 1, node_index] = -D / Delta_X

            # Handling the boundary between two regions
    if idx < len(regions_info) - 1:
        next_region = regions_info[idx + 1]
        next_length, next_material, next_nodes = next_region
        next_properties = material_properties[next_material]
        next_D = next_properties['D44']
        next_Removal = next_properties['Removal44']
        next_Delta_X = next_length / next_nodes

        # Modify the last node of current region considering the interface
        matrix_A4[node_index_end, node_index_end] /= 2
        matrix_A4[node_index_end, node_index_end] += (next_D / next_Delta_X) + (next_Removal * next_Delta_X / 2)
        if node_index_end + 1 < N:
            matrix_A4[node_index_end, node_index_end + 1] = -next_D / next_Delta_X
            matrix_A4[node_index_end + 1, node_index_end] = -next_D / next_Delta_X

    start_index1 += nodes


for k in  ['1','2','3','4']:
    for l in ['1','2','3','4']:
        start_index3 = 0

        for idx, region in enumerate(regions_info):
            length, material, nodes = region
            properties = material_properties[material]
            Delta_X = length / nodes
            nu_sigma_f12 = properties.get('nusigmaF'+k+'4',
                                          0.0)  # Default to 0 if not defined

            for i in range(nodes):
                node_index = start_index3 + i

                # Diagonal element for all nodes
                if idx < len(regions_info) - 1 and i == nodes - 1:
                    # Last node of the current region and first node of the next region
                    next_region = regions_info[idx + 1]
                    next_length, next_material, next_nodes = next_region
                    next_properties = material_properties[next_material]
                    next_Delta_X = next_length / next_nodes
                    next_nu_sigma_f12 = next_properties.get('nusigmaF'+k+'4', 0.0)

                    # Average properties at the interface
                    matrix_B[k+l][node_index, node_index] = ((Delta_X * nu_sigma_f12 * properties.get('Chi'+l+'4',0.0) + properties.get('SigmaS'+k+l,0.0) * Delta_X)/2) + ((next_Delta_X * next_nu_sigma_f12 * next_properties.get('Chi' + l + '4',0.0) + next_properties.get('SigmaS' + k + l,0.0) * next_Delta_X)/2)
                else:
                    matrix_B[k+l][node_index, node_index] = Delta_X * nu_sigma_f12 * properties.get('Chi'+l+'4',0.0) + properties.get('SigmaS'+k+l,0.0) * Delta_X

            start_index3 += nodes

# Non leakage probability


np.savetxt("matrix_A1.csv", matrix_A1, delimiter=",", fmt="%s")
# Section for calculation of K-effictive adn curve
flux1 = np.ones(N)
flux2 = np.ones(N)
flux3 = np.ones(N)
flux4 = np.ones(N)
k = 1.00
Source1 = matrix_B['11'].dot(flux1) + matrix_B['21'].dot(flux2) + matrix_B['31'].dot(flux3) + matrix_B['41'].dot(flux4)
Source2 = matrix_B['12'].dot(flux1) + matrix_B['22'].dot(flux2) + matrix_B['32'].dot(flux3) + matrix_B['42'].dot(flux4)
Source3 = matrix_B['13'].dot(flux1) + matrix_B['23'].dot(flux2) + matrix_B['33'].dot(flux3) + matrix_B['43'].dot(flux4)
Source4 = matrix_B['14'].dot(flux1) + matrix_B['24'].dot(flux2) + matrix_B['34'].dot(flux3) + matrix_B['44'].dot(flux4)


flxdiff1 = 1.0
flxdiff2 = 1.0
flxdiff3 = 1.0
flxdiff4 = 1.0
g = 1

while flxdiff1 > 0.0001 or flxdiff2 > 0.0001 or flxdiff3 > 0.0001 or flxdiff4 > 0.0001:
    oldflux1 = flux1.copy()
    oldflux2 = flux2.copy()
    oldflux3 = flux3.copy()
    oldflux4 = flux4.copy()
    oldk = k
    oldsource1 = Source1.copy()
    oldsource2 = Source2.copy()
    oldsource3 = Source3.copy()
    oldsource4 = Source4.copy()

    flux1 = (np.linalg.inv(matrix_A1).dot(oldsource1)) / oldk
    Source2 = matrix_B['12'].dot(flux1) + matrix_B['22'].dot(oldflux2) + matrix_B['32'].dot(oldflux3) + matrix_B['42'].dot(oldflux4)
    flux2 = (np.linalg.inv(matrix_A2).dot(oldsource2)) / oldk
    Source3 = matrix_B['13'].dot(flux1) + matrix_B['23'].dot(flux2) + matrix_B['33'].dot(oldflux3) + matrix_B['43'].dot(oldflux4)
    flux3 = (np.linalg.inv(matrix_A3).dot(oldsource3)) / oldk
    Source4 = matrix_B['14'].dot(flux1) + matrix_B['24'].dot(flux2) + matrix_B['34'].dot(flux3) + matrix_B['44'].dot(oldflux4)
    flux4 = (np.linalg.inv(matrix_A4).dot(oldsource4)) / oldk

    flxdiff1 = np.sqrt(np.sum((flux1 - oldflux1) ** 2))
    flxdiff2 = np.sqrt(np.sum((flux2 - oldflux2) ** 2))
    flxdiff3 = np.sqrt(np.sum((flux3 - oldflux3) ** 2))
    flxdiff4 = np.sqrt(np.sum((flux4 - oldflux4) ** 2))
    Source1 = matrix_B['11'].dot(flux1) + matrix_B['21'].dot(flux2) + matrix_B['31'].dot(flux3) + matrix_B['41'].dot(flux4)
    k = oldk * np.sum(Source1) / np.sum(oldsource1)
    g += 1
print(Source1)
    # Normalize the flux
S = 0
Pth = 3400 * (10 ** 6) * 1.602e19  # Convert MWth to eV/s
kappa = 200 * (10 ** 6)  # MeV per fission
dV = sum(lengths)

Sigma_f1 = 0.008476
Sigma_f2 = 0.18514
S = np.sum(Sigma_f1 * flux1 + Sigma_f2 * flux2)
C = 1#Pth / (kappa * S * dV ** 3)

normflux1 = C * flux1
normflux1 = [0,*normflux1,0]
normflux2 = C * flux2
normflux2 = [0,*normflux2,0]
normflux3 = C * flux3
normflux3 = [0,*normflux3,0]
normflux4 = C * flux4
normflux4 = [0,*normflux4,0]
print(g)
print(k)
start = 0
x_nodes = np.array(
    [])  # This will store all nodes' positions

# Create nodes for each region
for length, material, nodes in regions_info:
    # Generate evenly spaced nodes within this region
    if nodes > 1:  # More than one node, create nodes within the range
        region_nodes = np.linspace(start, start + length, nodes + 2)
    else:  # Handle case with only one node (potentially rare)
        region_nodes = np.array([start + length / 2])

        # Append the nodes to the full list, avoiding duplication of the last/first node
    if len(x_nodes) > 0:
        region_nodes = region_nodes[
                       1:]  # Skip the first node if not the first region
    x_nodes = np.append(x_nodes, region_nodes)

    # Update the start point for the next region
    start += length

# Plotting
# Calculate x_nodes and plot
x_nodes = np.array(
    [])  # Start with an empty array for all nodes
current_start = 0  # Start position for the first region
region_boundaries = [
    0]  # Start with the initial boundary

for length, material, nodes in regions_info:
    # Generate nodes for the current region
    current_nodes = np.linspace(current_start, current_start + length, nodes + 1)
    if len(x_nodes) > 0:
        current_nodes = current_nodes[
                        1:]  # Skip the first node if it's not the first region
    x_nodes = np.append(x_nodes, current_nodes)
    current_start += length  # Update the start position for the next region
    region_boundaries.append(
        current_start)  # Store the boundary of this region

# Create the plot
plt.figure(figsize=(12, 6))
x_nodes=[*x_nodes,sum(lengths)]
# Plot the nodes as points
#plt.scatter(x_nodes, [0] * len(x_nodes), color='blue', label='Nodes')

# Plot flux distributions
# Define `x_nodes` based on the full reactor length


# Plot the flux distributions based on the correctly defined `x_nodes_full`
plt.plot(x_nodes, normflux1, 'b',label='group 1')  # adjust indices as needed
plt.plot(x_nodes, normflux2, 'r', label='group2')
plt.plot(x_nodes, normflux3, 'g',label='group3')  # adjust indices as needed
plt.plot(x_nodes, normflux4, 'magenta', label='group4')
# Region boundaries adjustments
for boundary in region_boundaries[1:-1]:
    plt.axvline(x=boundary, color='black', linestyle='--',label='Boundary')

plt.title('Node Distribution and Neutron Flux Across Regions')
plt.xlabel('Position (cm)')
plt.ylabel('Normalized Neutron Flux')
plt.legend()
plt.grid(True)
plt.show()
