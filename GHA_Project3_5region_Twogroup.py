import numpy as np
import matplotlib.pyplot as plt
import math

### Values###

# Define material properties
material_properties = {
    "PWR": {
        "nusigmaF12": 0.008476,
        "nusigmaF22": 0.18514,
        "sigmaF12": 0.003320,
        "sigmaF22": 0.07537,
        "SigmaA12": 0.01207,
        "SigmaA22": 0.1210,
        "D12": 1.2627,
        "D22": 0.3543,
        "Removal12": 0.02619,
        "Removal22": 0.1210,
        "SigmaS12": 0.01412
    },
    "MOX": {
        "nusigmaF12": 0.00930,
        "nusigmaF22": 0.261,
        "sigmaF12": 0.0,
        "sigmaF22": 0.0,
        "SigmaA12": 0.0,
        "SigmaA22": 0.1637,
        "D12": 1.468,
        "D22": 0.296,
        "Removal12": 0.0285,
        "Removal22": 0.1637,
        "SigmaS12": 0.0178
    },
    "Reflector": {
        "nusigmaF12": 0.0,
        "nusigmaF22": 0.0,
        "sigmaF12": 0.0,
        "sigmaF22": 0.0,
        "SigmaA12": 0.0004,
        "SigmaA22": 0.0197,
        "D12": 1.13,
        "D22": 0.16,
        "Removal12": 0.0494,
        "Removal22": 0.0197,
        "SigmaS12": 0.0494
    },
    "FeulU": {
        "nusigmaF12": 0.0127,
        "nusigmaF22": 0.220,
        "sigmaF12": 0.0,
        "sigmaF22": 0.0,
        "SigmaA12": 0.0,
        "SigmaA22": 0.129,
        "D12": 1.4899,
        "D22": 0.323,
        "Removal12": 0.029,
        "Removal22": 0.129,
        "SigmaS12": 0.0178
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
Z = 20
Y = 6.614
regions_info = [(Y, 'Reflector', Z),(Y, 'FeulU', Z) , (Y, 'MOX', Z),(Y, 'FeulU', Z),(Y, 'Reflector', Z)]
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
matrix_B1 = np.zeros((N, N))
matrix_B2 = np.zeros((N, N))
matrix_B3 = np.zeros((N, N))

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
    Removal = properties['Removal12']
    D = properties['D12']
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
        next_D = next_properties['D12']
        next_Removal = next_properties['Removal12']
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
    Removal = properties['SigmaA22']
    D = properties['D22']

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
        next_D = next_properties['D22']
        next_Removal = next_properties['SigmaA22']
        next_Delta_X = next_length / next_nodes

        # Modify the last node of current region considering the interface
        matrix_A2[node_index_end, node_index_end] /= 2
        matrix_A2[node_index_end, node_index_end] += (next_D / next_Delta_X) + (next_Removal * next_Delta_X / 2)
        if node_index_end + 1 < N:
            matrix_A2[node_index_end, node_index_end + 1] = -next_D / next_Delta_X
            matrix_A2[node_index_end + 1, node_index_end] = -next_D / next_Delta_X

    start_index1 += nodes

# Matrix B1
start_index3 = 0
for idx, region in enumerate(regions_info):
    length, material, nodes = region
    properties = material_properties[material]
    Delta_X = length / nodes
    nu_sigma_f12 = properties.get('nusigmaF12',
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
            next_nu_sigma_f12 = next_properties.get('nusigmaF12', 0.0)

            # Average properties at the interface
            matrix_B1[node_index, node_index] = ((Delta_X * nu_sigma_f12) / 2) + (
                        (next_Delta_X * next_nu_sigma_f12) / 2)
        else:
            matrix_B1[node_index, node_index] = Delta_X * nu_sigma_f12

    start_index3 += nodes
# Matrix B2
start_index4 = 0
for idx, region in enumerate(regions_info):
    length, material, nodes = region
    properties = material_properties[material]
    Delta_X = length / nodes
    nu_sigma_f22 = properties.get('nusigmaF22',
                                  0.0)  # Default to 0 if not defined

    for i in range(nodes):
        node_index = start_index4 + i

        # Diagonal element for all nodes
        if idx < len(regions_info) - 1 and i == nodes - 1:
            # Last node of the current region and first node of the next region
            next_region = regions_info[idx + 1]
            next_length, next_material, next_nodes = next_region
            next_properties = material_properties[next_material]
            next_Delta_X = next_length / next_nodes
            next_nu_sigma_f22 = next_properties.get('nusigmaF22', 0.0)

            # Average properties at the interface
            matrix_B2[node_index, node_index] = ((Delta_X * nu_sigma_f22) / 2) + (
                        (next_Delta_X * next_nu_sigma_f22) / 2)
        else:
            matrix_B2[node_index, node_index] = Delta_X * nu_sigma_f22

    start_index4 += nodes
# matrix B3
start_index5 = 0
for idx, region in enumerate(regions_info):
    length, material, nodes = region
    properties = material_properties[material]
    Delta_X = length / nodes
    SigmaS12 = properties.get('SigmaS12',
                              0.0)  # Default to 0 if not defined

    for i in range(nodes):
        node_index = start_index5 + i

        # Diagonal element for all nodes
        if idx < len(regions_info) - 1 and i == nodes - 1:
            # Last node of the current region and first node of the next region
            next_region = regions_info[idx + 1]
            next_length, next_material, next_nodes = next_region
            next_properties = material_properties[next_material]
            next_Delta_X = next_length / next_nodes
            next_SigmaS12 = next_properties.get('SigmaS12', 0.0)

            # Average properties at the interface
            matrix_B3[node_index, node_index] = ((Delta_X * SigmaS12) / 2) + ((next_Delta_X * next_SigmaS12) / 2)
        else:
            matrix_B3[node_index, node_index] = Delta_X * SigmaS12

    start_index5 += nodes
# Non leakage probability
print(L1)
print(L2)

np.savetxt("matrix_A1.csv", matrix_A1, delimiter=",", fmt="%s")
# Section for calculation of K-effictive adn curve
flux1 = np.ones(N)
flux2 = np.ones(N)
k = 1.00
Source1 = matrix_B1.dot(flux1) + matrix_B2.dot(
    flux2)  # fast
Source2 = matrix_B3.dot(
    flux1)  # Thermal

flxdiff1 = 1.0
flxdiff2 = 1.0
g = 1

while flxdiff1 > 0.0001 or flxdiff2 > 0.0001:
    oldflux1 = flux1.copy()
    oldflux2 = flux2.copy()
    oldk = k
    oldsource1 = Source1.copy()
    oldsource2 = Source2.copy()

    flux1 = (np.linalg.inv(matrix_A1).dot(oldsource1)) / oldk
    flxdiff1 = np.sqrt(np.sum((flux1 - oldflux1) ** 2))

    Source2 = matrix_B3.dot(flux1)
    flux2 = np.linalg.inv(matrix_A2).dot(Source2)

    flxdiff2 = np.sqrt(np.sum((flux2 - oldflux2) ** 2))

    Source1 = (matrix_B1.dot(flux1) + matrix_B2.dot(flux2))

    k = oldk * np.sum(Source1) / np.sum(oldsource1)
    g += 1

    # Normalize the flux
S = 0
Pth = 3400 * (10 ** 6) * 1.602e19  # Convert MWth to eV/s
kappa = 200 * (10 ** 6)  # MeV per fission
dV = sum(lengths)
print(dV)
Sigma_f1 = 0.008476
Sigma_f2 = 0.18514
S = np.sum(Sigma_f1 * flux1 + Sigma_f2 * flux2)
C = Pth / (kappa * S * dV ** 3)

normflux1 = C * flux1
normflux1 = [0,*normflux1,0]
normflux2 = C * flux2
normflux2 = [0,*normflux2,0]

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
plt.plot(x_nodes, normflux1, 'b-+',
         label='Fast Flux')  # adjust indices as needed
plt.plot(x_nodes, normflux2, 'r-+', label='Thermal Flux')

# Region boundaries adjustments
for boundary in region_boundaries[1:-1]:
    plt.axvline(x=boundary, color='black', linestyle='--',label='Boundary')

plt.title('Node Distribution and Neutron Flux Across Regions')
plt.xlabel('Position (cm)')
plt.ylabel('Normalized Neutron Flux')
plt.legend()
plt.grid(True)
plt.show()
