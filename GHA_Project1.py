import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
import math

# For Deliverable A
S = 10 ** 12
W = 10
Delta_X = W / 4
N1 = 5
x1 = np.linspace(0, W, 500)

Sigma_a2_tr = 3.62e-2
Sigma_a2_abs = 0.1532
Dm = 1 / (3 * (Sigma_a2_abs + Sigma_a2_tr))
L = math.sqrt(Dm / Sigma_a2_abs)


# exponenent_1 = np.exp((2*W)/L)
# C1 = (-constant/(exponenent_1+1))
# C2= (constant-(constant/(exponenent_1+1)))
# flux = (C1*np.exp(-(x/L)))+(C2*np.exp(x/L))
def analytical_solution(x, L, W, S, Dm):
    constant = ((S * L) / (2 * Dm))
    num = np.exp(-(x / L)) - np.exp(((-((2 * W) - x)) / L))
    dom = 1 + np.exp(-(2 * W) / L)
    flux = constant * (num / dom)
    return flux


# Numerical analysis
# set up matrix
matrix_A = np.zeros((4, 4))
# populating the matrix manually
matrix_A[0, 0] = ((Dm / Delta_X ** 2) + (Sigma_a2_abs / 2))
matrix_A[0, 1] = -Dm / Delta_X ** 2
matrix_A[1, 0] = -Dm / Delta_X ** 2
matrix_A[1, 1] = Sigma_a2_abs + (2 * Dm / Delta_X ** 2)
matrix_A[1, 2] = -Dm / Delta_X ** 2
matrix_A[2, 1] = -Dm / Delta_X ** 2
matrix_A[2, 2] = Sigma_a2_abs + (2 * Dm / Delta_X ** 2)
matrix_A[2, 3] = -Dm / Delta_X ** 2
matrix_A[3, 2] = -Dm / Delta_X ** 2
matrix_A[3, 3] = Sigma_a2_abs + (2 * Dm / Delta_X ** 2)

print(matrix_A)
# Phi vector
Phi_vector = np.array([S / (2 * Delta_X), 0, 0, 0])
flux_num = np.linalg.inv(matrix_A).dot(Phi_vector)
flux_num = np.append(flux_num, 0)  # This line is to act as adding in the most right case which is 0

x2 = np.linspace(0, W, N1)
analytical_at_numerical_x = analytical_solution(x2, L, W, S, Dm)

# Calculate the difference
difference = analytical_at_numerical_x[:-1] - flux_num[:-1]  # Exclude the artificially added 0

# Plotting
plt.figure(figsize=(10, 7))
plt.title('Flux vs. Position for W= 10cm')
plt.xlabel('Position(cm)')
plt.ylabel('Flux')
# Original plot
plt.plot(x1, analytical_solution(x1, L, W, S, Dm), linestyle='-', color='blue', label="Analytical")
plt.plot(x2, flux_num, marker='o', linestyle='-', color='red', label="Numerical")

# Inset for difference
axins = inset_axes(plt.gca(), width="30%", height="30%", loc='upper right')
axins.plot(x2[:-1], difference, marker='o', linestyle='-', color='green', label="Difference")
mark_inset(plt.gca(), axins, loc1=2, loc2=3, fc="none", ec="0.5")
axins.set_title('Difference Plot')
axins.grid(True)

plt.xlabel('Position (cm)')
plt.ylabel('Flux')
plt.grid(True)
plt.legend()
plt.savefig("Deliverable_A_with_difference.svg", format='svg')
plt.show()

# Deliverable B
# Defining variables
N_values = [5, 10, 20, 30, 40, 50]
W_values = [5, 10, 50]


# Using the general equations from part A we can make an automated numerical analysis part
def solve_numerical(N, W, S, Sigma_a2_abs, Dm):
    # Calculate Delta_X based on the width W and the number of nodes N
    Delta_X = W / (N - 1)

    # Initialize the A matrix and the b vector (source term vector)
    matrix_A = np.zeros((N - 1, N - 1))  # Adjusted for the boundary conditions

    # Fill the A matrix
    for i in range(N - 1):
        if i == 0:
            # First row (left boundary condition)
            matrix_A[i, i] = ((Dm / Delta_X ** 2) + (Sigma_a2_abs / 2))
            matrix_A[i, i + 1] = -Dm / Delta_X ** 2
        elif i == N - 2:
            # Last row (right boundary condition before appending)
            matrix_A[i, i - 1] = -Dm / Delta_X ** 2
            matrix_A[i, i] = Sigma_a2_abs + (2 * Dm / Delta_X ** 2)
        else:
            # Internal nodes
            matrix_A[i, i - 1] = -Dm / Delta_X ** 2
            matrix_A[i, i] = Sigma_a2_abs + (2 * Dm / Delta_X ** 2)
            matrix_A[i, i + 1] = -Dm / Delta_X ** 2

    # Solve the system
    # Source term vector
    b_vector = np.zeros(N - 1)
    b_vector[0] = S / (2 * Delta_X)  # Adjusted source term for the first node, considering symmetry
    flux_num = np.linalg.inv(matrix_A).dot(b_vector)
    flux_num = np.append(flux_num, 0)  # This line is to act as adding in the most right case which is 0

    # Generate x values for plotting
    x_values = np.linspace(0, W, N)

    return x_values, flux_num


for W in W_values:
    fig, ax = plt.figure(figsize=(10, 6)), plt.gca()
    x_analytical = np.linspace(0, W, 500)
    analytical_flux = analytical_solution(x_analytical, math.sqrt(Dm / Sigma_a2_abs), W, S, Dm)
    ax.plot(x_analytical, analytical_flux, linestyle='-', color='black', label="Analytical")

    for N in N_values:
        x_values, flux_num = solve_numerical(N, W, S, Sigma_a2_abs, Dm)
        ax.plot(x_values, flux_num, marker='', linestyle='--', label=f"Numerical N={N}")

    # Calculate the middle range for the inset
    x_middle = W / 2.5
    zoom_width = W / 10  # Define zoom width as a fraction of the total width
    x1, x2 = x_middle - zoom_width / 2, x_middle + zoom_width / 2
    # Estimate y-range around the middle x-range for a dynamic y-axis zoom
    y_mid_values = analytical_flux[(x_analytical >= x1) & (x_analytical <= x2)]
    y1, y2 = min(y_mid_values), max(y_mid_values)
    y_margin = (y2 - y1) * 0.3  # Adding a margin to the y-range
    y1, y2 = y1 - y_margin, y2 + y_margin

    # Inset axes
    ax_inset = inset_axes(ax, width="30%", height="30%", loc=1)
    ax_inset.plot(x_analytical, analytical_flux, linestyle='-', color='black')
    for N in N_values:
        x_values, flux_num = solve_numerical(N, W, S, Sigma_a2_abs, Dm)
        ax_inset.plot(x_values, flux_num, marker='', linestyle='--')
    # Set dynamic inset zoom range
    ax_inset.set_xlim(x1, x2)
    ax_inset.set_ylim(y1, y2)
    mark_inset(ax, ax_inset, loc1=2, loc2=4, fc="none", ec="0.5")
    ax.set_title(f'Flux vs. Position for W= {W}cm')
    ax.set_xlabel('Position (cm)')
    ax.set_ylabel('Flux')
    ax.legend(loc='lower right', bbox_to_anchor=(1, 0))
    ax.grid(True)
    filename = f"Flux_vs_Position_W{W}_Nodes_with_Dynamic_Inset.svg"
    plt.savefig(filename, format='svg')
    plt.show()

# Deliverable C
f_values = {0.5, 1, 1.5}
W4 = 10


def analytical_solution_C(x, L, W, S, f):
    flux2 = ((-S * f * np.exp(W / L) + S * np.exp(2 * W / L)) / (np.exp(2 * W / L) - 1) * np.exp(-x / L)) + (
                (S * f * np.exp(W / L) - S) / (np.exp(2 * W / L) - 1) * np.exp(x / L))
    return flux2


def solve_numerical_C(N, W, S, Sigma_a2_abs, Dm, f):
    Delta_X = W / (N - 1)
    matrix_A = np.zeros((N, N))
    b_vector = np.zeros(N)

    for i in range(N):
        if i == 0:
            matrix_A[i, i] = 1
        elif i == N - 1:
            matrix_A[i, i] = 1
        else:
            matrix_A[i, i - 1] = -Dm / Delta_X ** 2
            matrix_A[i, i] = Sigma_a2_abs + (2 * Dm / Delta_X ** 2)
            matrix_A[i, i + 1] = -Dm / Delta_X ** 2

    b_vector[0] = S
    b_vector[-1] = f * S

    flux_num_C = np.linalg.solve(matrix_A, b_vector)

    x_values_C = np.linspace(0, W, N)
    return x_values_C, flux_num_C


# Loop through each f value and generate plots
for f in f_values:
    fig, ax = plt.subplots(figsize=(10, 6))

    # Plot the analytical solution
    x_analytical = np.linspace(0, W4, 500)
    analytical_flux_C = analytical_solution_C(x_analytical, L, W4, S, f)
    ax.plot(x_analytical, analytical_flux_C, linestyle='-', color='black', label=f"Analytical, f={f}")

    # Plot numerical solutions for different N values
    for N in N_values:
        x_values_C, flux_num_C = solve_numerical_C(N, W4, S, Sigma_a2_abs, Dm, f)
        ax.plot(x_values_C, flux_num_C, marker='o', linestyle='--', label=f"Numerical N={N}, f={f}")

    # Setup for dynamic zoom for the inset
    # Define middle range for the inset based on a target area or feature of interest
    target_x_range = (W4 * 0.4, W4 * 0.6)  # Adjust as necessary to focus on a specific part of the plot
    target_flux_C = analytical_flux_C[(x_analytical >= target_x_range[0]) & (x_analytical <= target_x_range[1])]
    # Calculate dynamic y-limits based on target range
    y3, y4 = target_flux_C.min(), target_flux_C.max()
    y_margin = (y4 - y3) * 0.3  # Add some margin to the y-limits
    y3 -= y_margin
    y4 += y_margin

    # Create inset axes and plot the same data
    ax_inset = inset_axes(ax, width="30%", height="30%", loc=1)
    ax_inset.plot(x_analytical, analytical_flux_C, linestyle='-', color='black')
    for N in N_values:
        x_values_C, flux_num_C = solve_numerical_C(N, W4, S, Sigma_a2_abs, Dm, f)
        ax_inset.plot(x_values_C, flux_num_C, marker='o', linestyle='--')

    # Apply dynamic zoom range to inset
    ax_inset.set_xlim(target_x_range)
    ax_inset.set_ylim(y3, y4)
    mark_inset(ax, ax_inset, loc1=2, loc2=4, fc="none", ec="0.5")

    ax.set_title(f'Flux vs. Position for W={W4}cm and f={f}')
    ax.set_xlabel('Position (cm)')
    ax.set_ylabel('Flux')
    ax.legend(loc='lower right', bbox_to_anchor=(1, 0))
    ax.grid(True)

    filename = f"Flux_vs_Position_W{W4}_f={f}_with_Dynamic_Inset.svg"
    plt.savefig(filename, format='svg')
    plt.show()
