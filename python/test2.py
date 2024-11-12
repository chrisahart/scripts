import numpy as np

import matplotlib.pyplot as plt

# Define the coefficients for each molecular orbital

coefficients = [
    [1.003, 0.132, 0.003, 0.000, 0.000, -0.000, 0.000, 0.000, 0.000, 0.000, -0.002, -0.008, 0.001, 0.000, 0.000, 0.001, 0.003, -0.001, 0.000, 0.000, 0.000, 0.000, -0.000, -0.000, -0.000, -0.000],
    [0.000, 0.000, 0.000, 0.163, 0.975, -0.091, -0.002, -0.013, 0.001, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, -0.003, -0.006, 0.001, 0.000, -0.000],
    [0.000, 0.000, 0.000, 0.786, -0.076, 0.600, -0.012, 0.001, -0.009, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.001, 0.000, -0.005, 0.001, -0.004, 0.002, -0.001],
    [0.000, 0.000, 0.000, -0.583, 0.171, 0.784, 0.009, -0.003, -0.012, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, -0.001, -0.002, 0.003, -0.003, -0.004, 0.000, -0.000],
    [0.030, -0.990, -0.002, 0.000, 0.000, -0.000, 0.000, 0.000, 0.000, 0.003, -0.032, -0.108, 0.018, 0.002, 0.000, -0.002, -0.006, 0.001, 0.000, 0.000, 0.001, 0.000, -0.000, -0.000, -0.000, -0.000],
    [0.000, 0.000, 0.000, 0.022, 0.133, -0.012, 0.042, 0.248, -0.023, 0.958, 0.061, -0.023, -0.174, -0.159, 0.029, 0.002, -0.001, -0.005, -0.005, 0.000, 0.000, 0.000, 0.001, 0.000, 0.000, 0.000]
]

# Define the grid for the 3D plot

x = np.linspace(-20, 20, 200)
y = np.linspace(-20, 20, 200)
z = np.linspace(-20, 20, 200)
x, y, z = np.meshgrid(x, y, z)

# Create a function to generate the wavefunction for a given MO

def generate_wavefunction(x, y, z, coefficients):
    wavefunction = np.zeros_like(x)
    for i, coeff in enumerate(coefficients):
        # For simplicity, we assume a Gaussian-like function for each AO

        wavefunction += coeff * np.exp(-0.5 * (x**2 + y**2 + z**2)) * np.cos(i * np.arctan2(y, x))
    return wavefunction

# Generate the wavefunctions for all molecular orbitals and sum them

total_wavefunction = np.zeros_like(x)
for i, coeffs in enumerate(coefficients):
    total_wavefunction += generate_wavefunction(x, y, z, coeffs)

# Average the total wavefunction along the x and y coordinates for each z coordinate

averaged_wavefunction = np.mean(total_wavefunction, axis=(1, 2))

# Plot the averaged wavefunction along the z-axis with a logarithmic y-axis

plt.figure(figsize=(10, 6))
plt.plot(z[0, 0, :], np.abs(averaged_wavefunction), label='Sum of All MOs')
plt.yscale('log')
plt.title('Sum of All Molecular Orbitals Averaged Along Y and Z')
plt.xlabel('z')
plt.ylabel('Averaged Wavefunction (log scale)')
plt.legend()
plt.show()