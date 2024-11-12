import numpy as np

import matplotlib.pyplot as plt

# O
exponents = [
    10.389228018317, 3.849621072005, 1.388401188741, 0.496955043655, 0.162491615040
]
coefficients = [
    [0.126240722900, 0.069215797900, -0.061302037200, -0.026862701100, 0.029845227500],
    [0.139933704300, 0.115634538900, -0.190087511700, -0.006283021000, 0.060939733900],
    [-0.434348231700, -0.322839719400, -0.377726982800, -0.224839187800, 0.732321580100],
    [-0.852791790900, -0.095944016600, -0.454266086000, 0.380324658600, 0.893564918400],
    [-0.242351537800, 1.102830348700, -0.257388983000, 1.054102919900, 0.152954188700]
]


# Hf
exponents = [
    2.671196497548, 2.597200237263, 1.091215248046, 0.482780902120, 0.184203719686, 0.060257634881
]
coefficients = [
    [0.788939427855, 0.166407997203, -0.808044747643, -0.301974408389, 0.200618001886, -0.047796823905, -0.576358120326, 0.136965103461],
    [-0.135944993309, 0.121000538909, 0.502710504170, 0.293385766135, -0.131814094303, 0.040716807355, 0.121047786291, -0.101127257056],
    [-0.911502384031, -0.450097538431, -0.571536407023, 0.044171452367, -0.119671510085, 0.281226003531, -0.260083157654, -0.239844496668],
    [-0.403765732131, -0.527116158052, 1.567869205746, 0.050912843752, -0.209054741982, 0.563159769901, -0.061407406614, -0.448266650992],
    [-0.116982787286, 0.367518491637, -3.196670125867, 0.006790589259, -0.050569947887, 0.537885417612, -0.642967474054, -0.402343230211],
    [-0.118019349325, 0.917540530320, 1.895484387467, -0.000081217891, -0.947169815010, 0.138650274135, 0.976981432228, 1.144416860469]
]



# Generate a range of r values
r = np.linspace(0, 20, 400)

# Define the Gaussian basis function

def gaussian_basis(r, alpha, N):
    return N * np.exp(-alpha * r**2)

# Compute the contracted basis function

phi = 0

for i, alpha in enumerate(exponents):
    N = (2 * alpha / np.pi) ** (3/4)
    for j, c in enumerate(coefficients[i]):
        phi += c * gaussian_basis(r, alpha, N)

print(phi)
print(np.sum((phi)**2))
phi = phi * 11/ np.sum((phi)**2)
print(np.sum((phi)**2))

# Plot the contracted basis function with a logarithmic y-axis

plt.figure(figsize=(10, 6))
plt.plot(r, phi, label='Hf DZVP-MOLOPT-SR-GTH Basis Function')
# plt.plot(r, np.abs(phi), label='Hf DZVP-MOLOPT-SR-GTH Basis Function')
plt.hlines(1e-5, 0, 100, 'r', alpha=0.5)
plt.hlines(1e-10, 0, 100, 'r', alpha=0.5)
plt.title('Hf DZVP-MOLOPT-SR-GTH Basis Function (Log Scale)')
plt.xlabel('r (distance from nucleus)')
plt.ylabel('φ(r)')
plt.yscale('log')  # Set the y-axis to a logarithmic scale

plt.legend()
plt.grid(True)
plt.show()