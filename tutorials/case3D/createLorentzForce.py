import numpy as np
import matplotlib.pyplot as plt
from scipy import special as scsp

def vector_field(x, y, x_center, y_center, k):
    """
    Compute the vector field at point (x, y) due to two vortices centered at (x_center, y_center).
    
    Parameters:
        x: X-coordinate of the point.
        y: Y-coordinate of the point.
        x_center: X-coordinate of the center of the vortices.
        y_center: Y-coordinate of the center of the vortices.
        k: Strength of the vortices.
    
    Returns:
        Tuple representing the vector field (u, v) at point (x, y).
    """
    # Calculate distance squared from the two vortex centers
    d1_squared = (x - x_center)**2 + (y - y_center/2)**2
    d2_squared = (x - x_center)**2 + (y - 3*y_center/2)**2
    
    # Compute the vector field components using the vortex flow formula
    u = np.exp(-2*(x)/delta)
    v = 0
    
    return u, v

# Generate a grid of points
x = np.linspace(0, 0.1, 20)
y = np.linspace(0, 0.1, 20)
X, Y = np.meshgrid(x, y)

# Parameters for the vortices
x_center = 0.1
y_center = 0.05
k = 0.01

mu0 = 1.25663706212e-6
f = 3e9
sigma = 1.4e6
alpha = 1

B0 = 1
delta = 1/np.sqrt(np.pi*mu0*sigma*f)*2000
print(delta)
F_x = np.exp(-2*x/delta)
# F_y = B0**2/2/mu0*alpha
F_Lz = sigma

F_Lr = -20*(np.exp(-2*(0.1-X)/delta))*np.exp(-(Y - 0.05)**2 / (2 * 0.02**2))
F_Lz = -10*(np.exp(-2*(0.1-X)/delta))*np.exp(-(Y - 0.05)**2 / (2 * 0.02**2))
Q_J = 200*(np.exp(-2*(0.1-X)/delta))*np.exp(-(Y - 0.05)**2 / (2 * 0.02**2))

# Flatten arrays
X_flat = X.flatten()
Y_flat = Y.flatten()

np.savetxt('lorentz-force.csv', np.column_stack((X_flat, Y_flat, Y_flat*0, F_Lr.flatten(), F_Lz.flatten(), Q_J.flatten())), delimiter=',')



# Plot the vector field
plt.figure(figsize=(8, 6))
plt.contourf(X, Y, Q_J)
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Vector Field with Two Vortices')
plt.grid(True)

plt.show()