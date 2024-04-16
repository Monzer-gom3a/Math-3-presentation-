import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

# Paramete rs
Lx = 1  # Length of the domain in x-direction
Ly = 1  # Length of the domain in y-direction
Nx = 50  # Number of grid points in x-direction
Ny = 50  # Number of grid points in y-direction
T = 2  # Total time
Nt = 200  # Number of time steps
alpha = 0.01  # Thermal diffusivity

# Grid spacin g
dx = Lx / Nx
dy = Ly / Ny
dt = T / Nt

# Create grid
x = np.linspace(0, Lx, Nx)
y = np.linspace(0, Ly, Ny)
X, Y = np.meshgrid(x, y)

# Initial condition
u0 = np.sin(np.pi * X) * np.sin(np.pi * Y)

# Boundary conditions
leftBC = np.zeros(Ny)
rightBC = np.zeros(Ny)
topBC = np.zeros(Nx)
bottomBC = np.zeros(Nx)

# Initialize solution matrix
u = np.zeros((Ny, Nx, Nt+1))
u[:, :, 0] = u0

# Finite difference scheme
for n in range(Nt):
    for i in range(1, Nx-1):
        for j in range(1, Ny-1):
            u[j, i, n+1] = u[j, i, n] + alpha * dt * (
                (u[j, i+1, n] - 2*u[j, i, n] + u[j, i-1, n]) / dx**2 +
                (u[j+1, i, n] - 2*u[j, i, n] + u[j-1, i, n]) / dy**2
            )

    # Apply boundary conditions
    u[:, 0, n+1] = leftBC
    u[:, -1, n+1] = rightBC
    u[0, :, n+1] = bottomBC
    u[-1, :, n+1] = topBC

# Plotting
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

def update_plot(frame, u, X, Y, dt):
    ax.clear()
    surf = ax.plot_surface(X, Y, u[:, :, frame], cmap='viridis')
    ax.set_title(f'Temperature distribution at time t = {frame * dt:.2f}')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('Temperature')
    ax.set_xlim(0, Lx)
    ax.set_ylim(0, Ly)
    ax.set_zlim(-1, 1)

# Create an animation
ani = FuncAnimation(fig, update_plot, frames=Nt, fargs=(u, X, Y, dt), interval=100)
plt.show()
