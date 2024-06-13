import numpy as np
import matplotlib.pyplot as plt
import casadi as ca

# Constants for the pendulum
g = 9.81  # acceleration due to gravity (m/s^2)
L = 1.0   # length of the pendulum (m)
m = 1.0   # mass of the pendulum (kg)

# Time and control parameters
dt = 0.1  # time step
N = 20    # number of control intervals
total_time = 20  # total time for simulation

# Define symbolic variables
theta = ca.MX.sym('theta')
omega = ca.MX.sym('omega')
torque = ca.MX.sym('torque')

# Dynamics model (nonlinear)
theta_dot = omega
omega_dot = -g/L * ca.sin(theta) + 1/(m*L**2) * torque

# Discrete dynamics using Euler integration
theta_next = theta + dt * theta_dot
omega_next = omega + dt * omega_dot

# Create a CasADi function for the dynamics
dynamics = ca.Function('dynamics', [theta, omega, torque], [theta_next, omega_next])

# Optimization variables
theta_opt = ca.MX.sym('theta_opt', N+1)
omega_opt = ca.MX.sym('omega_opt', N+1)
torque_opt = ca.MX.sym('torque_opt', N)

# Cost function and constraints
cost = 0
constraints = [theta_opt[0] - np.pi, omega_opt[0] - 0]  # Initial conditions

for i in range(N):
    # Cost function (minimize deviation from upright position and control effort)
    cost += (theta_opt[i] - 0)**2 + 0.01 * torque_opt[i]**2

    # System dynamics constraints
    if i < N:
        theta_next, omega_next = dynamics(theta_opt[i], omega_opt[i], torque_opt[i])
        constraints += [theta_opt[i+1] - theta_next, omega_opt[i+1] - omega_next]

# Problem definition
opt_problem = {'x': ca.vertcat(theta_opt, omega_opt, torque_opt), 'f': cost, 'g': ca.vertcat(*constraints)}

# Solver options
opts = {'ipopt': {'print_level': 0}, 'print_time': False}
solver = ca.nlpsol('solver', 'ipopt', opt_problem, opts)

# Initial guess and bounds for the solver
x0 = np.zeros(N+1 + N+1 + N)  # Properly sized initial guess
lbg = np.zeros(len(constraints))
ubg = np.zeros(len(constraints))

# Solve the problem
sol = solver(x0=x0, lbg=lbg, ubg=ubg)
x_opt = sol['x'].full().flatten()

# Extract the optimal values
theta_vals = x_opt[0:N+1]
omega_vals = x_opt[N+1:2*N+2]
torque_vals = x_opt[2*N+2:]

# Time vectors for plotting
t_theta_omega = np.linspace(0, total_time, N+1)  # Time vector for theta and omega
t_torque = np.linspace(0, total_time - dt, N)    # Time vector for torque

# Plot results
plt.figure(figsize=(12, 8))
plt.subplot(311)
plt.plot(t_theta_omega, theta_vals, label='Theta (angle)')
plt.title('Pendulum Angle (Theta)')
plt.legend()

plt.subplot(312)
plt.plot(t_theta_omega, omega_vals, label='Omega (angular velocity)')
plt.title('Angular Velocity')
plt.legend()

plt.subplot(313)
plt.step(t_torque, torque_vals, where='post', label='Torque (control input)')
plt.title('Control Input (Torque)')
plt.legend()

plt.xlabel('Time (s)')
plt.tight_layout()
plt.show()
