import gurobipy as gp
from gurobipy import GRB
import numpy as np
import matplotlib.pyplot as plt

# System parameters
R = 0.1  # Ohms (electric motor resistance)
L = 0.01  # Henry (electric motor inductance)
K = 0.01  # Nm/A (electric motor constant)
Km = 0.1  # Nm/V (electric motor torque constant)
J = 0.01  # kg m^2 (rotor inertia)
B = 0.1  # Nm s/rad (damping coefficient)
Tl = 0.05  # Nm (load torque)
Kp = 1.0  # Proportional gain

# MPC parameters
horizon = 200  # MPC horizon
dt = 0.1       # time step

# Sine wave parameters
amplitude = 0.5  # amplitude of sine wave
frequency = 0.1  # frequency of sine wave (rad/s)
phase_shift = np.pi / 2  # phase shift of sine wave (radians)

# Create sine wave reference function
def reference_position(t):
    return amplitude * np.sin(2 * np.pi * frequency * t + phase_shift)

# Create model
model = gp.Model("Electric Hydraulic Servo System")

# Decision variables
u = {}
for t in range(horizon):
    u[t] = model.addVar(lb=-10, ub=10, name=f"u[{t}]")  # voltage

theta = {}
omega = {}
for t in range(horizon + 1):
    theta[t] = model.addVar(lb=-gp.GRB.INFINITY, ub=gp.GRB.INFINITY, name=f"theta[{t}]")  # position
    omega[t] = model.addVar(lb=-gp.GRB.INFINITY, ub=gp.GRB.INFINITY, name=f"omega[{t}]")  # velocity

# Dynamics constraints
for t in range(horizon):
    model.addConstr(omega[t + 1] == omega[t] + (K * (u[t] - Km * omega[t]) - B * omega[t] - Tl) / J * dt)
    model.addConstr(theta[t + 1] == theta[t] + omega[t] * dt)

# Reference position constraint
for t in range(horizon + 1):
    model.addConstr(theta[t] == reference_position(t * dt))

# Objective (minimize position error)
error = gp.quicksum((reference_position(t * dt) - theta[t]) ** 2 for t in range(horizon + 1))
model.setObjective(error, GRB.MINIMIZE)

# Optimize
model.optimize()

# Extract results
positions = [theta[t].X for t in range(horizon + 1)]
voltages = [u[t].X for t in range(horizon)]

# Plot results
time = [t * dt for t in range(horizon + 1)]
reference_positions = [reference_position(t * dt) for t in range(horizon + 1)]

plt.figure(figsize=(10, 6))
plt.subplot(2, 1, 1)
plt.plot(time, positions, label='Actual Position')
plt.plot(time, reference_positions, 'r--', label='Reference Position')
plt.xlabel('Time (s)')
plt.ylabel('Position (m)')
plt.title('Hydraulic Cylinder Position')
plt.legend()

plt.subplot(2, 1, 2)
plt.step(time[:-1], voltages, where='post')
plt.xlabel('Time (s)')
plt.ylabel('Voltage (V)')
plt.title('Applied Voltage')

plt.tight_layout()
plt.show()
