import gurobipy as gp
from gurobipy import GRB

import numpy as np
import scipy.sparse as sp

# Create a new model
m = gp.Model("MPC")

# horizon
N = 10

# Create variables

x = np.array([])
z = np.array([])
u = np.array([])

for k in range(N+1):
    x = np.append(x,[m.addMVar(12, vtype=GRB.CONTINUOUS)])

for k in range(N+1):
    z = np.append(z,[m.addMVar(12, vtype=GRB.CONTINUOUS)])

for k in range(N):
    u = np.append(u,[m.addMVar(4, vtype=GRB.CONTINUOUS)])

# Bounds
u0 = 10.5916
umin = np.array([9.6, 9.6, 9.6, 9.6]) - u0
umax = np.array([13., 13., 13., 13.]) - u0
xmin = np.array([-np.pi/6,-np.pi/6,-np.inf,-np.inf,-np.inf,-1.,
-np.inf,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf])
xmax = np.array([ np.pi/6, np.pi/6, np.inf, np.inf, np.inf, np.inf,
np.inf, np.inf, np.inf, np.inf, np.inf, np.inf])

# Constants
x[0] = np.zeros(12)
xr = np.array([0.,0.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.])


Q = sp.diags([0., 0., 10., 10., 10., 10., 0., 0., 0., 5., 5., 5.])

R = 0.1*sp.eye(4)


A = sp.csc_matrix([
[1., 0., 0., 0., 0., 0., 0.1, 0., 0., 0., 0., 0. ],
[0., 1., 0., 0., 0., 0., 0., 0.1, 0., 0., 0., 0. ],
[0., 0., 1., 0., 0., 0., 0., 0., 0.1, 0., 0., 0. ],
[0.0488, 0., 0., 1., 0., 0., 0.0016, 0., 0., 0.0992, 0., 0. ],
[0., -0.0488, 0., 0., 1., 0., 0., -0.0016, 0., 0., 0.0992, 0. ],
[0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0.0992],
[0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0. ],
[0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0. ],
[0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0. ],
[0.9734, 0., 0., 0., 0., 0., 0.0488, 0., 0., 0.9846, 0., 0. ],
[0., -0.9734, 0., 0., 0., 0., 0., -0.0488, 0., 0., 0.9846, 0. ],
[0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.9846]
])

B = sp.csc_matrix([
[0., -0.0726, 0., 0.0726],
[-0.0726, 0., 0.0726, 0. ],
[-0.0152, 0.0152, -0.0152, 0.0152],
[-0., -0.0006, -0., 0.0006],
[0.0006, 0., -0.0006, 0.0000],
[0.0106, 0.0106, 0.0106, 0.0106],
[0, -1.4512, 0., 1.4512],
[-1.4512, 0., 1.4512, 0. ],
[-0.3049, 0.3049, -0.3049, 0.3049],
[-0., -0.0236, 0., 0.0236],
[0.0236, 0., -0.0236, 0. ],
[0.2107, 0.2107, 0.2107, 0.2107]])

# MPC Formulation
obj = 0
for k in range(N):
    diff = x[k] - xr
    m.addConstr(z[k] == diff)
    obj = obj + (z[k] @ Q @ z[k]) + (u[k] @ R @ u[k])
    m.addConstr(x[k+1] == A @ x[k] + B @ u[k])
    m.addConstr(umin <= u[k])
    m.addConstr(u[k] <= umax)
    m.addConstr(xmin <= x[k+1])
    m.addConstr(x[k+1] <= xmax)

diff = x[N] - xr
m.addConstr(z[N] == diff)
obj = obj + z[N] @ Q @ z[N]


m.setObjective(obj, GRB.MINIMIZE)
m.optimize()