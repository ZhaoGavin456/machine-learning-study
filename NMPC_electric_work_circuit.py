# main_nmpc.py

import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import scipy.io
from tqdm import tqdm  # Import tqdm for progress bar
from system_dynamics import whole_system_521
from cost_function import lift_and_tilt_cost
from plot_results import plot_results  # Import the new plotting function

# Load the MATLAB file
mat_data = scipy.io.loadmat('extracted_data_with_positions_521_long_cycle_1.mat')

# Extract data
orig_time = mat_data['time_extracted'].flatten()
lift_position = mat_data['lift_position'].flatten()
lift_force_extracted = mat_data['lift_force_extracted'].flatten()
tilt_position = mat_data['tilt_position'].flatten()
tilt_force_extracted = mat_data['tilt_force_extracted'].flatten()

# Simulation parameters
dt = 0.001 / 1
time = np.arange(0, orig_time[-1], dt)

# Interpolated reference and cyclical forces
d_ref_lift = np.interp(time, orig_time, lift_position)
cyc_lift_force = np.interp(time, orig_time, lift_force_extracted)
cyc_tilt_force = np.interp(time, orig_time, tilt_force_extracted)
d_ref_tilt = np.interp(time, orig_time, tilt_position)

# Initial conditions
SOC_initial = 0.8
x = [100, 5e5, 0, 0, 0, 0, 0, 0, 0, 0, SOC_initial]

# NMPC parameters
prediction_horizon = 20
control_horizon = 15

# Preallocate arrays for plotting
d_l = np.zeros_like(time)
p_1l_values = np.zeros_like(time)
p_2l_values = np.zeros_like(time)
v_l_values = np.zeros_like(time)
w_em_values = np.zeros_like(time)
p_p_values = np.zeros_like(time)
d_t = np.zeros_like(time)
p_1t_values = np.zeros_like(time)
p_2t_values = np.zeros_like(time)
v_t_values = np.zeros_like(time)
SOC_values = np.zeros_like(time)
eff_pump_pm = np.zeros_like(time)
eff_pump_pv = np.zeros_like(time)
Flow_value = np.zeros_like(time)
Torque_value = np.zeros_like(time)
Flow_value_after = np.zeros_like(time)
Torque_value_after = np.zeros_like(time)
Battery_current_value = np.zeros_like(time)
Motor_power_value = np.zeros_like(time)
cost_details_accumulated = {'position_error_lift': [], 'position_error_tilt': [], 'soc_change': [], 'total': []}

# Control inputs
U0 = np.zeros((control_horizon, 4))
U0[:, 0] = 70
U0[:, 1] = 0.9
U0[:, 2] = 0.02
U0[:, 3] = 0.045

lb = np.tile([0, 0.1, 0, 0], (control_horizon, 1))
ub = np.tile([500, 1, 1, 1], (control_horizon, 1))

# NMPC optimization function with reduced iterations
def nmpc_control(initial_state, target_state, control_horizon, weights):
    def objective_function(U):
        U = U.reshape((control_horizon, 4))
        state = initial_state.copy()
        total_cost = 0
        for i in range(control_horizon):
            state, _, _, _, _, _ = whole_system_521(state, U[i], dt, cyc_lift_force[i], cyc_tilt_force[i])
            cost, _ = lift_and_tilt_cost(U, state, dt, cyc_lift_force, cyc_tilt_force, d_ref_lift, d_ref_tilt, control_horizon)
            total_cost += cost
        return total_cost

    initial_guess = np.zeros(control_horizon * 4)
    bounds = [(0, 500), (0.1, 1), (0, 1), (0, 1)] * control_horizon
    options = {'maxiter': 1000}  # Set maximum iterations to 1000
    result = minimize(objective_function, initial_guess, bounds=bounds, method='SLSQP', options=options)
    return result.x.reshape((control_horizon, 4))

# Simulation loop with progress bar
control_inputs = np.zeros((4, len(time)))
for k in tqdm(range(len(time) - 1), desc='Simulating NMPC', ncols=100):
    if k % 10 == 0:
        future_cyc_lift_force = cyc_lift_force[k:k + prediction_horizon]
        future_cyc_tilt_force = cyc_tilt_force[k:k + prediction_horizon]
        future_d_ref_lift = d_ref_lift[k:k + prediction_horizon]
        future_d_ref_tilt = d_ref_tilt[k:k + prediction_horizon]
        U_opt = nmpc_control(x, [1, 1, 10, 100], control_horizon, {'lift': 1, 'tilt': 1, 'soc': 1})

    u = U_opt[0]
    control_inputs[:, k] = u
    x, eff, Flow, Torque, Battery_current, Motor_power = whole_system_521(x, u, dt, cyc_lift_force[k], cyc_tilt_force[k])

    w_em_values[k] = x[0]
    p_p_values[k] = x[1]
    p_1l_values[k] = x[2]
    p_2l_values[k] = x[3]
    d_l[k] = x[4]
    v_l_values[k] = x[5]
    p_1t_values[k] = x[6]
    p_2t_values[k] = x[7]
    d_t[k] = x[8]
    v_t_values[k] = x[9]
    SOC_values[k] = x[10]
    eff_pump_pm[k] = eff[0]
    eff_pump_pv[k] = eff[1]
    Flow_value[k] = Flow[0]
    Flow_value_after[k] = Flow[1]
    Torque_value[k] = Torque[0]
    Torque_value_after[k] = Torque[1]
    Battery_current_value[k] = Battery_current
    Motor_power_value[k] = Motor_power

# Calculate energy consumption
energy_consumption_kJ = np.cumsum(Motor_power_value) / 1000 * np.mean(np.diff(time))

# Plotting results using the new function
plot_results(time, cyc_lift_force, cyc_tilt_force, d_l, d_ref_lift, d_t, d_ref_tilt,
             control_inputs, w_em_values, Torque_value_after, eff_pump_pm, Motor_power_value, energy_consumption_kJ)
