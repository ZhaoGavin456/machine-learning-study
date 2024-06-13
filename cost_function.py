# cost_function.py

import numpy as np
from system_dynamics import whole_system_521

def lift_and_tilt_cost(U, x, dt, cyc_lift_forces, cyc_tilt_forces, d_refs_lift, d_refs_tilt, prediction_horizon):
    total_cost = 0
    cost_details = {'position_error_lift': [], 'position_error_tilt': [], 'soc_change': [], 'total': []}

    weight_position_error_lift = 60000
    weight_position_error_tilt = 20000
    weight_soc_change = 1

    x_pred = x
    for i in range(min(len(U), prediction_horizon)):
        idx_lift = min(i, len(cyc_lift_forces) - 1)
        idx_tilt = min(i, len(cyc_tilt_forces) - 1)
        x_pred, _, _, _, _, _ = whole_system_521(x_pred, U[i], dt, cyc_lift_forces[idx_lift], cyc_tilt_forces[idx_tilt])

        position_error_lift = ((x_pred[4] - d_refs_lift[idx_lift]) / 1.0) ** 2
        position_error_tilt = ((x_pred[8] - d_refs_tilt[idx_tilt]) / 1.0) ** 2
        soc_change = x_pred[-1] - x[10]

        weighted_position_error_lift = weight_position_error_lift * position_error_lift
        weighted_position_error_tilt = weight_position_error_tilt * position_error_tilt
        weighted_soc_change = weight_soc_change * soc_change

        current_total_cost = weighted_position_error_lift + weighted_position_error_tilt + weighted_soc_change
        total_cost += current_total_cost

        cost_details['position_error_lift'].append(weighted_position_error_lift)
        cost_details['position_error_tilt'].append(weighted_position_error_tilt)
        cost_details['soc_change'].append(weighted_soc_change)
        cost_details['total'].append(current_total_cost)

    return total_cost, cost_details
