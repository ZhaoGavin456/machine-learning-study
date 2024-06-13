# system_dynamics.py

import numpy as np

def whole_system_521(x, u, dt, cyc_lift_force, cyc_tilt_force):
    w_em, p_p, p_1l, p_2l, d_l, v_l, p_1t, p_2t, d_t, v_t, SOC = x
    u_1, u_2, u_3, u_4 = u

    J_em = 50e-3
    B_em = 0.1
    Dls = 71e-6 / (2 * np.pi) * 0.6
    vol_pump = (np.pi / 4) * 0.0254**2 * 20
    Eh = 11000e+5

    max_pressure = 400e5
    norm_pressure = p_p / max_pressure
    norm_torque = norm_pressure * u_2
    norm_torque = np.clip(norm_torque, 1e-6, 1)

    max_speed = 300
    norm_speed = w_em / max_speed
    norm_speed = np.clip(norm_speed, 0, 1)

    # Mechanical efficiency
    p_1_mech, p_2_mech, p_3_mech, p_4_mech, p_5_mech, p_6_mech = 0.8164, 0.4431, 0.0505, -0.3586, -0.0523, -0.0121
    eta_pm = p_1_mech + p_2_mech * norm_torque + p_3_mech * norm_speed + p_4_mech * norm_torque**2 + p_5_mech * norm_speed**2 + p_6_mech * norm_speed * norm_torque
    eta_pm = np.clip(eta_pm, 0.1, 1)

    # Volumetric efficiency
    p_1_vol, p_2_vol, p_3_vol, p_4_vol, p_5_vol, p_6_vol = 0.9744, -0.0031, -0.048, 0, -0.0479, -0.0283
    eta_pv = p_1_vol + p_2_vol * norm_pressure + p_3_vol * norm_torque + p_4_vol * norm_pressure**2 + p_5_vol * norm_torque**2 + p_6_vol * norm_pressure / norm_torque
    eta_pv = np.clip(eta_pv, 0.1, 1)

    n_lift, lift_m, lift_A, lift_B = 2, 30 * 2, np.pi * (101.6e-3)**2 / 4 * 2, (np.pi * (101.6e-3)**2 / 4 - np.pi * (57.2e-3)**2 / 4) * 2
    stroke_lift, lift_Bm = 783e-3, 28000
    vol_ch_lift = stroke_lift * lift_A

    n_tilt, tilt_m, tilt_A, tilt_B = 1, 30, np.pi * (114.3e-3)**2 / 4, np.pi * (114.3e-3)**2 / 4 - np.pi * (63.5e-3)**2 / 4
    storke_tilt, tilt_Bm = 529.3e-3, 28000
    vol_ch_tilt = storke_tilt * tilt_A

    # Update motor speed
    w_em += (1 / J_em * (u_1 - u_2 * Dls * p_p / eta_pm - B_em * w_em)) * dt
    w_em = np.clip(w_em, 64, 300)

    # Update pump pressure with checks for non-negative sqrt arguments
    delta_p_lift = np.maximum(p_p - p_1l, 0)  # Ensure non-negative
    delta_p_tilt = np.maximum(p_p - p_1t, 0)  # Ensure non-negative

    p_p += Eh / vol_pump * (
        u_2 * Dls * eta_pv * w_em -
        0.62 * np.pi * 5e-3 * 11e-3 * u_3 * np.sqrt(2 / 850 * delta_p_lift) -
        0.62 * np.pi * 8e-3 * 11e-3 * u_4 * np.sqrt(2 / 850 * delta_p_tilt)
    ) * dt

    p_p = np.clip(p_p, 51, 300e5)

    # Update lift dynamics
    p_1l += Eh / vol_ch_lift * (
        0.62 * np.pi * 5e-3 * 11e-3 * u_3 * np.sqrt(2 / 850 * delta_p_lift) -
        v_l * lift_A
    ) * dt

    delta_p_2_lift = np.maximum(p_2l, 0)  # Ensure non-negative
    p_2l += Eh / vol_ch_lift * (
        v_l * lift_B -
        0.62 * np.pi * 5e-3 * 11e-3 * u_3 * np.sqrt(2 / 850 * delta_p_2_lift)
    ) * dt

    d_l += v_l * dt
    v_l += (1 / lift_m * (p_1l * lift_A - p_2l * lift_B - cyc_lift_force - lift_Bm * v_l)) * dt

    # Update tilt dynamics
    delta_p_tilt_1 = np.maximum(p_p - p_1t, 0)  # Ensure non-negative
    p_1t += Eh / vol_ch_tilt * (
        0.62 * np.pi * 8e-3 * 11e-3 * u_4 * np.sqrt(2 / 850 * delta_p_tilt_1) -
        v_t * tilt_A
    ) * dt

    delta_p_2_tilt = np.maximum(p_2t, 0)  # Ensure non-negative
    p_2t += Eh / vol_ch_tilt * (
        v_t * tilt_B -
        0.62 * np.pi * 8e-3 * 11e-3 * u_4 * np.sqrt(2 / 850 * delta_p_2_tilt)
    ) * dt

    d_t += v_t * dt
    v_t += (1 / tilt_m * (p_1t * tilt_A - p_2t * tilt_B - cyc_tilt_force - tilt_Bm * v_t)) * dt

    # Ensure pressures are non-negative
    p_1l = max(p_1l, 0)
    p_2l = max(p_2l, 0)
    p_1t = max(p_1t, 0)
    p_2t = max(p_2t, 0)

    # Motor efficiency calculation
    motor_efficiency_coeff = [0.035473, 0.0010522, 5.66e-14, 0.016365]
    current_speed = w_em
    current_torque = u_1
    motor_loss = (motor_efficiency_coeff[0] * current_torque**2 +
                  motor_efficiency_coeff[1] * current_speed**2 +
                  motor_efficiency_coeff[2] +
                  motor_efficiency_coeff[3] * current_torque * current_speed)

    motor_efficiency = (current_torque * current_speed) / (current_torque * current_speed + motor_loss)
    motor_efficiency = np.clip(motor_efficiency, 1e-5, 1)

    # Electric motor power and battery current
    P_electric_motor = current_torque * current_speed / motor_efficiency
    V_nominal = 800
    R_ba = 1
    I_battery = (V_nominal - np.sqrt(V_nominal**2 - 4 * P_electric_motor * R_ba)) / (2 * R_ba)
    SOC += (-I_battery / (120 * 3600)) * dt

    return [w_em, p_p, p_1l, p_2l, d_l, v_l, p_1t, p_2t, d_t, v_t, SOC], [eta_pm, eta_pv, motor_efficiency], [0, 0, 0], [0, 0, 0], I_battery, P_electric_motor
