# plot_results.py

import matplotlib.pyplot as plt
import numpy as np

def plot_results(time, cyc_lift_force, cyc_tilt_force, d_l, d_ref_lift, d_t, d_ref_tilt,
                 control_inputs, w_em_values, Torque_value_after, eff_pump_pm, Motor_power_value, energy_consumption_kJ):
    plt.figure()
    plt.subplot(2, 1, 1)
    plt.plot(time, cyc_lift_force, 'b-', linewidth=2)
    plt.title('Cyclical Lift Force')
    plt.xlabel('Time (s)')
    plt.ylabel('Force (N)')

    plt.subplot(2, 1, 2)
    plt.plot(time, cyc_tilt_force, 'r-', linewidth=2)
    plt.title('Cyclical Tilt Force')
    plt.xlabel('Time (s)')
    plt.ylabel('Force (N)')

    plt.figure()
    plt.subplot(2, 2, 1)
    plt.plot(time, d_l, 'k--', linewidth=2)
    plt.plot(time, d_ref_lift, 'r-', linewidth=2)
    plt.title('Displacement of Lift')
    plt.xlabel('Time (s)')
    plt.ylabel('Displacement (m)')
    plt.legend(['Actual', 'Reference'])

    plt.subplot(2, 2, 2)
    plt.plot(time, d_t, 'k--', linewidth=2)
    plt.plot(time, d_ref_tilt, 'r-', linewidth=2)
    plt.title('Displacement of Tilt')
    plt.xlabel('Time (s)')
    plt.ylabel('Displacement (m)')
    plt.legend(['Actual', 'Reference'])

    plt.subplot(2, 2, 3)
    plt.plot(time, control_inputs[2], 'g-', linewidth=2)
    plt.title('Valve Opening for Lift')
    plt.xlabel('Time (s)')
    plt.ylabel('Opening')

    plt.subplot(2, 2, 4)
    plt.plot(time, control_inputs[3], 'y-', linewidth=2)
    plt.title('Valve Opening for Tilt')
    plt.xlabel('Time (s)')
    plt.ylabel('Opening')

    plt.figure()
    plt.subplot(5, 1, 1)
    plt.plot(time, w_em_values, 'b-', linewidth=2)
    plt.title('Motor Speed (w_em)')
    plt.xlabel('Time (s)')
    plt.ylabel('Speed (rad/s)')

    plt.subplot(5, 1, 2)
    plt.plot(time, Torque_value_after, 'b-', linewidth=2)
    plt.title('Motor Torque')
    plt.xlabel('Time (s)')
    plt.ylabel('Torque (Nm)')

    plt.subplot(5, 1, 3)
    plt.plot(time, eff_pump_pm, 'b-', linewidth=2)
    plt.title('Pump Mechanical Efficiency')
    plt.xlabel('Time (s)')
    plt.ylabel('Efficiency')

    plt.subplot(5, 1, 4)
    plt.plot(time, Motor_power_value / 1000, 'b-', linewidth=2)
    plt.title('Power Consumption (kW)')
    plt.xlabel('Time (s)')
    plt.ylabel('Power (kW)')

    plt.subplot(5, 1, 5)
    plt.plot(time, energy_consumption_kJ, 'b-', linewidth=2)
    plt.title('Energy Consumption (kJ)')
    plt.xlabel('Time (s)')
    plt.ylabel('Energy (kJ)')

    plt.show()
