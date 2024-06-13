import matplotlib.pyplot as plt
import numpy as np
import math
# time step
dt = 0.01
# parameters for system
r = 1.5
R_wh = 0.761
L = 3.34
kf = 24.67
Je = 8
mu = 0.0286
g = 9.8
eta = 0.55
n_cyl = 6
q_hv = 42.9e+06  # J/kg
C_fr1 = -3.8e-4
C_fr2 = 0.1155
C_fr3 = -4.8766
eta_ig = 0.4722
Vd = 6.7e-3  # m^3
rou_material = 1500  # kg/m^3
L_bucket = 2.95
R = 1.5
a1 = 1.7259e-04
a2 = 9.9186e-05
a3 = -1.6397e-04
b1 = 4.3442e-04
b2 = -3.1468e-04
b3 = -9.5534e-05
Tbp = 50
GR = 2.184
theta = 10 / 180 * math.pi
w_f = L_bucket
alpha_f = (4.23 / 12) * math.pi
ground_x = 4
ground_y = 0
gamma_f = 1500
# trained parameters for FEE model
c_f = 3.5020e+03
c_af = 1e+04
phi_f = 0.0198
delta_f = 0.4154
k_cf = 3.9762e+03
k_phif = 4.7337e+03
n_f = 1

# controller parameters
Kp = 50
Ki = 10


class State:

    def __init__(self, phi=0.0, omega_t=0.0, v_l=0.0, x_b=0.0, z_b=0.0, omega_e=90.0, area=0.0):
        self.phi = phi
        self.omega_t = omega_t
        self.v_l = v_l
        self.x_b = x_b
        self.z_b = z_b
        self.omega_e = omega_e
        self.area = area

    def update(self, u_t, u_p, u_f, dt, m_load, F_H, F_N):
        # u_t, u_p, u_f are control inputs
        # u_t is angular acceleration, u_p is linear acceleration, u_f is the injected fuel per cycle.
        # m_load is loaded material and bucket mass
        # F_H , F_N are resistance forces, F_H is tangential force, F_N is normal force.
        # bucket orientation angle
        self.phi = self.phi + self.omega_t * dt
        # bucket orientation angle constraints
        if self.phi < math.pi / 2:
            self.phi = math.pi / 2
        if self.phi > math.pi:
            self.phi = math.pi
        # rotational angular velocity
        self.omega_t = self.omega_t + u_t * dt
        # rotational angular velocity constraints
        if self.omega_t < -1:
            self.omega_t = -1
        if self.omega_t > 1:
            self.omega_t = 1
        # translational velocity vector of the bucket
        self.v_l = self.v_l + u_p * dt
        # translational velocity vector of the bucket constraints
        if self.v_l < 0:
            self.v_l = 0
        if self.v_l > 1:
            self.v_l = 1
        # bucket tip x position
        self.x_b = self.x_b + (
                (self.v_l + r * self.omega_t * math.cos(theta)) * math.cos(self.phi) + r * self.omega_t * math.sin(
            theta) * math.sin(math.pi - self.phi)) * dt
        # bucket tip x position constraints
        if self.x_b < 2:
            self.x_b = 2
        if self.x_b > 4:
            self.x_b = 4
        # bucket tip z position
        self.z_b = self.z_b + (
                (self.v_l + r * self.omega_t * math.cos(theta)) * math.sin(self.phi) + r * self.omega_t * math.sin(
            theta) * math.cos(math.pi - self.phi)) * dt
        # bucket tip z position constraints
        if self.z_b < 0:
            self.z_b = 0
        if self.z_b > 3:
            self.z_b = 3
        # engine speed
        Tig = (eta_ig * q_hv * n_cyl * u_f * 0.000001 / 4 / math.pi)
        Tfric = Vd * (100000) * (C_fr1 * pow(self.omega_e, 2) + C_fr2 * self.omega_e + C_fr3) / 4 / math.pi
        Te = Tig - Tfric
        # engine torque constraints
        if Te < 30:
            Te = 30
        if Te > 900:
            Te = 900
        Tp = pow((60 / 2 / math.pi), 2) * (a1 * pow(self.omega_e, 2))
        Tw = (F_H + m_load * g * (math.sin(math.pi - self.phi)) * (self.v_l + r * self.omega_t * math.cos(theta))
              + (F_N + m_load * g * (math.cos(math.pi - self.phi))) * (r * self.omega_t * math.sin(theta))) / (
                     eta * self.omega_e)
        # hydraulic pump torque constraints
        if Tw < 0:
            Tw = 0
        if Tw > 900:
            Tw = 900
        self.omega_e = self.omega_e + ((1 / Je) * (Te - Tp - Tw - Tbp)) * dt
        # engine speed constraints
        if self.omega_e < 90:
            self.omega_e = 90
        if self.omega_e > 220:
            self.omega_e = 220
        # cross-sectional area
        self.area = self.area + (
                ((-3 / 4 * pow(self.x_b, 2) + 3 * self.x_b) - self.z_b) * (self.v_l + r * self.omega_t) * (
            -math.cos(self.phi))) * dt
        # cross-sectional area constraints
        if self.area < 0:
            self.area = 0
        if self.area > 2:
            self.area = 2


class States:

    def __init__(self):
        self.phi = []
        self.omega_t = []
        self.v_l = []
        self.x_b = []
        self.z_b = []
        self.omega_e = []
        self.area = []

    def append(self, t, u_t, u_p, u_f, state):
        self.phi.append(state.phi)
        self.omega_t.append(state.omega_t)
        self.v_l.append(state.v_l)
        self.x_b.append(state.x_b)
        self.z_b.append(state.z_b)
        self.omega_e.append(state.omega_e)
        self.area.append(state.area)
        self.t.append(t)
        self.u_t.append(u_t)
        self.u_p.append(u_p)
        self.u_f.append(u_f)


def load_mass_cal(state):
    m_load = rou_material * state.area * L_bucket
    return m_load


def resistance_force(state, m_load):
    rho_f = alpha_f - math.atan((state.z_b) / (ground_x - state.x_b))  # atan vs arctan vs cot
    beta_f = math.pi / 2 - alpha_f
    dist_f = math.sqrt(pow((state.x_b - ground_x), 2) + pow((state.z_b - ground_y), 2))
    if dist_f > 1.44:
        dist_f = 1.44
    d_f = dist_f * math.sin(rho_f)
    denom = math.cos(rho_f + delta_f) + math.sin(rho_f + delta_f) * math.cot(beta_f + phi_f)
    N_gamma = (math.cot(beta_f) - math.tan(alpha_f)) * (
            math.cos(alpha_f) + math.sin(alpha_f) * math.cot(beta_f + phi_f)) / 2 / denom
    N_c = (1 + math.cot(beta_f) * math.cot(beta_f + phi_f)) / denom
    N_a = (1 - (math.cos(rho_f) / math.sin(rho_f)) * math.cot(beta_f + phi_f)) / denom
    N_q = (math.cos(alpha_f) + math.sin(alpha_f) * math.cot(beta_f + phi_f)) / denom
    F = pow(d_f, 2) * w_f * gamma_f * g * N_gamma + c_f * w_f * d_f * N_c + c_af * w_f * d_f * N_a + m_load * g * N_q
    F_B = (k_cf / 0.06 + k_phif) * pow(d_f, n_f)
    F_H = F_B + F * math.sin(delta_f) + c_af * dist_f * w_f
    F_N = F * math.cos(delta_f)
    return F_H, F_N


def pid_control(target, current):
    # p control

    u_t = Kp * (target - current)
    u_p = Kp * (target - current)
    u_f = Kp * (target - current)

    # control input constraints
    u_max = 1
    # u_t constraints
    if u_t > u_max:
        u_t = u_max
    if u_t < -u_max:
        u_t = -u_max
    # u_p constraints
    if u_p > u_max:
        u_p = u_max
    if u_p < -u_max:
        u_p = -u_max
    # u_f constraints
    if u_f < 15:
        u_f = 15
    if u_f > 115:
        u_f = 115

    return u_t, u_p, u_f


def main():
    #  target course
    T = 10.0  # max simulation time
    rt = np.arange(0, T, dt)
    rx = 0.02 * np.sin(rt)

    lastIndex = len(rx) - 1

    # initial state
    state = State(phi=0.0, omega_t=0.0, v_l=0.0, x_b=0.0, z_b=0.0, omega_e=90.0, area=0.0)

    time = 0.0
    ind = -1
    u_t = 0
    u_p = 0
    u_f = 15
    states = States()
    states.append(time, u_t, u_p, u_f, state)

    while time < T:
        # Calc control input
        u_t, u_p, u_f = pid_control(rx[ind], state.x_b)
        m_load = load_mass_cal(state)
        F_H, F_N = resistance_force(state, m_load)
        state.update(u_t, u_p, u_f, dt, m_load, F_H, F_N)
        time = time + dt
        states.append(time, u_t, u_p, u_f, state)
        ind = ind + 1

    ######################## plot #######################
    plt.figure()
    my_states_x_b = np.asarray(states.x_b)
    plt.plot(rt, rx)
    plt.plot(states.t, [ix for ix in states.x_b])
    plt.xlabel('time')
    plt.ylabel('lift position')
    plt.legend(['Ref', 'Actual'])
    plt.show()

    plt.figure()
    plt.plot(states.t, [iu for iu in states.u_f])
    plt.xlabel('time')
    plt.ylabel('lift input')
    plt.legend(['input'])
    plt.show()


if __name__ == '__main__':
    main()
