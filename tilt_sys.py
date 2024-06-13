import matplotlib.pyplot as plt
import numpy as np
import math

dt = 0.0001
# parameters for system
n_tilt = 1 # number of tilt cylinder;
m_t = 30*n_tilt
k_pist = 5000 # piston spring coefficient  
b_pist = 28000 # piston damping coefficient  
tilt_A = 1 * (math.pi*pow(165.1e-3,2)/4) # tilt piston A area, m^2;
tilt_B = 1 * (math.pi*pow(165.1e-3,2)/4 - math.pi*pow(88.9e-3,2)/4) # tilt piston B area, m^2;
tilt_ave = 0.5*(tilt_A+tilt_B) # average area
tilt_delta = 0.5*(tilt_A-tilt_B)  # area difference
A_p = tilt_ave # pressure area of the piston
A_spool = math.pi*5e-3*1e-3 # spool area
storke_tilt = 0.6 # piston stroke
beta = 1.8e9 # butk modutus
V_t = storke_tilt*tilt_ave # volume for tilt chamber
C_d = 0.6 # discharge coefficient
rho = 840 # density
P_s = 35e6 # supply pump pressure


# some combination of parameters for system
a1 = k_pist / m_t
a2 = b_pist / m_t
a3 = A_p / m_t
h1 = beta * (pow(tilt_A, 2) + pow(tilt_B, 2)) / V_t / tilt_ave
h2 = beta * (pow(tilt_A, 2) + pow(tilt_B, 2)) / V_t / tilt_ave * C_d * A_spool * math.sqrt(
    2 / rho / (pow(tilt_A, 3) + pow(tilt_B, 3)))

# controller parameters
Kp = 10
Ki = 1


class State:

    def __init__(self, x=0.0, v=0.0, p_t=0.0):
        self.x = x
        self.v = v
        self.p_t = p_t

    def update(self, u_t, dt):
        self.x = self.x + self.v * dt
        self.v = self.v + (-a1 * self.x - a2 * self.v + a3 * self.p_t) * dt
        self.p_t = self.p_t + (-h1 * self.v + h2 * u_t * math.sqrt(
            tilt_ave * P_s - np.sign(u_t) * (tilt_delta * P_s - tilt_ave * self.p_t))) * dt


class States:

    def __init__(self):
        self.x = []
        self.v = []
        self.p_t = []
        self.t = []
        self.ut = []

    def append(self, t, ut, state):
        self.x.append(state.x)
        self.v.append(state.v)
        self.p_t.append(state.p_t)
        self.t.append(t)
        self.ut.append(ut)




def pid_control(target, current, integral_error):
    # p control
    a_p = Kp * (target - current)
    # i control
    integral_error += (target - current)
    a_i = Ki * integral_error

    a = a_p + a_i

    # control input limits
    a_max = 1
    if a > a_max:
        a = a_max
    if a < -a_max:
        a = -a_max
    return a, integral_error


def main():
    #  target course
    T = 10.0  # max simutation time
    rt = np.arange(0, T, dt)
    rx = 0.02 * np.sin(rt)

    lastIndex = len(rx) - 1

    # initial state
    state = State(x=0.0, v=0.0, p_t=0.0)

    time = 0.0
    ind = -1
    integral_error = 0
    ut = 0
    states = States()
    states.append(time, ut, state)

    while time < T:
        # Calc control input

        ut, integral_error = pid_control(rx[ind], state.x, integral_error)

        state.update(ut, dt)

        time = time + dt
        states.append(time, ut, state)

        ind = ind + 1

    ######################## plot #######################
    plt.figure()
    my_states_x = np.asarray(states.x)
    plt.plot(rt, rx)
    plt.plot(states.t, [ix for ix in states.x])
    plt.xlabel('time')
    plt.ylabel('tilt position')
    plt.legend(['Ref', 'Actual'])

    plt.figure()
    plt.plot(states.t, [iu for iu in states.ut])
    plt.xlabel('time')
    plt.ylabel('tilt input')
    plt.legend(['input'])
    plt.show()


if __name__ == '__main__':
    main()
