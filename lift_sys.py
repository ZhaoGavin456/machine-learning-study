import matplotlib.pyplot as plt
import numpy as np
import math

dt = 0.0001
# parameters for system
n_lift = 2  # number of lift piston;
m_l = 30 * n_lift  # ----------- mass of lift piston
k_pist = 5000  # piston spring coefficient
b_pist = 28000 / 1  # piston damping coefficient
lift_A = n_lift * (math.pi * pow(133.4e-3, 2) / 4)  # lift piston A area, m^2;
lift_B = n_lift * (math.pi * pow(133.4e-3, 2) / 4 - math.pi * pow(76.2e-3, 2) / 4)  # lift piston B area, m^2;
lift_ave = 0.5 * (lift_A + lift_B)  # average area
lift_delta = 0.5 * (lift_A - lift_B)  # area difference
A_p = lift_ave  # pressure area of the piston
A_spool = math.pi * 5e-3 * 1e-3  # ----------- spool area
beta = 1.8e9  # bulk modulus
storke_lift = 0.845  # piston stroke
V_t = storke_lift * lift_ave  # ----------- volume for lift chamber
C_d = 0.6  # discharge coefficient
rho = 840  # density
P_s = 35e6  # supply pump pressure

# some combination of parameters for system
a1 = k_pist / m_l
a2 = b_pist / m_l
a3 = A_p / m_l
h1 = beta * (pow(lift_A, 2) + pow(lift_B, 2)) / V_t / lift_ave
h2 = beta * (pow(lift_A, 2) + pow(lift_B, 2)) / V_t / lift_ave * C_d * A_spool * math.sqrt(
    2 / rho / (pow(lift_A, 3) + pow(lift_B, 3)))

# controller parameters
Kp = 10
Ki = 1


class State:

    def __init__(self, x=0.0, v=0.0, p_l=0.0):
        self.x = x
        self.v = v
        self.p_l = p_l

    def update(self, u_l, dt):
        self.x = self.x + self.v * dt
        self.v = self.v + (-a1 * self.x - a2 * self.v + a3 * self.p_l) * dt
        self.p_l = self.p_l + (-h1 * self.v + h2 * u_l * math.sqrt(
            lift_ave * P_s - np.sign(u_l) * (lift_delta * P_s - lift_ave * self.p_l))) * dt


class States:

    def __init__(self):
        self.x = []
        self.v = []
        self.p_l = []
        self.t = []
        self.ul = []

    def append(self, t, ul, state):
        self.x.append(state.x)
        self.v.append(state.v)
        self.p_l.append(state.p_l)
        self.t.append(t)
        self.ul.append(ul)




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
    T = 10.0  # max simulation time
    rt = np.arange(0, T, dt)
    rx = 0.02 * np.sin(rt)

    lastIndex = len(rx) - 1

    # initial state
    state = State(x=0.0, v=0.0, p_l=0.0)

    time = 0.0
    ind = -1
    integral_error = 0
    ul = 0
    states = States()
    states.append(time, ul, state)

    while time < T:
        # Calc control input

        ul, integral_error = pid_control(rx[ind], state.x, integral_error)

        state.update(ul, dt)

        time = time + dt
        states.append(time, ul, state)

        ind = ind + 1

    ######################## plot #######################
    fg1 = plt.figure(1)
    my_states_x = np.asarray(states.x)
    plt.plot(rt, rx)
    plt.plot(states.t, [ix for ix in states.x])
    plt.xlabel('time')
    plt.ylabel('lift position')
    plt.legend(['Ref', 'Actual'])

    fg2 = plt.figure(2)
    plt.plot(states.t, [iu for iu in states.ul])
    plt.xlabel('time')
    plt.ylabel('lift input')
    plt.legend(['input'])
    plt.show()


if __name__ == '__main__':
    main()
