import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib as mpl
import json
import sys

n_rep = 100
n_var = 7
n_time = 5001
# first_inflammation_time = 544
first_inflammation_time = 0

var_names = ['$C_u$', '$C_b$', '$C_s$', '$\\phi_i$', '$\\phi_m$',
             '$\\phi_{{C_u}}$', '$\\phi_{{C_b}}$']

mpl.rc('mathtext', fontset='cm')

variables = [i for i in range(1, n_var + 1)]
times = [i for i in range(0, n_time)]

m_i_factor_min = 1.0
m_i_factor_max = 1000.0
j_phi_i_i_factor_min = 1.0
j_phi_i_i_factor_max = 1000.0
t_j_phi_i_lag_min = 0.0
t_j_phi_i_lag_max = 25.0
gamma_min = 0.0
gamma_max = 10.0


def get_params_scaled(rep):
    with open(f'{rep}/config.json', 'r') as config_file:
        data = json.load(config_file)
        return [(data['j_phi_i_i_factor'] - j_phi_i_i_factor_min) /
                (j_phi_i_i_factor_max - j_phi_i_i_factor_min),
                (data['m_i_factor'] - m_i_factor_min) /
                (m_i_factor_max - m_i_factor_min),
                (data['t_j_phi_i_lag'] - t_j_phi_i_lag_min) /
                (t_j_phi_i_lag_max - t_j_phi_i_lag_min),
                (data['gamma_ui'] - gamma_min) /
                (gamma_max - gamma_min)]


params_scaled = [get_params_scaled(rep) for rep in range(0, n_rep)]

fig, ax = plt.subplots()
plt.tight_layout()

fig.set_size_inches(16, 9, True)


def data_single(rep, var, time):
    try:
        if time < first_inflammation_time:
            rep = 'homeostasis'

        return np.genfromtxt(f'{rep}/output_{time:05d}.csv', skip_header=1,
                             usecols=(1, var + 1))

    except IOError:
        # ignore any io errors - they due to some simulations finishing earlier
        # than the maximum time
        return None


def plot_single(rep, var, time, solid=False):
    data = data_single(rep, var, time)

    if data is not None:
        line = None
        if solid:
            line, = ax.plot(data[:, 0], data[:, 1], color=(0.0, 0.0, 0.0, 1.0),
                            linewidth=1)
        else:
            j_phi_i_i_factor, m_i_factor, t_j_phi_i_lag, gamma = \
                params_scaled[rep]
            # r = 0.9 * j_phi_i_i_factor
            # g = 0.9 * m_i_factor
            # b = 0.9 * t_j_phi_i_lag
            r = 0.0
            g = 0.0
            b = gamma

            line, = ax.plot(data[:, 0], data[:, 1], color=(r, g, b, 1.0),
                            linewidth=1)

        return line


var = int(sys.argv[1])

ax.set_title(var_names[var - 1], fontsize=25)
ax.set_xlabel(r"$x$", fontsize=18)
lines = [plot_single(rep, var, 0) for rep in range(0, n_rep)]
data = [None for rep in range(0, n_rep)]


def animate(time):
    print(f"var: {var}, time: {time}")

    data_max = 0.0
    for rep in range(0, n_rep):
        data[rep] = data_single(rep, var, time)
        data_rep_max = data[rep][:, 1].max()
        if data_rep_max > data_max:
            data_max = data_rep_max

    ax.set_ylim(0.0, 1.01*data_max)

    for rep in range(0, n_rep):
        if data[rep] is not None:
            lines[rep].set_ydata(data[rep][:, 1])

            r, g, b = 0.0, 0.0, 0.0

            if time >= first_inflammation_time:
                j_phi_i_i_factor, m_i_factor, t_j_phi_i_lag, gamma = \
                    params_scaled[rep]
                # r = 0.9 * j_phi_i_i_factor
                # g = 0.9 * m_i_factor
                # b = 0.9 * t_j_phi_i_lag
                r = 0.0
                g = 0.0
                b = gamma

            lines[rep].set_color((r, g, b, 1.0))

    return lines


ani = animation.FuncAnimation(fig, animate, blit=False, interval=1000.0/60.0,
                              frames=range(0, n_time))

writer = animation.FFMpegWriter(fps=60)
ani.save(f"videos/{var}.mp4", dpi=160)

# plt.show()
