import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import json
import pandas as pd

n_rep = 100
n_var = 7
n_time = 1501
# first_inflammation_time = 544
first_inflammation_time = 0

var_names = ['$C_u$', '$C_b$', '$C_s$', '$\\phi_i$', '$\\phi_m$',
             '$\\phi_{{C_u}}$', '$\\phi_{{C_b}}$']

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
plt.subplots_adjust(bottom=0.10)

ax_var_slider = plt.axes([0.1, 0.05, 0.8, 0.015])
ax_time_slider = plt.axes([0.1, 0.025, 0.8, 0.015])

var_slider = Slider(ax_var_slider, 'Variable', valmin=1, valmax=n_var,
                    valstep=variables)

time_slider = Slider(ax_time_slider, 'Time', valmin=0, valmax=(n_time - 1),
                     valstep=times)

first_data_for_header = pd.read_csv('0/output_00000.csv', sep=' ',
                                    usecols=range(1, n_var + 2))
col_names = [n.replace('\\\\', '\\') for n in first_data_for_header.columns]


def data_single(rep, var, time):
    try:
        return np.genfromtxt(f'{rep}/output_{time:05d}.csv', skip_header=1,
                             usecols=(1, var + 1))

    except IOError:
        # ignore any io errors - they due to some simulations finishing earlier
        # than the maximum time
        return None


def plot_single(rep, var, time, solid=False):
    data = data_single(rep, var, time)

    ax.set_title(col_names[var])
    ax.set_xlabel(r'$x$')

    if data is not None:
        if solid:
            ax.plot(data[:, 0], data[:, 1], color=(0.0, 0.0, 0.0, 1.0),
                    linewidth=1)
        else:
            j_phi_i_i_factor, m_i_factor, t_j_phi_i_lag, gamma = \
                params_scaled[rep]
            r = 0.9 * j_phi_i_i_factor
            g = 0.9 * m_i_factor
            b = 0.9 * t_j_phi_i_lag
            # r = 0.0
            # g = 0.0
            # b = gamma

            ax.plot(data[:, 0], data[:, 1], color=(r, g, b, 1.0), linewidth=1)


def update(val):
    ax.cla()

    var = var_slider.val
    time = time_slider.val

    if time < first_inflammation_time:
        plot_single('homeostasis', var, time, solid=True)
    else:
        for rep in range(0, n_rep):
            plot_single(rep, var, time)


var_slider.on_changed(update)
time_slider.on_changed(update)

update(0)

plt.show()
