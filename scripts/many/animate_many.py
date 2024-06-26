import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib as mpl
import json
from count_outputs_since_inflammation import min_n_outputs

# n_rep = 600
# reps = range(n_rep)

reps = 138, 347, 68, 419
# reps are 0-indexed on disk but 1-indexed in R so decrement them
reps = [rep - 1 for rep in reps]

n_var = 7

n_time = min_n_outputs(reps)

var_names = ['$C_u$', '$C_b$', '$C_s$', '$\\phi_i$', '$\\phi_m$',
             '$\\phi_{{C_u}}$', '$\\phi_{{C_b}}$']

param_names = ['$J_{\\phi_i}^\\mathrm{I}$ (factor)',
               '$M^\\mathrm{I}$ (factor)',
               '$J_{\\phi_i}^\\mathrm{I}$ onset delay',
               '$\\gamma_{u,i}, \\gamma_{u,m}, \\dots$']

mpl.rc('mathtext', fontset='cm')

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


params_scaled = [get_params_scaled(rep) for rep in reps]

vars = [2, 7]
n_plot_var = len(vars)
n_vary_param = 4

fig, axs = plt.subplots(nrows=n_plot_var, ncols=n_vary_param, sharex=True)
plt.tight_layout()
plt.subplots_adjust(left=0.05, bottom=0.06, right=0.99, top=0.97,
                    wspace=0.2, hspace=0.075)

for ax in axs.flat:
    ax.margins(x=0)

fig.set_size_inches(16, 9, True)
time_annotation = fig.text(0.01, 0.98, 'Time = 0')


def data_single(rep, vars, time):
    try:
        return np.genfromtxt(f'{rep}/inflammation/output_{time:05d}.csv',
                             skip_header=1,
                             usecols=[1] + [(var + 1) for var in vars])

    except IOError:
        # This shouldn't happen, as we're meant to only be going up to the end
        # time of the simulation that finished first
        print("WARNING: TRYING TO READ A FILE THAT DOESN'T EXIST")


param_max_colours = [[0.0, 0.0, 1.0],
                     [1.0, 0.0, 0.0],
                     [0.0, 1.0, 0.0],
                     [1.0, 0.647, 0.0]]


def colour_from_param(param_scaled, max_colour):
    colour = mpl.colors.rgb_to_hsv(max_colour)
    colour[2] *= param_scaled
    return mpl.colors.hsv_to_rgb(colour)


def plot_single(idx, rep, vars, time):
    data = data_single(rep, vars, time)

    lines = [[None for j in range(n_vary_param)] for i in range(n_plot_var)]

    if data is not None:
        for i in range(n_plot_var):
            for j in range(n_vary_param):
                (r, g, b) = colour_from_param(params_scaled[idx][j],
                                              param_max_colours[j])

                line, = axs[i, j].plot(data[:, 0], data[:, i + 1],
                                       color=(r, g, b, 0.4), linewidth=1.0)
                lines[i][j] = line

        return lines


for i in range(n_plot_var):
    axs[i, 0].set_ylabel(var_names[vars[i] - 1], fontsize=18)

for j in range(n_vary_param):
    axs[n_plot_var - 1, j].set_xlabel(r"$x$", fontsize=18)
    axs[0, j].set_title(param_names[j])

lines = [plot_single(idx, rep, vars, 0) for idx, rep in enumerate(reps)]
data = [None for rep in reps]

use_fixed_y_axis_limits = True
fixed_y_axis_max = [9.0, 1.8]


def animate(time):
    print(f"vars: {vars}, time: {time} / {n_time}")

    data_max = [0.0] * n_plot_var
    for idx, rep in enumerate(reps):
        data[idx] = data_single(rep, vars, time)

        if not use_fixed_y_axis_limits:
            data_rep_max = \
                [data[idx][:, i + 1].max() for i in range(n_plot_var)]

            for i in range(n_plot_var):
                if data_rep_max[i] > data_max[i]:
                    data_max[i] = data_rep_max[i]

    if use_fixed_y_axis_limits:
        y_axis_max = fixed_y_axis_max
    else:
        y_axis_max = data_max

    for i in range(n_plot_var):
        for j in range(n_vary_param):
            axs[i, j].set_ylim(0.0, 1.01*y_axis_max[i])

    for idx, rep in enumerate(reps):
        if data[idx] is not None:
            for i in range(n_plot_var):
                for j in range(n_vary_param):
                    (r, g, b) = colour_from_param(params_scaled[idx][j],
                                                  param_max_colours[j])

                    lines[idx][i][j].set_ydata(data[idx][:, i + 1])
                    lines[idx][i][j].set_color((r, g, b, 0.4))

    time_annotation.set_text(f'Time = {time}')

    return lines


ani = animation.FuncAnimation(fig, animate, blit=False, interval=1000.0/60.0,
                              frames=range(0, n_time))

writer = animation.FFMpegWriter(fps=60)
ani.save("videos/grid.mp4", dpi=160)

# plt.show()
