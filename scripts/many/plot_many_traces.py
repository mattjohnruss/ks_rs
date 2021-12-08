import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import json

n_rep = 16
n_var = 10

j_phi_i_i_factor_min = 1.0
j_phi_i_i_factor_max = 1000.0
m_i_factor_min = 1.0
m_i_factor_max = 1000.0
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

file_list = [f'{rep}/trace.csv' for rep in range(0, n_rep)]

# Read data into a bunch of DataFrames rather than one big numpy array because
# some aren't full size (as the simulations finish early). There's probably a
# better way to do this but it works...
dataframes = [pd.read_csv(file_path, sep=' ') for file_path in file_list]
data = pd.concat(dataframes, keys=range(0, n_rep), axis=1)


mpl.rc('mathtext', fontset='cm')

for var in range(1, n_var + 1):
    fig, ax = plt.subplots()
    ax.yaxis.set_major_formatter(mpl.ticker.StrMethodFormatter("{x:.2f}"))

    fig.set_figwidth(15)
    fig.set_figheight(10)

    ax.set_xlabel(r'$t$')

    # set the plot titles from the column names from the first data series -
    # they will be the same as all other series
    ax.set_title(data[0].columns[var])

    for rep in range(0, n_rep):
        j_phi_i_i_factor, m_i_factor, t_j_phi_i_lag, gamma = params_scaled[rep]
        # r = 0.9 * j_phi_i_i_factor
        # g = 0.9 * m_i_factor
        # b = 0.9 * t_j_phi_i_lag
        r = 0.0
        g = 0.0
        b = gamma

        ax.plot(data[rep]['t'], data[rep].iloc[:, var],
                color=(r, g, b, 0.8), linewidth=0.5)

    fig.savefig(f'plots/{var}.pdf', bbox_inches='tight')
