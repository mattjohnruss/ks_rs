import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from glob import glob
from sys import argv

file_list = np.sort(glob("./output_*.csv"))
n_files = len(file_list)
n_points = int(argv[1])
n_var = int(argv[2])
n_timesteps = int(argv[3])
dt = float(argv[4])

# initialise empty storage
data = np.zeros([n_files, n_points, n_var + 2])

# loop over files and read data
for (i, file_path) in enumerate(file_list):
    data[i] = np.genfromtxt(file_path, skip_header=1, usecols=(0, 1, 2, 3, 4,
                                                               5, 6, 7, 8))

# extract the times and positions from the data, and then remove them
t = data[:, 0, 0]
x = data[0, :, 1]

# remove the times and positions from the data
data = data[:, :, 2:]

aspect_ratio = x.max()/t.max()

mpl.rc('mathtext', fontset='cm')

var_names = ['$C_u$', '$C_b$', '$C_s$', '$\\phi_i$', '$\\phi_m$',
             '$\\phi_{{C_u}}$', '$\\phi_{{C_b}}$']
simple_var_names = ['c_u', 'c_b', 'c_s', 'phi_i', 'phi_m', 'phi_c_u',
                    'phi_c_b']

for (i, var_name) in enumerate(var_names):
    fig = plt.figure()

    ax = fig.add_subplot(111)

    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$t$')

    # hack to reduce colobar range of phi_c_b
    # vmax = None
    # if i == 6:
    #     vmax = 20

    # im = ax.imshow(data[0:n_timesteps, :, i], interpolation='none',
    #                origin='lower', extent=[x.min(), x.max(), t.min(),
    #                t.max()],
    #                vmax=vmax, aspect=aspect_ratio, cmap='inferno')

    im = ax.imshow(data[0:n_timesteps, :, i], interpolation='none',
                   origin='lower', extent=[x.min(), x.max(), t.min(), t.max()],
                   aspect=aspect_ratio, cmap='inferno')

    # divider = make_axes_locatable(ax)
    # cax = divider.append_axes("right", size="5%", pad=0.01)
    # plt.colorbar(im, cax=cax)

    # hack to use 2 digits after decimal points, except for phi_c_b, to make
    # all the labels the same resulting length
    # format = '%.2f'
    # if i == 6:
    #     format = '%.1f'

    # plt.colorbar(im, format=format)
    plt.colorbar(im)

    ax.set_title(var_name)

    plt.savefig(f'{simple_var_names[i]}.pdf', bbox_inches='tight')
