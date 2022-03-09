import os


def n_outputs_one(i, subdir):
    path = os.path.join(".", str(i), subdir)
    return len([name for name in os.listdir(path) if
                os.path.isfile(os.path.join(path, name))])


def min_n_outputs(n_rep, subdir="inflammation"):
    len_all = [n_outputs_one(i, subdir) for i in range(n_rep)]
    return min(len_all)
