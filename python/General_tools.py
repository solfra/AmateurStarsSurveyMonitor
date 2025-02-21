import numpy as np

def init_dict_stat(stat_name):
    stat = {}
    for name in stat_name :
        stat[name] = np.array([])
    return stat