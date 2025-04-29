import numpy as np

def init_dict_stat(stat_name):
    stat = {}
    for name in stat_name :
        stat[name] = np.array([])
    return stat

def create_position_array(x,y):
    pos = []
    for i in range(len(x)):
        pos.append( [x[i] , y[i]] ) 
    return pos