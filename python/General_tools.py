import numpy as np

def init_dict_stat(stat_name):
    """
    Init a dictionary with an empty numpy array
    """
    stat = {}
    for name in stat_name :
        stat[name] = np.array([])
    return stat

def create_position_array(x,y):
    """
    Take 2 array X and Y and return a tuple of position (X,Y)
    """
    pos = []
    for i in range(len(x)):
        pos.append( [x[i] , y[i]] ) 
    return pos