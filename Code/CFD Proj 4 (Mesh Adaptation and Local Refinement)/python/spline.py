import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import splprep, splev

def snap2splines(point,airfoil):
    distances = np.sqrt(np.sum((airfoil - [point[0], point[1]])**2, axis=1))
    min_index = np.argsort(distances)[:2]

    if min(min_index) == 0 or min(min_index) == 201 or min(min_index) == 402:
        interval = np.array([airfoil[min(min_index)],airfoil[min(min_index)+1],airfoil[min(min_index)+2],airfoil[min(min_index)+3]])
    elif max(min_index) == 200 or max(min_index) == 401 or max(min_index) == 602:
        interval = np.array([airfoil[max(min_index)-3],airfoil[max(min_index)-2],airfoil[max(min_index)-1],airfoil[max(min_index)]])
    else:
        interval = np.array([airfoil[min(min_index)-1],airfoil[min(min_index)],airfoil[max(min_index)],airfoil[max(min_index)+1]])

    N = 10
    for i in range(5):
        tck, u = splprep(interval.T, u=None, s=0.0) 
        u_new = np.linspace(u.min(), u.max(), N)
        x_new, y_new = splev(u_new, tck, der=0)
        cstack = np.column_stack((x_new, y_new))
        distances = np.sqrt(np.sum((cstack - [point[0], point[1]])**2, axis=1))
        min_index = np.argsort(distances)[:2]
        if min(min_index) == 0:
            s = 1
        elif max(min_index) == N-1:
            s = -1
        else: 
            s = 0
        interval = np.array([cstack[min(min_index)-1+s],cstack[min(min_index)+s],cstack[max(min_index)+s],cstack[max(min_index)+1+s]])

    return interval[len(interval)//2]