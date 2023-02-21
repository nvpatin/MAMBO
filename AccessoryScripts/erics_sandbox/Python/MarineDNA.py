import numpy as np
import random as ran

def ranRelPct(mat):
    def drawOneSample(x, total):
        a = x + 1
        return ran.betavariate(a, total - a)
    p = [[drawOneSample(count, sum(row)) for count in row] for row in mat]
    return(np.array(p))

def harmonizeColumnSigns(mat_list):
    for i in range(1, len(mat_list)):
        for col in range(mat_list[i].shape[1]):
            if np.sign(mat_list[i][0, col]) != np.sign(mat_list[0][0, col]):
                mat_list[i][:, col] *= -1
    return(mat_list)