import numpy as np
import random as ran

def ranRelPct(mat):
    def drawOneSample(x, total):
        a = x + 1
        return ran.betavariate(a, total - a)
    p = [[drawOneSample(count, sum(row)) for count in row] for row in mat]
    return(np.array(p))