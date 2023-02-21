import numpy as np
import random as ran
import pandas as pd

def ranRelPct(df):
    def betaRow(row):
      ran_row = np.random.beta(row + 1, row.sum() - row + 1)
      return ran_row
    result = df.apply(lambda x: betaRow(x), axis = 1, result_type = 'expand') 
    result.index = df.index.values
    result.columns = df.columns.values
    return(result)

def harmonizeColumnSigns(mat_list):
    for i in range(1, len(mat_list)):
        for col in range(mat_list[i].shape[1]):
            if np.sign(mat_list[i][0, col]) != np.sign(mat_list[0][0, col]):
                mat_list[i][:, col] *= -1
    return(mat_list)