import numpy as np
import random as ran
import pandas as pd
from sklearn.decomposition import PCA

def ranRelPct(df, asLogOdds = True):
    def betaRow(row):
        ran_row = np.random.beta(row + 1, row.sum() - row + 1)
        return ran_row
    result = df.apply(lambda x: betaRow(x), axis = 1, result_type = 'expand') 
    result.index = df.index.values
    result.columns = df.columns.values
    if asLogOdds:
        result = np.log(result / (1 - result))
    return(result)

def samplePCA(df, num_pcs):
    ran_df = ranRelPct(df, asLogOdds = True)
    pca = PCA(n_components = num_pcs)
    pca_fit = pca.fit(ran_df)
    scores = pca_fit.transform(ran_df)
    loadings = np.transpose(pca_fit.components_)
    pca_results = {
        "df": ran_df,
        "scores": scores,
        "loadings": loadings
    }
    return(pca_results)

def harmonizeColumnSigns(mat_list):
    for i in range(1, len(mat_list)):
        for col in range(mat_list[i].shape[1]):
            if np.sign(mat_list[i][0, col]) != np.sign(mat_list[0][0, col]):
                mat_list[i][:, col] *= -1
    return(mat_list)