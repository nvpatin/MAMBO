import numpy as np
import random as ran
import pandas as pd
from sklearn.decomposition import PCA
from scipy.stats import rankdata


# Return data frame of a draw of relative percent of occurrence from a beta distribution
# fit to observed occurrence counts
#   df: data frame where rows = samples and columns = ASVs
def ranRelPct(df, asLogOdds = True):
    # function to return a random draw from a beta distribution for a row
    def betaRow(row):
        ran_row = np.random.beta(row + 1, row.sum() - row + 1)
        ran_row = ran_row / ran_row.sum()
        return ran_row
    # apply function to every row to draw sample of relative percent occurrence
    result = df.apply(betaRow, axis = 1, result_type = 'expand')
    # assign row and column names
    result.index = df.index.values
    result.columns = df.columns.values
    # convert to log-odds if requested
    if asLogOdds:
        result = np.log(result / (1 - result))
    return result


# Return data frame of a draw of relative percent of occurrence from a beta distribution
# fit to observed occurrence counts
#   df: data frame where rows = ASVs and columns = samples
def ranRelPct_cupy(df, asLogOdds = True):
    for i in range(df.shape[1]):
        col = df.iloc[:,i]
        a = col + 1
        b = col.sum() - col + 1
        beta_dist = cupy.random.beta(a,b)
        beta_dist /= col.sum()
        df.iloc[:,i] = beta_dist
    # convert to log-odds if requested
    if asLogOdds:
        result = np.log(df / (1 - df))
    return df.transpose()


# Draws sample of relative percent of occurrence and conducts PCA
#   df: data frame where rows = samples and columns = ASVs
# Returns a dictionary containing :
#   df: data frame containing a random draw of the original count data frame as log-odds
#   scores: array of PCA scores
#   loadings: array of PCA loadings
def samplePCA(df, num_pcs = None):
    max_pcs = min(df.shape[0] - 1, df.shape[1] - 1)
    if num_pcs is None:
        num_pcs = max_pcs
    elif num_pcs > max_pcs:
        num_pcs = max_pcs
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
    return pca_results


# If the sign of the first element in a column in matrices after the first is different than the first,
# multiply all values in that column by -1
def harmonizeColumnSigns(mat_list):
    for i in range(1, len(mat_list)):
        for col in range(mat_list[i].shape[1]):
            if np.sign(mat_list[i][0, col]) != np.sign(mat_list[0][0, col]):
                mat_list[i][:, col] *= -1
    return mat_list


# Sorts PCA loadings from a list 
def sortLoadings(loading_list, pc, asvs, asRanks = False):
    # Harmonize signs across replicates
    harm_loadings = harmonizeColumnSigns(loading_list)
    # Create 3 dimensional array and select component 'pc'
    loadings = np.stack(harm_loadings, axis = 2)[:, pc, :]
    # Convert to ranks if 'asRanks == True'
    if asRanks:
        loadings = np.array([rankdata(loadings[:, i]) for i in range(loadings.shape[1])]).transpose()
    # Get sorted order based on median for each ASV 
    row_sort = np.apply_along_axis(np.median, 1, loadings).ravel().argsort()[::-1]
    # Sort based on median, add ASV names (also sorted) and return data frame
    df = pd.DataFrame(loadings[row_sort, :])
    df.index = asvs[row_sort]
    return df