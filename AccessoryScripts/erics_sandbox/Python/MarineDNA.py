# Return 3-d array of Beta alpha and beta parameters for a data frame of raw read counts
# Output dimensions are (parameters, ASVs, samples)
def betaParams(raw_counts):
    def colParams(col):
        return (col + 1, col.sum() - col + 1)
    return np.stack(
        [colParams(raw_counts.iloc[:, i]) for i in range(raw_counts.shape[1])], 
        axis = 2,
        dtype = np.dtype(np.int64, metadata = {"asvs": asvs1.index, "samples": asvs1.columns})
    )
    
# Return random draw from Beta distribution of read percentages.
# Output is a data frame with rows = samples ad columns = ASVs
def ranRelPct_new(beta_params, asLogOdds = True):
    # function to draw from distribution and set percentages to sum to 1 for each sample
    def betaCol(col):
        beta_dist = np.random.beta(col[0, :], col[1, :])
        return beta_dist / beta_dist.sum()
    # pre-allocate result array (NOTE: transposes data structure)
    result = np.empty([beta_params.shape[2], beta_params.shape[1]])
    # draw from for each sample
    for i in range(result.shape[0]):
        result[i, :] = betaCol(beta_params[:, :, i])
    if asLogOdds:
        result = np.log(result / (1 - result))
    return pd.DataFrame(result, index = beta_params.dtype.metadata["samples"], columns = beta_params.dtype.metadata["asvs"])


# Return data frame of a draw of relative percent of occurrence from a beta distribution
# fit to observed occurrence counts
#   df: data frame where rows = ASVs and columns = samples
def ranRelPct(df, asLogOdds = True):
    import numpy as np
    
    result = df.copy()
    for i in range(df.shape[1]):
        col = df.iloc[:,i]
        a = col + 1
        b = col.sum() - col + 1
        beta_dist = np.random.beta(a,b)
        beta_dist /= beta_dist.sum()
        result.iloc[:,i] = beta_dist
    # convert to log-odds if requested
    if asLogOdds:
        result = np.log(result / (1 - result))
    return result.transpose()


# Does PCA
#   df: data frame where rows = samples and columns = ASVs
#   num_pcs: number of components to return. if None, return maximum number
# Returns a dictionary containing :
#   scores: array of PCA scores
#   loadings: array of PCA loadings
def doPCA(df, num_pcs = None):
    import numpy as np
    from sklearn.decomposition import PCA
    
    max_pcs = min(df.shape[0] - 1, df.shape[1] - 1)
    if num_pcs is None:
        num_pcs = max_pcs
    elif num_pcs > max_pcs:
        num_pcs = max_pcs
    pca = PCA(n_components = num_pcs)
    pca_fit = pca.fit(df)
    pca_results = {
        "scores": pca_fit.transform(df),
        "loadings": np.transpose(pca_fit.components_)
    }
    return pca_results


# If the sign of the first element in a column in matrices after the first is different than the first,
# multiply all values in that column by -1
def harmonizeColumnSigns(mat_list):
    for i in range(1, len(mat_list)):
        for col in range(mat_list[i].shape[1]):
            if cp.sign(mat_list[i][0, col]) != cp.sign(mat_list[0][0, col]):
                mat_list[i][:, col] *= -1
    return mat_list


# If the sign of the first element in a column in matrices after the first is different than the first,
# multiply all values in that column by -1
def harmonizeColumnSigns_std(mat_list):   
    import numpy as np
    
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
    loadings = cp.stack(harm_loadings, axis = 2)[:, pc, :]
    # Convert to ranks if 'asRanks == True'
    if asRanks:
        loadings = cp.array([rankdata(loadings[:, i]) for i in range(loadings.shape[1])]).transpose()
    # Get sorted order based on median for each ASV 
    row_sort = cp.apply_along_axis(np.median, 1, loadings).ravel().argsort()[::-1]
    # Sort based on median, add ASV names (also sorted) and return data frame
    df = cudf.DataFrame(loadings[row_sort, :])
    df.index = asvs[row_sort]
    return df


# Sorts PCA loadings from a list 
def sortLoadings_std(loading_list, pc, asvs, asRanks = False):
    from scipy.stats import rankdata  
    import numpy as np
    import pandas as pd
    
    # Harmonize signs across replicates
    harm_loadings = harmonizeColumnSigns_std(loading_list)
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

# Returns a vector of logicals identifying values that are outside of the inter-quartile range
def isOutlierIQR(x):
    import numpy as np
    
    quarts = np.quantile(x, [0.25, 0.75])
    iqr = np.diff(quarts)
    lower = quarts[0] - 1.5 * iqr
    upper = quarts[1] + 1.5 * iqr
    thresh = [lower, upper]
    return np.concatenate([xi <= thresh[0] or xi >= thresh[1] for xi in x])

# Returns a vector of logicals identifying values that have Z-scores larger than the threshold
def isOutlierZ(x, thresh = 3):
    import scipy
    
    z = scipy.stats.zscore(x)
    return abs(z) >= thresh


def loadingSummary(loadings, pc, asv_labels, z_thresh = 3, pct_outlier_thresh = 0.95):
    import numpy as np
    import pandas
    
    metrics = {
        'loadings': sortLoadings_std(loadings, 0, asv_labels),
        'ranks': sortLoadings_std(loadings, 0, asv_labels, True)
    }
    metrics['outlier_IQR'] = metrics['loadings'].apply(isOutlierIQR)
    metrics['outlier_Z'] = metrics['loadings'].apply(isOutlierZ, thresh = z_thresh)
    
    smry = pandas.DataFrame()
    
    smry['mean_loading'] = metrics['loadings'].apply(np.mean, axis = 1)
    smry['median_loading'] = metrics['loadings'].apply(np.median, axis = 1)

    smry['mean_rank'] = metrics['ranks'].apply(np.mean, axis = 1)
    smry['median_rank'] = metrics['ranks'].apply(np.median, axis = 1)

    smry['pct_outlier_IQR'] = metrics['outlier_IQR'].apply(np.mean, axis = 1)
    smry['pct_outlier_Z'] = metrics['outlier_Z'].apply(np.mean, axis = 1)

    return {
        'metrics': metrics,
        'smry': smry, 
        'positive': smry[(smry['pct_outlier_Z'] > pct_outlier_thresh) & (smry['mean_loading'] > 0)],
        'negative': smry[(smry['pct_outlier_Z'] > pct_outlier_thresh) & (smry['mean_loading'] < 0)].sort_values(by = 'mean_loading')
    }


def plotScoreDistribution(scores, x = 0, y = 1):
    import pandas as pd
    import numpy as np
    import plotly.offline as pyo
    import plotly.graph_objs as go
    
    score_list = md.harmonizeColumnSigns_std(scores)
    score_arr = np.stack(score_list, axis = 2)

    median_score = pd.DataFrame([[np.median(score_arr[row, col, :]) for col in range(score_arr.shape[1])] for row in range(score_arr.shape[0])])
    min_score = pd.DataFrame([[np.min(score_arr[row, col, :]) for col in range(score_arr.shape[1])] for row in range(score_arr.shape[0])])
    max_score = pd.DataFrame([[np.max(score_arr[row, col, :]) for col in range(score_arr.shape[1])] for row in range(score_arr.shape[0])])

    medians = go.Scatter(
        x = median_score.iloc[:, 0],
        y = median_score.iloc[:, 1],
        mode = 'markers'
    )

    horiz_lines = [
        dict(
            type = 'line',
            x0 = min_score.iloc[i, 0],
            y0 = median_score.iloc[i, 1],
            x1 = max_score.iloc[i, 0],
            y1 = median_score.iloc[i, 1],
            line = dict(
                color = 'grey',
                width = 1
            )
        )
        for i in range(median_score.shape[0])
    ]

    vert_lines = [
        dict(
            type = 'line',
            x0 = median_score.iloc[i, 0],
            y0 = min_score.iloc[i, 1],
            x1 = median_score.iloc[i, 0],
            y1 = max_score.iloc[i, 1],
            line = dict(
                color = 'grey',
                width = 1
            )
        )
        for i in range(median_score.shape[0])
    ]

    go.Figure(
        medians, 
        go.Layout(shapes = horiz_lines + vert_lines, autosize = False, width = 1000, height = 1000)
    ).show()

# Does hierarchical clustering on data frame where rows are samples and columns are ASVs
# Returns array of cluster labels for rows
def doClustering(df, num_clusts, num_pcs = None):
    from sklearn.cluster import AgglomerativeClustering
    
    agg_clust = AgglomerativeClustering(n_clusters = num_clusts, metric = "euclidean", linkage = "ward")
    labels = agg_clust.fit_predict(df)
    return labels.astype(str)


# Function to create n_rep draws of df and assign n_clusters and returns percent of draws that had same relative 
# cluster assignments
def pctSame(df, n_clust, n_rep):
    import pandas as pd
    import numpy as np
    import itertools
    
    if n_clust >= df.shape[0]:
        return 100
    
    def isSameCluster(pws, df, col):
        return df.iloc[pws[0], col] == df.iloc[pws[1], col]
    
    def maxSame(row):
        return row.value_counts().max()
    
    # cluster a random sample of logit(relative percentages)
    cluster_samples = [doClustering(ranRelPct(df), n_clust) for i in range(n_rep)]
    cluster_samples = pd.DataFrame(cluster_samples).transpose()
    # unique pairs of rows
    pws_rows = itertools.combinations(range(cluster_samples.shape[0]), 2)
    # identify pairs of samples that are in the same cluster (True) or in different clusters (False)
    pws_same = pd.DataFrame([[isSameCluster(pair, cluster_samples, col) for col in range(cluster_samples.shape[1])] for pair in pws_rows])
    # get the maximum number replicates that have the same value (True or False) for each sample
    num_same = [maxSame(pws_same.iloc[row, :]) for row in range(pws_same.shape[0])]
    # convert to percentage with maximum of same value across all replicates
    return np.sum(num_same) * 100 / (pws_same.shape[0] * pws_same.shape[1])