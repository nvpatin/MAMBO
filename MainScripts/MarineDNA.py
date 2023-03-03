# Return data frame of a draw of relative percent of occurrence from a beta distribution
# fit to observed occurrence counts
#   df: data frame where rows = ASVs and columns = samples
def ranRelPct(df, asLogOdds = True):
    import pandas as pd
    import numpy as np
    
    def betaCol(col):
        beta_dist = np.random.beta(col + 1, col.sum() - col + 1)
        return beta_dist / beta_dist.sum()
    result = np.empty([df.shape[0], df.shape[1]])
    for i in range(result.shape[1]):
        result[:,i] = betaCol(df.iloc[:,i])
    if asLogOdds:
        result = np.log(result / (1 - result))
    return pd.DataFrame(result, index = df.index, columns = df.columns).transpose()

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