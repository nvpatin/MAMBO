import pandas as pd
import MarineDNA as md
import plotly.express as px
from sklearn.cluster import AgglomerativeClustering

file1 = "../../../Data/Flyer2018_16S_table_counts.tsv"
asvs1 = pd.read_csv(file1, index_col = 0, sep = "\t").transpose()
pca = [md.samplePCA(asvs1, 3) for i in range(10)]
test_df = pca[1]["df"]
test_scores = pd.DataFrame(pca[1]["scores"])

agg_clust = AgglomerativeClustering(n_clusters = 3, metric = "euclidean", linkage = "ward")
labels = agg_clust.fit_predict(test_df)
labels = labels.astype(str)

test_scores['labels'] = labels
fig = px.scatter(
    test_scores,
    x = 0,
    y = 1,
    color = "labels"
)
fig.show()

def samplePCA(x, n_clust):
    agg_clust = AgglomerativeClustering(n_clusters = n_clust, metric = "euclidean", linkage = "ward")
    labels = agg_clust.fit_predict(x)
    return labels.astype(str)