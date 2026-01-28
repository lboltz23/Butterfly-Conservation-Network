import pandas as pd
from scipy.spatial.distance import pdist, squareform
from skbio.stats.distance import DistanceMatrix, permanova

traits = pd.read_csv("trait_incidence_with_communities.csv")

exclude = {
    "scientific_name",
    "leiden_community",
    "louvain_cluster",
    "gb_threatened_binary",
    "ie_threatened_binary"
}

trait_cols = [c for c in traits.columns if c not in exclude]
X = traits[trait_cols].values

dist = squareform(pdist(X, metric="jaccard"))
dm = DistanceMatrix(dist)

result = permanova(
    dm,
    traits["leiden_community"],
    permutations=999
)

print(result)