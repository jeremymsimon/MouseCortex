import pandas as pd
from sklearn.manifold import TSNE
import matplotlib as mpl
import matplotlib.pyplot as plt

#Import data table from PCA dimensionality reduction, columns are top X significant PCs, rows are cell IDs
e14_PCA_reduced = pd.read_table('DropSeqAnalysis_e14_WTonly_batchCorrected_WTonly_PCAdimReduced_top89PCs.txt',index_col=0)

#Run t-SNE
tsne_structure = TSNE(n_components=2, random_state=123, perplexity=50, n_iter=5000,learning_rate=750)
e14_tsne = pd.DataFrame(tsne_structure.fit_transform(e14_PCA_reduced), index=e14_PCA_reduced.index)

#Write out results table
e14_tsne.to_csv('DropSeqAnalysis_e14_WTonly_batchCorrected_WTonly_PCAdimReduced_top89PCs_perp50_learn750_tSNE.txt',sep='\t')

