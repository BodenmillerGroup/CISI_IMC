# Import libraries (from original script)
import numpy as np
from scipy.stats import entropy
import argparse


'''
Compute average gene abundances

inputs:
    X: proteins x cells (/pixels?) -> numpy array
    savepath: path where output will be stored
outputs:
    U: a dictionary of gene modules (cells/pixels x proteins)
    W: the module activity levels in each cell (/pixel?) of training data
       (modules x pixels)
'''

def gene_abundances(X, savepath):
	# X = np.load(args.datapath)
	# ### old, entropy-based way
	# AllGenes = np.load(args.all_genes)
	# f = open(args.selected_genes)
	# Genes = [line.strip() for line in f]
	# f.close()
	# gidx = [np.where(AllGenes == g)[0][0] for g in Genes]
	# X = X[gidx]
	# effective fraction (shannon diversity) for cells expressing gene
	# Xent = np.array([np.exp(entropy(x)) for x in X])/X.shape[1]

	# Average nonzero expression level of each gene
	Xavg = np.array([np.average(x[x>0]) for x in X])
	# # correct for shallow sequencing
	# RA = Xent/(1-np.exp(-Xavg/Xavg.max()*1.25)) # corresponds to a probability of ~29% of observing 0 counts for a gene expressed at Xavg.max() level
    # Same as Xavg?
    RA = np.average(X > 0, axis=1)

    # Save outputs
    path = Path(savepath)
    path.mkdir(parents=True, exist_ok=True)
    np.savetxt(os.path.join(savepath, 'relative_abundance.csv'), RA, delimiter=',')
    np.savetxt(os.path.join(savepath, 'average_on_expression_level.csv'), Xavg, delimiter=',')

    return RA, Xavg
