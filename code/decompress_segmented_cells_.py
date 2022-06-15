# Import libraries (from original script)
import numpy as np
import argparse, os
from optimize_dictionary import sparse_decode_blocks
import pandas as pd
import scipy.spatial as sp, scipy.cluster.hierarchy as hc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.style.use('ggplot')
import seaborn as sns
from scipy.spatial import distance
from scipy.stats import entropy
from scipy.sparse import load_npz
import imageio

# Import libraries (additionaly added)
from utils import check_file
from pathlib import Path
import errno


# Find genes that are comeasured but are not coexpressed (correlation<threshold)
# and correct their expressions
def select_and_correct_comeasured(x, xc, phi, phi_corr, training_corr,
                                  phi_thresh=0.6, train_thresh=0.1):
	# find comeasured genes that are not coexpressed
	comeasured = []
	for i in range(phi_corr.shape[0]):
		xs = np.argsort(-phi_corr[i])
		for j in xs:
			if phi_corr[i, j] < phi_thresh:
				break
			comeasured.append((phi_corr[i, j], i, j))
	corrected_pairs = []
	for c in sorted(comeasured, reverse=True):
		if training_corr[c[1], c[2]] < train_thresh:
			x, both_nz = correct_coexpression(x, xc, phi, c[1], c[2])
			corrected_pairs.append((c[1], c[2], both_nz))
	return x, corrected_pairs


# Set closest not coexpressed genes to 0?
def correct_coexpression(x, xc, phi, i, j):
	# pick the gene with nearest expression pattern in scRNA
	thresh_i = np.percentile(x[i], 99.9) / 100
	thresh_j = np.percentile(x[j], 99.9) / 100
	both_nz = (x[i] > thresh_i)*(x[j] > thresh_j)
	dist = distance.cdist([phi[:, i], phi[:, j]], xc[:, both_nz].T, 'correlation')
	i_closer = np.where(both_nz)[0][dist[0] < dist[1]]
	j_closer = np.where(both_nz)[0][dist[0] > dist[1]]
	x[i, j_closer] = 0
	x[j, i_closer] = 0
	return x, both_nz


# If genes are in specified percentile, then there are set to 0
def threshold_genes(X, p=99.9, t0=0.05):
	ptile = np.percentile(X, p, axis=1)*t0
	X.T[X.T < ptile] = 0
	return X



def decompress_cells(basepath, npypath=None, trainpath, outpath, Y_name,
					 X_name, modules='gene_modules.npy', phi_name='version_2.txt',
					 use_test = False, correct_comeasured = False,
					 save_output = True, method = "average_intensity",
					 region_mask = None, rois=[1], proteins=[]):
# 	parser = argparse.ArgumentParser()
# 	parser.add_argument('--basepath', help='Path to directory with subdirs for parsed images in each tissue',default="output")
# 	parser.add_argument('--npypath', help='Path to directory for compressed and individual cell expression',default="preparation/npy_decoding")
# 	parser.add_argument('--tissues', help='Comma-separated list of tissue numbers to include',default="2")
# 	parser.add_argument('--trainpath', help='Path to directory with composition (phi) and gene module (U) matrices',default='preparation/training')
# 	parser.add_argument('--modules', help='Filename of gene module dictionary', default='gene_modules')
# 	parser.add_argument('--method', help='Signal integration method',default='average_intensity')
# 	parser.add_argument('--region-mask', help='Mask for cells in a region of interest',default=None)
# 	parser.add_argument('--correct-comeasured',dest='correct_comeasured', action='store_true')
# 	parser.add_argument('--use-test-idx',dest='use_test', action='store_true')
# 	parser.add_argument('--save-output',dest='save_output', action='store_true')
# 	parser.set_defaults(correct_comeasured=False)
# 	parser.set_defaults(use_test=False)
# 	parser.set_defaults(save_output=True)
# 	args,_ = parser.parse_known_args()

	os.environ['KMP_WARNINGS'] = '0'
	# f = open('%s/%s/selected_genes.txt' % (basepath,configpath))
	# AllGenes = np.array([line.strip() for line in f])
	# f.close()

	# Get protein names
	AllGenes = proteins

	# Check that phi/A exists,if header exists and read into numpy array
	if utils.check_file(os.path.join(trainpath, phi_name), ['.txt']):
		has_header = open(Path(os.path.join(trainpath, phi_name)), 'r').read().find('channels')
		if has_header!=-1:
			phi = np.loadtxt(os.path.join(trainpath, phi_name),
					   skiprows=1, usecols=list(range(2, len(AllGenes)+2)))
		else:
			phi = np.loadtxt(os.path.join(trainpath, phi_name),
					   usecols=list(range(2, len(AllGenes)+2)))

    else:
        # If file is not found, throw error
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
                                os.path.join(trainpath, phi_name))



	# phi = np.load('%s/%s/phi.npy' % (basepath,dictpath))
	# Calculate correlation of phi?
	phi_corr = (np.einsum('ij,ik->jk', phi, phi)/phi.sum(0)).T - np.eye(phi.shape[1])

    if utils.check_file(os.path.join(trainpath, 'relative_abundance.npy'), ['.npy']):
		abundance_train = np.load(os.path.join(trainpath, 'relative_abundance.npy'))
    if utils.check_file(os.path.join(trainpath, 'average_on_expression_level.npy'), ['.npy']):
		abundance_train = np.load(os.path.join(trainpath, 'average_on_expression_level.npy'))

	# train_corr = np.load('%s/correlations.npy' % args.trainpath)
	# train_corr = distance.squareform(train_corr)

	# Load dictionary U
	if utils.check_file(os.path.join(trainpath, modules), ['.npy']):
		U = np.load(os.path.join(trainpath, modules))

	# U = np.load('%s/%s/gene_modules.npy' % (basepath,dictpath))
	# Inliers = np.load('%s/%s.inliers.npy' % (args.trainpath, args.modules),allow_pickle=True)

	if use_test:
		if utils.check_file(os.path.join(basepath, 'test_idx.npy'), ['.npy']):
			test_idx = np.load(os.path.join(basepath, 'test_idx.npy'))
		# test_idx = np.load('%s/%s.test_idx.npy' % (basepath, dictpath))
	else:
		test_idx = np.arange(1e7, dtype=np.int)

	idx_offset = 0
# 	T = [int(t) for t in args.tissues.split(',')]
	X0 = []
	X1 = []
	X2 = []
	C = []

    # Decompress and analize each tissue
    for r in rois:
        # If npypath, then X_name and Y_name should be numpy arrays
        if npypath is None:
            composite_measurements = Y_name[r]
            composite_measurements = np.float32(composite_measurements)
            direct_measurements = X_name[r]
            direct_measurements = np.float32(direct_measurements)
        # If npypath is set, than X_name and Y_name should be the names of the files
        # in npypath that should be read in
        else:
    		if utils.check_file(os.path.join(npypath, (Y_name + r + '.csv')), ['.csv']):
    			composite_measurements = np.genfromtxt(os.path.join(npypath,
                                                           (Y_name + r + '.csv')),
                                              delimiter=',', skip_header=1))
                composite_measurements = np.float32(composite_measurements)
            if utils.check_file(os.path.join(npypath, (X_name + r + '.csv')), ['.csv']):
    			direct_measurements = np.genfromtxt(os.path.join(npypath,
                                                           (X_name + r + '.csv')),
                                              delimiter=',', skip_header=1))
                direct_measurements = np.float32(direct_measurements)
    	direct_labels = AllGenes
        # TODO: Change if we need this part
		if region_mask is not None:
			CellMasks = load_npz('%s/output/%s/segmented/cell_masks.size_threshold.npz' % (basepath, npypath))
			region_mask = imageio.imread('%s/output/stitched_aligned_filtered/%s' % (basepath, npypath))
			cidx = CellMasks.dot(region_mask.flatten())
			cidx = np.where(cidx)[0]
			composite_measurements = composite_measurements[:, cidx]
			direct_measurements = direct_measurements[:, cidx]
        # Decode W blockwise?
		W = sparse_decode_blocks(composite_measurements, phi.dot(U).astype(np.float32),
                           0.1)
		print('Average W entropy: %.2f' % np.average([np.exp(entropy(abs(w)))
                                                for w in W.T if w.sum()]))
        # Compute predicted X_hat
        xhat = U.dot(W)
        # Set nan and values smaller than 0 to 0 in X_hat (not possible values)
		xhat[np.isnan(xhat)] = 0
		xhat[xhat < 0] = 0
        # If method == 'integrated_intensity_binary', then set genes in 99th
        # percentile to 0
		if method == 'integrated_intensity_binary':
			xhat = threshold_genes(xhat)
        # If save_output is set, then save W and X_hat
		if save_output:
			path = Path(outpath)
            path.mkdir(parents=True, exist_ok=True)
			np.savetxt(os.path.join(outpath, ('decoded_module_activity_%s.csv' % r),
                           W, delimiter=','))
            np.savetxt(os.path.join(outpath, ('decoded_expression%s.csv' % r),
                           xhat, delimiter=','))
        # Correct with phi_corr and training_corr
		if correct_comeasured:
 			xhat, cp = select_and_correct_comeasured(xhat, composite_measurements,
                                             phi, phi_corr, train_corr)
            path = Path(outpath)
            path.mkdir(parents=True, exist_ok=True)
			np.savetxt(os.path.join(outpath, ('segmented.segmentation.adjusted_%s.csv' % r),
                           xhat, delimiter=','))
		# # exclude outliers for evaluation
		# composite_measurements = composite_measurements[:,Inliers[ti]]
		# direct_measurements = direct_measurements[:,Inliers[ti]]
		# xhat = xhat[:,Inliers[ti]]

		# Only evaluate on test data
		test_idx_t = test_idx[(test_idx >= idx_offset)*
                        (test_idx < idx_offset+composite_measurements.shape[1])] - idx_offset
		direct_measurements = direct_measurements[:, test_idx_t]
		xhat = xhat[:, test_idx_t]
		idx_offset += composite_measurements.shape[1]
		X1.append(xhat)
# 		idx = [np.where(AllGenes == l)[0][0] for l in direct_labels]
# 		xhat = xhat[idx]
		n = direct_measurements.shape[0]+1
        # Create and save figures comparing results
		if save_output:
			fig, axes = plt.subplots(max(2, int(np.floor(np.sqrt(n)))),
                            int(np.ceil(np.sqrt(n))), figsize=(15, 15))
			plt.subplots_adjust(wspace=0.5, hspace=0.5)
			axes = axes.flatten()
			plt.rcParams["axes.labelsize"] = 2.5
		for i in range(n-1):
			corr = 1-distance.correlation(direct_measurements[i], xhat[i])
			if save_output:
				sns_plt = sns.scatterplot(direct_measurements[i], xhat[i], ax=axes[i])
				_= sns_plt.set(xlabel='Direct Intensity', ylabel='Recovered Intensity',
                   title='%s; r=%.4f' % (direct_labels[i], corr))
			print(r, direct_labels[i], corr)
			C.append((r, direct_labels[i], corr))
		if save_output:
			corr = 1-distance.correlation(direct_measurements.flatten(), xhat.flatten())
			sns_plt = sns.scatterplot(direct_measurements.flatten(), xhat.flatten(),
                             ax=axes[-1])
			_= sns_plt.set(xlabel='Direct Intensity', ylabel='Recovered Intensity',
                  title='%s; r=%.4f' % ('all points', corr))
			print(r, 'all points', corr)
# 			_=plt.tight_layout()
			if correct_comeasured:
				path = os.path.join(outpath, ('scatter.segmented.heuristic_correction_%s.png' % r))
			else:
				path = os.path.join(outpath, ('scatter.segmented%s_pergene_%s.png' % r))
			fig.savefig(path)
			plt.close()
		X0.append(direct_measurements.flatten())
		X2.append(xhat.flatten())

    # Analyse all tissues
	X0 = np.hstack(X0)
	X1 = np.hstack(X1)
	X2 = np.hstack(X2)
	corr = 1-distance.correlation(X0, X2)
	corr_gene_avg = np.average([c[2] for c in C])
	if save_output:
		sns_plt = sns.scatterplot(X0, X2)
		_= sns_plt.set(xlabel='Direct Intensity (arb. units)',
                 ylabel='Recovered Intensity (arb. units)',
                 title='%s; r=%.4f' % ('all tissues, all points', corr))
	print('all tissues, all points',corr)
	print('average gene corr', np.average(corr_gene_avg))
	avg_ng = np.average([np.exp(entropy(x)) for x in X1.T if x.sum()>0])
	print('average genes / cell:', avg_ng)

    # Save overall comparison outputs
	if save_output:
		fig = sns_plt.get_figure()
		if correct_comeasured:
			path = os.path.join(outpath, 'scatter.segmented.heuristic_correction.pdf')
		else:
            path = os.path.join(outpath, 'scatter.segmented.pdf')
		fig.savefig(path)
		plt.close()

		f = open(os.path.join(outpath, 'summary.segmentation.txt', 'w'))
		_=f.write('\t'.join(['Overall correlation', 'Average gene correlation',
                       'Genes per cell']) + '\n')
		_=f.write('\t'.join([str(x) for x in [corr, corr_gene_avg, avg_ng]]) + '\n')
		_=f.write('\n')
		_=f.write('\t'.join(['Gene', 'Tissue', 'Correlation across segmented cells',
                       'abundance_train', 'abundance_decode']) + '\n')
		for c in C:
 			el = abundance_train[np.where(AllGenes == c[1])[0][0]]
 			ra = abundance_decode[np.where(AllGenes == c[1])[0][0]]
 			_=f.write('\t'.join([c[1], c[0], str(c[2]), str(el), str(ra)]) + '\n')
		f.close()
