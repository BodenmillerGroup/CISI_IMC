import numpy as np
import argparse
from collections import defaultdict

if __name__ == '__main__':
# 	parser = argparse.ArgumentParser()
# 	parser.add_argument('--composites', help='Path to compositions (CSV channel, gene)')
# 	parser.add_argument('--genes', help='Path to gene labels')
# 	parser.add_argument('--savepath', help='Path to save')
    
	savepath = "dictionary/10_ms_2_max"
	A_path = savepath+"/version_48.csv"
    
	genes = "config/selected_genes.txt"
	
	f = open(genes)
	Genes = [line.strip() for line in f]
	f.close()
	f = open(A_path)
	header = f.readline()
	Phi_dict = defaultdict(list)
	for line in f:
		ls = line.strip().split(',')
		ls = ' '.join(ls).split()
		Phi_dict[int(ls[0])].extend(ls[1:])
	Phi = np.zeros((len(Phi_dict),len(Genes)))
	for k,v in Phi_dict.items():
		gi = [Genes.index(g) for g in v]
		Phi[k,gi] = 1
	np.save('%s/phi.npy' % savepath,Phi.astype(np.float32))