import pandas as pd
import numpy as np
from pylab import *
import os
import random
import time
import json
import copy
import sys
import scipy.io as sio

mer_len = int(sys.argv[1])
bases = ['A','T','C','G']
dna_dict = dict(zip(list('ATCG'),range(4)))
watsoncrick = {'N':'N','.':'.','C':'G','G':'C','A':'T','T':'A','*':'*'}
def add_base(li):
		"""Used in make_mer_list to add one more base to list"""
		new_li = []
		for s in li:
			for b in bases:
				new_li.append(s+b)
		return new_li

def make_mer_list(mer_len):
	"""Makes a list of all n-mers"""
	li = bases
	for i in range(mer_len-1):
		li = add_base(li)
	return li

data = sio.loadmat('../data/Reads.mat')
# A3SS
A3SS_data = data['A3SS']
# Only look at SA_1 usage:
A3SS_data = np.array(A3SS_data[:,235].todense()).reshape(-1)/np.array(A3SS_data.sum(axis=1),dtype=np.float64).reshape(-1)
# Get minigenes with reads
A3SS_nn = find(pd.notnull(A3SS_data))
A3SS_data = A3SS_data[A3SS_nn]
A3SS_seqs = pd.read_csv('../data/A3SS_Seqs.csv',index_col=0).Seq[A3SS_nn].values


seqs = pd.Series(A3SS_seqs,index=range(len(A3SS_seqs))).str.slice(-25)
data = A3SS_data

mer_list = make_mer_list(mer_len)

mer_counts = dict(zip(mer_list,np.zeros(len(mer_list))))
for seq in seqs:
    for b in range(25-mer_len+1):
        mer_counts[seq[b:b+mer_len]] += 1 
mer_counts = pd.Series(mer_counts)
comb_dict = {}
for mer1 in mer_list:
	print mer1,
	comb_dict[mer1] = {}
	mer1_found = seqs[seqs.str.contains(mer1)]
	for mer2 in mer_list:
		mer_found_inds = mer1_found[mer1_found.str.contains(mer1+'(.*)'+mer2)].index.values
		true_list = (seqs=='sdgsdg').values
		true_list[mer_found_inds] = True
		found = data[true_list].mean()
		not_found = data[true_list==False].mean()
		comb_dict[mer1][mer2] = log2(found/(1-found))-log2(not_found/(1-not_found))
df = pd.DataFrame(comb_dict)
if not os.path.exists('../results/N5_Combinatorial_Motif_Effects/A3SS/'):
	os.makedirs('../results/N5_Combinatorial_Motif_Effects/A3SS/')
df.to_pickle('../results/N5_Combinatorial_Motif_Effects/A3SS/comb_effects_mer_R2.df')






