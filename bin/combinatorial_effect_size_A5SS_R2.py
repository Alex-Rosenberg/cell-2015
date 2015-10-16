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

# A5SS
A5SS_data = data['A5SS']
A5SS_data = np.array(A5SS_data.todense())
# Get minigenes with reads
A5SS_nn = find(A5SS_data.sum(axis=1))
A5SS_data = A5SS_data[A5SS_nn]
A5SS_data = A5SS_data/A5SS_data.sum(axis=1)[:,newaxis]
A5SS_seqs = pd.read_csv('../data/A5SS_Seqs.csv',index_col=0).Seq[A5SS_nn]


seqs = pd.Series(A3SS_seqs,index=range(len(A3SS_seqs))).str.slice(50,75)
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
if not os.path.exists('../results/N5_Combinatorial_Motif_Effects/A5SS/'):
	os.makedirs('../results/N5_Combinatorial_Motif_Effects/A5SS/')
df.to_pickle('../results/N5_Combinatorial_Motif_Effects/A5SS/comb_effects_mer_R2.df')






