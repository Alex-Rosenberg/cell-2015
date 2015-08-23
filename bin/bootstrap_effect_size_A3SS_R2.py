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

seqs = pd.Series(A3SS_seqs).str.slice(-25)

mer6_list = make_mer_list(6)
c = 0
vals = arange(len(seqs),dtype=uint32)
num_bootstraps = 200
rand_inds = np.array([np.random.choice(vals,len(vals)) for j in range(num_bootstraps)])
vals = A3SS_data
effects_array = np.zeros((4**6,num_bootstraps))
effects_mean = np.zeros(4**6)
for mer in mer6_list:
	if(c%10)==0:
		print mer,c
	c += 1
	mer_found = seqs.str.contains(mer).values
	cur_mer = []
	for i in range(num_bootstraps):
		cur_mer_found = mer_found[rand_inds[i]]
		found_mean = mean(vals[rand_inds[i]][cur_mer_found])
		absent_mean = mean(vals[rand_inds[i]][cur_mer_found==False])		
		effects_array[c-1,i] = log2(found_mean/(1.-found_mean)) - log2(absent_mean/(1.-absent_mean))
	found_mean = mean(vals[mer_found])
	absent_mean = mean(vals[mer_found==False])
	effects_mean[c-1] = log2(found_mean/(1.-found_mean)) - log2(absent_mean/(1.-absent_mean))
	print mer
if not os.path.exists('../results/N4_Motif_Effect_Sizes/A3SS/'):
	os.makedirs('../results/N4_Motif_Effect_Sizes/A3SS/')
np.save('../results/N4_Motif_Effect_Sizes/A3SS/boot_strapped_effects_R2', effects_array)
np.save('../results/N4_Motif_Effect_Sizes/A3SS/mean_effects_sdpos_R2', effects_mean)




