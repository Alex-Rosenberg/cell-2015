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
# A5SS
A5SS_data = data['A5SS']
A5SS_data = np.array(A5SS_data.todense())
# Get minigenes with reads
A5SS_data = A5SS_data[:,44]/(A5SS_data[:,0]+ A5SS_data[:,44])
A5SS_nn = pd.notnull(A5SS_data)
A5SS_data = A5SS_data[A5SS_nn]
A5SS_seqs = pd.read_csv('../data/A5SS_Seqs.csv',index_col=0).Seq[A5SS_nn].values

seqs = pd.Series(A5SS_seqs).str.slice(50,75)

mer6_list = make_mer_list(6)
c = 0
vals = arange(len(seqs),dtype=uint32)
num_bootstraps = 200
rand_inds = np.array([np.random.choice(vals,len(vals)) for j in range(num_bootstraps)])
vals = A5SS_data
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
if not os.path.exists('../results/N4_Motif_Effect_Sizes/A5SS/'):
	os.makedirs('../results/N4_Motif_Effect_Sizes/A5SS/')
np.save('../results/N4_Motif_Effect_Sizes/A5SS/boot_strapped_effects_R2', effects_array)
np.save('../results/N4_Motif_Effect_Sizes/A5SS/mean_effects_sdpos_R2', effects_mean)





