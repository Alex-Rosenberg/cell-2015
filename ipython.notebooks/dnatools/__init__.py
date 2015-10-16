import numpy as np
import scipy
import scipy.sparse
from pylab import *

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

def reverse_complement(seq):
    outseq = ''
    for s in seq:
        outseq = watsoncrick[s] + outseq
    return outseq

def hamdist(str1, str2):
   diffs = 0
   for ch1, ch2 in zip(str1, str2):
       if ch1 != ch2:
           diffs += 1
   return diffs
        
def get_snp_pos(ref,mut):
    for p in range(len(ref)):
        if(ref[p]!=mut[p]):
            break
    return p

def make_mer_matrix_no_pos(seqs,mer_len):
    mer_dict = dict(zip(make_mer_list(mer_len),range(4**mer_len)))
    rows,cols = [],[]
    r = 0
    for i in xrange(len(seqs)):
        cur_seq = seqs[i]
        for b in range(len(cur_seq)-mer_len+1):
            rows.append(r)
            cols.append(mer_dict[cur_seq[b:b+mer_len]])
        if(r%10000)==0:
            print r,
        r+=1
    vals = np.ones_like(cols)
    rows.append(r-1)
    cols.append(4**mer_len-1)
    vals = np.append(vals,0)
    X = scipy.sparse.csr_matrix((vals,(rows,cols)),dtype=np.float64)
    return X

def find_seq_diff_pos(wt_seq,mut_seq):
    """ Function to find the actual position of mutations based
    on the WT and MUT seq"""
    muts = np.zeros(len(wt_seq))
    for i in range(len(wt_seq)):
        muts[i] = (wt_seq[i]!=mut_seq[i])
    return find(muts)[0]

def get_snp_pos(ref,mut):
    for p in range(len(ref)):
        if(ref[p]!=mut[p]):
            break
    return p
