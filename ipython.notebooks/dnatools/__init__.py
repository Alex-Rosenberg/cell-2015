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
