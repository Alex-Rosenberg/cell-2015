# A3SS_Seqs.csv
Contains the sequences for the alternative 3' library

# A5SS_Seqs.csv
Contains the sequences for the alternative 5' library

# Reads.mat
Contains spliced reads at each position in the intron for each plasmid.
Can be loaded in python with the following:
import scipy.io as sio
data = sio.loadmat('../data/Reads.mat')
A5SS_data = data['A5SS']
A3SS_data = data['A3SS']
