import pandas as pd
import numpy as np
from pylab import *
import os
import scipy.stats
import random
import time
import json
import scipy.io as sio
import scipy
import scipy.ndimage
import MLR_8regions

partition_num = int(sys.argv[1])
reg_lambda = float(sys.argv[2])
maxfun = int(sys.argv[3])
mer_len=6

data = sio.loadmat('../results/N7_A5SS_Model_Predictions/padded_6mer_matrix.mat')
A6pos  = data['X']
Y = scipy.matrix(data['Y'])
partition = np.genfromtxt('../results/N7_A5SS_Model_Predictions/Ten_fold_partition.txt')

training_inds = find((partition!=partition_num) & (partition!=((partition_num-1)%10)))
test_inds = find(partition==((partition_num-1)%10))

w,wfull,w0,f,d = MLR_8regions.MLR(A6pos[training_inds[:1000],:],Y[training_inds[:1000],:],reg_lambda,mer_len,maxfun=maxfun)
B = MLR_8regions.get_energy(A6pos,wfull,w0)[test_inds,:];

# Calculate kl-divergence
Ytest = Y[test_inds,:]
vv = scipy.sum(np.multiply(B,Ytest),axis=1)
Ypl = log(Ytest)
Ypl[Ytest==0]=0
Ypl = np.multiply(Ypl,Ytest)
vv = vv-scipy.sum(Ypl,axis=1)
V = -sum(vv)/float(len(Ytest))

f = open('../results/N7_A5SS_Model_Predictions/Reg8_lambdas.txt','a')
f.write(str(partition_num)+'\t'+str(reg_lambda)+'\t'+str(V)+'\n')
f.close()
