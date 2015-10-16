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

training = (partition!=partition_num)
if(partition_num==-1):
	test = (partition>-1)
else:
	test = (partition==partition_num)

training_inds = find(training)
test_inds = find(test)

w,wfull,w0,f,d = MLR_8regions.MLR(A6pos[training_inds,:],Y[training_inds,:],reg_lambda,mer_len,maxfun=maxfun)
Ypred = exp(MLR_8regions.get_energy(A6pos,wfull,w0))

MLR_data = {'Mer_scores':w,'w0':w0,'Prediction':Ypred,'Training':training,'Test':test,'Data':Y,'Reads':reads,'NotNull':not_null}
if not os.path.exists('../results/N7_A5SS_Model_Predictions/Partition'+str(partition_num)):
    os.makedirs('../results/N7_A5SS_Model_Predictions/Partition'+str(partition_num))
sio.savemat('../results/N7_A5SS_Model_Predictions/Partition'+str(partition_num)+ '/Training_data.mat',MLR_data)
