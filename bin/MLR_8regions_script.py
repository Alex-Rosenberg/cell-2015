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

reg_lambda = float(sys.argv[1])
mer_len = int(sys.argv[2])
maxfun = int(sys.argv[3])
output_path = sys.argv[4]


for i in range(5,len(sys.argv)):
	full_filepath = sys.argv[i]
	filename = '.'.join(full_filepath.split('.')[:-1]).split('/')[-1]
	if(i==5):
		Y = sio.loadmat(full_filepath)[filename].todense()
	else:
		Y += sio.loadmat(full_filepath)[filename].todense()

GT_inds = sio.loadmat('/net/shendure/vol7/abros/nobackup/Splicing_Project/Library_Files/GT_GC_Positions.mat')
Y = np.multiply(Y,GT_inds['inds'].todense())
A6pos  = sio.loadmat('../results/N7_A5SS_Model_Predictions/padded_6mer_matrix.mat')['X']
training = np.loadtxt('/net/shendure/vol7/abros/nobackup/Splicing_Project/Library_Files/training.csv')==1
print shape(training)
training_inds = find(training==1)
Y = scipy.hstack((Y[:,:80],Y[:,-1:]))
reads = float64(sum(Y,axis=1))
Y = Y/reads

not_null = np.array(reads>0).reshape(len(training))
print shape(not_null)
temp_training = (training==True)
print shape(temp_training)
temp_training[not_null==False] = False
temp_training_inds = find(temp_training)[:10000]
print 'Starting MLR...'

w,wfull,w0,f,d = MLR_8regions.MLR(A6pos[temp_training_inds,:],Y[temp_training_inds,:],reg_lambda,mer_len,maxfun=maxfun)

Ypred = exp(MLR_8regions.get_energy(A6pos,wfull,w0))

MLR_data = {'mer_scores':Wss,'w0':w0,'w':w,'Prediction':Ypred,'Training':training,'Data':Y,'Reads':reads,'NotNull':not_null}
if not os.path.exists(output_path):
    os.makedirs(output_path)
sio.savemat(output_path + 'Training_data.mat',MLR_data)
