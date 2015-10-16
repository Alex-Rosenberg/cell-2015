import scipy
from pylab import *

def get_energy(X,w,w0):
    B = X*w+w0;
    mx = amax(B);
    Z = scipy.matrix(mx+log(scipy.sum(exp(B-mx),axis=1)));
    try:
        B = B-Z
    except:
        B = B-Z.transpose()
    return B
	
def value(w,*args):
	Xpos=args[0]
	Y=args[1]
	Gobs=args[2]
	Gobs0 = args[3]
	reg_lambda = args[4]
	mer_len = args[5]
	Wss = np.zeros(150*4**mer_len)
	Wmers = w[:-shape(Y)[1]].reshape(4**mer_len,8)
	Wss[:4**mer_len*39] = np.tile(Wmers[:,0],39)
	Wss[4**mer_len*39:4**mer_len*70] = np.tile(Wmers[:,1],31)
	for i in range(4):
		Wss[4**mer_len*(76+i):4**mer_len*(77+i)] = Wmers[:,i+2]
	Wss[4**mer_len*85:4**mer_len*114] = np.tile(Wmers[:,6],29)
	Wss[4**mer_len*114:] = np.tile(Wmers[:,7],36)
	w0 = w[-shape(Y)[1]:]
	Wfull = scipy.zeros((shape(Xpos)[1],shape(Y)[1]));
	for i in [0]+range(6,33)+[44]+range(49,75)+[79]:
		Wfull[(i)*4**mer_len:(i+150)*4**mer_len,i] = Wss
	
	B = get_energy(Xpos,Wfull,w0);
	vv = scipy.sum(np.multiply(B,Y),axis=1)
	Ypl = log(Y)
	Ypl[Y==0]=0
	Ypl = np.multiply(Ypl,Y)
	vv = vv-scipy.sum(Ypl,axis=1)
	V = sum(vv)/float(shape(Xpos)[0])# - sum(w[:-shape(Y)[1]]**2)*reg_lambda/2.
	rV = - sum(abs(w[:-shape(Y)[1]])*reg_lambda);
	V += rV
	Gpred,Gpred0 = grad_pred(Xpos,exp(B),mer_len)
	Gw = (Gobs - Gpred)
	G0 = (Gobs0 - Gpred0)
	G = np.concatenate((Gw.reshape(size(Gw),1),G0.reshape(size(G0),1)))
	G = G/float(shape(Xpos)[0])
	w = w.reshape(size(w),1)
	rG = - sign(w)*reg_lambda;
	rG[-shape(Y)[1]:] = 0;
	G += rG
	G = np.reshape(G,size(G))
	print -V
	return -V, -np.array(G)
	
def grad_obs(Xpos,Y,mer_len):
	Gobs_ss_big = (Y.transpose()*Xpos).transpose()
	Gobs = scipy.zeros((150*4**mer_len,1))
	for i in [0]+range(6,33)+[44]+range(49,75)+[79]:
		Gobs += Gobs_ss_big[(i)*4**mer_len:(i+150)*4**mer_len,i]
	Gobs = Gobs.reshape((150,4**mer_len)).transpose()
	Gw = np.zeros((4**mer_len,8))
	Gw[:,0] = sum(Gobs[:,:39],axis=1)
	Gw[:,1] = sum(Gobs[:,44:70],axis=1)
	Gw[:,2] = Gobs[:,76]
	Gw[:,3] = Gobs[:,77]
	Gw[:,4] = Gobs[:,78]
	Gw[:,5] = Gobs[:,79]
	Gw[:,6] = sum(Gobs[:,85:114],axis=1)
	Gw[:,7] = sum(Gobs[:,114:],axis=1)
	Gw = Gw.reshape(-1)
	Gobs0 = scipy.sum(Y.transpose(),axis=1)
	return np.array(Gw),np.array(Gobs0)
	
def grad_pred(Xpos,B,mer_len):
	Gpred_ss_big = (B.transpose()*Xpos).transpose()
	Gpred = scipy.zeros((150*4**mer_len,1))
	for i in [0]+range(6,33)+[44]+range(49,75)+[79]:
		Gpred += Gpred_ss_big[(i)*4**mer_len:(i+150)*4**mer_len,i]
	Gpred = Gpred.reshape((150,4**mer_len)).transpose()
	Gw = np.zeros((4**mer_len,8))
	Gw[:,0] = sum(Gpred[:,:39],axis=1)
	Gw[:,1] = sum(Gpred[:,44:70],axis=1)
	Gw[:,2] = Gpred[:,76]
	Gw[:,3] = Gpred[:,77]
	Gw[:,4] = Gpred[:,78]
	Gw[:,5] = Gpred[:,79]
	Gw[:,6] = sum(Gpred[:,85:114],axis=1)
	Gw[:,7] = sum(Gpred[:,114:],axis=1)
	Gw = Gw.reshape(-1)
	Gpred0 = scipy.sum(B.transpose(),axis=1)
	return np.array(Gw),np.array(Gpred0)
	
def MLR(Xpos,Y,reg_lambda,mer_len,w_init=None,maxfun=200):
	Gobs,Gobs0 = grad_obs(Xpos,Y,mer_len)
	w_init = scipy.zeros(4**mer_len*8+shape(Y)[1])
	w,f,d = scipy.optimize.fmin_l_bfgs_b(value,w_init,args=(Xpos,Y,np.array(Gobs),np.array(Gobs0),reg_lambda,mer_len),maxfun=maxfun)
	if False: #(not(d['warnflag']==0)):
		print 'Failed to converge...bummer'
		for i in d.keys():
			print i,d[i]
	
	Wss = np.zeros(150*4**mer_len)
	Wmers = w[:-shape(Y)[1]].reshape(4**mer_len,8)
	Wss[:4**mer_len*39] = np.tile(Wmers[:,0],39)
	Wss[4**mer_len*39:4**mer_len*70] = np.tile(Wmers[:,1],31)
	for i in range(4):
		Wss[4**mer_len*(76+i):4**mer_len*(77+i)] = Wmers[:,i+2]
	Wss[4**mer_len*85:4**mer_len*114] = np.tile(Wmers[:,6],29)
	Wss[4**mer_len*114:] = np.tile(Wmers[:,7],36)
	w0 = w[-shape(Y)[1]:]
	Wfull = scipy.zeros((shape(Xpos)[1],shape(Y)[1]));
	for i in [0]+range(6,33)+[44]+range(49,75)+[79]:
		Wfull[(i)*4**mer_len:(i+150)*4**mer_len,i] = Wss

	return w,Wfull,w0,f,d
