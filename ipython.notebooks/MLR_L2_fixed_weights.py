import scipy
from pylab import *

def predict(X,w):
    B = X*w[:-1,:]+w[-1,:];
    mx = amax(B);
    Z = scipy.matrix(mx+log(scipy.sum(exp(B-mx),axis=1)));
    try:
        B = B-Z
    except:
        B = B-Z.transpose()
    return exp(B)

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
	X=args[0]
	Y=args[1]
	Gobs=args[2]
	reg_lambda = args[3]
	fixed_weights = args[4]
	wfull = scipy.reshape(w,((shape(X)[1]+1),shape(Y)[1]))
	B = get_energy(X,wfull[:-1,:],wfull[-1,:])
	G_pred = scipy.hstack(((exp(B).transpose()*X),scipy.sum(exp(B).transpose(),axis=1)))
	vv = scipy.sum(np.multiply(B,Y),axis=1)
	Ypl = log(Y)
	Ypl[Y==0]=0
	Ypl = np.multiply(Ypl,Y)
	vv = vv-scipy.sum(Ypl,axis=1)
	V = sum(vv)/float(shape(X)[0]) - sum(w[:shape(X)[1]*shape(Y)[1]]**2)*reg_lambda/2.
	G = np.array(Gobs-G_pred).transpose()
	G[:-1,:] -= reg_lambda*w.reshape(shape(G))[:-1,:]
	if not(fixed_weights is None):
		G[fixed_weights,:] = 0.
	G = np.reshape(G,size(G))/float(shape(X)[0])
	print -V,
	return -V, -np.array(G)
	
def grad_obs(X,Y):
    Gobs = scipy.hstack(((Y.transpose()*X),scipy.sum(Y.transpose(),axis=1)))
    return Gobs
	
def MLR(X,Y,reg_lambda,w_init=None,maxfun=200,fixed_weights=None):
	""" This function does multinomial logistic regression with L1 regularization
	X: nxd with n training cases and d features
	Y: nxk with n training cases and k classes
	fixed_weights: indices of weights that should not change when supplying w_init
	"""
	if w_init is None:
		w_init = scipy.zeros((shape(X)[1]+1)*shape(Y)[1])
	else:
		w_init = w_init.reshape((shape(X)[1]+1)*shape(Y)[1])
	Gobs = grad_obs(X,Y)
	w,f,d = scipy.optimize.fmin_l_bfgs_b(value,w_init,args=(X,Y,np.array(Gobs),reg_lambda,fixed_weights),maxfun=maxfun)
	if (not(d['warnflag']==0)):
		print 'Failed to converge...bummer'

	w = scipy.reshape(w,((shape(X)[1]+1),shape(Y)[1]))
	return w,f,d
