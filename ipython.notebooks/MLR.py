import scipy
from pylab import *

class MLR:

	def __init__(self,verbose=True):
		self.W = None
		self.b = None
		self.verbose = verbose

	def predict(self,X):
		if self.W is None:
			print 'Train model first.'
		else:
			B = X*self.W+self.b;
			mx = amax(B);
			Z = scipy.matrix(mx+log(scipy.sum(exp(B-mx),axis=1)));
			try:
				B = B-Z
			except:
				B = B-Z.transpose()
			return np.array(exp(B))

	def get_energy(self,X,W,b):
		B = X*W+b;
		mx = amax(B);
		Z = scipy.matrix(mx+log(scipy.sum(exp(B-mx),axis=1)));
		try:
			B = B-Z
		except:
			B = B-Z.transpose()
		return B

	def get_regularization_loss_grad(self,w_vector,X,Y,reg_type,reg_lambda):
		
		if(reg_type=='L1'):
			V = - sum(abs(w_vector[:shape(X)[1]*shape(Y)[1]]))*reg_lambda
			G = np.zeros_like(w_vector)
			G[:shape(X)[1]*shape(Y)[1]] -= reg_lambda*sign(w_vector[:shape(X)[1]*shape(Y)[1]])
		
		elif(reg_type=='L2'):
			V = - sum((w_vector[:shape(X)[1]*shape(Y)[1]])**2)*reg_lambda/2.
			G = np.zeros_like(w_vector)
			G[:shape(X)[1]*shape(Y)[1]] -= reg_lambda*w_vector[:shape(X)[1]*shape(Y)[1]]
		
		return V,G

	def get_loss_grad(self,w_vector,*args):
		X=args[0]
		Y=args[1]
		Gobs=args[2]
		reg_type = args[3]
		reg_lambda = args[4]
		wfull = scipy.reshape(w_vector,((shape(X)[1]+1),shape(Y)[1]))
		B = self.get_energy(X,wfull[:-1,:],wfull[-1,:])
		G_pred = scipy.hstack(((exp(B).transpose()*X),scipy.sum(exp(B).transpose(),axis=1)))
		vv = scipy.sum(np.multiply(B,Y),axis=1)
		Ypl = log(Y)
		Ypl[Y==0]=0
		Ypl = np.multiply(Ypl,Y)
		vv = vv-scipy.sum(Ypl,axis=1)
		V = sum(vv)/float(shape(X)[0])
		G = np.array(Gobs-G_pred).transpose()
		G = np.reshape(G,size(G))/float(shape(X)[0])
		
		V_reg, G_reg = self.get_regularization_loss_grad(w_vector,X,Y,reg_type,reg_lambda)
		V += V_reg
		G += G_reg
		if self.verbose:
			print -V,
		return -V, -np.array(G)

	def grad_obs(self,X,Y):
		Gobs = scipy.hstack(((Y.transpose()*X),scipy.sum(Y.transpose(),axis=1)))
		return Gobs

	def fit(self,X,Y,reg_type='L2',reg_lambda=0.0001,w_init=None,maxit=500):
		if self.W is None:
			self.W_vector = scipy.zeros((shape(X)[1]+1)*shape(Y)[1])
		else:
			W_bias = scipy.vstack((self.W,self.b[newaxis,:]))
			self.W_vector = W_bias.reshape((shape(X)[1]+1)*shape(Y)[1])
		Gobs = self.grad_obs(X,Y)
		self.W_vector,f,d = scipy.optimize.fmin_l_bfgs_b(self.get_loss_grad,self.W_vector,args=(X,Y,np.array(Gobs),reg_type,reg_lambda),maxfun=maxit)
		self.loss_val = f
		if (not(d['warnflag']==0)):
			print 'Failed to converge...'
			W_bias = scipy.reshape(self.W_vector,((shape(X)[1]+1),shape(Y)[1]))
		else:
			W_bias = scipy.reshape(self.W_vector,((shape(X)[1]+1),shape(Y)[1]))
	
		self.W, self.b = W_bias[:-1,:],W_bias[-1,:]



