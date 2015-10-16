import pandas as pd
import numpy as np
import scipy
import scipy.sparse
import scipy.stats
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

def make_mer_effect_fig(resultsdir,event_type,region,figname=None,fsize=14,figdir=None,savefig=False):
    bootstrap = np.load(resultsdir+event_type+'/boot_strapped_effects_'+region+'.npy')
    mean_vals = np.load(resultsdir+event_type+'/mean_effects_sdpos_'+region+'.npy')
    mer6_list = make_mer_list(6)
    bins = arange(-5,5,0.001)
    mer6_hists = np.array([(histogram(bootstrap[i,:],bins=bins)[0]) for i in range(len(bootstrap))])
    mer_cdf = cumsum(mer6_hists,axis=1)
    num_mers = 5
    fig = figure(figsize=(8,3.5))
    ax = fig.add_axes((0,0,0.4,1))
    mer_cdf = cumsum(mer6_hists,axis=1)
    maxes = bins[argmax(mer6_hists,axis=1)]
    mer6_arr = np.array(mer6_list)
    inds = argsort(maxes)
    c_low = bins[argmax((mer_cdf<5)[:,:-1]-(mer_cdf<5)[:,1:],axis=1)]
    c_high = bins[argmax((mer_cdf>=195)[:,:-1]-(mer_cdf>=195)[:,1:],axis=1)]
    xmin = -5#int(floor(min(c_low)))
    xmax = 4#int(ceil(max(c_high)))
    #ax.set_title('6-mers in Mut Region 1',fontsize=fsize)
    ax.plot([-100],[-100],marker='o',markersize=5,color='b',linewidth=0,label='Effect Size')
    ax.plot([-100,-101],[-100,-101],color='g',linewidth=1,label='95% CI')
    ax.errorbar(maxes[inds],arange(mer6_hists.shape[0]),xerr=[maxes[inds]-c_low[inds],c_high[inds]-maxes[inds]],\
             color='g',capsize=0,linewidth=1,fmt='o',markersize=2,markerfacecolor='b',markeredgecolor='None',alpha=0.05)
    ax.scatter(maxes[inds],arange(mer6_hists.shape[0]),marker='o',s=2,color='b',edgecolor='None')
    #ax.scatter(maxes[inds[:4**6]],arange(4**6),s=10)
    ax.xaxis.set_ticks(range(xmin,xmax+1));
    plot([xmin,xmax],[-30,-30],'k',linewidth=1)
    plot([0,0],[-30,4**6+30],'k',linewidth=1)
    ax.set_xlim(xmin,xmax)
    ax.yaxis.tick_right()
    ax.yaxis.set_ticks(arange(0,4**6,20));
    box('off')
    ax.tick_params(size=0,labelsize=fsize)
    leg = legend(bbox_to_anchor=(0.6,0.65),prop={'size':fsize},numpoints=1,fancybox=True)
    ax.legend_.get_frame().set_alpha(0.0)
    ax.yaxis.set_ticklabels(np.array(mer6_list)[inds[::20]],fontsize=4);
    ax.set_ylim(-50,4**6+50)
    ax.set_ylabel('6nt Sequences (4096)',fontsize=fsize)
    #ax.set_xlabel('Effect Size ($\Delta$ $log_2$odds ratio)',fontsize=fsize)
    ax.yaxis.set_label_position("right")
    ax.annotate(
        '', xy=(0, 100), xycoords = 'data',
        xytext = (xmax, 100), textcoords = 'data',
        arrowprops = {'arrowstyle':'<-'})
    ax.annotate(
        'Enhancers', xy=(0, 10), xycoords = 'data',
        xytext = (xmax, 10), textcoords = 'offset points',fontsize=fsize)

    ax.annotate(
        '', xy=(xmin, 3500), xycoords = 'data',
        xytext = (0, 3500), textcoords = 'data',
        arrowprops = {'arrowstyle':'->'})
    ax.annotate(
        'Silencers', xy=(0, 3500), xycoords = 'data',
        xytext = (xmin, 5), textcoords = 'offset points',fontsize=fsize,ha='right')
    ax.text(5,-800,'Effect Size ($\Delta$ $log_2$odds ratio)',fontsize=fsize,ha='center')
    
    
    ax = fig.add_axes((0.5,0.0,0.4,0.45))
    mer_cdf = cumsum(mer6_hists,axis=1)
    maxes = bins[argmax(mer6_hists,axis=1)]
    mer6_arr = np.array(mer6_list)
    inds = argsort(maxes)[:num_mers]
    c_low = bins[argmax((mer_cdf<5)[:,:-1]-(mer_cdf<5)[:,1:],axis=1)]
    c_high = bins[argmax((mer_cdf>=195)[:,:-1]-(mer_cdf>=195)[:,1:],axis=1)]
    ax.plot([-100],[-100],marker='o',markersize=5,color='b',linewidth=0,label='Effect Size')
    ax.plot([-100,-101],[-100,-101],color='g',linewidth=1,label='95% CI')
    ax.errorbar(maxes[inds],arange(num_mers),xerr=[maxes[inds]-c_low[inds],c_high[inds]-maxes[inds]],\
             color='g',capsize=0,linewidth=2,fmt='o',markersize=5,markerfacecolor='b',markeredgecolor='None')
    #ax.scatter(maxes[inds[:4**6]],arange(4**6),s=10)
    ax.xaxis.set_ticks(range(xmin,xmax+1));
    ax.plot([xmin,xmax],[-1,-1],'k',linewidth=1)
    ax.plot([0,0],[-1,num_mers-0.5],'k',linewidth=1)
    ax.set_xlim(xmin,xmax)
    ax.yaxis.tick_right()
    ax.yaxis.set_ticks(arange(0,num_mers));
    ax.tick_params(size=0,labelsize=fsize)
    ax.yaxis.set_ticklabels(np.array(mer6_list)[inds],fontsize=fsize);
    box('off')
    ax.set_ylim(-2,num_mers)
    #ax.set_xlabel('Effect Size ($\Delta$ $log_2$odds ratio)',fontsize=fsize)
    ax.yaxis.set_label_position("right")


    ax = fig.add_axes((0.5,0.5,0.4,0.45))
    mer_cdf = cumsum(mer6_hists,axis=1)
    maxes = bins[argmax(mer6_hists,axis=1)]
    mer6_arr = np.array(mer6_list)
    inds = argsort(maxes)[-num_mers:]
    c_low = bins[argmax((mer_cdf<5)[:,:-1]-(mer_cdf<5)[:,1:],axis=1)]
    c_high = bins[argmax((mer_cdf>=195)[:,:-1]-(mer_cdf>=195)[:,1:],axis=1)]
    ax.plot([-100],[-100],marker='o',markersize=5,color='b',linewidth=0,label='Effect Size')
    ax.plot([-100,-101],[-100,-101],color='g',linewidth=1,label='95% CI')
    ax.errorbar(maxes[inds],arange(num_mers),xerr=[maxes[inds]-c_low[inds],c_high[inds]-maxes[inds]],\
             color='g',capsize=0,linewidth=2,fmt='o',markersize=5,markerfacecolor='b',markeredgecolor='None')
    #ax.scatter(maxes[inds[:4**6]],arange(4**6),s=10)
    ax.xaxis.set_ticklabels([])

    ax.plot([0,0],[-0.5,num_mers-0.5],'k',linewidth=1)
    ax.set_xlim(xmin,xmax)
    ax.yaxis.tick_right()
    ax.yaxis.set_ticks(arange(0,num_mers));
    ax.tick_params(size=0,labelsize=fsize)
    ax.yaxis.set_ticklabels(np.array(mer6_list)[inds],fontsize=fsize);
    box('off')
    ax.set_ylim(-2,num_mers)
    ax.yaxis.set_label_position("right")
    ax.text(xmin+(xmax-xmin)*1.18,-1.75,'...',fontsize=fsize,rotation=90,ha='center',va='center')
    ax.text(-0.07,-1.75,'...',fontsize=fsize,rotation=90,ha='center',va='center')
    if figname and savefig:
        fig.savefig(figdir+figname+'.png', bbox_inches='tight',dpi=200)
        fig.savefig(figdir+figname+'.eps', bbox_inches='tight',dpi=200)
        fig.savefig(figdir+figname+'.pdf', bbox_inches='tight',dpi=200)
