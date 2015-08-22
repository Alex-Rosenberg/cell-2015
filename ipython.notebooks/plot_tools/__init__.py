from pylab import *

def simpleaxis(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

def plot_splicing_histogram(figname,title,data,
							y_max=26,
							y_step=5,
							smallplot=False,
							ylabel=False,
							mean_pos=None,
							smallplotpos=None,
							scale='T',
							figdir=None,
							save_plot=False):
	""" Plots histograms for splice site usage """
	if(scale=='M'):
		scalefactor=1000000.
	else:
		scalefactor=1000.

	fsize = 14
	fig = figure(figsize=(2.5,2))
	ax = fig.add_subplot(111)
	f,i = histogram(data,bins=arange(0,1.01,0.05));
	mean_splice = data.mean()
	ax.bar(i[:-1],(f)/scalefactor,width=0.05,color='gray')
	ax.axis([0,1.01,0,max(f)*1.1/scalefactor])
	ax.tick_params(labelsize=fsize)
	simpleaxis(ax)

	if(ylabel):
		if(scale=='M'):
			ax.set_ylabel('Plasmids (Millions)',fontsize=fsize)
		else:
			ax.set_ylabel('Plasmids (Thousands)',fontsize=fsize)
	ax.set_xlabel('Splicing Fraction at '+title,fontsize=fsize)
	ax.set_title(title ,fontsize=fsize)

	if(smallplot):
		rect = Rectangle((0, 0), 1.0, y_max, fc='None',ec='r',clip_on=False,linewidth=2)
		ax.add_patch(rect)
		if(smallplotpos):
			ax2 = axes([smallplotpos, 0.38, .49, .47])
		else:
			ax2 = axes([0.4, 0.38, .49, .47])
		ax2.bar(i[:-1],(f)/scalefactor,width=0.05,color='gray')
		#setp(ax2, xticks=[], yticks=[])
		ax2.axis([0,1,0,y_max]);
		ax2.set_yticks(arange(0,y_max,y_step))
		ax2.set_xticks(arange(0,1.01,0.5))
		ax2.tick_params(labelsize=fsize)
		plt.setp(ax2.spines.values(), color='r');
		[ax2.spines[i].set_linewidth(1.5) for i in ax2.spines]
		ax2.plot([mean_splice,mean_splice],[0,y_max],linewidth=2)
		if mean_pos:
			ax2.text(mean_pos,y_max*0.9,'%2.1f' %(mean_splice*100)+'%',fontsize=fsize,color='b',ha='left',va='top')
		else:
			ax2.text(mean_splice+0.1,y_max*0.9,'%2.1f' %(mean_splice*100)+'%',fontsize=fsize,color='b',ha='left',va='top')
	else:
		ax.plot([mean_splice,mean_splice],[0,max(f)/scalefactor],linewidth=2)
		ax.text(mean_splice*1.2,max(f)/scalefactor*0.9,'%2.1f' %(mean_splice*100)+'%',fontsize=fsize,color='b',ha='left',va='top')
	if save_plot and not (figdir is None):
		figname = 'Histogram_'+figname
		fig.savefig(figdir+figname+'.pdf',dpi=300,bbox_inches='tight')
		fig.savefig(figdir+figname+'.eps',dpi=300,bbox_inches='tight')
		fig.savefig(figdir+figname+'.png',dpi=300,bbox_inches='tight')
