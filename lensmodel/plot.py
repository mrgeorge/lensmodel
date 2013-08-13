#! env python

import matplotlib
matplotlib.use('Agg') # must appear before importing pyplot to get plots w/o GUI
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal
import scipy.stats
import fit

plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':20})
plt.rc('text', usetex=True)
plt.rc('axes',linewidth=1.5)


def contourPlot(xvals,yvals,smooth=0,percentiles=[0.68,0.95,0.99],colors=["red","green","blue"],xlabel=None,ylabel=None,xlim=None,ylim=None,filename=None,showPlot=False):
# make a 2d contour plot of parameter posteriors

    n2dbins=300

    # if it's a single ndarray wrapped in a list, convert to ndarray to use full color list
    if((type(xvals) is list) & (len(xvals) ==1)):
	xvals=xvals[0]
	yvals=yvals[0]
	
    if(type(xvals) is list):
	for ii in range(len(xvals)):
	    zz,xx,yy=np.histogram2d(xvals[ii],yvals[ii],bins=n2dbins)
	    xxbin=xx[1]-xx[0]
	    yybin=yy[1]-yy[0]
	    xx=xx[1:]+0.5*xxbin
	    yy=yy[1:]+0.5*yybin

	    if(smooth > 0):
		kernSize=int(10*smooth)
		sx,sy=scipy.mgrid[-kernSize:kernSize+1, -kernSize:kernSize+1]
		kern=np.exp(-(sx**2 + sy**2)/(2.*smooth**2))
		zz=scipy.signal.convolve2d(zz,kern/np.sum(kern),mode='same')
	
	    hist,bins=np.histogram(zz.flatten(),bins=1000)
	    sortzz=np.sort(zz.flatten())
	    cumhist=np.cumsum(sortzz)*1./np.sum(zz)
	    levels=np.array([sortzz[(cumhist>(1-pct)).nonzero()[0][0]] for pct in percentiles])

	    plt.contour(xx,yy,zz.T,levels=levels,colors=colors[ii])
    else: #we just have single ndarrays for xvals and yvals
        zz,xx,yy=np.histogram2d(xvals,yvals,bins=n2dbins)
	xxbin=xx[1]-xx[0]
	yybin=yy[1]-yy[0]
	xx=xx[1:]+0.5*xxbin
	yy=yy[1:]+0.5*yybin

	if(smooth > 0):
	    kernSize=int(10*smooth)
	    sx,sy=scipy.mgrid[-kernSize:kernSize+1, -kernSize:kernSize+1]
	    kern=np.exp(-(sx**2 + sy**2)/(2.*smooth**2))
	    zz=scipy.signal.convolve2d(zz,kern/np.sum(kern),mode='same')
	
	    hist,bins=np.histogram(zz.flatten(),bins=1000)
	    sortzz=np.sort(zz.flatten())
	    cumhist=np.cumsum(sortzz)*1./np.sum(zz)
	    levels=np.array([sortzz[(cumhist>(1-pct)).nonzero()[0][0]] for pct in percentiles])

	    plt.contour(xx,yy,zz.T,levels=levels,colors=colors)

    if(xlabel is not None):
	plt.xlabel(xlabel)
    if(ylabel is not None):
	plt.ylabel(ylabel)
    if(xlim is not None):
	plt.xlim(xlim)
    if(ylim is not None):
	plt.ylim(ylim)

    if(filename):
	plt.savefig(filename)
    if(showPlot):
	plt.show()
    
def contourPlotAll(chains,lnprobs=None,inputPars=None,showMax=True,showPeakKDE=True,show68=True,smooth=0,percentiles=[0.68,0.95,0.99],colors=["red","green","blue"],labels=None,figsize=(8,6),filename=None,showPlot=False):
# make a grid of contour plots for each pair of parameters
# chain is actually a list of 1 or more chains from emcee sampler

    nChains=len(chains)
    nPars=chains[0].shape[1]

    fig,axarr=plt.subplots(nPars,nPars,figsize=figsize)
    fig.subplots_adjust(hspace=0,wspace=0)

    if(labels is None):
	labels=np.repeat(None,nPars)

    # find max and min for all pars across chains
    limArr=np.tile((np.Inf,-np.Inf),nPars).reshape(nPars,2)
    for ch in chains:
	for par in range(nPars):
	    lo,hi=np.min(ch[:,par]), np.max(ch[:,par])
	    if(lo < limArr[par,0]):
		limArr[par,0]=lo.copy()
	    if(hi > limArr[par,1]):
		limArr[par,1]=hi.copy()

    # handle colors
    if(len(colors) == len(chains)):
	histColors=colors
	contourColors=colors
    if((nChains == 1) & (len(colors) == len(percentiles))):
	histColors=colors[0]
	contourColors=colors
	    
    # Get max posterior and width
    if((showMax) & (lnprobs is not None)):
        maxProbs=np.array([fit.getMaxProb(ch,lnp) for ch,lnp in zip(chains,lnprobs)])
    if((showPeakKDE) & (lnprobs is not None)):
        peakKDE=np.array([fit.getPeakKDE(ch,fit.getMaxProb(ch,lnp)) for ch,lnp in zip(chains,lnprobs)])
    if(show68):
        ranges=np.array([fit.get68(ch,opt="lowhigh") for ch in chains])
        
    # fill plot panels
    for row in range(nPars):
	for col in range(nPars):
	    fig.sca(axarr[row,col])

	    # setup axis labels
	    if(row == nPars-1):
		xlabel=labels[col]
		plt.setp(axarr[row,col].get_xticklabels(), rotation="vertical", fontsize="xx-small")
            else:
		xlabel=None
		plt.setp(axarr[row,col].get_xticklabels(),visible=False)
	    if(col == 0):
		ylabel=labels[row]
		plt.setp(axarr[row,col].get_yticklabels(), fontsize="xx-small")
            else:
		ylabel=None
		plt.setp(axarr[row,col].get_yticklabels(),visible=False)
		    
	    xarrs=[chain[:,col] for chain in chains]
	    yarrs=[chain[:,row] for chain in chains]
	    xlim=limArr[col]
	    ylim=limArr[row]
	    if(row == col):
                #		histvals=axarr[row,col].hist(xarrs,bins=50,range=xlim,histtype="step",color=histColors)
                xKDE=np.linspace(xlim[0],xlim[1],num=100)
                for ii in range(nChains):
                    kern=scipy.stats.gaussian_kde(xarrs[ii])
                    yKDE=kern(xKDE)
                    axarr[row,col].plot(xKDE,yKDE,color=histColors[ii])
                    if(showMax):
                        # add vertical lines marking the maximum posterior value
                        plt.plot(np.repeat(maxProbs[ii][col],2),np.array([0,kern(maxProbs[ii][col])]),color=histColors[ii],ls="-.")
                    if(showPeakKDE):
                        # add vertical lines marking the maximum posterior density value
                        plt.plot(np.repeat(peakKDE[ii][col],2),np.array([0,kern(peakKDE[ii][col])]),color=histColors[ii],ls=":")
                    if(show68):
                        # fill band marking 68% width
                        plt.fill_between(xKDE,yKDE,where=((xKDE > ranges[ii][0][col]) & (xKDE < ranges[ii][1][col])),color=histColors[ii],alpha=0.5)
                if(inputPars is not None):
                    # add vertical lines marking the input value
                    plt.plot(np.repeat(inputPars[col],2),np.array(plt.gca().get_ylim()),color="yellow",ls="--")

		if(xlabel is not None):
		    axarr[row,col].set_xlabel(xlabel)
		if(ylabel is not None):
		    axarr[row,col].set_ylabel(ylabel)
		axarr[row,col].set_xlim(xlim)
		plt.setp(axarr[row,col].get_yticklabels(),visible=False)
	    elif(col < row):
		contourPlot(xarrs,yarrs,smooth=smooth,percentiles=percentiles,colors=contourColors,xlabel=xlabel,ylabel=ylabel)
		axarr[row,col].set_xlim(xlim)
		axarr[row,col].set_ylim(ylim)
                if(inputPars is not None):
                    # add lines marking the input values
                    plt.plot(np.repeat(inputPars[col],2),ylim,color="yellow",ls="--")
                    plt.plot(xlim,np.repeat(inputPars[row],2),color="yellow",ls="--")
	    else:
		axarr[row,col].axis("off")

    fig.subplots_adjust(bottom=0.15)
    if(filename):
	fig.savefig(filename)
    if(showPlot):
	fig.show()

def plotProfile(xarr,yarr,xlim=None,ylim=None,xscale="log",yscale="log",xlabel="R (kpc)",ylabel=None,color="black",ls="-",lw=3,label=None):

    plt.plot(xarr,yarr,color=color,ls=ls,label=label,lw=lw)

    if(xlim is not None):
        plt.xlim(xlim)
    if(ylim is not None):
        plt.ylim(ylim)
    if(xscale is not None):
        if(xscale=="log"):
            plt.xscale(xscale,basex=10)
        else:
            plt.xscale(xscale)
    if(yscale is not None):
        if(yscale=="log"):
            plt.yscale(yscale,basey=10)
        else:
            plt.yscale(yscale)
    if(xlabel is not None):
        plt.xlabel(xlabel)
    if(ylabel is not None):
        plt.ylabel(ylabel)
