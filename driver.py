#! env python

import lensmodel
import numpy as np
import matplotlib.pyplot as plt


plotDir="../plots/"
# Load data to get binning and error bars
dataDir="../data/"
inFile=dataDir+"rebin.matt.10x.zebracut.src.sm_11_11.5.big.out"
x1,y1,e1=lensmodel.io.readData(inFile)
redshift=0.1
xshear,yshear,errshear=lensmodel.io.convertCosmo(x1,y1,e1,redshift=redshift)

# Use only bins at R<1Mpc
fitSel=(xshear < 1000.)


# Input model pars
logMstars=11.2
logRstars=1.
logMhalo=13.5
conc=5.
innerSlopeGNFW=1.
nuDutton=0.
AGnedin=-1.
wGnedin=-1.
inputPars=[logMstars, logRstars, logMhalo, conc, innerSlopeGNFW, nuDutton, AGnedin, wGnedin]
cenType="hernquist"
odType="critical"
delta=200.

#priors=[logMstars,logRstars,(logMhalo,1.),[2.,10.],innerSlopeGNFW,[-0.3,1.2],-1.,-1.]
priors=[logMstars,logRstars,[12.,14.],[2.,10.],[0.,1.5],[-0.5,1.2],-1.,-1.]
freePars=inputPars[2:6]
labels=np.array(["log(Mh)","conc","inner slope","nuDutton"])

yshear=lensmodel.profiles.deltaSigma(inputPars,xshear)
xmag=xshear
ymag=lensmodel.profiles.sigma(inputPars,xmag)
errmag=1.*errshear*(ymag/yshear) # say S/N (mag) = S/N (shear) at all radii

# Run chains
nWalkers=2000
nBurn=50
nSteps=250
nThreads=8
seed=None

chains,lnprobs=lensmodel.fit.fitObs(priors,xshear,yshear,errshear,xmag,ymag,errmag,redshift=redshift,cenType=cenType,delta=delta,odType=odType,nWalkers=nWalkers,nBurn=nBurn,nSteps=nSteps,nThreads=nThreads,seed=seed)

smooth=3
lensmodel.plot.contourPlotAll(chains,lnprobs=lnprobs,inputPars=freePars,smooth=smooth,labels=labels,showPlot=False,filename=plotDir+"contours.pdf")

