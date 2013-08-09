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
xshear_sdss,yshear_sdss,errshear_sdss=lensmodel.io.convertCosmo(x1,y1,e1,redshift=redshift)

shearSN_sdss=np.median(yshear_sdss/errshear_sdss)

n_sdss=1.2 # N sources / sq. arcmin
n_des=12.
n_lsst=37.

A_sdss=0.9243 # Survey area in 10**4 deg**2
A_des=0.5
A_lsst=1.8

shearSN_des=shearSN_sdss * np.sqrt(n_des/n_sdss * A_des/A_sdss)
shearSN_lsst=shearSN_sdss * np.sqrt(n_lsst/n_sdss * A_lsst/A_sdss)

# Input model pars
logMstars=11.
logRstars=0.5
logMhalo=12.0
conc=5.
innerSlopeGNFW=1.
nuDutton=0.
AGnedin=-1.
wGnedin=-1.
inputPars=[logMstars, logRstars, logMhalo, conc, innerSlopeGNFW, nuDutton, AGnedin, wGnedin]
allLabels=np.array(["log(SM)","log(Reff)","log(Mh)","conc","inner slope","nuDutton","AGnedin","wGnedin"])
cenType="hernquist"
odType="critical"
delta=200.

priors=[[10.5,11.5],logRstars,[11.5,12.5],5.,1.,[-0.2,1.0],-1.,-1.]
ind=[0,2,5]
freePars=[inputPars[ii] for ii in ind]
labels=allLabels[ind]


# Now redo the radial binning with similar scale as sdss
Rmin=40. # kpc - note this is the geometric bin center, so measurements would need to go to radii 0.5*bin width smaller
Rmax=1000. # kpc - Use only bins at R<1Mpc
dlog10R=np.median(np.log10(xshear_sdss[1:]/xshear_sdss[:-1])) # log bin spacing
nBins=np.floor(np.log10(Rmax/Rmin)/dlog10R)+1
xshear=Rmin*10.**(np.arange(nBins) * dlog10R)

# Evaluate the model at these points - this will be the data vector
yshear=lensmodel.profiles.deltaSigma(inputPars,xshear)

# Scale error bars
shearSN=shearSN_lsst
errshear=yshear / shearSN

# Get magnification signal
magSN=0.5*shearSN # assume shear and mag have same radial dependence
xmag=xshear
ymag=lensmodel.profiles.sigma(inputPars,xmag)
errmag=ymag / magSN

# Run chains
nWalkers=2000
#nWalkers=100
nBurn=100
nSteps=1000
nThreads=40
#nThreads=8
seed=None

chains,lnprobs=lensmodel.fit.fitObs(priors,xshear,yshear,errshear,xmag,ymag,errmag,redshift=redshift,cenType=cenType,delta=delta,odType=odType,nWalkers=nWalkers,nBurn=nBurn,nSteps=nSteps,nThreads=nThreads,seed=seed)


suffix="alt3par_lsst_40_long"

lensmodel.io.writeRec(lensmodel.io.chainToRec(chains[0],lnprobs[0],labels=labels),dataDir+"chainM_{}.fits.gz".format(suffix))
lensmodel.io.writeRec(lensmodel.io.chainToRec(chains[1],lnprobs[1],labels=labels),dataDir+"chainS_{}.fits.gz".format(suffix))
lensmodel.io.writeRec(lensmodel.io.chainToRec(chains[2],lnprobs[2],labels=labels),dataDir+"chainSM_{}.fits.gz".format(suffix))


smooth=3
lensmodel.plot.contourPlotAll(chains,lnprobs=lnprobs,inputPars=freePars,smooth=smooth,labels=labels,showPlot=False,filename=plotDir+"contours_{}.pdf".format(suffix))

