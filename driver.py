#! env python

import lensmodel
import numpy as np

dataDir="../data/"
inFile=dataDir+"rebin.matt.10x.zebracut.src.sm_11_11.5.big.out"
x1,y1,e1=lensmodel.io.readData(inFile)
redshift=0.1
xshear,yshear,errshear=lensmodel.io.convertCosmo(x1,y1,e1,redshift=redshift)

fitSel=(xshear < 1000.)


# Input model pars

logMstars=11.2
logRstars=1.
logMhalo=np.log10(5.e13)
conc=5.
innerSlopeGNFW=1.
nuDutton=0.
AGnedin=-1.
wGnedin=-1.
inputPars=[logMstars, logRstars, logMhalo, conc, innerSlopeGNFW, nuDutton, AGnedin, wGnedin]

priors=[logMstars,logRstars,(logMhalo,1.),[2.,10.],innerSlopeGNFW,[-0.3,1.2],-1.,-1.]

yshear=lensmodel.profiles.deltaSigma(inputPars,xshear)
ymag=lensmodel.profiles.sigma(inputPars,xmag)
xmag=xshear
errmag=2.*errshear*(ymag/yshear) # say S/N (mag) = 0.5 * S/N (shear) at all radii


# Before running chains (which are slow), try plotting chisq as a function of each free par while holding others fixed.
# See what radii would be needed for good discrimination of nu par.


samplerSM=lensmodel.fit.runMCMC(priors, xshear[fitSel], yshear[fitSel], errshear[fitSel], xmag[fitSel], ymag[fitSel], errmag[fitSel],redshift=redshift,cenType="hernquist",delta=200.,odType="critical",nWalkers=2000,nBurn=50,nSteps=250,seed=None)
samplerS=lensmodel.fit.runMCMC(priors, xshear[fitSel], yshear[fitSel], errshear[fitSel], None, None, None,redshift=redshift,cenType="hernquist",delta=200.,odType="critical",nWalkers=2000,nBurn=50,nSteps=250,seed=None)
samplerM=lensmodel.fit.runMCMC(priors, None, None, None, xmag[fitSel], ymag[fitSel], errmag[fitSel],redshift=redshift,cenType="hernquist",delta=200.,odType="critical",nWalkers=2000,nBurn=50,nSteps=250,seed=None)
