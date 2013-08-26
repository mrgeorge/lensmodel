#! env python

import numpy as np
import emcee
import profiles
import scipy.stats

def simData():
    """Generate fake lensing profiles with some S/N."""
    pass


####
# Priors
####

# These generate functions used in lnProb to compute chisq_prior
def makeFlatPrior(range):
    return lambda x: priorFlat(x, range)

def priorFlat(arg, range):
    if((arg >= range[0]) & (arg < range[1])):
	return 0
    else:
	return np.Inf

def makeGaussPrior(mean, sigma):
    return lambda x: ((x-mean)/sigma)**2

def makeGaussTruncPrior(mean, sigma, range):
    return lambda x: ((x-mean)/sigma)**2 + priorFlat(x, range)
    
def interpretPriors(priors):
    """Generate functions to evaluate priors and fix variables.

    Inputs:
        priors - a list or tuple with an entry for each of M parameters in the full model
            None - leave this variable completely free
            float - fix the variable to this value
            list[a,b] - flat prior between a and b
            tuple(a,b) - gaussian prior with mean a and stddev b
            tuple(a,b,c,d) - gaussian prior with mean a and stddev b, truncated at c and d

    Returns:
        priorFuncs - ndarray of length N (# of free parameters to fit)
        fixed - ndarray of length M; None for free pars, float for fixed pars
        guess - ndarray of length N with initial guesses
        guessScale - ndarray of length N with scale for range of initial guesses
    """

    # Define initial guess and range for emcee
    # Recall pars=[logMstars, logRstars, logMhalo, conc, innerSlopeGNFW, nuDutton, AGnedin, wGnedin]
    guess=np.array([     10.,0.5,13.,5.0,1.0,0.8,1.6,1.0])
    guessScale=np.array([ 1.,0.5, 1.,2.0,0.2,0.3,0.3,0.3])
    nPars=len(guess) # Number of pars in FULL MODEL (some may get fixed and not be sent to emcee)

    fixed=np.repeat(None, nPars)
    priorFuncs=np.repeat(None, nPars)
    if(priors is not None):
	for ii in xrange(nPars):
	    prior=priors[ii]
	    # note: each of the assignments below needs to *copy* aspects of prior to avoid pointer overwriting
	    if(prior is not None):
		if((type(prior) is int) | (type(prior) is float)):
		# entry will be removed from list of pars and guess but value is still sent to evauluate function
		    fixVal=np.copy(prior)
		    fixed[ii]=fixVal
		elif(type(prior) is list):
		    priorRange=np.copy(prior)
		    priorFuncs[ii]=makeFlatPrior(priorRange)
		elif(type(prior) is tuple):
                    if(len(prior)==2):
                        priorMean=np.copy(prior[0])
                        priorSigma=np.copy(prior[1])
                        priorFuncs[ii]=makeGaussPrior(priorMean,priorSigma)
                    elif(len(prior)==4):
                        priorMean=np.copy(prior[0])
                        priorSigma=np.copy(prior[1])
                        priorRange=np.copy(prior[2:])
                        priorFuncs[ii]=makeGaussTruncPrior(priorMean,priorSigma,priorRange)
                    else:
                        raise ValueError(prior)
                else:
                    raise ValueError(ii,prior,type(prior))


    # remove fixed entries from list of pars to fit
    delarr=np.array([])
    for ii in xrange(nPars):
	if(fixed[ii] is not None):
            delarr=np.append(delarr,ii)
    if(len(delarr) > 0):
        guess=np.delete(guess,delarr)
        guessScale=np.delete(guessScale,delarr)
        priorFuncs=np.delete(priorFuncs,delarr)

    return (priorFuncs,fixed,guess,guessScale)

def testPriors():
    return [10.,0.5,(13.,1),[3.,8.],1.0,[-0.2,1.0],-1.,-1.]

####
# Evaluate Likelihood
####
def lnProb(pars, priors, xshear, yshear, errshear, xmag, ymag, errmag, redshift, cenType, delta, odType):
    """Return ln(P(model|data)) = -0.5*chisq to evaluate likelihood surface.

    Take model parameters, priors, and data and compute chisq=sum[((model-data)/error)**2].
    Inputs:
        pars - ndarray of N model parameters to be fit (N<=M)
        priorFuncs - ndarray of N functions
        fixed - ndarray of M full model parameters, either None (for free pars) or the value to fix to
            (priorFuncs and fixed are generally from interpretPriors)
        xshear - radial coordinates of shear profiles (None to ignore)
        yshear - delta sigma values
        errshear - errors on yshear
        xmag - radial coordinates of magnification profiles (None to ignore)
        ymag - sigma values
        errmag - errors on ymag

    Returns:
        lnP - a float (this is what emcee needs)
    """

    priorFuncs,fixed,guess,guessScale = interpretPriors(priors)

    # First evaluate the prior to see if this set of pars should be ignored
    chisq_prior=0.
    if(priorFuncs is not None):
	for ii in range(len(priorFuncs)):
	    func=priorFuncs[ii]
	    if(func is not None):
		chisq_prior+=func(pars[ii])
    if(chisq_prior == np.Inf):
	return -np.Inf

    # re-insert any fixed parameters into pars array
    nPars=len(fixed)
    fullPars=np.zeros(nPars)
    parsInd=0
    for ii in xrange(nPars):
	if(fixed[ii] is None):
	    fullPars[ii]=pars[parsInd]
	    parsInd+=1
        else:
            fullPars[ii]=fixed[ii]

    if((xshear is None) & (xmag is None)): # no data, only priors
	chisq_like=0.
    else:
	if((xshear is not None) & (xmag is None)): # Shear only
            model=profiles.deltaSigma(fullPars, xshear, redshift=redshift, cenType=cenType, delta=delta, odType=odType)
	    data=yshear
	    error=errshear
        elif((xshear is None) & (xmag is not None)): # Magnification only
	    model=profiles.sigma(fullPars, xmag, redshift=redshift, cenType=cenType, delta=delta, odType=odType)
	    data=ymag
	    error=errmag
        elif((xshear is not None) & (xmag is not None)): # Shear + Magnification
            shearmodel=profiles.deltaSigma(fullPars,xshear, redshift=redshift, cenType=cenType, delta=delta, odType=odType)
            magmodel=profiles.sigma(fullPars,xmag, redshift=redshift, cenType=cenType, delta=delta, odType=odType)
	    model=np.concatenate([shearmodel,magmodel])
	    data=np.concatenate([yshear,ymag])
	    error=np.concatenate([errshear,errmag])

	chisq_like=np.sum(((model-data)/error)**2)

    return -0.5*(chisq_like+chisq_prior)

####
# MCMC
####
def runMCMC(priors, xshear, yshear, errshear, xmag, ymag, errmag,redshift=0.,cenType="hernquist",delta=200.,odType="critical",nWalkers=2000,nBurn=50,nSteps=250,nThreads=1,seed=None):
    """Call emcee and return sampler.

    See interpretPriors for format of priors array.
    """
    
    # SETUP PARS and PRIORS
    priorFuncs,fixed,guess,guessScale = interpretPriors(priors)
    nPars=len(guess) # Number of free parameters to be fitted

    # RUN MCMC
    np.random.seed(seed)
    walkerStart=np.array([np.random.randn(nWalkers)*guessScale[ii]+guess[ii] for ii in xrange(nPars)]).T
    sampler=emcee.EnsembleSampler(nWalkers,nPars,lnProb,args=[priors, xshear, yshear, errshear, xmag, ymag, errmag, redshift, cenType, delta, odType],threads=nThreads)
    print "emcee burnin"
    pos, prob, state = sampler.run_mcmc(walkerStart,nBurn)
    sampler.reset()

    print "emcee running"
    sampler.run_mcmc(pos, nSteps)

    return sampler

def fitObs(priors,xshear,yshear,errshear,xmag,ymag,errmag,**kwargs):
    """Wrapper to runMCMC to compare chains with shear, magnification, and combined observables.
    Passes kwargs to vmapFit
    """

    print "Magnification"
    samplerM=runMCMC(priors,None,None,None,xmag,ymag,errmag,**kwargs)
    print("Mean acceptance fraction: {0:.3f}".format(np.mean(samplerM.acceptance_fraction)))
    print "Shear"
    samplerS=runMCMC(priors,xshear,yshear,errshear,None,None,None,**kwargs)
    print("Mean acceptance fraction: {0:.3f}".format(np.mean(samplerS.acceptance_fraction)))
    print "Combined"
    samplerSM=runMCMC(priors,xshear,yshear,errshear,xmag,ymag,errmag,**kwargs)
    print("Mean acceptance fraction: {0:.3f}".format(np.mean(samplerSM.acceptance_fraction)))
    
    flatchainM=samplerM.flatchain
    flatlnprobM=samplerM.flatlnprobability
    flatchainS=samplerS.flatchain
    flatlnprobS=samplerS.flatlnprobability
    flatchainSM=samplerSM.flatchain
    flatlnprobSM=samplerSM.flatlnprobability
    
    goodM=(flatlnprobM > -np.Inf)
    goodS=(flatlnprobS > -np.Inf)
    goodSM=(flatlnprobSM > -np.Inf)

    print "M finite fraction: {0:.3f}".format(float(len(goodM.nonzero()[0]))/len(flatlnprobM))
    print "S finite fraction: {0:.3f}".format(float(len(goodS.nonzero()[0]))/len(flatlnprobS))
    print "SM finite fraction: {0:.3f}".format(float(len(goodSM.nonzero()[0]))/len(flatlnprobSM))
    
    chains=[flatchainM[goodM], flatchainS[goodS], flatchainSM[goodSM]]
    lnprobs=[flatlnprobM[goodM],flatlnprobS[goodS],flatlnprobSM[goodSM]]
    return (chains,lnprobs)

####
# Chain statistics
####
def getMaxProb(chain,lnprob):
    maxP=(lnprob == np.max(lnprob)).nonzero()[0][0]
    return chain[maxP]

def getPeakKDE(chain,guess):
    if(len(chain.shape)==1):
        nPars=1
        kern=scipy.stats.gaussian_kde(chain)
        peakKDE=scipy.optimize.fmin(lambda x: -kern(x), guess,disp=False)
        return peakKDE
    else:
        nPars=chain.shape[1]
        peakKDE=np.zeros(nPars)
        for ii in range(nPars):
            kern=scipy.stats.gaussian_kde(chain[:,ii])
            peakKDE[ii]=scipy.optimize.fmin(lambda x: -kern(x), guess[ii],disp=False)
        return peakKDE

def getMedPost(chain):
# this is not a good estimator when posteriors are flat
    return np.median(chain,axis=0)

def get68(chain,opt="hw"):
# get half-width of 68% confidence range
# for a gaussian distribution, this is 1-sigma
    nSteps=len(chain)
    chainSort=np.sort(chain,axis=0)
    low68=chainSort[0.16*nSteps]
    high68=chainSort[0.84*nSteps]
    hw68=0.5*(high68-low68)
    if(opt=="hw"):
        return hw68
    elif(opt=="low"):
        return low68
    elif(opt=="high"):
        return high68
    elif(opt=="lowhigh"):
        return (low68,high68)
