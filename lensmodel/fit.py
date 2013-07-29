#! env python

import numpy as np
import emcee
import profiles

def simData():
    """Generate fake lensing profiles with some S/N."""
    pass


####
# Priors
####
def makeFlatPrior(range):
    return lambda x: priorFlat(x, range)

def priorFlat(arg, range):
    if((arg >= range[0]) & (arg < range[1])):
	return 0
    else:
	return np.Inf

def makeGaussPrior(mean, sigma):
    return lambda x: ((x-mean)/sigma)**2

def interpretPriors(priors):
    """Generate functions to evaluate priors and fix variables.

    Inputs:
        priors - a list or tuple with an entry for each of M parameters in the full model
            None - leave this variable completely free
            float - fix the variable to this value
            list[a,b] - flat prior between a and b
            tuple(a,b) - gaussian prior with mean a and stddev b

    Returns:
        priorFuncs - ndarray of length N (# of free parameters to fit)
        fixed - ndarray of length M; None for free pars, float for fixed pars
        guess - ndarray of length N with initial guesses
        guessScale - ndarray of length N with scale for range of initial guesses
    """

    guess=np.array([10.,0.1,100.,0.,0.]) # TO DO!
    guessScale=np.array([10.,0.3,50.,0.02,0.02]) # TO DO!
    nPars=len(guess)
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
		    priorMean=np.copy(prior[0])
		    priorSigma=np.copy(prior[1])
		    priorFuncs[ii]=makeGaussPrior(priorMean,priorSigma)

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

####
# Evaluate Likelihood
####
def lnProb(modelPars, priorFuncs, fixed, xshear, yshear, errshear, xmag, ymag, errmag):
    """Return ln(P(model|data)) = -0.5*chisq to evaluate likelihood surface.

    Take model parameters, priors, and data and compute chisq=sum[((model-data)/error)**2].
    Inputs:
        modelPars - ndarray of N model parameters to be fit (N<=M)
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
            model=profiles.deltaSigma(fullPars, xshear)
	    data=yshear
	    error=errshear
        elif((xshear is None) & (xmag is not None)): # Magnification only
	    model=profiles.sigma(fullPars, xmag)
	    data=ymag
	    error=errmag
        elif((xshear is not None) & (xmag is not None)): # Shear + Magnification
            shearmodel=profiles.deltaSigma(fullPars,xshear)
            magmodel=profiles.sigma(fullPars,xmag)
	    model=np.concatenate([shearmodel,magmodel])
	    data=np.concatenate([yshear,ymag])
	    error=np.concatenate([errshear,errmag])

	chisq_like=np.sum(((model-data)/error)**2)

    return -0.5*(chisq_like+chisq_prior)

####
# MCMC
####
def runMCMC(priors, xshear, yshear, errshear, xmag, ymag, errmag,nWalkers=2000,nBurn=50,nSteps=250,seed=None):
    """Call emcee and return sampler.

    See interpretPriors for format of priors array.
    """
    
    # SETUP PARS and PRIORS
    priorFuncs,fixed,guess,guessScale = interpretPriors(priors)
    nPars=len(guess)

    # RUN MCMC
    np.random.seed(seed)
    walkerStart=np.array([np.random.randn(nWalkers)*guessScale[ii]+guess[ii] for ii in xrange(nPars)]).T
    sampler=emcee.EnsembleSampler(nWalkers,nPars,lnProb,args=[priorFuncs, fixed, xshear, yshear, errshear, xmag, ymag, errmag])
    print "emcee burnin"
    pos, prob, state = sampler.run_mcmc(walkerStart,nBurn)
    sampler.reset()

    print "emcee running"
    sampler.run_mcmc(pos, nSteps)

    return sampler
