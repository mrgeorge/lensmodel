#! env python

import lensmodel
import numpy as np
import sys

def main(survey, target, Rmin, magFrac, concPriorType, nThreads=8):
    """Forecast model constraints for a given lensing experiment

    Inputs:
        survey - string (sdss, des, lsst)
        target - string (galaxy, cluster)
        Rmin - minimum bin radius in kpc
        magFrac - ratio of errshear/errmag = (sigma_shear / sqrt(n_source_shear)) / (sigma_mag / sqrt(n_source_mag))
        concPriorType - how should concentration be treated in the fit (low, true, high, free)
            low=conc is fixed to conc - 1sigma
            true=conc is fixed to true value
            high=conc is fixed to conc + 1sigma
            free=conc is free with flat prior conc+/-2sigma
        nThreads - number of threads available (default 8)
    Outputs:
        Writes chains to dataDir and contour plot to plotDir
    
    """

    suffix="{}_{}_{}_{}_c{}".format(survey,target,Rmin,magFrac,concPriorType)
    dataDir="/data/mgeorge/sdsslens/data/"
    plotDir="/data/mgeorge/sdsslens/plots/"

    if(survey=="sdss"):
        n_source=1.2 # N shear sources / sq. arcmin
        z_source=1. # effective source redshift
        A_survey=9243. # survey area in deg**2
    elif(survey=="des"):
        n_source=12.
        z_source=1.
        A_survey=5000.
    elif(survey=="lsst"):
        n_source=37.
        z_source=1.
        A_survey=18000.
    else:
        raise ValueError(survey)

    # Input model pars
    if(target=="galaxy"):
        n_lens=1. # N lens galaxies / sq. deg
        z_lens=0.1
        logMstars=10.5
        logRstars=0.5
        logMhalo=12.0
        conc=10.
    elif(target=="cluster"):
        n_lens=1.e-3
        z_lens=0.1
        logMstars=11.
        logRstars=1.
        logMhalo=14.0
        conc=5.
    else:
        raise ValueError(target)
    innerSlopeGNFW=1.
    nuDutton=0.
    AGnedin=-1.
    wGnedin=-1.
    inputPars=[logMstars, logRstars, logMhalo, conc, innerSlopeGNFW, nuDutton, AGnedin, wGnedin]
    allLabels=np.array(["log(SM)","log(Reff)","log(Mh)","conc","inner slope","nuDutton","AGnedin","wGnedin"])
    cenType="hernquist"
    odType="critical"
    delta=200.


    # Choose free parameters and set priors
    freeInd=[0,2,5] # list of free pars - default is [logMstars, logMhalo, nuDutton] but conc can be added
    sigmaLogC=0.2 # base 10 lognormal scatter in c-M
    if(concPriorType=="low"):
        concPrior=10.**(np.log10(conc)-sigmaLogC)
    elif(concPriorType=="true"):
        concPrior=conc
    elif(concPriorType=="high"):
        concPrior=10.**(np.log10(conc)+sigmaLogC)
    elif(concPriorType=="free"):
        concPrior=[10.**(np.log10(conc) - 2.*sigmaLogC), 10.**(np.log10(conc) + 2.*sigmaLogC)]
        freeInd=[0,2,3,5]
    else:
        raise ValueError(concPriorType)
    freePars=[inputPars[ii] for ii in freeInd]
    labels=allLabels[freeInd]

    priors=[[logMstars-0.5,logMstars+0.5],logRstars,[logMhalo-0.5,logMhalo+0.5],concPrior,innerSlopeGNFW,[-0.2,1.0],AGnedin,wGnedin]


    # Setup radial binning
    Rmax=2000. # kpc - upper end of largest bin
    dlog10R=0.15
    nBins=np.floor(np.log10(Rmax/Rmin)/dlog10R)+1
    xbins=Rmin*10.**(np.arange(nBins) * dlog10R)
    xshear=10.**(np.log10(xbins[:-1])+0.5*dlog10R)

    # Evaluate the shear model at these points - this will be the shear data vector
    yshear=lensmodel.profiles.deltaSigma(inputPars,xshear)

    # Get magnification signal
    xmag=xshear
    ymag=lensmodel.profiles.sigma(inputPars,xmag)

    # Estimate errors starting with shape noise for shear and using magFrac to scale to magnification errors
    sigma_shear=0.4 # shape noise
    Sigma_crit=1./lensmodel.profiles.cosmo.sigmacritinv(z_lens, z_source) # in Msun/pc**2
    xarea_kpc2=np.pi*(xbins[1:]**2 - xbins[:-1]**2)
    xarea_arcmin2=xarea_kpc2/(lensmodel.profiles.cosmo.Da(0.,z_lens) * 1000. * np.deg2rad(1./60.))**2 # denom is (kpc/arcmin)**2
    errshear = sigma_shear / (Sigma_crit * np.sqrt(n_source * xarea_arcmin2 * n_lens * A_survey))

    # Scale error bars for magnification
    errmag=errshear / magFrac


    # Run chains
    nWalkers=200
    nBurn=10
    nSteps=30
    seed=None

    chains,lnprobs=lensmodel.fit.fitObs(priors,xshear,yshear,errshear,xmag,ymag,errmag,redshift=z_lens,cenType=cenType,delta=delta,odType=odType,nWalkers=nWalkers,nBurn=nBurn,nSteps=nSteps,nThreads=nThreads,seed=seed)


    # Save and plot outputs
    lensmodel.io.writeRec(lensmodel.io.chainToRec(chains[0],lnprobs[0],labels=labels),dataDir+"chainM_{}.fits.gz".format(suffix))
    lensmodel.io.writeRec(lensmodel.io.chainToRec(chains[1],lnprobs[1],labels=labels),dataDir+"chainS_{}.fits.gz".format(suffix))
    lensmodel.io.writeRec(lensmodel.io.chainToRec(chains[2],lnprobs[2],labels=labels),dataDir+"chainSM_{}.fits.gz".format(suffix))

    smooth=3
    lensmodel.plot.contourPlotAll(chains,lnprobs=lnprobs,inputPars=freePars,smooth=smooth,labels=labels,showPlot=False,filename=plotDir+"contours_{}.pdf".format(suffix))


if __name__ == "__main__":
# if called from command line

    if(len(sys.argv)!=7):
        print """Calling sequence: driver.py survey target Rmin magFrac concPriorType
        Inputs:
            survey - string (sdss, des, lsst)
            target - string (galaxy, cluster)
            Rmin - minimum bin radius in kpc
            magFrac - ratio of errshear/errmag = (sigma_shear / sqrt(n_source_shear)) / (sigma_mag / sqrt(n_source_mag))
            concPriorType - how should concentration be treated in the fit (low, true, high, free)
                low=conc is fixed to conc - 1sigma
                true=conc is fixed to true value
                high=conc is fixed to conc + 1sigma
                free=conc is free with flat prior conc+/-2sigma
            nThreads - number of threads available
            """
        raise ValueError(sys.argv)

    survey=sys.argv[1]
    target=sys.argv[2]
    Rmin=float(sys.argv[3])
    magFrac=float(sys.argv[4])
    concPriorType=sys.argv[5]
    nThreads=int(sys.argv[6])

    main(survey,target,Rmin,magFrac,concPriorType,nThreads=nThreads)
    
