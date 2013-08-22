#! env python

import lensmodel
import numpy as np
import sys

def getSurveyPars(survey):
    if(survey=="sdss"):
        n_source=1.2 # N shear sources / sq. arcmin (Reyes 2012)
        z_source=1. # effective source redshift [actually mean=0.42, med=0.39 from N11]
        A_survey=9243. # survey area in deg**2 (Reyes 2012 DR8 - they used only 7131 for DR7 lenses)
    elif(survey=="des"):
        n_source=12. # DES proposal https://www.darkenergysurvey.org/reports/proposal-standalone.pdf
        z_source=1. # DES proposal uses zmed=0.68
        A_survey=5000. # DES proposal
    elif(survey=="lsst"):
        n_source=37. # Chang 2013
        z_source=1. # Chang 2013 zmed = 0.82 for fiducial case 1. in Table 2
        A_survey=18000. # Chang 2013
    else:
        raise ValueError(survey)

    return (n_source,z_source,A_survey)

def getModelPars(target):
    odType="critical"
    delta=200.

    if(target=="galaxy"):
        n_lens=(0.5*75086 / 7131.) # N lens galaxies / sq. deg - Schulz 2010 faint sample [half the full sample 75086, DR7 area 7131]
        z_lens=0.1 # Actually 0.11 from Schulz 2010 faint sample
        logMstars=np.log10(7.e10) # Schulz 2010 faint sample
        logRstars=np.log10(lensmodel.profiles.RdevTorHern(4.3/(1.+z_lens))) # Schulz 2010 faint sample, convert comoving Rdev to physical rHern
        logMhalo=np.log10(7.e12) # Schulz 2010 faint sample
        conc=9. # Schulz 2010 faint sample
    elif(target=="cluster"):
        n_lens=(1711+787+272+47.)/6670. # Top 4 richness bins from Johnston 2007. N from Table 1, Area: SDSS DR4 - Sheldon 2009 says they use a somewhat smaller area, but doesn't seem to specify
        z_lens=0.1 # ? probably more like 0.2 or 0.3
        logMstars=12. # rough estimate from Johnston 2007
        logRstars=1. # ?
        logMhalo=np.log10(95.96e12) # 4th highest richness bin from Johnston 2007
        conc=5.82 # 4th highester richness bin from Johnston 2007
    elif(target[:2]=="sm"): # e.g. sm9.5
        logMstars=float(target[2:]) # parse target name to get log10 stellar mass
        z_lens=0.1
        dz_lens=0.05 # redshift range for lens sample centered on z_lens
        dlog10SM=0.5 # log10 SM bin width centered on logMstars
        n_lens=lensmodel.profiles.liwhiteSMF(logMstars)*lensmodel.profiles.cosmo.V(z_lens-0.5*dz_lens,z_lens+0.5*dz_lens) * dlog10SM * (np.pi/180.)**2/(4.*np.pi) # Li & White give SMF in N/Mpc**3/dex, this converts to N/deg**2 given survey depth dz_lens
        logRstars=0.5 # TO DO - use mass-size relation?
        logMvir=lensmodel.profiles.stellarMassToHaloMass(logMstars, z_lens) # Behroozi gives virial mass
        cvir=lensmodel.profiles.concentration(10.**logMvir, z_lens, type="KlypinVir")
        logMhalo=lensmodel.profiles.convertHaloMass(logMvir,cvir,z_lens,delta,odType) # convert to specified overdensity using Hu & Kravtsov 2003
        conc=lensmodel.profiles.concentration(10.**logMhalo, z_lens, type="Klypin200")
        
    else:
        raise ValueError(target)
    innerSlopeGNFW=1.
    nuDutton=0.
    AGnedin=-1.
    wGnedin=-1.
    inputPars=[logMstars, logRstars, logMhalo, conc, innerSlopeGNFW, nuDutton, AGnedin, wGnedin]
    allLabels=np.array(["logSM","logrh","logMh","conc","inner slope","nuDutton","AGnedin","wGnedin"])
    allPlotLabels=np.array([r"log($M_{\star}$)",r"log($r_{\star}$)",r"log($M_h$)",r"$c$",r"$\beta$",r"$\nu$",r"$A_G$",r"$w_G$"])
    cenType="hernquist"

    return (n_lens,z_lens,inputPars,allLabels,allPlotLabels,cenType,odType,delta)

def getRadialBins(Rmin,Rmax,dlog10R):
    nBins=np.floor(np.log10(Rmax/Rmin)/dlog10R)+1
    xbins=Rmin*10.**(np.arange(nBins) * dlog10R)
    xshear=10.**(np.log10(xbins[:-1])+0.5*dlog10R)

    return (xbins,xshear)

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

    # Survey pars
    n_source,z_source,A_survey=getSurveyPars(survey)

    # Input model pars
    n_lens,z_lens,inputPars,allLabels,allPlotLabels,cenType,odType,delta=getModelPars(target)
    logMstars, logRstars, logMhalo, conc, innerSlopeGNFW, nuDutton, AGnedin, wGnedin=inputPars

    # Choose free parameters and set priors
    freeInd=[0,2,5] # list of free pars - default is [logMstars, logMhalo, nuDutton] but conc can be added
    sigmaLogC=0.2 # base 10 lognormal scatter in c-M
    if(concPriorType=="low"):
        concPrior=float(10.**(np.log10(conc)-sigmaLogC))
    elif(concPriorType=="true"):
        concPrior=float(conc)
    elif(concPriorType=="high"):
        concPrior=float(10.**(np.log10(conc)+sigmaLogC))
    elif(concPriorType=="free"):
        concPrior=[10.**(np.log10(conc) - 2.*sigmaLogC), 10.**(np.log10(conc) + 2.*sigmaLogC)]
        freeInd=[0,2,3,5]
    else:
        raise ValueError(concPriorType)
    freePars=[inputPars[ii] for ii in freeInd]
    labels=allLabels[freeInd]
    plotLabels=allPlotLabels[freeInd]

    priors=[[logMstars-0.5,logMstars+0.5],float(logRstars),[logMhalo-0.5,logMhalo+0.5],concPrior,float(innerSlopeGNFW),[-0.2,1.0],float(AGnedin),float(wGnedin)]


    # Setup radial binning
    Rmax=2000. # kpc - upper end of largest bin
    dlog10R=0.15
    xbins,xshear=getRadialBins(Rmin,Rmax,dlog10R)

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
    errshear = sigma_shear * Sigma_crit / np.sqrt(n_source * xarea_arcmin2 * n_lens * A_survey)

    # Scale error bars for magnification
    errmag=errshear / magFrac


    # Run chains
    nWalkers=2000
    nBurn=100
    nSteps=2000
    seed=7

    chains,lnprobs=lensmodel.fit.fitObs(priors,xshear,yshear,errshear,xmag,ymag,errmag,redshift=z_lens,cenType=cenType,delta=delta,odType=odType,nWalkers=nWalkers,nBurn=nBurn,nSteps=nSteps,nThreads=nThreads,seed=seed)


    # Save and plot outputs
    lensmodel.io.writeRec(lensmodel.io.chainToRec(chains[0],lnprobs[0],labels=labels),dataDir+"chainM_{}.fits.gz".format(suffix))
    lensmodel.io.writeRec(lensmodel.io.chainToRec(chains[1],lnprobs[1],labels=labels),dataDir+"chainS_{}.fits.gz".format(suffix))
    lensmodel.io.writeRec(lensmodel.io.chainToRec(chains[2],lnprobs[2],labels=labels),dataDir+"chainSM_{}.fits.gz".format(suffix))

    smooth=3
    lensmodel.plot.contourPlotAll(chains,lnprobs=lnprobs,inputPars=freePars,smooth=smooth,percentiles=[0.68,0.95],labels=plotLabels,showPlot=False,filename=plotDir+"contours_{}.pdf".format(suffix))


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
    
