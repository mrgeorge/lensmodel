#! env python

# Plot shear and mag profiles for the fiducial galaxy and cluster model with error bars

import matplotlib.pyplot as plt
import driver
import lensmodel
import numpy as np

survey=["sdss","lsst"]
xoffs=[0.96,1.04]
surveyColors=["gray","slategray"]
elw=3

target=["galaxy","cluster"]
targetLWs=[2,4]

Rmin=40.
magFrac=0.5


plotDir="/data/mgeorge/sdsslens/plots/"
plt.clf()
figsize=(7,4)
logSigLim=(1.,2.e3)

plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':12})
plt.rc('text', usetex=True)
plt.rc('axes',linewidth=1.5)

fig,axarr=plt.subplots(1,2,figsize=figsize)
fig.subplots_adjust(hspace=0,wspace=0.45)

starsColor="royalblue"
haloColor="red"
totalColor="black"
starsLS="--"
haloLS=":"
totalLS="-"

# Setup radial binning
Rmax=2000. # kpc - upper end of largest bin
xlim=(0.5*Rmin,2.*Rmax)
dlog10R=0.15
xbins,xshear=driver.getRadialBins(Rmin,Rmax,dlog10R)
dlog10Rfine=0.01
xbinsfine,xshearfine=driver.getRadialBins(xlim[0],xlim[1],dlog10Rfine)
xmag=xshear
xmagfine=xshearfine
yshearstars=np.zeros(len(xshearfine))
yshearhalo=np.zeros_like(yshearstars)
ymagstars=np.zeros_like(yshearstars)
ymaghalo=np.zeros_like(yshearstars)

for tt,lw in zip(target,targetLWs):
    # Input model pars
    n_lens,z_lens,inputPars,allLabels,allPlotLabels,cenType,odType,delta=driver.getModelPars(tt)
    logMstars, logRstars, logMhalo, conc, innerSlopeGNFW, nuDutton, AGnedin, wGnedin=inputPars

    # Evaluate the shear model at these points - this will be the shear data vector
    yshear=lensmodel.profiles.deltaSigma(inputPars,xshear)
    yshearfine=lensmodel.profiles.deltaSigma(inputPars,xshearfine,deltaSigmaStarsOut=yshearstars,deltaSigmaHaloOut=yshearhalo)
    
    # Get magnification signal
    ymag=lensmodel.profiles.sigma(inputPars,xmag)
    ymagfine=lensmodel.profiles.sigma(inputPars,xmagfine,sigmaStarsOut=ymagstars,sigmaHaloOut=ymaghalo)

    for ss,xoff,scol in zip(survey,xoffs,surveyColors):
        # Survey pars
        n_source,z_source,A_survey=driver.getSurveyPars(ss)

        # Estimate errors starting with shape noise for shear and using magFrac to scale to magnification errors
        sigma_shear=0.4 # shape noise
        Sigma_crit=1./lensmodel.profiles.cosmo.sigmacritinv(z_lens, z_source) # in Msun/pc**2
        xarea_kpc2=np.pi*(xbins[1:]**2 - xbins[:-1]**2)
        xarea_arcmin2=xarea_kpc2/(lensmodel.profiles.cosmo.Da(0.,z_lens) * 1000. * np.deg2rad(1./60.))**2 # denom is (kpc/arcmin)**2
        errshear = sigma_shear * Sigma_crit / np.sqrt(n_source * xarea_arcmin2 * n_lens * A_survey)
    
        # Scale error bars for magnification
        errmag=errshear / magFrac


        #PLOTS

        # Magnification
        plt.sca(axarr[0])
        lensmodel.plot.plotProfile(xmagfine,ymagfine,xlim=xlim,ylim=logSigLim,xlabel=r"$R$ (kpc)",ylabel=r"$\Sigma~($M$_{\odot}~$pc$^{-2})$",lw=lw,ls=totalLS)
        lensmodel.plot.plotProfile(xmagfine,ymagstars,color=starsColor,ls=starsLS,lw=lw)
        lensmodel.plot.plotProfile(xmagfine,ymaghalo,color=haloColor,ls=haloLS,lw=lw)
        plt.errorbar(xmag*xoff,ymag,errmag,ecolor=scol,fmt=None,elinewidth=elw)

        # Shear
        plt.sca(axarr[1])
        lensmodel.plot.plotProfile(xshearfine,yshearfine,xlim=xlim,ylim=logSigLim,xlabel=r"$R$ (kpc)",ylabel=r"$\Delta\Sigma~($M$_{\odot}~$pc$^{-2})$",lw=lw,ls=totalLS)
        lensmodel.plot.plotProfile(xshearfine,yshearstars,color=starsColor,ls=starsLS,lw=lw)
        lensmodel.plot.plotProfile(xshearfine,yshearhalo,color=haloColor,ls=haloLS,lw=lw)
        plt.errorbar(xshear*xoff,yshear,errshear,ecolor=scol,fmt=None,elinewidth=elw)

plt.show()
